#![allow(clippy::too_many_arguments)]

use extendr_api::prelude::*;

mod conversions;

fn ensure_logger() {
    use std::sync::Once;
    static INIT: Once = Once::new();
    INIT.call_once(|| {
        let _ = env_logger::try_init();
    });
}

/// Convert R nullable floats to Option<f64>, mapping NA to None.
fn rfloat_to_options(vals: &[Rfloat]) -> Vec<Option<f64>> {
    vals.iter()
        .map(|v| if v.is_na() { None } else { Some(v.inner()) })
        .collect()
}

/// Clamp intercept matrix diagonal to >= 1.0 (GenomicSEM convention).
fn clamp_i_ld_diagonal(i_mat: &mut faer::Mat<f64>) {
    for i in 0..i_mat.nrows() {
        if i_mat[(i, i)] < 1.0 {
            i_mat[(i, i)] = 1.0;
        }
    }
}

/// Run LDSC pipeline.
///
/// @param trait_files Character vector of .sumstats.gz file paths
/// @param sample_prev Numeric vector of sample prevalences (NA for continuous)
/// @param pop_prev Numeric vector of population prevalences (NA for continuous)
/// @param ld_dir Path to LD score directory
/// @param wld_dir Path to weight LD score directory (or same as ld_dir)
/// @param n_blocks Number of jackknife blocks (default 200)
/// @param chr Number of chromosomes (default 22)
/// @param chisq_max Maximum chi-square filter (NaN = auto)
/// @param stand Standardize output
/// @param select Variable selection method
/// @return JSON string with S, V, I matrices, N vector, and m
#[extendr]
fn ldsc_rust(
    trait_files: Vec<String>,
    sample_prev: Vec<Rfloat>,
    pop_prev: Vec<Rfloat>,
    ld_dir: &str,
    wld_dir: &str,
    n_blocks: i32,
    chr: i32,
    chisq_max: Rfloat,
    stand: bool,
    select: &str,
) -> String {
    ensure_logger();
    let n_blocks = n_blocks as usize;
    let chr_count = chr as usize;

    let sp = rfloat_to_options(&sample_prev);
    let pp = rfloat_to_options(&pop_prev);

    let trait_data = match gsem::io::gwas_reader::load_trait_data(&trait_files) {
        Ok(d) => d,
        Err(e) => return format!("{{\"error\": \"{e}\"}}"),
    };

    let ld_path = std::path::Path::new(ld_dir);
    let wld_path = std::path::Path::new(wld_dir);
    // Parse `select` parameter: "FALSE" → all, "ODD" → odd, "EVEN" → even, or comma-separated
    let chromosomes: Vec<usize> = match select.to_uppercase().as_str() {
        "FALSE" | "" => (1..=chr_count).collect(),
        "ODD" => (1..=chr_count).filter(|c| c % 2 == 1).collect(),
        "EVEN" => (1..=chr_count).filter(|c| c % 2 == 0).collect(),
        other => other
            .split(',')
            .filter_map(|s| s.trim().parse::<usize>().ok())
            .collect(),
    };
    let ld_data = match gsem::io::ld_reader::read_ld_scores(ld_path, wld_path, &chromosomes) {
        Ok(d) => d,
        Err(e) => return format!("{{\"error\": \"{e}\"}}"),
    };

    let ld_snps: Vec<String> = ld_data.records.iter().map(|r| r.snp.clone()).collect();
    let ld_scores: Vec<f64> = ld_data.records.iter().map(|r| r.l2).collect();

    let chisq_max_opt = if chisq_max.is_na() || chisq_max.inner().is_nan() {
        None
    } else {
        Some(chisq_max.inner())
    };

    let config = gsem_ldsc::LdscConfig {
        n_blocks,
        chisq_max: chisq_max_opt,
    };

    let n_pairs = trait_data.len() * (trait_data.len() + 1) / 2;
    let pb = indicatif::ProgressBar::new(n_pairs as u64);
    pb.set_style(
        indicatif::ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} trait pairs ({eta})",
        )
        .unwrap_or_else(|_| indicatif::ProgressStyle::default_bar()),
    );

    match gsem_ldsc::ldsc(
        &trait_data,
        &sp,
        &pp,
        &ld_scores,
        &ld_data.w_ld,
        &ld_snps,
        ld_data.total_m,
        &config,
        Some(&|| pb.inc(1)),
    ) {
        Ok(result) => {
            pb.finish_with_message("complete");
            if stand {
                conversions::ldsc_result_to_json_stand(&result)
            } else {
                conversions::ldsc_result_to_json(&result)
            }
        }
        Err(e) => {
            pb.finish();
            format!("{{\"error\": \"{e}\"}}")
        }
    }
}

/// Fit a user-specified SEM model.
///
/// @param covstruc_json LDSC result as JSON string
/// @param model lavaan-style model syntax
/// @param estimation Estimation method: "DWLS" or "ML"
/// @param std_lv Standardize latent variables
/// @param fix_resid Fix residual variances to be positive
/// @param imp_cov Return model-implied covariance matrix
/// @param q_factor Compute Q_Factor heterogeneity statistic
/// @return JSON string with converged, objective, and parameter estimates
#[extendr]
fn usermodel_rust(
    covstruc_json: &str,
    model: &str,
    estimation: &str,
    std_lv: bool,
    fix_resid: bool,
    imp_cov: bool,
    q_factor: bool,
    toler: Rfloat,
    cfi_calc: bool,
) -> String {
    ensure_logger();
    let ldsc_result = match conversions::json_to_ldsc_result(covstruc_json) {
        Some(r) => r,
        None => return "{\"error\": \"failed to parse covstruc JSON\"}".to_string(),
    };

    let grad_tol = if toler.is_na() || toler.inner().is_nan() || toler.inner() <= 0.0 {
        None
    } else {
        Some(toler.inner())
    };

    let mut pt = match gsem_sem::syntax::parse_model(model, std_lv) {
        Ok(pt) => pt,
        Err(e) => return format!("{{\"error\": \"model parse error: {e}\"}}"),
    };

    let k = ldsc_result.s.nrows();
    // Extract observed variable names from the model partable.
    // Observed variables = RHS of =~ (loadings), excluding latent factors and SNP.
    let obs_names: Vec<String> = {
        let latents: std::collections::HashSet<String> = pt.rows.iter()
            .filter(|r| r.op == gsem_sem::syntax::Op::Loading)
            .map(|r| r.lhs.clone())
            .collect();
        let mut names = Vec::new();
        for row in &pt.rows {
            if row.op == gsem_sem::syntax::Op::Loading {
                let name = &row.rhs;
                if !latents.contains(name) && name != "SNP" && !names.contains(name) {
                    names.push(name.clone());
                }
            }
        }
        if names.is_empty() {
            // No loadings found — try self-covariances (models without =~)
            for row in &pt.rows {
                if row.op == gsem_sem::syntax::Op::Covariance
                    && row.lhs == row.rhs
                    && !latents.contains(&row.lhs)
                    && row.lhs != "SNP"
                    && !names.contains(&row.lhs)
                {
                    names.push(row.lhs.clone());
                }
            }
        }
        if names.len() != k {
            return format!(
                "{{\"error\": \"model has {} observed variables but S matrix is {}x{}: obs=[{}]\"}}",
                names.len(), k, k, names.join(", ")
            );
        }
        names
    };
    let mut sem_model = gsem_sem::model::Model::from_partable(&pt, &obs_names);

    let kstar = k * (k + 1) / 2;
    let v_diag: Vec<f64> = (0..kstar).map(|i| ldsc_result.v[(i, i)]).collect();

    let fit = if estimation.to_uppercase() == "ML" {
        gsem_sem::estimator::fit_ml(&mut sem_model, &ldsc_result.s, 1000, grad_tol)
    } else {
        gsem_sem::estimator::fit_dwls(&mut sem_model, &ldsc_result.s, &v_diag, 1000, grad_tol)
    };

    // If model failed to converge and fix_resid is true, add lower bounds and retry
    let (fit, final_pt) = if !fit.converged && fix_resid {
        for row in &mut pt.rows {
            if row.op == gsem_sem::syntax::Op::Covariance
                && row.lhs == row.rhs
                && row.free > 0
                && row.lower_bound.is_none()
            {
                row.lower_bound = Some(0.0001);
            }
        }
        let mut model_retry = gsem_sem::model::Model::from_partable(&pt, &obs_names);
        let fit2 = if estimation.to_uppercase() == "ML" {
            gsem_sem::estimator::fit_ml(&mut model_retry, &ldsc_result.s, 1000, grad_tol)
        } else {
            gsem_sem::estimator::fit_dwls(&mut model_retry, &ldsc_result.s, &v_diag, 1000, grad_tol)
        };
        (fit2, pt.clone())
    } else {
        (fit, pt.clone())
    };

    let params_json: Vec<String> = final_pt
        .rows
        .iter()
        .enumerate()
        .filter(|(_, row)| row.free > 0)
        .map(|(i, row)| {
            let est = fit.params.get(i).copied().unwrap_or(0.0);
            format!(
                "{{\"lhs\":\"{}\",\"op\":\"{}\",\"rhs\":\"{}\",\"est\":{:.6}}}",
                row.lhs, row.op, row.rhs, est
            )
        })
        .collect();

    // Build the fitted model for implied_cov / fit indices / Q_Factor
    let mut fitted_model = gsem_sem::model::Model::from_partable(&final_pt, &obs_names);
    fitted_model.set_param_vec(&fit.params);
    let sigma_hat = fitted_model.implied_cov();

    // Compute fit indices
    let n_free = fitted_model.n_free();
    let df = kstar.saturating_sub(n_free);

    // Null model for CFI (variances only, no covariances)
    let (q_null, df_null) = if cfi_calc {
        let null_model_str: String = obs_names
            .iter()
            .map(|v| format!("{v} ~~ {v}"))
            .collect::<Vec<_>>()
            .join("\n");
        if let Ok(null_pt) = gsem_sem::syntax::parse_model(&null_model_str, false) {
            let mut null_model = gsem_sem::model::Model::from_partable(&null_pt, &obs_names);
            if estimation.to_uppercase() == "ML" {
                gsem_sem::estimator::fit_ml(&mut null_model, &ldsc_result.s, 1000, None);
            } else {
                gsem_sem::estimator::fit_dwls(&mut null_model, &ldsc_result.s, &v_diag, 1000, None);
            }
            let null_sigma = null_model.implied_cov();
            let null_df = kstar - k;
            let null_stats = gsem_sem::fit_indices::compute_fit(
                &ldsc_result.s, &null_sigma, &ldsc_result.v, null_df, k, None, None,
            );
            (Some(null_stats.chisq), Some(null_df))
        } else {
            (None, None)
        }
    } else {
        (None, None)
    };

    let model_fit = gsem_sem::fit_indices::compute_fit(
        &ldsc_result.s, &sigma_hat, &ldsc_result.v, df, n_free, q_null, df_null,
    );

    let fit_json = if cfi_calc {
        format!(
            ",\"chisq\":{:.4},\"df\":{},\"p_chisq\":{:.6e},\"aic\":{:.4},\"cfi\":{:.4},\"srmr\":{:.4}",
            model_fit.chisq, model_fit.df, model_fit.p_chisq, model_fit.aic, model_fit.cfi, model_fit.srmr
        )
    } else {
        format!(
            ",\"chisq\":{:.4},\"df\":{},\"p_chisq\":{:.6e},\"aic\":{:.4},\"srmr\":{:.4}",
            model_fit.chisq, model_fit.df, model_fit.p_chisq, model_fit.aic, model_fit.srmr
        )
    };

    // Model-implied covariance
    let imp_cov_json = if imp_cov {
        let rows: Vec<String> = (0..k)
            .map(|i| {
                let vals: Vec<String> = (0..k).map(|j| format!("{:.6}", sigma_hat[(i, j)])).collect();
                format!("[{}]", vals.join(","))
            })
            .collect();
        format!(",\"implied_cov\":[{}]", rows.join(","))
    } else {
        String::new()
    };

    // Q_Factor heterogeneity
    let q_factor_json = if q_factor {
        let factor_inds = gsem_sem::q_factor::factor_indicators(&final_pt, &obs_names);
        let q_results = gsem_sem::q_factor::compute_q_factor(
            &ldsc_result.s,
            &sigma_hat,
            &ldsc_result.v,
            &factor_inds,
        ).expect("q_factor: matrices must be square");
        let entries: Vec<String> = q_results
            .iter()
            .map(|q| {
                format!(
                    "{{\"factor1\":\"{}\",\"factor2\":\"{}\",\"Q_chisq\":{:.4},\"Q_df\":{},\"Q_p\":{:.6e}}}",
                    q.factor1, q.factor2, q.q_chisq, q.q_df, q.q_p
                )
            })
            .collect();
        format!(",\"Q_Factor\":[{}]", entries.join(","))
    } else {
        String::new()
    };

    format!(
        "{{\"converged\":{},\"objective\":{:.6},\"parameters\":[{}]{}{}{}}}",
        fit.converged,
        fit.objective,
        params_json.join(","),
        fit_json,
        imp_cov_json,
        q_factor_json,
    )
}

/// Munge GWAS summary statistics files.
///
/// @param files Character vector of GWAS file paths
/// @param hm3 Path to HapMap3 SNP list
/// @param trait_names Character vector of trait names
/// @param info_filter INFO score filter threshold
/// @param maf_filter MAF filter threshold
/// @param out_dir Output directory
/// @param n_override Sample size override (NaN = no override)
/// @param column_names_json Column name overrides as JSON string
/// @return Character vector of output file paths
#[extendr]
fn munge_rust(
    files: Vec<String>,
    hm3: &str,
    trait_names: Vec<String>,
    info_filter: f64,
    maf_filter: f64,
    out_dir: &str,
    n_override: Rfloat,
    column_names_json: &str,
) -> Vec<String> {
    ensure_logger();
    let hm3_path = std::path::Path::new(hm3);
    let reference = match gsem::munge::read_reference(hm3_path) {
        Ok(r) => r,
        Err(e) => {
            log::error!("Failed to read reference: {e}");
            return vec![];
        }
    };

    let n_opt = if n_override.is_na() || n_override.inner().is_nan() {
        None
    } else {
        Some(n_override.inner())
    };

    let col_overrides: Option<std::collections::HashMap<String, String>> =
        if column_names_json.is_empty() || column_names_json == "{}" {
            None
        } else {
            serde_json::from_str(column_names_json).ok()
        };

    let config = gsem::munge::MungeConfig {
        info_filter,
        maf_filter,
        n_override: n_opt,
        column_overrides: col_overrides,
    };

    let mut output_paths = Vec::new();
    for (i, file) in files.iter().enumerate() {
        let name = trait_names.get(i).map(|s| s.as_str()).unwrap_or("trait");
        let out_path = std::path::Path::new(out_dir).join(format!("{name}.sumstats.gz"));

        match gsem::munge::munge_and_write(
            std::path::Path::new(file),
            &reference,
            &config,
            &out_path,
        ) {
            Ok(()) => output_paths.push(out_path.to_string_lossy().to_string()),
            Err(e) => log::error!("Failed to munge {file}: {e}"),
        }
    }
    output_paths
}

/// Fit common factor model.
#[extendr]
fn commonfactor_rust(covstruc_json: &str, estimation: &str) -> String {
    ensure_logger();
    let ldsc_result = match conversions::json_to_ldsc_result(covstruc_json) {
        Some(r) => r,
        None => return "{\"error\": \"failed to parse covstruc JSON\"}".to_string(),
    };

    match gsem_sem::commonfactor::run_commonfactor(&ldsc_result.s, &ldsc_result.v, estimation) {
        Ok(result) => {
            let params: Vec<String> = result.parameters.iter().map(|p| {
                format!(
                    "{{\"lhs\":\"{}\",\"op\":\"{}\",\"rhs\":\"{}\",\"est\":{:.6},\"se\":{:.6},\"z\":{:.4},\"p\":{:.6e}}}",
                    p.lhs, p.op, p.rhs, p.est, p.se, p.z, p.p
                )
            }).collect();
            format!(
                "{{\"parameters\":[{}],\"chisq\":{:.4},\"df\":{},\"p_chisq\":{:.6e},\"aic\":{:.4},\"cfi\":{:.4},\"srmr\":{:.4}}}",
                params.join(","), result.fit.chisq, result.fit.df, result.fit.p_chisq,
                result.fit.aic, result.fit.cfi, result.fit.srmr
            )
        }
        Err(e) => format!("{{\"error\": \"{e}\"}}"),
    }
}

/// Merge GWAS summary statistics.
#[extendr]
fn sumstats_rust(
    files: Vec<String>,
    ref_dir: &str,
    trait_names: Vec<String>,
    info_filter: f64,
    maf_filter: f64,
    keep_indel: bool,
    out: &str,
    se_logit: Vec<i32>,
    ols: Vec<i32>,
    linprob: Vec<i32>,
    n_overrides_json: &str,
    keep_ambig: bool,
    betas_json: &str,
    direct_filter: bool,
) -> String {
    ensure_logger();
    let se_logit_bool: Vec<bool> = se_logit.iter().map(|&v| v != 0).collect();
    let ols_bool: Vec<bool> = ols.iter().map(|&v| v != 0).collect();
    let linprob_bool: Vec<bool> = linprob.iter().map(|&v| v != 0).collect();

    let n_overrides: Vec<Option<f64>> = if n_overrides_json.is_empty() || n_overrides_json == "[]" {
        vec![None; files.len()]
    } else {
        let parsed: Vec<serde_json::Value> =
            serde_json::from_str(n_overrides_json).unwrap_or_default();
        parsed
            .iter()
            .map(|v| v.as_f64())
            .collect()
    };

    // Parse betas: JSON object like {"trait1":"BETA_COL"} or empty
    let beta_overrides: Vec<Option<String>> = if betas_json.is_empty() || betas_json == "{}" {
        Vec::new()
    } else {
        // Parse as map of trait_name -> column_name, then align to trait order
        let betas_map: std::collections::HashMap<String, String> =
            serde_json::from_str(betas_json).unwrap_or_default();
        trait_names
            .iter()
            .map(|name| betas_map.get(name).cloned())
            .collect()
    };

    let config = gsem::sumstats::SumstatsConfig {
        info_filter,
        maf_filter,
        keep_indel,
        keep_ambig,
        se_logit: se_logit_bool,
        ols: ols_bool,
        linprob: linprob_bool,
        n_overrides,
        beta_overrides,
        direct_filter,
    };
    let file_refs: Vec<&std::path::Path> =
        files.iter().map(|p| std::path::Path::new(p.as_str())).collect();
    let out_path = std::path::Path::new(out);
    match gsem::sumstats::merge_sumstats(
        &file_refs,
        std::path::Path::new(ref_dir),
        &trait_names,
        &config,
        out_path,
    ) {
        Ok(n) => format!("{{\"path\":\"{}\",\"n_snps\":{}}}", out, n),
        Err(e) => format!("{{\"error\": \"{e}\"}}"),
    }
}

/// Run common factor GWAS.
#[extendr]
fn commonfactor_gwas_rust(
    covstruc_json: &str,
    sumstats_path: &str,
    gc: &str,
    estimation: &str,
    snp_se: Rfloat,
    smooth_check: bool,
    twas: bool,
) -> String {

    let ldsc_result = match conversions::json_to_ldsc_result(covstruc_json) {
        Some(r) => r,
        None => return "{\"error\": \"failed to parse covstruc JSON\"}".to_string(),
    };

    let k = ldsc_result.s.nrows();
    let gc_mode: gsem::gwas::gc_correction::GcMode =
        gc.parse().unwrap_or(gsem::gwas::gc_correction::GcMode::Standard);

    let mut i_ld = ldsc_result.i_mat.to_owned();
    clamp_i_ld_diagonal(&mut i_ld);

    let snp_se_opt = if snp_se.is_na() || snp_se.inner().is_nan() {
        None
    } else {
        Some(snp_se.inner())
    };

    let cf_config = gsem::gwas::common_factor::CommonFactorGwasConfig {
        estimation: estimation.to_string(),
        gc: gc_mode,
        snp_se: snp_se_opt,
        smooth_check,
    };

    if twas {
        let twas_data = match gsem::io::twas_reader::read_twas_sumstats(
            std::path::Path::new(sumstats_path),
        ) {
            Ok(d) => d,
            Err(e) => return format!("{{\"error\": \"{e}\"}}"),
        };

        let beta_gene: Vec<Vec<f64>> = twas_data.genes.iter().map(|g| g.beta.clone()).collect();
        let se_gene: Vec<Vec<f64>> = twas_data.genes.iter().map(|g| g.se.clone()).collect();
        let var_gene: Vec<f64> = twas_data.genes.iter().map(|g| g.hsq).collect();

        let results = gsem::gwas::common_factor::run_common_factor_gwas(
            &twas_data.trait_names,
            &ldsc_result.s,
            &ldsc_result.v,
            &i_ld,
            &beta_gene,
            &se_gene,
            &var_gene,
            &cf_config,
        );

        return conversions::twas_results_to_json(&results, &twas_data.genes);
    }

    let merged = match gsem::io::sumstats_reader::read_merged_sumstats(std::path::Path::new(
        sumstats_path,
    )) {
        Ok(m) => m,
        Err(e) => return format!("{{\"error\": \"{e}\"}}"),
    };

    let beta_snp: Vec<Vec<f64>> = merged.snps.iter().map(|s| s.beta.clone()).collect();
    let se_snp: Vec<Vec<f64>> = merged.snps.iter().map(|s| s.se.clone()).collect();
    let var_snp: Vec<f64> = merged
        .snps
        .iter()
        .map(|s| 2.0 * s.maf * (1.0 - s.maf))
        .collect();

    let results = gsem::gwas::common_factor::run_common_factor_gwas(
        &merged.trait_names,
        &ldsc_result.s,
        &ldsc_result.v,
        &i_ld,
        &beta_snp,
        &se_snp,
        &var_snp,
        &cf_config,
    );

    conversions::snp_results_to_json(&results, &merged.snps)
}

/// Run user-specified GWAS.
#[extendr]
fn user_gwas_rust(
    covstruc_json: &str,
    sumstats_path: &str,
    model: &str,
    estimation: &str,
    gc: &str,
    sub: &str,
    snp_se: Rfloat,
    smooth_check: bool,
    std_lv: bool,
    fix_measurement: bool,
    q_snp: bool,
    twas: bool,
) -> String {

    let ldsc_result = match conversions::json_to_ldsc_result(covstruc_json) {
        Some(r) => r,
        None => return "{\"error\": \"failed to parse covstruc JSON\"}".to_string(),
    };

    let k = ldsc_result.s.nrows();
    let gc_mode: gsem::gwas::gc_correction::GcMode =
        gc.parse().unwrap_or(gsem::gwas::gc_correction::GcMode::Standard);

    let mut i_ld = ldsc_result.i_mat.to_owned();
    clamp_i_ld_diagonal(&mut i_ld);

    let snp_se_opt = if snp_se.is_na() || snp_se.inner().is_nan() {
        None
    } else {
        Some(snp_se.inner())
    };

    // TWAS mode: read gene-level data instead of SNP-level
    if twas {
        let twas_data = match gsem::io::twas_reader::read_twas_sumstats(
            std::path::Path::new(sumstats_path),
        ) {
            Ok(d) => d,
            Err(e) => return format!("{{\"error\": \"{e}\"}}"),
        };

        let beta_gene: Vec<Vec<f64>> = twas_data.genes.iter().map(|g| g.beta.clone()).collect();
        let se_gene: Vec<Vec<f64>> = twas_data.genes.iter().map(|g| g.se.clone()).collect();
        let var_gene: Vec<f64> = twas_data.genes.iter().map(|g| g.hsq).collect();

        // Replace "SNP" with "Gene" in model syntax
        let twas_model = model.replace("SNP", "Gene");

        let config = gsem::gwas::user_gwas::UserGwasConfig {
            model: twas_model,
            estimation: estimation.to_string(),
            gc: gc_mode,
            std_lv,
            smooth_check,
            snp_se: snp_se_opt,
            fix_measurement,
            q_snp,
            snp_label: "Gene".to_string(),
            ..Default::default()
        };

        let mut results = gsem::gwas::user_gwas::run_user_gwas(
            &config,
            &ldsc_result.s,
            &ldsc_result.v,
            &i_ld,
            &beta_gene,
            &se_gene,
            &var_gene,
        );

        if !sub.is_empty() {
            filter_results_by_sub(&mut results, sub);
        }

        return conversions::twas_results_to_json(&results, &twas_data.genes);
    }

    // Standard SNP mode
    let merged = match gsem::io::sumstats_reader::read_merged_sumstats(std::path::Path::new(
        sumstats_path,
    )) {
        Ok(m) => m,
        Err(e) => return format!("{{\"error\": \"{e}\"}}"),
    };

    let beta_snp: Vec<Vec<f64>> = merged.snps.iter().map(|s| s.beta.clone()).collect();
    let se_snp: Vec<Vec<f64>> = merged.snps.iter().map(|s| s.se.clone()).collect();
    let var_snp: Vec<f64> = merged
        .snps
        .iter()
        .map(|s| 2.0 * s.maf * (1.0 - s.maf))
        .collect();

    let config = gsem::gwas::user_gwas::UserGwasConfig {
        model: model.to_string(),
        estimation: estimation.to_string(),
        gc: gc_mode,
        std_lv,
        smooth_check,
        snp_se: snp_se_opt,
        fix_measurement,
        q_snp,
        ..Default::default()
    };

    let mut results = gsem::gwas::user_gwas::run_user_gwas(
        &config,
        &ldsc_result.s,
        &ldsc_result.v,
        &i_ld,
        &beta_snp,
        &se_snp,
        &var_snp,
    );

    if !sub.is_empty() {
        filter_results_by_sub(&mut results, sub);
    }

    conversions::snp_results_to_json(&results, &merged.snps)
}

fn filter_results_by_sub(results: &mut [gsem::gwas::user_gwas::SnpResult], sub: &str) {
    let patterns: Vec<String> = sub
        .split(',')
        .map(|s| s.trim().replace(' ', ""))
        .filter(|s| !s.is_empty())
        .collect();
    if !patterns.is_empty() {
        for snp_result in results.iter_mut() {
            snp_result.params.retain(|p| {
                let key = format!("{}{}{}", p.lhs, p.op, p.rhs).replace(' ', "");
                patterns.iter().any(|pat| key == *pat)
            });
        }
    }
}

/// Parallel analysis to determine number of factors.
#[extendr]
fn pa_ldsc_rust(s_json: &str, v_json: &str, n_sim: i32, percentile: Rfloat, diag_only: bool) -> String {
    let s_mat = match conversions::json_to_mat(s_json) {
        Some(m) => m,
        None => return "{\"error\": \"failed to parse S matrix JSON\"}".to_string(),
    };

    let v_mat = match conversions::json_to_mat(v_json) {
        Some(m) => m,
        None => return "{\"error\": \"failed to parse V matrix JSON\"}".to_string(),
    };

    let p = if percentile.is_na() || percentile.inner().is_nan() {
        0.95
    } else {
        percentile.inner()
    };

    let result =
        gsem::stats::parallel_analysis::parallel_analysis(&s_mat, &v_mat, n_sim as usize, p, diag_only);

    let obs: Vec<String> = result.observed.iter().map(|v| format!("{v:.6}")).collect();
    let sim: Vec<String> = result
        .simulated_95
        .iter()
        .map(|v| format!("{v:.6}"))
        .collect();
    format!(
        "{{\"observed\":[{}],\"simulated_95\":[{}],\"n_factors\":{}}}",
        obs.join(","),
        sim.join(","),
        result.n_factors
    )
}

/// Auto-generate model syntax from factor loadings.
#[extendr]
fn write_model_rust(
    loadings_flat: Vec<f64>,
    n_rows: i32,
    n_cols: i32,
    names: Vec<String>,
    cutoff: f64,
    fix_resid: bool,
    bifactor: bool,
    mustload: bool,
    common: bool,
) -> String {
    let loadings = faer::Mat::from_fn(n_rows as usize, n_cols as usize, |i, j| {
        loadings_flat[i * n_cols as usize + j]
    });
    gsem_sem::write_model::write_model(&loadings, &names, cutoff, fix_resid, bifactor, mustload, common)
}

/// Compute model-implied genetic correlation matrix.
#[extendr]
fn rgmodel_rust(
    covstruc_json: &str,
    estimation: &str,
    model: &str,
    std_lv: bool,
    sub: &str,
) -> String {
    ensure_logger();
    let ldsc_result = match conversions::json_to_ldsc_result(covstruc_json) {
        Some(r) => r,
        None => return "{\"error\": \"failed to parse covstruc JSON\"}".to_string(),
    };

    let model_opt = if model.is_empty() { None } else { Some(model) };

    let result = if !sub.is_empty() {
        // Parse sub as comma-separated phenotype names, map to indices
        let k = ldsc_result.s.nrows();
        let all_names: Vec<String> = (0..k).map(|i| format!("V{}", i + 1)).collect();
        let sub_names: Vec<&str> = sub.split(',').map(|s| s.trim()).filter(|s| !s.is_empty()).collect();
        let sub_indices: Vec<usize> = sub_names.iter()
            .filter_map(|name| all_names.iter().position(|n| n == name))
            .collect();
        gsem_sem::rgmodel::run_rgmodel_sub(
            &ldsc_result.s, &ldsc_result.v, estimation, model_opt, std_lv, &sub_indices,
        )
    } else {
        gsem_sem::rgmodel::run_rgmodel_with_model(
            &ldsc_result.s, &ldsc_result.v, estimation, model_opt, std_lv,
        )
    };

    match result {
        Ok(result) => {
            let k = result.r.nrows();
            let r_rows: Vec<String> = (0..k)
                .map(|i| {
                    let row: Vec<String> =
                        (0..k).map(|j| format!("{:.6}", result.r[(i, j)])).collect();
                    format!("[{}]", row.join(","))
                })
                .collect();
            let kstar = k * (k + 1) / 2;
            let v_rows: Vec<String> = (0..kstar)
                .map(|i| {
                    let row: Vec<String> = (0..kstar)
                        .map(|j| format!("{:.6}", result.v_r[(i, j)]))
                        .collect();
                    format!("[{}]", row.join(","))
                })
                .collect();
            format!(
                "{{\"R\":[{}],\"V_R\":[{}]}}",
                r_rows.join(","),
                v_rows.join(",")
            )
        }
        Err(e) => format!("{{\"error\": \"{e}\"}}"),
    }
}

/// Run HDL analysis.
#[extendr]
fn hdl_rust(
    trait_files: Vec<String>,
    sample_prev: Vec<Rfloat>,
    pop_prev: Vec<Rfloat>,
    ld_path: &str,
    n_ref: f64,
    method: &str,
) -> String {
    ensure_logger();
    use gsem_ldsc::hdl::{HdlConfig, HdlMethod, HdlTraitData, LdPiece};

    let sp = rfloat_to_options(&sample_prev);
    let pp = rfloat_to_options(&pop_prev);

    // Read trait files
    let trait_data: Vec<HdlTraitData> = match gsem::io::gwas_reader::load_trait_data(&trait_files) {
        Ok(d) => d.into_iter().map(HdlTraitData::from).collect(),
        Err(e) => return format!("{{\"error\": \"{e}\"}}"),
    };

    let hdl_method = match method.to_lowercase().as_str() {
        "jackknife" => HdlMethod::Jackknife,
        _ => HdlMethod::Piecewise,
    };

    let config = HdlConfig {
        method: hdl_method,
        n_ref,
    };

    // Load LD pieces from text format directory
    let ld_dir = std::path::Path::new(ld_path);
    let ld_pieces = match gsem::io::hdl_reader::load_hdl_pieces(ld_dir) {
        Ok(p) => p,
        Err(e) => return format!("{{\"error\": \"{e}\"}}"),
    };

    match gsem_ldsc::hdl::hdl(&trait_data, &sp, &pp, &ld_pieces, &config) {
        Ok(result) => {
            let ldsc_compat = result.to_ldsc_result();
            conversions::ldsc_result_to_json(&ldsc_compat)
        }
        Err(e) => format!("{{\"error\": \"{e}\"}}"),
    }
}

/// Run stratified LD Score Regression.
#[extendr]
fn s_ldsc_rust(
    trait_files: Vec<String>,
    sample_prev: Vec<Rfloat>,
    pop_prev: Vec<Rfloat>,
    ld_dir: &str,
    wld_dir: &str,
    frq_dir: &str,
    n_blocks: i32,
    exclude_cont: bool,
) -> String {
    ensure_logger();
    let sp = rfloat_to_options(&sample_prev);
    let pp = rfloat_to_options(&pop_prev);

    // Read trait files
    let trait_data = match gsem::io::gwas_reader::load_trait_data(&trait_files) {
        Ok(d) => d,
        Err(e) => return format!("{{\"error\": \"{e}\"}}"),
    };

    let ld = std::path::Path::new(ld_dir);
    let wld = std::path::Path::new(wld_dir);
    let chromosomes: Vec<usize> = (1..=22).collect();

    // Read annotation LD scores
    let mut annot_data = match gsem_ldsc::annot_reader::read_annot_ld_scores(ld, wld, &chromosomes) {
        Ok(d) => d,
        Err(e) => return format!("{{\"error\": \"{e}\"}}"),
    };

    // Filter by frq files (MAF 0.05–0.95) if frq_dir is provided
    if !frq_dir.is_empty() {
        let frq_path = std::path::Path::new(frq_dir);
        match gsem_ldsc::annot_reader::read_frq_files(frq_path, &chromosomes) {
            Ok(frq_snps) => {
                gsem_ldsc::annot_reader::filter_annot_by_frq(&mut annot_data, &frq_snps);
            }
            Err(e) => log::warn!("Failed to read frq files: {e}"),
        }
    }

    // Filter out continuous annotations if exclude_cont is true
    if exclude_cont {
        let n_snps = annot_data.annot_ld.nrows();
        let n_annot = annot_data.annot_ld.ncols();
        let mut keep: Vec<bool> = vec![true; n_annot];
        for j in 0..n_annot {
            for i in 0..n_snps {
                let v = annot_data.annot_ld[(i, j)];
                if v != 0.0 && v != 1.0 {
                    keep[j] = false;
                    break;
                }
            }
        }
        let kept_indices: Vec<usize> = keep.iter().enumerate()
            .filter(|&(_, &k)| k).map(|(i, _)| i).collect();
        if kept_indices.len() < n_annot {
            let new_annot_ld = faer::Mat::from_fn(n_snps, kept_indices.len(), |i, j| {
                annot_data.annot_ld[(i, kept_indices[j])]
            });
            let new_names: Vec<String> = kept_indices.iter()
                .map(|&i| annot_data.annotation_names[i].clone()).collect();
            let new_m: Vec<f64> = kept_indices.iter()
                .map(|&i| annot_data.m_annot[i]).collect();
            annot_data.annot_ld = new_annot_ld;
            annot_data.annotation_names = new_names;
            annot_data.m_annot = new_m;
        }
    }

    let config = gsem_ldsc::stratified::StratifiedLdscConfig {
        n_blocks: n_blocks as usize,
        rm_flank: false,
        flank_kb: 500,
    };

    match gsem_ldsc::stratified::s_ldsc(
        &trait_data,
        &sp,
        &pp,
        &annot_data.annot_ld,
        &annot_data.w_ld,
        &annot_data.snps,
        &annot_data.annotation_names,
        &annot_data.m_annot,
        &config,
        Some(&annot_data.chr),
        Some(&annot_data.bp),
    ) {
        Ok(result) => match result.to_json_string() {
            Ok(json) => json,
            Err(e) => format!("{{\"error\": \"{e}\"}}"),
        },
        Err(e) => format!("{{\"error\": \"{e}\"}}"),
    }
}

/// Enrichment analysis using stratified LDSC results.
#[extendr]
fn enrich_rust(
    s_baseline_json: &str,
    s_annot_json: &str,
    v_annot_json: &str,
    annotation_names: Vec<String>,
    m_annot: Vec<f64>,
    m_total: f64,
) -> String {
    ensure_logger();
    let s_baseline = match conversions::json_to_mat(s_baseline_json) {
        Some(m) => m,
        None => return "{\"error\": \"failed to parse S_baseline JSON\"}".to_string(),
    };

    // Parse arrays of matrices from JSON: [[mat1_rows], [mat2_rows], ...]
    let s_annot_outer: Vec<Vec<Vec<f64>>> = match serde_json::from_str(s_annot_json) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\": \"failed to parse S_annot JSON: {e}\"}}"),
    };
    let s_annot: Vec<faer::Mat<f64>> = s_annot_outer
        .iter()
        .map(|rows| {
            let nr = rows.len();
            let nc = if nr > 0 { rows[0].len() } else { 0 };
            faer::Mat::from_fn(nr, nc, |i, j| rows[i][j])
        })
        .collect();

    let v_annot_outer: Vec<Vec<Vec<f64>>> = match serde_json::from_str(v_annot_json) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\": \"failed to parse V_annot JSON: {e}\"}}"),
    };
    let v_annot: Vec<faer::Mat<f64>> = v_annot_outer
        .iter()
        .map(|rows| {
            let nr = rows.len();
            let nc = if nr > 0 { rows[0].len() } else { 0 };
            faer::Mat::from_fn(nr, nc, |i, j| rows[i][j])
        })
        .collect();

    let result = gsem::stats::enrich::enrichment_test(
        &s_baseline,
        &s_annot,
        &v_annot,
        &annotation_names,
        &m_annot,
        m_total,
    );

    let entries: Vec<String> = result
        .annotations
        .iter()
        .enumerate()
        .map(|(i, name)| {
            format!(
                "{{\"annotation\":\"{}\",\"enrichment\":{:.6},\"se\":{:.6},\"p\":{:.6e}}}",
                name, result.enrichment[i], result.se[i], result.p[i]
            )
        })
        .collect();
    format!("[{}]", entries.join(","))
}

/// Simulate GWAS summary statistics.
#[extendr]
fn sim_ldsc_rust(
    s_json: &str,
    n_per_trait: Vec<f64>,
    ld_scores: Vec<f64>,
    m: f64,
    int_json: &str,
    r_pheno_json: &str,
    n_overlap: f64,
) -> String {
    ensure_logger();
    let s_mat = match conversions::json_to_mat(s_json) {
        Some(m) => m,
        None => return "{\"error\": \"failed to parse S matrix JSON\"}".to_string(),
    };

    let intercepts = if int_json.is_empty() || int_json == "null" {
        None
    } else {
        conversions::json_to_mat(int_json)
    };

    let r_pheno = if r_pheno_json.is_empty() || r_pheno_json == "null" {
        None
    } else {
        conversions::json_to_mat(r_pheno_json)
    };

    let config = gsem::stats::simulation::SimConfig {
        intercepts,
        r_pheno,
        n_overlap,
    };

    let result = gsem::stats::simulation::simulate_sumstats(
        &s_mat,
        &n_per_trait,
        &ld_scores,
        m,
        &config,
    );

    // result is Vec<Vec<f64>> (k traits x n_snps)
    let trait_json: Vec<String> = result
        .iter()
        .map(|trait_z| {
            let vals: Vec<String> = trait_z.iter().map(|v| format!("{v:.6}")).collect();
            format!("[{}]", vals.join(","))
        })
        .collect();
    format!("[{}]", trait_json.join(","))
}

/// Run multi-SNP analysis.
#[extendr]
fn multi_snp_rust(
    covstruc_json: &str,
    model: &str,
    estimation: &str,
    beta_json: &str,
    se_json: &str,
    var_snp: Vec<f64>,
    ld_matrix_json: &str,
    snp_names: Vec<String>,
    snp_se: Rfloat,
) -> String {
    ensure_logger();
    let ldsc_result = match conversions::json_to_ldsc_result(covstruc_json) {
        Some(r) => r,
        None => return "{\"error\": \"failed to parse covstruc JSON\"}".to_string(),
    };

    let beta: Vec<Vec<f64>> = match serde_json::from_str(beta_json) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\": \"failed to parse beta JSON: {e}\"}}"),
    };
    let se: Vec<Vec<f64>> = match serde_json::from_str(se_json) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\": \"failed to parse se JSON: {e}\"}}"),
    };
    let ld_matrix = match conversions::json_to_mat(ld_matrix_json) {
        Some(m) => m,
        None => return "{\"error\": \"failed to parse LD matrix JSON\"}".to_string(),
    };

    let snp_var_se = if snp_se.is_na() || snp_se.inner().is_nan() {
        None
    } else {
        Some(snp_se.inner())
    };

    let config = gsem::gwas::multi_snp::MultiSnpConfig {
        model: model.to_string(),
        estimation: estimation.to_string(),
        max_iter: 1000,
        snp_var_se,
    };

    let result = gsem::gwas::multi_snp::run_multi_snp(
        &config,
        &ldsc_result.s,
        &ldsc_result.v,
        &beta,
        &se,
        &var_snp,
        &ld_matrix,
        &snp_names,
    );

    let params: Vec<String> = result
        .params
        .iter()
        .map(|p| {
            format!(
                "{{\"lhs\":\"{}\",\"op\":\"{}\",\"rhs\":\"{}\",\"est\":{:.6},\"se\":{:.6},\"z\":{:.4},\"p\":{:.6e}}}",
                p.lhs, p.op, p.rhs, p.est, p.se, p.z_stat, p.p_value
            )
        })
        .collect();
    format!(
        "{{\"converged\":{},\"chisq\":{:.4},\"df\":{},\"params\":[{}]}}",
        result.converged, result.chisq, result.chisq_df, params.join(",")
    )
}

/// Run multi-gene analysis (reuses multi-SNP engine).
#[extendr]
fn multi_gene_rust(
    covstruc_json: &str,
    model: &str,
    estimation: &str,
    beta_json: &str,
    se_json: &str,
    var_gene: Vec<f64>,
    ld_matrix_json: &str,
    gene_names: Vec<String>,
    snp_se: Rfloat,
) -> String {
    ensure_logger();
    // multiGene is the same as multiSNP but with gene-level data
    multi_snp_rust(
        covstruc_json,
        model,
        estimation,
        beta_json,
        se_json,
        var_gene,
        ld_matrix_json,
        gene_names,
        snp_se,
    )
}

/// Run Generalized Least Squares regression.
#[extendr]
fn summary_gls_rust(
    x_json: &str,
    y: Vec<f64>,
    v_json: &str,
    intercept: bool,
) -> String {
    ensure_logger();
    let mut x_mat = match conversions::json_to_mat(x_json) {
        Some(m) => m,
        None => return "{\"error\": \"failed to parse X matrix JSON\"}".to_string(),
    };
    let v_mat = match conversions::json_to_mat(v_json) {
        Some(m) => m,
        None => return "{\"error\": \"failed to parse V matrix JSON\"}".to_string(),
    };

    // Add intercept column if requested
    if intercept {
        let n = x_mat.nrows();
        let p = x_mat.ncols();
        let mut x_new = faer::Mat::zeros(n, p + 1);
        for i in 0..n {
            x_new[(i, 0)] = 1.0;
            for j in 0..p {
                x_new[(i, j + 1)] = x_mat[(i, j)];
            }
        }
        x_mat = x_new;
    }

    match gsem::stats::gls::summary_gls(&x_mat, &y, &v_mat) {
        Some(result) => {
            let entries: Vec<String> = result
                .beta
                .iter()
                .enumerate()
                .map(|(i, b)| {
                    format!(
                        "{{\"beta\":{:.6},\"se\":{:.6},\"z\":{:.4},\"p\":{:.6e}}}",
                        b, result.se[i], result.z[i], result.p[i]
                    )
                })
                .collect();
            format!("[{}]", entries.join(","))
        }
        None => "{\"error\": \"GLS failed (singular matrix?)\"}".to_string(),
    }
}

extendr_module! {
    mod gsemr;
    fn ldsc_rust;
    fn usermodel_rust;
    fn munge_rust;
    fn commonfactor_rust;
    fn sumstats_rust;
    fn commonfactor_gwas_rust;
    fn user_gwas_rust;
    fn pa_ldsc_rust;
    fn write_model_rust;
    fn rgmodel_rust;
    fn hdl_rust;
    fn s_ldsc_rust;
    fn enrich_rust;
    fn sim_ldsc_rust;
    fn multi_snp_rust;
    fn multi_gene_rust;
    fn summary_gls_rust;
}
