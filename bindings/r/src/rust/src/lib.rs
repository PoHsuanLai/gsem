use extendr_api::prelude::*;

mod conversions;

/// Run LDSC pipeline.
///
/// @param trait_files Character vector of .sumstats.gz file paths
/// @param sample_prev Numeric vector of sample prevalences (NA for continuous)
/// @param pop_prev Numeric vector of population prevalences (NA for continuous)
/// @param ld_dir Path to LD score directory
/// @param wld_dir Path to weight LD score directory (or same as ld_dir)
/// @param n_blocks Number of jackknife blocks (default 200)
/// @return JSON string with S, V, I matrices, N vector, and m
#[extendr]
fn ldsc_rust(
    trait_files: Vec<String>,
    sample_prev: Vec<Rfloat>,
    pop_prev: Vec<Rfloat>,
    ld_dir: &str,
    wld_dir: &str,
    n_blocks: i32,
) -> String {
    let n_blocks = n_blocks as usize;

    let sp: Vec<Option<f64>> = sample_prev
        .iter()
        .map(|v| if v.is_na() { None } else { Some(v.inner()) })
        .collect();
    let pp: Vec<Option<f64>> = pop_prev
        .iter()
        .map(|v| if v.is_na() { None } else { Some(v.inner()) })
        .collect();

    let mut trait_data = Vec::new();
    for path in &trait_files {
        let path = std::path::Path::new(path);
        match gsem::io::gwas_reader::read_sumstats(path) {
            Ok(records) => {
                trait_data.push(gsem_ldsc::TraitSumstats {
                    snp: records.iter().map(|r| r.snp.clone()).collect(),
                    z: records.iter().map(|r| r.z).collect(),
                    n: records.iter().map(|r| r.n).collect(),
                    a1: records.iter().map(|r| r.a1.clone()).collect(),
                    a2: records.iter().map(|r| r.a2.clone()).collect(),
                });
            }
            Err(e) => return format!("{{\"error\": \"{e}\"}}"),
        }
    }

    let ld_path = std::path::Path::new(ld_dir);
    let wld_path = std::path::Path::new(wld_dir);
    let chromosomes: Vec<usize> = (1..=22).collect();
    let ld_data = match gsem::io::ld_reader::read_ld_scores(ld_path, wld_path, &chromosomes) {
        Ok(d) => d,
        Err(e) => return format!("{{\"error\": \"{e}\"}}"),
    };

    let ld_snps: Vec<String> = ld_data.records.iter().map(|r| r.snp.clone()).collect();
    let ld_scores: Vec<f64> = ld_data.records.iter().map(|r| r.l2).collect();

    let config = gsem_ldsc::LdscConfig {
        n_blocks,
        chisq_max: None,
    };

    match gsem_ldsc::ldsc(
        &trait_data,
        &sp,
        &pp,
        &ld_scores,
        &ld_data.w_ld,
        &ld_snps,
        ld_data.total_m,
        &config,
    ) {
        Ok(result) => conversions::ldsc_result_to_json(&result),
        Err(e) => format!("{{\"error\": \"{e}\"}}"),
    }
}

/// Fit a user-specified SEM model.
///
/// @param covstruc_json LDSC result as JSON string
/// @param model lavaan-style model syntax
/// @param estimation Estimation method: "DWLS" or "ML"
/// @return JSON string with converged, objective, and parameter estimates
#[extendr]
fn usermodel_rust(covstruc_json: &str, model: &str, estimation: &str) -> String {
    let ldsc_result = match conversions::json_to_ldsc_result(covstruc_json) {
        Some(r) => r,
        None => return "{\"error\": \"failed to parse covstruc JSON\"}".to_string(),
    };

    let pt = match gsem_sem::syntax::parse_model(model, false) {
        Ok(pt) => pt,
        Err(e) => return format!("{{\"error\": \"model parse error: {e}\"}}"),
    };

    let k = ldsc_result.s.nrows();
    let obs_names: Vec<String> = (0..k).map(|i| format!("V{}", i + 1)).collect();
    let mut sem_model = gsem_sem::model::Model::from_partable(&pt, &obs_names);

    let kstar = k * (k + 1) / 2;
    let v_diag: Vec<f64> = (0..kstar).map(|i| ldsc_result.v[(i, i)]).collect();

    let fit = if estimation.to_uppercase() == "ML" {
        gsem_sem::estimator::fit_ml(&mut sem_model, &ldsc_result.s, 1000)
    } else {
        gsem_sem::estimator::fit_dwls(&mut sem_model, &ldsc_result.s, &v_diag, 1000)
    };

    let params_json: Vec<String> = pt
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

    format!(
        "{{\"converged\":{},\"objective\":{:.6},\"parameters\":[{}]}}",
        fit.converged,
        fit.objective,
        params_json.join(",")
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
/// @return Character vector of output file paths
#[extendr]
fn munge_rust(
    files: Vec<String>,
    hm3: &str,
    trait_names: Vec<String>,
    info_filter: f64,
    maf_filter: f64,
    out_dir: &str,
) -> Vec<String> {
    let hm3_path = std::path::Path::new(hm3);
    let reference = match gsem::munge::read_reference(hm3_path) {
        Ok(r) => r,
        Err(e) => {
            log::error!("Failed to read reference: {e}");
            return vec![];
        }
    };

    let config = gsem::munge::MungeConfig {
        info_filter,
        maf_filter,
        n_override: None,
        column_overrides: None,
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
) -> String {
    let config = gsem::sumstats::SumstatsConfig {
        info_filter,
        maf_filter,
        keep_indel,
        ..Default::default()
    };
    let file_refs: Vec<&std::path::Path> = files.iter().map(|p| std::path::Path::new(p.as_str())).collect();
    let out_path = std::path::Path::new(out);
    match gsem::sumstats::merge_sumstats(&file_refs, std::path::Path::new(ref_dir), &trait_names, &config, out_path) {
        Ok(n) => format!("{{\"path\":\"{}\",\"n_snps\":{}}}", out, n),
        Err(e) => format!("{{\"error\": \"{e}\"}}"),
    }
}

/// Run common factor GWAS.
#[extendr]
fn commonfactor_gwas_rust(covstruc_json: &str, sumstats_path: &str, gc: &str) -> String {
    let ldsc_result = match conversions::json_to_ldsc_result(covstruc_json) {
        Some(r) => r,
        None => return "{\"error\": \"failed to parse covstruc JSON\"}".to_string(),
    };

    let merged = match gsem::io::sumstats_reader::read_merged_sumstats(std::path::Path::new(sumstats_path)) {
        Ok(m) => m,
        Err(e) => return format!("{{\"error\": \"{e}\"}}"),
    };

    let k = ldsc_result.s.nrows();
    let gc_mode: gsem::gwas::gc_correction::GcMode = gc.parse().unwrap_or(gsem::gwas::gc_correction::GcMode::Standard);
    let beta_snp: Vec<Vec<f64>> = merged.snps.iter().map(|s| s.beta.clone()).collect();
    let se_snp: Vec<Vec<f64>> = merged.snps.iter().map(|s| s.se.clone()).collect();
    let var_snp: Vec<f64> = merged.snps.iter().map(|s| 2.0 * s.maf * (1.0 - s.maf)).collect();

    let mut i_ld = ldsc_result.i_mat.to_owned();
    for i in 0..k {
        if i_ld[(i, i)] < 1.0 { i_ld[(i, i)] = 1.0; }
    }

    let results = gsem::gwas::common_factor::run_common_factor_gwas(
        &merged.trait_names, &ldsc_result.s, &ldsc_result.v, &i_ld,
        &beta_snp, &se_snp, &var_snp, gc_mode,
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
) -> String {
    let ldsc_result = match conversions::json_to_ldsc_result(covstruc_json) {
        Some(r) => r,
        None => return "{\"error\": \"failed to parse covstruc JSON\"}".to_string(),
    };

    let merged = match gsem::io::sumstats_reader::read_merged_sumstats(std::path::Path::new(sumstats_path)) {
        Ok(m) => m,
        Err(e) => return format!("{{\"error\": \"{e}\"}}"),
    };

    let k = ldsc_result.s.nrows();
    let gc_mode: gsem::gwas::gc_correction::GcMode = gc.parse().unwrap_or(gsem::gwas::gc_correction::GcMode::Standard);
    let beta_snp: Vec<Vec<f64>> = merged.snps.iter().map(|s| s.beta.clone()).collect();
    let se_snp: Vec<Vec<f64>> = merged.snps.iter().map(|s| s.se.clone()).collect();
    let var_snp: Vec<f64> = merged.snps.iter().map(|s| 2.0 * s.maf * (1.0 - s.maf)).collect();

    let mut i_ld = ldsc_result.i_mat.to_owned();
    for i in 0..k {
        if i_ld[(i, i)] < 1.0 { i_ld[(i, i)] = 1.0; }
    }

    let config = gsem::gwas::user_gwas::UserGwasConfig {
        model: model.to_string(),
        estimation: estimation.to_string(),
        gc: gc_mode,
        ..Default::default()
    };

    let results = gsem::gwas::user_gwas::run_user_gwas(
        &config, &ldsc_result.s, &ldsc_result.v, &i_ld,
        &beta_snp, &se_snp, &var_snp,
    );

    conversions::snp_results_to_json(&results, &merged.snps)
}

/// Parallel analysis to determine number of factors.
#[extendr]
fn pa_ldsc_rust(covstruc_json: &str, n_sim: i32) -> String {
    let ldsc_result = match conversions::json_to_ldsc_result(covstruc_json) {
        Some(r) => r,
        None => return "{\"error\": \"failed to parse covstruc JSON\"}".to_string(),
    };

    let result = gsem::stats::parallel_analysis::parallel_analysis(
        &ldsc_result.s, &ldsc_result.v, n_sim as usize,
    );

    let obs: Vec<String> = result.observed.iter().map(|v| format!("{v:.6}")).collect();
    let sim: Vec<String> = result.simulated_95.iter().map(|v| format!("{v:.6}")).collect();
    format!(
        "{{\"observed\":[{}],\"simulated_95\":[{}],\"n_factors\":{}}}",
        obs.join(","), sim.join(","), result.n_factors
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
) -> String {
    let loadings = faer::Mat::from_fn(n_rows as usize, n_cols as usize, |i, j| {
        loadings_flat[i * n_cols as usize + j]
    });
    gsem_sem::write_model::write_model(&loadings, &names, cutoff, fix_resid, bifactor)
}

/// Compute model-implied genetic correlation matrix.
#[extendr]
fn rgmodel_rust(covstruc_json: &str, estimation: &str) -> String {
    let ldsc_result = match conversions::json_to_ldsc_result(covstruc_json) {
        Some(r) => r,
        None => return "{\"error\": \"failed to parse covstruc JSON\"}".to_string(),
    };

    match gsem_sem::rgmodel::run_rgmodel(&ldsc_result.s, &ldsc_result.v, estimation) {
        Ok(result) => {
            let k = result.r.nrows();
            let r_rows: Vec<String> = (0..k).map(|i| {
                let row: Vec<String> = (0..k).map(|j| format!("{:.6}", result.r[(i, j)])).collect();
                format!("[{}]", row.join(","))
            }).collect();
            let kstar = k * (k + 1) / 2;
            let v_rows: Vec<String> = (0..kstar).map(|i| {
                let row: Vec<String> = (0..kstar).map(|j| format!("{:.6}", result.v_r[(i, j)])).collect();
                format!("[{}]", row.join(","))
            }).collect();
            format!("{{\"R\":[{}],\"V_R\":[{}]}}", r_rows.join(","), v_rows.join(","))
        }
        Err(e) => format!("{{\"error\": \"{e}\"}}"),
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
}
