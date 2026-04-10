#![allow(clippy::too_many_arguments)]

use extendr_api::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

mod conversions;

/// Custom log backend that routes messages to R's message stream.
struct RLogger;

static RLOGGER: RLogger = RLogger;

impl log::Log for RLogger {
    fn enabled(&self, metadata: &log::Metadata) -> bool {
        metadata.level() <= log::Level::Info
    }

    fn log(&self, record: &log::Record) {
        if self.enabled(record.metadata()) {
            reprintln!("[{}] {}", record.level(), record.args());
        }
    }

    fn flush(&self) {}
}

fn ensure_logger() {
    use std::sync::Once;
    static INIT: Once = Once::new();
    INIT.call_once(|| {
        let _ = log::set_logger(&RLOGGER).map(|()| log::set_max_level(log::LevelFilter::Info));
    });
}

/// Progress counter for reporting to R console at reasonable intervals.
struct RProgress {
    counter: AtomicUsize,
    total: usize,
    interval: usize,
}

impl RProgress {
    fn new(total: usize) -> Self {
        let interval = if total <= 20 { 1 } else { (total / 20).max(1) };
        Self {
            counter: AtomicUsize::new(0),
            total,
            interval,
        }
    }

    fn callback(&self) -> impl Fn() + Sync + '_ {
        move || {
            let done = self.counter.fetch_add(1, Ordering::Relaxed) + 1;
            if done == self.total || done.is_multiple_of(self.interval) {
                reprintln!("Progress: {done}/{}", self.total);
            }
        }
    }
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

/// Convert an R integer to the `num_threads` value the in-process API
/// expects: NA or non-positive → `None` (use rayon default = all cores),
/// positive → `Some(n)`. R wrappers should pass `NA_integer_` when the
/// caller wants the default and an explicit positive value otherwise.
fn rint_to_num_threads(n: Rint) -> Option<usize> {
    if n.is_na() {
        return None;
    }
    let v = n.inner();
    if v > 0 { Some(v as usize) } else { None }
}

/// Run LDSC pipeline.
///
/// Returns a named R list with `s`, `v`, `i_mat`, `n_vec`, `m` — plus
/// `s_stand` and `v_stand` when `stand=TRUE`.
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
    num_threads: Rint,
) -> List {
    ensure_logger();
    let n_blocks = n_blocks as usize;
    let chr_count = chr as usize;

    let sp = rfloat_to_options(&sample_prev);
    let pp = rfloat_to_options(&pop_prev);

    let trait_data = match gsem::io::gwas_reader::load_trait_data(&trait_files) {
        Ok(d) => d,
        Err(e) => return conversions::error_list(e.to_string()),
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
        Err(e) => return conversions::error_list(e.to_string()),
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
        num_threads: rint_to_num_threads(num_threads),
    };

    let n_pairs = trait_data.len() * (trait_data.len() + 1) / 2;
    let progress = RProgress::new(n_pairs);
    let cb = progress.callback();

    match gsem_ldsc::ldsc(
        &trait_data,
        &sp,
        &pp,
        &ld_scores,
        &ld_data.w_ld,
        &ld_snps,
        ld_data.total_m,
        &config,
        Some(&cb),
    ) {
        Ok(result) => {
            if stand {
                conversions::ldsc_result_to_list_stand(&result)
            } else {
                conversions::ldsc_result_to_list(&result)
            }
        }
        Err(e) => conversions::error_list(e.to_string()),
    }
}

/// Fit a user-specified SEM model.
///
/// Takes `covstruc` as a named R list (`s`, `v`, `i_mat`, `n_vec`, `m`) and
/// returns a named R list containing `converged`, `objective`, `parameters`
/// (columnar list), fit indices, and optional `implied_cov` / `Q_Factor`.
#[extendr]
fn usermodel_rust(
    covstruc: List,
    model: &str,
    estimation: &str,
    std_lv: bool,
    fix_resid: bool,
    imp_cov: bool,
    q_factor: bool,
    toler: Rfloat,
    cfi_calc: bool,
) -> List {
    ensure_logger();
    let ldsc_result = match conversions::list_to_ldsc_result(&covstruc) {
        Ok(r) => r,
        Err(e) => return conversions::error_list(e.to_string()),
    };

    let grad_tol = if toler.is_na() || toler.inner().is_nan() || toler.inner() <= 0.0 {
        None
    } else {
        Some(toler.inner())
    };

    let mut pt = match gsem_sem::syntax::parse_model(model, std_lv) {
        Ok(pt) => pt,
        Err(e) => return conversions::error_list(format!("model parse error: {e}")),
    };

    let k = ldsc_result.s.nrows();
    // Extract observed variable names from the model partable.
    // Observed variables = RHS of =~ (loadings), excluding latent factors and SNP.
    let obs_names: Vec<String> = {
        let latents: std::collections::HashSet<String> = pt
            .rows
            .iter()
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
            return conversions::error_list(format!(
                "model has {} observed variables but S matrix is {}x{}: obs=[{}]",
                names.len(),
                k,
                k,
                names.join(", ")
            ));
        }
        names
    };
    let mut sem_model = gsem_sem::model::Model::from_partable(&pt, &obs_names);

    let kstar = k * (k + 1) / 2;
    let v_diag: Vec<f64> = (0..kstar).map(|i| ldsc_result.v[(i, i)]).collect();

    let est_method = gsem_sem::EstimationMethod::from_str_lossy(estimation);
    let fit = match est_method {
        gsem_sem::EstimationMethod::Ml => {
            gsem_sem::estimator::fit_ml(&mut sem_model, &ldsc_result.s, 1000, grad_tol)
        }
        gsem_sem::EstimationMethod::Dwls => {
            gsem_sem::estimator::fit_dwls(&mut sem_model, &ldsc_result.s, &v_diag, 1000, grad_tol)
        }
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
        let fit2 = match est_method {
            gsem_sem::EstimationMethod::Ml => {
                gsem_sem::estimator::fit_ml(&mut model_retry, &ldsc_result.s, 1000, grad_tol)
            }
            gsem_sem::EstimationMethod::Dwls => gsem_sem::estimator::fit_dwls(
                &mut model_retry,
                &ldsc_result.s,
                &v_diag,
                1000,
                grad_tol,
            ),
        };
        (fit2, pt.clone())
    } else {
        (fit, pt.clone())
    };

    // Build a columnar parameter list from the free rows. We stash the
    // point estimate into each row so we can reuse the `ParamEstimate`
    // builder; se/z/p are unavailable for this path, so we fill them with
    // NaN (which maps to R NA via the default Doubles path).
    let free_rows: Vec<_> = final_pt.rows.iter().filter(|r| r.free > 0).collect();
    let (lhs_col, op_col, rhs_col, est_col): (Vec<String>, Vec<String>, Vec<String>, Vec<f64>) = {
        let mut lhs = Vec::with_capacity(free_rows.len());
        let mut op = Vec::with_capacity(free_rows.len());
        let mut rhs = Vec::with_capacity(free_rows.len());
        let mut est = Vec::with_capacity(free_rows.len());
        for (i, row) in free_rows.iter().enumerate() {
            lhs.push(row.lhs.clone());
            op.push(row.op.to_string());
            rhs.push(row.rhs.clone());
            est.push(fit.params.get(i).copied().unwrap_or(0.0));
        }
        (lhs, op, rhs, est)
    };
    let parameters = list!(lhs = lhs_col, op = op_col, rhs = rhs_col, est = est_col);

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
            match est_method {
                gsem_sem::EstimationMethod::Ml => {
                    gsem_sem::estimator::fit_ml(&mut null_model, &ldsc_result.s, 1000, None);
                }
                gsem_sem::EstimationMethod::Dwls => {
                    gsem_sem::estimator::fit_dwls(
                        &mut null_model,
                        &ldsc_result.s,
                        &v_diag,
                        1000,
                        None,
                    );
                }
            }
            let null_sigma = null_model.implied_cov();
            let null_df = kstar - k;
            let null_stats = gsem_sem::fit_indices::compute_fit(
                &ldsc_result.s,
                &null_sigma,
                &ldsc_result.v,
                null_df,
                k,
                None,
                None,
            );
            (Some(null_stats.chisq), Some(null_df))
        } else {
            (None, None)
        }
    } else {
        (None, None)
    };

    let model_fit = gsem_sem::fit_indices::compute_fit(
        &ldsc_result.s,
        &sigma_hat,
        &ldsc_result.v,
        df,
        n_free,
        q_null,
        df_null,
    );

    // Assemble the output list. We build a Vec of (name, Robj) pairs and
    // let `List::from_pairs` turn them into a named list — this keeps the
    // conditional CFI / implied_cov / Q_Factor fields readable.
    let mut pairs: Vec<(&'static str, Robj)> = vec![
        ("converged", fit.converged.into_robj()),
        ("objective", fit.objective.into_robj()),
        ("parameters", parameters.into_robj()),
        ("chisq", model_fit.chisq.into_robj()),
        ("df", (model_fit.df as i32).into_robj()),
        ("p_chisq", model_fit.p_chisq.into_robj()),
        ("aic", model_fit.aic.into_robj()),
        ("srmr", model_fit.srmr.into_robj()),
    ];
    if cfi_calc {
        pairs.push(("cfi", model_fit.cfi.into_robj()));
    }

    if imp_cov {
        pairs.push((
            "implied_cov",
            conversions::mat_to_rmatrix(&sigma_hat).into_robj(),
        ));
    }

    if q_factor {
        let factor_inds = gsem_sem::q_factor::factor_indicators(&final_pt, &obs_names);
        let q_results = gsem_sem::q_factor::compute_q_factor(
            &ldsc_result.s,
            &sigma_hat,
            &ldsc_result.v,
            &factor_inds,
        )
        .expect("q_factor: matrices must be square");
        let mut factor1 = Vec::with_capacity(q_results.len());
        let mut factor2 = Vec::with_capacity(q_results.len());
        let mut q_chisq = Vec::with_capacity(q_results.len());
        let mut q_df = Vec::with_capacity(q_results.len());
        let mut q_p = Vec::with_capacity(q_results.len());
        for q in &q_results {
            factor1.push(q.factor1.clone());
            factor2.push(q.factor2.clone());
            q_chisq.push(q.q_chisq);
            q_df.push(q.q_df as i32);
            q_p.push(q.q_p);
        }
        let q_list = list!(
            factor1 = factor1,
            factor2 = factor2,
            Q_chisq = q_chisq,
            Q_df = q_df,
            Q_p = q_p
        );
        pairs.push(("Q_Factor", q_list.into_robj()));
    }

    List::from_pairs(pairs)
}

/// Munge GWAS summary statistics files.
///
/// `column_names_json` is a small JSON object (trait column name overrides);
/// this boundary is kept as a JSON string because it's low-volume and
/// infrequent.
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
        parse_simple_json_map(column_names_json);

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

/// Parse the `{"key":"value","k2":"v2"}` mini-format that the R `munge` and
/// `sumstats` wrappers emit for column-override options. Returns `None` for
/// empty/`"{}"` inputs or on any parse error.
///
/// This replaces a `serde_json::from_str` call for the few config-string
/// boundaries we keep as JSON. It only needs to handle the exact shape the
/// R wrappers produce — not arbitrary JSON.
fn parse_simple_json_map(s: &str) -> Option<std::collections::HashMap<String, String>> {
    let trimmed = s.trim();
    if trimmed.is_empty() || trimmed == "{}" {
        return None;
    }
    let inner = trimmed.strip_prefix('{')?.strip_suffix('}')?;
    let mut map = std::collections::HashMap::new();
    // Split on top-level commas (no nesting expected from the R side).
    for pair in inner.split(',') {
        let pair = pair.trim();
        if pair.is_empty() {
            continue;
        }
        let (k, v) = pair.split_once(':')?;
        let k = unquote(k.trim())?;
        let v = unquote(v.trim())?;
        map.insert(k, v);
    }
    if map.is_empty() { None } else { Some(map) }
}

/// Parse a comma-separated list of optional f64 values in the format
/// emitted by R's `jsonlite::toJSON(as.list(N), auto_unbox=TRUE)`. Accepts
/// `null` for missing entries.
fn parse_simple_json_f64_array(s: &str) -> Vec<Option<f64>> {
    let trimmed = s.trim();
    if trimmed.is_empty() || trimmed == "[]" {
        return Vec::new();
    }
    let inner = trimmed
        .strip_prefix('[')
        .and_then(|s| s.strip_suffix(']'))
        .unwrap_or(trimmed);
    inner
        .split(',')
        .map(|v| {
            let v = v.trim();
            if v.eq_ignore_ascii_case("null") || v.is_empty() {
                None
            } else {
                v.parse::<f64>().ok()
            }
        })
        .collect()
}

/// Strip surrounding double quotes from a JSON string token.
fn unquote(s: &str) -> Option<String> {
    s.strip_prefix('"')
        .and_then(|s| s.strip_suffix('"'))
        .map(|s| s.to_string())
}

/// Fit common factor model.
#[extendr]
fn commonfactor_rust(covstruc: List, estimation: &str) -> List {
    ensure_logger();
    let ldsc_result = match conversions::list_to_ldsc_result(&covstruc) {
        Ok(r) => r,
        Err(e) => return conversions::error_list(e.to_string()),
    };

    match gsem_sem::commonfactor::run_commonfactor(
        &ldsc_result.s,
        &ldsc_result.v,
        gsem_sem::EstimationMethod::from_str_lossy(estimation),
    ) {
        Ok(result) => {
            list!(
                parameters = conversions::param_estimates_to_list(&result.parameters),
                chisq = result.fit.chisq,
                df = result.fit.df as i32,
                p_chisq = result.fit.p_chisq,
                aic = result.fit.aic,
                cfi = result.fit.cfi,
                srmr = result.fit.srmr
            )
        }
        Err(e) => conversions::error_list(e.to_string()),
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
    num_threads: Rint,
) -> List {
    ensure_logger();
    let se_logit_bool: Vec<bool> = se_logit.iter().map(|&v| v != 0).collect();
    let ols_bool: Vec<bool> = ols.iter().map(|&v| v != 0).collect();
    let linprob_bool: Vec<bool> = linprob.iter().map(|&v| v != 0).collect();

    let n_overrides: Vec<Option<f64>> = if n_overrides_json.trim().is_empty() {
        vec![None; files.len()]
    } else {
        let parsed = parse_simple_json_f64_array(n_overrides_json);
        if parsed.is_empty() {
            vec![None; files.len()]
        } else {
            parsed
        }
    };

    // Parse betas: JSON object like {"trait1":"BETA_COL"} or empty
    let beta_overrides: Vec<Option<String>> = match parse_simple_json_map(betas_json) {
        Some(map) => trait_names
            .iter()
            .map(|name| map.get(name).cloned())
            .collect(),
        None => Vec::new(),
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
        num_threads: rint_to_num_threads(num_threads),
    };
    let file_refs: Vec<&std::path::Path> = files
        .iter()
        .map(|p| std::path::Path::new(p.as_str()))
        .collect();
    let out_path = std::path::Path::new(out);
    match gsem::sumstats::merge_sumstats(
        &file_refs,
        std::path::Path::new(ref_dir),
        &trait_names,
        &config,
        out_path,
    ) {
        Ok(n) => list!(path = out.to_string(), n_snps = n as i32),
        Err(e) => conversions::error_list(e.to_string()),
    }
}

/// Run common factor GWAS.
#[extendr]
fn commonfactor_gwas_rust(
    covstruc: List,
    sumstats_path: &str,
    gc: &str,
    estimation: &str,
    snp_se: Rfloat,
    smooth_check: bool,
    twas: bool,
    identification: &str,
    num_threads: Rint,
) -> List {
    let ldsc_result = match conversions::list_to_ldsc_result(&covstruc) {
        Ok(r) => r,
        Err(e) => return conversions::error_list(e.to_string()),
    };

    let gc_mode: gsem::gwas::gc_correction::GcMode = gc
        .parse()
        .unwrap_or(gsem::gwas::gc_correction::GcMode::Standard);

    let mut i_ld = ldsc_result.i_mat.to_owned();
    clamp_i_ld_diagonal(&mut i_ld);

    let snp_se_opt = if snp_se.is_na() || snp_se.inner().is_nan() {
        None
    } else {
        Some(snp_se.inner())
    };

    let cf_config = gsem::gwas::common_factor::CommonFactorGwasConfig {
        estimation: gsem_sem::EstimationMethod::from_str_lossy(estimation),
        gc: gc_mode,
        snp_se: snp_se_opt,
        smooth_check,
        identification: gsem::gwas::common_factor::Identification::from_str_lossy(identification),
        num_threads: rint_to_num_threads(num_threads),
        ..Default::default()
    };

    if twas {
        let twas_data =
            match gsem::io::twas_reader::read_twas_sumstats(std::path::Path::new(sumstats_path)) {
                Ok(d) => d,
                Err(e) => return conversions::error_list(e.to_string()),
            };

        let beta_gene: Vec<&[f64]> = twas_data.genes.iter().map(|g| g.beta.as_slice()).collect();
        let se_gene: Vec<&[f64]> = twas_data.genes.iter().map(|g| g.se.as_slice()).collect();
        let var_gene: Vec<f64> = twas_data.genes.iter().map(|g| g.hsq).collect();

        let n_genes = var_gene.len();
        let progress = RProgress::new(n_genes);
        let cb = progress.callback();
        let results = gsem::gwas::common_factor::run_common_factor_gwas(
            &twas_data.trait_names,
            &ldsc_result.s,
            &ldsc_result.v,
            &i_ld,
            &beta_gene,
            &se_gene,
            &var_gene,
            &cf_config,
            Some(&cb),
        );

        return conversions::twas_results_to_list(&results, &twas_data.genes);
    }

    let merged = match gsem::io::sumstats_reader::read_merged_sumstats(std::path::Path::new(
        sumstats_path,
    )) {
        Ok(m) => m,
        Err(e) => return conversions::error_list(e.to_string()),
    };

    // Filter out SNPs with zero MAF (var_snp=0 causes singular matrices)
    let valid_idx: Vec<usize> = merged
        .snps
        .iter()
        .enumerate()
        .filter(|(_, s)| {
            let v = 2.0 * s.maf * (1.0 - s.maf);
            v > 1e-10
        })
        .map(|(i, _)| i)
        .collect();

    if valid_idx.is_empty() {
        return conversions::error_list("no SNPs with valid MAF (all MAF=0)");
    }

    let filtered_snps: Vec<_> = valid_idx.iter().map(|&i| &merged.snps[i]).collect();
    let beta_snp: Vec<&[f64]> = filtered_snps.iter().map(|s| s.beta.as_slice()).collect();
    let se_snp: Vec<&[f64]> = filtered_snps.iter().map(|s| s.se.as_slice()).collect();
    let var_snp: Vec<f64> = filtered_snps
        .iter()
        .map(|s| 2.0 * s.maf * (1.0 - s.maf))
        .collect();

    // Remap snp indices: create a new merged.snps with only valid entries
    let valid_merged_snps: Vec<_> = valid_idx.iter().map(|&i| merged.snps[i].clone()).collect();

    let n_snps = var_snp.len();
    let progress = RProgress::new(n_snps);
    let cb = progress.callback();
    let results = gsem::gwas::common_factor::run_common_factor_gwas(
        &merged.trait_names,
        &ldsc_result.s,
        &ldsc_result.v,
        &i_ld,
        &beta_snp,
        &se_snp,
        &var_snp,
        &cf_config,
        Some(&cb),
    );

    conversions::snp_results_to_list(&results, &valid_merged_snps)
}

/// Run user-specified GWAS.
#[extendr]
fn user_gwas_rust(
    covstruc: List,
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
    num_threads: Rint,
) -> List {
    let nt = rint_to_num_threads(num_threads);

    let ldsc_result = match conversions::list_to_ldsc_result(&covstruc) {
        Ok(r) => r,
        Err(e) => return conversions::error_list(e.to_string()),
    };

    let gc_mode: gsem::gwas::gc_correction::GcMode = gc
        .parse()
        .unwrap_or(gsem::gwas::gc_correction::GcMode::Standard);

    let mut i_ld = ldsc_result.i_mat.to_owned();
    clamp_i_ld_diagonal(&mut i_ld);

    let snp_se_opt = if snp_se.is_na() || snp_se.inner().is_nan() {
        None
    } else {
        Some(snp_se.inner())
    };

    // TWAS mode: read gene-level data instead of SNP-level
    if twas {
        let twas_data =
            match gsem::io::twas_reader::read_twas_sumstats(std::path::Path::new(sumstats_path)) {
                Ok(d) => d,
                Err(e) => return conversions::error_list(e.to_string()),
            };

        let beta_gene: Vec<&[f64]> = twas_data.genes.iter().map(|g| g.beta.as_slice()).collect();
        let se_gene: Vec<&[f64]> = twas_data.genes.iter().map(|g| g.se.as_slice()).collect();
        let var_gene: Vec<f64> = twas_data.genes.iter().map(|g| g.hsq).collect();

        // Replace "SNP" with "Gene" in model syntax
        let twas_model_str = model.replace("SNP", "Gene");
        let twas_pt = match gsem_sem::syntax::parse_model(&twas_model_str, std_lv) {
            Ok(pt) => pt,
            Err(e) => return conversions::error_list(format!("model parse error: {e}")),
        };

        let config = gsem::gwas::user_gwas::UserGwasConfig {
            model: twas_pt,
            estimation: gsem_sem::EstimationMethod::from_str_lossy(estimation),
            gc: gc_mode,
            smooth_check,
            snp_se: snp_se_opt,
            fix_measurement,
            q_snp,
            variant_label: gsem::gwas::user_gwas::VariantLabel::Gene,
            max_iter: 500,
            num_threads: nt,
        };

        let n_genes = var_gene.len();
        let progress = RProgress::new(n_genes);
        let cb = progress.callback();
        let mut results = gsem::gwas::user_gwas::run_user_gwas(
            &config,
            &ldsc_result.s,
            &ldsc_result.v,
            &i_ld,
            &beta_gene,
            &se_gene,
            &var_gene,
            Some(&cb),
        );

        if !sub.is_empty() {
            filter_results_by_sub(&mut results, sub);
        }

        return conversions::twas_results_to_list(&results, &twas_data.genes);
    }

    // Standard SNP mode
    let merged = match gsem::io::sumstats_reader::read_merged_sumstats(std::path::Path::new(
        sumstats_path,
    )) {
        Ok(m) => m,
        Err(e) => return conversions::error_list(e.to_string()),
    };

    let beta_snp: Vec<&[f64]> = merged.snps.iter().map(|s| s.beta.as_slice()).collect();
    let se_snp: Vec<&[f64]> = merged.snps.iter().map(|s| s.se.as_slice()).collect();
    let var_snp: Vec<f64> = merged
        .snps
        .iter()
        .map(|s| 2.0 * s.maf * (1.0 - s.maf))
        .collect();

    let snp_pt = match gsem_sem::syntax::parse_model(model, std_lv) {
        Ok(pt) => pt,
        Err(e) => return conversions::error_list(format!("model parse error: {e}")),
    };
    let config = gsem::gwas::user_gwas::UserGwasConfig {
        model: snp_pt,
        estimation: gsem_sem::EstimationMethod::from_str_lossy(estimation),
        gc: gc_mode,
        max_iter: 500,
        smooth_check,
        snp_se: snp_se_opt,
        variant_label: gsem::gwas::user_gwas::VariantLabel::Snp,
        fix_measurement,
        q_snp,
        num_threads: nt,
    };

    let n_snps = var_snp.len();
    let progress = RProgress::new(n_snps);
    let cb = progress.callback();
    let mut results = gsem::gwas::user_gwas::run_user_gwas(
        &config,
        &ldsc_result.s,
        &ldsc_result.v,
        &i_ld,
        &beta_snp,
        &se_snp,
        &var_snp,
        Some(&cb),
    );

    if !sub.is_empty() {
        filter_results_by_sub(&mut results, sub);
    }

    conversions::snp_results_to_list(&results, &merged.snps)
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
                patterns.contains(&key)
            });
        }
    }
}

/// Parallel analysis to determine number of factors.
#[extendr]
fn pa_ldsc_rust(
    s: RMatrix<f64>,
    v: RMatrix<f64>,
    n_sim: i32,
    percentile: Rfloat,
    diag_only: bool,
    num_threads: Rint,
) -> List {
    let s_mat = conversions::rmatrix_to_mat(&s);
    let v_mat = conversions::rmatrix_to_mat(&v);

    let p = if percentile.is_na() || percentile.inner().is_nan() {
        0.95
    } else {
        percentile.inner()
    };

    let n = n_sim as usize;
    let progress = RProgress::new(n);
    let cb = progress.callback();
    let result = gsem::stats::parallel_analysis::parallel_analysis(
        &s_mat,
        &v_mat,
        n,
        p,
        diag_only,
        rint_to_num_threads(num_threads),
        Some(&cb),
    );

    list!(
        observed = result.observed,
        simulated_95 = result.simulated_95,
        n_factors = result.n_factors as i32
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
    gsem_sem::write_model::write_model(
        &loadings, &names, cutoff, fix_resid, bifactor, mustload, common,
    )
}

/// Compute model-implied genetic correlation matrix.
#[extendr]
fn rgmodel_rust(covstruc: List, estimation: &str, model: &str, std_lv: bool, sub: &str) -> List {
    ensure_logger();
    let ldsc_result = match conversions::list_to_ldsc_result(&covstruc) {
        Ok(r) => r,
        Err(e) => return conversions::error_list(e.to_string()),
    };

    let model_opt = if model.is_empty() { None } else { Some(model) };

    let result = if !sub.is_empty() {
        // Parse sub as comma-separated phenotype names, map to indices
        let k = ldsc_result.s.nrows();
        let all_names: Vec<String> = (0..k).map(|i| format!("V{}", i + 1)).collect();
        let sub_names: Vec<&str> = sub
            .split(',')
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .collect();
        let sub_indices: Vec<usize> = sub_names
            .iter()
            .filter_map(|name| all_names.iter().position(|n| n == name))
            .collect();
        let est_method = gsem_sem::EstimationMethod::from_str_lossy(estimation);
        gsem_sem::rgmodel::run_rgmodel_sub(
            &ldsc_result.s,
            &ldsc_result.v,
            est_method,
            model_opt,
            std_lv,
            &sub_indices,
        )
    } else {
        let est_method = gsem_sem::EstimationMethod::from_str_lossy(estimation);
        gsem_sem::rgmodel::run_rgmodel_with_model(
            &ldsc_result.s,
            &ldsc_result.v,
            est_method,
            model_opt,
            std_lv,
        )
    };

    match result {
        Ok(result) => list!(
            R = conversions::mat_to_rmatrix(&result.r),
            V_R = conversions::mat_to_rmatrix(&result.v_r)
        ),
        Err(e) => conversions::error_list(e.to_string()),
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
) -> List {
    ensure_logger();
    use gsem_ldsc::hdl::{HdlConfig, HdlMethod, HdlTraitData};

    let sp = rfloat_to_options(&sample_prev);
    let pp = rfloat_to_options(&pop_prev);

    // Read trait files
    let trait_data: Vec<HdlTraitData> = match gsem::io::gwas_reader::load_trait_data(&trait_files) {
        Ok(d) => d.into_iter().map(HdlTraitData::from).collect(),
        Err(e) => return conversions::error_list(e.to_string()),
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
        Err(e) => return conversions::error_list(e.to_string()),
    };

    match gsem_ldsc::hdl::hdl(&trait_data, &sp, &pp, &ld_pieces, &config) {
        Ok(result) => {
            let ldsc_compat = result.to_ldsc_result();
            conversions::ldsc_result_to_list(&ldsc_compat)
        }
        Err(e) => conversions::error_list(e.to_string()),
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
) -> List {
    ensure_logger();
    let sp = rfloat_to_options(&sample_prev);
    let pp = rfloat_to_options(&pop_prev);

    // Read trait files
    let trait_data = match gsem::io::gwas_reader::load_trait_data(&trait_files) {
        Ok(d) => d,
        Err(e) => return conversions::error_list(e.to_string()),
    };

    let ld = std::path::Path::new(ld_dir);
    let wld = std::path::Path::new(wld_dir);
    let chromosomes: Vec<usize> = (1..=22).collect();

    // Read annotation LD scores
    let mut annot_data = match gsem_ldsc::annot_reader::read_annot_ld_scores(ld, wld, &chromosomes)
    {
        Ok(d) => d,
        Err(e) => return conversions::error_list(e.to_string()),
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
        for (j, k) in keep.iter_mut().enumerate() {
            for i in 0..n_snps {
                let v = annot_data.annot_ld[(i, j)];
                if v != 0.0 && v != 1.0 {
                    *k = false;
                    break;
                }
            }
        }
        let kept_indices: Vec<usize> = keep
            .iter()
            .enumerate()
            .filter(|&(_, &k)| k)
            .map(|(i, _)| i)
            .collect();
        if kept_indices.len() < n_annot {
            let new_annot_ld = faer::Mat::from_fn(n_snps, kept_indices.len(), |i, j| {
                annot_data.annot_ld[(i, kept_indices[j])]
            });
            let new_names: Vec<String> = kept_indices
                .iter()
                .map(|&i| annot_data.annotation_names[i].clone())
                .collect();
            let new_m: Vec<f64> = kept_indices
                .iter()
                .map(|&i| annot_data.m_annot[i])
                .collect();
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
        Ok(result) => {
            let s_annot_list =
                List::from_values(result.s_annot.iter().map(conversions::mat_to_rmatrix));
            let v_annot_list =
                List::from_values(result.v_annot.iter().map(conversions::mat_to_rmatrix));
            list!(
                annotations = result.annotations,
                m_annot = result.m_annot,
                m_total = result.m_total,
                prop = result.prop,
                I = conversions::mat_to_rmatrix(&result.i_mat),
                S_annot = s_annot_list,
                V_annot = v_annot_list
            )
        }
        Err(e) => conversions::error_list(e.to_string()),
    }
}

/// Enrichment analysis using stratified LDSC results.
#[extendr]
fn enrich_rust(
    s_baseline: RMatrix<f64>,
    s_annot: List,
    v_annot: List,
    annotation_names: Vec<String>,
    m_annot: Vec<f64>,
    m_total: f64,
) -> List {
    ensure_logger();
    let s_baseline_mat = conversions::rmatrix_to_mat(&s_baseline);

    let parse_list_of_matrices = |list: &List, label: &str| -> Result<Vec<faer::Mat<f64>>> {
        list.values()
            .enumerate()
            .map(|(i, obj)| {
                conversions::robj_to_mat(&obj)
                    .map_err(|e| Error::Other(format!("{label}[{}]: {e}", i + 1)))
            })
            .collect()
    };

    let s_annot_mats = match parse_list_of_matrices(&s_annot, "S_annot") {
        Ok(v) => v,
        Err(e) => return conversions::error_list(e.to_string()),
    };
    let v_annot_mats = match parse_list_of_matrices(&v_annot, "V_annot") {
        Ok(v) => v,
        Err(e) => return conversions::error_list(e.to_string()),
    };

    let result = gsem::stats::enrich::enrichment_test(
        &s_baseline_mat,
        &s_annot_mats,
        &v_annot_mats,
        &annotation_names,
        &m_annot,
        m_total,
    );

    list!(
        annotation = result.annotations,
        enrichment = result.enrichment,
        se = result.se,
        p = result.p
    )
}

/// Simulate GWAS summary statistics.
///
/// Returns an `RMatrix<f64>` of simulated Z-scores with shape `k × n_snps`.
#[extendr]
fn sim_ldsc_rust(
    s: RMatrix<f64>,
    n_per_trait: Vec<f64>,
    ld_scores: Vec<f64>,
    m: f64,
    int_mat: Nullable<RMatrix<f64>>,
    r_pheno: Nullable<RMatrix<f64>>,
    n_overlap: f64,
) -> Robj {
    ensure_logger();
    let s_mat = conversions::rmatrix_to_mat(&s);

    let intercepts = match int_mat {
        Nullable::NotNull(m) => Some(conversions::rmatrix_to_mat(&m)),
        Nullable::Null => None,
    };
    let r_pheno_mat = match r_pheno {
        Nullable::NotNull(m) => Some(conversions::rmatrix_to_mat(&m)),
        Nullable::Null => None,
    };

    let config = gsem::stats::simulation::SimConfig {
        intercepts,
        r_pheno: r_pheno_mat,
        n_overlap,
    };

    let result =
        gsem::stats::simulation::simulate_sumstats(&s_mat, &n_per_trait, &ld_scores, m, &config);

    // Result is Vec<Vec<f64>> with shape (k × n_snps). Pack it into an
    // `RMatrix<f64>` directly rather than a nested list.
    let k = result.len();
    let n_snps = if k == 0 { 0 } else { result[0].len() };
    conversions::mat_to_rmatrix(&faer::Mat::from_fn(k, n_snps, |i, j| result[i][j])).into_robj()
}

/// Run multi-SNP analysis.
#[extendr]
fn multi_snp_rust(
    covstruc: List,
    model: &str,
    estimation: &str,
    beta: RMatrix<f64>,
    se: RMatrix<f64>,
    var_snp: Vec<f64>,
    ld_matrix: RMatrix<f64>,
    snp_names: Vec<String>,
    snp_se: Rfloat,
) -> List {
    ensure_logger();
    let ldsc_result = match conversions::list_to_ldsc_result(&covstruc) {
        Ok(r) => r,
        Err(e) => return conversions::error_list(e.to_string()),
    };

    // beta / se are passed as n_snps × k matrices. Unpack column-major into
    // the `Vec<Vec<f64>>` layout the engine expects.
    let rmatrix_to_rows = |r: &RMatrix<f64>| -> Vec<Vec<f64>> {
        let nrows = r.nrows();
        let ncols = r.ncols();
        let data = r.data();
        (0..nrows)
            .map(|i| (0..ncols).map(|j| data[i + j * nrows]).collect())
            .collect()
    };
    let beta_rows = rmatrix_to_rows(&beta);
    let se_rows = rmatrix_to_rows(&se);
    let ld_mat = conversions::rmatrix_to_mat(&ld_matrix);

    let snp_var_se = if snp_se.is_na() || snp_se.inner().is_nan() {
        None
    } else {
        Some(snp_se.inner())
    };

    let pt = match gsem_sem::syntax::parse_model(model, false) {
        Ok(pt) => pt,
        Err(e) => return conversions::error_list(format!("model parse error: {e}")),
    };
    let config = gsem::gwas::multi_snp::MultiSnpConfig {
        model: pt,
        estimation: gsem_sem::EstimationMethod::from_str_lossy(estimation),
        max_iter: 1000,
        snp_var_se,
    };

    let result = gsem::gwas::multi_snp::run_multi_snp(
        &config,
        &ldsc_result.s,
        &ldsc_result.v,
        &beta_rows,
        &se_rows,
        &var_snp,
        &ld_mat,
        &snp_names,
    );

    list!(
        converged = result.converged,
        chisq = result.chisq,
        df = result.chisq_df as i32,
        params = conversions::param_results_to_list(&result.params)
    )
}

/// Run multi-gene analysis (reuses multi-SNP engine).
#[extendr]
fn multi_gene_rust(
    covstruc: List,
    model: &str,
    estimation: &str,
    beta: RMatrix<f64>,
    se: RMatrix<f64>,
    var_gene: Vec<f64>,
    ld_matrix: RMatrix<f64>,
    gene_names: Vec<String>,
    snp_se: Rfloat,
) -> List {
    ensure_logger();
    multi_snp_rust(
        covstruc, model, estimation, beta, se, var_gene, ld_matrix, gene_names, snp_se,
    )
}

/// Run Generalized Least Squares regression.
#[extendr]
fn summary_gls_rust(x: RMatrix<f64>, y: Vec<f64>, v: RMatrix<f64>, intercept: bool) -> List {
    ensure_logger();
    let mut x_mat = conversions::rmatrix_to_mat(&x);
    let v_mat = conversions::rmatrix_to_mat(&v);

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
        Some(result) => list!(
            beta = result.beta,
            se = result.se,
            z = result.z,
            p = result.p
        ),
        None => conversions::error_list("GLS failed (singular matrix?)"),
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
