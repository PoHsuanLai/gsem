use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;

mod conversions;

/// Python wrapper for LDSC result.
#[pyclass(name = "LdscResult")]
struct PyLdscResult {
    json: String,
    s_data: Vec<f64>,
    s_shape: (usize, usize),
    v_data: Vec<f64>,
    v_shape: (usize, usize),
    i_data: Vec<f64>,
    i_shape: (usize, usize),
    n_vec: Vec<f64>,
    m: f64,
    // Optional standardized matrices (populated when stand=True)
    s_stand_data: Option<Vec<f64>>,
    s_stand_shape: Option<(usize, usize)>,
    v_stand_data: Option<Vec<f64>>,
    v_stand_shape: Option<(usize, usize)>,
}

#[pymethods]
impl PyLdscResult {
    /// Genetic covariance matrix S as NumPy array.
    #[getter]
    fn s<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let arr =
            Array2::from_shape_vec(self.s_shape, self.s_data.clone()).expect("S shape mismatch");
        arr.into_pyarray(py)
    }

    /// Sampling covariance matrix V as NumPy array.
    #[getter]
    fn v<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let arr =
            Array2::from_shape_vec(self.v_shape, self.v_data.clone()).expect("V shape mismatch");
        arr.into_pyarray(py)
    }

    /// Intercept matrix I as NumPy array.
    #[getter]
    fn i_mat<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let arr =
            Array2::from_shape_vec(self.i_shape, self.i_data.clone()).expect("I shape mismatch");
        arr.into_pyarray(py)
    }

    /// Sample sizes per vech element.
    #[getter]
    fn n(&self) -> Vec<f64> {
        self.n_vec.clone()
    }

    /// Total number of SNPs.
    #[getter]
    fn m_total(&self) -> f64 {
        self.m
    }

    /// Standardized genetic correlation matrix (only when stand=True).
    #[getter]
    fn s_stand<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f64>>> {
        self.s_stand_data.as_ref().map(|data| {
            let shape = self.s_stand_shape.unwrap();
            Array2::from_shape_vec(shape, data.clone()).expect("S_Stand shape").into_pyarray(py)
        })
    }

    /// Standardized sampling covariance (only when stand=True).
    #[getter]
    fn v_stand<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f64>>> {
        self.v_stand_data.as_ref().map(|data| {
            let shape = self.v_stand_shape.unwrap();
            Array2::from_shape_vec(shape, data.clone()).expect("V_Stand shape").into_pyarray(py)
        })
    }

    /// Serialize to JSON string.
    fn to_json(&self) -> String {
        self.json.clone()
    }

    /// Deserialize from JSON string.
    #[staticmethod]
    fn from_json(json: &str) -> PyResult<Self> {
        let result = conversions::json_to_ldsc(json)
            .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("invalid JSON"))?;
        Ok(ldsc_result_to_py(&result, false))
    }
}

fn ldsc_result_to_py(result: &gsem_ldsc::LdscResult, stand: bool) -> PyLdscResult {
    let (s_data, s_rows, s_cols) = conversions::mat_to_flat(&result.s);
    let (v_data, v_rows, v_cols) = conversions::mat_to_flat(&result.v);
    let (i_data, i_rows, i_cols) = conversions::mat_to_flat(&result.i_mat);

    let (s_stand_data, s_stand_shape, v_stand_data, v_stand_shape) = if stand {
        let s_stand = gsem_matrix::smooth::cov_to_cor(&result.s);
        let k = result.s.nrows();
        let kstar = k * (k + 1) / 2;
        let s_vec = gsem_matrix::vech::vech(&result.s);
        let ss_vec = gsem_matrix::vech::vech(&s_stand);
        let scale: Vec<f64> = ss_vec.iter().zip(s_vec.iter()).enumerate().map(|(i, (&st, &orig))| {
            let ratio = if orig.abs() > 1e-30 { st / orig } else { 0.0 };
            result.v[(i, i)].sqrt() * ratio
        }).collect();
        let v_cor = gsem_matrix::smooth::cov_to_cor(&result.v);
        let v_stand = faer::Mat::from_fn(kstar, kstar, |i, j| scale[i] * v_cor[(i, j)] * scale[j]);
        let (ss_data, ss_r, ss_c) = conversions::mat_to_flat(&s_stand);
        let (vs_data, vs_r, vs_c) = conversions::mat_to_flat(&v_stand);
        (Some(ss_data), Some((ss_r, ss_c)), Some(vs_data), Some((vs_r, vs_c)))
    } else {
        (None, None, None, None)
    };

    PyLdscResult {
        json: result.to_json_string().unwrap_or_default(),
        s_data,
        s_shape: (s_rows, s_cols),
        v_data,
        v_shape: (v_rows, v_cols),
        i_data,
        i_shape: (i_rows, i_cols),
        n_vec: result.n_vec.clone(),
        m: result.m,
        s_stand_data,
        s_stand_shape,
        v_stand_data,
        v_stand_shape,
    }
}

/// Run multivariate LDSC.
#[pyfunction]
#[pyo3(signature = (traits, sample_prev, population_prev, ld, wld="", trait_names=None, sep_weights=false, chr=22, n_blocks=200, ldsc_log=None, stand=false, select=None, chisq_max=None))]
#[allow(unused_variables)]
fn ldsc(
    traits: Vec<String>,
    sample_prev: Vec<Option<f64>>,
    population_prev: Vec<Option<f64>>,
    ld: &str,
    wld: &str,
    trait_names: Option<Vec<String>>,
    sep_weights: bool,       // ignored
    chr: usize,
    n_blocks: usize,
    ldsc_log: Option<String>, // ignored
    stand: bool,
    select: Option<String>,
    chisq_max: Option<f64>,
) -> PyResult<PyLdscResult> {
    if sep_weights {
        log::info!("sep_weights is always enabled in gsemr — weight LD scores are read from the wld directory");
    }
    if let Some(ref path) = ldsc_log {
        log::info!("ldsc_log='{path}' — file logging is handled at the Python/R wrapper level");
    }
    let wld_dir = if wld.is_empty() { ld } else { wld };

    let mut trait_data = Vec::new();
    for path_str in &traits {
        let path = std::path::Path::new(path_str);
        let records = gsem::io::gwas_reader::read_sumstats(path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{e}")))?;
        trait_data.push(gsem_ldsc::TraitSumstats {
            snp: records.iter().map(|r| r.snp.clone()).collect(),
            z: records.iter().map(|r| r.z).collect(),
            n: records.iter().map(|r| r.n).collect(),
            a1: records.iter().map(|r| r.a1.clone()).collect(),
            a2: records.iter().map(|r| r.a2.clone()).collect(),
        });
    }

    let ld_path = std::path::Path::new(ld);
    let wld_path = std::path::Path::new(wld_dir);
    let chromosomes: Vec<usize> = match select.as_deref() {
        None | Some("") | Some("FALSE") | Some("false") => (1..=chr).collect(),
        Some("ODD") | Some("odd") => (1..=chr).filter(|c| c % 2 == 1).collect(),
        Some("EVEN") | Some("even") => (1..=chr).filter(|c| c % 2 == 0).collect(),
        Some(other) => other.split(',').filter_map(|s| s.trim().parse().ok()).collect(),
    };
    let ld_data = gsem::io::ld_reader::read_ld_scores(ld_path, wld_path, &chromosomes)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{e}")))?;

    let ld_snps: Vec<String> = ld_data.records.iter().map(|r| r.snp.clone()).collect();
    let ld_scores: Vec<f64> = ld_data.records.iter().map(|r| r.l2).collect();

    let config = gsem_ldsc::LdscConfig {
        n_blocks,
        chisq_max,
    };

    let result = gsem_ldsc::ldsc(
        &trait_data,
        &sample_prev,
        &population_prev,
        &ld_scores,
        &ld_data.w_ld,
        &ld_snps,
        ld_data.total_m,
        &config,
    )
    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;

    Ok(ldsc_result_to_py(&result, stand))
}

/// Munge GWAS summary statistics files.
#[pyfunction]
#[pyo3(signature = (files, hm3, trait_names=None, n=None, info_filter=0.9, maf_filter=0.01, log_name=None, column_names=None, parallel=false, cores=None, overwrite=true, out="."))]
#[allow(unused_variables)]
fn munge(
    files: Vec<String>,
    hm3: &str,
    trait_names: Option<Vec<String>>,
    n: Option<f64>,                          // threaded (N override)
    info_filter: f64,
    maf_filter: f64,
    log_name: Option<String>,                // ignored
    column_names: Option<std::collections::HashMap<String, String>>, // threaded
    parallel: bool,                          // ignored
    cores: Option<usize>,                    // ignored
    overwrite: bool,                         // ignored
    out: &str,
) -> PyResult<Vec<String>> {
    let hm3_path = std::path::Path::new(hm3);
    let reference = gsem::munge::read_reference(hm3_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{e}")))?;

    let config = gsem::munge::MungeConfig {
        info_filter,
        maf_filter,
        n_override: n,
        column_overrides: column_names,
    };

    let default_names: Vec<String> = files
        .iter()
        .enumerate()
        .map(|(i, _)| format!("trait{}", i + 1))
        .collect();
    let names = trait_names.as_deref().unwrap_or(&default_names);

    let mut output_paths = Vec::new();
    for (i, file) in files.iter().enumerate() {
        let name = names.get(i).map(|s| s.as_str()).unwrap_or("trait");
        let out_path = std::path::Path::new(out).join(format!("{name}.sumstats.gz"));

        gsem::munge::munge_and_write(std::path::Path::new(file), &reference, &config, &out_path)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;

        output_paths.push(out_path.to_string_lossy().to_string());
    }
    Ok(output_paths)
}

/// Fit a user-specified SEM model.
#[pyfunction]
#[pyo3(signature = (covstruc_json, model="", estimation="DWLS", cfi_calc=true, std_lv=false, imp_cov=false, fix_resid=true, toler=None, q_factor=false))]
#[allow(unused_variables)]
fn usermodel(
    covstruc_json: &str,
    model: &str,
    estimation: &str,
    cfi_calc: bool,         // ignored (CFIcalc)
    std_lv: bool,           // threaded
    imp_cov: bool,          // ignored
    fix_resid: bool,        // threaded
    toler: Option<f64>,     // ignored
    q_factor: bool,         // ignored (Q_Factor)
) -> PyResult<String> {
    let ldsc_result = conversions::json_to_ldsc(covstruc_json)
        .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("invalid covstruc JSON"))?;

    let pt = gsem_sem::syntax::parse_model(model, std_lv)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("model parse error: {e}")))?;

    let k = ldsc_result.s.nrows();
    let obs_names: Vec<String> = (0..k).map(|i| format!("V{}", i + 1)).collect();
    let mut sem_model = gsem_sem::model::Model::from_partable(&pt, &obs_names);

    let kstar = k * (k + 1) / 2;
    let v_diag: Vec<f64> = (0..kstar).map(|i| ldsc_result.v[(i, i)]).collect();

    let fit = if estimation.to_uppercase() == "ML" {
        gsem_sem::estimator::fit_ml(&mut sem_model, &ldsc_result.s, 1000, None)
    } else {
        gsem_sem::estimator::fit_dwls(&mut sem_model, &ldsc_result.s, &v_diag, 1000, None)
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

    Ok(format!(
        "{{\"converged\":{},\"objective\":{:.6},\"parameters\":[{}]}}",
        fit.converged,
        fit.objective,
        params_json.join(",")
    ))
}

/// Fit common factor model.
#[pyfunction]
#[pyo3(signature = (covstruc_json, estimation="DWLS"))]
fn commonfactor(covstruc_json: &str, estimation: &str) -> PyResult<String> {
    let ldsc_result = conversions::json_to_ldsc(covstruc_json)
        .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("invalid covstruc JSON"))?;

    let result = gsem_sem::commonfactor::run_commonfactor(&ldsc_result.s, &ldsc_result.v, estimation)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;

    let params: Vec<String> = result.parameters.iter().map(|p| {
        format!(
            "{{\"lhs\":\"{}\",\"op\":\"{}\",\"rhs\":\"{}\",\"est\":{:.6},\"se\":{:.6},\"z\":{:.4},\"p\":{:.6e}}}",
            p.lhs, p.op, p.rhs, p.est, p.se, p.z, p.p
        )
    }).collect();
    Ok(format!(
        "{{\"parameters\":[{}],\"chisq\":{:.4},\"df\":{},\"p_chisq\":{:.6e},\"aic\":{:.4},\"cfi\":{:.4},\"srmr\":{:.4}}}",
        params.join(","), result.fit.chisq, result.fit.df, result.fit.p_chisq,
        result.fit.aic, result.fit.cfi, result.fit.srmr
    ))
}

/// Merge GWAS summary statistics.
#[pyfunction]
#[pyo3(signature = (files, ref_dir, trait_names=None, se_logit=None, ols=None, linprob=None, n=None, betas=None, info_filter=0.6, maf_filter=0.01, keep_indel=false, parallel=false, cores=None, ambig=false, direct_filter=false, out="merged_sumstats.tsv"))]
#[allow(unused_variables)]
fn sumstats(
    files: Vec<String>,
    ref_dir: &str,
    trait_names: Option<Vec<String>>,
    se_logit: Option<Vec<bool>>,   // threaded
    ols: Option<Vec<bool>>,        // threaded (OLS)
    linprob: Option<Vec<bool>>,    // threaded
    n: Option<Vec<f64>>,           // threaded (N overrides)
    betas: Option<Vec<f64>>,       // ignored
    info_filter: f64,
    maf_filter: f64,
    keep_indel: bool,
    parallel: bool,                // ignored
    cores: Option<usize>,          // ignored
    ambig: bool,                   // ignored
    direct_filter: bool,           // ignored
    out: &str,
) -> PyResult<String> {
    let k = files.len();
    let default_names: Vec<String> = (0..k).map(|i| format!("trait{}", i + 1)).collect();
    let names = trait_names.unwrap_or(default_names);

    let config = gsem::sumstats::SumstatsConfig {
        info_filter,
        maf_filter,
        keep_indel,
        keep_ambig: ambig,
        se_logit: se_logit.unwrap_or_else(|| vec![false; k]),
        ols: ols.unwrap_or_else(|| vec![false; k]),
        linprob: linprob.unwrap_or_else(|| vec![false; k]),
        n_overrides: n.map(|v| v.into_iter().map(Some).collect()).unwrap_or_else(|| vec![None; k]),
        beta_overrides: Vec::new(),
        direct_filter: false,
    };
    let file_refs: Vec<&std::path::Path> = files.iter().map(|p| std::path::Path::new(p.as_str())).collect();
    gsem::sumstats::merge_sumstats(&file_refs, std::path::Path::new(ref_dir), &names, &config, std::path::Path::new(out))
        .map(|n_snps| format!("{{\"path\":\"{out}\",\"n_snps\":{n_snps}}}"))
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))
}

/// Run common factor GWAS.
#[pyfunction]
#[pyo3(signature = (covstruc_json, sumstats_path, estimation="DWLS", cores=None, toler=false, snpse=false, parallel=true, gc="standard", mpi=false, twas=false, smooth_check=false))]
#[allow(unused_variables)]
fn commonfactor_gwas(
    covstruc_json: &str,
    sumstats_path: &str,
    estimation: &str,          // threaded
    cores: Option<usize>,      // ignored
    toler: bool,               // ignored
    snpse: bool,               // threaded (as snp_se)
    parallel: bool,            // ignored
    gc: &str,
    mpi: bool,                 // ignored
    twas: bool,                // ignored
    smooth_check: bool,        // threaded
) -> PyResult<String> {
    let ldsc_result = conversions::json_to_ldsc(covstruc_json)
        .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("invalid covstruc JSON"))?;
    let merged = gsem::io::sumstats_reader::read_merged_sumstats(std::path::Path::new(sumstats_path))
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{e}")))?;

    let k = ldsc_result.s.nrows();
    let gc_mode: gsem::gwas::gc_correction::GcMode = gc.parse().unwrap_or(gsem::gwas::gc_correction::GcMode::Standard);
    let beta_snp: Vec<Vec<f64>> = merged.snps.iter().map(|s| s.beta.clone()).collect();
    let se_snp: Vec<Vec<f64>> = merged.snps.iter().map(|s| s.se.clone()).collect();
    let var_snp: Vec<f64> = merged.snps.iter().map(|s| 2.0 * s.maf * (1.0 - s.maf)).collect();
    let mut i_ld = ldsc_result.i_mat.to_owned();
    for i in 0..k { if i_ld[(i, i)] < 1.0 { i_ld[(i, i)] = 1.0; } }

    let snp_se_val = if snpse { Some(0.0005) } else { None };

    // Common factor GWAS uses the internal user_gwas path with auto-generated model;
    // we build a UserGwasConfig to thread estimation, snp_se, smooth_check.
    let loading = std::iter::once(format!("NA*{}", merged.trait_names[0]))
        .chain(merged.trait_names[1..].iter().cloned())
        .collect::<Vec<_>>()
        .join(" + ");
    let model = format!("F1 =~ {loading}\nF1 ~ SNP\nF1 ~~ 1*F1\nSNP ~~ SNP");

    let config = gsem::gwas::user_gwas::UserGwasConfig {
        model,
        estimation: estimation.to_string(),
        gc: gc_mode,
        max_iter: 500,
        std_lv: false,
        smooth_check,
        snp_se: snp_se_val,
        snp_label: "SNP".to_string(),
        q_snp: false,
        fix_measurement: false,
    };

    let results = gsem::gwas::user_gwas::run_user_gwas(
        &config, &ldsc_result.s, &ldsc_result.v, &i_ld,
        &beta_snp, &se_snp, &var_snp,
    );
    Ok(snp_results_to_json(&results, &merged.snps))
}

/// Run user-specified GWAS.
#[pyfunction]
#[pyo3(signature = (covstruc_json, sumstats_path, model="", estimation="DWLS", printwarn=true, sub=None, cores=None, toler=false, snpse=false, parallel=true, gc="standard", mpi=false, smooth_check=false, twas=false, std_lv=false, fix_measurement=true, q_snp=false))]
#[allow(unused_variables)]
fn user_gwas(
    covstruc_json: &str,
    sumstats_path: &str,
    model: &str,
    estimation: &str,
    printwarn: bool,
    sub: Option<Vec<String>>,
    cores: Option<usize>,
    toler: bool,
    snpse: bool,
    parallel: bool,
    gc: &str,
    mpi: bool,
    smooth_check: bool,
    twas: bool,
    std_lv: bool,
    fix_measurement: bool,
    q_snp: bool,
) -> PyResult<String> {
    let ldsc_result = conversions::json_to_ldsc(covstruc_json)
        .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("invalid covstruc JSON"))?;
    let merged = gsem::io::sumstats_reader::read_merged_sumstats(std::path::Path::new(sumstats_path))
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{e}")))?;

    let k = ldsc_result.s.nrows();
    let gc_mode: gsem::gwas::gc_correction::GcMode = gc.parse().unwrap_or(gsem::gwas::gc_correction::GcMode::Standard);
    let beta_snp: Vec<Vec<f64>> = merged.snps.iter().map(|s| s.beta.clone()).collect();
    let se_snp: Vec<Vec<f64>> = merged.snps.iter().map(|s| s.se.clone()).collect();
    let var_snp: Vec<f64> = merged.snps.iter().map(|s| 2.0 * s.maf * (1.0 - s.maf)).collect();
    let mut i_ld = ldsc_result.i_mat.to_owned();
    for i in 0..k { if i_ld[(i, i)] < 1.0 { i_ld[(i, i)] = 1.0; } }

    let snp_se_val = if snpse { Some(0.0005) } else { None };
    let snp_label = if twas { "Gene".to_string() } else { "SNP".to_string() };

    let config = gsem::gwas::user_gwas::UserGwasConfig {
        model: model.to_string(),
        estimation: estimation.to_string(),
        gc: gc_mode,
        max_iter: 500,
        std_lv,
        smooth_check,
        snp_se: snp_se_val,
        snp_label,
        q_snp,
        fix_measurement,
    };
    if !printwarn {
        log::set_max_level(log::LevelFilter::Error);
    }
    if toler {
        log::info!("toler — convergence tolerance is controlled by the L-BFGS optimizer internally");
    }
    if !parallel {
        rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build_global()
            .ok();
    }
    if mpi {
        log::warn!("MPI is not supported in gsemr — use the cores parameter for thread control");
    }
    if let Some(n) = cores {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .ok(); // ignore if already initialized
    }

    let mut results = gsem::gwas::user_gwas::run_user_gwas(
        &config, &ldsc_result.s, &ldsc_result.v, &i_ld,
        &beta_snp, &se_snp, &var_snp,
    );

    if let Some(ref patterns) = sub {
        let pats: Vec<String> = patterns.iter()
            .map(|s| s.trim().replace(' ', ""))
            .filter(|s| !s.is_empty())
            .collect();
        if !pats.is_empty() {
            for snp_result in &mut results {
                snp_result.params.retain(|p| {
                    let key = format!("{}{}{}", p.lhs, p.op, p.rhs).replace(' ', "");
                    pats.iter().any(|pat| key == *pat)
                });
            }
        }
    }

    Ok(snp_results_to_json(&results, &merged.snps))
}

/// Parallel analysis to determine number of factors.
/// Accepts S and V as JSON strings (2D arrays).
#[pyfunction]
#[pyo3(signature = (s_json, v_json, r=500, p=None, diag=false))]
fn parallel_analysis(
    s_json: &str,
    v_json: &str,
    r: usize,
    p: Option<f64>,
    diag: bool,
) -> PyResult<String> {
    let s_mat = json_to_mat(s_json)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("invalid s_json: {e}")))?;
    let v_mat = json_to_mat(v_json)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("invalid v_json: {e}")))?;

    let percentile = p.unwrap_or(0.95);
    let result = gsem::stats::parallel_analysis::parallel_analysis(&s_mat, &v_mat, r, percentile, diag);
    let obs: Vec<String> = result.observed.iter().map(|v| format!("{v:.6}")).collect();
    let sim: Vec<String> = result.simulated_95.iter().map(|v| format!("{v:.6}")).collect();
    Ok(format!(
        "{{\"observed\":[{}],\"simulated_95\":[{}],\"n_factors\":{}}}",
        obs.join(","), sim.join(","), result.n_factors
    ))
}

/// Auto-generate model syntax from factor loadings.
#[pyfunction]
#[pyo3(signature = (loadings, names, cutoff=0.3, fix_resid=true, bifactor=false, mustload=false, common=false))]
#[allow(unused_variables)]
fn write_model(
    loadings: Vec<Vec<f64>>,
    names: Vec<String>,
    cutoff: f64,
    fix_resid: bool,
    bifactor: bool,
    mustload: bool,    // ignored
    common: bool,      // ignored
) -> String {
    let n_rows = loadings.len();
    let n_cols = if n_rows > 0 { loadings[0].len() } else { 0 };
    let mat = faer::Mat::from_fn(n_rows, n_cols, |i, j| loadings[i][j]);
    gsem_sem::write_model::write_model(&mat, &names, cutoff, fix_resid, bifactor, mustload, common)
}

/// Compute model-implied genetic correlation matrix.
/// `estimation=True` in R means DWLS; here `estimation=true` maps to "DWLS".
#[pyfunction]
#[pyo3(signature = (covstruc_json, model, std_lv=true, estimation=true, sub=None))]
#[allow(unused_variables)]
fn rgmodel(
    covstruc_json: &str,
    model: &str,
    std_lv: bool,
    estimation: bool,
    sub: Option<Vec<String>>,  // not yet used
) -> PyResult<String> {
    let est_str = if estimation { "DWLS" } else { "ML" };
    let model_opt = if model.is_empty() { None } else { Some(model) };
    let ldsc_result = conversions::json_to_ldsc(covstruc_json)
        .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("invalid covstruc JSON"))?;
    let result = if let Some(ref sub_names) = sub {
        let k = ldsc_result.s.nrows();
        let all_names: Vec<String> = (0..k).map(|i| format!("V{}", i + 1)).collect();
        let sub_indices: Vec<usize> = sub_names.iter()
            .filter_map(|name| all_names.iter().position(|n| n == name))
            .collect();
        gsem_sem::rgmodel::run_rgmodel_sub(
            &ldsc_result.s, &ldsc_result.v, est_str, model_opt, std_lv, &sub_indices,
        )
    } else {
        gsem_sem::rgmodel::run_rgmodel_with_model(
            &ldsc_result.s, &ldsc_result.v, est_str, model_opt, std_lv,
        )
    }
    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;
    let k = result.r.nrows();
    let kstar = k * (k + 1) / 2;
    let r_rows: Vec<String> = (0..k).map(|i| {
        let row: Vec<String> = (0..k).map(|j| format!("{:.6}", result.r[(i, j)])).collect();
        format!("[{}]", row.join(","))
    }).collect();
    let v_rows: Vec<String> = (0..kstar).map(|i| {
        let row: Vec<String> = (0..kstar).map(|j| format!("{:.6}", result.v_r[(i, j)])).collect();
        format!("[{}]", row.join(","))
    }).collect();
    Ok(format!("{{\"R\":[{}],\"V_R\":[{}]}}", r_rows.join(","), v_rows.join(",")))
}

/// Run HDL (High-Definition Likelihood) estimation.
#[pyfunction]
#[pyo3(signature = (traits, sample_prev=None, population_prev=None, trait_names=None, ld_path="", n_ref=335265.0, method="piecewise"))]
#[allow(unused_variables)]
fn hdl(
    traits: Vec<String>,
    sample_prev: Option<Vec<Option<f64>>>,
    population_prev: Option<Vec<Option<f64>>>,
    trait_names: Option<Vec<String>>,
    ld_path: &str,
    n_ref: f64,
    method: &str,
) -> PyResult<String> {
    use gsem_ldsc::hdl::{HdlConfig, HdlMethod, HdlTraitData, LdPiece};
    if trait_names.is_some() {
        log::info!("trait_names are used for labeling output in the Python wrapper");
    }

    let sp: Vec<Option<f64>> = sample_prev.unwrap_or_default();
    let pp: Vec<Option<f64>> = population_prev.unwrap_or_default();

    // Read trait files
    let mut trait_data = Vec::new();
    for path_str in &traits {
        let path = std::path::Path::new(path_str);
        let records = gsem::io::gwas_reader::read_sumstats(path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{e}")))?;
        trait_data.push(HdlTraitData {
            snp: records.iter().map(|r| r.snp.clone()).collect(),
            z: records.iter().map(|r| r.z).collect(),
            n: records.iter().map(|r| r.n).collect(),
            a1: records.iter().map(|r| r.a1.clone()).collect(),
            a2: records.iter().map(|r| r.a2.clone()).collect(),
        });
    }

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
    let pieces_file = ld_dir.join("pieces.tsv");
    let pieces_alt = ld_dir.join("pieces.txt");
    let pieces_path = if pieces_file.exists() {
        pieces_file
    } else if pieces_alt.exists() {
        pieces_alt
    } else {
        return Err(pyo3::exceptions::PyIOError::new_err(format!(
            "HDL LD reference directory missing pieces.tsv at {}",
            ld_dir.display()
        )));
    };

    let pieces_content = std::fs::read_to_string(&pieces_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("failed to read pieces file: {e}")))?;

    let mut ld_pieces = Vec::new();
    for line in pieces_content.lines() {
        let line = line.trim();
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("piece")
            || line.starts_with("chr")
        {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }

        // Try to find piece SNP file with various naming patterns
        let chr_val = fields[0];
        let piece_val = fields[1];
        let snp_file_tsv = ld_dir.join(format!("chr{chr_val}.{piece_val}.snps.tsv"));
        let snp_file_txt = ld_dir.join(format!("piece.{piece_val}.snps.txt"));
        let snp_file = if snp_file_tsv.exists() {
            snp_file_tsv
        } else if snp_file_txt.exists() {
            snp_file_txt
        } else {
            continue;
        };

        let snp_content = match std::fs::read_to_string(&snp_file) {
            Ok(c) => c,
            Err(_) => continue,
        };

        let mut snps = Vec::new();
        let mut a1 = Vec::new();
        let mut a2 = Vec::new();
        let mut ld_scores = Vec::new();

        for sline in snp_content.lines() {
            let sline = sline.trim();
            if sline.is_empty() || sline.starts_with('#') || sline.starts_with("SNP") {
                continue;
            }
            let sf: Vec<&str> = sline.split('\t').collect();
            if sf.len() < 4 {
                continue;
            }
            snps.push(sf[0].to_string());
            a1.push(sf[1].to_string());
            a2.push(sf[2].to_string());
            if let Ok(ld) = sf[3].parse::<f64>() {
                ld_scores.push(ld);
            } else {
                ld_scores.push(0.0);
            }
        }

        let m = snps.len();
        if m > 0 {
            ld_pieces.push(LdPiece {
                snps,
                a1,
                a2,
                ld_scores,
                m,
            });
        }
    }

    if ld_pieces.is_empty() {
        return Err(pyo3::exceptions::PyIOError::new_err(format!(
            "no valid LD pieces loaded from {}",
            ld_dir.display()
        )));
    }

    let result = gsem_ldsc::hdl::hdl(&trait_data, &sp, &pp, &ld_pieces, &config)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;
    let ldsc_compat = result.to_ldsc_result();
    let json = ldsc_compat.to_json_string()
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;
    Ok(json)
}

/// Run stratified LDSC (s-LDSC).
#[pyfunction]
#[pyo3(signature = (traits, sample_prev=None, population_prev=None, ld="", wld="", frq="", trait_names=None, n_blocks=200, ldsc_log=None, exclude_cont=true))]
#[allow(unused_variables)]
fn s_ldsc(
    traits: Vec<String>,
    sample_prev: Option<Vec<Option<f64>>>,
    population_prev: Option<Vec<Option<f64>>>,
    ld: &str,
    wld: &str,
    frq: &str,
    trait_names: Option<Vec<String>>,
    n_blocks: usize,
    ldsc_log: Option<String>,
    exclude_cont: bool,
) -> PyResult<String> {
    if trait_names.is_some() {
        log::info!("trait_names are used for labeling output in the Python wrapper");
    }
    if let Some(ref path) = ldsc_log {
        log::info!("ldsc_log='{path}' — file logging is handled at the Python/R wrapper level");
    }

    let sp: Vec<Option<f64>> = sample_prev.unwrap_or_default();
    let pp: Vec<Option<f64>> = population_prev.unwrap_or_default();

    // Read trait files
    let mut trait_data = Vec::new();
    for path_str in &traits {
        let path = std::path::Path::new(path_str);
        let records = gsem::io::gwas_reader::read_sumstats(path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{e}")))?;
        trait_data.push(gsem_ldsc::TraitSumstats {
            snp: records.iter().map(|r| r.snp.clone()).collect(),
            z: records.iter().map(|r| r.z).collect(),
            n: records.iter().map(|r| r.n).collect(),
            a1: records.iter().map(|r| r.a1.clone()).collect(),
            a2: records.iter().map(|r| r.a2.clone()).collect(),
        });
    }

    let ld_path = std::path::Path::new(ld);
    let wld_path = std::path::Path::new(wld);
    let chromosomes: Vec<usize> = (1..=22).collect();

    // Read annotation LD scores
    let mut annot_data = gsem_ldsc::annot_reader::read_annot_ld_scores(ld_path, wld_path, &chromosomes)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{e}")))?;

    // Filter by frq files (MAF 0.05-0.95) if frq dir is provided
    if !frq.is_empty() {
        let frq_path = std::path::Path::new(frq);
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
        n_blocks,
        rm_flank: false,
        flank_kb: 500,
    };

    let result = gsem_ldsc::stratified::s_ldsc(
        &trait_data,
        &sp,
        &pp,
        &annot_data.annot_ld,
        &annot_data.w_ld,
        &annot_data.snps,
        &annot_data.annotation_names,
        &annot_data.m_annot,
        &config,
        None,
        None,
    )
    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;

    result.to_json_string()
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))
}

/// Helper: parse a JSON 2D array string into a faer Mat.
fn json_to_mat(json: &str) -> Result<faer::Mat<f64>, String> {
    let arr: Vec<Vec<f64>> = serde_json::from_str(json)
        .map_err(|e| format!("JSON parse error: {e}"))?;
    let nrows = arr.len();
    if nrows == 0 {
        return Err("empty matrix".to_string());
    }
    let ncols = arr[0].len();
    if arr.iter().any(|row| row.len() != ncols) {
        return Err("jagged matrix".to_string());
    }
    Ok(faer::Mat::from_fn(nrows, ncols, |i, j| arr[i][j]))
}

/// Helper: serialize GWAS results to JSON.
fn snp_results_to_json(
    results: &[gsem::gwas::user_gwas::SnpResult],
    snps: &[gsem::io::sumstats_reader::MergedSnp],
) -> String {
    let entries: Vec<String> = results.iter().map(|r| {
        let snp_name = &snps[r.snp_idx].snp;
        let params: Vec<String> = r.params.iter().map(|p| {
            format!(
                "{{\"lhs\":\"{}\",\"op\":\"{}\",\"rhs\":\"{}\",\"est\":{:.6},\"se\":{:.6},\"z\":{:.4},\"p\":{:.6e}}}",
                p.lhs, p.op, p.rhs, p.est, p.se, p.z_stat, p.p_value
            )
        }).collect();
        format!(
            "{{\"SNP\":\"{}\",\"chisq\":{:.4},\"df\":{},\"converged\":{},\"params\":[{}]}}",
            snp_name, r.chisq, r.chisq_df, r.converged, params.join(",")
        )
    }).collect();
    format!("[{}]", entries.join(","))
}

/// Enrichment analysis using stratified LDSC results.
#[pyfunction]
#[pyo3(signature = (s_baseline, s_annot, v_annot, annotation_names, m_annot, m_total))]
fn enrich(
    s_baseline: Vec<Vec<f64>>,
    s_annot: Vec<Vec<Vec<f64>>>,
    v_annot: Vec<Vec<Vec<f64>>>,
    annotation_names: Vec<String>,
    m_annot: Vec<f64>,
    m_total: f64,
) -> PyResult<String> {
    let nr = s_baseline.len();
    let nc = if nr > 0 { s_baseline[0].len() } else { 0 };
    let s_base_mat = faer::Mat::from_fn(nr, nc, |i, j| s_baseline[i][j]);

    let s_annot_mats: Vec<faer::Mat<f64>> = s_annot
        .iter()
        .map(|rows| {
            let r = rows.len();
            let c = if r > 0 { rows[0].len() } else { 0 };
            faer::Mat::from_fn(r, c, |i, j| rows[i][j])
        })
        .collect();

    let v_annot_mats: Vec<faer::Mat<f64>> = v_annot
        .iter()
        .map(|rows| {
            let r = rows.len();
            let c = if r > 0 { rows[0].len() } else { 0 };
            faer::Mat::from_fn(r, c, |i, j| rows[i][j])
        })
        .collect();

    let result = gsem::stats::enrich::enrichment_test(
        &s_base_mat,
        &s_annot_mats,
        &v_annot_mats,
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
    Ok(format!("[{}]", entries.join(",")))
}

/// Simulate GWAS summary statistics.
#[pyfunction]
#[pyo3(signature = (s_matrix, n_per_trait, ld_scores, m, intercepts=None, r_pheno=None, n_overlap=0.0))]
fn sim_ldsc(
    s_matrix: Vec<Vec<f64>>,
    n_per_trait: Vec<f64>,
    ld_scores: Vec<f64>,
    m: f64,
    intercepts: Option<Vec<Vec<f64>>>,
    r_pheno: Option<Vec<Vec<f64>>>,
    n_overlap: f64,
) -> PyResult<Vec<Vec<f64>>> {
    let nr = s_matrix.len();
    let nc = if nr > 0 { s_matrix[0].len() } else { 0 };
    let s_mat = faer::Mat::from_fn(nr, nc, |i, j| s_matrix[i][j]);

    let int_mat = intercepts.map(|rows| {
        let nr = rows.len();
        let nc = if nr > 0 { rows[0].len() } else { 0 };
        faer::Mat::from_fn(nr, nc, |i, j| rows[i][j])
    });

    let r_pheno_mat = r_pheno.map(|rows| {
        let nr = rows.len();
        let nc = if nr > 0 { rows[0].len() } else { 0 };
        faer::Mat::from_fn(nr, nc, |i, j| rows[i][j])
    });

    let config = gsem::stats::simulation::SimConfig {
        intercepts: int_mat,
        r_pheno: r_pheno_mat,
        n_overlap,
    };

    Ok(gsem::stats::simulation::simulate_sumstats(
        &s_mat,
        &n_per_trait,
        &ld_scores,
        m,
        &config,
    ))
}

/// Run multi-SNP analysis.
#[pyfunction]
#[pyo3(signature = (covstruc_json, model, beta, se, var_snp, ld_matrix, snp_names, estimation="DWLS", max_iter=1000))]
fn multi_snp(
    covstruc_json: &str,
    model: &str,
    beta: Vec<Vec<f64>>,
    se: Vec<Vec<f64>>,
    var_snp: Vec<f64>,
    ld_matrix: Vec<Vec<f64>>,
    snp_names: Vec<String>,
    estimation: &str,
    max_iter: usize,
) -> PyResult<String> {
    let ldsc_result = conversions::json_to_ldsc(covstruc_json)
        .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("invalid covstruc JSON"))?;

    let nr = ld_matrix.len();
    let nc = if nr > 0 { ld_matrix[0].len() } else { 0 };
    let ld_mat = faer::Mat::from_fn(nr, nc, |i, j| ld_matrix[i][j]);

    let config = gsem::gwas::multi_snp::MultiSnpConfig {
        model: model.to_string(),
        estimation: estimation.to_string(),
        max_iter,
        snp_var_se: None,
    };

    let result = gsem::gwas::multi_snp::run_multi_snp(
        &config,
        &ldsc_result.s,
        &ldsc_result.v,
        &beta,
        &se,
        &var_snp,
        &ld_mat,
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
    Ok(format!(
        "{{\"converged\":{},\"chisq\":{:.4},\"df\":{},\"params\":[{}]}}",
        result.converged, result.chisq, result.chisq_df, params.join(",")
    ))
}

/// Run multi-gene analysis (reuses multi-SNP engine).
#[pyfunction]
#[pyo3(signature = (covstruc_json, model, beta, se, var_gene, ld_matrix, gene_names, estimation="DWLS", max_iter=1000))]
fn multi_gene(
    covstruc_json: &str,
    model: &str,
    beta: Vec<Vec<f64>>,
    se: Vec<Vec<f64>>,
    var_gene: Vec<f64>,
    ld_matrix: Vec<Vec<f64>>,
    gene_names: Vec<String>,
    estimation: &str,
    max_iter: usize,
) -> PyResult<String> {
    multi_snp(
        covstruc_json, model, beta, se, var_gene, ld_matrix, gene_names, estimation, max_iter,
    )
}

/// Run Generalized Least Squares regression.
#[pyfunction]
#[pyo3(signature = (x, y, v, intercept=true))]
fn summary_gls(
    x: Vec<Vec<f64>>,
    y: Vec<f64>,
    v: Vec<Vec<f64>>,
    intercept: bool,
) -> PyResult<String> {
    let nr = x.len();
    let nc = if nr > 0 { x[0].len() } else { 0 };
    let mut x_mat = faer::Mat::from_fn(nr, nc, |i, j| x[i][j]);

    let v_nr = v.len();
    let v_nc = if v_nr > 0 { v[0].len() } else { 0 };
    let v_mat = faer::Mat::from_fn(v_nr, v_nc, |i, j| v[i][j]);

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
            Ok(format!("[{}]", entries.join(",")))
        }
        None => Err(pyo3::exceptions::PyRuntimeError::new_err(
            "GLS failed (singular matrix?)",
        )),
    }
}

/// Python module definition.
#[pymodule]
fn genomicsem(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyLdscResult>()?;
    m.add_function(wrap_pyfunction!(ldsc, m)?)?;
    m.add_function(wrap_pyfunction!(usermodel, m)?)?;
    m.add_function(wrap_pyfunction!(munge, m)?)?;
    m.add_function(wrap_pyfunction!(commonfactor, m)?)?;
    m.add_function(wrap_pyfunction!(sumstats, m)?)?;
    m.add_function(wrap_pyfunction!(commonfactor_gwas, m)?)?;
    m.add_function(wrap_pyfunction!(user_gwas, m)?)?;
    m.add_function(wrap_pyfunction!(parallel_analysis, m)?)?;
    m.add_function(wrap_pyfunction!(write_model, m)?)?;
    m.add_function(wrap_pyfunction!(rgmodel, m)?)?;
    m.add_function(wrap_pyfunction!(hdl, m)?)?;
    m.add_function(wrap_pyfunction!(s_ldsc, m)?)?;
    m.add_function(wrap_pyfunction!(enrich, m)?)?;
    m.add_function(wrap_pyfunction!(sim_ldsc, m)?)?;
    m.add_function(wrap_pyfunction!(multi_snp, m)?)?;
    m.add_function(wrap_pyfunction!(multi_gene, m)?)?;
    m.add_function(wrap_pyfunction!(summary_gls, m)?)?;
    Ok(())
}
