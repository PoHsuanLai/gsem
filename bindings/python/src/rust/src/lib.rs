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

    /// Serialize to JSON string.
    fn to_json(&self) -> String {
        self.json.clone()
    }

    /// Deserialize from JSON string.
    #[staticmethod]
    fn from_json(json: &str) -> PyResult<Self> {
        let result = conversions::json_to_ldsc(json)
            .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("invalid JSON"))?;
        Ok(ldsc_result_to_py(&result))
    }
}

fn ldsc_result_to_py(result: &gsem_ldsc::LdscResult) -> PyLdscResult {
    let (s_data, s_rows, s_cols) = conversions::mat_to_flat(&result.s);
    let (v_data, v_rows, v_cols) = conversions::mat_to_flat(&result.v);
    let (i_data, i_rows, i_cols) = conversions::mat_to_flat(&result.i_mat);

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
    }
}

/// Run multivariate LDSC.
#[pyfunction]
#[pyo3(signature = (traits, sample_prev, pop_prev, ld, wld="", n_blocks=200))]
fn ldsc(
    traits: Vec<String>,
    sample_prev: Vec<Option<f64>>,
    pop_prev: Vec<Option<f64>>,
    ld: &str,
    wld: &str,
    n_blocks: usize,
) -> PyResult<PyLdscResult> {
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
    let chromosomes: Vec<usize> = (1..=22).collect();
    let ld_data = gsem::io::ld_reader::read_ld_scores(ld_path, wld_path, &chromosomes)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{e}")))?;

    let ld_snps: Vec<String> = ld_data.records.iter().map(|r| r.snp.clone()).collect();
    let ld_scores: Vec<f64> = ld_data.records.iter().map(|r| r.l2).collect();

    let config = gsem_ldsc::LdscConfig {
        n_blocks,
        chisq_max: None,
    };

    let result = gsem_ldsc::ldsc(
        &trait_data,
        &sample_prev,
        &pop_prev,
        &ld_scores,
        &ld_data.w_ld,
        &ld_snps,
        ld_data.total_m,
        &config,
    )
    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;

    Ok(ldsc_result_to_py(&result))
}

/// Fit a user-specified SEM model.
#[pyfunction]
#[pyo3(signature = (covstruc_json, model, estimation="DWLS"))]
fn usermodel(covstruc_json: &str, model: &str, estimation: &str) -> PyResult<String> {
    let ldsc_result = conversions::json_to_ldsc(covstruc_json)
        .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("invalid covstruc JSON"))?;

    let pt = gsem_sem::syntax::parse_model(model, false)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("model parse error: {e}")))?;

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

    Ok(format!(
        "{{\"converged\":{},\"objective\":{:.6},\"parameters\":[{}]}}",
        fit.converged,
        fit.objective,
        params_json.join(",")
    ))
}

/// Munge GWAS summary statistics files.
#[pyfunction]
#[pyo3(signature = (files, hm3, trait_names, info_filter=0.9, maf_filter=0.01, out_dir="."))]
fn munge(
    files: Vec<String>,
    hm3: &str,
    trait_names: Vec<String>,
    info_filter: f64,
    maf_filter: f64,
    out_dir: &str,
) -> PyResult<Vec<String>> {
    let hm3_path = std::path::Path::new(hm3);
    let reference = gsem::munge::read_reference(hm3_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{e}")))?;

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

        gsem::munge::munge_and_write(std::path::Path::new(file), &reference, &config, &out_path)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;

        output_paths.push(out_path.to_string_lossy().to_string());
    }
    Ok(output_paths)
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
#[pyo3(signature = (files, ref_dir, trait_names, info_filter=0.6, maf_filter=0.01, keep_indel=false, out="merged_sumstats.tsv"))]
fn sumstats(
    files: Vec<String>,
    ref_dir: &str,
    trait_names: Vec<String>,
    info_filter: f64,
    maf_filter: f64,
    keep_indel: bool,
    out: &str,
) -> PyResult<String> {
    let config = gsem::sumstats::SumstatsConfig {
        info_filter,
        maf_filter,
        keep_indel,
        ..Default::default()
    };
    let file_refs: Vec<&std::path::Path> = files.iter().map(|p| std::path::Path::new(p.as_str())).collect();
    gsem::sumstats::merge_sumstats(&file_refs, std::path::Path::new(ref_dir), &trait_names, &config, std::path::Path::new(out))
        .map(|n| format!("{{\"path\":\"{out}\",\"n_snps\":{n}}}"))
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))
}

/// Run common factor GWAS.
#[pyfunction]
#[pyo3(signature = (covstruc_json, sumstats_path, gc="standard"))]
fn commonfactor_gwas(covstruc_json: &str, sumstats_path: &str, gc: &str) -> PyResult<String> {
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

    let results = gsem::gwas::common_factor::run_common_factor_gwas(
        &merged.trait_names, &ldsc_result.s, &ldsc_result.v, &i_ld,
        &beta_snp, &se_snp, &var_snp, gc_mode,
    );
    Ok(snp_results_to_json(&results, &merged.snps))
}

/// Run user-specified GWAS.
#[pyfunction]
#[pyo3(signature = (covstruc_json, sumstats_path, model, estimation="DWLS", gc="standard"))]
fn user_gwas(covstruc_json: &str, sumstats_path: &str, model: &str, estimation: &str, gc: &str) -> PyResult<String> {
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
    Ok(snp_results_to_json(&results, &merged.snps))
}

/// Parallel analysis to determine number of factors.
#[pyfunction]
#[pyo3(signature = (covstruc_json, n_sim=500))]
fn parallel_analysis(covstruc_json: &str, n_sim: usize) -> PyResult<String> {
    let ldsc_result = conversions::json_to_ldsc(covstruc_json)
        .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("invalid covstruc JSON"))?;
    let result = gsem::stats::parallel_analysis::parallel_analysis(
        &ldsc_result.s, &ldsc_result.v, n_sim,
    );
    let obs: Vec<String> = result.observed.iter().map(|v| format!("{v:.6}")).collect();
    let sim: Vec<String> = result.simulated_95.iter().map(|v| format!("{v:.6}")).collect();
    Ok(format!(
        "{{\"observed\":[{}],\"simulated_95\":[{}],\"n_factors\":{}}}",
        obs.join(","), sim.join(","), result.n_factors
    ))
}

/// Auto-generate model syntax from factor loadings.
#[pyfunction]
#[pyo3(signature = (loadings, names, cutoff=0.3, fix_resid=false, bifactor=false))]
fn write_model(
    loadings: Vec<Vec<f64>>,
    names: Vec<String>,
    cutoff: f64,
    fix_resid: bool,
    bifactor: bool,
) -> String {
    let n_rows = loadings.len();
    let n_cols = if n_rows > 0 { loadings[0].len() } else { 0 };
    let mat = faer::Mat::from_fn(n_rows, n_cols, |i, j| loadings[i][j]);
    gsem_sem::write_model::write_model(&mat, &names, cutoff, fix_resid, bifactor)
}

/// Compute model-implied genetic correlation matrix.
#[pyfunction]
#[pyo3(signature = (covstruc_json, estimation="DWLS"))]
fn rgmodel(covstruc_json: &str, estimation: &str) -> PyResult<String> {
    let ldsc_result = conversions::json_to_ldsc(covstruc_json)
        .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("invalid covstruc JSON"))?;
    let result = gsem_sem::rgmodel::run_rgmodel(&ldsc_result.s, &ldsc_result.v, estimation)
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
    Ok(())
}
