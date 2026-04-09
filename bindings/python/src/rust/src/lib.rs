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

/// Python module definition.
#[pymodule]
fn genomicsem(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyLdscResult>()?;
    m.add_function(wrap_pyfunction!(ldsc, m)?)?;
    m.add_function(wrap_pyfunction!(usermodel, m)?)?;
    m.add_function(wrap_pyfunction!(munge, m)?)?;
    Ok(())
}
