#![allow(clippy::too_many_arguments)]

use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::sync::atomic::{AtomicBool, Ordering};

mod conversions;

use conversions::{mat_to_pyarray, pyany_to_ldsc_result, pyarray_to_mat};

/// One-shot session flag for the commonfactor_gwas compatibility warning.
static COMMONFACTOR_GWAS_WARNED: AtomicBool = AtomicBool::new(false);

/// Emit a first-use warning that gsem.commonfactor_gwas uses a different
/// parameterization than R's GenomicSEM::commonfactorGWAS. Fires once per
/// process. Suppress by setting the GSEMR_COMMONFACTOR_GWAS_QUIET env var
/// to any non-empty value.
fn warn_commonfactor_gwas_semantics(py: Python<'_>) {
    if std::env::var("GSEMR_COMMONFACTOR_GWAS_QUIET")
        .map(|v| !v.is_empty())
        .unwrap_or(false)
    {
        return;
    }
    // Compare-and-set: return early on subsequent calls within the same
    // process. Ordering::Relaxed is sufficient because we don't need to
    // synchronize any other memory with this flag.
    if COMMONFACTOR_GWAS_WARNED
        .compare_exchange(false, true, Ordering::Relaxed, Ordering::Relaxed)
        .is_err()
    {
        return;
    }

    let msg = "gsem.commonfactor_gwas fits a single-factor model using \
fixed-variance identification (F1 ~~ 1*F1, loadings free) with the \
fix_measurement baseline optimization. This is numerically stable and \
matches GenomicSEM::userGWAS on the equivalent model, but it does NOT \
numerically match GenomicSEM::commonfactorGWAS, which uses marker-\
indicator identification and refits the full model per SNP. On real \
GWAS data the two can disagree in sign and magnitude. If you need bit-\
for-bit parity with R GenomicSEM::commonfactorGWAS, there is no exact \
replacement currently. If you want stable single-factor GWAS with the \
same model R userGWAS would fit, this function is the right call. See \
ARCHITECTURE.md section 3.3 for the full rationale. Suppress this \
warning by setting GSEMR_COMMONFACTOR_GWAS_QUIET=1 in the environment.";

    // Route through Python's warnings module so standard warning filters apply.
    if let Ok(warnings) = py.import("warnings") {
        let _ = warnings.call_method1("warn", (msg,));
    }
}

/// Clamp intercept matrix diagonal to >= 1.0 (GenomicSEM convention).
fn clamp_i_ld_diagonal(i_mat: &mut faer::Mat<f64>) {
    for i in 0..i_mat.nrows() {
        if i_mat[(i, i)] < 1.0 {
            i_mat[(i, i)] = 1.0;
        }
    }
}

// ---------------------------------------------------------------------------
// Two-table columnar result dicts for per-SNP / per-gene GWAS.
// ---------------------------------------------------------------------------
//
// Return shape:
//   {"snps":   {SNP, CHR, BP, MAF, A1, A2, chisq, df, converged[, Q_*]},
//    "params": {SNP, lhs, op, rhs, est, se, z, p}}
// TWAS uses "genes" with columns Gene/Panel/HSQ instead of "snps"; the
// params table's identifier column is "Gene". Both inner dicts are
// columnar and feed directly into `pandas.DataFrame(result["..."])`.
//
// The flat params scatter runs under rayon inside `py.detach()` so the
// GIL is released while Rust touches only Rust-owned buffers;
// `PyDict::new` / `set_item` calls happen only after reacquiring the
// GIL. Mirrors `bindings/r/src/rust/src/conversions.rs`.

/// Below this many total params, the serial flatten beats the rayon setup.
const PARAMS_PARALLEL_THRESHOLD: usize = 50_000;

/// Columnar buffers for the per-SNP table.
struct SnpColsPy {
    snp: Vec<String>,
    chr: Vec<Option<i64>>,
    bp: Vec<Option<i64>>,
    maf: Vec<f64>,
    a1: Vec<String>,
    a2: Vec<String>,
    chisq: Vec<f64>,
    df: Vec<i64>,
    converged: Vec<bool>,
}

/// Columnar buffers for the per-gene (TWAS) table.
struct GeneColsPy {
    gene: Vec<String>,
    panel: Vec<String>,
    hsq: Vec<f64>,
    chisq: Vec<f64>,
    df: Vec<i64>,
    converged: Vec<bool>,
}

/// Flat params table, sized to `total_params` so rayon can fill
/// disjoint windows.
struct ParamColsPy {
    id: Vec<String>,
    lhs: Vec<String>,
    op: Vec<String>,
    rhs: Vec<String>,
    est: Vec<f64>,
    se: Vec<f64>,
    z: Vec<f64>,
    p: Vec<f64>,
}

impl ParamColsPy {
    fn zeroed(total: usize) -> Self {
        Self {
            id: vec![String::new(); total],
            lhs: vec![String::new(); total],
            op: vec![String::new(); total],
            rhs: vec![String::new(); total],
            est: vec![0.0_f64; total],
            se: vec![0.0_f64; total],
            z: vec![0.0_f64; total],
            p: vec![0.0_f64; total],
        }
    }
}

/// `offsets[i]..offsets[i+1]` is result `i`'s row range in the flat
/// params table; `offsets[n] == total_params`.
fn params_offsets_py(results: &[gsem::gwas::user_gwas::SnpResult]) -> Vec<usize> {
    let n = results.len();
    let mut offsets = Vec::with_capacity(n + 1);
    let mut acc = 0usize;
    offsets.push(0);
    for r in results {
        acc += r.params.len();
        offsets.push(acc);
    }
    offsets
}

/// Write one result's params into its window of every column buffer.
#[inline]
fn fill_params_window_py(
    snp_name: &str,
    r: &gsem::gwas::user_gwas::SnpResult,
    id_chunk: &mut [String],
    lhs_chunk: &mut [String],
    op_chunk: &mut [String],
    rhs_chunk: &mut [String],
    est_chunk: &mut [f64],
    se_chunk: &mut [f64],
    z_chunk: &mut [f64],
    p_chunk: &mut [f64],
) {
    for (k, pr) in r.params.iter().enumerate() {
        id_chunk[k] = snp_name.to_string();
        lhs_chunk[k] = pr.lhs.clone();
        op_chunk[k] = pr.op.to_string();
        rhs_chunk[k] = pr.rhs.clone();
        est_chunk[k] = pr.est;
        se_chunk[k] = pr.se;
        z_chunk[k] = pr.z_stat;
        p_chunk[k] = pr.p_value;
    }
}

/// Carve the eight column buffers into non-aliasing per-result windows
/// via repeated `split_at_mut`.
fn split_param_windows_py<'a>(
    offsets: &[usize],
    id: &'a mut [String],
    lhs: &'a mut [String],
    op: &'a mut [String],
    rhs: &'a mut [String],
    est: &'a mut [f64],
    se: &'a mut [f64],
    z: &'a mut [f64],
    p: &'a mut [f64],
) -> Vec<(
    &'a mut [String],
    &'a mut [String],
    &'a mut [String],
    &'a mut [String],
    &'a mut [f64],
    &'a mut [f64],
    &'a mut [f64],
    &'a mut [f64],
)> {
    let n = offsets.len().saturating_sub(1);
    let mut out = Vec::with_capacity(n);
    let mut id_rest = id;
    let mut lhs_rest = lhs;
    let mut op_rest = op;
    let mut rhs_rest = rhs;
    let mut est_rest = est;
    let mut se_rest = se;
    let mut z_rest = z;
    let mut p_rest = p;
    for i in 0..n {
        let take = offsets[i + 1] - offsets[i];
        let (id_head, id_tail) = id_rest.split_at_mut(take);
        let (lhs_head, lhs_tail) = lhs_rest.split_at_mut(take);
        let (op_head, op_tail) = op_rest.split_at_mut(take);
        let (rhs_head, rhs_tail) = rhs_rest.split_at_mut(take);
        let (est_head, est_tail) = est_rest.split_at_mut(take);
        let (se_head, se_tail) = se_rest.split_at_mut(take);
        let (z_head, z_tail) = z_rest.split_at_mut(take);
        let (p_head, p_tail) = p_rest.split_at_mut(take);
        out.push((
            id_head, lhs_head, op_head, rhs_head, est_head, se_head, z_head, p_head,
        ));
        id_rest = id_tail;
        lhs_rest = lhs_tail;
        op_rest = op_tail;
        rhs_rest = rhs_tail;
        est_rest = est_tail;
        se_rest = se_tail;
        z_rest = z_tail;
        p_rest = p_tail;
    }
    out
}

/// Scatter-fill the params table, serial below
/// [`PARAMS_PARALLEL_THRESHOLD`] and rayon above.
fn fill_params_table_py(
    results: &[gsem::gwas::user_gwas::SnpResult],
    ids_per_snp: &[String],
    offsets: &[usize],
    pc: &mut ParamColsPy,
) {
    use rayon::prelude::*;

    let total = *offsets.last().unwrap_or(&0);
    if total < PARAMS_PARALLEL_THRESHOLD {
        for (i, r) in results.iter().enumerate() {
            let start = offsets[i];
            let end = offsets[i + 1];
            fill_params_window_py(
                &ids_per_snp[i],
                r,
                &mut pc.id[start..end],
                &mut pc.lhs[start..end],
                &mut pc.op[start..end],
                &mut pc.rhs[start..end],
                &mut pc.est[start..end],
                &mut pc.se[start..end],
                &mut pc.z[start..end],
                &mut pc.p[start..end],
            );
        }
        return;
    }

    let windows = split_param_windows_py(
        offsets, &mut pc.id, &mut pc.lhs, &mut pc.op, &mut pc.rhs, &mut pc.est, &mut pc.se,
        &mut pc.z, &mut pc.p,
    );
    windows
        .into_par_iter()
        .zip(results.par_iter())
        .zip(ids_per_snp.par_iter())
        .for_each(
            |(((id_c, lhs_c, op_c, rhs_c, est_c, se_c, z_c, p_c), r), snp_name)| {
                fill_params_window_py(
                    snp_name, r, id_c, lhs_c, op_c, rhs_c, est_c, se_c, z_c, p_c,
                );
            },
        );
}

/// Extract `Q_SNP` columns into three parallel vectors, `None`-filled
/// where missing. Returns `None` when no result carries any Q_SNP info.
fn q_snp_columns_py(
    results: &[gsem::gwas::user_gwas::SnpResult],
) -> Option<(Vec<Option<f64>>, Vec<Option<i64>>, Vec<Option<f64>>)> {
    let any = results
        .iter()
        .any(|r| r.q_snp.is_some() || r.q_snp_df.is_some() || r.q_snp_p.is_some());
    if !any {
        return None;
    }
    let chisq: Vec<Option<f64>> = results.iter().map(|r| r.q_snp).collect();
    let df: Vec<Option<i64>> = results
        .iter()
        .map(|r| r.q_snp_df.map(|d| d as i64))
        .collect();
    let pval: Vec<Option<f64>> = results.iter().map(|r| r.q_snp_p).collect();
    Some((chisq, df, pval))
}

/// Build the `{"snps": ..., "params": ...}` dict for `user_gwas` /
/// `commonfactor_gwas`. Column buffers are populated on rayon threads
/// inside `py.detach()`; `PyDict::new` / `set_item` calls happen after
/// the GIL is reacquired.
///
/// `index_map` lets callers that filter the merged file (e.g. MAF=0
/// drop in `commonfactor_gwas`) pass the original `MergedSumstats`
/// plus a `Vec<usize>` of surviving rows rather than a cloned subset.
/// Pass `None` for the identity mapping.
fn snp_results_to_dict<'py>(
    py: Python<'py>,
    results: &[gsem::gwas::user_gwas::SnpResult],
    merged: &gsem::io::sumstats_reader::MergedSumstats,
    index_map: Option<&[usize]>,
) -> PyResult<Bound<'py, PyDict>> {
    let (sc, pc, q_opt) = py.detach(|| {
        let n = results.len();

        let mut sc = SnpColsPy {
            snp: Vec::with_capacity(n),
            chr: Vec::with_capacity(n),
            bp: Vec::with_capacity(n),
            maf: Vec::with_capacity(n),
            a1: Vec::with_capacity(n),
            a2: Vec::with_capacity(n),
            chisq: Vec::with_capacity(n),
            df: Vec::with_capacity(n),
            converged: Vec::with_capacity(n),
        };
        for r in results {
            let idx = match index_map {
                Some(m) => m[r.snp_idx],
                None => r.snp_idx,
            };
            sc.snp.push(merged.snp[idx].clone());
            sc.chr.push(
                merged
                    .chr
                    .as_ref()
                    .and_then(|c| c.get(idx).copied())
                    .map(|c| c as i64),
            );
            sc.bp.push(
                merged
                    .bp
                    .as_ref()
                    .and_then(|b| b.get(idx).copied())
                    .map(|b| b as i64),
            );
            sc.maf.push(merged.maf[idx]);
            sc.a1.push(merged.a1_string(idx));
            sc.a2.push(merged.a2_string(idx));
            sc.chisq.push(r.chisq);
            sc.df.push(r.chisq_df as i64);
            sc.converged.push(r.converged);
        }

        // Reuse `sc.snp` for params-table ids — no separate clone.
        let offsets = params_offsets_py(results);
        let mut pc = ParamColsPy::zeroed(*offsets.last().unwrap_or(&0));
        fill_params_table_py(results, &sc.snp, &offsets, &mut pc);

        (sc, pc, q_snp_columns_py(results))
    });

    let snps_dict = PyDict::new(py);
    snps_dict.set_item("SNP", sc.snp)?;
    snps_dict.set_item("CHR", sc.chr)?;
    snps_dict.set_item("BP", sc.bp)?;
    snps_dict.set_item("MAF", sc.maf)?;
    snps_dict.set_item("A1", sc.a1)?;
    snps_dict.set_item("A2", sc.a2)?;
    snps_dict.set_item("chisq", sc.chisq)?;
    snps_dict.set_item("df", sc.df)?;
    snps_dict.set_item("converged", sc.converged)?;
    if let Some((qc, qd, qp)) = q_opt {
        snps_dict.set_item("Q_chisq", qc)?;
        snps_dict.set_item("Q_df", qd)?;
        snps_dict.set_item("Q_pval", qp)?;
    }

    let params_dict = PyDict::new(py);
    params_dict.set_item("SNP", pc.id)?;
    params_dict.set_item("lhs", pc.lhs)?;
    params_dict.set_item("op", pc.op)?;
    params_dict.set_item("rhs", pc.rhs)?;
    params_dict.set_item("est", pc.est)?;
    params_dict.set_item("se", pc.se)?;
    params_dict.set_item("z", pc.z)?;
    params_dict.set_item("p", pc.p)?;

    let out = PyDict::new(py);
    out.set_item("snps", snps_dict)?;
    out.set_item("params", params_dict)?;
    Ok(out)
}

/// TWAS variant of [`snp_results_to_dict`]. Returns
/// `{"genes": {...}, "params": {...}}` with `Gene` as the identifier
/// column in both inner tables. Called from `user_gwas(twas = True)`.
fn twas_results_to_dict<'py>(
    py: Python<'py>,
    results: &[gsem::gwas::user_gwas::SnpResult],
    genes: &[gsem::io::twas_reader::TwasGene],
) -> PyResult<Bound<'py, PyDict>> {
    let (gc, pc, q_opt) = py.detach(|| {
        let n = results.len();

        let mut gc = GeneColsPy {
            gene: Vec::with_capacity(n),
            panel: Vec::with_capacity(n),
            hsq: Vec::with_capacity(n),
            chisq: Vec::with_capacity(n),
            df: Vec::with_capacity(n),
            converged: Vec::with_capacity(n),
        };
        let mut ids_per_snp: Vec<String> = Vec::with_capacity(n);
        for r in results {
            let g = &genes[r.snp_idx];
            gc.gene.push(g.gene.clone());
            gc.panel.push(g.panel.clone());
            gc.hsq.push(g.hsq);
            gc.chisq.push(r.chisq);
            gc.df.push(r.chisq_df as i64);
            gc.converged.push(r.converged);
            ids_per_snp.push(g.gene.clone());
        }

        let offsets = params_offsets_py(results);
        let mut pc = ParamColsPy::zeroed(*offsets.last().unwrap_or(&0));
        fill_params_table_py(results, &ids_per_snp, &offsets, &mut pc);

        (gc, pc, q_snp_columns_py(results))
    });

    let genes_dict = PyDict::new(py);
    genes_dict.set_item("Gene", gc.gene)?;
    genes_dict.set_item("Panel", gc.panel)?;
    genes_dict.set_item("HSQ", gc.hsq)?;
    genes_dict.set_item("chisq", gc.chisq)?;
    genes_dict.set_item("df", gc.df)?;
    genes_dict.set_item("converged", gc.converged)?;
    if let Some((qc, qd, qp)) = q_opt {
        genes_dict.set_item("Q_chisq", qc)?;
        genes_dict.set_item("Q_df", qd)?;
        genes_dict.set_item("Q_pval", qp)?;
    }

    let params_dict = PyDict::new(py);
    params_dict.set_item("Gene", pc.id)?;
    params_dict.set_item("lhs", pc.lhs)?;
    params_dict.set_item("op", pc.op)?;
    params_dict.set_item("rhs", pc.rhs)?;
    params_dict.set_item("est", pc.est)?;
    params_dict.set_item("se", pc.se)?;
    params_dict.set_item("z", pc.z)?;
    params_dict.set_item("p", pc.p)?;

    let out = PyDict::new(py);
    out.set_item("genes", genes_dict)?;
    out.set_item("params", params_dict)?;
    Ok(out)
}

/// Build a columnar dict from a slice of `SnpParamResult`s — used by
/// `multi_snp` / `multi_gene` where there's a single parameter table
/// rather than per-SNP grouping.
fn param_results_to_dict<'py>(
    py: Python<'py>,
    params: &[gsem::gwas::user_gwas::SnpParamResult],
) -> PyResult<Bound<'py, PyDict>> {
    let n = params.len();
    let mut lhs: Vec<String> = Vec::with_capacity(n);
    let mut op: Vec<String> = Vec::with_capacity(n);
    let mut rhs: Vec<String> = Vec::with_capacity(n);
    let mut est: Vec<f64> = Vec::with_capacity(n);
    let mut se: Vec<f64> = Vec::with_capacity(n);
    let mut z: Vec<f64> = Vec::with_capacity(n);
    let mut p: Vec<f64> = Vec::with_capacity(n);
    for r in params {
        lhs.push(r.lhs.clone());
        op.push(r.op.to_string());
        rhs.push(r.rhs.clone());
        est.push(r.est);
        se.push(r.se);
        z.push(r.z_stat);
        p.push(r.p_value);
    }
    let out = PyDict::new(py);
    out.set_item("lhs", lhs)?;
    out.set_item("op", op)?;
    out.set_item("rhs", rhs)?;
    out.set_item("est", est)?;
    out.set_item("se", se)?;
    out.set_item("z", z)?;
    out.set_item("p", p)?;
    Ok(out)
}

/// Columnar dict for SEM-fit `ParamEstimate`s (used by `commonfactor` /
/// `usermodel` whose engine returns an SE-aware parameter table).
fn param_estimates_to_dict<'py>(
    py: Python<'py>,
    params: &[gsem_sem::ParamEstimate],
) -> PyResult<Bound<'py, PyDict>> {
    let n = params.len();
    let mut lhs: Vec<String> = Vec::with_capacity(n);
    let mut op: Vec<String> = Vec::with_capacity(n);
    let mut rhs: Vec<String> = Vec::with_capacity(n);
    let mut est: Vec<f64> = Vec::with_capacity(n);
    let mut se: Vec<f64> = Vec::with_capacity(n);
    let mut z: Vec<f64> = Vec::with_capacity(n);
    let mut p: Vec<f64> = Vec::with_capacity(n);
    for r in params {
        lhs.push(r.lhs.clone());
        op.push(r.op.to_string());
        rhs.push(r.rhs.clone());
        est.push(r.est);
        se.push(r.se);
        z.push(r.z);
        p.push(r.p);
    }
    let out = PyDict::new(py);
    out.set_item("lhs", lhs)?;
    out.set_item("op", op)?;
    out.set_item("rhs", rhs)?;
    out.set_item("est", est)?;
    out.set_item("se", se)?;
    out.set_item("z", z)?;
    out.set_item("p", p)?;
    Ok(out)
}

/// Python wrapper for LDSC result.
#[pyclass(name = "LdscResult")]
struct PyLdscResult {
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
    fn s<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let arr = Array2::from_shape_vec(self.s_shape, self.s_data.clone())
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("S shape: {e}")))?;
        Ok(arr.into_pyarray(py))
    }

    /// Sampling covariance matrix V as NumPy array.
    #[getter]
    fn v<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let arr = Array2::from_shape_vec(self.v_shape, self.v_data.clone())
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("V shape: {e}")))?;
        Ok(arr.into_pyarray(py))
    }

    /// Intercept matrix I as NumPy array.
    #[getter]
    fn i_mat<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let arr = Array2::from_shape_vec(self.i_shape, self.i_data.clone())
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("I shape: {e}")))?;
        Ok(arr.into_pyarray(py))
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
    fn s_stand<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray2<f64>>>> {
        match (&self.s_stand_data, self.s_stand_shape) {
            (Some(data), Some(shape)) => {
                let arr = Array2::from_shape_vec(shape, data.clone()).map_err(|e| {
                    pyo3::exceptions::PyValueError::new_err(format!("S_Stand shape: {e}"))
                })?;
                Ok(Some(arr.into_pyarray(py)))
            }
            _ => Ok(None),
        }
    }

    /// Standardized sampling covariance (only when stand=True).
    #[getter]
    fn v_stand<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray2<f64>>>> {
        match (&self.v_stand_data, self.v_stand_shape) {
            (Some(data), Some(shape)) => {
                let arr = Array2::from_shape_vec(shape, data.clone()).map_err(|e| {
                    pyo3::exceptions::PyValueError::new_err(format!("V_Stand shape: {e}"))
                })?;
                Ok(Some(arr.into_pyarray(py)))
            }
            _ => Ok(None),
        }
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
        let s_vec = gsem_matrix::vech::vech(&result.s).expect("S must be square");
        let ss_vec = gsem_matrix::vech::vech(&s_stand).expect("S_Stand must be square");
        let scale: Vec<f64> = ss_vec
            .iter()
            .zip(s_vec.iter())
            .enumerate()
            .map(|(i, (&st, &orig))| {
                let ratio = if orig.abs() > 1e-30 { st / orig } else { 0.0 };
                result.v[(i, i)].sqrt() * ratio
            })
            .collect();
        let v_cor = gsem_matrix::smooth::cov_to_cor(&result.v);
        let v_stand = faer::Mat::from_fn(kstar, kstar, |i, j| scale[i] * v_cor[(i, j)] * scale[j]);
        let (ss_data, ss_r, ss_c) = conversions::mat_to_flat(&s_stand);
        let (vs_data, vs_r, vs_c) = conversions::mat_to_flat(&v_stand);
        (
            Some(ss_data),
            Some((ss_r, ss_c)),
            Some(vs_data),
            Some((vs_r, vs_c)),
        )
    } else {
        (None, None, None, None)
    };

    PyLdscResult {
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
///
/// Estimates a genetic covariance matrix (and sampling-covariance matrix)
/// across a set of munged GWAS summary statistics using LD score
/// regression with a block jackknife.
///
/// Examples
/// --------
/// >>> import genomicsem as gs
/// >>> cov = gs.ldsc(
/// ...     traits=["T1.sumstats.gz", "T2.sumstats.gz"],
/// ...     sample_prev=[None, None],
/// ...     population_prev=[None, None],
/// ...     ld="eur_w_ld_chr/",
/// ...     trait_names=["T1", "T2"],
/// ... )
/// >>> cov.s, cov.v  # 2x2 genetic covariance, 3x3 sampling covariance
#[pyfunction]
#[pyo3(signature = (traits, sample_prev, population_prev, ld, wld="", trait_names=None, sep_weights=false, chr=22, n_blocks=200, ldsc_log=None, stand=false, select=None, chisq_max=None, parallel=true, cores=None))]
#[allow(unused_variables)]
fn ldsc(
    py: Python<'_>,
    traits: Vec<String>,
    sample_prev: Vec<Option<f64>>,
    population_prev: Vec<Option<f64>>,
    ld: &str,
    wld: &str,
    trait_names: Option<Vec<String>>,
    sep_weights: bool, // ignored
    chr: usize,
    n_blocks: usize,
    ldsc_log: Option<String>, // ignored
    stand: bool,
    select: Option<String>,
    chisq_max: Option<f64>,
    parallel: bool,
    cores: Option<usize>,
) -> PyResult<PyLdscResult> {
    if sep_weights {
        log::info!(
            "sep_weights is always enabled in gsemr — weight LD scores are read from the wld directory"
        );
    }
    if let Some(ref path) = ldsc_log {
        log::info!("ldsc_log='{path}' — file logging is handled at the Python/R wrapper level");
    }
    let wld_dir = if wld.is_empty() { ld } else { wld };

    let trait_data = gsem::io::gwas_reader::load_trait_data(&traits)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{e}")))?;

    let ld_path = std::path::Path::new(ld);
    let wld_path = std::path::Path::new(wld_dir);
    let chromosomes: Vec<usize> = match select.as_deref() {
        None | Some("") | Some("FALSE") | Some("false") => (1..=chr).collect(),
        Some("ODD") | Some("odd") => (1..=chr).filter(|c| c % 2 == 1).collect(),
        Some("EVEN") | Some("even") => (1..=chr).filter(|c| c % 2 == 0).collect(),
        Some(other) => other
            .split(',')
            .filter_map(|s| s.trim().parse().ok())
            .collect(),
    };
    let ld_data = gsem::io::ld_reader::read_ld_scores(ld_path, wld_path, &chromosomes)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{e}")))?;

    let ld_snps: Vec<String> = ld_data.records.iter().map(|r| r.snp.clone()).collect();
    let ld_scores: Vec<f64> = ld_data.records.iter().map(|r| r.l2).collect();

    let config = gsem_ldsc::LdscConfig {
        n_blocks,
        chisq_max,
        num_threads: if parallel { cores } else { Some(1) },
    };

    let n_pairs = trait_data.len() * (trait_data.len() + 1) / 2;
    let counter = std::sync::atomic::AtomicUsize::new(0);
    let interval = if n_pairs <= 20 {
        1
    } else {
        (n_pairs / 20).max(1)
    };
    let cb = || {
        let done = counter.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
        if done == n_pairs || done.is_multiple_of(interval) {
            log::info!("LDSC progress: {done}/{n_pairs} trait pairs");
        }
    };

    // Release the GIL while rayon runs the parallel jackknife. Without this,
    // worker threads that try to log via pyo3_log block on the GIL while the
    // main thread sits inside this FFI call holding the GIL — classic deadlock.
    let result = py
        .detach(|| {
            gsem_ldsc::ldsc(
                &trait_data,
                &sample_prev,
                &population_prev,
                &ld_scores,
                &ld_data.w_ld,
                &ld_snps,
                ld_data.total_m,
                &config,
                Some(&cb),
            )
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;
    Ok(ldsc_result_to_py(&result, stand))
}

/// Munge GWAS summary statistics files.
///
/// Reads one or more raw GWAS files, harmonises columns against an HM3
/// SNPlist, applies QC filters, and writes `<trait>.sumstats.gz` files
/// ready for `ldsc()`.
///
/// Examples
/// --------
/// >>> import genomicsem as gs
/// >>> gs.munge(
/// ...     files=["raw_T1.txt.gz", "raw_T2.txt.gz"],
/// ...     hm3="w_hm3.snplist",
/// ...     trait_names=["T1", "T2"],
/// ...     n=1e5,
/// ... )
/// ['./T1.sumstats.gz', './T2.sumstats.gz']
#[pyfunction]
#[pyo3(signature = (files, hm3, trait_names=None, n=None, info_filter=0.9, maf_filter=0.01, log_name=None, column_names=None, parallel=false, cores=None, overwrite=true, out="."))]
#[allow(unused_variables)]
fn munge(
    files: Vec<String>,
    hm3: &str,
    trait_names: Option<Vec<String>>,
    n: Option<f64>, // threaded (N override)
    info_filter: f64,
    maf_filter: f64,
    log_name: Option<String>, // ignored
    column_names: Option<std::collections::HashMap<String, String>>, // threaded
    parallel: bool,           // ignored
    cores: Option<usize>,     // ignored
    overwrite: bool,          // ignored
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
///
/// `covstruc` accepts either a `PyLdscResult` returned from `ldsc()` or a
/// dict with `s`/`v`/`i_mat`/`n_vec`/`m` (or their uppercase aliases)
/// entries.
///
/// Examples
/// --------
/// >>> import numpy as np, genomicsem as gs
/// >>> covstruc = {
/// ...     "s": np.array([[0.60, 0.42, 0.35],
/// ...                    [0.42, 0.50, 0.30],
/// ...                    [0.35, 0.30, 0.40]]),
/// ...     "v": np.eye(6) * 1e-3,
/// ...     "i_mat": np.eye(3),
/// ...     "n_vec": [1e5, 1e5, 1e5],
/// ...     "m": 1e6,
/// ... }
/// >>> fit = gs.usermodel(
/// ...     covstruc,
/// ...     model="F1 =~ NA*V1 + V2 + V3\nF1 ~~ 1*F1",
/// ...     estimation="DWLS",
/// ... )
/// >>> fit["converged"], fit["parameters"]
#[pyfunction]
#[pyo3(signature = (covstruc, model="", estimation="DWLS", cfi_calc=true, std_lv=false, imp_cov=false, fix_resid=true, toler=None, q_factor=false))]
#[allow(unused_variables)]
fn usermodel<'py>(
    py: Python<'py>,
    covstruc: &Bound<'py, PyAny>,
    model: &str,
    estimation: &str,
    cfi_calc: bool,     // ignored (CFIcalc)
    std_lv: bool,       // threaded
    imp_cov: bool,      // ignored
    fix_resid: bool,    // threaded
    toler: Option<f64>, // ignored
    q_factor: bool,     // ignored (Q_Factor)
) -> PyResult<Bound<'py, PyDict>> {
    let ldsc_result = pyany_to_ldsc_result(covstruc)?;

    let pt = gsem_sem::syntax::parse_model(model, std_lv)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("model parse error: {e}")))?;

    let k = ldsc_result.s.nrows();
    let obs_names: Vec<String> = (0..k).map(|i| format!("V{}", i + 1)).collect();
    let mut sem_model = gsem_sem::model::Model::from_partable(&pt, &obs_names);

    let kstar = k * (k + 1) / 2;
    let v_diag: Vec<f64> = (0..kstar).map(|i| ldsc_result.v[(i, i)]).collect();

    let est_method = gsem_sem::EstimationMethod::from_str_lossy(estimation);
    let fit = match est_method {
        gsem_sem::EstimationMethod::Ml => {
            gsem_sem::estimator::fit_ml(&mut sem_model, &ldsc_result.s, 1000, None)
        }
        gsem_sem::EstimationMethod::Dwls => {
            gsem_sem::estimator::fit_dwls(&mut sem_model, &ldsc_result.s, &v_diag, 1000, None)
        }
    };

    // Build a columnar parameter dict: {lhs, op, rhs, est}. SE/Z/p are not
    // computed on this code path (matches the R binding's usermodel shape).
    let free_rows: Vec<_> = pt.rows.iter().filter(|r| r.free > 0).collect();
    let mut lhs: Vec<String> = Vec::with_capacity(free_rows.len());
    let mut op: Vec<String> = Vec::with_capacity(free_rows.len());
    let mut rhs: Vec<String> = Vec::with_capacity(free_rows.len());
    let mut est: Vec<f64> = Vec::with_capacity(free_rows.len());
    for (i, row) in free_rows.iter().enumerate() {
        lhs.push(row.lhs.clone());
        op.push(row.op.to_string());
        rhs.push(row.rhs.clone());
        est.push(fit.params.get(i).copied().unwrap_or(0.0));
    }
    let params = PyDict::new(py);
    params.set_item("lhs", lhs)?;
    params.set_item("op", op)?;
    params.set_item("rhs", rhs)?;
    params.set_item("est", est)?;

    let out = PyDict::new(py);
    out.set_item("converged", fit.converged)?;
    out.set_item("objective", fit.objective)?;
    out.set_item("parameters", params)?;
    Ok(out)
}

/// Fit common factor model.
///
/// Auto-generates a 1-factor CFA from `covstruc` and fits it via DWLS
/// (default) or ML. Returns parameter estimates plus fit indices.
///
/// Examples
/// --------
/// >>> import numpy as np, genomicsem as gs
/// >>> covstruc = {
/// ...     "s": np.array([[0.60, 0.42, 0.35],
/// ...                    [0.42, 0.50, 0.30],
/// ...                    [0.35, 0.30, 0.40]]),
/// ...     "v": np.eye(6) * 1e-3,
/// ...     "i_mat": np.eye(3),
/// ...     "n_vec": [1e5, 1e5, 1e5],
/// ...     "m": 1e6,
/// ... }
/// >>> cf = gs.commonfactor(covstruc, estimation="DWLS")
/// >>> cf["chisq"], cf["df"], cf["cfi"]
#[pyfunction]
#[pyo3(signature = (covstruc, estimation="DWLS"))]
fn commonfactor<'py>(
    py: Python<'py>,
    covstruc: &Bound<'py, PyAny>,
    estimation: &str,
) -> PyResult<Bound<'py, PyDict>> {
    let ldsc_result = pyany_to_ldsc_result(covstruc)?;

    let result = gsem_sem::commonfactor::run_commonfactor(
        &ldsc_result.s,
        &ldsc_result.v,
        gsem_sem::EstimationMethod::from_str_lossy(estimation),
    )
    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;

    let out = PyDict::new(py);
    out.set_item(
        "parameters",
        param_estimates_to_dict(py, &result.parameters)?,
    )?;
    out.set_item("chisq", result.fit.chisq)?;
    out.set_item("df", result.fit.df as i64)?;
    out.set_item("p_chisq", result.fit.p_chisq)?;
    out.set_item("aic", result.fit.aic)?;
    out.set_item("cfi", result.fit.cfi)?;
    out.set_item("srmr", result.fit.srmr)?;
    Ok(out)
}

/// Merge GWAS summary statistics.
///
/// Returns a dict `{"path": out, "n_snps": n}` describing the written file.
///
/// `parallel`/`cores` control the local rayon pool that fans the
/// reference + per-trait reads across worker threads. `parallel=False`
/// forces a single-threaded run; otherwise `cores` caps the pool size
/// (default = rayon default, typically one thread per logical core).
///
/// Examples
/// --------
/// >>> import genomicsem as gs
/// >>> gs.sumstats(
/// ...     files=["T1.sumstats.gz", "T2.sumstats.gz"],
/// ...     ref_dir="eur_w_ld_chr/",
/// ...     trait_names=["T1", "T2"],
/// ...     out="merged_sumstats.tsv",
/// ... )
/// {'path': 'merged_sumstats.tsv', 'n_snps': ...}
#[pyfunction]
#[pyo3(signature = (files, ref_dir, trait_names=None, se_logit=None, ols=None, linprob=None, n=None, betas=None, info_filter=0.6, maf_filter=0.01, keep_indel=false, parallel=true, cores=None, ambig=false, direct_filter=false, out="merged_sumstats.tsv"))]
#[allow(unused_variables)]
fn sumstats<'py>(
    py: Python<'py>,
    files: Vec<String>,
    ref_dir: &str,
    trait_names: Option<Vec<String>>,
    se_logit: Option<Vec<bool>>, // threaded
    ols: Option<Vec<bool>>,      // threaded (OLS)
    linprob: Option<Vec<bool>>,  // threaded
    n: Option<Vec<f64>>,         // threaded (N overrides)
    betas: Option<Vec<f64>>,     // ignored
    info_filter: f64,
    maf_filter: f64,
    keep_indel: bool,
    parallel: bool,
    cores: Option<usize>,
    ambig: bool,          // ignored
    direct_filter: bool,  // ignored
    out: &str,
) -> PyResult<Bound<'py, PyDict>> {
    let k = files.len();
    let default_names: Vec<String> = (0..k).map(|i| format!("trait{}", i + 1)).collect();
    let names = trait_names.unwrap_or(default_names);

    // Mirror the R wrapper's resolution: cores > 0 wins, else
    // parallel=False means one thread, else None = rayon default.
    let num_threads: Option<usize> = if let Some(n) = cores.filter(|&n| n > 0) {
        Some(n)
    } else if !parallel {
        Some(1)
    } else {
        None
    };

    let config = gsem::sumstats::SumstatsConfig {
        info_filter,
        maf_filter,
        keep_indel,
        keep_ambig: ambig,
        se_logit: se_logit.unwrap_or_else(|| vec![false; k]),
        ols: ols.unwrap_or_else(|| vec![false; k]),
        linprob: linprob.unwrap_or_else(|| vec![false; k]),
        n_overrides: n
            .map(|v| v.into_iter().map(Some).collect())
            .unwrap_or_else(|| vec![None; k]),
        beta_overrides: Vec::new(),
        direct_filter: false,
        num_threads,
    };
    let file_refs: Vec<&std::path::Path> = files
        .iter()
        .map(|p| std::path::Path::new(p.as_str()))
        .collect();
    let n_snps = gsem::sumstats::merge_sumstats(
        &file_refs,
        std::path::Path::new(ref_dir),
        &names,
        &config,
        std::path::Path::new(out),
    )
    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;

    let d = PyDict::new(py);
    d.set_item("path", out)?;
    d.set_item("n_snps", n_snps as i64)?;
    Ok(d)
}

/// Run common factor GWAS.
///
/// Fits a 1-factor model with the SNP regressed on the common factor
/// at every SNP in a merged sumstats file. Returns a dict with `snps`
/// (per-SNP metadata + fit stats) and `params` (flat parameter table).
///
/// Warnings
/// --------
/// `gsem.commonfactor_gwas` uses an auto-generated fixed-variance
/// parameterisation that differs from R's `GenomicSEM::commonfactorGWAS`,
/// which uses a marker-indicator approach. Results are not bit-identical.
/// Set `GSEMR_COMMONFACTOR_GWAS_QUIET=1` in the environment to suppress
/// the first-use warning.
///
/// Examples
/// --------
/// >>> import genomicsem as gs
/// >>> covstruc = gs.ldsc(
/// ...     traits=["T1.sumstats.gz", "T2.sumstats.gz", "T3.sumstats.gz"],
/// ...     sample_prev=[None]*3, population_prev=[None]*3,
/// ...     ld="eur_w_ld_chr/",
/// ... )
/// >>> res = gs.commonfactor_gwas(
/// ...     covstruc=covstruc,
/// ...     snps="merged_sumstats.tsv",
/// ...     estimation="DWLS",
/// ... )
/// >>> res["snps"][:5]  # head of per-SNP table
#[pyfunction]
#[pyo3(signature = (
    covstruc,
    sumstats_path,
    estimation="DWLS",
    cores=None,
    toler=false,
    snpse=false,
    parallel=true,
    gc="standard",
    mpi=false,
    twas=false,
    smooth_check=false,
    identification="fixed_variance",
))]
#[allow(unused_variables)]
fn commonfactor_gwas<'py>(
    py: Python<'py>,
    covstruc: &Bound<'py, PyAny>,
    sumstats_path: &str,
    estimation: &str,
    cores: Option<usize>,
    toler: bool, // ignored: convergence is controlled internally
    snpse: bool,
    parallel: bool,
    gc: &str,
    mpi: bool, // ignored: not applicable to the Rust backend
    twas: bool,
    smooth_check: bool,
    identification: &str,
) -> PyResult<Bound<'py, PyDict>> {
    warn_commonfactor_gwas_semantics(py);

    if twas {
        return Err(pyo3::exceptions::PyNotImplementedError::new_err(
            "twas=True is not supported by the Python commonfactor_gwas binding yet; \
             use the R binding or open an issue if you need this",
        ));
    }

    let ldsc_result = pyany_to_ldsc_result(covstruc)?;
    let merged =
        gsem::io::sumstats_reader::read_merged_sumstats(std::path::Path::new(sumstats_path))
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{e}")))?;

    let gc_mode: gsem::gwas::gc_correction::GcMode = gc
        .parse()
        .unwrap_or(gsem::gwas::gc_correction::GcMode::Standard);
    let mut i_ld = ldsc_result.i_mat.to_owned();
    clamp_i_ld_diagonal(&mut i_ld);

    // Filter SNPs with zero MAF: var_snp=0 → singular per-SNP matrices.
    // Work against the SoA `maf` column directly to avoid any
    // per-SNP cloning.
    let valid_idx: Vec<usize> = (0..merged.len())
        .filter(|&i| {
            let m = merged.maf[i];
            2.0 * m * (1.0 - m) > 1e-10
        })
        .collect();
    if valid_idx.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "no SNPs with valid MAF (all MAF=0)",
        ));
    }
    // Build ref vecs directly from `valid_idx` — no cloned subset.
    let beta_snp: Vec<&[f64]> = valid_idx.iter().map(|&i| merged.beta_row(i)).collect();
    let se_snp: Vec<&[f64]> = valid_idx.iter().map(|&i| merged.se_row(i)).collect();
    let var_snp: Vec<f64> = valid_idx
        .iter()
        .map(|&i| {
            let m = merged.maf[i];
            2.0 * m * (1.0 - m)
        })
        .collect();

    let cf_config = gsem::gwas::common_factor::CommonFactorGwasConfig {
        estimation: gsem_sem::EstimationMethod::from_str_lossy(estimation),
        gc: gc_mode,
        snp_se: if snpse { Some(0.0005) } else { None },
        smooth_check,
        identification: gsem::gwas::common_factor::Identification::from_str_lossy(identification),
        num_threads: if parallel { cores } else { Some(1) },
        ..Default::default()
    };

    // Release the GIL while rayon runs the per-SNP fits.
    let results = py.detach(|| {
        gsem::gwas::common_factor::run_common_factor_gwas(
            &merged.trait_names,
            &ldsc_result.s,
            &ldsc_result.v,
            &i_ld,
            &beta_snp,
            &se_snp,
            &var_snp,
            &cf_config,
            None,
        )
    });

    // Drop the ref vecs before building the result dict.
    drop(beta_snp);
    drop(se_snp);
    drop(var_snp);

    snp_results_to_dict(py, &results, &merged, Some(&valid_idx))
}

/// Run user-specified GWAS.
///
/// Fits a user-specified SEM at every SNP. Returns a dict with `snps`
/// (per-SNP metadata + fit stats) and `params` (flat parameter table).
/// Filter by SNP and operator to extract per-SNP regression effects.
///
/// Examples
/// --------
/// >>> import genomicsem as gs
/// >>> res = gs.user_gwas(
/// ...     covstruc=covstruc,
/// ...     snps="merged_sumstats.tsv",
/// ...     model="F1 =~ NA*V1 + V2 + V3\nF1 ~ SNP\nF1 ~~ 1*F1",
/// ... )
/// >>> # every SNP's F1 ~ SNP regression effect:
/// >>> [(s, e) for s, o, r, e in zip(res["params"]["SNP"],
/// ...                               res["params"]["op"],
/// ...                               res["params"]["rhs"],
/// ...                               res["params"]["est"])
/// ...  if o == "~" and r == "SNP"][:5]
#[pyfunction]
#[pyo3(signature = (covstruc, sumstats_path, model="", estimation="DWLS", printwarn=true, sub=None, cores=None, toler=false, snpse=false, parallel=true, gc="standard", mpi=false, smooth_check=false, twas=false, std_lv=false, fix_measurement=true, q_snp=false))]
#[allow(unused_variables)]
fn user_gwas<'py>(
    py: Python<'py>,
    covstruc: &Bound<'py, PyAny>,
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
) -> PyResult<Bound<'py, PyDict>> {
    let ldsc_result = pyany_to_ldsc_result(covstruc)?;
    let k = ldsc_result.s.nrows();
    let gc_mode: gsem::gwas::gc_correction::GcMode = gc
        .parse()
        .unwrap_or(gsem::gwas::gc_correction::GcMode::Standard);
    let mut i_ld = ldsc_result.i_mat.to_owned();
    clamp_i_ld_diagonal(&mut i_ld);

    let snp_se_val = if snpse { Some(0.0005) } else { None };

    // TWAS renames the observed variant from `SNP` to `Gene` throughout
    // the model; see crates/gsem/src/main.rs::run_twas_gwas.
    let model_str = if twas {
        model.replace("SNP", "Gene")
    } else {
        model.to_string()
    };
    let pt = gsem_sem::syntax::parse_model(&model_str, std_lv)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("model parse error: {e}")))?;

    let variant_label = if twas {
        gsem::gwas::user_gwas::VariantLabel::Gene
    } else {
        gsem::gwas::user_gwas::VariantLabel::Snp
    };
    let config = gsem::gwas::user_gwas::UserGwasConfig {
        model: pt,
        estimation: gsem_sem::EstimationMethod::from_str_lossy(estimation),
        gc: gc_mode,
        max_iter: 500,
        smooth_check,
        snp_se: snp_se_val,
        variant_label,
        q_snp,
        fix_measurement,
        num_threads: if parallel { cores } else { Some(1) },
    };

    if !printwarn {
        log::set_max_level(log::LevelFilter::Error);
    }
    if toler {
        log::info!(
            "toler — convergence tolerance is controlled by the L-BFGS optimizer internally"
        );
    }
    if mpi {
        log::warn!("MPI is not supported in gsemr — use the cores parameter for thread control");
    }

    // Apply the `sub` parameter-name filter in place.
    let apply_sub_filter = |results: &mut Vec<gsem::gwas::user_gwas::SnpResult>| {
        if let Some(ref patterns) = sub {
            let pats: Vec<String> = patterns
                .iter()
                .map(|s| s.trim().replace(' ', ""))
                .filter(|s| !s.is_empty())
                .collect();
            if !pats.is_empty() {
                for snp_result in results.iter_mut() {
                    snp_result.params.retain(|p| {
                        let key = format!("{}{}{}", p.lhs, p.op, p.rhs).replace(' ', "");
                        pats.contains(&key)
                    });
                }
            }
        }
    };

    if twas {
        let twas_data = gsem::io::twas_reader::read_twas_sumstats(std::path::Path::new(
            sumstats_path,
        ))
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{e}")))?;

        if twas_data.trait_names.len() != k {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "trait count mismatch: LDSC has {k} traits, TWAS sumstats has {} (beta.* columns)",
                twas_data.trait_names.len()
            )));
        }

        let beta_gene: Vec<&[f64]> = twas_data.genes.iter().map(|g| g.beta.as_slice()).collect();
        let se_gene: Vec<&[f64]> = twas_data.genes.iter().map(|g| g.se.as_slice()).collect();
        // In TWAS mode var_snp = HSQ (heritability of expression).
        let var_gene: Vec<f64> = twas_data.genes.iter().map(|g| g.hsq).collect();

        let mut results = py.detach(|| {
            gsem::gwas::user_gwas::run_user_gwas(
                &config,
                &ldsc_result.s,
                &ldsc_result.v,
                &i_ld,
                &beta_gene,
                &se_gene,
                &var_gene,
                None,
            )
        });
        apply_sub_filter(&mut results);

        twas_results_to_dict(py, &results, &twas_data.genes)
    } else {
        let merged =
            gsem::io::sumstats_reader::read_merged_sumstats(std::path::Path::new(sumstats_path))
                .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{e}")))?;

        let beta_snp = merged.beta_rows();
        let se_snp = merged.se_rows();
        let var_snp = merged.var_snp();

        let mut results = py.detach(|| {
            gsem::gwas::user_gwas::run_user_gwas(
                &config,
                &ldsc_result.s,
                &ldsc_result.v,
                &i_ld,
                &beta_snp,
                &se_snp,
                &var_snp,
                None,
            )
        });
        apply_sub_filter(&mut results);

        drop(beta_snp);
        drop(se_snp);
        drop(var_snp);

        snp_results_to_dict(py, &results, &merged, None)
    }
}

/// Parallel analysis to determine number of factors.
///
/// Compares observed eigenvalues of the genetic covariance matrix to
/// the 95th percentile of eigenvalues simulated from sampling
/// covariances, returning the recommended number of factors.
///
/// Examples
/// --------
/// >>> import numpy as np, genomicsem as gs
/// >>> s = np.array([[0.60, 0.42, 0.35],
/// ...               [0.42, 0.50, 0.30],
/// ...               [0.35, 0.30, 0.40]])
/// >>> v = np.eye(6) * 1e-3
/// >>> pa = gs.parallel_analysis(s, v, r=500)
/// >>> pa["n_factors"], pa["observed"]
#[pyfunction]
#[pyo3(signature = (s, v, r=500, p=None, diag=false, parallel=true, cores=None))]
fn parallel_analysis<'py>(
    py: Python<'py>,
    s: PyReadonlyArray2<'py, f64>,
    v: PyReadonlyArray2<'py, f64>,
    r: usize,
    p: Option<f64>,
    diag: bool,
    parallel: bool,
    cores: Option<usize>,
) -> PyResult<Bound<'py, PyDict>> {
    let s_mat = pyarray_to_mat(&s);
    let v_mat = pyarray_to_mat(&v);

    let percentile = p.unwrap_or(0.95);
    let num_threads = if parallel { cores } else { Some(1) };
    // Defensive: release the GIL while rayon runs simulations. Today the
    // par_iter has no log calls, but adding one later would deadlock.
    let result = py.detach(|| {
        gsem::stats::parallel_analysis::parallel_analysis(
            &s_mat,
            &v_mat,
            r,
            percentile,
            diag,
            num_threads,
            None,
        )
    });

    let out = PyDict::new(py);
    out.set_item("observed", result.observed)?;
    out.set_item("simulated_95", result.simulated_95)?;
    out.set_item("n_factors", result.n_factors as i64)?;
    Ok(out)
}

/// Auto-generate model syntax from factor loadings.
///
/// Converts a factor loading matrix into lavaan-style model syntax,
/// dropping loadings below `cutoff` and optionally fixing residual
/// variances to be positive.
///
/// Examples
/// --------
/// >>> import genomicsem as gs
/// >>> syntax = gs.write_model(
/// ...     loadings=[[0.80], [0.60], [0.70]],
/// ...     names=["V1", "V2", "V3"],
/// ...     cutoff=0.3,
/// ... )
/// >>> print(syntax)
/// F1 =~ ...
#[pyfunction]
#[pyo3(signature = (loadings, names, cutoff=0.3, fix_resid=true, bifactor=false, mustload=false, common=false))]
#[allow(unused_variables)]
fn write_model(
    loadings: Vec<Vec<f64>>,
    names: Vec<String>,
    cutoff: f64,
    fix_resid: bool,
    bifactor: bool,
    mustload: bool, // ignored
    common: bool,   // ignored
) -> String {
    let n_rows = loadings.len();
    let n_cols = if n_rows > 0 { loadings[0].len() } else { 0 };
    let mat = faer::Mat::from_fn(n_rows, n_cols, |i, j| loadings[i][j]);
    gsem_sem::write_model::write_model(&mat, &names, cutoff, fix_resid, bifactor, mustload, common)
}

/// Compute model-implied genetic correlation matrix.
/// `estimation=True` means DWLS; `estimation=False` means ML.
///
/// Examples
/// --------
/// >>> import numpy as np, genomicsem as gs
/// >>> covstruc = {
/// ...     "s": np.array([[0.60, 0.42, 0.35],
/// ...                    [0.42, 0.50, 0.30],
/// ...                    [0.35, 0.30, 0.40]]),
/// ...     "v": np.eye(6) * 1e-3,
/// ...     "i_mat": np.eye(3),
/// ...     "n_vec": [1e5, 1e5, 1e5],
/// ...     "m": 1e6,
/// ... }
/// >>> rg = gs.rgmodel(covstruc, model="", std_lv=True, estimation=True)
/// >>> rg["rg"]  # model-implied genetic correlation matrix
#[pyfunction]
#[pyo3(signature = (covstruc, model, std_lv=true, estimation=true, sub=None))]
#[allow(unused_variables)]
fn rgmodel<'py>(
    py: Python<'py>,
    covstruc: &Bound<'py, PyAny>,
    model: &str,
    std_lv: bool,
    estimation: bool,
    sub: Option<Vec<String>>,
) -> PyResult<Bound<'py, PyDict>> {
    let est_method = if estimation {
        gsem_sem::EstimationMethod::Dwls
    } else {
        gsem_sem::EstimationMethod::Ml
    };
    let model_opt = if model.is_empty() { None } else { Some(model) };
    let ldsc_result = pyany_to_ldsc_result(covstruc)?;

    let result = if let Some(ref sub_names) = sub {
        let k = ldsc_result.s.nrows();
        let all_names: Vec<String> = (0..k).map(|i| format!("V{}", i + 1)).collect();
        let sub_indices: Vec<usize> = sub_names
            .iter()
            .filter_map(|name| all_names.iter().position(|n| n == name))
            .collect();
        gsem_sem::rgmodel::run_rgmodel_sub(
            &ldsc_result.s,
            &ldsc_result.v,
            est_method,
            model_opt,
            std_lv,
            &sub_indices,
        )
    } else {
        gsem_sem::rgmodel::run_rgmodel_with_model(
            &ldsc_result.s,
            &ldsc_result.v,
            est_method,
            model_opt,
            std_lv,
        )
    }
    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;

    let out = PyDict::new(py);
    out.set_item("R", mat_to_pyarray(py, &result.r))?;
    out.set_item("V_R", mat_to_pyarray(py, &result.v_r))?;
    Ok(out)
}

/// Run HDL (High-Definition Likelihood) estimation.
///
/// Returns a `PyLdscResult` with the same shape as `ldsc()`'s output, so
/// downstream functions (`usermodel`, `commonfactor`, `rgmodel`, ...) can
/// consume it directly.
///
/// Examples
/// --------
/// >>> import genomicsem as gs
/// >>> cov = gs.hdl(
/// ...     traits=["T1.sumstats.gz", "T2.sumstats.gz"],
/// ...     trait_names=["T1", "T2"],
/// ...     ld_path="UKB_imputed_SVD_eigen99_extraction/",
/// ...     n_ref=335265.0,
/// ...     method="piecewise",
/// ... )
/// >>> cov.s  # 2x2 genetic covariance
#[pyfunction]
#[pyo3(signature = (traits, sample_prev=None, population_prev=None, trait_names=None, ld_path="", n_ref=335265.0, method="piecewise"))]
#[allow(unused_variables)]
fn hdl(
    py: Python<'_>,
    traits: Vec<String>,
    sample_prev: Option<Vec<Option<f64>>>,
    population_prev: Option<Vec<Option<f64>>>,
    trait_names: Option<Vec<String>>,
    ld_path: &str,
    n_ref: f64,
    method: &str,
) -> PyResult<PyLdscResult> {
    use gsem_ldsc::hdl::{HdlConfig, HdlMethod, HdlTraitData};
    if trait_names.is_some() {
        log::info!("trait_names are used for labeling output in the Python wrapper");
    }

    let sp: Vec<Option<f64>> = sample_prev.unwrap_or_default();
    let pp: Vec<Option<f64>> = population_prev.unwrap_or_default();

    let trait_data: Vec<HdlTraitData> = gsem::io::gwas_reader::load_trait_data(&traits)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{e}")))?
        .into_iter()
        .map(HdlTraitData::from)
        .collect();

    let hdl_method = match method.to_lowercase().as_str() {
        "jackknife" => HdlMethod::Jackknife,
        _ => HdlMethod::Piecewise,
    };

    let config = HdlConfig {
        method: hdl_method,
        n_ref,
    };

    let ld_dir = std::path::Path::new(ld_path);
    let ld_pieces = gsem::io::hdl_reader::load_hdl_pieces(ld_dir)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{e}")))?;

    // Defensive: release the GIL while HDL runs (may be parallelized).
    let result = py
        .detach(|| gsem_ldsc::hdl::hdl(&trait_data, &sp, &pp, &ld_pieces, &config))
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;
    Ok(ldsc_result_to_py(&result.to_ldsc_result(), false))
}

/// Run stratified LDSC (s-LDSC).
///
/// Partitions heritability across functional annotations. Returns per-annotation
/// S / V matrices and the `tau` table.
///
/// Examples
/// --------
/// >>> import genomicsem as gs
/// >>> res = gs.s_ldsc(
/// ...     traits=["T1.sumstats.gz"],
/// ...     ld="baselineLD_v2.2/",
/// ...     wld="weights_hm3_no_hla/",
/// ...     frq="1000G_Phase3_frq/",
/// ...     trait_names=["T1"],
/// ... )
/// >>> res["tau"]
#[pyfunction]
#[pyo3(signature = (traits, sample_prev=None, population_prev=None, ld="", wld="", frq="", trait_names=None, n_blocks=200, ldsc_log=None, exclude_cont=true))]
#[allow(unused_variables)]
fn s_ldsc<'py>(
    py: Python<'py>,
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
) -> PyResult<Bound<'py, PyDict>> {
    if trait_names.is_some() {
        log::info!("trait_names are used for labeling output in the Python wrapper");
    }
    if let Some(ref path) = ldsc_log {
        log::info!("ldsc_log='{path}' — file logging is handled at the Python/R wrapper level");
    }

    let sp: Vec<Option<f64>> = sample_prev.unwrap_or_default();
    let pp: Vec<Option<f64>> = population_prev.unwrap_or_default();

    // Read trait files
    let trait_data = gsem::io::gwas_reader::load_trait_data(&traits)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("{e}")))?;

    let ld_path = std::path::Path::new(ld);
    let wld_path = std::path::Path::new(wld);
    let chromosomes: Vec<usize> = (1..=22).collect();

    // Read annotation LD scores
    let mut annot_data =
        gsem_ldsc::annot_reader::read_annot_ld_scores(ld_path, wld_path, &chromosomes)
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
        n_blocks,
        rm_flank: false,
        flank_kb: 500,
    };

    // Defensive: release the GIL while stratified LDSC runs (uses rayon).
    let result = py
        .detach(|| {
            gsem_ldsc::stratified::s_ldsc(
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
            )
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;

    let s_list = PyList::empty(py);
    for m in &result.s_annot {
        s_list.append(mat_to_pyarray(py, m))?;
    }
    let v_list = PyList::empty(py);
    for m in &result.v_annot {
        v_list.append(mat_to_pyarray(py, m))?;
    }

    let out = PyDict::new(py);
    out.set_item("annotations", result.annotations)?;
    out.set_item("m_annot", result.m_annot)?;
    out.set_item("m_total", result.m_total)?;
    out.set_item("prop", result.prop)?;
    out.set_item("I", mat_to_pyarray(py, &result.i_mat))?;
    out.set_item("S_annot", s_list)?;
    out.set_item("V_annot", v_list)?;
    Ok(out)
}

/// Enrichment analysis using stratified LDSC results.
///
/// `s_annot` / `v_annot` are lists of NumPy 2D arrays (one matrix per
/// annotation). Returns a columnar dict.
///
/// Examples
/// --------
/// >>> import genomicsem as gs
/// >>> sldsc = gs.s_ldsc(traits=["T1.sumstats.gz"], ld="baselineLD/", wld="weights/", frq="1000G/")
/// >>> enr = gs.enrich(
/// ...     s_annot=sldsc["s_annot"],
/// ...     v_annot=sldsc["v_annot"],
/// ...     model="F1 =~ NA*V1",
/// ... )
/// >>> enr["enrichment"]
#[pyfunction]
#[pyo3(signature = (s_baseline, s_annot, v_annot, annotation_names, m_annot, m_total))]
fn enrich<'py>(
    py: Python<'py>,
    s_baseline: PyReadonlyArray2<'py, f64>,
    s_annot: &Bound<'py, PyAny>,
    v_annot: &Bound<'py, PyAny>,
    annotation_names: Vec<String>,
    m_annot: Vec<f64>,
    m_total: f64,
) -> PyResult<Bound<'py, PyDict>> {
    let s_base_mat = pyarray_to_mat(&s_baseline);

    let extract_mats = |obj: &Bound<'_, PyAny>, label: &str| -> PyResult<Vec<faer::Mat<f64>>> {
        let mut mats = Vec::new();
        for (i, item) in obj.try_iter()?.enumerate() {
            let item = item?;
            let arr: PyReadonlyArray2<'_, f64> = item.extract().map_err(|e| {
                pyo3::exceptions::PyValueError::new_err(format!("{label}[{i}]: {e}"))
            })?;
            mats.push(pyarray_to_mat(&arr));
        }
        Ok(mats)
    };
    let s_annot_mats = extract_mats(s_annot, "s_annot")?;
    let v_annot_mats = extract_mats(v_annot, "v_annot")?;

    let result = gsem::stats::enrich::enrichment_test(
        &s_base_mat,
        &s_annot_mats,
        &v_annot_mats,
        &annotation_names,
        &m_annot,
        m_total,
    );

    let out = PyDict::new(py);
    out.set_item("annotation", result.annotations)?;
    out.set_item("enrichment", result.enrichment)?;
    out.set_item("se", result.se)?;
    out.set_item("p", result.p)?;
    Ok(out)
}

/// Simulate GWAS summary statistics.
///
/// Returns a NumPy 2D array of simulated Z-scores with shape `k × n_snps`.
///
/// Examples
/// --------
/// >>> import numpy as np, genomicsem as gs
/// >>> covmat = np.array([[0.60, 0.42],
/// ...                    [0.42, 0.50]])
/// >>> z = gs.sim_ldsc(
/// ...     covmat=covmat,
/// ...     n=[1e5, 1e5],
/// ...     ld="eur_w_ld_chr/",
/// ... )
/// >>> z.shape  # (2, n_snps)
#[pyfunction]
#[pyo3(signature = (s_matrix, n_per_trait, ld_scores, m, intercepts=None, r_pheno=None, n_overlap=0.0))]
fn sim_ldsc<'py>(
    py: Python<'py>,
    s_matrix: PyReadonlyArray2<'py, f64>,
    n_per_trait: Vec<f64>,
    ld_scores: Vec<f64>,
    m: f64,
    intercepts: Option<PyReadonlyArray2<'py, f64>>,
    r_pheno: Option<PyReadonlyArray2<'py, f64>>,
    n_overlap: f64,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let s_mat = pyarray_to_mat(&s_matrix);

    let config = gsem::stats::simulation::SimConfig {
        intercepts: intercepts.as_ref().map(pyarray_to_mat),
        r_pheno: r_pheno.as_ref().map(pyarray_to_mat),
        n_overlap,
    };

    let result =
        gsem::stats::simulation::simulate_sumstats(&s_mat, &n_per_trait, &ld_scores, m, &config);
    // result is Vec<Vec<f64>> with shape (k × n_snps).
    let k = result.len();
    let n_snps = if k == 0 { 0 } else { result[0].len() };
    let packed = faer::Mat::from_fn(k, n_snps, |i, j| result[i][j]);
    Ok(mat_to_pyarray(py, &packed))
}

/// Run multi-SNP analysis.
///
/// Fits a SEM that regresses the common factor on multiple SNPs
/// simultaneously, returning a single parameter table.
///
/// Examples
/// --------
/// >>> import genomicsem as gs
/// >>> res = gs.multi_snp(
/// ...     covstruc=covstruc,
/// ...     snps="merged_sumstats.tsv",
/// ...     model="F1 =~ NA*V1 + V2 + V3\nF1 ~ rs1 + rs2\nF1 ~~ 1*F1",
/// ...     snp_list=["rs1", "rs2"],
/// ... )
/// >>> res["parameters"]
#[pyfunction]
#[pyo3(signature = (covstruc, model, beta, se, var_snp, ld_matrix, snp_names, estimation="DWLS", max_iter=1000))]
fn multi_snp<'py>(
    py: Python<'py>,
    covstruc: &Bound<'py, PyAny>,
    model: &str,
    beta: PyReadonlyArray2<'py, f64>,
    se: PyReadonlyArray2<'py, f64>,
    var_snp: Vec<f64>,
    ld_matrix: PyReadonlyArray2<'py, f64>,
    snp_names: Vec<String>,
    estimation: &str,
    max_iter: usize,
) -> PyResult<Bound<'py, PyDict>> {
    let ldsc_result = pyany_to_ldsc_result(covstruc)?;
    let ld_mat = pyarray_to_mat(&ld_matrix);

    // beta / se are passed as n_snps × k 2D arrays; unpack into the
    // Vec<Vec<f64>> layout the engine expects.
    let pyarray_to_rows = |arr: &PyReadonlyArray2<'_, f64>| -> Vec<Vec<f64>> {
        let view = arr.as_array();
        let (n, k) = (view.nrows(), view.ncols());
        (0..n)
            .map(|i| (0..k).map(|j| view[[i, j]]).collect())
            .collect()
    };
    let beta_rows = pyarray_to_rows(&beta);
    let se_rows = pyarray_to_rows(&se);

    let pt = gsem_sem::syntax::parse_model(model, false)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("model parse error: {e}")))?;
    let config = gsem::gwas::multi_snp::MultiSnpConfig {
        model: pt,
        estimation: gsem_sem::EstimationMethod::from_str_lossy(estimation),
        max_iter,
        snp_var_se: None,
    };

    // `run_multi_snp` now takes borrowed rows.
    let beta_refs: Vec<&[f64]> = beta_rows.iter().map(Vec::as_slice).collect();
    let se_refs: Vec<&[f64]> = se_rows.iter().map(Vec::as_slice).collect();
    let result = gsem::gwas::multi_snp::run_multi_snp(
        &config,
        &ldsc_result.s,
        &ldsc_result.v,
        &beta_refs,
        &se_refs,
        &var_snp,
        &ld_mat,
        &snp_names,
    );

    let out = PyDict::new(py);
    out.set_item("converged", result.converged)?;
    out.set_item("chisq", result.chisq)?;
    out.set_item("df", result.chisq_df as i64)?;
    out.set_item("params", param_results_to_dict(py, &result.params)?)?;
    Ok(out)
}

/// Run multi-gene analysis (reuses multi-SNP engine).
///
/// TWAS variant of `multi_snp` — feeds gene-level summary stats (FUSION
/// panels) through the same per-locus SEM fitting engine.
///
/// Examples
/// --------
/// >>> import genomicsem as gs
/// >>> res = gs.multi_gene(
/// ...     covstruc=covstruc,
/// ...     genes="gene_level.tsv",
/// ...     model="F1 =~ NA*V1 + V2 + V3\nF1 ~ Gene\nF1 ~~ 1*F1",
/// ...     gene_list=["ENSG1", "ENSG2"],
/// ... )
/// >>> res["parameters"]
#[pyfunction]
#[pyo3(signature = (covstruc, model, beta, se, var_gene, ld_matrix, gene_names, estimation="DWLS", max_iter=1000))]
fn multi_gene<'py>(
    py: Python<'py>,
    covstruc: &Bound<'py, PyAny>,
    model: &str,
    beta: PyReadonlyArray2<'py, f64>,
    se: PyReadonlyArray2<'py, f64>,
    var_gene: Vec<f64>,
    ld_matrix: PyReadonlyArray2<'py, f64>,
    gene_names: Vec<String>,
    estimation: &str,
    max_iter: usize,
) -> PyResult<Bound<'py, PyDict>> {
    multi_snp(
        py, covstruc, model, beta, se, var_gene, ld_matrix, gene_names, estimation, max_iter,
    )
}

/// Run Generalized Least Squares regression.
///
/// Returns a columnar dict `{beta, se, z, p}`.
///
/// Examples
/// --------
/// >>> import numpy as np, genomicsem as gs
/// >>> y = np.array([0.50, 0.30, 0.20, 0.10])
/// >>> v_y = np.eye(4) * 0.01
/// >>> x = np.array([[1, 0],
/// ...               [2, 1],
/// ...               [3, 0],
/// ...               [4, 1]], dtype=float)
/// >>> fit = gs.summary_gls(y=y, v_y=v_y, predictors=x, intercept=True)
/// >>> fit["beta"], fit["p"]
#[pyfunction]
#[pyo3(signature = (x, y, v, intercept=true))]
fn summary_gls<'py>(
    py: Python<'py>,
    x: PyReadonlyArray2<'py, f64>,
    y: Vec<f64>,
    v: PyReadonlyArray2<'py, f64>,
    intercept: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let mut x_mat = pyarray_to_mat(&x);
    let v_mat = pyarray_to_mat(&v);

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
            let out = PyDict::new(py);
            out.set_item("beta", result.beta)?;
            out.set_item("se", result.se)?;
            out.set_item("z", result.z)?;
            out.set_item("p", result.p)?;
            Ok(out)
        }
        None => Err(pyo3::exceptions::PyRuntimeError::new_err(
            "GLS failed (singular matrix?)",
        )),
    }
}

/// Python module definition.
#[pymodule]
fn genomicsem(m: &Bound<'_, PyModule>) -> PyResult<()> {
    pyo3_log::init();
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
