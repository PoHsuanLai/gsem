//! Zero-copy conversions between extendr native types and the internal gsem
//! matrix/struct types used by the R binding.
//!
//! The R binding used to route every matrix and result struct across the
//! boundary as a JSON string; this module replaces that with direct
//! `RMatrix<f64>` / `List` construction so no serialization happens.

use extendr_api::prelude::*;
use faer::Mat;
use gsem::gwas::user_gwas::{SnpParamResult, SnpResult};
use gsem::io::sumstats_reader::MergedSnp;
use gsem::io::twas_reader::TwasGene;
use gsem_ldsc::LdscResult;

// ---------------------------------------------------------------------------
// Matrix <-> RMatrix<f64>
// ---------------------------------------------------------------------------

/// Convert a faer `Mat<f64>` into an extendr column-major `RMatrix<f64>`.
pub fn mat_to_rmatrix(m: &Mat<f64>) -> RMatrix<f64> {
    RMatrix::new_matrix(m.nrows(), m.ncols(), |i, j| m[(i, j)])
}

/// Convert an `RMatrix<f64>` (column-major) into a faer `Mat<f64>`.
pub fn rmatrix_to_mat(r: &RMatrix<f64>) -> Mat<f64> {
    let nrows = r.nrows();
    let ncols = r.ncols();
    let data = r.data();
    Mat::from_fn(nrows, ncols, |i, j| data[i + j * nrows])
}

/// Try to convert an arbitrary `Robj` into a faer `Mat<f64>`. Accepts any
/// numeric matrix the `RMatrix<f64>::try_from` path accepts.
pub fn robj_to_mat(obj: &Robj) -> Result<Mat<f64>> {
    let m: RMatrix<f64> = obj.try_into()?;
    Ok(rmatrix_to_mat(&m))
}

// ---------------------------------------------------------------------------
// LdscResult <-> List
// ---------------------------------------------------------------------------

/// Build a named R list from an `LdscResult`.
pub fn ldsc_result_to_list(result: &LdscResult) -> List {
    list!(
        s = mat_to_rmatrix(&result.s),
        v = mat_to_rmatrix(&result.v),
        i_mat = mat_to_rmatrix(&result.i_mat),
        n_vec = result.n_vec.clone(),
        m = result.m
    )
}

/// Build a named R list from an `LdscResult`, adding the standardized
/// `s_stand` / `v_stand` matrices (same computation as the old
/// `ldsc_result_to_json_stand` path).
pub fn ldsc_result_to_list_stand(result: &LdscResult) -> List {
    let s_stand = gsem_matrix::smooth::cov_to_cor(&result.s);

    let k = result.s.nrows();
    let kstar = k * (k + 1) / 2;
    let s_vec = gsem_matrix::vech::vech(&result.s).expect("S must be square");
    let s_stand_vec = gsem_matrix::vech::vech(&s_stand).expect("S_Stand must be square");
    let scale: Vec<f64> = s_stand_vec
        .iter()
        .zip(s_vec.iter())
        .enumerate()
        .map(|(i, (&st, &orig))| {
            let ratio = if orig.abs() > 1e-30 { st / orig } else { 0.0 };
            result.v[(i, i)].sqrt() * ratio
        })
        .collect();
    let v_cor = gsem_matrix::smooth::cov_to_cor(&result.v);
    let v_stand = Mat::from_fn(kstar, kstar, |i, j| scale[i] * v_cor[(i, j)] * scale[j]);

    list!(
        s = mat_to_rmatrix(&result.s),
        v = mat_to_rmatrix(&result.v),
        i_mat = mat_to_rmatrix(&result.i_mat),
        n_vec = result.n_vec.clone(),
        m = result.m,
        s_stand = mat_to_rmatrix(&s_stand),
        v_stand = mat_to_rmatrix(&v_stand)
    )
}

/// Look up a field in a named list, accepting either of two possible names.
fn list_field(list: &List, primary: &str, alt: &str) -> Option<Robj> {
    for (name, obj) in list.iter() {
        if name == primary || name == alt {
            return Some(obj);
        }
    }
    None
}

/// Parse an R list with fields `s`, `v`, `i_mat`, `n_vec`, `m` into an
/// `LdscResult`. Accepts the capitalized variants (`S`, `V`, `I`, `N`, `m`)
/// as a convenience, so R callers can pass `ldsc()` output directly.
pub fn list_to_ldsc_result(list: &List) -> Result<LdscResult> {
    let s = list_field(list, "s", "S")
        .ok_or_else(|| Error::Other("covstruc list missing field 's' (or 'S')".into()))?;
    let v = list_field(list, "v", "V")
        .ok_or_else(|| Error::Other("covstruc list missing field 'v' (or 'V')".into()))?;
    let i = list_field(list, "i_mat", "I")
        .ok_or_else(|| Error::Other("covstruc list missing field 'i_mat' (or 'I')".into()))?;

    let s_mat = robj_to_mat(&s).map_err(|e| Error::Other(format!("covstruc$s: {e}")))?;
    let v_mat = robj_to_mat(&v).map_err(|e| Error::Other(format!("covstruc$v: {e}")))?;
    let i_mat = robj_to_mat(&i).map_err(|e| Error::Other(format!("covstruc$i_mat: {e}")))?;

    let n_vec: Vec<f64> = match list_field(list, "n_vec", "N") {
        Some(n) if !n.is_null() => {
            let doubles: Doubles = n
                .try_into()
                .map_err(|e| Error::Other(format!("covstruc$n_vec: {e}")))?;
            doubles.iter().map(|x| x.inner()).collect()
        }
        _ => Vec::new(),
    };

    let m: f64 = match list_field(list, "m", "M") {
        Some(m_obj) if !m_obj.is_null() => m_obj
            .as_real()
            .or_else(|| m_obj.as_integer().map(|v| v as f64))
            .ok_or_else(|| Error::Other("covstruc$m must be numeric".into()))?,
        _ => 0.0,
    };

    Ok(LdscResult {
        s: s_mat,
        v: v_mat,
        i_mat,
        n_vec,
        m,
    })
}

// ---------------------------------------------------------------------------
// Per-SNP / per-gene result tables as columnar lists
// ---------------------------------------------------------------------------

/// Temporary buffer for building a columnar parameter table.
struct ParamColumns {
    snp_idx: Vec<i32>, // 1-based R index of the containing SNP/gene row
    lhs: Vec<String>,
    op: Vec<String>,
    rhs: Vec<String>,
    est: Vec<f64>,
    se: Vec<f64>,
    z: Vec<f64>,
    p: Vec<f64>,
}

impl ParamColumns {
    fn with_capacity(cap: usize) -> Self {
        Self {
            snp_idx: Vec::with_capacity(cap),
            lhs: Vec::with_capacity(cap),
            op: Vec::with_capacity(cap),
            rhs: Vec::with_capacity(cap),
            est: Vec::with_capacity(cap),
            se: Vec::with_capacity(cap),
            z: Vec::with_capacity(cap),
            p: Vec::with_capacity(cap),
        }
    }

    fn push(&mut self, idx_1based: i32, p: &SnpParamResult) {
        self.snp_idx.push(idx_1based);
        self.lhs.push(p.lhs.clone());
        self.op.push(p.op.to_string());
        self.rhs.push(p.rhs.clone());
        self.est.push(p.est);
        self.se.push(p.se);
        self.z.push(p.z_stat);
        self.p.push(p.p_value);
    }

    fn into_list(self) -> List {
        list!(
            snp_idx = self.snp_idx,
            lhs = self.lhs,
            op = self.op,
            rhs = self.rhs,
            est = self.est,
            se = self.se,
            z = self.z,
            p = self.p
        )
    }
}

/// If any result reports Q_SNP statistics, return three parallel vectors (one
/// per result row) with NA for rows where the value is missing. Returns
/// `None` when no row has any Q_SNP information.
fn q_snp_columns(results: &[SnpResult]) -> Option<(Vec<Rfloat>, Vec<Rint>, Vec<Rfloat>)> {
    let any = results
        .iter()
        .any(|r| r.q_snp.is_some() || r.q_snp_df.is_some() || r.q_snp_p.is_some());
    if !any {
        return None;
    }
    let chisq: Vec<Rfloat> = results
        .iter()
        .map(|r| r.q_snp.map(Rfloat::from).unwrap_or(Rfloat::na()))
        .collect();
    let df: Vec<Rint> = results
        .iter()
        .map(|r| {
            r.q_snp_df
                .map(|d| Rint::from(d as i32))
                .unwrap_or(Rint::na())
        })
        .collect();
    let pval: Vec<Rfloat> = results
        .iter()
        .map(|r| r.q_snp_p.map(Rfloat::from).unwrap_or(Rfloat::na()))
        .collect();
    Some((chisq, df, pval))
}

/// Return value for per-SNP GWAS functions, structured as a named list of
/// equal-length column vectors. `params` is a nested named list with the
/// per-parameter rows (SNP, lhs/op/rhs, est/se/z/p); the row's parent SNP is
/// identified by the 1-based `snp_idx` column.
pub fn snp_results_to_list(results: &[SnpResult], snps: &[MergedSnp]) -> List {
    let n = results.len();
    let mut snp_name: Vec<String> = Vec::with_capacity(n);
    let mut chisq: Vec<f64> = Vec::with_capacity(n);
    let mut df: Vec<i32> = Vec::with_capacity(n);
    let mut converged: Vec<bool> = Vec::with_capacity(n);

    let mut params = ParamColumns::with_capacity(n * 4);

    for (row_idx, r) in results.iter().enumerate() {
        snp_name.push(snps[r.snp_idx].snp.clone());
        chisq.push(r.chisq);
        df.push(r.chisq_df as i32);
        converged.push(r.converged);

        let r_idx = (row_idx + 1) as i32;
        for p in &r.params {
            params.push(r_idx, p);
        }
    }

    let params_list = params.into_list();

    if let Some((q_chisq, q_df, q_pval)) = q_snp_columns(results) {
        list!(
            SNP = snp_name,
            chisq = chisq,
            df = df,
            converged = converged,
            Q_chisq = q_chisq,
            Q_df = q_df,
            Q_pval = q_pval,
            params = params_list
        )
    } else {
        list!(
            SNP = snp_name,
            chisq = chisq,
            df = df,
            converged = converged,
            params = params_list
        )
    }
}

/// TWAS variant: same layout as `snp_results_to_list`, but the identifier
/// columns are `Gene`/`Panel`/`HSQ` rather than `SNP`.
pub fn twas_results_to_list(results: &[SnpResult], genes: &[TwasGene]) -> List {
    let n = results.len();
    let mut gene_name: Vec<String> = Vec::with_capacity(n);
    let mut panel: Vec<String> = Vec::with_capacity(n);
    let mut hsq: Vec<f64> = Vec::with_capacity(n);
    let mut chisq: Vec<f64> = Vec::with_capacity(n);
    let mut df: Vec<i32> = Vec::with_capacity(n);
    let mut converged: Vec<bool> = Vec::with_capacity(n);

    let mut params = ParamColumns::with_capacity(n * 4);

    for (row_idx, r) in results.iter().enumerate() {
        let g = &genes[r.snp_idx];
        gene_name.push(g.gene.clone());
        panel.push(g.panel.clone());
        hsq.push(g.hsq);
        chisq.push(r.chisq);
        df.push(r.chisq_df as i32);
        converged.push(r.converged);

        let r_idx = (row_idx + 1) as i32;
        for p in &r.params {
            params.push(r_idx, p);
        }
    }

    let params_list = params.into_list();

    if let Some((q_chisq, q_df, q_pval)) = q_snp_columns(results) {
        list!(
            Gene = gene_name,
            Panel = panel,
            HSQ = hsq,
            chisq = chisq,
            df = df,
            converged = converged,
            Q_chisq = q_chisq,
            Q_df = q_df,
            Q_pval = q_pval,
            params = params_list
        )
    } else {
        list!(
            Gene = gene_name,
            Panel = panel,
            HSQ = hsq,
            chisq = chisq,
            df = df,
            converged = converged,
            params = params_list
        )
    }
}

/// Build a columnar list from a flat `Vec<SnpParamResult>`. Used for the
/// parameter tables returned by `usermodel` / `multi_snp` (single set of
/// rows, no SNP grouping).
pub fn param_results_to_list(params: &[SnpParamResult]) -> List {
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
    list!(
        lhs = lhs,
        op = op,
        rhs = rhs,
        est = est,
        se = se,
        z = z,
        p = p
    )
}

/// Build a columnar list from `ParamEstimate`s — used for the SEM-fit
/// parameter tables returned by `run_commonfactor`.
pub fn param_estimates_to_list(params: &[gsem_sem::ParamEstimate]) -> List {
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
    list!(
        lhs = lhs,
        op = op,
        rhs = rhs,
        est = est,
        se = se,
        z = z,
        p = p
    )
}

/// Helper: build a single-element error list `list(error = "...")` for the
/// early-return paths in each `#[extendr]` function.
pub fn error_list(msg: impl Into<String>) -> List {
    list!(error = msg.into())
}
