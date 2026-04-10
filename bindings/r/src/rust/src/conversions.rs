//! Zero-copy conversions between extendr native types and the internal gsem
//! matrix/struct types used by the R binding.
//!
//! The R binding used to route every matrix and result struct across the
//! boundary as a JSON string; this module replaces that with direct
//! `RMatrix<f64>` / `List` construction so no serialization happens.

use extendr_api::prelude::*;
use faer::Mat;
use gsem::gwas::user_gwas::{SnpParamResult, SnpResult};
use gsem::io::sumstats_reader::MergedSumstats;
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
// Per-SNP / per-gene result tables as two-table normalized lists
// ---------------------------------------------------------------------------

use rayon::prelude::*;

/// Below this many total params, the serial flatten beats the rayon setup.
const PARAMS_PARALLEL_THRESHOLD: usize = 50_000;

/// Columnar buffers for the per-SNP table.
struct SnpCols {
    snp: Vec<String>,
    chr: Vec<Rint>,
    bp: Vec<f64>,
    maf: Vec<f64>,
    a1: Vec<String>,
    a2: Vec<String>,
    chisq: Vec<f64>,
    df: Vec<i32>,
    converged: Vec<bool>,
}

impl SnpCols {
    fn with_capacity(n: usize) -> Self {
        Self {
            snp: Vec::with_capacity(n),
            chr: Vec::with_capacity(n),
            bp: Vec::with_capacity(n),
            maf: Vec::with_capacity(n),
            a1: Vec::with_capacity(n),
            a2: Vec::with_capacity(n),
            chisq: Vec::with_capacity(n),
            df: Vec::with_capacity(n),
            converged: Vec::with_capacity(n),
        }
    }
}

/// Columnar buffers for the per-gene (TWAS) table.
struct GeneCols {
    gene: Vec<String>,
    panel: Vec<String>,
    hsq: Vec<f64>,
    chisq: Vec<f64>,
    df: Vec<i32>,
    converged: Vec<bool>,
}

impl GeneCols {
    fn with_capacity(n: usize) -> Self {
        Self {
            gene: Vec::with_capacity(n),
            panel: Vec::with_capacity(n),
            hsq: Vec::with_capacity(n),
            chisq: Vec::with_capacity(n),
            df: Vec::with_capacity(n),
            converged: Vec::with_capacity(n),
        }
    }
}

/// Flat params table, sized to `total_params` so rayon can fill
/// disjoint windows.
struct ParamCols {
    id: Vec<String>,
    lhs: Vec<String>,
    op: Vec<String>,
    rhs: Vec<String>,
    est: Vec<f64>,
    se: Vec<f64>,
    z: Vec<f64>,
    p: Vec<f64>,
}

impl ParamCols {
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
fn params_offsets(results: &[SnpResult]) -> Vec<usize> {
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
fn fill_params_window(
    snp_name: &str,
    r: &SnpResult,
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
fn split_param_windows<'a>(
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
fn fill_params_table(
    results: &[SnpResult],
    ids_per_snp: &[String],
    offsets: &[usize],
    pc: &mut ParamCols,
) {
    let total = *offsets.last().unwrap_or(&0);
    if total < PARAMS_PARALLEL_THRESHOLD {
        for (i, r) in results.iter().enumerate() {
            let start = offsets[i];
            let end = offsets[i + 1];
            fill_params_window(
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

    let windows = split_param_windows(
        offsets, &mut pc.id, &mut pc.lhs, &mut pc.op, &mut pc.rhs, &mut pc.est, &mut pc.se,
        &mut pc.z, &mut pc.p,
    );

    windows
        .into_par_iter()
        .zip(results.par_iter())
        .zip(ids_per_snp.par_iter())
        .for_each(
            |(((id_c, lhs_c, op_c, rhs_c, est_c, se_c, z_c, p_c), r), snp_name)| {
                fill_params_window(
                    snp_name, r, id_c, lhs_c, op_c, rhs_c, est_c, se_c, z_c, p_c,
                );
            },
        );
}

/// Extract `Q_SNP` columns into three parallel vectors, NA-filled where
/// missing. Returns `None` when no result carries any Q_SNP information.
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

/// Build the inner `snps` list, with or without `Q_SNP` columns.
fn build_snps_list(sc: SnpCols, q: Option<(Vec<Rfloat>, Vec<Rint>, Vec<Rfloat>)>) -> List {
    match q {
        Some((qc, qd, qp)) => list!(
            SNP = sc.snp,
            CHR = sc.chr,
            BP = sc.bp,
            MAF = sc.maf,
            A1 = sc.a1,
            A2 = sc.a2,
            chisq = sc.chisq,
            df = sc.df,
            converged = sc.converged,
            Q_chisq = qc,
            Q_df = qd,
            Q_pval = qp
        ),
        None => list!(
            SNP = sc.snp,
            CHR = sc.chr,
            BP = sc.bp,
            MAF = sc.maf,
            A1 = sc.a1,
            A2 = sc.a2,
            chisq = sc.chisq,
            df = sc.df,
            converged = sc.converged
        ),
    }
}

/// Build the inner `genes` list (TWAS), with or without `Q_SNP` columns.
fn build_genes_list(gc: GeneCols, q: Option<(Vec<Rfloat>, Vec<Rint>, Vec<Rfloat>)>) -> List {
    match q {
        Some((qc, qd, qp)) => list!(
            Gene = gc.gene,
            Panel = gc.panel,
            HSQ = gc.hsq,
            chisq = gc.chisq,
            df = gc.df,
            converged = gc.converged,
            Q_chisq = qc,
            Q_df = qd,
            Q_pval = qp
        ),
        None => list!(
            Gene = gc.gene,
            Panel = gc.panel,
            HSQ = gc.hsq,
            chisq = gc.chisq,
            df = gc.df,
            converged = gc.converged
        ),
    }
}

/// Build the inner `params` list. The id column is `SNP` or `Gene`
/// depending on `id_label`; extendr's `list!` macro binds names at
/// compile time, so the label is dispatched here rather than passed
/// through.
fn build_params_list(pc: ParamCols, id_label: &str) -> List {
    if id_label == "Gene" {
        list!(
            Gene = pc.id,
            lhs = pc.lhs,
            op = pc.op,
            rhs = pc.rhs,
            est = pc.est,
            se = pc.se,
            z = pc.z,
            p = pc.p
        )
    } else {
        list!(
            SNP = pc.id,
            lhs = pc.lhs,
            op = pc.op,
            rhs = pc.rhs,
            est = pc.est,
            se = pc.se,
            z = pc.z,
            p = pc.p
        )
    }
}

/// Build the two-element `list(snps = ..., params = ...)` return value
/// for `userGWAS` / `commonfactorGWAS`. Both inner slots are columnar
/// and get wrapped in `data.frame` by the R shim via its
/// `structure(..., class = "data.frame")` fast path.
///
/// `index_map` lets callers that filter the merged file (e.g. MAF=0
/// drop in `common_factor_gwas_rust`) pass the original
/// `MergedSumstats` plus a `Vec<usize>` of surviving rows rather than
/// a cloned subset — saves ~1.3 GB transient on the full bench.
/// Pass `None` for the identity mapping.
pub fn snp_results_to_list(
    results: &[SnpResult],
    merged: &MergedSumstats,
    index_map: Option<&[usize]>,
) -> List {
    let n = results.len();

    let mut sc = SnpCols::with_capacity(n);
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
                .map(|c| Rint::from(c as i32))
                .unwrap_or(Rint::na()),
        );
        sc.bp.push(
            merged
                .bp
                .as_ref()
                .and_then(|b| b.get(idx).copied())
                .map(|b| b as f64)
                .unwrap_or(f64::NAN),
        );
        sc.maf.push(merged.maf[idx]);
        sc.a1.push(merged.a1_string(idx));
        sc.a2.push(merged.a2_string(idx));
        sc.chisq.push(r.chisq);
        sc.df.push(r.chisq_df as i32);
        sc.converged.push(r.converged);
    }

    // `fill_params_table` takes `&[String]`, so reuse `sc.snp`
    // directly instead of cloning into a separate `ids_per_snp`
    // buffer — saves another ~80 MB transient at 4.94M SNPs.
    let offsets = params_offsets(results);
    let mut pc = ParamCols::zeroed(*offsets.last().unwrap_or(&0));
    fill_params_table(results, &sc.snp, &offsets, &mut pc);

    let snps_list = build_snps_list(sc, q_snp_columns(results));
    let params_list = build_params_list(pc, "SNP");
    list!(snps = snps_list, params = params_list)
}

/// TWAS variant of [`snp_results_to_list`]. Returns
/// `list(genes = ..., params = ...)` with `Gene` as the identifier
/// column in both inner tables.
pub fn twas_results_to_list(results: &[SnpResult], genes: &[TwasGene]) -> List {
    let n = results.len();

    let mut gc = GeneCols::with_capacity(n);
    let mut ids_per_snp: Vec<String> = Vec::with_capacity(n);
    for r in results {
        let g = &genes[r.snp_idx];
        gc.gene.push(g.gene.clone());
        gc.panel.push(g.panel.clone());
        gc.hsq.push(g.hsq);
        gc.chisq.push(r.chisq);
        gc.df.push(r.chisq_df as i32);
        gc.converged.push(r.converged);
        ids_per_snp.push(g.gene.clone());
    }

    let offsets = params_offsets(results);
    let mut pc = ParamCols::zeroed(*offsets.last().unwrap_or(&0));
    fill_params_table(results, &ids_per_snp, &offsets, &mut pc);

    let genes_list = build_genes_list(gc, q_snp_columns(results));
    let params_list = build_params_list(pc, "Gene");

    list!(genes = genes_list, params = params_list)
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
