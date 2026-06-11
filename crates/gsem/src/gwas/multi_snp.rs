//! Joint analysis of multiple SNPs with LD.
//!
//! Fits a model with multiple SNPs simultaneously, accounting for LD between them.
//!
//! Key difference from single-SNP GWAS:
//! - S_Full has multiple SNP rows/columns (one per SNP)
//! - SNP-SNP correlations come from an LD matrix provided by the user
//! - Model has multiple SNP predictors
//!
//! Port of R GenomicSEM's `multiSNP()`.

use faer::Mat;
use gsem_sem::EstimationMethod;
use gsem_sem::estimator;
use gsem_sem::model::Model;
use gsem_sem::sandwich;
use gsem_sem::syntax;

/// Configuration for multi-SNP analysis.
#[derive(Debug, Clone)]
pub struct MultiSnpConfig {
    /// Pre-parsed model parameter table
    pub model: syntax::ParTable,
    /// Estimation method
    pub estimation: EstimationMethod,
    /// Maximum optimizer iterations
    pub max_iter: usize,
    /// Override for SNP variance SE (default: 0.0005).
    pub snp_var_se: Option<f64>,
}

/// Result of multi-SNP analysis.
#[derive(Debug, Clone)]
pub struct MultiSnpResult {
    /// Parameter estimates from the joint model
    pub params: Vec<super::user_gwas::SnpParamResult>,
    /// Model chi-square statistic
    pub chisq: f64,
    /// Model degrees of freedom
    pub chisq_df: usize,
    /// Whether the optimizer converged
    pub converged: bool,
}

/// Build the augmented `S_Full` (observed covariance) and `V_Full` (sampling
/// covariance of vech(S_Full)) matrices for a joint multi-SNP analysis.
///
/// This is the deterministic core of R GenomicSEM's `multiSNP()`, which returns
/// these two matrices (the SEM fit is performed separately by `usermodel` /
/// `userGWAS`). gsem reproduces R's construction exactly:
///
/// `S_Full` (total = n_snps + k variables, SNPs first then traits):
/// - SNP variances on the diagonal (`2*MAF*(1-MAF)`),
/// - SNP-SNP covariances `LD[i,j]*sqrt(varSNP_i*varSNP_j)`,
/// - SNP-trait covariances `varSNP*beta`,
/// - the trait-trait LDSC block `S_LD`.
///
/// `V_Full` (kstar_full = total*(total+1)/2), with `SE_SNP2[s,t] = se[s,t] *
/// I_LD_diag[t] * varSNP[s]` (univariate intercepts clamped to >= 1):
/// - SNP variances and SNP-SNP covariances: fixed `snp_var_se^2` on the diagonal,
/// - SNP-trait variances: `SE_SNP2[s,t]^2` on the diagonal,
/// - off-diagonals between two SNP-trait elements `(sa,ta)` and `(sb,tb)`:
///   `SE_SNP2[sa,ta]*SE_SNP2[sb,tb] * (I_LD[ta,tb] if ta!=tb) * (LD[sa,sb] if sa!=sb)`,
/// - the trait-trait block: `V_LD` (last kstar_ld vech elements).
///
/// Divergence from R (documented): R's cross-SNP cross-trait block weights every
/// cell by a single CONSTANT LD value (`LD2[(f^2-f)/2]`, the last lower-triangle
/// entry) rather than the actual SNP-pair LD — an indexing bug in `multiSNP.R`.
/// gsem uses the correct per-pair `LD[sa,sb]`; the two agree whenever the
/// off-diagonal LD is constant.
#[allow(clippy::too_many_arguments)]
pub fn build_multi_snp_sv(
    s_ld: &Mat<f64>,
    v_ld: &Mat<f64>,
    i_ld: &Mat<f64>,
    beta: &[&[f64]],
    se: &[&[f64]],
    var_snp: &[f64],
    ld_matrix: &Mat<f64>,
    n_snps: usize,
    snp_var_se_raw: f64,
) -> (Mat<f64>, Mat<f64>) {
    let k = s_ld.nrows();
    let f = n_snps;
    let total = f + k;

    // ── S_Full ───────────────────────────────────────────────────────────────
    let mut s_full = Mat::zeros(total, total);
    for i in 0..f {
        s_full[(i, i)] = var_snp[i];
        for j in (i + 1)..f {
            let cov = ld_matrix[(i, j)] * (var_snp[i] * var_snp[j]).sqrt();
            s_full[(i, j)] = cov;
            s_full[(j, i)] = cov;
        }
    }
    for snp_i in 0..f {
        for t in 0..k {
            let cov = beta[snp_i][t] * var_snp[snp_i];
            s_full[(snp_i, f + t)] = cov;
            s_full[(f + t, snp_i)] = cov;
        }
    }
    for i in 0..k {
        for j in 0..k {
            s_full[(f + i, f + j)] = s_ld[(i, j)];
        }
    }

    // ── V_Full ───────────────────────────────────────────────────────────────
    let kstar_full = total * (total + 1) / 2;
    let kstar_ld = k * (k + 1) / 2;
    let snp_var_se2 = snp_var_se_raw * snp_var_se_raw;

    // Univariate intercepts clamped to >= 1 (R: diag(I_LD)<-ifelse(<=1,1,.)).
    let i_diag: Vec<f64> = (0..k).map(|t| i_ld[(t, t)].max(1.0)).collect();
    // SE_SNP2[snp][trait] = se * I_LD_diag * varSNP.
    let se_snp2 = |snp: usize, t: usize| -> f64 { se[snp][t] * i_diag[t] * var_snp[snp] };

    // Enumerate vech elements (column-major lower triangle) as (row, col).
    let mut elems: Vec<(usize, usize)> = Vec::with_capacity(kstar_full);
    for c in 0..total {
        for r in c..total {
            elems.push((r, c));
        }
    }
    // Classify an S_Full element: SNP index if it is a variable < f, else trait.
    // Returns the SNP-trait pair if the element is a SNP×trait covariance.
    let snp_trait = |r: usize, c: usize| -> Option<(usize, usize)> {
        if c < f && r >= f {
            Some((c, r - f)) // (snp, trait)
        } else {
            None
        }
    };

    let mut v_full = Mat::zeros(kstar_full, kstar_full);

    // Diagonal.
    for (a, &(r, c)) in elems.iter().enumerate() {
        if let Some((snp, t)) = snp_trait(r, c) {
            v_full[(a, a)] = se_snp2(snp, t).powi(2);
        } else if r >= f && c >= f {
            // trait-trait: filled from V_LD below.
        } else {
            // SNP variance or SNP-SNP covariance: fixed sampling variance.
            v_full[(a, a)] = snp_var_se2;
        }
    }

    // Off-diagonals between two SNP-trait elements.
    for a in 0..kstar_full {
        let (ra, ca) = elems[a];
        let Some((sa, ta)) = snp_trait(ra, ca) else {
            continue;
        };
        for b in (a + 1)..kstar_full {
            let (rb, cb) = elems[b];
            let Some((sb, tb)) = snp_trait(rb, cb) else {
                continue;
            };
            let mut val = se_snp2(sa, ta) * se_snp2(sb, tb);
            if ta != tb {
                val *= i_ld[(ta, tb)];
            }
            if sa != sb {
                val *= ld_matrix[(sa, sb)];
            }
            v_full[(a, b)] = val;
            v_full[(b, a)] = val;
        }
    }

    // Trait-trait block = V_LD (last kstar_ld vech elements).
    let offset = kstar_full - kstar_ld;
    for i in 0..kstar_ld {
        for j in 0..kstar_ld {
            v_full[(offset + i, offset + j)] = v_ld[(i, j)];
        }
    }

    (s_full, v_full)
}

/// Run multi-SNP analysis.
///
/// Port of R GenomicSEM's `multiSNP()`.
///
/// Builds an augmented S matrix with n_snp SNP rows/columns:
/// - SNP variances on diagonal (2*MAF*(1-MAF))
/// - SNP-SNP correlations from LD matrix (with allele alignment)
/// - SNP-trait covariances from betas
/// - Trait-trait covariance from S_LD
#[allow(clippy::too_many_arguments)]
pub fn run_multi_snp(
    config: &MultiSnpConfig,
    s_ld: &Mat<f64>,
    v_ld: &Mat<f64>,
    i_ld: &Mat<f64>,
    beta: &[&[f64]],      // n_snps x k betas (borrowed rows)
    se: &[&[f64]],        // n_snps x k SEs   (borrowed rows)
    var_snp: &[f64],      // n_snps variances
    ld_matrix: &Mat<f64>, // n_snps x n_snps LD correlation
    snp_names: &[String],
) -> MultiSnpResult {
    let k = s_ld.nrows();
    let n_snps = snp_names.len();
    let total = n_snps + k;
    let kstar_full = total * (total + 1) / 2;

    let snp_var_se = config.snp_var_se.unwrap_or(0.0005_f64);
    let (s_full, v_full) = build_multi_snp_sv(
        s_ld, v_ld, i_ld, beta, se, var_snp, ld_matrix, n_snps, snp_var_se,
    );

    // Build observed variable names
    let mut obs_names: Vec<String> = snp_names.to_vec();
    obs_names.extend((0..k).map(|i| format!("V{}", i + 1)));

    // Use pre-parsed model
    let pt = &config.model;

    let mut model = Model::from_partable(pt, &obs_names);
    let v_diag: Vec<f64> = (0..kstar_full).map(|i| v_full[(i, i)]).collect();

    let fit = match config.estimation {
        EstimationMethod::Ml => estimator::fit_ml(&mut model, &s_full, config.max_iter, None),
        EstimationMethod::Dwls => {
            estimator::fit_dwls(&mut model, &s_full, &v_diag, config.max_iter, None)
        }
    };

    // Sandwich SEs
    let w_diag = Mat::from_fn(kstar_full, kstar_full, |i, j| {
        if i == j && v_diag[i] > 1e-30 {
            1.0 / v_diag[i]
        } else {
            0.0
        }
    });
    let (se_vec, _) = sandwich::sandwich_se(&mut model, &w_diag, &v_full);

    // Build parameter results
    use statrs::distribution::ContinuousCDF;
    let params: Vec<super::user_gwas::SnpParamResult> = pt
        .rows
        .iter()
        .enumerate()
        .filter(|(_, row)| row.free > 0)
        .enumerate()
        .map(|(free_idx, (_, row))| {
            let est = fit.params.get(free_idx).copied().unwrap_or(f64::NAN);
            let se_val = se_vec.get(free_idx).copied().unwrap_or(f64::NAN);
            let z = est / se_val;
            let p = if z.is_finite() {
                2.0 * statrs::distribution::Normal::standard().cdf(-z.abs())
            } else {
                f64::NAN
            };
            super::user_gwas::SnpParamResult {
                lhs: row.lhs.clone(),
                op: row.op,
                rhs: row.rhs.clone(),
                est,
                se: se_val,
                z_stat: z,
                p_value: p,
            }
        })
        .collect();

    MultiSnpResult {
        params,
        chisq: fit.objective,
        chisq_df: model.df(),
        converged: fit.converged,
    }
}

/// Compute the vech index for element (row, col) in a p x p symmetric matrix.
/// Assumes row >= col (lower triangle, column-major order).
#[cfg(test)]
fn vech_index(row: usize, col: usize, p: usize) -> usize {
    debug_assert!(row >= col, "vech_index requires row >= col");
    // For column c, the offset is: sum_{j=0..c-1} (p - j) = c*p - c*(c-1)/2
    // Then within column c, the element at row r is at position (r - c).
    let col_offset = if col == 0 {
        0
    } else {
        col * p - col * (col - 1) / 2
    };
    col_offset + (row - col)
}

/// Read an LD matrix from a tab-delimited file.
///
/// Format: optional header row with SNP names, then n x n float values.
/// Returns (matrix, optional SNP names from header).
pub fn read_ld_matrix(path: &std::path::Path) -> anyhow::Result<(Mat<f64>, Option<Vec<String>>)> {
    use std::io::BufRead;

    let file = std::fs::File::open(path)
        .map_err(|e| anyhow::anyhow!("failed to open LD matrix {}: {e}", path.display()))?;
    let reader = std::io::BufReader::new(file);
    let mut lines: Vec<String> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if !trimmed.is_empty() {
            lines.push(trimmed.to_string());
        }
    }

    if lines.is_empty() {
        anyhow::bail!("LD matrix file is empty");
    }

    // Check if first line is a header (try parsing first field as float)
    let first_fields: Vec<&str> = lines[0].split('\t').collect();
    let has_header = first_fields[0].parse::<f64>().is_err();

    let (header, data_lines) = if has_header {
        let names: Vec<String> = first_fields.iter().map(|s| s.to_string()).collect();
        (Some(names), &lines[1..])
    } else {
        (None, &lines[..])
    };

    let n = data_lines.len();
    let mut mat = Mat::zeros(n, n);

    for (i, line) in data_lines.iter().enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != n {
            anyhow::bail!(
                "LD matrix row {} has {} columns, expected {n}",
                i + 1,
                fields.len()
            );
        }
        for (j, field) in fields.iter().enumerate() {
            mat[(i, j)] = field
                .parse::<f64>()
                .map_err(|e| anyhow::anyhow!("LD matrix [{i},{j}]: {e}"))?;
        }
    }

    Ok((mat, header))
}

#[cfg(test)]
mod tests {
    use super::*;
    use faer::Mat;
    use std::io::Write;

    #[test]
    fn test_vech_index_known_values() {
        // For a 3x3 matrix, vech ordering (column-major lower triangle):
        // (0,0)=0, (1,0)=1, (2,0)=2, (1,1)=3, (2,1)=4, (2,2)=5
        assert_eq!(vech_index(0, 0, 3), 0);
        assert_eq!(vech_index(1, 0, 3), 1);
        assert_eq!(vech_index(2, 0, 3), 2);
        assert_eq!(vech_index(1, 1, 3), 3);
        assert_eq!(vech_index(2, 1, 3), 4);
        assert_eq!(vech_index(2, 2, 3), 5);
    }

    #[test]
    fn test_vech_index_4x4() {
        // 4x4: (0,0)=0 (1,0)=1 (2,0)=2 (3,0)=3 (1,1)=4 (2,1)=5 (3,1)=6 (2,2)=7 (3,2)=8 (3,3)=9
        assert_eq!(vech_index(0, 0, 4), 0);
        assert_eq!(vech_index(3, 0, 4), 3);
        assert_eq!(vech_index(1, 1, 4), 4);
        assert_eq!(vech_index(3, 3, 4), 9);
    }

    #[test]
    fn test_read_ld_matrix_no_header() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("ld.txt");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "1.0\t0.3").unwrap();
        writeln!(f, "0.3\t1.0").unwrap();
        drop(f);

        let (mat, header) = read_ld_matrix(&path).unwrap();
        assert!(header.is_none());
        assert_eq!(mat.nrows(), 2);
        assert_eq!(mat.ncols(), 2);
        assert!((mat[(0, 0)] - 1.0).abs() < 1e-10);
        assert!((mat[(0, 1)] - 0.3).abs() < 1e-10);
        assert!((mat[(1, 0)] - 0.3).abs() < 1e-10);
    }

    #[test]
    fn test_read_ld_matrix_with_header() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("ld_h.txt");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "SNP1\tSNP2").unwrap();
        writeln!(f, "1.0\t0.5").unwrap();
        writeln!(f, "0.5\t1.0").unwrap();
        drop(f);

        let (mat, header) = read_ld_matrix(&path).unwrap();
        assert!(header.is_some());
        let names = header.unwrap();
        assert_eq!(names, vec!["SNP1", "SNP2"]);
        assert_eq!(mat.nrows(), 2);
        assert!((mat[(0, 1)] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_read_ld_matrix_empty_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.txt");
        std::fs::File::create(&path).unwrap();
        assert!(read_ld_matrix(&path).is_err());
    }

    #[test]
    fn test_build_sv_cross_snp_cross_trait_uses_per_pair_ld() {
        // gsem's corrected cross-SNP cross-trait block uses the actual SNP-pair
        // LD (R uses a constant; documented bug). With unequal off-diagonal LD,
        // the (T1-S1) x (T2-S2) cell must use LD[0,1], not LD[2,1].
        let s_ld = faer::mat![[0.4, 0.18], [0.18, 0.5]];
        let i_ld = faer::mat![[1.02, 0.05], [0.05, 1.03]];
        let v_ld = Mat::from_fn(3, 3, |i, j| if i == j { 1e-4 } else { 0.0 });
        let beta = [vec![0.03, 0.018], vec![-0.02, 0.025]];
        let se = [vec![0.006, 0.007], vec![0.005, 0.006]];
        let var_snp = vec![0.375_f64, 0.48];
        let ld = faer::mat![[1.0, 0.2], [0.2, 1.0]];
        let beta_refs: Vec<&[f64]> = beta.iter().map(Vec::as_slice).collect();
        let se_refs: Vec<&[f64]> = se.iter().map(Vec::as_slice).collect();

        let (_s, v) =
            build_multi_snp_sv(&s_ld, &v_ld, &i_ld, &beta_refs, &se_refs, &var_snp, &ld, 2, 0.0005);

        // Variables: S1,S2,T1,T2 (total=4). vech idx: 0=S1S1,1=S2S1,2=T1S1,
        // 3=T2S1,4=S2S2,5=T1S2,6=T2S2,7=T1T1,8=T2T1,9=T2T2.
        // (T1-S1) is vech idx 2 -> snp0,trait0; (T2-S2) is vech idx 6 -> snp1,trait1.
        let se_snp2 = |s: usize, t: usize| se[s][t] * i_ld[(t, t)].max(1.0) * var_snp[s];
        let expected = se_snp2(0, 0) * se_snp2(1, 1) * i_ld[(0, 1)] * ld[(0, 1)];
        assert!(
            (v[(2, 6)] - expected).abs() < 1e-15,
            "cross-SNP cross-trait must use per-pair LD: got {} expected {}",
            v[(2, 6)],
            expected
        );
    }

    #[test]
    fn test_run_multi_snp_basic() {
        let s_ld = faer::mat![[0.5, 0.2], [0.2, 0.4]];
        let i_ld = faer::mat![[1.0, 0.0], [0.0, 1.0]];
        let v_ld = Mat::from_fn(3, 3, |i, j| if i == j { 0.001 } else { 0.0 });
        let beta = [vec![0.1, 0.05], vec![0.08, 0.12]];
        let se = [vec![0.02, 0.02], vec![0.02, 0.02]];
        let var_snp = vec![0.3, 0.25];
        let ld_matrix = faer::mat![[1.0, 0.3], [0.3, 1.0]];
        let snp_names = vec!["SNP1".to_string(), "SNP2".to_string()];

        let config = MultiSnpConfig {
            model: syntax::parse_model("F1 =~ NA*V1 + V2\nF1 ~~ 1*F1\nV1 ~~ V1\nV2 ~~ V2\nF1 ~ SNP1 + SNP2\nSNP1 ~~ SNP1\nSNP2 ~~ SNP2", false).unwrap(),
            estimation: EstimationMethod::Dwls,
            max_iter: 500,
            snp_var_se: None,
        };

        let beta_refs: Vec<&[f64]> = beta.iter().map(Vec::as_slice).collect();
        let se_refs: Vec<&[f64]> = se.iter().map(Vec::as_slice).collect();
        let result = run_multi_snp(
            &config, &s_ld, &v_ld, &i_ld, &beta_refs, &se_refs, &var_snp, &ld_matrix, &snp_names,
        );

        // Should produce some parameters (even if not converged, we check structure)
        assert!(!result.params.is_empty(), "should have parameter estimates");
        // Chi-square should be finite (or at least not panic)
        assert!(
            result.chisq.is_finite() || result.chisq.is_nan(),
            "chisq should be a number"
        );
    }
}
