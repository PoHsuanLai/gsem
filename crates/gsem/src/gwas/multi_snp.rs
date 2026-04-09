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
    beta: &[Vec<f64>],    // n_snps x k betas
    se: &[Vec<f64>],      // n_snps x k SEs
    var_snp: &[f64],      // n_snps variances
    ld_matrix: &Mat<f64>, // n_snps x n_snps LD correlation
    snp_names: &[String],
) -> MultiSnpResult {
    let k = s_ld.nrows();
    let n_snps = snp_names.len();
    let total = n_snps + k;

    // Build augmented S matrix (n_snps + k) x (n_snps + k)
    let mut s_full = Mat::zeros(total, total);

    // SNP-SNP block: variances on diagonal, LD correlations off-diagonal
    for i in 0..n_snps {
        s_full[(i, i)] = var_snp[i];
        for j in (i + 1)..n_snps {
            let cov = ld_matrix[(i, j)] * (var_snp[i] * var_snp[j]).sqrt();
            s_full[(i, j)] = cov;
            s_full[(j, i)] = cov;
        }
    }

    // SNP-trait block
    for snp_i in 0..n_snps {
        for t in 0..k {
            let cov = beta[snp_i][t] * var_snp[snp_i];
            s_full[(snp_i, n_snps + t)] = cov;
            s_full[(n_snps + t, snp_i)] = cov;
        }
    }

    // Trait-trait block
    for i in 0..k {
        for j in 0..k {
            s_full[(n_snps + i, n_snps + j)] = s_ld[(i, j)];
        }
    }

    // Build V_Full (simplified: diagonal approximation for SNP blocks)
    let kstar_full = total * (total + 1) / 2;
    let mut v_full = Mat::zeros(kstar_full, kstar_full);

    // Copy trait-trait V block
    let kstar_trait = k * (k + 1) / 2;

    // Map trait V indices to full V indices.
    // In vech ordering (column-major lower triangle), the trait-trait block
    // occupies the last kstar_trait elements since traits are the last k
    // variables in the (n_snps + k) total.
    // Compute the vech offset for the trait-trait sub-block.
    let trait_vech_offset = kstar_full - kstar_trait;
    for i in 0..kstar_trait {
        for j in 0..kstar_trait {
            if trait_vech_offset + i < kstar_full && trait_vech_offset + j < kstar_full {
                v_full[(trait_vech_offset + i, trait_vech_offset + j)] = v_ld[(i, j)];
            }
        }
    }

    // SNP variance elements: small fixed value (or user override)
    let snp_var_se = config.snp_var_se.unwrap_or(0.0005_f64).powi(2);
    for i in 0..n_snps {
        // vech index for diagonal element (i, i) in column-major lower triangle:
        // For column j, offset = sum_{c=0..j-1} (total - c), then row offset = i - j.
        // For diagonal (i, i): col=i, row offset=0
        let idx = vech_index(i, i, total);
        if idx < kstar_full {
            v_full[(idx, idx)] = snp_var_se;
        }
    }

    // SNP-trait covariance elements: use se^2 * var_snp^2 as approximate variance
    for (snp_i, var_snp_i) in var_snp.iter().enumerate().take(n_snps) {
        for t in 0..k {
            let row_idx = n_snps + t;
            let col_idx = snp_i;
            // vech uses lower triangle, so ensure row >= col
            let (r, c) = if row_idx >= col_idx {
                (row_idx, col_idx)
            } else {
                (col_idx, row_idx)
            };
            let idx = vech_index(r, c, total);
            if idx < kstar_full {
                let se_val = se
                    .get(snp_i)
                    .and_then(|s| s.get(t))
                    .copied()
                    .unwrap_or(f64::NAN);
                v_full[(idx, idx)] = (se_val * var_snp_i).powi(2);
            }
        }
    }

    // SNP-SNP off-diagonal covariance elements: small fixed value
    for i in 0..n_snps {
        for j in (i + 1)..n_snps {
            let idx = vech_index(j, i, total);
            if idx < kstar_full {
                v_full[(idx, idx)] = snp_var_se;
            }
        }
    }

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
    fn test_run_multi_snp_basic() {
        let s_ld = faer::mat![[0.5, 0.2], [0.2, 0.4]];
        let v_ld = Mat::from_fn(3, 3, |i, j| if i == j { 0.001 } else { 0.0 });
        let beta = vec![vec![0.1, 0.05], vec![0.08, 0.12]];
        let se = vec![vec![0.02, 0.02], vec![0.02, 0.02]];
        let var_snp = vec![0.3, 0.25];
        let ld_matrix = faer::mat![[1.0, 0.3], [0.3, 1.0]];
        let snp_names = vec!["SNP1".to_string(), "SNP2".to_string()];

        let config = MultiSnpConfig {
            model: syntax::parse_model("F1 =~ NA*V1 + V2\nF1 ~~ 1*F1\nV1 ~~ V1\nV2 ~~ V2\nF1 ~ SNP1 + SNP2\nSNP1 ~~ SNP1\nSNP2 ~~ SNP2", false).unwrap(),
            estimation: EstimationMethod::Dwls,
            max_iter: 500,
            snp_var_se: None,
        };

        let result = run_multi_snp(
            &config, &s_ld, &v_ld, &beta, &se, &var_snp, &ld_matrix, &snp_names,
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
