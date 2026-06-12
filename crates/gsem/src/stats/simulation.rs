use faer::{Mat, Side};
use rand::RngExt;
use rand_distr::StandardNormal;

/// Configuration for simLDSC.
pub struct SimConfig {
    /// LDSC intercept matrix (k x k). Defaults to identity if None.
    pub intercepts: Option<Mat<f64>>,
    /// Phenotypic correlation matrix (k x k). Defaults to None (no environmental correlation).
    pub r_pheno: Option<Mat<f64>>,
    /// Sample overlap proportion (0 to 1). Only used when r_pheno is set.
    pub n_overlap: f64,
}

impl Default for SimConfig {
    fn default() -> Self {
        Self {
            intercepts: None,
            r_pheno: None,
            n_overlap: 0.0,
        }
    }
}

/// Per-SNP variance-covariance matrix of Z statistics for a single LD score.
///
/// This is the deterministic core of GenomicSEM's `simLDSC()` (the matrix the
/// random Z vector is drawn from). For an LD score `ld`:
///   - `Sigma[i,i] = int_i + (N_i * S_ii / M) * ld`
///   - `Sigma[i,j] = (sqrt(N_i N_j) * S_ij / M) * ld + rPheno_ij * N_ij/sqrt(N_i N_j)`
///
/// `intercepts` carries the per-trait LDSC intercepts on its diagonal (off-diagonal
/// ignored). The phenotypic-overlap term uses the scalar `n_overlap`, which equals
/// R's `N_ij/sqrt(N_i N_j)` when the sample-overlap matrix is constructed that way.
pub fn per_snp_z_cov(
    s: &Mat<f64>,
    n_per_trait: &[f64],
    ld: f64,
    m: f64,
    intercepts: &Mat<f64>,
    r_pheno: Option<&Mat<f64>>,
    n_overlap: f64,
) -> Mat<f64> {
    let k = s.nrows();
    Mat::from_fn(k, k, |i, j| {
        let sqrt_nn = (n_per_trait[i] * n_per_trait[j]).sqrt();
        let genetic = s[(i, j)] / m * ld * sqrt_nn;
        if i == j {
            intercepts[(i, i)] + genetic
        } else {
            let env = r_pheno.map_or(0.0, |rp| rp[(i, j)] * n_overlap);
            genetic + env
        }
    })
}

/// Simulate GWAS summary statistics under a specified genetic model.
///
/// Port of GenomicSEM's `simLDSC()`. For each SNP the per-SNP Z covariance is
/// built by [`per_snp_z_cov`]; an independent multivariate-normal Z vector is
/// then drawn from it. (The random draw is platform-dependent; the deterministic
/// `per_snp_z_cov` is the R-equivalent piece.)
pub fn simulate_sumstats(
    s: &Mat<f64>,
    n_per_trait: &[f64],
    ld_scores: &[f64],
    m: f64,
    config: &SimConfig,
) -> Vec<Vec<f64>> {
    let k = s.nrows();
    let n_snps = ld_scores.len();
    let mut rng = rand::rng();

    // Intercept matrix (default: identity)
    let int_mat = config
        .intercepts
        .as_ref()
        .cloned()
        .unwrap_or_else(|| Mat::<f64>::identity(k, k));

    let mut z_all = vec![vec![0.0; n_snps]; k];
    let mut ind = vec![0.0; k];

    for s_idx in 0..n_snps {
        let ld = ld_scores[s_idx];

        // Build per-SNP covariance (deterministic R-equivalent construction).
        let sigma_z = per_snp_z_cov(
            s,
            n_per_trait,
            ld,
            m,
            &int_mat,
            config.r_pheno.as_ref(),
            config.n_overlap,
        );

        // Cholesky of per-SNP covariance
        let chol = cholesky_or_sqrt(&sigma_z, k);

        // Generate independent standard normals
        for v in ind.iter_mut() {
            *v = rng.sample(StandardNormal);
        }

        // Correlate using Cholesky factor
        for j in 0..k {
            z_all[j][s_idx] = (0..k).map(|l| chol[(j, l)] * ind[l]).sum();
        }
    }

    z_all
}

/// Cholesky decomposition with eigendecomposition fallback.
fn cholesky_or_sqrt(mat: &Mat<f64>, k: usize) -> Mat<f64> {
    if let Ok(llt) = mat.llt(Side::Lower) {
        let l = llt.L();
        return Mat::from_fn(k, k, |i, j| l[(i, j)]);
    }
    // Eigendecomposition square root fallback
    let Ok(eigen) = mat.self_adjoint_eigen(Side::Lower) else {
        return Mat::zeros(k, k);
    };
    let u = eigen.U();
    let sv = eigen.S().column_vector();
    Mat::from_fn(k, k, |i, j| {
        (0..k)
            .map(|l| u[(i, l)] * sv[l].max(0.0).sqrt() * u[(j, l)])
            .sum()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simulate_produces_correct_dimensions() {
        let s = faer::mat![[0.3, 0.1], [0.1, 0.4]];
        let n = vec![50000.0, 50000.0];
        let ld = vec![10.0; 100];
        let z = simulate_sumstats(&s, &n, &ld, 1000000.0, &SimConfig::default());
        assert_eq!(z.len(), 2);
        assert_eq!(z[0].len(), 100);
        assert_eq!(z[1].len(), 100);
    }

    #[test]
    fn test_simulate_z_finite() {
        let s = faer::mat![[0.5]];
        let n = vec![10000.0];
        let ld = vec![5.0; 50];
        let z = simulate_sumstats(&s, &n, &ld, 500000.0, &SimConfig::default());
        assert!(z[0].iter().all(|z| z.is_finite()));
    }

    #[test]
    fn test_simulate_with_intercept() {
        let s = faer::mat![[0.3, 0.1], [0.1, 0.4]];
        let n = vec![50000.0, 50000.0];
        let ld = vec![10.0; 50];
        let int = faer::mat![[1.1, 0.05], [0.05, 1.2]];
        let config = SimConfig {
            intercepts: Some(int),
            ..Default::default()
        };
        let z = simulate_sumstats(&s, &n, &ld, 1000000.0, &config);
        assert_eq!(z.len(), 2);
        assert!(z[0].iter().all(|z| z.is_finite()));
    }

    #[test]
    fn test_simulate_with_rpheno_overlap() {
        let s = faer::mat![[0.3, 0.1], [0.1, 0.4]];
        let n = vec![50000.0, 50000.0];
        let ld = vec![10.0; 50];
        let r_pheno = faer::mat![[1.0, 0.5], [0.5, 1.0]];
        let config = SimConfig {
            r_pheno: Some(r_pheno),
            n_overlap: 0.8,
            ..Default::default()
        };
        let z = simulate_sumstats(&s, &n, &ld, 1000000.0, &config);
        assert_eq!(z.len(), 2);
        assert!(z[0].iter().all(|z| z.is_finite()));
    }
}
