use faer::{Mat, Side};
use rand::Rng;
use rand_distr::StandardNormal;

/// Simulate GWAS summary statistics under a specified genetic model.
///
/// Port of GenomicSEM's `simLDSC()`.
///
/// Generates multivariate normal Z-statistics with covariance structure
/// determined by the genetic model and LD scores.
pub fn simulate_sumstats(
    s: &Mat<f64>,
    n_per_trait: &[f64],
    ld_scores: &[f64],
    m: f64,
) -> Vec<Vec<f64>> {
    let k = s.nrows();
    let n_snps = ld_scores.len();
    let mut rng = rand::rng();

    // Cholesky of S for generating correlated Z-scores
    let chol = match s.llt(Side::Lower) {
        Ok(llt) => {
            let l = llt.L();
            Mat::from_fn(k, k, |i, j| l[(i, j)])
        }
        Err(_) => {
            // If not PD, use eigendecomposition
            let eigen = s.self_adjoint_eigen(Side::Lower).unwrap();
            let u = eigen.U();
            let sv = eigen.S().column_vector();
            let sqrt_s = Mat::from_fn(k, k, |i, j| {
                let mut val = 0.0;
                for l in 0..k {
                    val += u[(i, l)] * sv[l].max(0.0).sqrt() * u[(j, l)];
                }
                val
            });
            sqrt_s
        }
    };

    // Generate Z-scores per SNP
    let mut z_all = vec![vec![0.0; n_snps]; k];

    for s_idx in 0..n_snps {
        // Expected chi2 contribution from this SNP's LD score
        let ld = ld_scores[s_idx];

        // Generate independent standard normals
        let ind: Vec<f64> = (0..k).map(|_| rng.sample(StandardNormal)).collect();

        // Correlate using Cholesky factor and scale by LD
        for j in 0..k {
            let mut z = 0.0;
            for l in 0..k {
                z += chol[(j, l)] * ind[l];
            }
            // Scale: Z ~ N(0, 1 + h2/M * LD * N)
            let scale = (1.0 + s[(j, j)] / m * ld * n_per_trait[j]).sqrt();
            z_all[j][s_idx] = z * scale;
        }
    }

    z_all
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simulate_produces_correct_dimensions() {
        let s = faer::mat![[0.3, 0.1], [0.1, 0.4]];
        let n = vec![50000.0, 50000.0];
        let ld = vec![10.0; 100];
        let z = simulate_sumstats(&s, &n, &ld, 1000000.0);
        assert_eq!(z.len(), 2);
        assert_eq!(z[0].len(), 100);
        assert_eq!(z[1].len(), 100);
    }

    #[test]
    fn test_simulate_z_finite() {
        let s = faer::mat![[0.5]];
        let n = vec![10000.0];
        let ld = vec![5.0; 50];
        let z = simulate_sumstats(&s, &n, &ld, 500000.0);
        assert!(z[0].iter().all(|z| z.is_finite()));
    }
}
