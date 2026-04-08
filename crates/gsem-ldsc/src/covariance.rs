use anyhow::Result;

use crate::jackknife;
use crate::regression;
use crate::weights;

/// Result of genetic covariance estimation between two traits.
#[derive(Debug, Clone)]
pub struct GcovResult {
    /// Genetic covariance estimate
    pub gcov: f64,
    /// Cross-trait LDSC intercept
    pub intercept: f64,
    /// Geometric mean of sample sizes
    pub mean_n: f64,
    /// Jackknife SE
    pub se: f64,
    /// Jackknife pseudo-values (used for V matrix)
    pub pseudo_values: Vec<f64>,
}

/// Estimate genetic covariance between two traits via LDSC regression.
///
/// Regresses Z1*Z2 product on LD scores:
///   Z1*Z2 = gcov_coef * LD + intercept
///   gcov = gcov_coef * M / sqrt(N1 * N2)
pub fn estimate_gcov(
    z1: &[f64],
    z2: &[f64],
    n1: &[f64],
    n2: &[f64],
    ld: &[f64],
    w_ld: &[f64],
    m: f64,
    n_blocks: usize,
) -> Result<GcovResult> {
    let n_snps = z1.len();
    assert_eq!(n_snps, z2.len());

    // Compute Z1*Z2 product (cross-trait chi statistic)
    let zz: Vec<f64> = z1.iter().zip(z2.iter()).map(|(a, b)| a * b).collect();

    // Mean sample sizes
    let mean_n1 = n1.iter().sum::<f64>() / n_snps as f64;
    let mean_n2 = n2.iter().sum::<f64>() / n_snps as f64;
    let mean_n = (mean_n1 * mean_n2).sqrt();

    // Compute chi-squared for each trait (for weight computation)
    let chi2_1: Vec<f64> = z1.iter().map(|z| z * z).collect();
    let chi2_2: Vec<f64> = z2.iter().map(|z| z * z).collect();

    // Average weights from both traits
    let weights = weights::compute_gcov_weights(&chi2_1, &chi2_2, ld, w_ld, mean_n1, mean_n2, m);

    // Run block regression
    let reg_result = regression::block_regression(ld, &zz, &weights, n_blocks);

    // gcov = slope * M / sqrt(N1 * N2)
    let gcov = reg_result.coef[0] * m / mean_n;
    let intercept = reg_result.coef[1];

    // Jackknife
    let jk = jackknife::jackknife(&reg_result);

    // Scale pseudo-values
    let scale = m / mean_n;
    let pseudo_values: Vec<f64> = jk.pseudo_slope.iter().map(|p| p * scale).collect();

    Ok(GcovResult {
        gcov,
        intercept,
        mean_n,
        se: jk.se[0] * scale,
        pseudo_values,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gcov_independent_traits() {
        // Two independent traits should have gcov ≈ 0
        let n_snps = 500;
        // Alternating signs, no correlation
        let z1: Vec<f64> = (0..n_snps)
            .map(|i| if i % 2 == 0 { 1.0 } else { -1.0 })
            .collect();
        let z2: Vec<f64> = (0..n_snps)
            .map(|i| if i % 3 == 0 { 1.0 } else { -1.0 })
            .collect();
        let n1 = vec![50000.0; n_snps];
        let n2 = vec![50000.0; n_snps];
        let ld: Vec<f64> = (0..n_snps).map(|i| 10.0 + (i as f64) * 0.02).collect();
        let w_ld = ld.clone();

        let result = estimate_gcov(&z1, &z2, &n1, &n2, &ld, &w_ld, 1000000.0, 50).unwrap();
        assert!(
            result.gcov.abs() < 0.5,
            "gcov should be near 0, got {}",
            result.gcov
        );
    }
}
