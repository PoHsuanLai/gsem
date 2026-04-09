use anyhow::Result;

use crate::jackknife;
use crate::regression;
use crate::weights;

/// Result of heritability estimation for a single trait.
#[derive(Debug, Clone)]
pub struct H2Result {
    /// Heritability estimate (h2)
    pub h2: f64,
    /// LDSC intercept
    pub intercept: f64,
    /// Mean sample size
    pub mean_n: f64,
    /// Jackknife SE for h2
    pub se: f64,
    /// Jackknife pseudo-values for the h2 coefficient (used for V matrix)
    pub pseudo_values: Vec<f64>,
}

/// Estimate heritability for a single trait via LDSC regression.
///
/// Regresses chi-squared statistics on LD scores:
///   chi2 = h2_coef * LD + intercept
///   h2 = h2_coef * M / N_bar
pub fn estimate_h2(
    z: &[f64],
    n: &[f64],
    ld: &[f64],
    w_ld: &[f64],
    m: f64,
    n_blocks: usize,
) -> Result<H2Result> {
    let n_snps = z.len();
    assert_eq!(n_snps, n.len());
    assert_eq!(n_snps, ld.len());
    assert_eq!(n_snps, w_ld.len());

    // Compute chi-squared values
    let chi2: Vec<f64> = z.iter().map(|zi| zi * zi).collect();

    // Mean sample size
    let mean_n = n.iter().sum::<f64>() / n_snps as f64;

    // Compute weights
    let weights = weights::compute_h2_weights(&chi2, ld, w_ld, mean_n, m);

    // Run block regression
    let reg_result = regression::block_regression(ld, &chi2, &weights, n_blocks);

    // h2 = slope * M / N_bar
    let h2 = reg_result.coef[0] * m / mean_n;
    let intercept = reg_result.coef[1];

    // Jackknife for SE and pseudo-values
    let jk = jackknife::jackknife(&reg_result);

    let scale = m / mean_n;

    // Return RAW pseudo-values (unscaled regression slope).
    // The V matrix normalization is applied later in lib.rs, matching R's approach:
    //   v.out = cov(V.hold) / crossprod(N.vec * sqrt(n.blocks) / m)
    let pseudo_values = jk.pseudo_slope;

    Ok(H2Result {
        h2,
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
    fn test_h2_null() {
        // Under null (no heritability), chi2 ≈ 1 for all SNPs
        let n_snps = 1000;
        let z: Vec<f64> = (0..n_snps).map(|_| 1.0).collect();
        let n = vec![50000.0; n_snps];
        let ld: Vec<f64> = (0..n_snps).map(|i| 10.0 + (i as f64) * 0.01).collect();
        let w_ld = ld.clone();

        let result = estimate_h2(&z, &n, &ld, &w_ld, 1000000.0, 50).unwrap();
        // h2 should be close to 0 under null
        assert!(
            result.h2.abs() < 0.1,
            "h2 should be near 0, got {}",
            result.h2
        );
        // Intercept should be close to 1
        assert!(
            (result.intercept - 1.0).abs() < 0.5,
            "intercept should be near 1, got {}",
            result.intercept
        );
    }
}
