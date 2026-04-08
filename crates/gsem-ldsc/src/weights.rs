/// Compute LDSC regression weights for the heritability path (j==k).
///
/// From ldsc.R:
/// - tot_agg = M * (mean(chi2) - 1) / mean(L2 * N), clamped to [0, 1]
/// - c = tot_agg * N / M
/// - het_w = 1 / (2 * (1 + c * ld)^2)
/// - oc_w = 1 / max(w_ld, 1)
/// - w = het_w * oc_w, normalized: sqrt(w) / sum(sqrt(w))
pub fn compute_h2_weights(
    chi2: &[f64],
    ld: &[f64],
    w_ld: &[f64],
    n_bar: f64,
    m: f64,
) -> Vec<f64> {
    let n_snps = chi2.len();
    assert_eq!(n_snps, ld.len());
    assert_eq!(n_snps, w_ld.len());

    // Compute total aggregated heritability
    let mean_chi2: f64 = chi2.iter().sum::<f64>() / n_snps as f64;
    let mean_ld_n: f64 = ld.iter().zip(std::iter::repeat(n_bar)).map(|(l, n)| l * n).sum::<f64>()
        / n_snps as f64;

    let tot_agg = if mean_ld_n > 0.0 {
        ((m * (mean_chi2 - 1.0)) / mean_ld_n).clamp(0.0, 1.0)
    } else {
        0.0
    };

    let c = tot_agg * n_bar / m;

    let mut w = Vec::with_capacity(n_snps);
    for i in 0..n_snps {
        let het_w = 1.0 / (2.0 * (1.0 + c * ld[i]).powi(2));
        let oc_w = 1.0 / w_ld[i].max(1.0);
        w.push(het_w * oc_w);
    }

    // Normalize: sqrt(w) / sum(sqrt(w))
    let sqrt_w: Vec<f64> = w.iter().map(|x| x.sqrt()).collect();
    let sum_sqrt: f64 = sqrt_w.iter().sum();
    if sum_sqrt > 0.0 {
        sqrt_w.iter().map(|x| x / sum_sqrt).collect()
    } else {
        vec![1.0 / n_snps as f64; n_snps]
    }
}

/// Compute LDSC regression weights for the covariance path (j!=k).
///
/// Average of weights computed from each trait's chi-squared values.
pub fn compute_gcov_weights(
    chi2_j: &[f64],
    chi2_k: &[f64],
    ld: &[f64],
    w_ld: &[f64],
    n_j: f64,
    n_k: f64,
    m: f64,
) -> Vec<f64> {
    let w_j = compute_h2_weights(chi2_j, ld, w_ld, n_j, m);
    let w_k = compute_h2_weights(chi2_k, ld, w_ld, n_k, m);

    // Average weights
    let mut w_avg: Vec<f64> = w_j.iter().zip(w_k.iter()).map(|(a, b)| (a + b) / 2.0).collect();

    // Re-normalize
    let sum: f64 = w_avg.iter().sum();
    if sum > 0.0 {
        for x in &mut w_avg {
            *x /= sum;
        }
    }

    w_avg
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_weights_sum_to_one() {
        let chi2 = vec![1.5, 2.0, 1.2, 3.0, 1.8];
        let ld = vec![10.0, 20.0, 15.0, 25.0, 12.0];
        let w_ld = vec![10.0, 20.0, 15.0, 25.0, 12.0];
        let w = compute_h2_weights(&chi2, &ld, &w_ld, 50000.0, 1000000.0);
        let sum: f64 = w.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_weights_positive() {
        let chi2 = vec![1.0; 10];
        let ld = vec![5.0; 10];
        let w_ld = vec![5.0; 10];
        let w = compute_h2_weights(&chi2, &ld, &w_ld, 10000.0, 500000.0);
        assert!(w.iter().all(|&x| x > 0.0));
    }
}
