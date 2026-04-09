/// Compute LDSC regression weights for the heritability path (j==k).
///
/// From ldsc.R:
/// - `tot_agg = M * (mean(chi2) - 1) / mean(L2 * N)`, clamped to 0..1
/// - c[i] = tot_agg * N[i] / M   (per-SNP)
/// - het_w[i] = 1 / (2 * (1 + c[i] * ld[i])^2)
/// - oc_w[i] = 1 / max(w_ld[i], 1)
/// - w = het_w * oc_w, normalized: sqrt(w) / sum(sqrt(w))
pub fn compute_h2_weights(chi2: &[f64], ld: &[f64], w_ld: &[f64], n: &[f64], m: f64) -> Vec<f64> {
    let n_snps = chi2.len();
    assert_eq!(n_snps, ld.len());
    assert_eq!(n_snps, w_ld.len());
    assert_eq!(n_snps, n.len());

    // Compute total aggregated heritability
    // R: tot.agg <- (M.tot*(mean(merged$chi1)-1))/mean(merged$L2*merged$N)
    let mean_chi2: f64 = chi2.iter().sum::<f64>() / n_snps as f64;
    let mean_ld_n: f64 = ld
        .iter()
        .zip(n.iter())
        .map(|(&l, &ni)| l * ni)
        .sum::<f64>()
        / n_snps as f64;

    let tot_agg = if mean_ld_n > 0.0 {
        ((m * (mean_chi2 - 1.0)) / mean_ld_n).clamp(0.0, 1.0)
    } else {
        0.0
    };

    // R: merged$c <- tot.agg*merged$N/M.tot  (per-SNP)
    // R: merged$het.w <- 1/(2*(1+(merged$c*merged$ld))^2)
    let mut sqrt_w: Vec<f64> = ld
        .iter()
        .zip(w_ld.iter())
        .zip(n.iter())
        .map(|((&l, &wl), &ni)| {
            let c_i = tot_agg * ni / m;
            let het_w = 1.0 / (2.0 * (1.0 + c_i * l.max(1.0)).powi(2));
            let oc_w = 1.0 / wl.max(1.0);
            (het_w * oc_w).sqrt()
        })
        .collect();

    let sum_sqrt: f64 = sqrt_w.iter().sum();
    if sum_sqrt > 0.0 {
        for x in &mut sqrt_w {
            *x /= sum_sqrt;
        }
        sqrt_w
    } else {
        vec![1.0 / n_snps as f64; n_snps]
    }
}

/// Compute LDSC regression weights for the covariance path (j!=k).
///
/// R computes initial.w from trait 1 and initial.w2 from trait 2,
/// then: weights_cov = (initial.w + initial.w2) / sum(initial.w + initial.w2)
/// Due to R's `$` partial matching, both LD and chi use weights_cov.
pub fn compute_gcov_weights(
    chi2_j: &[f64],
    chi2_k: &[f64],
    ld: &[f64],
    w_ld: &[f64],
    n_j: &[f64],
    n_k: &[f64],
    m: f64,
) -> Vec<f64> {
    let w_j = compute_h2_weights(chi2_j, ld, w_ld, n_j, m);
    let w_k = compute_h2_weights(chi2_k, ld, w_ld, n_k, m);

    // R: merged$weights_cov <- (initial.w + initial.w2) / sum(initial.w + initial.w2)
    // Note: compute_h2_weights returns normalized sqrt(w), i.e. initial.w/sum(initial.w)
    // We need unnormalized initial.w + initial.w2, then renormalize.
    // Since w_j = initial.w_j / sum(initial.w_j) and w_k = initial.w_k / sum(initial.w_k),
    // averaging and renormalizing is equivalent to R's formula.
    let mut w_avg: Vec<f64> = w_j
        .iter()
        .zip(w_k.iter())
        .map(|(a, b)| a + b)
        .collect();

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
        let n = vec![50000.0; 5];
        let w = compute_h2_weights(&chi2, &ld, &w_ld, &n, 1000000.0);
        let sum: f64 = w.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_weights_positive() {
        let chi2 = vec![1.0; 10];
        let ld = vec![5.0; 10];
        let w_ld = vec![5.0; 10];
        let n = vec![10000.0; 10];
        let w = compute_h2_weights(&chi2, &ld, &w_ld, &n, 500000.0);
        assert!(w.iter().all(|&x| x > 0.0));
    }
}
