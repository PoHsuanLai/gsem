use statrs::distribution::{Continuous, ContinuousCDF, Normal};

use crate::LdscResult;

/// Apply liability scale conversion to S and V matrices.
///
/// For binary traits with known sample prevalence and population prevalence:
///   conversion_factor = (K^2 * (1-K)^2) / (P * (1-P) * dnorm(qnorm(1-K))^2)
/// where K = population prevalence, P = sample prevalence.
///
/// S_liab[j,k] = S[j,k] * sqrt(cf_j) * sqrt(cf_k)
/// V_liab = V * outer(scale_ratios, scale_ratios)
pub fn apply_liability_scale(
    result: &mut LdscResult,
    sample_prev: &[Option<f64>],
    pop_prev: &[Option<f64>],
) {
    let k = result.s.nrows();
    if sample_prev.len() != k || pop_prev.len() != k {
        return;
    }

    let normal = Normal::standard();

    // Compute conversion factors per trait (1.0 for continuous traits)
    let cf: Vec<f64> = (0..k)
        .map(|i| {
            match (sample_prev[i], pop_prev[i]) {
                (Some(sp), Some(pp)) if sp > 0.0 && sp < 1.0 && pp > 0.0 && pp < 1.0 => {
                    let threshold = normal.inverse_cdf(1.0 - pp);
                    let z_density = normal.pdf(threshold);
                    (pp.powi(2) * (1.0 - pp).powi(2)) / (sp * (1.0 - sp) * z_density.powi(2))
                }
                _ => 1.0, // continuous trait, no conversion
            }
        })
        .collect();

    // Scale S matrix
    for j in 0..k {
        for jj in 0..k {
            result.s[(j, jj)] *= cf[j].sqrt() * cf[jj].sqrt();
        }
    }

    // Compute scale ratios for V matrix
    // scale_ratio[vech_idx] = sqrt(cf_j) * sqrt(cf_jj)
    let kstar = result.v.nrows();
    let mut scale_ratios = Vec::with_capacity(kstar);
    for j in 0..k {
        for jj in j..k {
            scale_ratios.push(cf[j].sqrt() * cf[jj].sqrt());
        }
    }

    // Scale V matrix: V_liab[i,j] = V[i,j] * scale_ratios[i] * scale_ratios[j]
    for i in 0..kstar {
        for j in 0..kstar {
            result.v[(i, j)] *= scale_ratios[i] * scale_ratios[j];
        }
    }
}

/// Compute the liability conversion factor for a single binary trait.
pub fn liability_conversion_factor(sample_prev: f64, pop_prev: f64) -> f64 {
    let normal = Normal::standard();
    let threshold = normal.inverse_cdf(1.0 - pop_prev);
    let z_density = normal.pdf(threshold);
    (pop_prev.powi(2) * (1.0 - pop_prev).powi(2))
        / (sample_prev * (1.0 - sample_prev) * z_density.powi(2))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_liability_factor_reasonable() {
        // For a rare disease with K=0.01, P=0.5 (case-control 50/50)
        let cf = liability_conversion_factor(0.5, 0.01);
        // Conversion factor should be positive and finite
        assert!(
            cf > 0.0 && cf.is_finite(),
            "cf should be positive and finite, got {cf}"
        );
    }

    #[test]
    fn test_no_conversion_for_continuous() {
        let mut result = LdscResult {
            s: faer::mat![[0.5]],
            v: faer::mat![[0.01]],
            i_mat: faer::mat![[1.0]],
            n_vec: vec![50000.0],
            m: 1000000.0,
        };
        let original_s = result.s[(0, 0)];
        apply_liability_scale(&mut result, &[None], &[None]);
        assert!((result.s[(0, 0)] - original_s).abs() < 1e-15);
    }
}
