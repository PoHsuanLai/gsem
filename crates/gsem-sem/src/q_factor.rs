use faer::{Mat, Side};
use gsem_matrix::error::MatrixError;
use gsem_matrix::vech;
use statrs::distribution::{ChiSquared, ContinuousCDF};

use crate::syntax::{Op, ParTable};

/// Result of Q_Factor test for a pair of factors.
#[derive(Debug, Clone)]
pub struct QFactorResult {
    /// Name of the first factor
    pub factor1: String,
    /// Name of the second factor
    pub factor2: String,
    /// Q chi-square test statistic
    pub q_chisq: f64,
    /// Degrees of freedom for the Q test
    pub q_df: usize,
    /// P-value for the Q test
    pub q_p: f64,
}

/// Extract factor-to-indicator mappings from a parameter table.
///
/// Returns a vector of `(factor_name, indicator_obs_indices)` where
/// each indicator index refers to its position in `obs_names`.
pub fn factor_indicators(pt: &ParTable, obs_names: &[String]) -> Vec<(String, Vec<usize>)> {
    let mut factors: Vec<(String, Vec<usize>)> = Vec::new();
    let mut seen_factors = std::collections::HashSet::new();

    for row in &pt.rows {
        if row.op != Op::Loading {
            continue;
        }
        let factor_name = &row.lhs;
        let indicator_name = &row.rhs;

        // Find indicator index in obs_names
        let Some(obs_idx) = obs_names.iter().position(|n| n == indicator_name) else {
            continue;
        };

        if !seen_factors.contains(factor_name) {
            seen_factors.insert(factor_name.clone());
            factors.push((factor_name.clone(), Vec::new()));
        }

        if let Some(entry) = factors
            .iter_mut()
            .find(|(name, _)| name == factor_name)
            .filter(|entry| !entry.1.contains(&obs_idx))
        {
            entry.1.push(obs_idx);
        }
    }

    factors
}

/// Compute Q_Factor heterogeneity statistics for all factor pairs.
///
/// Tests whether cross-factor indicator relationships are adequately
/// captured by the factor model.
///
/// For each factor pair (Fi, Fj):
/// 1. Identify indicators that load on Fi and indicators that load on Fj
/// 2. Compute residuals for cross-factor indicator pairs
/// 3. Q = eta' * V_sub^{-1} * eta (chi-square test)
pub fn compute_q_factor(
    s_obs: &Mat<f64>,
    sigma_hat: &Mat<f64>,
    v: &Mat<f64>,
    factor_inds: &[(String, Vec<usize>)],
) -> Result<Vec<QFactorResult>, MatrixError> {
    let k = s_obs.nrows();
    let mut results = Vec::new();

    // Build vech index lookup
    let vech_idx_pairs = vech::vech_indices(k);

    let n_factors = factor_inds.len();
    for fi in 0..n_factors {
        for fj in (fi + 1)..n_factors {
            let (ref name_i, ref inds_i) = factor_inds[fi];
            let (ref name_j, ref inds_j) = factor_inds[fj];

            // Cross-factor pairs: indicators of Fi x indicators of Fj
            // We want (row, col) pairs in vech ordering (row >= col)
            let mut cross_pairs: Vec<(usize, usize)> = Vec::new();
            for &a in inds_i {
                for &b in inds_j {
                    if a != b {
                        let pair = (a.max(b), a.min(b));
                        cross_pairs.push(pair);
                    }
                }
            }
            cross_pairs.sort();
            cross_pairs.dedup();

            if cross_pairs.is_empty() {
                continue;
            }

            // Map cross-factor pairs to vech indices
            let vech_indices: Vec<usize> = cross_pairs
                .iter()
                .map(|&(r, c)| {
                    vech_idx_pairs
                        .iter()
                        .position(|&(vi, vj)| vi == r && vj == c)
                        .expect("cross-factor pair must exist in vech indices")
                })
                .collect();

            let n_cross = vech_indices.len();

            // Residual vector: eta = vech(S - Sigma)[cross indices]
            let s_vec = vech::vech(s_obs)?;
            let sig_vec = vech::vech(sigma_hat)?;
            let eta: Vec<f64> = vech_indices
                .iter()
                .map(|&idx| s_vec[idx] - sig_vec[idx])
                .collect();

            // V submatrix for cross-factor indices
            let v_sub = Mat::from_fn(n_cross, n_cross, |i, j| {
                v[(vech_indices[i], vech_indices[j])]
            });

            // Q = eta' * V_sub^{-1} * eta
            let q_chisq = compute_q_from_eigen(&eta, &v_sub);
            let q_df = n_cross;

            let q_p = if q_df > 0 && q_chisq.is_finite() {
                ChiSquared::new(q_df as f64)
                    .map(|chi2| 1.0 - chi2.cdf(q_chisq))
                    .unwrap_or(1.0)
            } else {
                1.0
            };

            results.push(QFactorResult {
                factor1: name_i.clone(),
                factor2: name_j.clone(),
                q_chisq,
                q_df,
                q_p,
            });
        }
    }

    Ok(results)
}

/// Compute Q = eta' * V^{-1} * eta via eigendecomposition for numerical stability.
fn compute_q_from_eigen(eta: &[f64], v: &Mat<f64>) -> f64 {
    let n = eta.len();
    let Ok(eigen) = v.self_adjoint_eigen(Side::Lower) else {
        return 0.0;
    };

    let u = eigen.U();
    let s = eigen.S().column_vector();

    let mut q = 0.0;
    for i in 0..n {
        let eigenval = s[i];
        if eigenval > 1e-10 {
            let mut dot = 0.0;
            for j in 0..n {
                dot += u[(j, i)] * eta[j];
            }
            q += dot * dot / eigenval;
        }
    }

    q.max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_factor_indicators() {
        use crate::syntax::parse_model;

        let model = "F1 =~ V1 + V2 + V3\nF2 =~ V4 + V5 + V6";
        let pt = parse_model(model, false).unwrap();
        let obs_names: Vec<String> = (1..=6).map(|i| format!("V{i}")).collect();

        let fi = factor_indicators(&pt, &obs_names);
        assert_eq!(fi.len(), 2);
        assert_eq!(fi[0].0, "F1");
        assert_eq!(fi[0].1, vec![0, 1, 2]);
        assert_eq!(fi[1].0, "F2");
        assert_eq!(fi[1].1, vec![3, 4, 5]);
    }

    #[test]
    fn test_factor_indicators_cross_loading() {
        use crate::syntax::parse_model;

        // V3 loads on both factors
        let model = "F1 =~ V1 + V2 + V3\nF2 =~ V3 + V4 + V5";
        let pt = parse_model(model, false).unwrap();
        let obs_names: Vec<String> = (1..=5).map(|i| format!("V{i}")).collect();

        let fi = factor_indicators(&pt, &obs_names);
        assert_eq!(fi.len(), 2);
        assert_eq!(fi[0].1, vec![0, 1, 2]); // V1, V2, V3
        assert_eq!(fi[1].1, vec![2, 3, 4]); // V3, V4, V5
    }

    #[test]
    fn test_compute_q_factor_zero_residual() {
        // When implied == observed, Q should be 0
        let s = faer::mat![
            [1.0, 0.5, 0.3, 0.1],
            [0.5, 1.0, 0.2, 0.15],
            [0.3, 0.2, 1.0, 0.4],
            [0.1, 0.15, 0.4, 1.0],
        ];
        let sigma = s.clone();
        let kstar = 4 * 5 / 2; // 10
        let v = Mat::from_fn(kstar, kstar, |i, j| if i == j { 0.01 } else { 0.0 });

        let fi = vec![
            ("F1".to_string(), vec![0, 1]),
            ("F2".to_string(), vec![2, 3]),
        ];

        let results = compute_q_factor(&s, &sigma, &v, &fi).unwrap();
        assert_eq!(results.len(), 1);
        assert!(
            results[0].q_chisq.abs() < 1e-10,
            "Q should be 0 when residuals are 0, got {}",
            results[0].q_chisq
        );
        assert_eq!(results[0].factor1, "F1");
        assert_eq!(results[0].factor2, "F2");
    }

    #[test]
    fn test_compute_q_factor_nonzero_residual() {
        // Observed and implied differ in cross-factor entries
        let s = faer::mat![
            [1.0, 0.5, 0.3, 0.1],
            [0.5, 1.0, 0.2, 0.15],
            [0.3, 0.2, 1.0, 0.4],
            [0.1, 0.15, 0.4, 1.0],
        ];
        let mut sigma = s.clone();
        // Perturb cross-factor elements
        sigma[(0, 2)] = 0.1;
        sigma[(2, 0)] = 0.1;
        sigma[(0, 3)] = 0.0;
        sigma[(3, 0)] = 0.0;
        sigma[(1, 2)] = 0.1;
        sigma[(2, 1)] = 0.1;
        sigma[(1, 3)] = 0.05;
        sigma[(3, 1)] = 0.05;

        let kstar = 4 * 5 / 2;
        let v = Mat::from_fn(kstar, kstar, |i, j| if i == j { 0.01 } else { 0.0 });

        let fi = vec![
            ("F1".to_string(), vec![0, 1]),
            ("F2".to_string(), vec![2, 3]),
        ];

        let results = compute_q_factor(&s, &sigma, &v, &fi).unwrap();
        assert_eq!(results.len(), 1);
        assert!(
            results[0].q_chisq > 0.0,
            "Q should be positive with nonzero residuals"
        );
        assert_eq!(results[0].q_df, 4); // 2x2 cross-factor pairs
        assert!(results[0].q_p >= 0.0 && results[0].q_p <= 1.0);
    }
}
