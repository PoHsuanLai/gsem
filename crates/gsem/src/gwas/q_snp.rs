use faer::{Mat, Side};
use gsem_matrix::error::MatrixError;
use statrs::distribution::{ChiSquared, ContinuousCDF};

/// Compute the Q_SNP heterogeneity statistic for one factor.
///
/// Tests whether a SNP's effects on a factor's indicators are homogeneous with
/// the common-pathway model (i.e. fully mediated by the factor). This mirrors R
/// GenomicSEM's `userGWAS_main.R` Q_SNP block exactly:
///
/// ```text
/// eta       = residual[SNP, indicators]          # length m (SNP→indicator residual cov)
/// V_SNP_i   = v_full[indicators, indicators]     # m×m SNP sampling-variance block
/// Q_SNP     = eta' · P · Eig^{-1} · P' · eta      # P, Eig = eigvecs/vals of V_SNP_i
/// df        = m - 1
/// ```
///
/// Inputs (all in `s_full` / `v_full` coordinates, where the SNP is variable 0
/// and the `k` phenotypes are variables `1..=k`):
/// * `residual`   = `S_full - Sigma_hat`, the `(k+1)×(k+1)` residual covariance.
/// * `v_full`     = the `(k+1)(k+2)/2`-sized sampling covariance; its `[1..=k, 1..=k]`
///   block is R's `V_SNP` (the SNP sampling-variance matrix indexed by phenotype).
/// * `indicators` = 0-based phenotype indices of this factor's indicators (so
///   phenotype `p` lives at `s_full` index `p + 1` and `v_full` block row `p`).
///
/// Returns `(Q_SNP, df, p_value)`.
pub fn compute_q_snp(
    residual: &Mat<f64>,
    v_full: &Mat<f64>,
    indicators: &[usize],
) -> Result<(f64, usize, f64), MatrixError> {
    let m = indicators.len();
    if m == 0 {
        return Ok((f64::NAN, 0, 1.0));
    }

    // SNP is row 0; phenotype p is at index p+1 in both s_full and the v_full
    // V_SNP block (which occupies v_full[1..=k, 1..=k]).
    let eta: Vec<f64> = indicators.iter().map(|&p| residual[(0, p + 1)]).collect();
    let v_snp_i = Mat::from_fn(m, m, |i, j| v_full[(indicators[i] + 1, indicators[j] + 1)]);

    let Ok(eigen) = v_snp_i.self_adjoint_eigen(Side::Lower) else {
        return Ok((f64::NAN, m.saturating_sub(1), 1.0));
    };
    let u = eigen.U();
    let s_diag = eigen.S().column_vector();

    // Q = sum_i (P_i · eta)^2 / lambda_i over all m eigencomponents.
    let mut q = 0.0_f64;
    for i in 0..m {
        let lambda = s_diag[i];
        if lambda.abs() < 1e-12 {
            // A non-positive / zero eigenvalue means V_SNP_i is singular; R would
            // error on solve(). Bail to NaN rather than fabricate a value.
            return Ok((f64::NAN, m.saturating_sub(1), 1.0));
        }
        let dot: f64 = (0..m).map(|j| u[(j, i)] * eta[j]).sum();
        q += dot * dot / lambda;
    }

    let df = m - 1;
    let p_val = if df > 0 && q.is_finite() {
        ChiSquared::new(df as f64)
            .map(|chi2| 1.0 - chi2.cdf(q))
            .unwrap_or(1.0)
    } else {
        1.0
    };

    Ok((q, df, p_val))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_q_snp_identity_v() {
        // With V_SNP_i = I, Q = eta'eta = sum of squared SNP→indicator residuals.
        // residual is (k+1)×(k+1); SNP is row/col 0, k=3 indicators at 1..=3.
        // Put eta = [0.2, -0.3, 0.5] in the SNP row.
        let residual = faer::mat![
            [0.0, 0.2, -0.3, 0.5],
            [0.2, 0.0, 0.0, 0.0],
            [-0.3, 0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0, 0.0],
        ];
        // v_full's [1..=3, 1..=3] block = identity.
        let kstar = (3 + 1) * (3 + 2) / 2; // 10
        let v_full = Mat::from_fn(kstar, kstar, |i, j| {
            if (1..=3).contains(&i) && i == j {
                1.0
            } else if i == j {
                1e-3 // other diagonal entries, irrelevant to Q
            } else {
                0.0
            }
        });
        let (q, df, _p) = compute_q_snp(&residual, &v_full, &[0, 1, 2]).unwrap();
        let expected = 0.2_f64.powi(2) + 0.3_f64.powi(2) + 0.5_f64.powi(2);
        assert!((q - expected).abs() < 1e-12, "Q={q} expected={expected}");
        assert_eq!(df, 2, "df should be #indicators - 1");
    }

    #[test]
    fn test_q_snp_empty_indicators() {
        let residual = Mat::<f64>::zeros(2, 2);
        let v_full = Mat::<f64>::identity(3, 3);
        let (q, df, p) = compute_q_snp(&residual, &v_full, &[]).unwrap();
        assert!(q.is_nan());
        assert_eq!(df, 0);
        assert_eq!(p, 1.0);
    }
}
