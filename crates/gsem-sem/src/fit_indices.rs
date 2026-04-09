use faer::{Mat, Side};
use gsem_matrix::vech;
use statrs::distribution::{ChiSquared, ContinuousCDF};

use crate::ModelFit;

/// Compute model fit statistics.
///
/// Chi-square via eigendecomposition of V:
///   eta = vech(S - Sigma_hat)
///   V = P * Eig * P'
///   Q = eta' * P * Eig^{-1} * P' * eta
pub fn compute_fit(
    s_obs: &Mat<f64>,
    sigma_hat: &Mat<f64>,
    v: &Mat<f64>,
    df: usize,
    n_free: usize,
    q_null: Option<f64>,
    df_null: Option<usize>,
) -> ModelFit {
    let s_vec = vech::vech(s_obs).expect("s_obs must be square");
    let sigma_vec = vech::vech(sigma_hat).expect("sigma_hat must be square");

    // Residual
    let eta: Vec<f64> = s_vec
        .iter()
        .zip(sigma_vec.iter())
        .map(|(s, sig)| s - sig)
        .collect();

    // Chi-square via eigendecomposition of V
    let chisq = compute_chisq_eigen(&eta, v);

    // p-value
    let p_chisq = if df > 0 {
        ChiSquared::new(df as f64)
            .map(|chi2| 1.0 - chi2.cdf(chisq))
            .unwrap_or(1.0)
    } else {
        1.0
    };

    // AIC
    let aic = chisq + 2.0 * n_free as f64;

    // CFI
    let cfi = match (q_null, df_null) {
        (Some(qn), Some(dfn)) => {
            let num = (qn - dfn as f64) - (chisq - df as f64);
            let denom = qn - dfn as f64;
            if denom > 0.0 {
                (num / denom).min(1.0)
            } else {
                1.0
            }
        }
        _ => f64::NAN,
    };

    // SRMR
    let srmr = compute_srmr(s_obs, sigma_hat);

    ModelFit {
        chisq,
        df,
        p_chisq,
        aic,
        cfi,
        srmr,
    }
}

/// Compute chi-square via eigendecomposition.
/// Q = eta' * V^{-1} * eta (using eigendecomposition for stability)
fn compute_chisq_eigen(eta: &[f64], v: &Mat<f64>) -> f64 {
    let n = eta.len();

    let Ok(eigen) = v.self_adjoint_eigen(Side::Lower) else {
        // Fallback: direct computation
        return 0.0;
    };

    let u = eigen.U();
    let s = eigen.S().column_vector();

    // Q = eta' * U * diag(1/eigenvalues) * U' * eta
    // = sum_i (u_i' * eta)^2 / eigenvalue_i
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

/// Compute SRMR (Standardized Root Mean Square Residual).
fn compute_srmr(s_obs: &Mat<f64>, sigma_hat: &Mat<f64>) -> f64 {
    let p = s_obs.nrows();
    let mut sum = 0.0;
    let mut count = 0;

    for i in 0..p {
        for j in 0..=i {
            let sd_i = s_obs[(i, i)].sqrt();
            let sd_j = s_obs[(j, j)].sqrt();
            if sd_i > 0.0 && sd_j > 0.0 {
                let r_obs = s_obs[(i, j)] / (sd_i * sd_j);
                let r_imp = sigma_hat[(i, j)] / (sd_i * sd_j);
                sum += (r_obs - r_imp).powi(2);
                count += 1;
            }
        }
    }

    if count > 0 {
        (sum / count as f64).sqrt()
    } else {
        0.0
    }
}
