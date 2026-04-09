//! Model-implied genetic correlation matrix.
//!
//! Wrapper around [`crate::commonfactor`] or a user-specified model that additionally computes:
//! - R: the model-implied genetic correlation matrix (standardized implied covariance)
//! - V_R: the sampling covariance of R (via delta method / numerical Jacobian)
//!
//! Port of R GenomicSEM's `rgmodel()`.

use anyhow::Result;
use faer::Mat;
use gsem_matrix::vech;

/// Result of rgmodel.
#[derive(Debug, Clone)]
pub struct RgModelResult {
    /// Model-implied genetic correlation matrix (k x k)
    pub r: Mat<f64>,
    /// Sampling covariance of vech(R) (kstar x kstar)
    pub v_r: Mat<f64>,
    /// The underlying SEM result
    pub sem_result: crate::SemResult,
}

/// Compute model-implied genetic correlation matrix and its sampling covariance.
///
/// Port of R GenomicSEM's `rgmodel()`.
///
/// If `model` is `Some(model_str)`, fit the user-specified model.
/// If `model` is `None`, fit a common factor model (original behavior).
pub fn run_rgmodel(
    s: &Mat<f64>,
    v: &Mat<f64>,
    estimation: &str,
    model: Option<&str>,
    std_lv: bool,
) -> Result<RgModelResult> {
    // Fit either user-specified or common factor model
    let sem_result = if let Some(model_str) = model.filter(|m| !m.is_empty()) {
        fit_user_model(s, v, estimation, model_str, std_lv)?
    } else {
        crate::commonfactor::run_commonfactor(s, v, estimation)?
    };

    let k = s.nrows();
    let sigma = &sem_result.implied_cov;

    // Compute R = D^{-1/2} * Sigma * D^{-1/2} where D = diag(Sigma)
    let sds: Vec<f64> = (0..k).map(|i| sigma[(i, i)].sqrt().max(1e-10)).collect();
    let r = Mat::from_fn(k, k, |i, j| sigma[(i, j)] / (sds[i] * sds[j]));

    // Compute V_R via delta method
    // The Jacobian of the standardization: d(r_ij)/d(sigma_ab)
    // Use numerical differentiation on the vech elements
    let kstar = k * (k + 1) / 2;
    let sigma_vec = vech::vech(sigma);

    // Numerical Jacobian: perturb each element of vech(Sigma), compute change in vech(R)
    let eps = 1e-7;
    let r_vec = vech::vech(&r);
    let mut jacobian = Mat::zeros(kstar, kstar); // d(vech(R)) / d(vech(Sigma))

    for col in 0..kstar {
        let mut perturbed = sigma_vec.clone();
        perturbed[col] += eps;
        let sigma_plus = vech::vech_reverse(&perturbed, k);
        let sds_plus: Vec<f64> = (0..k)
            .map(|i| sigma_plus[(i, i)].sqrt().max(1e-10))
            .collect();
        let r_plus = Mat::from_fn(k, k, |i, j| {
            sigma_plus[(i, j)] / (sds_plus[i] * sds_plus[j])
        });
        let r_plus_vec = vech::vech(&r_plus);

        for row in 0..kstar {
            jacobian[(row, col)] = (r_plus_vec[row] - r_vec[row]) / eps;
        }
    }

    // V_R = J * V * J'
    let jv = &jacobian * v;
    let v_r = &jv * jacobian.transpose();

    Ok(RgModelResult { r, v_r, sem_result })
}

/// Fit a user-specified SEM model (same logic as commonfactor but with arbitrary model syntax).
fn fit_user_model(
    s: &Mat<f64>,
    v: &Mat<f64>,
    estimation: &str,
    model_str: &str,
    std_lv: bool,
) -> Result<crate::SemResult> {
    use crate::ParamEstimate;
    use crate::estimator;
    use crate::fit_indices;
    use crate::model::Model;
    use crate::sandwich;
    use crate::syntax::parse_model;
    use statrs::distribution::{ChiSquared, ContinuousCDF};

    let k = s.nrows();
    let obs_names: Vec<String> = (0..k).map(|i| format!("V{}", i + 1)).collect();

    // Parse user model
    let pt = parse_model(model_str, std_lv).map_err(|e| anyhow::anyhow!("{e}"))?;
    let mut model = Model::from_partable(&pt, &obs_names);

    // Fit
    let kstar = k * (k + 1) / 2;
    let v_diag: Vec<f64> = (0..kstar).map(|i| v[(i, i)]).collect();

    let fit = if estimation.to_uppercase() == "ML" {
        estimator::fit_ml(&mut model, s, 1000)
    } else {
        estimator::fit_dwls(&mut model, s, &v_diag, 1000)
    };

    if !fit.converged {
        log::warn!("User model did not converge");
    }

    // Sandwich SEs
    let w = Mat::from_fn(kstar, kstar, |i, j| {
        if i == j {
            if v_diag[i] > 1e-30 {
                1.0 / v_diag[i]
            } else {
                0.0
            }
        } else {
            0.0
        }
    });
    let (se_vec, _ohtt) = sandwich::sandwich_se(&mut model, &w, v);

    // Build parameter estimates
    let mut parameters = Vec::new();
    let mut free_idx = 0;
    for row in &pt.rows {
        if row.free > 0 {
            let est = fit.params.get(free_idx).copied().unwrap_or(0.0);
            let se = se_vec.get(free_idx).copied().unwrap_or(0.0);
            let z = if se > 0.0 { est / se } else { 0.0 };
            let p = if z.abs() > 0.0 {
                let chi2 = ChiSquared::new(1.0).unwrap();
                2.0 * (1.0 - chi2.cdf(z * z))
            } else {
                1.0
            };
            parameters.push(ParamEstimate {
                lhs: row.lhs.clone(),
                op: row.op.to_string(),
                rhs: row.rhs.clone(),
                est,
                se,
                z,
                p,
            });
            free_idx += 1;
        }
    }

    // Fit indices
    let null_model_str: String = obs_names
        .iter()
        .map(|v| format!("{v} ~~ {v}"))
        .collect::<Vec<_>>()
        .join("\n");
    let null_pt = parse_model(&null_model_str, false).map_err(|e| anyhow::anyhow!("{e}"))?;
    let mut null_model = Model::from_partable(&null_pt, &obs_names);
    let null_fit = if estimation.to_uppercase() == "ML" {
        estimator::fit_ml(&mut null_model, s, 1000)
    } else {
        estimator::fit_dwls(&mut null_model, s, &v_diag, 1000)
    };
    let _ = null_fit;
    let null_sigma = null_model.implied_cov();
    let null_df = kstar - k;
    let null_fit_stats = fit_indices::compute_fit(s, &null_sigma, v, null_df, k, None, None);

    let sigma_hat = model.implied_cov();
    let n_free = model.n_free();
    let df = kstar.saturating_sub(n_free);
    let model_fit = fit_indices::compute_fit(
        s,
        &sigma_hat,
        v,
        df,
        n_free,
        Some(null_fit_stats.chisq),
        Some(null_df),
    );

    Ok(crate::SemResult {
        parameters,
        fit: model_fit,
        implied_cov: sigma_hat,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rgmodel_diagonal_ones() {
        // For a 2x2 case, the diagonal of R should be 1.0
        let s = faer::mat![[0.5, 0.1], [0.1, 0.3]];
        let kstar = 3;
        let v = Mat::<f64>::identity(kstar, kstar) * 0.001;

        let result = run_rgmodel(&s, &v, "DWLS", None, false).unwrap();
        assert!((result.r[(0, 0)] - 1.0).abs() < 1e-6);
        assert!((result.r[(1, 1)] - 1.0).abs() < 1e-6);
        assert_eq!(result.v_r.nrows(), kstar);
        assert_eq!(result.v_r.ncols(), kstar);
    }

    #[test]
    fn test_rgmodel_with_user_model() {
        // 3x3 case with explicit 1-factor model
        let s = faer::mat![[0.60, 0.42, 0.35], [0.42, 0.50, 0.30], [0.35, 0.30, 0.40]];
        let kstar = 6;
        let v = Mat::<f64>::identity(kstar, kstar) * 0.001;

        let model = "F1 =~ NA*V1 + V2 + V3\nF1 ~~ 1*F1\nV1 ~~ V1\nV2 ~~ V2\nV3 ~~ V3";
        let result = run_rgmodel(&s, &v, "DWLS", Some(model), false).unwrap();
        assert!((result.r[(0, 0)] - 1.0).abs() < 1e-6);
        assert!((result.r[(1, 1)] - 1.0).abs() < 1e-6);
        assert!((result.r[(2, 2)] - 1.0).abs() < 1e-6);
    }
}
