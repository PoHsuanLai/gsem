//! Common factor model fitting.
//!
//! Auto-generates and fits a 1-factor CFA model:
//!   F1 =~ NA*V1 + V2 + ... + Vk
//!   F1 ~~ 1*F1
//!   V1 ~~ V1; V2 ~~ V2; ...; Vk ~~ Vk

use anyhow::Result;
use faer::Mat;
use statrs::distribution::{ChiSquared, ContinuousCDF};

use crate::estimator;
use crate::fit_indices;
use crate::model::Model;
use crate::sandwich;
use crate::syntax::parse_model;
use crate::{EstimationMethod, ParamEstimate, SemResult};

/// Run common factor analysis on LDSC output.
///
/// Equivalent to R GenomicSEM's `commonfactor()`.
pub fn run_commonfactor(s: &Mat<f64>, v: &Mat<f64>, estimation: EstimationMethod) -> Result<SemResult> {
    let k = s.nrows();
    let obs_names: Vec<String> = (0..k).map(|i| format!("V{}", i + 1)).collect();

    // Auto-generate model: F1 =~ NA*V1 + V2 + ... + Vk; F1 ~~ 1*F1; residual variances
    let loading = std::iter::once(format!("NA*{}", obs_names[0]))
        .chain(obs_names[1..].iter().cloned())
        .collect::<Vec<_>>()
        .join(" + ");
    let residuals: String = obs_names
        .iter()
        .map(|v| format!("{v} ~~ {v}"))
        .collect::<Vec<_>>()
        .join("\n");
    let model_str = format!("F1 =~ {loading}\nF1 ~~ 1*F1\n{residuals}");

    // Parse and build model
    let pt = parse_model(&model_str, false).map_err(|e| anyhow::anyhow!("{e}"))?;
    let mut model = Model::from_partable(&pt, &obs_names);

    // Fit
    let kstar = k * (k + 1) / 2;
    let v_diag: Vec<f64> = (0..kstar).map(|i| v[(i, i)]).collect();

    let fit = match estimation {
        EstimationMethod::Ml => estimator::fit_ml(&mut model, s, 1000, None),
        EstimationMethod::Dwls => estimator::fit_dwls(&mut model, s, &v_diag, 1000, None),
    };

    if !fit.converged {
        log::warn!("Common factor model did not converge");
    }

    // Compute implied cov BEFORE sandwich_se (which may mutate model internals)
    let sigma_hat = model.implied_cov();

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
                ChiSquared::new(1.0)
                    .map(|chi2| 2.0 * (1.0 - chi2.cdf(z * z)))
                    .unwrap_or(1.0)
            } else {
                1.0
            };
            parameters.push(ParamEstimate {
                lhs: row.lhs.clone(),
                op: row.op,
                rhs: row.rhs.clone(),
                est,
                se,
                z,
                p,
            });
            free_idx += 1;
        }
    }

    // Fit indices: need null model for CFI
    // Null model: all variances free, all covariances fixed to 0
    let null_model_str: String = obs_names
        .iter()
        .map(|v| format!("{v} ~~ {v}"))
        .collect::<Vec<_>>()
        .join("\n");
    let null_pt = parse_model(&null_model_str, false).map_err(|e| anyhow::anyhow!("{e}"))?;
    let mut null_model = Model::from_partable(&null_pt, &obs_names);
    let null_fit = match estimation {
        EstimationMethod::Ml => estimator::fit_ml(&mut null_model, s, 1000, None),
        EstimationMethod::Dwls => estimator::fit_dwls(&mut null_model, s, &v_diag, 1000, None),
    };
    let _ = null_fit; // used only to set parameter values in null_model
    let null_sigma = null_model.implied_cov();
    let null_df = kstar - k; // null model has k free params (variances only)
    let null_fit_stats = fit_indices::compute_fit(s, &null_sigma, v, null_df, k, None, None);

    // Model fit
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

    Ok(SemResult {
        parameters,
        fit: model_fit,
        implied_cov: sigma_hat,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use faer::Mat;

    fn simple_s_and_v() -> (Mat<f64>, Mat<f64>) {
        let s = faer::mat![[0.60, 0.42, 0.35], [0.42, 0.50, 0.30], [0.35, 0.30, 0.40],];
        // V is kstar x kstar where kstar = k*(k+1)/2 = 6
        let v = Mat::from_fn(6, 6, |i, j| if i == j { 0.001 } else { 0.0 });
        (s, v)
    }

    #[test]
    fn test_commonfactor_dwls_converges() {
        let (s, v) = simple_s_and_v();
        let result = run_commonfactor(&s, &v, crate::EstimationMethod::Dwls).expect("DWLS should not error");
        // Should have parameters: 3 loadings + 3 residual variances = 6
        assert_eq!(result.parameters.len(), 6, "expected 6 free parameters");
        // Check that estimates are finite
        for p in &result.parameters {
            assert!(
                p.est.is_finite(),
                "estimate for {} {} {} should be finite",
                p.lhs,
                p.op,
                p.rhs
            );
            assert!(
                p.se.is_finite(),
                "SE for {} {} {} should be finite",
                p.lhs,
                p.op,
                p.rhs
            );
        }
        // Fit indices should be finite
        assert!(result.fit.chisq.is_finite(), "chi-square should be finite");
        assert!(result.fit.cfi.is_finite(), "CFI should be finite");
    }

    #[test]
    fn test_commonfactor_ml_converges() {
        let (s, v) = simple_s_and_v();
        let result = run_commonfactor(&s, &v, crate::EstimationMethod::Ml).expect("ML should not error");
        assert_eq!(result.parameters.len(), 6, "expected 6 free parameters");
        for p in &result.parameters {
            assert!(p.est.is_finite(), "estimate should be finite");
        }
    }

    #[test]
    fn test_commonfactor_parameter_structure() {
        let (s, v) = simple_s_and_v();
        let result = run_commonfactor(&s, &v, crate::EstimationMethod::Dwls).unwrap();

        // First 3 params should be loadings (F1 =~ V1, V2, V3)
        let loadings: Vec<_> = result.parameters.iter().filter(|p| p.op == crate::syntax::Op::Loading).collect();
        assert_eq!(loadings.len(), 3, "should have 3 loadings");

        // Last 3 params should be residual variances (V1 ~~ V1, etc.)
        let resid: Vec<_> = result.parameters.iter().filter(|p| p.op == crate::syntax::Op::Covariance).collect();
        assert_eq!(resid.len(), 3, "should have 3 residual variances");

        // Loadings should be positive for this well-behaved covariance matrix
        for l in &loadings {
            assert!(
                l.est > 0.0,
                "loading {} should be positive, got {}",
                l.rhs,
                l.est
            );
        }

        // Residual variances should be positive
        for r in &resid {
            assert!(
                r.est > 0.0,
                "residual variance {} should be positive, got {}",
                r.rhs,
                r.est
            );
        }
    }

    #[test]
    fn test_commonfactor_implied_cov_shape() {
        let (s, v) = simple_s_and_v();
        let result = run_commonfactor(&s, &v, crate::EstimationMethod::Dwls).unwrap();
        assert_eq!(result.implied_cov.nrows(), 3);
        assert_eq!(result.implied_cov.ncols(), 3);
    }
}
