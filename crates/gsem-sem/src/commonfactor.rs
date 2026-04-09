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
use crate::{ParamEstimate, SemResult};

/// Run common factor analysis on LDSC output.
///
/// Equivalent to R GenomicSEM's `commonfactor()`.
pub fn run_commonfactor(s: &Mat<f64>, v: &Mat<f64>, estimation: &str) -> Result<SemResult> {
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

    let fit = if estimation.to_uppercase() == "ML" {
        estimator::fit_ml(&mut model, s, 1000)
    } else {
        estimator::fit_dwls(&mut model, s, &v_diag, 1000)
    };

    if !fit.converged {
        log::warn!("Common factor model did not converge");
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

    // Fit indices: need null model for CFI
    // Null model: all variances free, all covariances fixed to 0
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
    let _ = null_fit; // used only to set parameter values in null_model
    let null_sigma = null_model.implied_cov();
    let null_df = kstar - k; // null model has k free params (variances only)
    let null_fit_stats = fit_indices::compute_fit(s, &null_sigma, v, null_df, k, None, None);

    // Model fit
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

    Ok(SemResult {
        parameters,
        fit: model_fit,
        implied_cov: sigma_hat,
    })
}
