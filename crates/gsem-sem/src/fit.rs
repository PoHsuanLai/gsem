//! Fit a parsed model (`ParTable`) to an S/V covariance structure via
//! DWLS/ML with sandwich standard errors, returning parameter estimates.
//!
//! This is the shared fitting path used by `usermodel`, `rgmodel`, and
//! `enrich` (model-based functional enrichment).

use faer::Mat;
use statrs::distribution::{ChiSquared, ContinuousCDF};

use crate::syntax::ParTable;
use crate::{ParamEstimate, SemResult, estimator, fit_indices, model::Model, sandwich};

/// Fit a parsed parameter table to `s` (genetic covariance) with sampling
/// covariance `v`, returning the parameter estimates with sandwich SEs.
pub fn fit_partable(
    pt: &ParTable,
    obs_names: &[String],
    s: &Mat<f64>,
    v: &Mat<f64>,
    estimation: crate::EstimationMethod,
) -> SemResult {
    let k = s.nrows();
    let kstar = k * (k + 1) / 2;
    let mut model = Model::from_partable(pt, obs_names);

    let v_diag: Vec<f64> = (0..kstar).map(|i| v[(i, i)]).collect();
    let fit = match estimation {
        crate::EstimationMethod::Ml => estimator::fit_ml(&mut model, s, 1000, None),
        crate::EstimationMethod::Dwls => estimator::fit_dwls(&mut model, s, &v_diag, 1000, None),
    };

    let w = Mat::from_fn(kstar, kstar, |i, j| {
        if i == j && v_diag[i] > 1e-30 {
            1.0 / v_diag[i]
        } else {
            0.0
        }
    });
    let (se_vec, _) = sandwich::sandwich_se(&mut model, &w, v);

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

    let sigma_hat = model.implied_cov();
    let n_free = model.n_free();
    let df = kstar.saturating_sub(n_free);
    let model_fit = fit_indices::compute_fit(s, &sigma_hat, v, df, n_free, None, None);

    SemResult {
        parameters,
        fit: model_fit,
        implied_cov: sigma_hat,
    }
}
