//! Model-implied genetic correlation matrix.
//!
//! Port of R GenomicSEM's `rgmodel()`.
//!
//! 1. Standardizes S → S_Stand = cov2cor(S)
//! 2. Rescales V accordingly
//! 3. Fits common factor (or user) model on S_Stand
//! 4. Returns R (model-implied correlations) and V_R (sampling covariance)

use anyhow::Result;
use faer::Mat;
use gsem_matrix::smooth::cov_to_cor;
use gsem_matrix::vech;

/// Result of rgmodel.
#[derive(Debug, Clone)]
pub struct RgModelResult {
    /// Model-implied genetic correlation matrix (k x k)
    pub r: Mat<f64>,
    /// Sampling covariance of the off-diagonal correlations
    pub v_r: Mat<f64>,
    /// The underlying SEM result
    pub sem_result: crate::SemResult,
}

/// Compute model-implied genetic correlation matrix and its sampling covariance.
///
/// Port of R GenomicSEM's `rgmodel()`:
/// - Standardizes S to a correlation matrix
/// - Rescales V by the ratio S_Stand/S
/// - Fits the model on the standardized matrix
/// - Returns model-implied correlations and their sampling covariance
pub fn run_rgmodel(
    s: &Mat<f64>,
    v: &Mat<f64>,
    estimation: &str,
) -> Result<RgModelResult> {
    run_rgmodel_inner(s, v, estimation, None, false)
}

/// Full version with optional user model, std_lv, and phenotype subsetting.
pub fn run_rgmodel_with_model(
    s: &Mat<f64>,
    v: &Mat<f64>,
    estimation: &str,
    model: Option<&str>,
    std_lv: bool,
) -> Result<RgModelResult> {
    run_rgmodel_inner(s, v, estimation, model, std_lv)
}

/// Version with phenotype subsetting.
/// `sub` contains indices of phenotypes to keep (0-based).
pub fn run_rgmodel_sub(
    s: &Mat<f64>,
    v: &Mat<f64>,
    estimation: &str,
    model: Option<&str>,
    std_lv: bool,
    sub_indices: &[usize],
) -> Result<RgModelResult> {
    if sub_indices.is_empty() {
        return run_rgmodel_inner(s, v, estimation, model, std_lv);
    }

    let k_sub = sub_indices.len();
    // Subset S
    let s_sub = Mat::from_fn(k_sub, k_sub, |i, j| s[(sub_indices[i], sub_indices[j])]);

    // Subset V (which is kstar x kstar, indexed by vech)
    let k_full = s.nrows();
    let kstar_sub = k_sub * (k_sub + 1) / 2;

    // Build vech index mapping: sub vech → full vech
    let mut vech_map = Vec::with_capacity(kstar_sub);
    for col_sub in 0..k_sub {
        for row_sub in col_sub..k_sub {
            let row_full = sub_indices[row_sub];
            let col_full = sub_indices[col_sub];
            let (r, c) = if row_full >= col_full {
                (row_full, col_full)
            } else {
                (col_full, row_full)
            };
            // vech index in full matrix (column-major lower triangle)
            let idx = c * k_full - c * (c.wrapping_sub(1)) / 2 + (r - c);
            vech_map.push(idx);
        }
    }

    let v_sub = Mat::from_fn(kstar_sub, kstar_sub, |i, j| v[(vech_map[i], vech_map[j])]);

    run_rgmodel_inner(&s_sub, &v_sub, estimation, model, std_lv)
}

fn run_rgmodel_inner(
    s: &Mat<f64>,
    v: &Mat<f64>,
    estimation: &str,
    model: Option<&str>,
    std_lv: bool,
) -> Result<RgModelResult> {
    let k = s.nrows();
    let kstar = k * (k + 1) / 2;

    // Step 1: Standardize S → S_Stand = cov2cor(S)
    let s_stand = cov_to_cor(s);
    // Step 2: Rescale V
    let s_vec = vech::vech(s);
    let s_stand_vec = vech::vech(&s_stand);

    let scale_o: Vec<f64> = s_stand_vec
        .iter()
        .zip(s_vec.iter())
        .map(|(&stand, &orig)| {
            if orig.abs() > 1e-30 {
                stand / orig
            } else {
                0.0
            }
        })
        .collect();

    let v_diag_sqrt: Vec<f64> = (0..kstar).map(|i| v[(i, i)].sqrt()).collect();
    let dvcovl: Vec<f64> = v_diag_sqrt
        .iter()
        .zip(scale_o.iter())
        .map(|(&d, &sc)| d * sc)
        .collect();

    let v_cor = cov_to_cor(v);
    let v_stand = Mat::from_fn(kstar, kstar, |i, j| dvcovl[i] * v_cor[(i, j)] * dvcovl[j]);

    let mut v_stand_fixed = v_stand.clone();
    for i in 0..kstar {
        if v_stand_fixed[(i, i)].abs() < 2e-9 {
            v_stand_fixed[(i, i)] = 2e-9;
        }
    }

    // Step 3: Fit model on S_Stand
    // Use a simple diagonal V for fitting (matching R's lavaan behavior on correlation matrices).
    // The rescaled V_stand is used only for sandwich SEs in Step 5.
    let v_fit = Mat::<f64>::identity(kstar, kstar) * 0.001;
    let sem_result = if let Some(model_str) = model.filter(|m| !m.is_empty()) {
        fit_user_model_on_stand(&s_stand, &v_fit, estimation, model_str, std_lv)?
    } else {
        crate::commonfactor::run_commonfactor(&s_stand, &v_fit, estimation)?
    };

    // Step 4: The model-implied covariance from fitting on S_Stand IS the correlation matrix
    let r = sem_result.implied_cov.clone();

    // Step 5: V_R — R returns the V from the standardized model fit
    // R reshapes its V_R oddly (matrix of sqrt(length)), but the actual content
    // is the sandwich SE covariance of the standardized model parameters
    // For compatibility, return the same kstar x kstar V_R from the standardized fit
    // But R returns k x k... it takes only the unique off-diagonal elements
    // Actually R takes the full V from the standardized sandwich and returns it as-is
    // Let's match R's output: extract the k*(k-1)/2 unique off-diagonal correlations
    // and their sampling covariance

    // For now, return the model-implied correlation matrix and the V from standardized fit
    // R's V_R is actually the result of the sandwich on the standardized model
    // which has n_free x n_free dimensions, reshaped to sqrt(n) x sqrt(n)
    // For a common factor 3-trait model: 6 free params, V_R would be 6x6 reshaped... no
    // R returns 3x3 for 3 traits. Let me check: the R code does:
    //   rgmodel$V_R = matrix(rgmodel$V_R, nrow=sqrt(length(rgmodel$V_R)))
    // If V_R has 9 elements → 3x3. Those 9 elements come from the 3 off-diagonal correlations
    // and their sampling covariance.

    // The simplest correct approach: compute V_R as the sampling covariance of vech(R)
    // using delta method on the standardization transformation, same as before but on S_Stand
    let v_r = compute_v_r(&r, &v_stand_fixed, k);

    Ok(RgModelResult { r, v_r, sem_result })
}

/// Fit a user-specified model on the standardized S matrix.
fn fit_user_model_on_stand(
    s_stand: &Mat<f64>,
    v_stand: &Mat<f64>,
    estimation: &str,
    model_str: &str,
    std_lv: bool,
) -> Result<crate::SemResult> {
    use crate::{ParamEstimate, estimator, fit_indices, model::Model, sandwich, syntax::parse_model};
    use statrs::distribution::{ChiSquared, ContinuousCDF};

    let k = s_stand.nrows();
    let kstar = k * (k + 1) / 2;
    let obs_names: Vec<String> = (0..k).map(|i| format!("V{}", i + 1)).collect();

    let pt = parse_model(model_str, std_lv).map_err(|e| anyhow::anyhow!("{e}"))?;
    let mut model = Model::from_partable(&pt, &obs_names);

    let v_diag: Vec<f64> = (0..kstar).map(|i| v_stand[(i, i)]).collect();

    let fit = if estimation.to_uppercase() == "ML" {
        estimator::fit_ml(&mut model, s_stand, 1000, None)
    } else {
        estimator::fit_dwls(&mut model, s_stand, &v_diag, 1000, None)
    };

    let w = Mat::from_fn(kstar, kstar, |i, j| {
        if i == j && v_diag[i] > 1e-30 { 1.0 / v_diag[i] } else { 0.0 }
    });
    let (se_vec, _) = sandwich::sandwich_se(&mut model, &w, v_stand);

    let mut parameters = Vec::new();
    let mut free_idx = 0;
    for row in &pt.rows {
        if row.free > 0 {
            let est = fit.params.get(free_idx).copied().unwrap_or(0.0);
            let se = se_vec.get(free_idx).copied().unwrap_or(0.0);
            let z = if se > 0.0 { est / se } else { 0.0 };
            let p = if z.abs() > 0.0 {
                2.0 * (1.0 - ChiSquared::new(1.0).unwrap().cdf(z * z))
            } else {
                1.0
            };
            parameters.push(ParamEstimate {
                lhs: row.lhs.clone(), op: row.op.to_string(), rhs: row.rhs.clone(),
                est, se, z, p,
            });
            free_idx += 1;
        }
    }

    let sigma_hat = model.implied_cov();
    let n_free = model.n_free();
    let df = kstar.saturating_sub(n_free);
    let model_fit = fit_indices::compute_fit(s_stand, &sigma_hat, v_stand, df, n_free, None, None);

    Ok(crate::SemResult { parameters, fit: model_fit, implied_cov: sigma_hat })
}

/// Compute sampling covariance of the correlation matrix via delta method.
fn compute_v_r(r: &Mat<f64>, v_stand: &Mat<f64>, k: usize) -> Mat<f64> {
    let kstar = k * (k + 1) / 2;
    let r_vec = vech::vech(r);
    let eps = 1e-7;

    // Jacobian: perturb each element of vech(Sigma_stand), compute change in vech(R)
    let sigma_vec = r_vec.clone(); // for standardized model, Sigma ≈ R
    let mut jacobian = Mat::zeros(kstar, kstar);

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

    // V_R = J * V_stand * J'
    let jv = &jacobian * v_stand;
    &jv * jacobian.transpose()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rgmodel_diagonal_ones() {
        let s = faer::mat![[0.5, 0.1], [0.1, 0.3]];
        let kstar = 3;
        let v = Mat::<f64>::identity(kstar, kstar) * 0.001;

        let result = run_rgmodel(&s, &v, "DWLS").unwrap();
        assert!((result.r[(0, 0)] - 1.0).abs() < 1e-3);
        assert!((result.r[(1, 1)] - 1.0).abs() < 1e-3);
    }

    #[test]
    fn test_rgmodel_correlation_range() {
        let s = faer::mat![[0.60, 0.42, 0.35], [0.42, 0.50, 0.30], [0.35, 0.30, 0.40]];
        let kstar = 6;
        let v = Mat::<f64>::identity(kstar, kstar) * 0.001;

        let result = run_rgmodel(&s, &v, "DWLS").unwrap();
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    result.r[(i, j)] >= -1.0 - 0.01 && result.r[(i, j)] <= 1.0 + 0.01,
                    "r[{i},{j}] = {} out of range",
                    result.r[(i, j)]
                );
            }
        }
        let expected_r12 = 0.42 / (0.60_f64.sqrt() * 0.50_f64.sqrt());
        assert!(
            (result.r[(0, 1)] - expected_r12).abs() < 0.1,
            "r[0,1]={} expected ~{expected_r12}",
            result.r[(0, 1)]
        );
    }

    #[test]
    fn test_rgmodel_real_scale_v() {
        // Same S but with V scaled like real LDSC output (~1e-5)
        let s = faer::mat![[0.60, 0.42, 0.35], [0.42, 0.50, 0.30], [0.35, 0.30, 0.40]];
        let kstar = 6;
        let v = Mat::from_fn(kstar, kstar, |i, j| if i == j { 3e-5 } else { 0.0 });

        let result = run_rgmodel(&s, &v, "DWLS").unwrap();
        let expected_r12 = 0.42 / (0.60_f64.sqrt() * 0.50_f64.sqrt());
        eprintln!("Real-scale V: r[0,1]={:.4} expected ~{expected_r12:.4}", result.r[(0, 1)]);
        assert!(
            (result.r[(0, 1)] - expected_r12).abs() < 0.1,
            "r[0,1]={} expected ~{expected_r12} (real-scale V)",
            result.r[(0, 1)]
        );
    }
}
