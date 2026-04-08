use faer::Mat;
use faer::prelude::Solve;

use gsem_matrix::vech;

use crate::model::Model;

/// Result of SEM estimation.
#[derive(Debug, Clone)]
pub struct FitResult {
    /// Optimized parameter values
    pub params: Vec<f64>,
    /// Final objective function value
    pub objective: f64,
    /// Did the optimizer converge?
    pub converged: bool,
}

/// Fit a model using Diagonally Weighted Least Squares (DWLS).
///
/// Minimizes: F = (s - sigma(theta))' W (s - sigma(theta))
/// where s = vech(S), sigma(theta) = vech(implied_cov), W = diag(V)^{-1}
pub fn fit_dwls(model: &mut Model, s_obs: &Mat<f64>, v_diag: &[f64], max_iter: usize) -> FitResult {
    let s_vec = vech::vech(s_obs);
    // Weight matrix: W = 1/diag(V)
    let w: Vec<f64> = v_diag
        .iter()
        .map(|&v| if v > 1e-30 { 1.0 / v } else { 0.0 })
        .collect();

    let mut params = model.get_param_vec();

    // Simple gradient descent with line search
    let mut best_obj = dwls_objective(model, &s_vec, &w);
    let mut converged = false;

    for _iter in 0..max_iter {
        // Compute gradient via central differences
        let grad = numerical_gradient(model, &params, &s_vec, &w);

        // Line search
        let mut step = 0.1;
        let mut improved = false;

        for _ in 0..20 {
            let trial: Vec<f64> = params
                .iter()
                .zip(grad.iter())
                .enumerate()
                .map(|(i, (&p, &g))| {
                    let new_val = p - step * g;
                    // Apply lower bounds
                    if let Some(lb) = model.lower_bounds[i] {
                        new_val.max(lb)
                    } else {
                        new_val
                    }
                })
                .collect();

            model.set_param_vec(&trial);
            let obj = dwls_objective(model, &s_vec, &w);

            if obj < best_obj {
                params = trial;
                best_obj = obj;
                improved = true;
                break;
            }
            step *= 0.5;
        }

        if !improved {
            // Try smaller step sizes
            step = 0.001;
            let trial: Vec<f64> = params
                .iter()
                .zip(grad.iter())
                .enumerate()
                .map(|(i, (&p, &g))| {
                    let new_val = p - step * g;
                    if let Some(lb) = model.lower_bounds[i] {
                        new_val.max(lb)
                    } else {
                        new_val
                    }
                })
                .collect();

            model.set_param_vec(&trial);
            let obj = dwls_objective(model, &s_vec, &w);
            if obj < best_obj {
                params = trial;
                best_obj = obj;
            }
        }

        // Check convergence
        let grad_norm: f64 = grad.iter().map(|g| g * g).sum::<f64>().sqrt();
        if grad_norm < 1e-6 {
            converged = true;
            break;
        }
    }

    model.set_param_vec(&params);

    FitResult {
        params,
        objective: best_obj,
        converged,
    }
}

/// Fit a model using Maximum Likelihood.
///
/// Minimizes: F = log|Sigma| + tr(S * Sigma^{-1}) - log|S| - p
pub fn fit_ml(model: &mut Model, s_obs: &Mat<f64>, max_iter: usize) -> FitResult {
    let _s_vec = vech::vech(s_obs);
    let mut params = model.get_param_vec();

    let mut best_obj = ml_objective(model, s_obs);
    let mut converged = false;

    for _ in 0..max_iter {
        let grad = numerical_gradient_ml(model, &params, s_obs);

        let mut step = 0.05;
        let mut improved = false;

        for _ in 0..20 {
            let trial: Vec<f64> = params
                .iter()
                .zip(grad.iter())
                .enumerate()
                .map(|(i, (&p_val, &g))| {
                    let new_val = p_val - step * g;
                    if let Some(lb) = model.lower_bounds[i] {
                        new_val.max(lb)
                    } else {
                        new_val
                    }
                })
                .collect();

            model.set_param_vec(&trial);
            let obj = ml_objective(model, s_obs);

            if obj < best_obj && obj.is_finite() {
                params = trial;
                best_obj = obj;
                improved = true;
                break;
            }
            step *= 0.5;
        }

        if !improved {
            break;
        }

        let grad_norm: f64 = grad.iter().map(|g| g * g).sum::<f64>().sqrt();
        if grad_norm < 1e-6 {
            converged = true;
            break;
        }
    }

    model.set_param_vec(&params);
    FitResult {
        params,
        objective: best_obj,
        converged,
    }
}

/// DWLS objective function value.
fn dwls_objective(model: &Model, s_vec: &[f64], w: &[f64]) -> f64 {
    let sigma = model.implied_cov();
    let sigma_vec = vech::vech(&sigma);
    let mut obj = 0.0;
    for i in 0..s_vec.len() {
        let r = s_vec[i] - sigma_vec[i];
        obj += w[i] * r * r;
    }
    obj
}

/// ML objective function value.
fn ml_objective(model: &Model, s_obs: &Mat<f64>) -> f64 {
    let sigma = model.implied_cov();
    let p = s_obs.nrows();

    // log|Sigma|
    let chol = sigma.llt(faer::Side::Lower);
    let log_det_sigma = match chol {
        Ok(llt) => {
            let l = llt.L();
            let mut ld = 0.0;
            for i in 0..p {
                ld += l[(i, i)].ln();
            }
            2.0 * ld
        }
        Err(_) => return f64::INFINITY,
    };

    // tr(S * Sigma^{-1})
    let sigma_inv = match sigma.llt(faer::Side::Lower) {
        Ok(llt) => {
            let mut eye = Mat::identity(p, p);
            llt.solve_in_place(eye.as_mut());
            eye
        }
        Err(_) => return f64::INFINITY,
    };

    let s_sigma_inv = s_obs * &sigma_inv;
    let trace = (0..p).map(|i| s_sigma_inv[(i, i)]).sum::<f64>();

    // log|S|
    let log_det_s = match s_obs.llt(faer::Side::Lower) {
        Ok(llt) => {
            let l = llt.L();
            let mut ld = 0.0;
            for i in 0..p {
                ld += l[(i, i)].ln();
            }
            2.0 * ld
        }
        Err(_) => 0.0,
    };

    log_det_sigma + trace - log_det_s - p as f64
}

/// Numerical gradient for DWLS objective via central differences.
fn numerical_gradient(model: &mut Model, params: &[f64], s_vec: &[f64], w: &[f64]) -> Vec<f64> {
    let eps = 1e-7;
    let n = params.len();
    let mut grad = vec![0.0; n];

    for i in 0..n {
        let mut p_plus = params.to_vec();
        let mut p_minus = params.to_vec();
        p_plus[i] += eps;
        p_minus[i] -= eps;

        model.set_param_vec(&p_plus);
        let f_plus = dwls_objective(model, s_vec, w);

        model.set_param_vec(&p_minus);
        let f_minus = dwls_objective(model, s_vec, w);

        grad[i] = (f_plus - f_minus) / (2.0 * eps);
    }

    model.set_param_vec(params);
    grad
}

/// Numerical gradient for ML objective.
fn numerical_gradient_ml(model: &mut Model, params: &[f64], s_obs: &Mat<f64>) -> Vec<f64> {
    let eps = 1e-7;
    let n = params.len();
    let mut grad = vec![0.0; n];

    for i in 0..n {
        let mut p_plus = params.to_vec();
        let mut p_minus = params.to_vec();
        p_plus[i] += eps;
        p_minus[i] -= eps;

        model.set_param_vec(&p_plus);
        let f_plus = ml_objective(model, s_obs);

        model.set_param_vec(&p_minus);
        let f_minus = ml_objective(model, s_obs);

        grad[i] = (f_plus - f_minus) / (2.0 * eps);
    }

    model.set_param_vec(params);
    grad
}
