use std::collections::VecDeque;

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
/// Fit a model using Diagonally Weighted Least Squares (DWLS).
///
/// Minimizes: F = (s - sigma(theta))' W (s - sigma(theta))
/// where s = vech(S), sigma(theta) = vech(implied_cov), W = diag(V)^{-1}
///
/// `grad_tol`: gradient norm convergence threshold (default 1e-6).
pub fn fit_dwls(
    model: &mut Model,
    s_obs: &Mat<f64>,
    v_diag: &[f64],
    max_iter: usize,
    grad_tol: Option<f64>,
) -> FitResult {
    let s_vec = vech::vech(s_obs).expect("s_obs must be square");
    let w: Vec<f64> = v_diag
        .iter()
        .map(|&v| if v > 1e-30 { 1.0 / v } else { 0.0 })
        .collect();

    let lower_bounds: Vec<Option<f64>> = model.lower_bounds.clone();

    let obj_fn = |model: &Model| -> f64 { dwls_objective(model, &s_vec, &w) };

    let grad_fn = |model: &mut Model, params: &[f64]| -> Vec<f64> {
        numerical_gradient(model, params, |m| dwls_objective(m, &s_vec, &w))
    };

    lbfgs_minimize(model, max_iter, &lower_bounds, obj_fn, grad_fn, grad_tol)
}

/// Fit a model using Maximum Likelihood.
///
/// Minimizes: F = log|Sigma| + tr(S * Sigma^{-1}) - log|S| - p
///
/// `grad_tol`: gradient norm convergence threshold (default 1e-6).
pub fn fit_ml(
    model: &mut Model,
    s_obs: &Mat<f64>,
    max_iter: usize,
    grad_tol: Option<f64>,
) -> FitResult {
    let lower_bounds: Vec<Option<f64>> = model.lower_bounds.clone();
    let s_owned = s_obs.to_owned();

    let obj_fn = |model: &Model| -> f64 { ml_objective(model, &s_owned) };

    let grad_fn = |model: &mut Model, params: &[f64]| -> Vec<f64> {
        numerical_gradient(model, params, |m| ml_objective(m, &s_owned))
    };

    lbfgs_minimize(model, max_iter, &lower_bounds, obj_fn, grad_fn, grad_tol)
}

/// Apply a parameter step with projection onto lower bounds.
fn apply_step(
    params: &[f64],
    direction: &[f64],
    step: f64,
    bounds: &[Option<f64>],
    out: &mut [f64],
) {
    for (((o, &p), &d), bound) in out
        .iter_mut()
        .zip(params)
        .zip(direction)
        .zip(bounds)
    {
        *o = p + step * d;
        if let Some(lb) = bound {
            *o = o.max(*lb);
        }
    }
}

/// L-BFGS optimizer with projected lower bounds.
///
/// Uses the two-loop recursion for approximate inverse Hessian direction
/// and backtracking line search with Armijo condition.
fn lbfgs_minimize<F, G>(
    model: &mut Model,
    max_iter: usize,
    lower_bounds: &[Option<f64>],
    obj_fn: F,
    mut grad_fn: G,
    tolerance: Option<f64>,
) -> FitResult
where
    F: Fn(&Model) -> f64,
    G: FnMut(&mut Model, &[f64]) -> Vec<f64>,
{
    let gtol = tolerance.unwrap_or(1e-6);
    let memory_size = 10;
    let c1 = 1e-4; // Armijo condition constant
    let max_ls = 30; // Max line search steps

    let mut params = model.get_param_vec();
    let n = params.len();

    if n == 0 {
        return FitResult {
            params,
            objective: obj_fn(model),
            converged: true,
        };
    }

    let mut obj = obj_fn(model);
    let mut grad = grad_fn(model, &params);

    // L-BFGS storage: circular buffer of {s, y, rho} pairs
    let mut s_history: VecDeque<Vec<f64>> = VecDeque::with_capacity(memory_size);
    let mut y_history: VecDeque<Vec<f64>> = VecDeque::with_capacity(memory_size);
    let mut rho_history: VecDeque<f64> = VecDeque::with_capacity(memory_size);
    let mut converged = false;

    for _iter in 0..max_iter {
        // Check gradient convergence
        let grad_norm: f64 = grad.iter().map(|g| g * g).sum::<f64>().sqrt();
        if grad_norm < gtol {
            converged = true;
            break;
        }

        // Compute search direction via L-BFGS two-loop recursion
        let direction = lbfgs_direction(&grad, &s_history, &y_history, &rho_history);

        // Ensure descent direction, then compute directional derivative
        let dg: f64 = direction.iter().zip(grad.iter()).map(|(d, g)| d * g).sum();
        let (direction, dg) = if dg >= 0.0 {
            // Fall back to steepest descent
            let dir: Vec<_> = grad.iter().map(|g| -g).collect();
            let dg = dir.iter().zip(grad.iter()).map(|(d, g)| d * g).sum();
            (dir, dg)
        } else {
            (direction, dg)
        };

        // Backtracking line search with Armijo condition
        let mut step = 1.0;
        let mut new_params = vec![0.0; n];
        let mut new_obj = f64::INFINITY;
        let mut ls_success = false;

        for _ in 0..max_ls {
            apply_step(&params, &direction, step, lower_bounds, &mut new_params);

            model.set_param_vec(&new_params);
            new_obj = obj_fn(model);

            // Armijo sufficient decrease condition
            if new_obj.is_finite() && new_obj <= obj + c1 * step * dg {
                ls_success = true;
                break;
            }
            step *= 0.5;
        }

        if !ls_success {
            // Try a very small step as last resort
            apply_step(&params, &direction, 1e-4, lower_bounds, &mut new_params);
            model.set_param_vec(&new_params);
            new_obj = obj_fn(model);

            if !new_obj.is_finite() || new_obj >= obj {
                // Cannot make progress — converged or stuck
                model.set_param_vec(&params);
                break;
            }
        }

        // Compute new gradient
        let new_grad = grad_fn(model, &new_params);

        // Update L-BFGS history
        let s_k: Vec<f64> = new_params
            .iter()
            .zip(params.iter())
            .map(|(a, b)| a - b)
            .collect();
        let y_k: Vec<f64> = new_grad
            .iter()
            .zip(grad.iter())
            .map(|(a, b)| a - b)
            .collect();
        let sy: f64 = s_k.iter().zip(y_k.iter()).map(|(s, y)| s * y).sum();

        if sy > 1e-10 {
            // Only update if curvature condition holds
            if s_history.len() >= memory_size {
                s_history.pop_front();
                y_history.pop_front();
                rho_history.pop_front();
            }
            s_history.push_back(s_k);
            y_history.push_back(y_k);
            rho_history.push_back(1.0 / sy);
        }

        // Check for convergence: relative objective change
        let rel_change = (obj - new_obj).abs() / (obj.abs() + 1e-10);
        let param_change: f64 = new_params
            .iter()
            .zip(params.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();

        params = new_params;
        obj = new_obj;
        grad = new_grad;

        if rel_change < 1e-8 && param_change < 1e-6 {
            converged = true;
            break;
        }
    }

    model.set_param_vec(&params);
    FitResult {
        params,
        objective: obj,
        converged,
    }
}

/// L-BFGS two-loop recursion to compute search direction.
///
/// Returns: direction = -H_k * grad (approximate Newton direction)
fn lbfgs_direction(
    grad: &[f64],
    s_history: &VecDeque<Vec<f64>>,
    y_history: &VecDeque<Vec<f64>>,
    rho_history: &VecDeque<f64>,
) -> Vec<f64> {
    let n = grad.len();
    let m = s_history.len();

    if m == 0 {
        // No history: use steepest descent
        return grad.iter().map(|g| -g).collect();
    }

    // First loop: traverse history from most recent to oldest
    let mut q = grad.to_vec();
    let mut alpha = vec![0.0; m];

    for i in (0..m).rev() {
        alpha[i] = rho_history[i] * dot(&s_history[i], &q);
        for j in 0..n {
            q[j] -= alpha[i] * y_history[i][j];
        }
    }

    // Initial Hessian approximation: H0 = gamma * I
    // gamma = s_{k-1}' y_{k-1} / y_{k-1}' y_{k-1}
    let last = m - 1;
    let sy: f64 = dot(&s_history[last], &y_history[last]);
    let yy: f64 = dot(&y_history[last], &y_history[last]);
    let gamma = if yy > 1e-30 { sy / yy } else { 1.0 };

    let mut r: Vec<f64> = q.iter().map(|&qi| gamma * qi).collect();

    // Second loop: traverse from oldest to most recent
    for i in 0..m {
        let beta = rho_history[i] * dot(&y_history[i], &r);
        for j in 0..n {
            r[j] += (alpha[i] - beta) * s_history[i][j];
        }
    }

    // Negate for descent direction
    r.iter_mut().for_each(|x| *x = -*x);
    r
}

/// Dot product of two vectors.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// DWLS objective function value.
fn dwls_objective(model: &Model, s_vec: &[f64], w: &[f64]) -> f64 {
    let sigma = model.implied_cov();
    let sigma_vec = vech::vech(&sigma).expect("implied cov must be square");
    s_vec
        .iter()
        .zip(sigma_vec.iter())
        .zip(w.iter())
        .map(|((&s, &sig), &w)| {
            let r = s - sig;
            w * r * r
        })
        .sum()
}

/// Log-determinant of a symmetric positive-definite matrix via Cholesky.
fn log_det(mat: &Mat<f64>) -> Option<f64> {
    let llt = mat.llt(faer::Side::Lower).ok()?;
    let l = llt.L();
    let ld: f64 = (0..mat.nrows()).map(|i| l[(i, i)].ln()).sum();
    Some(2.0 * ld)
}

/// ML objective function value.
fn ml_objective(model: &Model, s_obs: &Mat<f64>) -> f64 {
    let sigma = model.implied_cov();
    let p = s_obs.nrows();

    let Some(log_det_sigma) = log_det(&sigma) else {
        return f64::INFINITY;
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

    let log_det_s = log_det(s_obs).unwrap_or(0.0);

    log_det_sigma + trace - log_det_s - p as f64
}

/// Numerical gradient via central differences.
fn numerical_gradient<F>(model: &mut Model, params: &[f64], obj_fn: F) -> Vec<f64>
where
    F: Fn(&Model) -> f64,
{
    let eps = 1e-7;
    let n = params.len();
    let mut grad = vec![0.0; n];
    let mut perturbed = params.to_vec();

    for i in 0..n {
        perturbed[i] += eps;
        model.set_param_vec(&perturbed);
        let f_plus = obj_fn(model);

        perturbed[i] -= 2.0 * eps;
        model.set_param_vec(&perturbed);
        let f_minus = obj_fn(model);

        perturbed[i] += eps; // restore
        grad[i] = (f_plus - f_minus) / (2.0 * eps);
    }

    model.set_param_vec(params);
    grad
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Model;
    use crate::syntax::parse_model;

    #[test]
    fn test_lbfgs_converges_simple_cfa() {
        // 1-factor CFA: F1 =~ V1 + V2 + V3
        // Must explicitly add residual variances (parser doesn't auto-add like lavaan)
        let s = faer::mat![[1.0, 0.6, 0.5], [0.6, 1.0, 0.4], [0.5, 0.4, 1.0],];
        let v_diag = vec![0.01; 6]; // kstar = 3*4/2 = 6

        let model_str = "F1 =~ NA*V1 + V2 + V3\nF1 ~~ 1*F1\nV1 ~~ V1\nV2 ~~ V2\nV3 ~~ V3";
        let pt = parse_model(model_str, false).unwrap();
        let obs_names: Vec<String> = vec!["V1", "V2", "V3"]
            .into_iter()
            .map(String::from)
            .collect();
        let mut model = Model::from_partable(&pt, &obs_names);

        let result = fit_dwls(&mut model, &s, &v_diag, 1000, None);
        assert!(result.converged, "L-BFGS should converge for simple CFA");
        assert!(
            result.objective < 0.1,
            "Objective should be small: {}",
            result.objective
        );
    }
}
