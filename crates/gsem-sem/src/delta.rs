use faer::Mat;
use gsem_matrix::error::MatrixError;
use gsem_matrix::vech;

use crate::model::Model;

/// Compute the delta (Jacobian) matrix: d(vech(Sigma)) / d(theta).
///
/// Uses central finite differences with step size eps.
/// Output: kstar × n_free matrix.
pub fn compute_delta(model: &mut Model, eps: f64) -> Result<Mat<f64>, MatrixError> {
    let params = model.get_param_vec();
    let n_free = params.len();
    let p = model.obs_names.len();
    let kstar = p * (p + 1) / 2;

    let mut delta = Mat::zeros(kstar, n_free);

    for j in 0..n_free {
        let mut p_plus = params.clone();
        let mut p_minus = params.clone();
        p_plus[j] += eps;
        p_minus[j] -= eps;

        model.set_param_vec(&p_plus);
        let sigma_plus = vech::vech(&model.implied_cov())?;

        model.set_param_vec(&p_minus);
        let sigma_minus = vech::vech(&model.implied_cov())?;

        for i in 0..kstar {
            delta[(i, j)] = (sigma_plus[i] - sigma_minus[i]) / (2.0 * eps);
        }
    }

    // Restore original parameters
    model.set_param_vec(&params);
    Ok(delta)
}
