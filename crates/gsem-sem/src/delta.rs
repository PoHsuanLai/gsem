use faer::Mat;
use gsem_matrix::error::MatrixError;

use crate::jacobian;
use crate::model::Model;

/// Compute the delta (Jacobian) matrix: d(vech(Sigma)) / d(theta).
///
/// Uses analytical derivatives for efficiency (single pass instead of 2×n_free
/// evaluations of implied_cov).
///
/// Output: kstar × n_free matrix.
pub fn compute_delta(model: &mut Model, _eps: f64) -> Result<Mat<f64>, MatrixError> {
    Ok(jacobian::analytical_jacobian(model))
}
