//! Model-implied genetic correlation matrix.
//!
//! Wrapper around [`commonfactor`] that additionally computes:
//! - R: the model-implied genetic correlation matrix (standardized implied covariance)
//! - V_R: the sampling covariance of R (via delta method / numerical Jacobian)
//!
//! Port of R GenomicSEM's `rgmodel()`.

use faer::Mat;
use anyhow::Result;
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
pub fn run_rgmodel(
    s: &Mat<f64>,
    v: &Mat<f64>,
    estimation: &str,
) -> Result<RgModelResult> {
    // Fit common factor model
    let sem_result = crate::commonfactor::run_commonfactor(s, v, estimation)?;

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rgmodel_diagonal_ones() {
        // For a 2x2 case, the diagonal of R should be 1.0
        let s = faer::mat![[0.5, 0.1], [0.1, 0.3]];
        let kstar = 3;
        let v = Mat::<f64>::identity(kstar, kstar) * 0.001;

        let result = run_rgmodel(&s, &v, "DWLS").unwrap();
        assert!((result.r[(0, 0)] - 1.0).abs() < 1e-6);
        assert!((result.r[(1, 1)] - 1.0).abs() < 1e-6);
        assert_eq!(result.v_r.nrows(), kstar);
        assert_eq!(result.v_r.ncols(), kstar);
    }
}
