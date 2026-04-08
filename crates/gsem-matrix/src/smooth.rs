use faer::{Mat, Side};

use crate::near_pd::nearest_pd_default;

/// Check if a matrix is positive definite via Cholesky decomposition.
pub fn is_pd(mat: &Mat<f64>) -> bool {
    mat.llt(Side::Lower).is_ok()
}

/// Check if a matrix is positive semi-definite via eigendecomposition.
/// Returns true if all eigenvalues >= -tol.
pub fn is_psd(mat: &Mat<f64>, tol: f64) -> bool {
    let Ok(eigen) = mat.self_adjoint_eigen(Side::Lower) else {
        return false;
    };
    let s_diag = eigen.S().column_vector();
    for i in 0..mat.nrows() {
        if s_diag[i] < -tol {
            return false;
        }
    }
    true
}

/// If the matrix is not positive definite, smooth it via nearest_pd.
/// Returns true if smoothing was applied.
pub fn smooth_if_needed(mat: &mut Mat<f64>) -> bool {
    if is_pd(mat) {
        return false;
    }

    match nearest_pd_default(mat) {
        Ok(smoothed) => {
            *mat = smoothed;
            true
        }
        Err(e) => {
            log::warn!("nearest_pd failed: {e}");
            false
        }
    }
}

/// Convert a covariance matrix to a correlation matrix.
/// `cor[i,j] = cov[i,j] / sqrt(cov[i,i] * cov[j,j])`
pub fn cov_to_cor(cov: &Mat<f64>) -> Mat<f64> {
    let n = cov.nrows();
    assert_eq!(n, cov.ncols());
    let sds: Vec<f64> = (0..n).map(|i| cov[(i, i)].sqrt()).collect();
    Mat::from_fn(n, n, |i, j| {
        if sds[i] > 0.0 && sds[j] > 0.0 {
            cov[(i, j)] / (sds[i] * sds[j])
        } else {
            0.0
        }
    })
}

/// Convert a correlation matrix back to a covariance matrix given standard deviations.
pub fn cor_to_cov(cor: &Mat<f64>, sds: &[f64]) -> Mat<f64> {
    let n = cor.nrows();
    assert_eq!(n, cor.ncols());
    assert_eq!(n, sds.len());
    Mat::from_fn(n, n, |i, j| cor[(i, j)] * sds[i] * sds[j])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_pd_positive() {
        let mat = faer::mat![[2.0, 1.0], [1.0, 2.0],];
        assert!(is_pd(&mat));
    }

    #[test]
    fn test_is_pd_negative() {
        let mat = faer::mat![[1.0, 2.0], [2.0, 1.0],];
        assert!(!is_pd(&mat));
    }

    #[test]
    fn test_smooth_if_needed_already_pd() {
        let mut mat = faer::mat![[2.0, 1.0], [1.0, 2.0],];
        let original = mat.clone();
        let smoothed = smooth_if_needed(&mut mat);
        assert!(!smoothed);
        for i in 0..2 {
            for j in 0..2 {
                assert!((mat[(i, j)] - original[(i, j)]).abs() < 1e-15);
            }
        }
    }

    #[test]
    fn test_smooth_if_needed_not_pd() {
        let mut mat = faer::mat![[1.0, 2.0], [2.0, 1.0],];
        let smoothed = smooth_if_needed(&mut mat);
        assert!(smoothed);
        assert!(is_pd(&mat));
    }

    #[test]
    fn test_cov_to_cor() {
        let cov = faer::mat![[4.0, 2.0], [2.0, 9.0],];
        let cor = cov_to_cor(&cov);
        assert!((cor[(0, 0)] - 1.0).abs() < 1e-12);
        assert!((cor[(1, 1)] - 1.0).abs() < 1e-12);
        // 2 / (2*3) = 1/3
        assert!((cor[(0, 1)] - 1.0 / 3.0).abs() < 1e-12);
        assert!((cor[(1, 0)] - 1.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_cor_to_cov_roundtrip() {
        let cov = faer::mat![[4.0, 2.0], [2.0, 9.0],];
        let sds: Vec<f64> = (0..2usize).map(|i| f64::sqrt(cov[(i, i)])).collect();
        let cor = cov_to_cor(&cov);
        let recovered = cor_to_cov(&cor, &sds);
        for i in 0..2 {
            for j in 0..2 {
                assert!((recovered[(i, j)] - cov[(i, j)]).abs() < 1e-12);
            }
        }
    }
}
