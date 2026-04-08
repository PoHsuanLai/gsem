use faer::{Mat, Side};

use crate::error::MatrixError;

/// Compute the nearest positive-definite matrix using the Higham (2002)
/// alternating projections algorithm.
///
/// This is a port of R's `Matrix::nearPD(x, corr=FALSE, keepDiag=keepDiag)`.
///
/// # Algorithm
/// Alternates between:
/// 1. Projecting onto the PSD cone (clamp negative eigenvalues to eps)
/// 2. Projecting onto the set of matrices with the original diagonal (if keep_diag)
///
/// # Arguments
/// * `mat` - A symmetric matrix (not necessarily PD)
/// * `keep_diag` - If true, preserve the original diagonal elements
/// * `max_iter` - Maximum iterations (default: 100)
/// * `tol` - Convergence tolerance (default: 1e-8)
pub fn nearest_pd(
    mat: &Mat<f64>,
    keep_diag: bool,
    max_iter: usize,
    tol: f64,
) -> Result<Mat<f64>, MatrixError> {
    let n = mat.nrows();
    if n != mat.ncols() {
        return Err(MatrixError::NotSquare {
            rows: n,
            cols: mat.ncols(),
        });
    }

    if n == 0 {
        return Ok(Mat::zeros(0, 0));
    }

    // Save original diagonal
    let orig_diag: Vec<f64> = (0..n).map(|i| mat[(i, i)]).collect();

    // Ensure symmetry: X = (mat + mat') / 2
    let mut x = Mat::from_fn(n, n, |i, j| (mat[(i, j)] + mat[(j, i)]) / 2.0);

    // Dykstra correction term
    let mut ds = Mat::zeros(n, n);

    let mut iter = 0;
    loop {
        if iter >= max_iter {
            return Err(MatrixError::NearPdNotConverged {
                iterations: max_iter,
            });
        }

        let x_old = x.clone();

        // Apply Dykstra correction: R = X - DS
        let r = &x - &ds;

        // Project onto PSD cone
        let eigen = r
            .self_adjoint_eigen(Side::Lower)
            .map_err(|_| MatrixError::SingularMatrix)?;

        let u = eigen.U();
        let s_diag = eigen.S().column_vector();

        // Clamp negative eigenvalues to a small positive value
        let eps = f64::EPSILON * n as f64;
        let max_eval = (0..n).map(|i| s_diag[i]).fold(f64::NEG_INFINITY, f64::max);
        let eig_threshold = eps * max_eval.max(1.0);

        // X_new = U * max(S, threshold) * U'
        let mut d_clamped = Mat::zeros(n, n);
        for i in 0..n {
            d_clamped[(i, i)] = s_diag[i].max(eig_threshold);
        }

        x = &u * &d_clamped * u.transpose();

        // Update Dykstra correction
        ds = &x - &r;

        // Ensure symmetry
        for i in 0..n {
            for j in (i + 1)..n {
                let avg = (x[(i, j)] + x[(j, i)]) / 2.0;
                x[(i, j)] = avg;
                x[(j, i)] = avg;
            }
        }

        // Project onto original diagonal constraint
        if keep_diag {
            for i in 0..n {
                x[(i, i)] = orig_diag[i];
            }
        }

        // Check convergence: relative change in Frobenius norm
        let mut norm_diff = 0.0;
        let mut norm_old = 0.0;
        for i in 0..n {
            for j in 0..n {
                let d = x[(i, j)] - x_old[(i, j)];
                norm_diff += d * d;
                norm_old += x_old[(i, j)] * x_old[(i, j)];
            }
        }

        if norm_old > 0.0 && (norm_diff / norm_old).sqrt() < tol {
            break;
        }

        iter += 1;
    }

    // Final PSD enforcement: one more eigenvalue clamp without Dykstra
    let eigen = x
        .self_adjoint_eigen(Side::Lower)
        .map_err(|_| MatrixError::SingularMatrix)?;

    let u = eigen.U();
    let s_diag = eigen.S().column_vector();
    let eps = f64::EPSILON * n as f64;
    let max_eval = (0..n).map(|i| s_diag[i]).fold(f64::NEG_INFINITY, f64::max);
    let eig_threshold = eps * max_eval.max(1.0);

    let mut d_clamped = Mat::zeros(n, n);
    for i in 0..n {
        d_clamped[(i, i)] = s_diag[i].max(eig_threshold);
    }

    let mut result = &u * &d_clamped * u.transpose();

    // Ensure symmetry
    for i in 0..n {
        for j in (i + 1)..n {
            let avg = (result[(i, j)] + result[(j, i)]) / 2.0;
            result[(i, j)] = avg;
            result[(j, i)] = avg;
        }
    }

    // Restore diagonal if needed
    if keep_diag {
        for i in 0..n {
            result[(i, i)] = orig_diag[i];
        }
    }

    Ok(result)
}

/// Convenience wrapper with default parameters matching R's nearPD behavior.
pub fn nearest_pd_default(mat: &Mat<f64>) -> Result<Mat<f64>, MatrixError> {
    nearest_pd(mat, false, 100, 1e-8)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_already_pd() {
        let mat = faer::mat![
            [2.0, 1.0],
            [1.0, 2.0],
        ];
        let result = nearest_pd_default(&mat).unwrap();
        for i in 0..2 {
            for j in 0..2 {
                assert!(
                    (result[(i, j)] - mat[(i, j)]).abs() < 1e-6,
                    "({i},{j}): {} vs {}",
                    result[(i, j)],
                    mat[(i, j)]
                );
            }
        }
    }

    #[test]
    fn test_not_pd() {
        // eigenvalues are 3 and -1, so not PD
        let mat = faer::mat![
            [1.0, 2.0],
            [2.0, 1.0],
        ];
        let result = nearest_pd_default(&mat).unwrap();
        assert!(result.llt(Side::Lower).is_ok(), "result should be PD");
        assert!((result[(0, 1)] - result[(1, 0)]).abs() < 1e-12);
    }

    #[test]
    fn test_3x3_not_pd() {
        let mat = faer::mat![
            [1.0, 0.9, 0.9],
            [0.9, 1.0, 0.9],
            [0.9, 0.9, 0.0],
        ];
        let result = nearest_pd_default(&mat).unwrap();
        assert!(result.llt(Side::Lower).is_ok(), "result should be PD");
        for i in 0..3 {
            for j in (i + 1)..3 {
                assert!(
                    (result[(i, j)] - result[(j, i)]).abs() < 1e-12,
                    "not symmetric at ({i},{j})"
                );
            }
        }
    }

    #[test]
    fn test_keep_diag() {
        let mat = faer::mat![
            [1.0, 2.0],
            [2.0, 1.0],
        ];
        let result = nearest_pd(&mat, true, 100, 1e-8).unwrap();
        assert!((result[(0, 0)] - 1.0).abs() < 1e-12);
        assert!((result[(1, 1)] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_identity_unchanged() {
        let mat = Mat::<f64>::identity(3, 3);
        let result = nearest_pd_default(&mat).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (result[(i, j)] - expected).abs() < 1e-6,
                    "({i},{j}): {} vs {expected}",
                    result[(i, j)]
                );
            }
        }
    }
}
