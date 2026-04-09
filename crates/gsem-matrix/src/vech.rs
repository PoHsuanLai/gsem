use faer::Mat;

use crate::error::MatrixError;

/// Extract the lower triangle (including diagonal) of a symmetric matrix,
/// column-major order. This matches R's `lowerTriangle(mat, diag=TRUE)`.
///
/// For a k×k matrix, returns a vector of length k*(k+1)/2.
/// Order: (0,0), (1,0), (2,0), ..., (k-1,0), (1,1), (2,1), ..., (k-1,k-1)
pub fn vech(mat: &Mat<f64>) -> Result<Vec<f64>, MatrixError> {
    let k = mat.nrows();
    if k != mat.ncols() {
        return Err(MatrixError::NotSquare {
            rows: k,
            cols: mat.ncols(),
        });
    }
    let n = k * (k + 1) / 2;
    let mut v = Vec::with_capacity(n);
    for j in 0..k {
        for i in j..k {
            v.push(mat[(i, j)]);
        }
    }
    Ok(v)
}

/// Reconstruct a symmetric matrix from its half-vectorization.
///
/// Given a vector of length k*(k+1)/2, reconstruct the k×k symmetric matrix.
pub fn vech_reverse(v: &[f64], k: usize) -> Result<Mat<f64>, MatrixError> {
    let expected = k * (k + 1) / 2;
    if v.len() != expected {
        return Err(MatrixError::DimensionMismatch {
            expected,
            got: v.len(),
        });
    }
    let mut mat = Mat::zeros(k, k);
    let mut idx = 0;
    for j in 0..k {
        for i in j..k {
            mat[(i, j)] = v[idx];
            mat[(j, i)] = v[idx];
            idx += 1;
        }
    }
    Ok(mat)
}

/// Return the vech indices for a k×k matrix.
///
/// Returns pairs (row, col) in the same order as `vech()`.
pub fn vech_indices(k: usize) -> Vec<(usize, usize)> {
    let mut indices = Vec::with_capacity(k * (k + 1) / 2);
    for j in 0..k {
        for i in j..k {
            indices.push((i, j));
        }
    }
    indices
}

/// Return the length of vech for a k×k matrix.
pub fn vech_size(k: usize) -> usize {
    k * (k + 1) / 2
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vech_2x2() {
        let mat = faer::mat![[1.0, 2.0], [2.0, 3.0],];
        let v = vech(&mat).unwrap();
        assert_eq!(v, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_vech_3x3() {
        let mat = faer::mat![[1.0, 4.0, 5.0], [4.0, 2.0, 6.0], [5.0, 6.0, 3.0],];
        let v = vech(&mat).unwrap();
        assert_eq!(v, vec![1.0, 4.0, 5.0, 2.0, 6.0, 3.0]);
    }

    #[test]
    fn test_vech_roundtrip() {
        let mat = faer::mat![[1.0, 0.5, 0.3], [0.5, 2.0, 0.7], [0.3, 0.7, 3.0],];
        let v = vech(&mat).unwrap();
        let reconstructed = vech_reverse(&v, 3).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                assert!((mat[(i, j)] - reconstructed[(i, j)]).abs() < 1e-15);
            }
        }
    }

    #[test]
    fn test_vech_size() {
        assert_eq!(vech_size(1), 1);
        assert_eq!(vech_size(2), 3);
        assert_eq!(vech_size(3), 6);
        assert_eq!(vech_size(10), 55);
    }

    #[test]
    fn test_vech_indices() {
        let idx = vech_indices(3);
        assert_eq!(idx, vec![(0, 0), (1, 0), (2, 0), (1, 1), (2, 1), (2, 2)]);
    }

    #[test]
    fn test_vech_not_square() {
        let mat = Mat::zeros(2, 3);
        assert!(vech(&mat).is_err());
    }

    #[test]
    fn test_vech_reverse_wrong_len() {
        assert!(vech_reverse(&[1.0, 2.0], 3).is_err());
    }
}
