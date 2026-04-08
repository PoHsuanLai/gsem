use faer::Mat;

/// Extract the lower triangle (including diagonal) of a symmetric matrix,
/// column-major order. This matches R's `lowerTriangle(mat, diag=TRUE)`.
///
/// For a k×k matrix, returns a vector of length k*(k+1)/2.
/// Order: (0,0), (1,0), (2,0), ..., (k-1,0), (1,1), (2,1), ..., (k-1,k-1)
pub fn vech(mat: &Mat<f64>) -> Vec<f64> {
    let k = mat.nrows();
    assert_eq!(k, mat.ncols(), "vech requires a square matrix");
    let n = k * (k + 1) / 2;
    let mut v = Vec::with_capacity(n);
    for j in 0..k {
        for i in j..k {
            v.push(mat[(i, j)]);
        }
    }
    v
}

/// Reconstruct a symmetric matrix from its half-vectorization.
///
/// Given a vector of length k*(k+1)/2, reconstruct the k×k symmetric matrix.
pub fn vech_reverse(v: &[f64], k: usize) -> Mat<f64> {
    assert_eq!(v.len(), k * (k + 1) / 2, "vector length must be k*(k+1)/2");
    let mut mat = Mat::zeros(k, k);
    let mut idx = 0;
    for j in 0..k {
        for i in j..k {
            mat[(i, j)] = v[idx];
            mat[(j, i)] = v[idx];
            idx += 1;
        }
    }
    mat
}

/// Return the (row, col) index pairs for vech ordering.
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
        let v = vech(&mat);
        assert_eq!(v, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_vech_3x3() {
        let mat = faer::mat![[1.0, 4.0, 5.0], [4.0, 2.0, 6.0], [5.0, 6.0, 3.0],];
        // Column-major lower tri: (0,0), (1,0), (2,0), (1,1), (2,1), (2,2)
        let v = vech(&mat);
        assert_eq!(v, vec![1.0, 4.0, 5.0, 2.0, 6.0, 3.0]);
    }

    #[test]
    fn test_vech_roundtrip() {
        let mat = faer::mat![[1.0, 0.5, 0.3], [0.5, 2.0, 0.7], [0.3, 0.7, 3.0],];
        let v = vech(&mat);
        let reconstructed = vech_reverse(&v, 3);
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
}
