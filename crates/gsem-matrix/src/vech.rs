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

/// Which set of vech positions `subset_sv` indexes into.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubsetType {
    /// Full lower triangle including the diagonal (R's `TYPE = "S"` / `"S_Stand"`).
    /// Positions are numbered `1..=k(k+1)/2` in column-major lower-tri order.
    WithDiagonal,
    /// Strict lower triangle, excluding the diagonal (R's `TYPE = "R"` for a
    /// correlation matrix). Positions are numbered `1..=k(k-1)/2`.
    OffDiagonal,
}

/// Result of [`subset_sv`]: the selected genetic-covariance entries and the
/// matching block of the sampling-covariance matrix.
#[derive(Debug, Clone)]
pub struct SubsetSV {
    /// The selected entries of `vech(S)` (in `index_vals` order's natural
    /// vech ordering — i.e. ascending vech position, as R returns).
    pub sub_s: Vec<f64>,
    /// The `V` submatrix at rows/cols `index_vals`.
    pub sub_v: Mat<f64>,
}

/// Port of R GenomicSEM's `subSV`: subset the genetic-covariance vector and the
/// sampling-covariance matrix by a set of **vech-position indices**.
///
/// `S` is a k×k symmetric genetic-covariance (or correlation) matrix and `V` is
/// its `kstar × kstar` sampling-covariance matrix, where `kstar` is the number
/// of vech positions for `subset_type`. `index_vals` are **1-based** vech
/// positions (matching R's `INDEXVALS`), numbered in column-major lower-tri
/// order. Returns the selected `vech(S)` entries and the `V[idx, idx]` block.
///
/// This differs from `rgmodel(sub=)`, which subsets by whole *traits*; `subSV`
/// selects arbitrary individual covariance elements.
pub fn subset_sv(
    s: &Mat<f64>,
    v: &Mat<f64>,
    index_vals: &[usize],
    subset_type: SubsetType,
) -> Result<SubsetSV, MatrixError> {
    let k = s.nrows();
    if k != s.ncols() {
        return Err(MatrixError::NotSquare {
            rows: k,
            cols: s.ncols(),
        });
    }

    // Enumerate the (row, col) of each numbered vech position for this type,
    // column-major lower triangle — exactly how R builds `Snum`.
    let positions: Vec<(usize, usize)> = match subset_type {
        SubsetType::WithDiagonal => {
            let mut p = Vec::with_capacity(k * (k + 1) / 2);
            for j in 0..k {
                for i in j..k {
                    p.push((i, j));
                }
            }
            p
        }
        SubsetType::OffDiagonal => {
            let mut p = Vec::with_capacity(k * (k - 1) / 2);
            for j in 0..k {
                for i in (j + 1)..k {
                    p.push((i, j));
                }
            }
            p
        }
    };
    let kstar = positions.len();

    if v.nrows() != kstar || v.ncols() != kstar {
        return Err(MatrixError::DimensionMismatch {
            expected: kstar,
            got: v.nrows().min(v.ncols()),
        });
    }

    // 1-based positions, ascending, deduplicated — `%in%` keeps vech order.
    let mut idx0: Vec<usize> = index_vals.iter().map(|&p| p.saturating_sub(1)).collect();
    idx0.sort_unstable();
    idx0.dedup();
    for &i in &idx0 {
        if i >= kstar {
            return Err(MatrixError::DimensionMismatch {
                expected: kstar,
                got: i + 1,
            });
        }
    }

    let sub_s: Vec<f64> = idx0
        .iter()
        .map(|&i| {
            let (r, c) = positions[i];
            s[(r, c)]
        })
        .collect();
    let m = idx0.len();
    let sub_v = Mat::from_fn(m, m, |i, j| v[(idx0[i], idx0[j])]);

    Ok(SubsetSV { sub_s, sub_v })
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

    #[test]
    fn test_subset_sv_with_diagonal() {
        // vech(S) for this 3x3 (col-major lower tri incl diag) is
        // [s00, s10, s20, s11, s21, s22] = [0.6, 0.42, 0.1, 0.5, 0.08, 0.45].
        let s = faer::mat![[0.60, 0.42, 0.10], [0.42, 0.50, 0.08], [0.10, 0.08, 0.45]];
        let v = Mat::from_fn(6, 6, |i, j| {
            if i == j {
                0.001 * (i as f64 + 1.0)
            } else {
                0.0
            }
        });
        // 1-based positions 1,3,6 -> vech entries [0.6, 0.1, 0.45].
        let out = subset_sv(&s, &v, &[1, 3, 6], SubsetType::WithDiagonal).unwrap();
        assert_eq!(out.sub_s, vec![0.60, 0.10, 0.45]);
        assert_eq!(out.sub_v.nrows(), 3);
        assert!((out.sub_v[(0, 0)] - 0.001).abs() < 1e-15); // V[1,1]
        assert!((out.sub_v[(1, 1)] - 0.003).abs() < 1e-15); // V[3,3]
        assert!((out.sub_v[(2, 2)] - 0.006).abs() < 1e-15); // V[6,6]
    }

    #[test]
    fn test_subset_sv_off_diagonal() {
        // Off-diagonal numbering: strict lower tri, col-major.
        // positions = (1,0),(2,0),(2,1) -> values [0.42, 0.10, 0.08].
        let s = faer::mat![[0.60, 0.42, 0.10], [0.42, 0.50, 0.08], [0.10, 0.08, 0.45]];
        let v = Mat::from_fn(3, 3, |i, j| if i == j { 0.01 } else { 0.0 });
        let out = subset_sv(&s, &v, &[1, 3], SubsetType::OffDiagonal).unwrap();
        assert_eq!(out.sub_s, vec![0.42, 0.08]);
        assert_eq!(out.sub_v.nrows(), 2);
    }

    #[test]
    fn test_subset_sv_out_of_range_errors() {
        let s = faer::mat![[1.0, 0.0], [0.0, 1.0]];
        let v = Mat::<f64>::identity(3, 3); // kstar = 3 for k=2 with diagonal
        // Position 4 is out of range (only 3 vech positions exist).
        assert!(subset_sv(&s, &v, &[4], SubsetType::WithDiagonal).is_err());
    }
}
