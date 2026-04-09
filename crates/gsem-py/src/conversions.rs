//! Conversion utilities between Python/NumPy and Rust types.
//!
//! - `faer::Mat<f64>` <-> flat row-major Vec (for NumPy array construction)
//! - `LdscResult` <-> JSON (matrices as row-major 2D arrays)
//!
//! The PyO3-exported functions in `py_exports.rs` use these for
//! converting between Rust and Python data structures.

use faer::Mat;
use gsem_ldsc::LdscResult;

/// Convert a faer Mat to a flat column-major Vec (NumPy default).
pub fn mat_to_numpy_flat(mat: &Mat<f64>) -> (Vec<f64>, usize, usize) {
    let nrows = mat.nrows();
    let ncols = mat.ncols();
    // NumPy default is row-major (C order), but we can also do column-major (F order)
    // Use row-major for Python convention
    let mut data = Vec::with_capacity(nrows * ncols);
    for i in 0..nrows {
        for j in 0..ncols {
            data.push(mat[(i, j)]);
        }
    }
    (data, nrows, ncols)
}

/// Convert a flat row-major Vec from NumPy to a faer Mat.
pub fn numpy_flat_to_mat(data: &[f64], nrows: usize, ncols: usize) -> Mat<f64> {
    Mat::from_fn(nrows, ncols, |i, j| data[i * ncols + j])
}

/// Convert LdscResult to a Python-friendly dict-like structure (via JSON).
pub fn ldsc_to_python_dict(result: &LdscResult) -> String {
    result.to_json_string().unwrap_or_default()
}

/// Convert JSON string back to LdscResult.
pub fn python_dict_to_ldsc(json: &str) -> Option<LdscResult> {
    LdscResult::from_json_string(json).ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_numpy_roundtrip() {
        let mat = faer::mat![[1.0, 2.0], [3.0, 4.0],];
        let (flat, nrows, ncols) = mat_to_numpy_flat(&mat);
        assert_eq!(flat, vec![1.0, 2.0, 3.0, 4.0]);
        let recovered = numpy_flat_to_mat(&flat, nrows, ncols);
        for i in 0..2 {
            for j in 0..2 {
                assert!((recovered[(i, j)] - mat[(i, j)]).abs() < 1e-15);
            }
        }
    }
}
