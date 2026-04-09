use faer::Mat;
use gsem_ldsc::LdscResult;

/// Convert a faer Mat to a flat row-major Vec for NumPy.
pub fn mat_to_flat(mat: &Mat<f64>) -> (Vec<f64>, usize, usize) {
    let nrows = mat.nrows();
    let ncols = mat.ncols();
    let mut data = Vec::with_capacity(nrows * ncols);
    for i in 0..nrows {
        for j in 0..ncols {
            data.push(mat[(i, j)]);
        }
    }
    (data, nrows, ncols)
}

/// Convert JSON string back to LdscResult.
pub fn json_to_ldsc(json: &str) -> Option<LdscResult> {
    LdscResult::from_json_string(json).ok()
}
