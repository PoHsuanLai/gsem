//! Conversion utilities between R data structures and Rust types.
//!
//! These functions handle the translation layer:
//! - R matrix (column-major numeric) <-> faer::Mat<f64>
//! - R list with $S, $V, $I, $N, $m <-> gsem_ldsc::LdscResult
//! - Rust SemResult -> R list with $results and $modelfit
//! - Rust Vec<SnpResult> -> R data.frame
//!
//! When extendr is enabled, these use extendr_api types (Robj, List, RMatrix).
//! For now, they work with the JSON serialization layer.

use faer::Mat;
use gsem_ldsc::LdscResult;
use gsem_sem::SemResult;

/// Convert a faer Mat to a row-major Vec<Vec<f64>> (for R matrix construction).
pub fn mat_to_r_matrix(mat: &Mat<f64>) -> Vec<Vec<f64>> {
    (0..mat.nrows())
        .map(|i| (0..mat.ncols()).map(|j| mat[(i, j)]).collect())
        .collect()
}

/// Convert a column-major flat vector (R's matrix storage) to a faer Mat.
pub fn r_matrix_to_mat(data: &[f64], nrow: usize, ncol: usize) -> Mat<f64> {
    // R stores matrices in column-major order
    Mat::from_fn(nrow, ncol, |i, j| data[j * nrow + i])
}

/// Convert LdscResult to a JSON string for R interchange.
pub fn ldsc_result_to_json(result: &LdscResult) -> String {
    result.to_json_string().unwrap_or_default()
}

/// Convert JSON string to LdscResult for R interchange.
pub fn json_to_ldsc_result(json: &str) -> Option<LdscResult> {
    LdscResult::from_json_string(json).ok()
}

/// Convert SemResult parameters to a flat table (for R data.frame).
pub fn sem_result_to_table(result: &SemResult) -> Vec<Vec<String>> {
    result
        .parameters
        .iter()
        .map(|p| {
            vec![
                p.lhs.clone(),
                p.op.clone(),
                p.rhs.clone(),
                format!("{:.6}", p.est),
                format!("{:.6}", p.se),
                format!("{:.4}", p.z),
                format!("{:.6e}", p.p),
            ]
        })
        .collect()
}

/// Column names for SEM result table.
pub fn sem_result_headers() -> Vec<&'static str> {
    vec!["lhs", "op", "rhs", "est", "se", "z", "p"]
}

/// Convert model fit to a flat row (for R data.frame).
pub fn model_fit_to_row(result: &SemResult) -> Vec<(String, f64)> {
    vec![
        ("chisq".to_string(), result.fit.chisq),
        ("df".to_string(), result.fit.df as f64),
        ("p_chisq".to_string(), result.fit.p_chisq),
        ("AIC".to_string(), result.fit.aic),
        ("CFI".to_string(), result.fit.cfi),
        ("SRMR".to_string(), result.fit.srmr),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_r_matrix_roundtrip() {
        let mat = faer::mat![
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
        ];
        let r_data: Vec<f64> = {
            // Column-major: col0=[1,4], col1=[2,5], col2=[3,6]
            vec![1.0, 4.0, 2.0, 5.0, 3.0, 6.0]
        };
        let recovered = r_matrix_to_mat(&r_data, 2, 3);
        for i in 0..2 {
            for j in 0..3 {
                assert!((recovered[(i, j)] - mat[(i, j)]).abs() < 1e-15);
            }
        }
    }

    #[test]
    fn test_ldsc_json_roundtrip() {
        let result = LdscResult {
            s: faer::mat![[0.3, 0.1], [0.1, 0.4]],
            v: faer::mat![[0.01, 0.002, 0.003], [0.002, 0.02, 0.004], [0.003, 0.004, 0.03]],
            i_mat: faer::mat![[1.05, 0.02], [0.02, 1.03]],
            n_vec: vec![50000.0, 50000.0, 50000.0],
            m: 1000000.0,
        };
        let json = ldsc_result_to_json(&result);
        let recovered = json_to_ldsc_result(&json).unwrap();
        assert!((recovered.s[(0, 0)] - 0.3).abs() < 1e-10);
        assert!((recovered.m - 1000000.0).abs() < 1e-10);
    }
}
