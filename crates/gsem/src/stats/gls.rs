use faer::Mat;
use faer::linalg::solvers::DenseSolveCore;

/// Result of Generalized Least Squares regression.
#[derive(Debug, Clone)]
pub struct GlsResult {
    /// Coefficient estimates
    pub beta: Vec<f64>,
    /// Standard errors
    pub se: Vec<f64>,
    /// Z-statistics
    pub z: Vec<f64>,
    /// P-values
    pub p: Vec<f64>,
}

/// Run GLS regression on genetic parameters.
///
/// Port of GenomicSEM's `summaryGLS()`.
///
/// beta_gls = (X' V^{-1} X)^{-1} X' V^{-1} y
/// SE = sqrt(diag((X' V^{-1} X)^{-1}))
pub fn summary_gls(x: &Mat<f64>, y: &[f64], v: &Mat<f64>) -> Option<GlsResult> {
    let n = x.nrows();
    let p = x.ncols();
    if n != y.len() || n != v.nrows() {
        return None;
    }

    // V^{-1}
    let v_inv = v.partial_piv_lu().inverse();

    // X' V^{-1}
    let xt_vinv = x.transpose() * &v_inv;

    // X' V^{-1} X
    let xt_vinv_x = &xt_vinv * x;

    // (X' V^{-1} X)^{-1}
    let bread = xt_vinv_x.partial_piv_lu().inverse();

    // X' V^{-1} y
    let y_mat = Mat::from_fn(n, 1, |i, _| y[i]);
    let xt_vinv_y = &xt_vinv * &y_mat;

    // beta = bread * X' V^{-1} y
    let beta_mat = &bread * &xt_vinv_y;

    let beta: Vec<f64> = (0..p).map(|i| beta_mat[(i, 0)]).collect();
    let se: Vec<f64> = (0..p)
        .map(|i| {
            let v = bread[(i, i)];
            if v > 0.0 { v.sqrt() } else { 0.0 }
        })
        .collect();
    let z: Vec<f64> = beta
        .iter()
        .zip(se.iter())
        .map(|(b, s)| if *s > 0.0 { b / s } else { 0.0 })
        .collect();
    let p_vals: Vec<f64> = z
        .iter()
        .map(|&zi| {
            use statrs::distribution::{ContinuousCDF, Normal};
            let n = Normal::standard();
            2.0 * n.cdf(-zi.abs())
        })
        .collect();

    Some(GlsResult {
        beta,
        se,
        z,
        p: p_vals,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gls_identity_weights() {
        // With V = I, GLS = OLS
        let x = faer::mat![[1.0, 1.0], [1.0, 2.0], [1.0, 3.0], [1.0, 4.0],];
        let y = vec![2.0, 4.0, 6.0, 8.0]; // y = 2x
        let v = Mat::<f64>::identity(4, 4);
        let result = summary_gls(&x, &y, &v).unwrap();
        assert_eq!(result.beta.len(), 2);
        // Slope should be ~2
        assert!((result.beta[1] - 2.0).abs() < 1e-10);
    }
}
