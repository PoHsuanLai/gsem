use faer::{Mat, Side};
use gsem_matrix::vech;
use statrs::distribution::{ChiSquared, ContinuousCDF};

/// Compute Q_SNP heterogeneity statistic for a factor.
///
/// Tests whether SNP effects are homogeneous across indicators of a factor.
///
/// Q_SNP = eta' * P * Eig^{-1} * P' * eta
/// where eta = vech(S_subset - Sigma_subset)
pub fn compute_q_snp(
    s_subset: &Mat<f64>,
    sigma_subset: &Mat<f64>,
    v_subset: &Mat<f64>,
) -> (f64, usize, f64) {
    let p = s_subset.nrows();
    let s_vec = vech::vech(s_subset);
    let sig_vec = vech::vech(sigma_subset);

    let eta: Vec<f64> = s_vec
        .iter()
        .zip(sig_vec.iter())
        .map(|(s, sig)| s - sig)
        .collect();

    let n = eta.len();

    let Ok(eigen) = v_subset.self_adjoint_eigen(Side::Lower) else {
        return (f64::NAN, 0, 1.0);
    };

    let u = eigen.U();
    let s_diag = eigen.S().column_vector();

    let q: f64 = (0..n)
        .filter_map(|i| {
            let eigenval = s_diag[i];
            (eigenval > 1e-10).then(|| {
                let dot: f64 = eta.iter().enumerate().map(|(j, &e)| u[(j, i)] * e).sum();
                dot * dot / eigenval
            })
        })
        .sum();

    let df = p.saturating_sub(1);
    let p_val = if df > 0 && q.is_finite() {
        ChiSquared::new(df as f64)
            .map(|chi2| 1.0 - chi2.cdf(q))
            .unwrap_or(1.0)
    } else {
        1.0
    };

    (q, df, p_val)
}
