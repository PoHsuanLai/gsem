use faer::{Mat, Side};
use gsem_matrix::vech;
use indicatif::{ProgressBar, ProgressStyle};
use rand::RngExt;
use rand_distr::StandardNormal;
use rayon::prelude::*;

/// Result of parallel analysis.
#[derive(Debug, Clone)]
pub struct PaResult {
    /// Observed eigenvalues
    pub observed: Vec<f64>,
    /// 95th percentile of simulated eigenvalues
    pub simulated_95: Vec<f64>,
    /// Suggested number of factors
    pub n_factors: usize,
}

/// Perform parallel analysis on an LDSC genetic correlation matrix.
///
/// Monte Carlo simulation: generates random correlation matrices from
/// N(null_vec, V) where null_vec is vech of the diagonal-only S,
/// then compares observed eigenvalues to the given percentile threshold.
///
/// `percentile`: quantile threshold (e.g. 0.95 for 95th percentile).
/// `diag_only`: if true, use only `diag(V)` for simulation (ignore off-diagonal sampling covariance).
///
/// Port of GenomicSEM's `paLDSC()`.
pub fn parallel_analysis(
    s: &Mat<f64>,
    v: &Mat<f64>,
    n_sim: usize,
    percentile: f64,
    diag_only: bool,
) -> PaResult {
    let k = s.nrows();
    let kstar = k * (k + 1) / 2;

    // Convert S to correlation matrix for eigenanalysis
    let sds: Vec<f64> = (0..k).map(|i| s[(i, i)].sqrt().max(1e-10)).collect();
    let cor = Mat::from_fn(k, k, |i, j| s[(i, j)] / (sds[i] * sds[j]));

    // Observed eigenvalues
    let observed = eigenvalues_sorted(&cor);

    // Null model: identity correlation matrix (no off-diagonal correlations)
    let s_null = Mat::<f64>::identity(k, k);
    let null_vec = vech::vech(&s_null).expect("identity matrix vech should not fail");

    // When diag_only=true, use only diagonal of V (ignore sampling covariance between elements)
    let v_sim = if diag_only {
        Mat::from_fn(kstar, kstar, |i, j| if i == j { v[(i, j)] } else { 0.0 })
    } else {
        v.to_owned()
    };

    // Cholesky of V for multivariate normal sampling: L such that V = L L'
    let chol_l = compute_cholesky_l(&v_sim, kstar);

    // Simulate null eigenvalue distributions (parallelized — each iteration is independent)
    let pb = ProgressBar::new(n_sim as u64);
    pb.set_style(
        ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} simulations ({eta})",
        )
        .unwrap(),
    );

    let all_eigs: Vec<Vec<f64>> = (0..n_sim)
        .into_par_iter()
        .map(|_| {
            let mut rng = rand::rng();
            let mut z = vec![0.0; kstar];
            let mut sample = vec![0.0; kstar];

            // Generate z ~ N(0, I)
            for zi in z.iter_mut() {
                *zi = rng.sample(StandardNormal);
            }

            // sample = null_vec + L * z (multivariate normal)
            for i in 0..kstar {
                let lz: f64 = (0..kstar).map(|j| chol_l[(i, j)] * z[j]).sum();
                sample[i] = null_vec[i] + lz;
            }

            let sim_mat = vech::vech_reverse(&sample, k)
                .expect("vech_reverse should not fail for valid kstar");
            let eigs = eigenvalues_sorted(&sim_mat);
            pb.inc(1);
            eigs
        })
        .collect();

    pb.finish_with_message("complete");

    // Transpose: from n_sim × k to k × n_sim
    let mut sim_eigenvalues: Vec<Vec<f64>> = vec![Vec::with_capacity(n_sim); k];
    for eigs in &all_eigs {
        for (f, &eig) in eigs.iter().enumerate() {
            sim_eigenvalues[f].push(eig);
        }
    }

    // Compute percentile for each eigenvalue using partial sort
    let percentile_idx = ((n_sim as f64) * percentile) as usize;
    let simulated_95: Vec<f64> = sim_eigenvalues
        .iter_mut()
        .map(|eigs| {
            let idx = percentile_idx.min(eigs.len().saturating_sub(1));
            eigs.select_nth_unstable_by(idx, |a, b| {
                a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
            });
            eigs[idx]
        })
        .collect();

    // Count factors: observed > simulated 95th percentile
    let n_factors = observed
        .iter()
        .zip(simulated_95.iter())
        .take_while(|(obs, sim)| obs > sim)
        .count()
        .max(1);

    PaResult {
        observed,
        simulated_95,
        n_factors,
    }
}

/// Compute Cholesky factor L of V, with fallback to nearPD then diagonal.
fn compute_cholesky_l(v: &Mat<f64>, kstar: usize) -> Mat<f64> {
    if let Ok(llt) = v.llt(Side::Lower) {
        let l_ref = llt.L();
        return Mat::from_fn(kstar, kstar, |i, j| l_ref[(i, j)]);
    }

    // V is not PD; apply nearPD then retry
    #[allow(clippy::collapsible_if)]
    if let Ok(v_pd) = gsem_matrix::near_pd::nearest_pd(v, false, 100, 1e-8) {
        if let Ok(llt) = v_pd.llt(Side::Lower) {
            let l_ref = llt.L();
            return Mat::from_fn(kstar, kstar, |i, j| l_ref[(i, j)]);
        }
    }

    // Last resort: use diagonal sqrt
    Mat::from_fn(kstar, kstar, |i, j| {
        if i == j {
            v[(i, i)].max(0.0).sqrt()
        } else {
            0.0
        }
    })
}

/// Get eigenvalues of a symmetric matrix, sorted descending.
fn eigenvalues_sorted(mat: &Mat<f64>) -> Vec<f64> {
    let Ok(eigen) = mat.self_adjoint_eigen(Side::Lower) else {
        return vec![0.0; mat.nrows()];
    };
    let s = eigen.S().column_vector();
    let mut eigs: Vec<f64> = (0..mat.nrows()).map(|i| s[i]).collect();
    eigs.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    eigs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parallel_analysis_identity() {
        // Identity correlation matrix: all eigenvalues = 1
        let s = Mat::<f64>::identity(3, 3);
        let v = Mat::from_fn(6, 6, |i, j| if i == j { 0.01 } else { 0.0 });
        let result = parallel_analysis(&s, &v, 100, 0.95, false);
        assert_eq!(result.observed.len(), 3);
        for &eig in &result.observed {
            assert!((eig - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_parallel_analysis_one_factor() {
        // Strong 1-factor structure
        let s = faer::mat![[1.0, 0.8, 0.8], [0.8, 1.0, 0.8], [0.8, 0.8, 1.0],];
        let v = Mat::from_fn(6, 6, |i, j| if i == j { 0.001 } else { 0.0 });
        let result = parallel_analysis(&s, &v, 200, 0.95, false);
        assert!(result.n_factors >= 1);
    }
}
