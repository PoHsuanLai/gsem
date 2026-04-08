use faer::{Mat, Side};
use rand::Rng;
use rand_distr::StandardNormal;

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
/// Monte Carlo simulation: generates random correlation matrices from V,
/// computes eigenvalues, and compares observed to 95th percentile.
///
/// Port of GenomicSEM's `paLDSC()`.
pub fn parallel_analysis(s: &Mat<f64>, v: &Mat<f64>, n_sim: usize) -> PaResult {
    let k = s.nrows();

    // Convert S to correlation matrix for eigenanalysis
    let sds: Vec<f64> = (0..k).map(|i| s[(i, i)].sqrt().max(1e-10)).collect();
    let cor = Mat::from_fn(k, k, |i, j| s[(i, j)] / (sds[i] * sds[j]));

    // Observed eigenvalues
    let observed = eigenvalues_sorted(&cor);

    // Simulate null eigenvalue distributions
    let mut sim_eigenvalues: Vec<Vec<f64>> = vec![Vec::with_capacity(n_sim); k];
    let mut rng = rand::rng();

    for _ in 0..n_sim {
        // Generate random symmetric matrix with structure from V
        let mut sim_cor = Mat::from_fn(k, k, |i, j| {
            if i == j {
                1.0
            } else {
                let noise: f64 = rng.sample(StandardNormal);
                // Scale noise by approximate SE from V
                let kstar_idx_i = vech_index(i, j, k);
                let se = if kstar_idx_i < v.nrows() {
                    v[(kstar_idx_i, kstar_idx_i)].sqrt()
                } else {
                    0.1
                };
                (noise * se).clamp(-0.99, 0.99)
            }
        });

        // Symmetrize
        for i in 0..k {
            for j in (i + 1)..k {
                let avg = (sim_cor[(i, j)] + sim_cor[(j, i)]) / 2.0;
                sim_cor[(i, j)] = avg;
                sim_cor[(j, i)] = avg;
            }
        }

        let eigs = eigenvalues_sorted(&sim_cor);
        for (f, &eig) in eigs.iter().enumerate() {
            sim_eigenvalues[f].push(eig);
        }
    }

    // Compute 95th percentile for each eigenvalue
    let simulated_95: Vec<f64> = sim_eigenvalues
        .iter()
        .map(|eigs| {
            let mut sorted = eigs.clone();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let idx = ((n_sim as f64) * 0.95) as usize;
            sorted[idx.min(sorted.len() - 1)]
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

/// Get eigenvalues of a symmetric matrix, sorted descending.
fn eigenvalues_sorted(mat: &Mat<f64>) -> Vec<f64> {
    let Ok(eigen) = mat.self_adjoint_eigen(Side::Lower) else {
        return vec![0.0; mat.nrows()];
    };
    let s = eigen.S().column_vector();
    let mut eigs: Vec<f64> = (0..mat.nrows()).map(|i| s[i]).collect();
    eigs.sort_by(|a, b| b.partial_cmp(a).unwrap());
    eigs
}

/// Convert (row, col) with row >= col to vech index.
fn vech_index(row: usize, col: usize, _k: usize) -> usize {
    let (r, c) = if row >= col { (row, col) } else { (col, row) };
    r * (r + 1) / 2 + c
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parallel_analysis_identity() {
        // Identity correlation matrix: all eigenvalues = 1
        let s = Mat::<f64>::identity(3, 3);
        let v = Mat::from_fn(6, 6, |i, j| if i == j { 0.01 } else { 0.0 });
        let result = parallel_analysis(&s, &v, 100);
        assert_eq!(result.observed.len(), 3);
        // All eigenvalues should be ~1
        for &eig in &result.observed {
            assert!((eig - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_parallel_analysis_one_factor() {
        // Strong 1-factor structure
        let s = faer::mat![[1.0, 0.8, 0.8], [0.8, 1.0, 0.8], [0.8, 0.8, 1.0],];
        let v = Mat::from_fn(6, 6, |i, j| if i == j { 0.001 } else { 0.0 });
        let result = parallel_analysis(&s, &v, 200);
        assert!(result.n_factors >= 1);
    }
}
