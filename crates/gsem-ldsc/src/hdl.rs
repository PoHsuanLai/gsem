//! High-Definition Likelihood (HDL) estimation of genetic covariance.
//!
//! A likelihood-based alternative to LDSC that uses LD eigenvalues
//! from reference panels to estimate genetic covariance more precisely.
//!
//! The piecewise method (default in GenomicSEM) partitions the genome
//! into LD blocks and optimizes the likelihood per block.

use std::collections::HashMap;

use anyhow::Result;
use faer::Mat;
use serde::{Deserialize, Serialize};

/// HDL method: piecewise (default in GenomicSEM) or jackknife.
#[derive(Debug, Clone, Copy)]
pub enum HdlMethod {
    Piecewise,
    Jackknife,
}

/// Result of HDL analysis. Same structure as LdscResult for compatibility.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HdlResult {
    /// Genetic covariance matrix (k x k)
    #[serde(serialize_with = "crate::ser_mat", deserialize_with = "crate::de_mat")]
    pub s: Mat<f64>,
    /// Sampling covariance of vech(S) (kstar x kstar)
    #[serde(serialize_with = "crate::ser_mat", deserialize_with = "crate::de_mat")]
    pub v: Mat<f64>,
    /// Intercept matrix (k x k)
    #[serde(serialize_with = "crate::ser_mat", deserialize_with = "crate::de_mat")]
    pub i_mat: Mat<f64>,
    /// Total number of SNPs used
    pub m: f64,
}

/// Per-piece LD reference data (pre-computed).
pub struct LdPiece {
    /// SNP names in this piece
    pub snps: Vec<String>,
    /// A1 alleles
    pub a1: Vec<String>,
    /// A2 alleles
    pub a2: Vec<String>,
    /// LD scores (eigenvalues of LD correlation matrix)
    pub ld_scores: Vec<f64>,
    /// Number of SNPs in piece
    pub m: usize,
}

/// HDL configuration.
#[derive(Debug, Clone)]
pub struct HdlConfig {
    pub method: HdlMethod,
    /// Reference panel sample size (default: 335265 for UKB)
    pub n_ref: f64,
}

impl Default for HdlConfig {
    fn default() -> Self {
        Self {
            method: HdlMethod::Piecewise,
            n_ref: 335265.0,
        }
    }
}

/// Input trait data for HDL.
pub struct HdlTraitData {
    pub snp: Vec<String>,
    pub z: Vec<f64>,
    pub n: Vec<f64>,
    pub a1: Vec<String>,
    pub a2: Vec<String>,
}

/// Run HDL analysis.
///
/// Port of R GenomicSEM's `hdl()`.
///
/// Uses a likelihood-based approach per LD block:
///   For h2: minimize -2*loglik where L(h2, int) = N(bstar | 0, Sigma(h2, int))
///   For gcov: condition on h2 estimates, minimize for h12
pub fn hdl(
    traits: &[HdlTraitData],
    sample_prev: &[Option<f64>],
    pop_prev: &[Option<f64>],
    ld_pieces: &[LdPiece],
    config: &HdlConfig,
) -> Result<HdlResult> {
    let k = traits.len();
    if k == 0 {
        anyhow::bail!("no traits provided");
    }

    let num_pieces = ld_pieces.len();
    let m_total: f64 = ld_pieces.iter().map(|p| p.m as f64).sum();
    let kstar = k * (k + 1) / 2;

    let mut s = Mat::zeros(k, k);
    let mut i_mat = Mat::zeros(k, k);

    // Collect per-piece estimates for V matrix construction
    let mut all_piece_estimates: Vec<Vec<f64>> = Vec::with_capacity(kstar);

    for j in 0..k {
        for d in j..k {
            let mut piece_estimates = Vec::with_capacity(num_pieces);
            let mut piece_intercepts = Vec::new();

            if j == d {
                // Univariate: h2 estimation
                for piece in ld_pieces {
                    let (h2, intercept) = estimate_h2_piece(&traits[j], piece, config.n_ref);
                    piece_estimates.push(h2);
                    piece_intercepts.push(intercept);
                }

                let h2_total: f64 = piece_estimates.iter().sum();
                let int_mean = if piece_intercepts.is_empty() {
                    1.0
                } else {
                    piece_intercepts.iter().sum::<f64>() / piece_intercepts.len() as f64
                };

                s[(j, j)] = h2_total;
                i_mat[(j, j)] = int_mean;
            } else {
                // Bivariate: genetic covariance
                // First need per-piece h2 for both traits
                let mut h2_j_pieces: Vec<(f64, f64)> = Vec::new();
                let mut h2_d_pieces: Vec<(f64, f64)> = Vec::new();

                for piece in ld_pieces {
                    h2_j_pieces.push(estimate_h2_piece(&traits[j], piece, config.n_ref));
                    h2_d_pieces.push(estimate_h2_piece(&traits[d], piece, config.n_ref));
                }

                for (p_idx, piece) in ld_pieces.iter().enumerate() {
                    let (gcov, intercept) = estimate_gcov_piece(
                        &traits[j],
                        &traits[d],
                        piece,
                        config.n_ref,
                        h2_j_pieces[p_idx],
                        h2_d_pieces[p_idx],
                    );
                    piece_estimates.push(gcov);
                    piece_intercepts.push(intercept);
                }

                let gcov_total: f64 = piece_estimates.iter().sum();
                let int_mean = if piece_intercepts.is_empty() {
                    0.0
                } else {
                    piece_intercepts.iter().sum::<f64>() / piece_intercepts.len() as f64
                };

                s[(j, d)] = gcov_total;
                s[(d, j)] = gcov_total;
                i_mat[(j, d)] = int_mean;
                i_mat[(d, j)] = int_mean;
            }

            all_piece_estimates.push(piece_estimates);
        }
    }

    // Construct V matrix from leave-one-piece-out pseudo-values (piecewise method)
    let v = construct_hdl_v(&all_piece_estimates, num_pieces);

    // Apply liability scale if prevalences provided
    let mut result = HdlResult {
        s,
        v,
        i_mat,
        m: m_total,
    };

    // Convert to LdscResult for liability scaling, then convert back
    let n_vec: Vec<f64> = vec![0.0; kstar];
    let mut ldsc_compat = crate::LdscResult {
        s: result.s.to_owned(),
        v: result.v.to_owned(),
        i_mat: result.i_mat.to_owned(),
        n_vec,
        m: result.m,
    };
    crate::liability::apply_liability_scale(&mut ldsc_compat, sample_prev, pop_prev);
    result.s = ldsc_compat.s;
    result.v = ldsc_compat.v;

    Ok(result)
}

/// Estimate h2 for a single LD piece.
/// Returns (h2_contribution, intercept).
fn estimate_h2_piece(trait_data: &HdlTraitData, piece: &LdPiece, n_ref: f64) -> (f64, f64) {
    let m = piece.m as f64;

    // Match trait SNPs to piece SNPs
    let bhat = extract_bhat(trait_data, piece);
    if bhat.is_empty() {
        return (0.0, 1.0);
    }

    let mean_n = trait_data.n.iter().sum::<f64>() / trait_data.n.len() as f64;

    // Initial OLS: a11 = bhat^2, regress on LD scores
    let a11: Vec<f64> = bhat.iter().map(|b| b * b).collect();
    let (slope_ols, int_ols) = simple_regression(&a11, &piece.ld_scores);

    let h2_init = (slope_ols * m).clamp(0.001, 0.999);
    let int_init = (int_ols * mean_n).clamp(0.5, 2.0);

    // Optimize likelihood
    optimize_h2_likelihood(&bhat, &piece.ld_scores, mean_n, m, n_ref, h2_init, int_init)
}

/// Estimate genetic covariance for a single LD piece.
fn estimate_gcov_piece(
    trait_j: &HdlTraitData,
    trait_d: &HdlTraitData,
    piece: &LdPiece,
    n_ref: f64,
    h2_j: (f64, f64), // (h2, intercept) for trait j
    h2_d: (f64, f64), // (h2, intercept) for trait d
) -> (f64, f64) {
    let m = piece.m as f64;

    let bhat_j = extract_bhat(trait_j, piece);
    let bhat_d = extract_bhat(trait_d, piece);

    if bhat_j.is_empty() || bhat_d.is_empty() {
        return (0.0, 0.0);
    }

    let mean_n_j = trait_j.n.iter().sum::<f64>() / trait_j.n.len() as f64;
    let mean_n_d = trait_d.n.iter().sum::<f64>() / trait_d.n.len() as f64;
    let n0 = mean_n_j.min(mean_n_d); // approximate overlap

    // Initial OLS: a12 = bhat_j * bhat_d, regress on LD scores
    let a12: Vec<f64> = bhat_j
        .iter()
        .zip(bhat_d.iter())
        .map(|(a, b)| a * b)
        .collect();
    let (slope_ols, int_ols) = simple_regression(&a12, &piece.ld_scores);

    let h12_init = (slope_ols * m).clamp(-0.999, 0.999);
    let int_init = int_ols * (mean_n_j * mean_n_d).sqrt();

    // Optimize covariance likelihood conditioned on h2 estimates
    optimize_gcov_likelihood(
        &bhat_j,
        &bhat_d,
        &piece.ld_scores,
        mean_n_j,
        mean_n_d,
        n0,
        m,
        n_ref,
        h2_j,
        h2_d,
        h12_init,
        int_init,
    )
}

/// Extract aligned bhat values for SNPs in a piece.
fn extract_bhat(trait_data: &HdlTraitData, piece: &LdPiece) -> Vec<f64> {
    let trait_map: HashMap<&str, usize> = trait_data
        .snp
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();

    let mut bhat = vec![0.0; piece.snps.len()];
    let mut found = 0;

    for (p_idx, snp) in piece.snps.iter().enumerate() {
        if let Some(&t_idx) = trait_map.get(snp.as_str()) {
            let z = trait_data.z[t_idx];
            let n = trait_data.n[t_idx];
            if n > 0.0 {
                let mut b = z / n.sqrt();
                // Allele alignment
                if trait_data.a2[t_idx].to_uppercase() != piece.a2[p_idx].to_uppercase() {
                    b = -b;
                }
                bhat[p_idx] = b;
                found += 1;
            }
        }
    }

    if found == 0 {
        Vec::new()
    } else {
        bhat
    }
}

/// Simple OLS regression: y = slope * x + intercept.
fn simple_regression(y: &[f64], x: &[f64]) -> (f64, f64) {
    let n = y.len() as f64;
    if n < 2.0 {
        return (0.0, 0.0);
    }

    let mean_x = x.iter().sum::<f64>() / n;
    let mean_y = y.iter().sum::<f64>() / n;

    let ss_xy: f64 = x
        .iter()
        .zip(y.iter())
        .map(|(&xi, &yi)| (xi - mean_x) * (yi - mean_y))
        .sum();
    let ss_xx: f64 = x.iter().map(|&xi| (xi - mean_x).powi(2)).sum();

    let slope = if ss_xx > 1e-30 { ss_xy / ss_xx } else { 0.0 };
    let intercept = mean_y - slope * mean_x;
    (slope, intercept)
}

/// Optimize h2 likelihood for a single piece.
/// Uses coordinate descent (simpler than L-BFGS-B for 2 parameters).
fn optimize_h2_likelihood(
    bhat: &[f64],
    ld_scores: &[f64],
    n: f64,
    m: f64,
    n_ref: f64,
    h2_init: f64,
    int_init: f64,
) -> (f64, f64) {
    let mut best_h2 = h2_init;
    let mut best_int = int_init;
    let mut best_ll = h2_neg_loglik(bhat, ld_scores, best_h2, best_int, n, m, n_ref);

    // Coordinate descent
    for _ in 0..50 {
        // Optimize h2
        for &delta in &[-0.1, -0.01, -0.001, 0.001, 0.01, 0.1] {
            let h2_try = (best_h2 + delta).clamp(0.0001, 0.9999);
            let ll = h2_neg_loglik(bhat, ld_scores, h2_try, best_int, n, m, n_ref);
            if ll < best_ll {
                best_h2 = h2_try;
                best_ll = ll;
            }
        }
        // Optimize intercept
        for &delta in &[-0.5, -0.1, -0.01, 0.01, 0.1, 0.5] {
            let int_try = (best_int + delta).clamp(0.1, 10.0);
            let ll = h2_neg_loglik(bhat, ld_scores, best_h2, int_try, n, m, n_ref);
            if ll < best_ll {
                best_int = int_try;
                best_ll = ll;
            }
        }
    }

    (best_h2, best_int)
}

/// Negative log-likelihood for h2 estimation.
/// L = sum(log(lamh2)) + sum(bstar^2 / lamh2)
fn h2_neg_loglik(
    bhat: &[f64],
    lam: &[f64], // LD scores
    h2: f64,
    intercept: f64,
    n: f64,
    m: f64,
    n_ref: f64,
) -> f64 {
    let floor = (-10.0_f64).exp();
    let mut ll = 0.0;

    for (&b, &l) in bhat.iter().zip(lam.iter()) {
        let lamh2 = ((h2 / m) * l * l - h2 * l / n_ref + intercept * l / n).max(floor);
        ll += lamh2.ln() + b * b / lamh2;
    }

    ll
}

/// Optimize genetic covariance likelihood for a single piece.
#[allow(clippy::too_many_arguments)]
fn optimize_gcov_likelihood(
    bhat_j: &[f64],
    bhat_d: &[f64],
    ld_scores: &[f64],
    n_j: f64,
    n_d: f64,
    n0: f64,
    m: f64,
    n_ref: f64,
    h2_j: (f64, f64),
    h2_d: (f64, f64),
    gcov_init: f64,
    int_init: f64,
) -> (f64, f64) {
    let mut best_gcov = gcov_init;
    let mut best_int = int_init;
    let mut best_ll = gcov_neg_loglik(
        bhat_j, bhat_d, ld_scores, best_gcov, best_int, n_j, n_d, n0, m, n_ref, h2_j, h2_d,
    );

    // Coordinate descent
    for _ in 0..50 {
        for &delta in &[-0.1, -0.01, -0.001, 0.001, 0.01, 0.1] {
            let gcov_try = (best_gcov + delta).clamp(-0.9999, 0.9999);
            let ll = gcov_neg_loglik(
                bhat_j, bhat_d, ld_scores, gcov_try, best_int, n_j, n_d, n0, m, n_ref, h2_j,
                h2_d,
            );
            if ll < best_ll {
                best_gcov = gcov_try;
                best_ll = ll;
            }
        }
        for &delta in &[-0.5, -0.1, -0.01, 0.01, 0.1, 0.5] {
            let int_try = (best_int + delta).clamp(-5.0, 5.0);
            let ll = gcov_neg_loglik(
                bhat_j, bhat_d, ld_scores, best_gcov, int_try, n_j, n_d, n0, m, n_ref, h2_j,
                h2_d,
            );
            if ll < best_ll {
                best_int = int_try;
                best_ll = ll;
            }
        }
    }

    (best_gcov, best_int)
}

/// Negative log-likelihood for genetic covariance estimation.
#[allow(clippy::too_many_arguments)]
fn gcov_neg_loglik(
    bhat_j: &[f64],
    bhat_d: &[f64],
    lam: &[f64],
    h12: f64,
    intercept: f64,
    n_j: f64,
    n_d: f64,
    n0: f64,
    m: f64,
    n_ref: f64,
    h2_j: (f64, f64),
    h2_d: (f64, f64),
) -> f64 {
    let floor = (-10.0_f64).exp();
    let p1 = if n_j > 0.0 { n0 / n_j } else { 0.0 };
    let p2 = if n_d > 0.0 { n0 / n_d } else { 0.0 };
    let mut ll = 0.0;

    for i in 0..lam.len() {
        let l = lam[i];
        let lam11 = ((h2_j.0 / m) * l * l - h2_j.0 * l / n_ref + h2_j.1 * l / n_j).max(floor);
        let lam22 = ((h2_d.0 / m) * l * l - h2_d.0 * l / n_ref + h2_d.1 * l / n_d).max(floor);

        let lam12 = if n0 > 0.0 {
            (h12 / m) * l * l + p1 * p2 * intercept * l / n0
        } else {
            (h12 / m) * l * l
        };

        let ustar = bhat_d[i] - (lam12 / lam11) * bhat_j[i];
        let lam22_1 = (lam22 - lam12 * lam12 / lam11).max(floor);

        ll += lam22_1.ln() + ustar * ustar / lam22_1;
    }

    ll
}

/// Construct V matrix from leave-one-piece-out pseudo-values.
fn construct_hdl_v(all_piece_estimates: &[Vec<f64>], num_pieces: usize) -> Mat<f64> {
    let kstar = all_piece_estimates.len();
    let np = num_pieces as f64;

    // Compute totals
    let totals: Vec<f64> = all_piece_estimates
        .iter()
        .map(|ests| ests.iter().sum::<f64>())
        .collect();

    // Leave-one-piece-out pseudo-values
    let mut pseudo_values: Vec<Vec<f64>> = vec![Vec::with_capacity(num_pieces); kstar];
    for p in 0..num_pieces {
        for (elem, ests) in all_piece_estimates.iter().enumerate() {
            let loo_sum = totals[elem] - ests[p];
            let pseudo = np * totals[elem] - (np - 1.0) * loo_sum;
            pseudo_values[elem].push(pseudo);
        }
    }

    // V = cov(pseudo) * (num_pieces - 1)
    let means: Vec<f64> = pseudo_values
        .iter()
        .map(|p| p.iter().sum::<f64>() / np)
        .collect();

    let mut v = Mat::zeros(kstar, kstar);
    for i in 0..kstar {
        for j in i..kstar {
            let cov: f64 = pseudo_values[i]
                .iter()
                .zip(pseudo_values[j].iter())
                .map(|(&a, &b)| (a - means[i]) * (b - means[j]))
                .sum::<f64>()
                / (np - 1.0);

            // Scale by (num_pieces - 1) matching R: V = cov(V.hold) * (num_pieces - 1)
            let val = cov * (np - 1.0);
            v[(i, j)] = val;
            v[(j, i)] = val;
        }
    }

    v
}

impl HdlResult {
    /// Serialize to JSON string.
    pub fn to_json_string(&self) -> Result<String> {
        Ok(serde_json::to_string_pretty(self)?)
    }

    /// Deserialize from JSON string.
    pub fn from_json_string(s: &str) -> Result<Self> {
        Ok(serde_json::from_str(s)?)
    }

    /// Convert to LdscResult for downstream SEM compatibility.
    pub fn to_ldsc_result(&self) -> crate::LdscResult {
        let k = self.s.nrows();
        let kstar = k * (k + 1) / 2;
        crate::LdscResult {
            s: self.s.to_owned(),
            v: self.v.to_owned(),
            i_mat: self.i_mat.to_owned(),
            n_vec: vec![0.0; kstar],
            m: self.m,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_regression() {
        let y = vec![3.0, 5.0, 7.0, 9.0];
        let x = vec![1.0, 2.0, 3.0, 4.0];
        let (slope, int) = simple_regression(&y, &x);
        assert!((slope - 2.0).abs() < 1e-10);
        assert!((int - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_h2_neg_loglik_finite() {
        let bhat = vec![0.01, -0.02, 0.005];
        let lam = vec![10.0, 20.0, 15.0];
        let ll = h2_neg_loglik(&bhat, &lam, 0.3, 1.0, 50000.0, 100.0, 335265.0);
        assert!(ll.is_finite());
    }

    #[test]
    fn test_construct_hdl_v_dimensions() {
        // 2 traits -> kstar = 3, with 5 pieces
        let estimates = vec![
            vec![0.1, 0.12, 0.09, 0.11, 0.08],
            vec![0.05, 0.06, 0.04, 0.055, 0.045],
            vec![0.2, 0.22, 0.19, 0.21, 0.18],
        ];
        let v = construct_hdl_v(&estimates, 5);
        assert_eq!(v.nrows(), 3);
        assert_eq!(v.ncols(), 3);
        // V should be symmetric
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (v[(i, j)] - v[(j, i)]).abs() < 1e-15,
                    "V should be symmetric"
                );
            }
        }
    }

    #[test]
    fn test_extract_bhat_empty_when_no_overlap() {
        let trait_data = HdlTraitData {
            snp: vec!["rs1".to_string(), "rs2".to_string()],
            z: vec![1.0, 2.0],
            n: vec![10000.0, 10000.0],
            a1: vec!["A".to_string(), "C".to_string()],
            a2: vec!["G".to_string(), "T".to_string()],
        };
        let piece = LdPiece {
            snps: vec!["rs99".to_string(), "rs100".to_string()],
            a1: vec!["A".to_string(), "C".to_string()],
            a2: vec!["G".to_string(), "T".to_string()],
            ld_scores: vec![10.0, 20.0],
            m: 2,
        };
        let bhat = extract_bhat(&trait_data, &piece);
        assert!(bhat.is_empty());
    }
}
