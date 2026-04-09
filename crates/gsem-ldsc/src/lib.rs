//! LD Score Regression engine for estimating genetic covariance.
//!
//! Standalone crate implementing multivariate LDSC with block jackknife,
//! producing the S (genetic covariance), V (sampling covariance), and
//! I (intercept) matrices needed for downstream SEM analysis.
//!
//! # Entry point
//!
//! Use [`ldsc()`] to run the full pipeline. Results are returned as
//! [`LdscResult`], which supports JSON serialization.

pub mod annot_reader;
pub mod covariance;
pub mod error;
pub mod hdl;
pub mod heritability;
pub mod jackknife;
pub mod liability;
pub mod regression;
pub mod stratified;
pub mod weights;

use std::collections::HashMap;

use anyhow::Result;
use faer::Mat;
use gsem_matrix::vech;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Result of multivariate LD Score Regression.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LdscResult {
    /// Genetic covariance matrix (k×k)
    #[serde(serialize_with = "ser_mat", deserialize_with = "de_mat")]
    pub s: Mat<f64>,
    /// Sampling covariance of vech(S) (kstar×kstar where kstar=k*(k+1)/2)
    #[serde(serialize_with = "ser_mat", deserialize_with = "de_mat")]
    pub v: Mat<f64>,
    /// LDSC intercept matrix (k×k)
    #[serde(serialize_with = "ser_mat", deserialize_with = "de_mat")]
    pub i_mat: Mat<f64>,
    /// Sample sizes per vech element (kstar×1)
    pub n_vec: Vec<f64>,
    /// Total number of SNPs used
    pub m: f64,
}

/// Serialize a `Mat<f64>` as a row-major 2D array: `[[1,2],[3,4]]`.
pub fn ser_mat<S: serde::Serializer>(mat: &Mat<f64>, s: S) -> std::result::Result<S::Ok, S::Error> {
    use serde::ser::SerializeSeq;
    let mut seq = s.serialize_seq(Some(mat.nrows()))?;
    for i in 0..mat.nrows() {
        let row: Vec<f64> = (0..mat.ncols()).map(|j| mat[(i, j)]).collect();
        seq.serialize_element(&row)?;
    }
    seq.end()
}

/// Deserialize a `Mat<f64>` from a row-major 2D array: `[[1,2],[3,4]]`.
pub fn de_mat<'de, D: serde::Deserializer<'de>>(d: D) -> std::result::Result<Mat<f64>, D::Error> {
    let rows: Vec<Vec<f64>> = Deserialize::deserialize(d)?;
    if rows.is_empty() {
        return Ok(Mat::zeros(0, 0));
    }
    let nrows = rows.len();
    let ncols = rows[0].len();
    Ok(Mat::from_fn(nrows, ncols, |i, j| rows[i][j]))
}

/// Input data for a single trait's summary statistics.
#[derive(Debug, Clone)]
pub struct TraitSumstats {
    pub snp: Vec<String>,
    pub z: Vec<f64>,
    pub n: Vec<f64>,
    pub a1: Vec<String>,
    pub a2: Vec<String>,
}

/// Merged data for a single trait pair, aligned by SNP.
#[derive(Debug, Clone)]
struct PairMergedData {
    /// Z-scores for trait j (n_snps)
    z_j: Vec<f64>,
    /// Z-scores for trait k (n_snps, same as z_j for diagonal)
    z_k: Vec<f64>,
    /// Sample sizes for trait j (n_snps)
    n_j: Vec<f64>,
    /// Sample sizes for trait k (n_snps)
    n_k: Vec<f64>,
    /// LD scores per SNP
    ld: Vec<f64>,
    /// Weight LD scores per SNP
    w_ld: Vec<f64>,
    /// Total M (number of SNPs in reference)
    m: f64,
    /// Number of SNPs after merging
    n_snps: usize,
}

/// Result from a single trait-pair regression (used internally for parallel collection).
struct PairResult {
    j: usize,
    jj: usize,
    value: f64,
    intercept: f64,
    mean_n: f64,
    pseudos: Vec<f64>,
}

/// Configuration for LDSC.
#[derive(Debug, Clone)]
pub struct LdscConfig {
    pub n_blocks: usize,
    pub chisq_max: Option<f64>,
}

impl Default for LdscConfig {
    fn default() -> Self {
        Self {
            n_blocks: 200,
            chisq_max: None,
        }
    }
}

/// Run multivariate LD Score Regression.
///
/// This is the main entry point, equivalent to R's `ldsc()`.
/// Merges SNPs per trait-pair (not globally) to match R's behavior.
#[allow(clippy::too_many_arguments)]
pub fn ldsc(
    traits: &[TraitSumstats],
    sample_prev: &[Option<f64>],
    pop_prev: &[Option<f64>],
    ld_scores: &[f64],
    w_ld_scores: &[f64],
    ld_snps: &[String],
    m_total: f64,
    config: &LdscConfig,
    on_pair_done: Option<&(dyn Fn() + Sync)>,
) -> Result<LdscResult> {
    let k = traits.len();
    if k == 0 {
        anyhow::bail!("no traits provided");
    }

    // Build trait SNP lookups: SNP -> index
    let trait_maps: Vec<HashMap<&str, usize>> = traits
        .iter()
        .map(|t| {
            t.snp
                .iter()
                .enumerate()
                .map(|(i, s)| (s.as_str(), i))
                .collect()
        })
        .collect();

    let kstar = vech::vech_size(k);

    // Build (j, jj) pairs in vech order for parallel processing
    let pairs: Vec<(usize, usize)> = (0..k).flat_map(|j| (j..k).map(move |jj| (j, jj))).collect();

    // Run all trait pair regressions in parallel, each with its own per-pair merge
    let pair_results: Vec<Result<PairResult>> = pairs
        .par_iter()
        .map(|&(j, jj)| {
            // Merge data for this specific pair
            let merged = merge_pair(
                traits,
                ld_snps,
                ld_scores,
                w_ld_scores,
                &trait_maps,
                j,
                jj,
                m_total,
                config,
            )?;
            let n_blocks = config.n_blocks.min(merged.n_snps);

            let pair_result = if j == jj {
                let result = heritability::estimate_h2(
                    &merged.z_j,
                    &merged.n_j,
                    &merged.ld,
                    &merged.w_ld,
                    merged.m,
                    n_blocks,
                )?;
                Ok(PairResult {
                    j,
                    jj,
                    value: result.h2,
                    intercept: result.intercept,
                    mean_n: result.mean_n,
                    pseudos: result.pseudo_values,
                })
            } else {
                let result = covariance::estimate_gcov(
                    &merged.z_j,
                    &merged.z_k,
                    &merged.n_j,
                    &merged.n_k,
                    &merged.ld,
                    &merged.w_ld,
                    merged.m,
                    n_blocks,
                )?;
                Ok(PairResult {
                    j,
                    jj,
                    value: result.gcov,
                    intercept: result.intercept,
                    mean_n: result.mean_n,
                    pseudos: result.pseudo_values,
                })
            };
            if let Some(cb) = on_pair_done {
                cb();
            }
            pair_result
        })
        .collect();

    // Collect results into matrices (sequential, preserves vech order)
    let mut s = Mat::zeros(k, k);
    let mut i_mat = Mat::zeros(k, k);
    let mut n_vec = vec![0.0; kstar];
    let mut all_pseudos: Vec<Vec<f64>> = Vec::with_capacity(kstar);

    for (vech_idx, pr) in pair_results.into_iter().enumerate() {
        let r = pr?;
        s[(r.j, r.jj)] = r.value;
        i_mat[(r.j, r.jj)] = r.intercept;
        if r.j != r.jj {
            s[(r.jj, r.j)] = r.value;
            i_mat[(r.jj, r.j)] = r.intercept;
        }
        n_vec[vech_idx] = r.mean_n;
        all_pseudos.push(r.pseudos);
    }

    // Construct V matrix from jackknife pseudo-values.
    // n_vec contains per-element mean N (N_bar for diag, sqrt(N1*N2) for off-diag).
    let n_blocks = config.n_blocks;
    let v = jackknife::construct_v_matrix(&all_pseudos, n_blocks, &n_vec, m_total);

    // Apply liability scale conversion if prevalences provided
    let mut result = LdscResult {
        s,
        v,
        i_mat,
        n_vec,
        m: m_total,
    };

    liability::apply_liability_scale(&mut result, sample_prev, pop_prev);

    Ok(result)
}

/// Merge two traits with LD scores on SNP ID (per-pair, matching R's behavior).
///
/// R merges per trait-pair, then sorts by CHR, BP (genomic position).
/// We iterate over LD SNPs in order (which are already in genomic order,
/// read chromosome by chromosome) to match R's block structure.
#[allow(clippy::too_many_arguments)]
fn merge_pair(
    traits: &[TraitSumstats],
    ld_snps: &[String],
    ld_scores: &[f64],
    w_ld_scores: &[f64],
    trait_maps: &[HashMap<&str, usize>],
    j: usize,
    jj: usize,
    m_total: f64,
    config: &LdscConfig,
) -> Result<PairMergedData> {
    // Iterate over LD SNPs in genomic order (chr1, chr2, ..., chr22).
    // This matches R's `merged[order(CHR, BP), ]` which is critical for
    // jackknife block structure and V matrix accuracy.
    let mut z_j = Vec::new();
    let mut z_k = Vec::new();
    let mut n_j = Vec::new();
    let mut n_k = Vec::new();
    let mut ld = Vec::new();
    let mut w_ld = Vec::new();

    for (ld_idx, snp) in ld_snps.iter().enumerate() {
        let snp_str = snp.as_str();

        // Must be in trait j
        let Some(&idx_j) = trait_maps[j].get(snp_str) else {
            continue;
        };

        // Must be in trait jj
        let Some(&idx_jj) = trait_maps[jj].get(snp_str) else {
            continue;
        };

        // Chi-square max filter (max across the two traits in this pair)
        let chi_j = traits[j].z[idx_j].powi(2);
        let chi_k = traits[jj].z[idx_jj].powi(2);
        let max_chi = chi_j.max(chi_k);

        let chisq_max = config.chisq_max.unwrap_or_else(|| {
            let max_n = traits[j].n[idx_j].max(traits[jj].n[idx_jj]);
            (0.001 * max_n).max(80.0)
        });

        if max_chi > chisq_max {
            continue;
        }

        z_j.push(traits[j].z[idx_j]);
        z_k.push(traits[jj].z[idx_jj]);
        n_j.push(traits[j].n[idx_j]);
        n_k.push(traits[jj].n[idx_jj]);
        ld.push(ld_scores[ld_idx]);
        w_ld.push(w_ld_scores[ld_idx]);
    }

    let n_snps = ld.len();
    if n_snps == 0 {
        anyhow::bail!("no common SNPs for trait pair ({j}, {jj})");
    }

    log::info!("Pair ({j},{jj}): {n_snps} SNPs after merge and chi2 filter");

    Ok(PairMergedData {
        z_j,
        z_k,
        n_j,
        n_k,
        ld,
        w_ld,
        m: m_total,
        n_snps,
    })
}

impl LdscResult {
    /// Serialize to JSON string.
    pub fn to_json_string(&self) -> Result<String> {
        Ok(serde_json::to_string_pretty(self)?)
    }

    /// Deserialize from JSON string.
    pub fn from_json_string(s: &str) -> Result<Self> {
        Ok(serde_json::from_str(s)?)
    }
}
