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

/// Merged data across all traits and LD scores, ready for LDSC regression.
#[derive(Debug, Clone)]
struct MergedData {
    /// Z-scores per trait, aligned by SNP (n_traits × n_snps)
    z: Vec<Vec<f64>>,
    /// Sample sizes per trait (n_traits × n_snps)
    n: Vec<Vec<f64>>,
    /// LD scores per SNP
    ld: Vec<f64>,
    /// Weight LD scores per SNP
    w_ld: Vec<f64>,
    /// Total M (number of SNPs in reference)
    m: f64,
    /// Number of SNPs after merging
    n_snps: usize,
    /// Number of traits
    _n_traits: usize,
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
) -> Result<LdscResult> {
    let k = traits.len();
    if k == 0 {
        anyhow::bail!("no traits provided");
    }

    // Merge traits with LD scores on SNP ID
    let merged = merge_data(traits, ld_scores, w_ld_scores, ld_snps, m_total, config)?;
    let n_snps = merged.n_snps;
    let n_blocks = config.n_blocks.min(n_snps);

    log::info!("LDSC: {k} traits, {n_snps} SNPs after merge, {n_blocks} jackknife blocks");

    let kstar = vech::vech_size(k);

    // Build (j, jj) pairs in vech order for parallel processing
    let pairs: Vec<(usize, usize)> = (0..k).flat_map(|j| (j..k).map(move |jj| (j, jj))).collect();

    // Run all trait pair regressions in parallel
    let pair_results: Vec<Result<PairResult>> = pairs
        .par_iter()
        .map(|&(j, jj)| {
            if j == jj {
                let result = heritability::estimate_h2(
                    &merged.z[j],
                    &merged.n[j],
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
                    &merged.z[j],
                    &merged.z[jj],
                    &merged.n[j],
                    &merged.n[jj],
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
            }
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
    let v = jackknife::construct_v_matrix(&all_pseudos, n_blocks, &n_vec, merged.m);

    // Apply liability scale conversion if prevalences provided
    let mut result = LdscResult {
        s,
        v,
        i_mat,
        n_vec,
        m: merged.m,
    };

    liability::apply_liability_scale(&mut result, sample_prev, pop_prev);

    Ok(result)
}

/// Merge trait summary statistics with LD scores on SNP ID.
fn merge_data(
    traits: &[TraitSumstats],
    ld_scores: &[f64],
    w_ld_scores: &[f64],
    ld_snps: &[String],
    m_total: f64,
    config: &LdscConfig,
) -> Result<MergedData> {
    let k = traits.len();

    // Build LD score lookup: SNP -> (index, l2, w_ld)
    let mut ld_map: HashMap<&str, (usize, f64, f64)> = HashMap::with_capacity(ld_snps.len());
    for (i, snp) in ld_snps.iter().enumerate() {
        ld_map.insert(snp.as_str(), (i, ld_scores[i], w_ld_scores[i]));
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

    // Find SNPs present in ALL traits AND in LD scores
    let mut common_snps: Vec<&str> = Vec::new();
    for snp in &traits[0].snp {
        let snp_str = snp.as_str();
        if ld_map.contains_key(snp_str) && trait_maps.iter().all(|m| m.contains_key(snp_str)) {
            common_snps.push(snp_str);
        }
    }

    if common_snps.is_empty() {
        anyhow::bail!("no common SNPs across traits and LD scores");
    }

    // Extract aligned data
    let mut z = vec![Vec::with_capacity(common_snps.len()); k];
    let mut n = vec![Vec::with_capacity(common_snps.len()); k];
    let mut ld = Vec::with_capacity(common_snps.len());
    let mut w_ld = Vec::with_capacity(common_snps.len());

    for &snp in &common_snps {
        let (_, l2, wl) = ld_map[snp];

        // Apply chi-square max filter
        let max_chi = traits
            .iter()
            .enumerate()
            .map(|(t, _)| {
                let idx = trait_maps[t][snp];
                traits[t].z[idx].powi(2)
            })
            .fold(0.0f64, f64::max);

        let chisq_max = config.chisq_max.unwrap_or_else(|| {
            let max_n: f64 = traits
                .iter()
                .enumerate()
                .map(|(t, _)| {
                    let idx = trait_maps[t][snp];
                    traits[t].n[idx]
                })
                .fold(0.0f64, f64::max);
            (0.001 * max_n).max(80.0)
        });

        if max_chi > chisq_max {
            continue;
        }

        for t in 0..k {
            let idx = trait_maps[t][snp];
            z[t].push(traits[t].z[idx]);
            n[t].push(traits[t].n[idx]);
        }
        ld.push(l2);
        w_ld.push(wl);
    }

    let n_snps = ld.len();
    log::info!(
        "Merged: {} common SNPs (from {} in trait 0, {} in LD scores)",
        n_snps,
        traits[0].snp.len(),
        ld_snps.len()
    );

    Ok(MergedData {
        z,
        n,
        ld,
        w_ld,
        m: m_total,
        n_snps,
        _n_traits: k,
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
