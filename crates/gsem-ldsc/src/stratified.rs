//! Stratified (partitioned) LD Score Regression.
//!
//! Extends standard LDSC to multiple annotation categories, estimating
//! per-annotation contributions to genetic covariance and heritability.
//!
//! The key difference from standard LDSC is that regression uses p+1
//! predictors (p annotation LD scores + intercept) instead of 2.

use std::collections::HashMap;

use anyhow::Result;
use faer::Mat;
use faer::prelude::Solve;

use crate::weights;

/// Result of stratified LD Score Regression.
#[derive(Debug, Clone)]
pub struct StratifiedLdscResult {
    /// Annotation names.
    pub annotations: Vec<String>,
    /// Per-annotation S matrices (k x k genetic covariance). Length = n_annot.
    pub s_annot: Vec<Mat<f64>>,
    /// Per-annotation V matrices (kstar x kstar sampling covariance). Length = n_annot.
    pub v_annot: Vec<Mat<f64>>,
    /// Intercept matrix (k x k), same as standard LDSC.
    pub i_mat: Mat<f64>,
    /// Per-annotation M counts (number of SNPs in each annotation).
    pub m_annot: Vec<f64>,
    /// Total M across all annotations.
    pub m_total: f64,
    /// Proportion of SNPs per annotation.
    pub prop: Vec<f64>,
}

/// Configuration for stratified LDSC.
#[derive(Debug, Clone)]
pub struct StratifiedLdscConfig {
    /// Number of jackknife blocks.
    pub n_blocks: usize,
    /// Remove flanking SNPs (within `flank_kb` kb) around each annotation's SNPs
    /// before running the regression for that annotation.
    pub rm_flank: bool,
    /// Flanking window size in kilobases (default 500).
    pub flank_kb: u64,
}

impl Default for StratifiedLdscConfig {
    fn default() -> Self {
        Self {
            n_blocks: 200,
            rm_flank: false,
            flank_kb: 500,
        }
    }
}

/// Run stratified LD Score Regression.
///
/// # Arguments
/// * `traits` - Per-trait summary statistics.
/// * `sample_prev` - Sample prevalences (None for continuous traits).
/// * `pop_prev` - Population prevalences (None for continuous traits).
/// * `annot_ld` - Annotation-specific LD scores (n_snps x n_annot).
/// * `w_ld_scores` - Weight LD scores per SNP.
/// * `ld_snps` - SNP names in LD reference.
/// * `annot_names` - Annotation names.
/// * `m_annot` - Per-annotation M counts.
/// * `config` - Configuration.
#[allow(clippy::too_many_arguments)]
pub fn s_ldsc(
    traits: &[crate::TraitSumstats],
    sample_prev: &[Option<f64>],
    pop_prev: &[Option<f64>],
    annot_ld: &Mat<f64>,
    w_ld_scores: &[f64],
    ld_snps: &[String],
    annot_names: &[String],
    m_annot: &[f64],
    config: &StratifiedLdscConfig,
    snp_chr: Option<&[u32]>,
    snp_bp: Option<&[u64]>,
) -> Result<StratifiedLdscResult> {
    let k = traits.len();
    let n_annot = annot_names.len();

    if k == 0 {
        anyhow::bail!("no traits provided");
    }
    if n_annot == 0 {
        anyhow::bail!("no annotations provided");
    }

    // Merge traits with LD scores on SNP ID
    let merged = merge_stratified_data(traits, annot_ld, w_ld_scores, ld_snps)?;
    let n_snps = merged.n_snps;
    let n_blocks = config.n_blocks.min(n_snps);
    let m_total: f64 = m_annot.iter().sum();

    log::info!("s_ldsc: {k} traits, {n_annot} annotations, {n_snps} SNPs, {n_blocks} blocks");

    // Pre-compute flanking masks if rm_flank is enabled
    let flank_masks: Vec<Vec<bool>> = if config.rm_flank {
        match (snp_chr, snp_bp) {
            (Some(chr), Some(bp)) => compute_flank_masks(
                &merged.annot_ld,
                &merged.ld_indices,
                chr,
                bp,
                n_annot,
                config.flank_kb,
            ),
            _ => {
                log::warn!("rm_flank=true but chr/bp data not available; skipping flanking removal");
                vec![vec![false; n_snps]; n_annot]
            }
        }
    } else {
        vec![vec![false; n_snps]; n_annot]
    };

    let kstar = k * (k + 1) / 2;
    let mut i_mat = Mat::zeros(k, k);

    // Pre-compute per-trait chi-squared (z^2) and mean N
    let chi2_per_trait: Vec<Vec<f64>> = (0..k)
        .map(|j| merged.z[j].iter().map(|z| z * z).collect())
        .collect();
    let mean_n_per_trait: Vec<f64> = (0..k)
        .map(|j| merged.n[j].iter().sum::<f64>() / n_snps as f64)
        .collect();

    // Pre-compute total LD score per SNP (sum across annotations)
    let total_ld: Vec<f64> = (0..n_snps)
        .map(|i| (0..n_annot).map(|a| merged.annot_ld[(i, a)]).sum())
        .collect();

    // Per-annotation S matrices
    let mut s_annot: Vec<Mat<f64>> = (0..n_annot).map(|_| Mat::zeros(k, k)).collect();

    // Collect pseudo-values for V: per-annotation, each has kstar vectors of length n_blocks
    let mut all_pseudos_per_annot: Vec<Vec<Vec<f64>>> =
        (0..n_annot).map(|_| Vec::with_capacity(kstar)).collect();
    let mut n_vec = vec![0.0f64; kstar];

    let mut vech_idx = 0;
    for j in 0..k {
        for jj in j..k {
            let (outcome, mean_n) = if j == jj {
                (chi2_per_trait[j].clone(), mean_n_per_trait[j])
            } else {
                let zz: Vec<f64> = merged.z[j]
                    .iter()
                    .zip(merged.z[jj].iter())
                    .map(|(a, b)| a * b)
                    .collect();
                (zz, (mean_n_per_trait[j] * mean_n_per_trait[jj]).sqrt())
            };

            n_vec[vech_idx] = mean_n;

            let w = if j == jj {
                weights::compute_h2_weights(
                    &chi2_per_trait[j],
                    &total_ld,
                    &merged.w_ld,
                    &merged.n[j],
                    m_total,
                )
            } else {
                weights::compute_gcov_weights(
                    &chi2_per_trait[j],
                    &chi2_per_trait[jj],
                    &total_ld,
                    &merged.w_ld,
                    &merged.n[j],
                    &merged.n[jj],
                    m_total,
                )
            };

            if config.rm_flank {
                // rm_flank: run a separate regression per annotation,
                // excluding flanking SNPs around that annotation's SNPs.
                // First run the full regression once for the intercept.
                let reg_full = block_regression_multi(&merged.annot_ld, &outcome, &w, n_blocks, n_annot);
                let intercept = reg_full.coef[n_annot];
                if j == jj {
                    i_mat[(j, j)] = intercept;
                } else {
                    i_mat[(j, jj)] = intercept;
                    i_mat[(jj, j)] = intercept;
                }

                for a in 0..n_annot {
                    let mask = &flank_masks[a];
                    // Build filtered data for this annotation
                    let kept: Vec<usize> = (0..n_snps).filter(|i| !mask[*i]).collect();
                    let n_kept = kept.len();
                    if n_kept < n_annot + 2 {
                        // Not enough SNPs after flanking removal
                        all_pseudos_per_annot[a].push(vec![0.0; n_blocks.min(n_kept.max(1))]);
                        continue;
                    }

                    let filt_annot = Mat::from_fn(n_kept, n_annot, |i, col| merged.annot_ld[(kept[i], col)]);
                    let filt_outcome: Vec<f64> = kept.iter().map(|&i| outcome[i]).collect();
                    let filt_w: Vec<f64> = kept.iter().map(|&i| w[i]).collect();
                    let filt_blocks = n_blocks.min(n_kept);

                    let reg_a = block_regression_multi(&filt_annot, &filt_outcome, &filt_w, filt_blocks, n_annot);
                    let tau = reg_a.coef[a];
                    let contribution = tau * m_annot[a] / mean_n;

                    if j == jj {
                        s_annot[a][(j, j)] = contribution;
                    } else {
                        s_annot[a][(j, jj)] = contribution;
                        s_annot[a][(jj, j)] = contribution;
                    }

                    let jk = jackknife_multi(&reg_a, n_annot);
                    all_pseudos_per_annot[a].push(jk.pseudo_coefs[a].clone());
                }
            } else {
                // Standard: single regression for all annotations
                let reg = block_regression_multi(&merged.annot_ld, &outcome, &w, n_blocks, n_annot);

                let intercept = reg.coef[n_annot];
                if j == jj {
                    i_mat[(j, j)] = intercept;
                } else {
                    i_mat[(j, jj)] = intercept;
                    i_mat[(jj, j)] = intercept;
                }

                for a in 0..n_annot {
                    let tau = reg.coef[a];
                    let contribution = tau * m_annot[a] / mean_n;

                    if j == jj {
                        s_annot[a][(j, j)] = contribution;
                    } else {
                        s_annot[a][(j, jj)] = contribution;
                        s_annot[a][(jj, j)] = contribution;
                    }
                }

                let jk = jackknife_multi(&reg, n_annot);
                for (a, pseudos) in all_pseudos_per_annot.iter_mut().enumerate().take(n_annot) {
                    pseudos.push(jk.pseudo_coefs[a].clone());
                }
            }

            vech_idx += 1;
        }
    }

    // Construct V matrices per annotation from jackknife pseudo-values
    let v_annot: Vec<Mat<f64>> = (0..n_annot)
        .map(|a| {
            crate::jackknife::construct_v_matrix(
                &all_pseudos_per_annot[a],
                n_blocks,
                &n_vec,
                m_total,
            )
        })
        .collect();

    // Apply liability scale if prevalences provided
    for a in 0..n_annot {
        let mut partial_result = crate::LdscResult {
            s: s_annot[a].to_owned(),
            v: v_annot[a].to_owned(),
            i_mat: i_mat.to_owned(),
            n_vec: n_vec.clone(),
            m: m_total,
        };
        crate::liability::apply_liability_scale(&mut partial_result, sample_prev, pop_prev);
        s_annot[a] = partial_result.s;
    }

    let prop: Vec<f64> = m_annot.iter().map(|&m| m / m_total).collect();

    Ok(StratifiedLdscResult {
        annotations: annot_names.to_vec(),
        s_annot,
        v_annot,
        i_mat,
        m_annot: m_annot.to_vec(),
        m_total,
        prop,
    })
}

// ---------------------------------------------------------------------------
// Internal: data merging
// ---------------------------------------------------------------------------

/// Merged data for stratified LDSC.
struct StratifiedMergedData {
    /// Z-scores per trait (k x n_snps).
    z: Vec<Vec<f64>>,
    /// Sample sizes per trait (k x n_snps).
    n: Vec<Vec<f64>>,
    /// Annotation LD scores (n_snps x n_annot).
    annot_ld: Mat<f64>,
    /// Weight LD scores (n_snps).
    w_ld: Vec<f64>,
    /// Number of SNPs after merging.
    n_snps: usize,
    /// Original LD reference indices for each merged SNP (for chr/bp lookup).
    ld_indices: Vec<usize>,
}

/// Merge trait data with annotation LD scores on SNP ID.
fn merge_stratified_data(
    traits: &[crate::TraitSumstats],
    annot_ld: &Mat<f64>,
    w_ld_scores: &[f64],
    ld_snps: &[String],
) -> Result<StratifiedMergedData> {
    let k = traits.len();
    let n_annot = annot_ld.ncols();

    // Build LD SNP lookup
    let mut ld_map: HashMap<&str, usize> = HashMap::with_capacity(ld_snps.len());
    for (i, snp) in ld_snps.iter().enumerate() {
        ld_map.insert(snp.as_str(), i);
    }

    // Build trait lookups
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

    // Find common SNPs present in all traits and in LD scores
    let mut common_indices: Vec<(usize, Vec<usize>)> = Vec::new();
    for snp in &traits[0].snp {
        let snp_str = snp.as_str();
        if let Some(&ld_idx) = ld_map.get(snp_str) {
            let trait_indices: Option<Vec<usize>> =
                trait_maps.iter().map(|m| m.get(snp_str).copied()).collect();
            if let Some(indices) = trait_indices {
                common_indices.push((ld_idx, indices));
            }
        }
    }

    if common_indices.is_empty() {
        anyhow::bail!("no common SNPs across traits and LD scores");
    }

    let n_snps = common_indices.len();
    let mut z = vec![Vec::with_capacity(n_snps); k];
    let mut n = vec![Vec::with_capacity(n_snps); k];
    let mut merged_annot = Mat::zeros(n_snps, n_annot);
    let mut merged_wld = Vec::with_capacity(n_snps);

    for (snp_i, (ld_idx, trait_indices)) in common_indices.iter().enumerate() {
        for t in 0..k {
            let ti = trait_indices[t];
            z[t].push(traits[t].z[ti]);
            n[t].push(traits[t].n[ti]);
        }
        for a in 0..n_annot {
            merged_annot[(snp_i, a)] = annot_ld[(*ld_idx, a)];
        }
        merged_wld.push(w_ld_scores[*ld_idx]);
    }

    log::info!(
        "s_ldsc merge: {} common SNPs (from {} in trait 0, {} in LD scores)",
        n_snps,
        traits[0].snp.len(),
        ld_snps.len()
    );

    let ld_indices: Vec<usize> = common_indices.iter().map(|(ld_idx, _)| *ld_idx).collect();

    Ok(StratifiedMergedData {
        z,
        n,
        annot_ld: merged_annot,
        w_ld: merged_wld,
        n_snps,
        ld_indices,
    })
}

// ---------------------------------------------------------------------------
// Internal: multi-predictor block regression
// ---------------------------------------------------------------------------

/// Result of multi-predictor block regression.
struct MultiRegressionResult {
    /// Coefficients: [annot_0, annot_1, ..., annot_{p-1}, intercept].
    coef: Vec<f64>,
    /// Per-block XtX matrices (each p x p).
    xtx_blocks: Vec<Mat<f64>>,
    /// Per-block Xty vectors (each length p).
    xty_blocks: Vec<Vec<f64>>,
    /// Total XtX (p x p).
    xtx_total: Mat<f64>,
    /// Total Xty (length p).
    xty_total: Vec<f64>,
}

/// Multi-predictor block regression: y ~ annot_ld_1 + ... + annot_ld_p + intercept.
fn block_regression_multi(
    annot_ld: &Mat<f64>,
    y: &[f64],
    w: &[f64],
    n_blocks: usize,
    n_annot: usize,
) -> MultiRegressionResult {
    let n = y.len();
    let p = n_annot + 1; // annotation predictors + intercept

    // Block boundaries
    let block_bounds: Vec<usize> = (0..=n_blocks).map(|i| (i * n) / n_blocks).collect();
    let actual_blocks = block_bounds.len() - 1;

    let mut xtx_blocks: Vec<Mat<f64>> = (0..actual_blocks).map(|_| Mat::zeros(p, p)).collect();
    let mut xty_blocks: Vec<Vec<f64>> = vec![vec![0.0; p]; actual_blocks];

    for b in 0..actual_blocks {
        let start = block_bounds[b];
        let end = block_bounds[b + 1];

        for i in start..end {
            let wi = w[i];
            let yi = y[i] * wi;

            // Accumulate XtX and Xty
            for a1 in 0..n_annot {
                let x_a1 = annot_ld[(i, a1)] * wi;
                xty_blocks[b][a1] += x_a1 * yi;

                for a2 in a1..n_annot {
                    let x_a2 = annot_ld[(i, a2)] * wi;
                    xtx_blocks[b][(a1, a2)] += x_a1 * x_a2;
                    if a1 != a2 {
                        xtx_blocks[b][(a2, a1)] += x_a1 * x_a2;
                    }
                }
                // Cross-term with intercept column
                let x_int = wi;
                xtx_blocks[b][(a1, n_annot)] += x_a1 * x_int;
                xtx_blocks[b][(n_annot, a1)] += x_a1 * x_int;
            }
            // Intercept x intercept
            let x_int = wi;
            xtx_blocks[b][(n_annot, n_annot)] += x_int * x_int;
            xty_blocks[b][n_annot] += x_int * yi;
        }
    }

    // Sum totals
    let mut xtx_total = Mat::zeros(p, p);
    let mut xty_total = vec![0.0; p];
    for b in 0..actual_blocks {
        for i in 0..p {
            xty_total[i] += xty_blocks[b][i];
            for j in 0..p {
                xtx_total[(i, j)] += xtx_blocks[b][(i, j)];
            }
        }
    }

    // Solve: coef = XtX^{-1} Xty
    let coef = solve_system(&xtx_total, &xty_total);

    MultiRegressionResult {
        coef,
        xtx_blocks,
        xty_blocks,
        xtx_total,
        xty_total,
    }
}

/// Solve a p x p linear system using LU decomposition.
fn solve_system(xtx: &Mat<f64>, xty: &[f64]) -> Vec<f64> {
    let p = xty.len();
    let rhs = Mat::from_fn(p, 1, |i, _| xty[i]);
    let solution = xtx.partial_piv_lu().solve(&rhs);
    (0..p).map(|i| solution[(i, 0)]).collect()
}

// ---------------------------------------------------------------------------
// Internal: jackknife for multi-predictor regression
// ---------------------------------------------------------------------------

/// Result of multi-predictor jackknife.
struct MultiJackknifeResult {
    /// Per-coefficient pseudo-values: n_annot vectors, each of length n_blocks.
    /// (We only track annotation coefficients, not the intercept.)
    pseudo_coefs: Vec<Vec<f64>>,
}

/// Delete-one-block jackknife for multi-predictor regression.
fn jackknife_multi(result: &MultiRegressionResult, n_annot: usize) -> MultiJackknifeResult {
    let n_blocks = result.xtx_blocks.len();
    let nb = n_blocks as f64;
    let p = n_annot + 1;

    let mut pseudo_coefs: Vec<Vec<f64>> =
        (0..n_annot).map(|_| Vec::with_capacity(n_blocks)).collect();

    for b in 0..n_blocks {
        // Leave-one-out: subtract block b from totals
        let mut xtx_loo = result.xtx_total.to_owned();
        let mut xty_loo = result.xty_total.clone();

        for i in 0..p {
            xty_loo[i] -= result.xty_blocks[b][i];
            for j in 0..p {
                xtx_loo[(i, j)] -= result.xtx_blocks[b][(i, j)];
            }
        }

        let coef_loo = solve_system(&xtx_loo, &xty_loo);

        for a in 0..n_annot {
            let pseudo = nb * result.coef[a] - (nb - 1.0) * coef_loo[a];
            pseudo_coefs[a].push(pseudo);
        }
    }

    MultiJackknifeResult { pseudo_coefs }
}

// ---------------------------------------------------------------------------
// JSON serialization helpers
// ---------------------------------------------------------------------------

impl StratifiedLdscResult {
    /// Serialize to a JSON string.
    ///
    /// Since `Mat<f64>` does not implement `Serialize`, this manually builds
    /// the JSON representation with matrices as nested arrays.
    pub fn to_json_string(&self) -> Result<String> {
        use std::fmt::Write;

        let mut out = String::from("{\n");

        // annotations
        writeln!(out, "  \"annotations\": {:?},", self.annotations)?;

        // m_annot
        writeln!(out, "  \"m_annot\": {:?},", self.m_annot)?;

        // m_total
        writeln!(out, "  \"m_total\": {},", self.m_total)?;

        // prop
        writeln!(out, "  \"prop\": {:?},", self.prop)?;

        // i_mat
        writeln!(out, "  \"i_mat\": {},", mat_to_json(&self.i_mat))?;

        // s_annot
        write!(out, "  \"s_annot\": [")?;
        for (i, s) in self.s_annot.iter().enumerate() {
            if i > 0 {
                write!(out, ", ")?;
            }
            write!(out, "{}", mat_to_json(s))?;
        }
        writeln!(out, "],")?;

        // v_annot
        write!(out, "  \"v_annot\": [")?;
        for (i, v) in self.v_annot.iter().enumerate() {
            if i > 0 {
                write!(out, ", ")?;
            }
            write!(out, "{}", mat_to_json(v))?;
        }
        writeln!(out, "]")?;

        out.push('}');
        Ok(out)
    }
}

/// Compute per-annotation flanking masks.
///
/// For each annotation, finds SNPs where the annotation value > 0 (annotation SNPs),
/// then marks all SNPs within `flank_kb` kilobases on the same chromosome as flanking.
/// Returns a Vec of bool masks (one per annotation), where `true` = exclude this SNP.
fn compute_flank_masks(
    annot_ld: &Mat<f64>,
    ld_indices: &[usize],
    chr: &[u32],
    bp: &[u64],
    n_annot: usize,
    flank_kb: u64,
) -> Vec<Vec<bool>> {
    let n_snps = ld_indices.len();
    let flank_bp = flank_kb * 1000;

    (0..n_annot)
        .map(|a| {
            // Find annotation SNPs (non-zero annotation value)
            let annot_snp_positions: Vec<(u32, u64)> = (0..n_snps)
                .filter(|&i| annot_ld[(i, a)] > 0.0)
                .map(|i| {
                    let orig_idx = ld_indices[i];
                    (chr[orig_idx], bp[orig_idx])
                })
                .collect();

            // Mark flanking SNPs
            let mut mask = vec![false; n_snps];
            for i in 0..n_snps {
                let orig_idx = ld_indices[i];
                let snp_chr = chr[orig_idx];
                let snp_bp = bp[orig_idx];

                // Check if this SNP is within flank_bp of any annotation SNP on same chr
                // (but is not itself an annotation SNP)
                if annot_ld[(i, a)] <= 0.0 {
                    for &(a_chr, a_bp) in &annot_snp_positions {
                        if snp_chr == a_chr && snp_bp.abs_diff(a_bp) <= flank_bp {
                            mask[i] = true;
                            break;
                        }
                    }
                }
            }
            mask
        })
        .collect()
}

/// Convert a faer Mat to a JSON-style nested array string.
fn mat_to_json(m: &Mat<f64>) -> String {
    let rows: Vec<String> = (0..m.nrows())
        .map(|i| {
            let cols: Vec<String> = (0..m.ncols()).map(|j| format!("{}", m[(i, j)])).collect();
            format!("[{}]", cols.join(", "))
        })
        .collect();
    format!("[{}]", rows.join(", "))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::TraitSumstats;

    /// Smoke test: single trait, two annotations, verifies output structure.
    #[test]
    fn test_s_ldsc_smoke() {
        let n_snps = 100;
        let n_annot = 2;

        // Create synthetic trait data
        let trait_data = TraitSumstats {
            snp: (0..n_snps).map(|i| format!("rs{i}")).collect(),
            z: (0..n_snps).map(|i| 0.5 + 0.01 * i as f64).collect(),
            n: vec![50000.0; n_snps],
            a1: vec!["A".to_string(); n_snps],
            a2: vec!["G".to_string(); n_snps],
        };

        // Create annotation LD scores
        let annot_ld = Mat::from_fn(n_snps, n_annot, |i, a| {
            if a == 0 {
                10.0 + 0.1 * i as f64
            } else {
                5.0 + 0.05 * i as f64
            }
        });

        let w_ld: Vec<f64> = vec![15.0; n_snps];
        let ld_snps: Vec<String> = (0..n_snps).map(|i| format!("rs{i}")).collect();
        let annot_names = vec!["annot1".to_string(), "annot2".to_string()];
        let m_annot = vec![500000.0, 300000.0];

        let config = StratifiedLdscConfig { n_blocks: 10, rm_flank: false, flank_kb: 500 };

        let result = s_ldsc(
            &[trait_data],
            &[None],
            &[None],
            &annot_ld,
            &w_ld,
            &ld_snps,
            &annot_names,
            &m_annot,
            &config,
            None,
            None,
        )
        .expect("s_ldsc should succeed");

        // Check structure
        assert_eq!(result.annotations.len(), 2);
        assert_eq!(result.s_annot.len(), 2);
        assert_eq!(result.v_annot.len(), 2);
        assert_eq!(result.i_mat.nrows(), 1);
        assert_eq!(result.i_mat.ncols(), 1);
        assert_eq!(result.m_annot.len(), 2);
        assert!((result.m_total - 800000.0).abs() < 1e-6);
        assert_eq!(result.prop.len(), 2);

        // Each S_annot should be 1x1
        for s in &result.s_annot {
            assert_eq!(s.nrows(), 1);
            assert_eq!(s.ncols(), 1);
        }

        // Each V_annot should be 1x1
        for v in &result.v_annot {
            assert_eq!(v.nrows(), 1);
            assert_eq!(v.ncols(), 1);
        }
    }

    #[test]
    fn test_solve_system() {
        // [2 1; 1 3] * [x; y] = [5; 7]
        let xtx = Mat::from_fn(2, 2, |i, j| [[2.0, 1.0], [1.0, 3.0]][i][j]);
        let xty = [5.0, 7.0];
        let sol = solve_system(&xtx, &xty);
        assert!((sol[0] - 1.6).abs() < 1e-10);
        assert!((sol[1] - 1.8).abs() < 1e-10);
    }

    #[test]
    fn test_block_regression_multi_simple() {
        // Simple test: 2 annotations with independent variation
        let n = 50;
        let n_annot = 2;
        let annot_ld = Mat::from_fn(n, n_annot, |i, a| {
            if a == 0 {
                (i + 1) as f64
            } else {
                // Use i^2 / n to create an independent predictor
                ((i as f64) * (i as f64)) / n as f64 + 1.0
            }
        });
        // y = 2 * annot_0 + 3 * annot_1 + 1 (intercept)
        let y: Vec<f64> = (0..n)
            .map(|i| 2.0 * annot_ld[(i, 0)] + 3.0 * annot_ld[(i, 1)] + 1.0)
            .collect();
        let w = vec![1.0; n];

        let result = block_regression_multi(&annot_ld, &y, &w, 5, n_annot);

        assert!(
            (result.coef[0] - 2.0).abs() < 1e-6,
            "coef[0] should be 2.0, got {}",
            result.coef[0]
        );
        assert!(
            (result.coef[1] - 3.0).abs() < 1e-6,
            "coef[1] should be 3.0, got {}",
            result.coef[1]
        );
        assert!(
            (result.coef[2] - 1.0).abs() < 1e-6,
            "intercept should be 1.0, got {}",
            result.coef[2]
        );
    }

    #[test]
    fn test_json_serialization() {
        let result = StratifiedLdscResult {
            annotations: vec!["a1".to_string(), "a2".to_string()],
            s_annot: vec![faer::mat![[0.5]], faer::mat![[0.3]]],
            v_annot: vec![faer::mat![[0.01]], faer::mat![[0.02]]],
            i_mat: faer::mat![[1.0]],
            m_annot: vec![500000.0, 300000.0],
            m_total: 800000.0,
            prop: vec![0.625, 0.375],
        };

        let json = result.to_json_string().expect("serialization should work");
        assert!(json.contains("\"annotations\""));
        assert!(json.contains("\"s_annot\""));
        assert!(json.contains("\"v_annot\""));
    }
}
