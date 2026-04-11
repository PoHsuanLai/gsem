use std::collections::HashMap;

use anyhow::Result;
use statrs::distribution::{ContinuousCDF, Normal};

use crate::io::gwas_reader::{GwasRecord, MungedRecord};

use super::allele::{AlleleMatch, alleles_match};
use super::{MungeConfig, RefSnp};

/// Run the full QC pipeline on GWAS records, producing munged output.
///
/// Pipeline steps (matching R GenomicSEM's munge_main.R):
/// 1. Handle special N columns (NEFFDIV2/NEFF_HALF doubling)
/// 2. Validate alleles are A/C/G/T
/// 3. Merge with HapMap3 reference
/// 4. Remove missing P/effect
/// 5. Detect OR vs beta (median effect ≈ 1 → OR)
/// 6. Allele flipping against reference
/// 7. Remove non-matching alleles
/// 8. Filter by INFO and MAF
/// 9. Compute Z-scores: sign(effect) * sqrt(qchisq(p, 1, lower=FALSE))
pub fn run_qc_pipeline(
    records: Vec<GwasRecord>,
    reference: &HashMap<String, RefSnp>,
    config: &MungeConfig,
) -> Result<Vec<MungedRecord>> {
    let normal = Normal::standard();
    let mut output = Vec::new();
    let mut n_no_ref = 0usize;
    let mut n_missing = 0usize;
    let mut n_allele_mismatch = 0usize;
    let mut n_info_filtered = 0usize;
    let mut n_maf_filtered = 0usize;
    let mut n_bad_alleles = 0usize;

    // Detect if effects are OR (median ≈ 1) or beta (median ≈ 0)
    let is_or = detect_or(&records);
    if is_or {
        log::info!("Detected OR format (median effect ≈ 1), converting to log(OR)");
    }

    for mut rec in records {
        // Apply N override
        if let Some(n_override) = config.n_override {
            rec.n = Some(n_override);
        }

        // Must have P and effect
        let (Some(mut effect), Some(p)) = (rec.effect, rec.p) else {
            n_missing += 1;
            continue;
        };

        // Validate P
        if !(0.0..=1.0).contains(&p) || p.is_nan() {
            n_missing += 1;
            continue;
        }

        // Must have N
        let n = rec.n.unwrap_or(0.0);
        if n <= 0.0 {
            n_missing += 1;
            continue;
        }

        // Validate alleles
        let (Some(a1_str), Some(a2_str)) = (rec.a1.as_ref(), rec.a2.as_ref()) else {
            n_bad_alleles += 1;
            continue;
        };
        let a1 = a1_str.to_uppercase();
        let a2 = a2_str.to_uppercase();
        if !is_valid_allele(&a1) || !is_valid_allele(&a2) {
            n_bad_alleles += 1;
            continue;
        }

        // Merge with reference
        let Some(ref_snp) = reference.get(&rec.snp) else {
            n_no_ref += 1;
            continue;
        };

        // INFO filter
        #[allow(clippy::collapsible_if)]
        if let Some(info) = rec.info {
            if info < config.info_filter {
                n_info_filtered += 1;
                continue;
            }
        }

        // MAF filter
        if let Some(mut maf) = rec.maf {
            // Ensure MAF <= 0.5
            if maf > 0.5 {
                maf = 1.0 - maf;
            }
            if maf < config.maf_filter {
                n_maf_filtered += 1;
                continue;
            }
        }

        // Convert OR to log(OR)
        if is_or {
            if effect <= 0.0 {
                n_missing += 1;
                continue;
            }
            effect = effect.ln();
        }

        // Allele matching and flipping
        let amatch = alleles_match(&a1, &a2, &ref_snp.a1, &ref_snp.a2);
        match amatch {
            AlleleMatch::Match => {}
            AlleleMatch::Flipped => {
                effect = -effect;
            }
            AlleleMatch::NoMatch => {
                n_allele_mismatch += 1;
                continue;
            }
        }

        // Compute Z-score: sign(effect) * sqrt(qchisq(p, df=1, lower.tail=FALSE))
        let z = compute_z(effect, p, &normal);

        output.push(MungedRecord {
            snp: rec.snp,
            n,
            z,
            a1: ref_snp.a1.clone(),
            a2: ref_snp.a2.clone(),
        });
    }

    log::info!(
        "QC summary: {} no ref, {} missing P/effect/N, {} bad alleles, {} allele mismatch, {} INFO filtered, {} MAF filtered",
        n_no_ref,
        n_missing,
        n_bad_alleles,
        n_allele_mismatch,
        n_info_filtered,
        n_maf_filtered
    );

    Ok(output)
}

/// Detect if effects are odds ratios (median ≈ 1) or betas (median ≈ 0).
fn detect_or(records: &[GwasRecord]) -> bool {
    let effects: Vec<f64> = records
        .iter()
        .filter_map(|r| r.effect)
        .filter(|e| e.is_finite())
        .collect();

    if effects.is_empty() {
        return false;
    }

    let mut sorted = effects;
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median = sorted[sorted.len() / 2];

    // R uses: if round(median, 0) == 1 → OR
    median.round() == 1.0
}

/// Compute Z-score: sign(effect) * sqrt(qchisq(p, 1, lower.tail=FALSE)).
///
/// For chi-square with 1 df the identity
/// `sqrt(qchisq(p, 1, lower.tail=FALSE)) = qnorm(1 - p/2)`
/// holds for all `p` in `[0, 1]`, and `Normal::inverse_cdf` is
/// numerically stable across the full range. The earlier chi-square
/// path hit a NaN regression in `statrs::distribution::ChiSquared`
/// around `p ≳ 0.97` which silently collapsed large-p SNPs to `Z = 0`.
///
/// When `effect` is exactly zero, the Z is zero — this mirrors R's
/// `sign(0) == 0` so records with a literal `Effect = 0` entry from
/// upstream pipelines produce `Z = 0` instead of an arbitrary-sign
/// near-zero value.
fn compute_z(effect: f64, p: f64, normal: &Normal) -> f64 {
    if effect == 0.0 {
        return 0.0;
    }
    if p == 0.0 {
        // Extreme: return large Z with correct sign.
        return if effect > 0.0 { 37.0 } else { -37.0 };
    }

    let z_abs = normal.inverse_cdf(1.0 - p / 2.0);
    if !z_abs.is_finite() {
        return if effect > 0.0 { 37.0 } else { -37.0 };
    }

    if effect > 0.0 { z_abs } else { -z_abs }
}

/// Check if an allele string is a valid single nucleotide (A, C, G, T).
fn is_valid_allele(a: &str) -> bool {
    matches!(a, "A" | "C" | "G" | "T")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_z_positive() {
        let normal = Normal::standard();
        let z = compute_z(0.5, 0.05, &normal);
        assert!(z > 0.0);
        assert!((z - 1.96).abs() < 0.01);
    }

    #[test]
    fn test_compute_z_negative() {
        let normal = Normal::standard();
        let z = compute_z(-0.5, 0.05, &normal);
        assert!(z < 0.0);
        assert!((z + 1.96).abs() < 0.01);
    }

    #[test]
    fn test_compute_z_large_p() {
        // Regression: statrs 0.18's ChiSquared::inverse_cdf returned NaN
        // for p ≳ 0.97, collapsing Z to 0. The Normal-based path must
        // survive the full range and match R's qnorm(1 - p/2).
        let normal = Normal::standard();
        for &p in &[0.9722_f64, 0.98, 0.99, 0.999] {
            let z = compute_z(-0.001, p, &normal);
            assert!(z.is_finite(), "z non-finite at p={p}");
            // R: qnorm(1 - p/2) for these p: 0.0348, 0.0251, 0.0125, 0.00125
            let expected = -match p {
                p if (p - 0.9722).abs() < 1e-9 => 0.034849,
                p if (p - 0.98).abs() < 1e-9 => 0.025069,
                p if (p - 0.99).abs() < 1e-9 => 0.012533,
                p if (p - 0.999).abs() < 1e-9 => 0.001253,
                _ => unreachable!(),
            };
            assert!(
                (z - expected).abs() < 1e-4,
                "p={p} z={z} expected≈{expected}"
            );
        }
    }

    #[test]
    fn test_compute_z_tiny_p() {
        let normal = Normal::standard();
        let z = compute_z(1.0, 1e-300, &normal);
        assert!(z > 30.0);
    }

    #[test]
    fn test_compute_z_zero_p() {
        let normal = Normal::standard();
        let z = compute_z(1.0, 0.0, &normal);
        assert_eq!(z, 37.0);
    }

    #[test]
    fn test_detect_or_beta() {
        let records: Vec<GwasRecord> = (0..100)
            .map(|_| GwasRecord {
                snp: "rs1".to_string(),
                effect: Some(0.05),
                p: Some(0.5),
                a1: None,
                a2: None,
                se: None,
                n: None,
                z: None,
                info: None,
                maf: None,
            })
            .collect();
        assert!(!detect_or(&records));
    }

    #[test]
    fn test_detect_or_or() {
        let records: Vec<GwasRecord> = (0..100)
            .map(|_| GwasRecord {
                snp: "rs1".to_string(),
                effect: Some(1.05),
                p: Some(0.5),
                a1: None,
                a2: None,
                se: None,
                n: None,
                z: None,
                info: None,
                maf: None,
            })
            .collect();
        assert!(detect_or(&records));
    }

    #[test]
    fn test_is_valid_allele() {
        assert!(is_valid_allele("A"));
        assert!(is_valid_allele("C"));
        assert!(is_valid_allele("G"));
        assert!(is_valid_allele("T"));
        assert!(!is_valid_allele("N"));
        assert!(!is_valid_allele("AC"));
        assert!(!is_valid_allele(""));
    }
}
