use std::collections::HashMap;

use anyhow::Result;
use statrs::distribution::{ChiSquared, ContinuousCDF};

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
    let chi2_dist = ChiSquared::new(1.0).unwrap();
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
        if let Some(info) = rec.info
            && info < config.info_filter
        {
            n_info_filtered += 1;
            continue;
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
        let z = compute_z(effect, p, &chi2_dist);

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
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = sorted[sorted.len() / 2];

    // R uses: if round(median, 0) == 1 → OR
    median.round() == 1.0
}

/// Compute Z-score: sign(effect) * sqrt(qchisq(p, 1, lower.tail=FALSE))
///
/// For very small p-values (< 1e-300), uses log-space approximation.
fn compute_z(effect: f64, p: f64, chi2_dist: &ChiSquared) -> f64 {
    if p == 0.0 {
        // Extreme: return large Z with correct sign
        return if effect >= 0.0 { 37.0 } else { -37.0 };
    }

    // qchisq(p, 1, lower.tail=FALSE) = inverse_survival(p) for chi2(1)
    let chi2_val = if p < 1e-300 {
        // Log-space approximation for tiny p-values
        // For chi2(1), survival function ≈ erfc(sqrt(x/2))/2
        // So for very small p: x ≈ 2 * (erfcinv(2p))^2 ≈ -2 * ln(p)
        -2.0 * p.ln()
    } else {
        // Standard inverse: qchisq(p, 1, lower.tail=FALSE) = inverse_cdf(1-p)
        chi2_dist.inverse_cdf(1.0 - p)
    };

    let z_abs = chi2_val.max(0.0).sqrt();
    if effect >= 0.0 { z_abs } else { -z_abs }
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
        let chi2 = ChiSquared::new(1.0).unwrap();
        let z = compute_z(0.5, 0.05, &chi2);
        assert!(z > 0.0);
        assert!((z - 1.96).abs() < 0.01);
    }

    #[test]
    fn test_compute_z_negative() {
        let chi2 = ChiSquared::new(1.0).unwrap();
        let z = compute_z(-0.5, 0.05, &chi2);
        assert!(z < 0.0);
        assert!((z + 1.96).abs() < 0.01);
    }

    #[test]
    fn test_compute_z_tiny_p() {
        let chi2 = ChiSquared::new(1.0).unwrap();
        let z = compute_z(1.0, 1e-300, &chi2);
        assert!(z > 30.0);
    }

    #[test]
    fn test_compute_z_zero_p() {
        let chi2 = ChiSquared::new(1.0).unwrap();
        let z = compute_z(1.0, 0.0, &chi2);
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
