//! Merge multiple GWAS summary statistics files for multivariate GWAS.
//!
//! Port of R GenomicSEM's `sumstats()`. Reads multiple GWAS files, applies QC
//! filters, aligns alleles against an LD score reference, and outputs a merged
//! tab-delimited file with columns: SNP, A1, A2, MAF, beta.T1, se.T1, ...

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};

use crate::io::gwas_reader;
use crate::munge::allele::{AlleleMatch, alleles_match};

/// Configuration for the sumstats merge pipeline.
#[derive(Debug, Clone)]
pub struct SumstatsConfig {
    pub info_filter: f64,
    pub maf_filter: f64,
    /// Per-trait sample size overrides (None = use file values)
    pub n_overrides: Vec<Option<f64>>,
    /// Per-trait: are betas on logistic scale? (need se.logit transform)
    pub se_logit: Vec<bool>,
    /// Per-trait: are effects from OLS linear probability model?
    pub ols: Vec<bool>,
    /// Per-trait: are effects linear probability betas?
    pub linprob: Vec<bool>,
    /// Keep indels (multi-character alleles)?
    pub keep_indel: bool,
    /// If false (default), remove ambiguous strand SNPs (A/T and C/G pairs).
    pub keep_ambig: bool,
}

impl Default for SumstatsConfig {
    fn default() -> Self {
        Self {
            info_filter: 0.6,
            maf_filter: 0.01,
            n_overrides: Vec::new(),
            se_logit: Vec::new(),
            ols: Vec::new(),
            linprob: Vec::new(),
            keep_indel: false,
            keep_ambig: false,
        }
    }
}

/// A single SNP's merged data across all traits.
struct MergedSnpData {
    snp: String,
    a1: String,
    a2: String,
    maf: f64,
    beta: Vec<f64>,
    se: Vec<f64>,
}

/// Run the sumstats merge pipeline.
///
/// 1. Read LD score reference for allele alignment and MAF
/// 2. Read each GWAS file, apply QC filters
/// 3. Inner join on SNP across all files
/// 4. Align alleles against reference (flip beta signs as needed)
/// 5. Write merged output
pub fn merge_sumstats(
    files: &[&Path],
    ref_file: &Path,
    trait_names: &[String],
    config: &SumstatsConfig,
    out: &Path,
) -> Result<usize> {
    let k = files.len();
    if k == 0 {
        anyhow::bail!("no input files provided");
    }

    // Read reference panel file (e.g., w_hm3.snplist or reference.1000G.maf.0.005.txt)
    // Matches R's behavior: ref is a single file with SNP, A1, A2, MAF columns
    let ref_snps = read_reference_file(ref_file, config.maf_filter)?;
    log::info!("Loaded {} reference SNPs", ref_snps.len());

    // Read and QC each GWAS file
    let mut trait_records: Vec<HashMap<String, QcRecord>> = Vec::with_capacity(k);
    for (i, file) in files.iter().enumerate() {
        let n_override = config.n_overrides.get(i).copied().flatten();
        let records = read_and_qc_gwas(file, &ref_snps, config, n_override, i)?;
        log::info!("  {}: {} SNPs after QC", trait_names[i], records.len());
        trait_records.push(records);
    }

    // Inner join: find SNPs present in ALL traits AND in reference
    let common_snps = find_common_snps(&trait_records, &ref_snps);
    log::info!("Common SNPs across all traits: {}", common_snps.len());

    if common_snps.is_empty() {
        anyhow::bail!("no common SNPs found across all traits and reference");
    }

    // Build merged data with allele alignment
    let merged = build_merged(&common_snps, &trait_records, &ref_snps, k);

    // Write output
    write_merged_sumstats(&merged, trait_names, out)?;
    log::info!(
        "Wrote {} SNPs x {} traits to {}",
        merged.len(),
        k,
        out.display()
    );

    Ok(merged.len())
}

/// Reference SNP with alleles and MAF.
struct RefSnpInfo {
    a1: String,
    a2: String,
    maf: Option<f64>,
}

/// Read reference panel file (like R's `fread(ref)`).
///
/// Expects a tab/whitespace-delimited file with at least SNP column.
/// Optionally A1, A2, MAF columns for allele alignment.
/// Filters by MAF if provided.
fn read_reference_file(path: &Path, maf_filter: f64) -> Result<HashMap<String, RefSnpInfo>> {
    let data = gwas_reader::read_gwas_file(path)
        .with_context(|| format!("failed to read reference file {}", path.display()))?;

    let mut map = HashMap::with_capacity(data.records.len());
    for rec in &data.records {
        // Filter by MAF if available (R: ref <- ref[ref$MAF >= maf.filter, ])
        if let Some(maf) = rec.maf {
            if maf < maf_filter {
                continue;
            }
        }

        map.insert(
            rec.snp.clone(),
            RefSnpInfo {
                a1: rec.a1.as_deref().unwrap_or("").to_uppercase(),
                a2: rec.a2.as_deref().unwrap_or("").to_uppercase(),
                maf: rec.maf,
            },
        );
    }

    Ok(map)
}

/// A QC-passed record from a single GWAS file.
struct QcRecord {
    a1: String,
    a2: String,
    beta: f64,
    se: f64,
    maf: Option<f64>,
}

/// Read a GWAS file and apply QC filters.
fn read_and_qc_gwas(
    path: &Path,
    ref_snps: &HashMap<String, RefSnpInfo>,
    config: &SumstatsConfig,
    n_override: Option<f64>,
    trait_idx: usize,
) -> Result<HashMap<String, QcRecord>> {
    let data = gwas_reader::read_gwas_file(path)?;
    let mut records = HashMap::new();

    for rec in data.records {
        // Must be in reference
        if !ref_snps.contains_key(&rec.snp) {
            continue;
        }

        // Must have effect and SE
        let (Some(effect), Some(se)) = (rec.effect, rec.se) else {
            // Try to derive from Z and N if available
            if let (Some(z), Some(n_val)) = (rec.z, rec.n.or(n_override)) {
                if n_val > 0.0 {
                    let beta = z / n_val.sqrt();
                    let se_val = 1.0 / n_val.sqrt();
                    let a1 = rec.a1.as_deref().unwrap_or("").to_uppercase();
                    let a2 = rec.a2.as_deref().unwrap_or("").to_uppercase();
                    // Accept even without A2 (common for .sumstats.gz from simLDSC)
                    if !a1.is_empty() {
                        if a2.is_empty()
                            || (valid_allele_pair(&a1, &a2, config.keep_indel)
                                && (config.keep_ambig || !is_ambiguous_snp(&a1, &a2)))
                        {
                            records.insert(
                                rec.snp,
                                QcRecord {
                                    a1,
                                    a2,
                                    beta,
                                    se: se_val,
                                    maf: rec.maf,
                                },
                            );
                        }
                    }
                }
                continue;
            }
            continue;
        };

        let se = if se <= 0.0 || !se.is_finite() {
            continue;
        } else {
            se
        };

        let effect = if !effect.is_finite() {
            continue;
        } else {
            effect
        };

        // A1 required, A2 optional (munged .sumstats.gz may lack A2)
        let a1 = match rec.a1.as_ref() {
            Some(a) if !a.is_empty() => a.to_uppercase(),
            _ => continue,
        };
        let a2 = rec.a2.as_deref().unwrap_or("").to_uppercase();

        // Allele validation (skip if A2 missing — trust the munged file)
        if !a2.is_empty() && !valid_allele_pair(&a1, &a2, config.keep_indel) {
            continue;
        }

        // Ambiguous strand SNP filtering (skip check if A2 missing)
        if !a2.is_empty() && !config.keep_ambig && is_ambiguous_snp(&a1, &a2) {
            continue;
        }

        // INFO filter
        if rec.info.is_some_and(|info| info < config.info_filter) {
            continue;
        }

        // MAF filter
        let maf = rec.maf.map(|m| if m > 0.5 { 1.0 - m } else { m });
        if maf.is_some_and(|maf_val| maf_val < config.maf_filter) {
            continue;
        }

        // Detect and convert OR to log(OR)
        // Skip OR auto-detection if se_logit or ols is set for this trait
        // (effects are already on the correct scale)
        let skip_or_detect = config.se_logit.get(trait_idx).copied().unwrap_or(false)
            || config.ols.get(trait_idx).copied().unwrap_or(false);
        let beta = if !skip_or_detect && is_or_value(effect) {
            if effect <= 0.0 {
                continue;
            }
            effect.ln()
        } else {
            effect
        };

        records.insert(
            rec.snp,
            QcRecord {
                a1,
                a2,
                beta,
                se,
                maf,
            },
        );
    }

    Ok(records)
}

/// Check if allele pair is valid.
fn valid_allele_pair(a1: &str, a2: &str, keep_indel: bool) -> bool {
    if keep_indel {
        !a1.is_empty() && !a2.is_empty()
    } else {
        is_single_nucleotide(a1) && is_single_nucleotide(a2)
    }
}

fn is_single_nucleotide(a: &str) -> bool {
    matches!(a, "A" | "C" | "G" | "T")
}

/// Check if a SNP has ambiguous strand (A/T or C/G pairs).
fn is_ambiguous_snp(a1: &str, a2: &str) -> bool {
    matches!((a1, a2), ("A", "T") | ("T", "A") | ("C", "G") | ("G", "C"))
}

/// Rough heuristic: if |effect| is close to 1, it's likely an OR.
fn is_or_value(effect: f64) -> bool {
    // Effects that round to 1 are likely OR
    effect.round() == 1.0 && effect > 0.0
}

/// Find SNPs present in ALL trait HashMaps and in the reference.
fn find_common_snps(
    trait_records: &[HashMap<String, QcRecord>],
    ref_snps: &HashMap<String, RefSnpInfo>,
) -> Vec<String> {
    if trait_records.is_empty() {
        return Vec::new();
    }

    // Start with reference SNPs and filter to those present in all traits
    ref_snps
        .keys()
        .filter(|snp| trait_records.iter().all(|m| m.contains_key(snp.as_str())))
        .cloned()
        .collect()
}

/// Build merged SNP data, aligning alleles across traits against the reference panel.
fn build_merged(
    common_snps: &[String],
    trait_records: &[HashMap<String, QcRecord>],
    ref_snps: &HashMap<String, RefSnpInfo>,
    k: usize,
) -> Vec<MergedSnpData> {
    let mut merged = Vec::with_capacity(common_snps.len());

    for snp in common_snps {
        let ref_info = &ref_snps[snp];
        let ref_a1 = &ref_info.a1;
        let ref_a2 = &ref_info.a2;

        let mut betas = Vec::with_capacity(k);
        let mut ses = Vec::with_capacity(k);
        let mut skip = false;
        let mut maf_val = ref_info.maf.unwrap_or(0.0);

        for records in trait_records.iter() {
            let rec = &records[snp];

            // Align alleles against reference panel
            if ref_a1.is_empty() && ref_a2.is_empty() {
                // No reference alleles — use as-is (trust the munged file)
                betas.push(rec.beta);
                ses.push(rec.se);
            } else if rec.a1.is_empty() && rec.a2.is_empty() {
                // No trait alleles (e.g., simLDSC output) — use as-is
                betas.push(rec.beta);
                ses.push(rec.se);
            } else {
                let amatch = alleles_match(&rec.a1, &rec.a2, ref_a1, ref_a2);
                match amatch {
                    AlleleMatch::Match => {
                        betas.push(rec.beta);
                        ses.push(rec.se);
                    }
                    AlleleMatch::Flipped => {
                        betas.push(-rec.beta);
                        ses.push(rec.se);
                    }
                    AlleleMatch::NoMatch => {
                        skip = true;
                        break;
                    }
                }
            }

            // Use MAF from reference or trait
            if let Some(m) = rec.maf.filter(|_| maf_val == 0.0) {
                maf_val = m;
            }
        }

        if skip {
            continue;
        }

        merged.push(MergedSnpData {
            snp: snp.clone(),
            a1: ref_a1.clone(),
            a2: ref_a2.clone(),
            maf: maf_val,
            beta: betas,
            se: ses,
        });
    }

    merged
}

/// Write merged sumstats to a tab-delimited file.
fn write_merged_sumstats(
    merged: &[MergedSnpData],
    trait_names: &[String],
    path: &Path,
) -> Result<()> {
    let file = File::create(path).with_context(|| format!("cannot create {}", path.display()))?;
    let mut w = BufWriter::new(file);

    // Header
    let mut header = String::from("SNP\tA1\tA2\tMAF");
    for name in trait_names {
        header.push_str(&format!("\tbeta.{name}\tse.{name}"));
    }
    writeln!(w, "{header}")?;

    // Data
    for snp in merged {
        write!(w, "{}\t{}\t{}\t{:.6}", snp.snp, snp.a1, snp.a2, snp.maf)?;
        for t in 0..trait_names.len() {
            write!(w, "\t{:.6}\t{:.6}", snp.beta[t], snp.se[t])?;
        }
        writeln!(w)?;
    }

    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_ambiguous_snp() {
        assert!(is_ambiguous_snp("A", "T"));
        assert!(is_ambiguous_snp("T", "A"));
        assert!(is_ambiguous_snp("C", "G"));
        assert!(is_ambiguous_snp("G", "C"));
        assert!(!is_ambiguous_snp("A", "G"));
        assert!(!is_ambiguous_snp("A", "C"));
        assert!(!is_ambiguous_snp("T", "G"));
        assert!(!is_ambiguous_snp("T", "C"));
    }

    #[test]
    fn test_valid_allele_pair() {
        assert!(valid_allele_pair("A", "G", false));
        assert!(valid_allele_pair("C", "T", false));
        assert!(!valid_allele_pair("AC", "G", false));
        assert!(valid_allele_pair("AC", "G", true)); // indels allowed
        assert!(!valid_allele_pair("", "G", false));
    }

    #[test]
    fn test_find_common_snps() {
        let mut m1 = HashMap::new();
        let mut m2 = HashMap::new();

        let rec = || QcRecord {
            a1: "A".into(),
            a2: "G".into(),
            beta: 0.1,
            se: 0.01,
            maf: Some(0.3),
        };

        m1.insert("rs1".into(), rec());
        m1.insert("rs2".into(), rec());
        m1.insert("rs3".into(), rec());

        m2.insert("rs1".into(), rec());
        m2.insert("rs3".into(), rec());
        m2.insert("rs4".into(), rec());

        let ref_info = || RefSnpInfo {
            a1: "A".into(),
            a2: "G".into(),
            maf: None,
        };
        let mut ref_snps = HashMap::new();
        ref_snps.insert("rs1".into(), ref_info());
        ref_snps.insert("rs2".into(), ref_info());
        ref_snps.insert("rs3".into(), ref_info());
        ref_snps.insert("rs4".into(), ref_info());

        let mut common = find_common_snps(&[m1, m2], &ref_snps);
        common.sort();
        assert_eq!(common, vec!["rs1", "rs3"]);
    }
}
