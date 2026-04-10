//! Merge multiple GWAS summary statistics files for multivariate GWAS.
//!
//! Port of R GenomicSEM's `sumstats()`. Reads multiple GWAS files, applies QC
//! filters, aligns alleles against an LD score reference, and outputs a merged
//! tab-delimited file with columns: SNP, A1, A2, MAF, beta.T1, se.T1, ...

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use rayon::prelude::*;

use crate::io::column_detect;
use crate::io::gwas_reader::open_file_reader;
use crate::munge::allele::{AlleleMatch, alleles_match};

/// Configuration for the sumstats merge pipeline.
#[derive(Debug, Clone)]
pub struct SumstatsConfig {
    /// Minimum imputation INFO score (SNPs below this are removed)
    pub info_filter: f64,
    /// Minimum minor allele frequency (SNPs below this are removed)
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
    /// Per-trait column name overrides for the effect (beta) column.
    /// If set, the i-th entry overrides auto-detection for trait i.
    pub beta_overrides: Vec<Option<String>>,
    /// If true, apply `maf_filter` directly to each GWAS file's MAF/FRQ column
    /// (in addition to the reference-based MAF filter).
    pub direct_filter: bool,
    /// Number of threads for the local rayon pool. `None` or `Some(0)` = rayon default.
    pub num_threads: Option<usize>,
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
            beta_overrides: Vec::new(),
            direct_filter: false,
            num_threads: None,
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

    // Build a local rayon pool so callers control parallelism per-invocation
    // and concurrent calls don't share global state.
    let mut builder = rayon::ThreadPoolBuilder::new();
    if let Some(n) = config.num_threads
        && n > 0
    {
        builder = builder.num_threads(n);
    }
    let pool = builder.build().expect("failed to build rayon thread pool");

    // Read the reference file and all GWAS files in parallel. The reference
    // is needed to filter GWAS rows, so we first load it, then fan out the
    // GWAS reads. Decompression and parsing dominate here, so running each
    // file on its own worker thread gives roughly a k× speedup until we hit
    // memory bandwidth.
    let (ref_snps, ref_order) = pool.install(|| read_reference_file(ref_file, config.maf_filter))?;
    log::info!("Loaded {} reference SNPs", ref_snps.len());

    // Read and QC each GWAS file in parallel.
    let trait_records: Vec<HashMap<String, QcRecord>> = pool.install(|| {
        files
            .par_iter()
            .enumerate()
            .map(|(i, file)| {
                let records = read_and_qc_gwas(file, &ref_snps, config, i)?;
                log::info!("  {}: {} SNPs after QC", trait_names[i], records.len());
                Ok(records)
            })
            .collect::<Result<Vec<_>>>()
    })?;

    // Inner join: find SNPs present in ALL traits AND in reference.
    // We iterate `ref_order` (reference-file order) rather than HashMap
    // keys so the output is deterministic and matches R's behavior, which
    // starts from the reference and successively inner-joins each trait.
    let common_snps = find_common_snps(&trait_records, &ref_order);
    log::info!("Common SNPs across all traits: {}", common_snps.len());

    if common_snps.is_empty() {
        anyhow::bail!("no common SNPs found across all traits and reference");
    }

    // Build merged data with allele alignment
    let merged = build_merged(&common_snps, &trait_records, &ref_snps, k);

    // Check for missing MAF: if all SNPs have MAF=0 it means neither the reference
    // nor GWAS files had MAF data. Error instead of silently producing MAF=0.
    let n_zero_maf = merged.iter().filter(|s| s.maf == 0.0).count();
    if n_zero_maf == merged.len() && !merged.is_empty() {
        anyhow::bail!(
            "All {} SNPs have MAF=0. Neither the reference file nor the GWAS files \
             contain a MAF column. A reference file with MAF data is required.",
            merged.len()
        );
    }

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

/// Iterator over the fields of a GWAS/reference line.
///
/// Tab-delimited is the common case; we fall back to whitespace if no tab is
/// present (the 1000G reference file ships space-separated). This enum
/// avoids the heap allocation of a `Box<dyn Iterator>` per line — at ~5M
/// lines per GWAS file the extra allocation and virtual dispatch are
/// measurable.
enum LineFields<'a> {
    Tab(std::str::Split<'a, char>),
    Whitespace(std::str::SplitWhitespace<'a>),
}

impl<'a> Iterator for LineFields<'a> {
    type Item = &'a str;

    #[inline]
    fn next(&mut self) -> Option<&'a str> {
        match self {
            LineFields::Tab(it) => it.next().map(|s| s.trim()),
            LineFields::Whitespace(it) => it.next(),
        }
    }
}

#[inline]
fn line_fields(line: &str) -> LineFields<'_> {
    if line.contains('\t') {
        LineFields::Tab(line.split('\t'))
    } else {
        LineFields::Whitespace(line.split_whitespace())
    }
}

/// Read reference panel file (like R's `fread(ref)`).
///
/// Streams the file line-by-line instead of materializing a
/// `Vec<GwasRecord>` first — important because a full 1000G reference is
/// ~10M rows and allocating one struct-with-Strings per row is a sizable
/// chunk of the total sumstats runtime.
///
/// Returns a `(HashMap<SNP, info>, Vec<SNP>)` pair: the map for lookups and
/// the vector for reference-file ordering. Downstream code iterates the
/// vector so the merged output is deterministic and matches R's behavior.
fn read_reference_file(
    path: &Path,
    maf_filter: f64,
) -> Result<(HashMap<String, RefSnpInfo>, Vec<String>)> {
    let mut reader = open_file_reader(path)
        .with_context(|| format!("failed to open reference file {}", path.display()))?;

    // Parse header.
    let mut header_line = String::new();
    reader
        .read_line(&mut header_line)
        .with_context(|| format!("failed to read reference header {}", path.display()))?;
    let headers: Vec<String> = if header_line.contains('\t') {
        header_line.split('\t').map(|s| s.trim().to_string()).collect()
    } else {
        header_line.split_whitespace().map(|s| s.to_string()).collect()
    };

    let detected = column_detect::detect_columns(&headers);
    let snp_idx = detected.get("SNP").context("SNP column not found in reference")?;
    let a1_idx = detected.get("A1");
    let a2_idx = detected.get("A2");
    let maf_idx = detected.get("MAF");

    // Reasonable starting capacity: the 1000G reference is ~10M SNPs, but
    // smaller reference files exist too. The allocator will grow as needed.
    let mut map: HashMap<String, RefSnpInfo> = HashMap::with_capacity(10_000_000);
    let mut order: Vec<String> = Vec::with_capacity(10_000_000);

    let mut line = String::new();
    loop {
        line.clear();
        let n = reader
            .read_line(&mut line)
            .with_context(|| format!("failed to read from {}", path.display()))?;
        if n == 0 {
            break;
        }
        // Strip trailing newline(s) without reallocating.
        let trimmed = line.trim_end_matches(['\n', '\r']);
        if trimmed.is_empty() {
            continue;
        }

        let mut snp: Option<&str> = None;
        let mut a1: &str = "";
        let mut a2: &str = "";
        let mut maf_val: Option<f64> = None;
        for (i, field) in line_fields(trimmed).enumerate() {
            if i == snp_idx {
                snp = Some(field);
            } else if a1_idx == Some(i) {
                a1 = field;
            } else if a2_idx == Some(i) {
                a2 = field;
            } else if maf_idx == Some(i) {
                maf_val = field.parse().ok();
            }
        }

        let Some(snp) = snp else { continue };

        // Filter by MAF if available (R: ref <- ref[ref$MAF >= maf.filter, ]).
        if let Some(maf) = maf_val
            && maf < maf_filter
        {
            continue;
        }

        // Skip duplicate rsIDs so the ordering vector and the map stay in sync.
        if map.contains_key(snp) {
            continue;
        }
        let snp_owned = snp.to_string();
        order.push(snp_owned.clone());
        map.insert(
            snp_owned,
            RefSnpInfo {
                a1: a1.to_ascii_uppercase(),
                a2: a2.to_ascii_uppercase(),
                maf: maf_val,
            },
        );
    }

    Ok((map, order))
}

/// A QC-passed record from a single GWAS file.
struct QcRecord {
    a1: String,
    a2: String,
    beta: f64,
    se: f64,
    maf: Option<f64>,
}

/// Indices of the columns the sumstats QC pipeline cares about in a GWAS
/// file. Computed once per file from the header, then used to drive the
/// per-line fused parse+QC loop.
struct GwasColumnPlan {
    snp: usize,
    a1: Option<usize>,
    a2: Option<usize>,
    effect: usize,
    se: Option<usize>,
    info: Option<usize>,
    maf: Option<usize>,
}

/// Read a GWAS file and apply QC filters in a single streaming pass.
///
/// The previous implementation went through `read_gwas_file_with_overrides`,
/// which materializes every row as a fully-populated `GwasRecord` (a struct
/// with three owned `String`s plus seven `Option<f64>` slots) before we
/// throw most of them away here. At ~5M rows per file that was a dominant
/// cost. This version parses the header once, then walks each line exactly
/// once, extracting only the columns sumstats needs and pushing straight
/// into the QC HashMap.
fn read_and_qc_gwas(
    path: &Path,
    ref_snps: &HashMap<String, RefSnpInfo>,
    config: &SumstatsConfig,
    trait_idx: usize,
) -> Result<HashMap<String, QcRecord>> {
    let mut reader =
        open_file_reader(path).with_context(|| format!("failed to open {}", path.display()))?;

    // Parse header.
    let mut header_line = String::new();
    reader
        .read_line(&mut header_line)
        .with_context(|| format!("failed to read header from {}", path.display()))?;
    let headers: Vec<String> = if header_line.contains('\t') {
        header_line.split('\t').map(|s| s.trim().to_string()).collect()
    } else {
        header_line.split_whitespace().map(|s| s.to_string()).collect()
    };

    // Build column overrides if a beta column name is specified for this trait.
    let col_overrides = config
        .beta_overrides
        .get(trait_idx)
        .and_then(|v| v.as_ref())
        .map(|beta_name| {
            let mut m = HashMap::new();
            m.insert("effect".to_string(), beta_name.clone());
            m
        });

    let detected = column_detect::detect_columns_with_overrides(&headers, col_overrides.as_ref());

    let plan = GwasColumnPlan {
        snp: detected
            .get("SNP")
            .with_context(|| format!("SNP column not found in {}", path.display()))?,
        a1: detected.get("A1"),
        a2: detected.get("A2"),
        effect: detected
            .get("effect")
            .with_context(|| format!("effect column not found in {}", path.display()))?,
        se: detected.get("SE"),
        info: detected.get("INFO"),
        maf: detected.get("MAF"),
    };

    let skip_or_detect = config.se_logit.get(trait_idx).copied().unwrap_or(false)
        || config.ols.get(trait_idx).copied().unwrap_or(false);

    // Preallocate assuming most reference SNPs will also be present in the
    // file — close enough that we avoid most resizes.
    let mut records: HashMap<String, QcRecord> = HashMap::with_capacity(ref_snps.len());

    let mut line = String::new();
    loop {
        line.clear();
        let n = reader
            .read_line(&mut line)
            .with_context(|| format!("failed to read from {}", path.display()))?;
        if n == 0 {
            break;
        }
        let trimmed = line.trim_end_matches(['\n', '\r']);
        if trimmed.is_empty() {
            continue;
        }

        let mut snp: Option<&str> = None;
        let mut a1_str: &str = "";
        let mut a2_str: &str = "";
        let mut effect: Option<f64> = None;
        let mut se: Option<f64> = None;
        let mut info: Option<f64> = None;
        let mut maf_raw: Option<f64> = None;

        for (i, field) in line_fields(trimmed).enumerate() {
            if i == plan.snp {
                snp = Some(field);
            } else if i == plan.effect {
                effect = field.parse().ok();
            } else if plan.a1 == Some(i) {
                a1_str = field;
            } else if plan.a2 == Some(i) {
                a2_str = field;
            } else if plan.se == Some(i) {
                se = field.parse().ok();
            } else if plan.info == Some(i) {
                info = field.parse().ok();
            } else if plan.maf == Some(i) {
                maf_raw = field.parse().ok();
            }
        }

        let Some(snp) = snp else { continue };

        // Must be in reference. Checking early lets us skip most of the
        // per-row parse work for rows that won't survive the join anyway.
        if !ref_snps.contains_key(snp) {
            continue;
        }

        // Must have effect and SE. R's sumstats stops when these are missing.
        let Some(effect) = effect else { continue };
        let Some(se) = se else { continue };
        if !effect.is_finite() || se <= 0.0 || !se.is_finite() {
            continue;
        }

        // A1 required, A2 optional (munged .sumstats.gz may lack A2).
        if a1_str.is_empty() {
            continue;
        }
        let a1 = a1_str.to_ascii_uppercase();
        let a2 = a2_str.to_ascii_uppercase();

        // Allele validation (skip if A2 missing — trust the munged file).
        if !a2.is_empty() && !valid_allele_pair(&a1, &a2, config.keep_indel) {
            continue;
        }

        // Ambiguous strand SNP filtering (skip check if A2 missing).
        if !a2.is_empty() && !config.keep_ambig && is_ambiguous_snp(&a1, &a2) {
            continue;
        }

        // INFO filter.
        if info.is_some_and(|v| v < config.info_filter) {
            continue;
        }

        // MAF filter.
        let maf = maf_raw.map(|m| if m > 0.5 { 1.0 - m } else { m });
        if maf.is_some_and(|maf_val| maf_val < config.maf_filter) {
            continue;
        }

        // Detect and convert OR to log(OR). Skip auto-detect if se_logit or
        // ols is set for this trait (effects are already on the correct scale).
        let beta = if !skip_or_detect && is_or_value(effect) {
            if effect <= 0.0 {
                continue;
            }
            effect.ln()
        } else {
            effect
        };

        records.insert(
            snp.to_string(),
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

/// Find SNPs present in ALL trait HashMaps, preserving the reference-file
/// order supplied in `ref_order`.
fn find_common_snps(
    trait_records: &[HashMap<String, QcRecord>],
    ref_order: &[String],
) -> Vec<String> {
    if trait_records.is_empty() {
        return Vec::new();
    }

    ref_order
        .iter()
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

        let ref_order: Vec<String> =
            ["rs1", "rs2", "rs3", "rs4"].iter().map(|s| s.to_string()).collect();

        let common = find_common_snps(&[m1, m2], &ref_order);
        assert_eq!(common, vec!["rs1".to_string(), "rs3".to_string()]);
    }
}
