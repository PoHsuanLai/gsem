use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::GzDecoder;

use super::column_detect::{self, DetectedColumns};

/// A single row of raw GWAS summary statistics.
#[derive(Debug, Clone)]
pub struct GwasRecord {
    /// SNP identifier (e.g., rsID)
    pub snp: String,
    /// Effect allele
    pub a1: Option<String>,
    /// Non-effect allele
    pub a2: Option<String>,
    /// Effect size (beta or OR)
    pub effect: Option<f64>,
    /// Standard error of the effect
    pub se: Option<f64>,
    /// P-value
    pub p: Option<f64>,
    /// Sample size
    pub n: Option<f64>,
    /// Z-statistic
    pub z: Option<f64>,
    /// Imputation quality score
    pub info: Option<f64>,
    /// Minor allele frequency
    pub maf: Option<f64>,
}

/// Parsed GWAS summary statistics file.
#[derive(Debug, Clone)]
pub struct GwasData {
    /// All parsed SNP records
    pub records: Vec<GwasRecord>,
    /// Auto-detected column name mappings
    pub detected: DetectedColumns,
}

/// A munged summary statistics record (SNP, N, Z, A1, A2).
#[derive(Debug, Clone)]
pub struct MungedRecord {
    /// SNP identifier
    pub snp: String,
    /// Sample size
    pub n: f64,
    /// Z-statistic
    pub z: f64,
    /// Effect allele
    pub a1: String,
    /// Non-effect allele
    pub a2: String,
}

/// Read a raw GWAS summary statistics file.
///
/// Supports tab/space-delimited, optionally gzipped (.gz extension).
/// Auto-detects columns using the alias map.
pub fn read_gwas_file(path: &Path) -> Result<GwasData> {
    read_gwas_file_with_overrides(path, None)
}

/// Read a raw GWAS summary statistics file with optional column name overrides.
///
/// Overrides map canonical names (e.g., "SNP", "P", "effect") to actual header names.
/// User overrides take priority over alias matching.
pub fn read_gwas_file_with_overrides(
    path: &Path,
    column_overrides: Option<&HashMap<String, String>>,
) -> Result<GwasData> {
    let reader = open_file_reader(path)?;
    let mut lines = reader.lines();

    // Read header line
    let header_line = lines.next().context("empty file")??;
    let headers: Vec<String> = split_delimited(&header_line)
        .into_iter()
        .map(|s| s.to_owned())
        .collect();
    let detected = column_detect::detect_columns_with_overrides(&headers, column_overrides);

    let snp_idx = detected.get("SNP").context("SNP column not found")?;

    // Determine N column source and whether doubling is needed.
    // NEFFDIV2 / NEFF_HALF columns store N_eff / 2 and must be doubled.
    let (n_col_idx, n_multiply) = if let Some(idx) = detected.get("N") {
        (Some(idx), 1.0)
    } else if let Some(idx) = detected.get("NEFFDIV2") {
        log::info!("Detected NEFFDIV2 column — N values will be doubled");
        (Some(idx), 2.0)
    } else if let Some(idx) = detected.get("NEFF_HALF") {
        log::info!("Detected NEFF_HALF column — N values will be doubled");
        (Some(idx), 2.0)
    } else {
        (None, 1.0)
    };

    let mut records = Vec::new();

    for line_result in lines {
        let line = line_result?;
        if line.is_empty() {
            continue;
        }
        let fields = split_delimited(&line);
        if fields.len() <= snp_idx {
            continue;
        }

        let n_raw: Option<f64> = n_col_idx
            .and_then(|i| fields.get(i))
            .and_then(|s| s.parse().ok());

        let record = GwasRecord {
            snp: fields[snp_idx].to_owned(),
            a1: detected
                .get("A1")
                .and_then(|i| fields.get(i))
                .map(|s| s.to_string()),
            a2: detected
                .get("A2")
                .and_then(|i| fields.get(i))
                .map(|s| s.to_string()),
            effect: detected
                .get("effect")
                .and_then(|i| fields.get(i))
                .and_then(|s| s.parse().ok()),
            se: detected
                .get("SE")
                .and_then(|i| fields.get(i))
                .and_then(|s| s.parse().ok()),
            p: detected
                .get("P")
                .and_then(|i| fields.get(i))
                .and_then(|s| s.parse().ok()),
            n: n_raw.map(|v| v * n_multiply),
            z: detected
                .get("Z")
                .and_then(|i| fields.get(i))
                .and_then(|s| s.parse().ok()),
            info: detected
                .get("INFO")
                .and_then(|i| fields.get(i))
                .and_then(|s| s.parse().ok()),
            maf: detected
                .get("MAF")
                .and_then(|i| fields.get(i))
                .and_then(|s| s.parse().ok()),
        };

        records.push(record);
    }

    Ok(GwasData { records, detected })
}

/// Read a munged .sumstats.gz file (columns: SNP, N, Z, A1, A2).
pub fn read_sumstats(path: &Path) -> Result<Vec<MungedRecord>> {
    let reader = open_file_reader(path)?;
    let mut lines = reader.lines();

    // Read and validate header
    let header_line = lines.next().context("empty sumstats file")??;
    let headers = split_delimited(&header_line);
    let upper: Vec<String> = headers.iter().map(|h| h.to_uppercase()).collect();

    let snp_idx = upper
        .iter()
        .position(|h| h == "SNP")
        .context("SNP column not found in sumstats")?;
    let n_idx = upper
        .iter()
        .position(|h| h == "N")
        .context("N column not found in sumstats")?;
    let z_idx = upper
        .iter()
        .position(|h| h == "Z")
        .context("Z column not found in sumstats")?;
    let a1_idx = upper.iter().position(|h| h == "A1");
    let a2_idx = upper.iter().position(|h| h == "A2");

    let mut required = vec![snp_idx, n_idx, z_idx];
    if let Some(i) = a1_idx {
        required.push(i);
    }
    if let Some(i) = a2_idx {
        required.push(i);
    }
    let max_idx = *required.iter().max().unwrap();

    let mut records = Vec::new();
    let mut dropped_na = 0usize;
    // Line 1 is the header; data rows start at line 2.
    for (line_no, line_result) in lines.enumerate().map(|(i, r)| (i + 2, r)) {
        let line = line_result
            .with_context(|| format!("{}:{}: I/O error reading line", path.display(), line_no))?;
        if line.is_empty() {
            continue;
        }
        let fields = split_delimited(&line);
        if fields.len() <= max_idx {
            continue;
        }
        // Match GenomicSEM::ldsc semantics: silently drop rows with missing
        // N or Z (its reader calls na.omit() on read_delim output). Literal
        // "NA"/"NaN"/empty strings are treated as missing; anything else
        // that fails to parse is a real error.
        let n_field = fields[n_idx];
        let z_field = fields[z_idx];
        if is_na_token(n_field) || is_na_token(z_field) {
            dropped_na += 1;
            continue;
        }
        let n: f64 = n_field.parse().with_context(|| {
            format!(
                "{}:{}: invalid N value {:?}",
                path.display(),
                line_no,
                n_field
            )
        })?;
        let z: f64 = z_field.parse().with_context(|| {
            format!(
                "{}:{}: invalid Z value {:?}",
                path.display(),
                line_no,
                z_field
            )
        })?;
        records.push(MungedRecord {
            snp: fields[snp_idx].to_owned(),
            n,
            z,
            a1: a1_idx.map(|i| fields[i].to_owned()).unwrap_or_default(),
            a2: a2_idx.map(|i| fields[i].to_owned()).unwrap_or_default(),
        });
    }

    if dropped_na > 0 {
        log::info!(
            "{}: dropped {} row(s) with missing N or Z",
            path.display(),
            dropped_na
        );
    }

    Ok(records)
}

/// Returns `true` if `s` represents a missing value in a sumstats file
/// (empty, or one of `NA`, `NaN`, `.`, case-insensitive — matching
/// `readr::read_delim` / `na.omit` semantics).
fn is_na_token(s: &str) -> bool {
    let t = s.trim();
    t.is_empty()
        || t == "."
        || t.eq_ignore_ascii_case("NA")
        || t.eq_ignore_ascii_case("NaN")
        || t.eq_ignore_ascii_case("N/A")
        || t.eq_ignore_ascii_case("NULL")
}

/// Convert munged records into an LDSC TraitSumstats struct.
pub fn records_to_trait_sumstats(records: Vec<MungedRecord>) -> gsem_ldsc::TraitSumstats {
    let (snp, z, n, a1, a2) = records.into_iter().fold(
        (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new()),
        |(mut snps, mut zs, mut ns, mut a1s, mut a2s), r| {
            snps.push(r.snp);
            zs.push(r.z);
            ns.push(r.n);
            a1s.push(r.a1);
            a2s.push(r.a2);
            (snps, zs, ns, a1s, a2s)
        },
    );
    gsem_ldsc::TraitSumstats { snp, z, n, a1, a2 }
}

/// Read multiple munged sumstats files and return LDSC trait data.
pub fn load_trait_data(paths: &[impl AsRef<Path>]) -> Result<Vec<gsem_ldsc::TraitSumstats>> {
    paths
        .iter()
        .map(|p| {
            let path = p.as_ref();
            let records = read_sumstats(path)
                .with_context(|| format!("failed to read sumstats file {}", path.display()))?;
            Ok(records_to_trait_sumstats(records))
        })
        .collect()
}

/// Open a file for reading, with automatic gzip detection based on .gz extension.
pub fn open_file_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    let is_gz = path
        .extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"));

    if is_gz {
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Split a line by tab or whitespace. Prefers tab if present.
fn split_delimited(line: &str) -> Vec<&str> {
    if line.contains('\t') {
        line.split('\t').map(|s| s.trim()).collect()
    } else {
        line.split_whitespace().collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_tmp(name: &str, contents: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir();
        let path = dir.join(format!("gsemr_test_{}_{}", std::process::id(), name));
        let mut f = File::create(&path).unwrap();
        f.write_all(contents.as_bytes()).unwrap();
        path
    }

    #[test]
    fn read_sumstats_drops_rows_with_na_n_or_z() {
        // Mirrors GenomicSEM::ldsc, which applies na.omit() after read_delim.
        let contents = "\
SNP\tN\tZ\tA1\tA2
rs1\t10000\t1.5\tA\tG
rs2\tNA\t2.0\tC\tT
rs3\t9500\tNA\tG\tA
rs4\t8000\t-0.7\tT\tC
rs5\t\t0.3\tA\tG
rs6\t.\t0.1\tC\tT
";
        let path = write_tmp("na_drop.sumstats", contents);
        let recs = read_sumstats(&path).expect("should not error on NA rows");
        assert_eq!(recs.len(), 2, "only rs1 and rs4 should survive");
        assert_eq!(recs[0].snp, "rs1");
        assert_eq!(recs[1].snp, "rs4");
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn read_gwas_file_tolerates_na_tokens_in_numeric_columns() {
        // Mirrors GenomicSEM::sumstats, which calls read.table(..., na.string =
        // c(".", NA, "NA", "")). Raw GWAS files with NA/./empty in numeric
        // columns should parse without error, producing `None` slots. Filtering
        // of NA rows happens downstream in munge::qc and sumstats::read_and_qc.
        let contents = "\
SNP\tA1\tA2\tBETA\tSE\tP\tN\tMAF
rs1\tA\tG\t0.05\t0.01\t0.001\t10000\t0.3
rs2\tA\tG\tNA\t0.01\t0.001\t10000\t0.3
rs3\tC\tT\t0.05\t.\t0.001\t10000\t0.3
rs4\tC\tT\t0.05\t0.01\tNA\t10000\t0.3
rs5\tG\tA\t0.05\t0.01\t0.01\t\t0.3
rs6\tT\tA\t0.05\t0.01\t0.01\t8000\t0.2
";
        let path = write_tmp("raw_na.gwas", contents);
        let data = read_gwas_file(&path).expect("raw reader must tolerate NA tokens");
        assert_eq!(data.records.len(), 6, "no rows should be dropped by reader");

        // rs2 has NA beta
        assert!(
            data.records[1].effect.is_none(),
            "rs2 effect should be None"
        );
        // rs3 has "." SE
        assert!(data.records[2].se.is_none(), "rs3 se should be None");
        // rs4 has NA P
        assert!(data.records[3].p.is_none(), "rs4 p should be None");
        // rs5 has empty N
        assert!(data.records[4].n.is_none(), "rs5 n should be None");
        // rs1 and rs6 are fully populated
        assert!(data.records[0].effect.is_some());
        assert!(data.records[0].n.is_some());
        assert!(data.records[5].effect.is_some());
        assert!(data.records[5].n.is_some());
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn read_sumstats_reports_file_and_line_on_bad_number() {
        let contents = "\
SNP\tN\tZ\tA1\tA2
rs1\t10000\t1.5\tA\tG
rs2\tnot_a_number\t2.0\tC\tT
";
        let path = write_tmp("bad_n.sumstats", contents);
        let err = read_sumstats(&path).expect_err("should error on garbage N");
        let msg = format!("{err:#}");
        assert!(msg.contains(":3:"), "expected line number in: {msg}");
        assert!(msg.contains("invalid N"), "expected 'invalid N' in: {msg}");
        assert!(
            msg.contains("not_a_number"),
            "expected offending value in: {msg}"
        );
        std::fs::remove_file(&path).ok();
    }
}
