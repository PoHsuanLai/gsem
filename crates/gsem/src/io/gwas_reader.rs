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
    pub snp: String,
    pub a1: Option<String>,
    pub a2: Option<String>,
    pub effect: Option<f64>,
    pub se: Option<f64>,
    pub p: Option<f64>,
    pub n: Option<f64>,
    pub z: Option<f64>,
    pub info: Option<f64>,
    pub maf: Option<f64>,
}

/// Parsed GWAS summary statistics file.
#[derive(Debug, Clone)]
pub struct GwasData {
    pub records: Vec<GwasRecord>,
    pub detected: DetectedColumns,
}

/// A munged summary statistics record (SNP, N, Z, A1, A2).
#[derive(Debug, Clone)]
pub struct MungedRecord {
    pub snp: String,
    pub n: f64,
    pub z: f64,
    pub a1: String,
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
    for line_result in lines {
        let line = line_result?;
        if line.is_empty() {
            continue;
        }
        let fields = split_delimited(&line);
        if fields.len() <= max_idx {
            continue;
        }
        let n: f64 = fields[n_idx].parse().context("invalid N")?;
        let z: f64 = fields[z_idx].parse().context("invalid Z")?;
        records.push(MungedRecord {
            snp: fields[snp_idx].to_owned(),
            n,
            z,
            a1: a1_idx.map(|i| fields[i].to_owned()).unwrap_or_default(),
            a2: a2_idx.map(|i| fields[i].to_owned()).unwrap_or_default(),
        });
    }

    Ok(records)
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
            let records = read_sumstats(p.as_ref())?;
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
