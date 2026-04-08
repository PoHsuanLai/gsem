use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::GzDecoder;

use super::column_detect::{DetectedColumns, detect_columns};

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
    let reader = open_file_reader(path)?;
    let mut lines = reader.lines();

    // Read header line
    let header_line = lines.next().context("empty file")??;
    let headers = split_delimited(&header_line);
    let detected = detect_columns(&headers);

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
            snp: fields[snp_idx].clone(),
            a1: detected.get("A1").and_then(|i| fields.get(i)).cloned(),
            a2: detected.get("A2").and_then(|i| fields.get(i)).cloned(),
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
    let a1_idx = upper
        .iter()
        .position(|h| h == "A1")
        .context("A1 column not found in sumstats")?;
    let a2_idx = upper
        .iter()
        .position(|h| h == "A2")
        .context("A2 column not found in sumstats")?;

    let mut records = Vec::new();
    for line_result in lines {
        let line = line_result?;
        if line.is_empty() {
            continue;
        }
        let fields = split_delimited(&line);
        let n: f64 = fields[n_idx].parse().context("invalid N")?;
        let z: f64 = fields[z_idx].parse().context("invalid Z")?;
        records.push(MungedRecord {
            snp: fields[snp_idx].clone(),
            n,
            z,
            a1: fields[a1_idx].clone(),
            a2: fields[a2_idx].clone(),
        });
    }

    Ok(records)
}

/// Open a file for reading, with automatic gzip detection based on .gz extension.
fn open_file_reader(path: &Path) -> Result<Box<dyn BufRead>> {
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
fn split_delimited(line: &str) -> Vec<String> {
    if line.contains('\t') {
        line.split('\t').map(|s| s.trim().to_string()).collect()
    } else {
        line.split_whitespace().map(|s| s.to_string()).collect()
    }
}
