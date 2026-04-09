use std::fs;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::GzDecoder;

/// A single LD score record.
#[derive(Debug, Clone)]
pub struct LdScoreRecord {
    pub chr: u8,
    pub snp: String,
    pub bp: u64,
    pub l2: f64,
}

/// LD scores and associated data loaded from reference panel files.
#[derive(Debug, Clone)]
pub struct LdScores {
    pub records: Vec<LdScoreRecord>,
    pub w_ld: Vec<f64>,
    pub m_per_chr: Vec<f64>,
    pub total_m: f64,
}

/// Read LD score files for chromosomes 1..n_chr.
///
/// Expects files at `{dir}/{chr}.l2.ldscore.gz` with columns: CHR, SNP, BP, L2
/// and M count files at `{dir}/{chr}.l2.M_5_50`.
/// Weight files at `{wld_dir}/{chr}.l2.ldscore.gz`.
pub fn read_ld_scores(ld_dir: &Path, wld_dir: &Path, n_chr: usize) -> Result<LdScores> {
    let mut records = Vec::new();
    let mut w_ld_all = Vec::new();
    let mut m_per_chr = Vec::with_capacity(n_chr);
    let mut total_m = 0.0;

    for chr in 1..=n_chr {
        let ld_path = ld_dir.join(format!("{chr}.l2.ldscore.gz"));
        let chr_records = read_ld_score_file(&ld_path)?;

        let wld_path = wld_dir.join(format!("{chr}.l2.ldscore.gz"));
        let wld_records = read_ld_score_file(&wld_path)?;

        let m_path = ld_dir.join(format!("{chr}.l2.M_5_50"));
        let m = read_m_file(&m_path)?;
        total_m += m;
        m_per_chr.push(m);

        w_ld_all.extend(wld_records.iter().map(|r| r.l2));
        records.extend(chr_records);
    }

    Ok(LdScores {
        records,
        w_ld: w_ld_all,
        m_per_chr,
        total_m,
    })
}

/// Read a single chromosome LD score file (.l2.ldscore.gz).
///
/// Detects L2 column by header name. Common formats:
/// - `CHR SNP BP L2` (4 columns)
/// - `CHR SNP BP CM MAF L2` (6 columns)
fn read_ld_score_file(path: &Path) -> Result<Vec<LdScoreRecord>> {
    let file =
        std::fs::File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    let reader = BufReader::new(GzDecoder::new(file));
    let mut records = Vec::new();

    let mut lines = reader.lines();
    let header = lines.next().context("empty LD score file")??;
    let header_fields: Vec<&str> = if header.contains('\t') {
        header.split('\t').collect()
    } else {
        header.split_whitespace().collect()
    };

    // Find column indices by name
    let chr_idx = header_fields
        .iter()
        .position(|&h| h == "CHR")
        .context("CHR column not found in LD score file")?;
    let snp_idx = header_fields
        .iter()
        .position(|&h| h == "SNP")
        .context("SNP column not found in LD score file")?;
    let bp_idx = header_fields
        .iter()
        .position(|&h| h == "BP")
        .context("BP column not found in LD score file")?;
    let l2_idx = header_fields
        .iter()
        .position(|&h| h == "L2")
        .context("L2 column not found in LD score file")?;

    let min_cols = [chr_idx, snp_idx, bp_idx, l2_idx]
        .into_iter()
        .max()
        .unwrap()
        + 1;

    for line_result in lines {
        let line = line_result?;
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = if line.contains('\t') {
            line.split('\t').collect()
        } else {
            line.split_whitespace().collect()
        };

        if fields.len() < min_cols {
            continue;
        }

        let chr: u8 = fields[chr_idx].parse().context("invalid CHR")?;
        let snp = fields[snp_idx].to_string();
        let bp: u64 = fields[bp_idx].parse().context("invalid BP")?;
        let l2: f64 = fields[l2_idx].parse().context("invalid L2")?;

        records.push(LdScoreRecord { chr, snp, bp, l2 });
    }

    Ok(records)
}

/// Read M count file (.l2.M_5_50). Contains a single line with one or more
/// whitespace-separated numbers; we sum them.
fn read_m_file(path: &Path) -> Result<f64> {
    let content =
        fs::read_to_string(path).with_context(|| format!("cannot open {}", path.display()))?;
    content
        .split_whitespace()
        .map(|s| {
            s.parse::<f64>()
                .with_context(|| format!("invalid M value: {s}"))
        })
        .try_fold(0.0, |acc, r| r.map(|v| acc + v))
}
