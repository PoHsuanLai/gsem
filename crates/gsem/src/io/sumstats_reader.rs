//! Reader for merged summary statistics files (output of sumstats() in R).
//!
//! Expected format: tab-delimited with columns:
//!   SNP, A1, A2, MAF, beta.Trait1, se.Trait1, beta.Trait2, se.Trait2, ...
//!
//! The beta.* and se.* columns are matched by prefix.

use std::io::BufRead;
use std::path::Path;

use anyhow::{Context, Result};

use super::gwas_reader::open_file_reader;

/// Per-SNP data from merged sumstats.
#[derive(Debug, Clone)]
pub struct MergedSnp {
    pub snp: String,
    pub a1: String,
    pub a2: String,
    pub maf: f64,
    /// beta per trait (length k)
    pub beta: Vec<f64>,
    /// SE per trait (length k)
    pub se: Vec<f64>,
}

/// Result of reading merged sumstats.
#[derive(Debug)]
pub struct MergedSumstats {
    pub snps: Vec<MergedSnp>,
    /// Trait names extracted from column headers (beta.X -> X)
    pub trait_names: Vec<String>,
}

/// Read a merged summary statistics file.
///
/// Columns are identified by name:
/// - `SNP` — SNP identifier
/// - `A1`, `A2` — alleles
/// - `MAF` — minor allele frequency
/// - `beta.*` — effect sizes per trait
/// - `se.*` — standard errors per trait
pub fn read_merged_sumstats(path: &Path) -> Result<MergedSumstats> {
    let reader = open_file_reader(path)?;
    let mut lines = reader.lines();

    let header_line = lines.next().context("empty file")??;
    let headers: Vec<String> = header_line.split('\t').map(|s| s.trim().to_string()).collect();

    // Find fixed columns
    let snp_idx = col_idx(&headers, "SNP")?;
    let a1_idx = col_idx(&headers, "A1")?;
    let a2_idx = col_idx(&headers, "A2")?;
    let maf_idx = col_idx(&headers, "MAF")?;

    // Find beta.* and se.* columns (case-insensitive prefix match)
    let mut beta_cols: Vec<(usize, String)> = Vec::new();
    let mut se_cols: Vec<(usize, String)> = Vec::new();

    for (i, h) in headers.iter().enumerate() {
        let lower = h.to_lowercase();
        if let Some(name) = lower.strip_prefix("beta.") {
            beta_cols.push((i, name.to_string()));
        } else if let Some(name) = lower.strip_prefix("se.") {
            se_cols.push((i, name.to_string()));
        }
    }

    if beta_cols.is_empty() {
        anyhow::bail!("no beta.* columns found in merged sumstats");
    }
    if beta_cols.len() != se_cols.len() {
        anyhow::bail!(
            "mismatched beta ({}) and se ({}) column counts",
            beta_cols.len(),
            se_cols.len()
        );
    }

    let trait_names: Vec<String> = beta_cols.iter().map(|(_, name)| name.clone()).collect();
    let k = trait_names.len();

    let beta_indices: Vec<usize> = beta_cols.iter().map(|(i, _)| *i).collect();
    let se_indices: Vec<usize> = se_cols.iter().map(|(i, _)| *i).collect();

    let mut snps = Vec::new();

    for line_result in lines {
        let line = line_result?;
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let max_idx = *[snp_idx, a1_idx, a2_idx, maf_idx]
            .iter()
            .chain(beta_indices.iter())
            .chain(se_indices.iter())
            .max()
            .unwrap();
        if fields.len() <= max_idx {
            continue;
        }

        let maf: f64 = match fields[maf_idx].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };

        let mut beta = Vec::with_capacity(k);
        let mut se = Vec::with_capacity(k);
        let mut skip = false;

        for t in 0..k {
            match (
                fields[beta_indices[t]].parse::<f64>(),
                fields[se_indices[t]].parse::<f64>(),
            ) {
                (Ok(b), Ok(s)) if b.is_finite() && s.is_finite() => {
                    beta.push(b);
                    se.push(s);
                }
                _ => {
                    skip = true;
                    break;
                }
            }
        }

        if skip {
            continue;
        }

        snps.push(MergedSnp {
            snp: fields[snp_idx].to_string(),
            a1: fields[a1_idx].to_uppercase(),
            a2: fields[a2_idx].to_uppercase(),
            maf,
            beta,
            se,
        });
    }

    log::info!(
        "Read {} SNPs with {} traits from {}",
        snps.len(),
        k,
        path.display()
    );

    Ok(MergedSumstats { snps, trait_names })
}

fn col_idx(headers: &[String], name: &str) -> Result<usize> {
    let upper = name.to_uppercase();
    headers
        .iter()
        .position(|h| h.to_uppercase() == upper)
        .with_context(|| format!("{name} column not found in merged sumstats"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_read_merged_sumstats() {
        let dir = std::env::temp_dir().join("gsem_test_merged");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("test_merged.tsv");

        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "SNP\tA1\tA2\tMAF\tbeta.T1\tse.T1\tbeta.T2\tse.T2").unwrap();
        writeln!(f, "rs1\tA\tG\t0.3\t0.05\t0.01\t-0.03\t0.02").unwrap();
        writeln!(f, "rs2\tC\tT\t0.1\t0.10\t0.02\t0.08\t0.03").unwrap();

        let result = read_merged_sumstats(&path).unwrap();
        assert_eq!(result.snps.len(), 2);
        assert_eq!(result.trait_names, vec!["t1", "t2"]);
        assert_eq!(result.snps[0].snp, "rs1");
        assert!((result.snps[0].beta[0] - 0.05).abs() < 1e-10);
        assert!((result.snps[1].se[1] - 0.03).abs() < 1e-10);

        std::fs::remove_dir_all(&dir).ok();
    }
}
