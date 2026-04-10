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
use super::{col_idx, parse_beta_se_columns, parse_beta_se_values};

/// Per-SNP data from merged sumstats.
#[derive(Debug, Clone)]
pub struct MergedSnp {
    pub snp: String,
    pub a1: String,
    pub a2: String,
    pub maf: f64,
    /// Chromosome (if CHR column present in the merged file).
    pub chr: Option<u8>,
    /// Base-pair position (if BP column present in the merged file).
    pub bp: Option<u64>,
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
    let headers: Vec<String> = header_line
        .split('\t')
        .map(|s| s.trim().to_string())
        .collect();

    let snp_idx = col_idx(&headers, "SNP", "merged sumstats")?;
    let a1_idx = col_idx(&headers, "A1", "merged sumstats")?;
    let a2_idx = col_idx(&headers, "A2", "merged sumstats")?;
    let maf_idx = col_idx(&headers, "MAF", "merged sumstats")?;
    // CHR and BP are optional (needed only for Manhattan plots)
    let chr_idx = headers.iter().position(|h| h.eq_ignore_ascii_case("CHR"));
    let bp_idx = headers
        .iter()
        .position(|h| h.eq_ignore_ascii_case("BP") || h.eq_ignore_ascii_case("POS"));

    let bs = parse_beta_se_columns(&headers, "merged sumstats")?;

    let max_idx = *[snp_idx, a1_idx, a2_idx, maf_idx]
        .iter()
        .chain(chr_idx.iter())
        .chain(bp_idx.iter())
        .chain(bs.beta_indices.iter())
        .chain(bs.se_indices.iter())
        .max()
        .unwrap();

    let mut snps = Vec::new();

    for line_result in lines {
        let line = line_result?;
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() <= max_idx {
            continue;
        }

        let maf: f64 = match fields[maf_idx].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };

        let Some((beta, se)) = parse_beta_se_values(&fields, &bs.beta_indices, &bs.se_indices)
        else {
            continue;
        };

        let chr = chr_idx.and_then(|i| fields.get(i).and_then(|s| s.trim().parse().ok()));
        let bp = bp_idx.and_then(|i| fields.get(i).and_then(|s| s.trim().parse().ok()));

        snps.push(MergedSnp {
            snp: fields[snp_idx].to_string(),
            a1: fields[a1_idx].to_uppercase(),
            a2: fields[a2_idx].to_uppercase(),
            maf,
            chr,
            bp,
            beta,
            se,
        });
    }

    log::info!(
        "Read {} SNPs with {} traits from {}",
        snps.len(),
        bs.trait_names.len(),
        path.display()
    );

    Ok(MergedSumstats {
        snps,
        trait_names: bs.trait_names,
    })
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
