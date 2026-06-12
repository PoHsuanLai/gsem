//! Reader for TWAS merged summary statistics files.
//!
//! Expected format: tab-delimited with columns:
//!   Gene, Panel, HSQ, beta.Trait1, se.Trait1, beta.Trait2, se.Trait2, ...
//!
//! The beta.* and se.* columns are matched by prefix, same as sumstats_reader.

use std::io::BufRead;
use std::path::Path;

use anyhow::{Context, Result};

use super::gwas_reader::open_file_reader;
use super::{col_idx, parse_beta_se_columns, parse_beta_se_values};

/// Per-gene data from TWAS merged file.
#[derive(Debug, Clone)]
pub struct TwasGene {
    pub gene: String,
    pub panel: String,
    pub hsq: f64,
    /// beta per trait (length k)
    pub beta: Vec<f64>,
    /// SE per trait (length k)
    pub se: Vec<f64>,
}

/// Result of reading TWAS merged sumstats.
#[derive(Debug)]
pub struct TwasSumstats {
    pub genes: Vec<TwasGene>,
    /// Trait names extracted from column headers (beta.X -> X)
    pub trait_names: Vec<String>,
}

/// Read a TWAS merged summary statistics file.
///
/// Columns are identified by name:
/// - `Gene` — gene identifier
/// - `Panel` — expression panel name
/// - `HSQ` — heritability of expression
/// - `beta.*` — effect sizes per trait
/// - `se.*` — standard errors per trait
pub fn read_twas_sumstats(path: &Path) -> Result<TwasSumstats> {
    let reader = open_file_reader(path)?;
    let mut lines = reader.lines();

    let header_line = lines.next().context("empty file")??;
    let headers: Vec<String> = header_line
        .split('\t')
        .map(|s| s.trim().to_string())
        .collect();

    let gene_idx = col_idx(&headers, "Gene", "TWAS sumstats")?;
    let panel_idx = col_idx(&headers, "Panel", "TWAS sumstats")?;
    let hsq_idx = col_idx(&headers, "HSQ", "TWAS sumstats")?;

    let bs = parse_beta_se_columns(&headers, "TWAS sumstats")?;

    let max_idx = *[gene_idx, panel_idx, hsq_idx]
        .iter()
        .chain(bs.beta_indices.iter())
        .chain(bs.se_indices.iter())
        .max()
        .unwrap();

    let mut genes = Vec::new();

    for line_result in lines {
        let line = line_result?;
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() <= max_idx {
            continue;
        }

        let hsq: f64 = match fields[hsq_idx].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };

        let Some((beta, se)) = parse_beta_se_values(&fields, &bs.beta_indices, &bs.se_indices)
        else {
            continue;
        };

        genes.push(TwasGene {
            gene: fields[gene_idx].to_string(),
            panel: fields[panel_idx].to_string(),
            hsq,
            beta,
            se,
        });
    }

    log::info!(
        "Read {} genes with {} traits from {}",
        genes.len(),
        bs.trait_names.len(),
        path.display()
    );

    Ok(TwasSumstats {
        genes,
        trait_names: bs.trait_names,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_read_twas_sumstats() {
        let dir = std::env::temp_dir().join("gsem_test_twas");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("test_twas.tsv");

        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "Gene\tPanel\tHSQ\tbeta.T1\tse.T1\tbeta.T2\tse.T2").unwrap();
        writeln!(f, "BRCA1\tGTEx_Brain\t0.15\t0.05\t0.01\t-0.03\t0.02").unwrap();
        writeln!(f, "TP53\tGTEx_Liver\t0.22\t0.10\t0.02\t0.08\t0.03").unwrap();

        let result = read_twas_sumstats(&path).unwrap();
        assert_eq!(result.genes.len(), 2);
        assert_eq!(result.trait_names, vec!["t1", "t2"]);
        assert_eq!(result.genes[0].gene, "BRCA1");
        assert_eq!(result.genes[0].panel, "GTEx_Brain");
        assert!((result.genes[0].hsq - 0.15).abs() < 1e-10);
        assert!((result.genes[0].beta[0] - 0.05).abs() < 1e-10);
        assert!((result.genes[1].se[1] - 0.03).abs() < 1e-10);

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_read_twas_sumstats_skips_malformed_rows() {
        // Exercises every skip branch in the row loop: blank line, a row too
        // short to reach max_idx, an unparseable HSQ, and an unparseable
        // beta/se. Only the two well-formed genes should survive.
        let dir = std::env::temp_dir().join("gsem_test_twas_malformed");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("malformed.tsv");

        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "Gene\tPanel\tHSQ\tbeta.T1\tse.T1").unwrap();
        writeln!(f, "GOOD1\tP\t0.10\t0.05\t0.01").unwrap();
        writeln!(f).unwrap(); // blank line -> skipped
        writeln!(f, "SHORT\tP\t0.10").unwrap(); // too few fields -> skipped
        writeln!(f, "BADHSQ\tP\tNA\t0.05\t0.01").unwrap(); // HSQ unparseable -> skipped
        writeln!(f, "BADBETA\tP\t0.10\tNA\t0.01").unwrap(); // beta unparseable -> skipped
        writeln!(f, "GOOD2\tP\t0.20\t0.07\t0.02").unwrap();

        let result = read_twas_sumstats(&path).unwrap();
        let names: Vec<&str> = result.genes.iter().map(|g| g.gene.as_str()).collect();
        assert_eq!(names, vec!["GOOD1", "GOOD2"], "only well-formed rows kept");

        std::fs::remove_dir_all(&dir).ok();
    }
}
