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

    // Find fixed columns
    let gene_idx = col_idx(&headers, "Gene")?;
    let panel_idx = col_idx(&headers, "Panel")?;
    let hsq_idx = col_idx(&headers, "HSQ")?;

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
        anyhow::bail!("no beta.* columns found in TWAS sumstats");
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

    let mut genes = Vec::new();

    for line_result in lines {
        let line = line_result?;
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let max_idx = *[gene_idx, panel_idx, hsq_idx]
            .iter()
            .chain(beta_indices.iter())
            .chain(se_indices.iter())
            .max()
            .unwrap();
        if fields.len() <= max_idx {
            continue;
        }

        let hsq: f64 = match fields[hsq_idx].parse() {
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
        k,
        path.display()
    );

    Ok(TwasSumstats { genes, trait_names })
}

fn col_idx(headers: &[String], name: &str) -> Result<usize> {
    let upper = name.to_uppercase();
    headers
        .iter()
        .position(|h| h.to_uppercase() == upper)
        .with_context(|| format!("{name} column not found in TWAS sumstats"))
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
}
