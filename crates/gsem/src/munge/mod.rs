pub mod allele;
pub mod qc;

use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};

use crate::io::gwas_reader::{self, MungedRecord};
use crate::io::writer;

/// Configuration for the munge pipeline.
#[derive(Debug, Clone)]
pub struct MungeConfig {
    pub info_filter: f64,
    pub maf_filter: f64,
    pub n_override: Option<f64>,
    pub column_overrides: Option<HashMap<String, String>>,
}

impl Default for MungeConfig {
    fn default() -> Self {
        Self {
            info_filter: 0.9,
            maf_filter: 0.01,
            n_override: None,
            column_overrides: None,
        }
    }
}

/// Reference panel record (e.g., HapMap3 SNP list).
#[derive(Debug, Clone)]
pub struct RefSnp {
    pub snp: String,
    pub a1: String,
    pub a2: String,
}

/// Read a HapMap3 reference SNP list.
///
/// Expects a file with at least columns: SNP, A1, A2.
pub fn read_reference(path: &Path) -> Result<HashMap<String, RefSnp>> {
    let data = gwas_reader::read_gwas_file(path)?;
    let _snp_idx = data
        .detected
        .get("SNP")
        .context("SNP column not found in reference")?;
    let _a1_idx = data
        .detected
        .get("A1")
        .context("A1 column not found in reference")?;
    let _a2_idx = data
        .detected
        .get("A2")
        .context("A2 column not found in reference")?;

    let mut map = HashMap::new();
    for rec in &data.records {
        let a1 = rec.a1.as_deref().unwrap_or("").to_uppercase();
        let a2 = rec.a2.as_deref().unwrap_or("").to_uppercase();
        map.insert(
            rec.snp.clone(),
            RefSnp {
                snp: rec.snp.clone(),
                a1,
                a2,
            },
        );
    }
    Ok(map)
}

/// Run the full munge pipeline on a single GWAS summary statistics file.
///
/// Returns the munged records (SNP, N, Z, A1, A2).
pub fn munge_file(
    gwas_path: &Path,
    reference: &HashMap<String, RefSnp>,
    config: &MungeConfig,
) -> Result<Vec<MungedRecord>> {
    let data = gwas_reader::read_gwas_file_with_overrides(gwas_path, config.column_overrides.as_ref())?;
    let n_input = data.records.len();
    log::info!("Read {n_input} SNPs from {}", gwas_path.display());

    let records = qc::run_qc_pipeline(data.records, reference, config)?;
    let n_output = records.len();
    log::info!(
        "After QC: {n_output} SNPs remain ({} removed)",
        n_input - n_output
    );

    Ok(records)
}

/// Munge a file and write the output to .sumstats.gz.
pub fn munge_and_write(
    gwas_path: &Path,
    reference: &HashMap<String, RefSnp>,
    config: &MungeConfig,
    output_path: &Path,
) -> Result<()> {
    let records = munge_file(gwas_path, reference, config)?;
    writer::write_sumstats_gz(&records, output_path)?;
    log::info!("Wrote {} SNPs to {}", records.len(), output_path.display());
    Ok(())
}
