use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::Compression;
use flate2::write::GzEncoder;

use super::gwas_reader::MungedRecord;

/// Write munged summary statistics to a .sumstats.gz file.
///
/// Output format: tab-delimited with columns SNP, N, Z, A1, A2.
pub fn write_sumstats_gz(records: &[MungedRecord], path: &Path) -> Result<()> {
    let file = File::create(path).with_context(|| format!("cannot create {}", path.display()))?;
    let gz = GzEncoder::new(file, Compression::default());
    let mut writer = BufWriter::new(gz);

    writeln!(writer, "SNP\tN\tZ\tA1\tA2")?;
    for rec in records {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            rec.snp, rec.n, rec.z, rec.a1, rec.a2
        )?;
    }

    writer.flush()?;
    Ok(())
}

/// Write a TSV file from rows of string vectors.
pub fn write_tsv(headers: &[&str], rows: &[Vec<String>], path: &Path) -> Result<()> {
    let file = File::create(path).with_context(|| format!("cannot create {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "{}", headers.join("\t"))?;
    for row in rows {
        writeln!(writer, "{}", row.join("\t"))?;
    }

    writer.flush()?;
    Ok(())
}
