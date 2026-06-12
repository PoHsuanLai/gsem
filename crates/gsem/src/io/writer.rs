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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::gwas_reader::read_sumstats;
    use std::io::{BufRead, BufReader};

    fn temp_path(name: &str) -> std::path::PathBuf {
        let mut p = std::env::temp_dir();
        // Vary by name + pid to avoid collisions across parallel test runs.
        p.push(format!("gsem_writer_test_{}_{}", std::process::id(), name));
        p
    }

    #[test]
    fn test_write_sumstats_gz_roundtrips() {
        let records = vec![
            MungedRecord {
                snp: "rs1".into(),
                n: 1000.0,
                z: 1.5,
                a1: "A".into(),
                a2: "G".into(),
            },
            MungedRecord {
                snp: "rs2".into(),
                n: 2000.0,
                z: -2.25,
                a1: "C".into(),
                a2: "T".into(),
            },
        ];
        let path = temp_path("roundtrip.sumstats.gz");
        write_sumstats_gz(&records, &path).unwrap();

        let read_back = read_sumstats(&path).unwrap();
        std::fs::remove_file(&path).ok();

        assert_eq!(read_back.len(), records.len());
        for (orig, got) in records.iter().zip(read_back.iter()) {
            assert_eq!(orig.snp, got.snp);
            assert_eq!(orig.a1, got.a1);
            assert_eq!(orig.a2, got.a2);
            assert!((orig.n - got.n).abs() < 1e-9);
            assert!((orig.z - got.z).abs() < 1e-9);
        }
    }

    #[test]
    fn test_write_tsv_writes_header_and_rows() {
        let headers = ["SNP", "BETA", "SE"];
        let rows = vec![
            vec!["rs1".to_string(), "0.1".to_string(), "0.05".to_string()],
            vec!["rs2".to_string(), "-0.2".to_string(), "0.06".to_string()],
        ];
        let path = temp_path("out.tsv");
        write_tsv(&headers, &rows, &path).unwrap();

        let f = File::open(&path).unwrap();
        let lines: Vec<String> = BufReader::new(f).lines().map(|l| l.unwrap()).collect();
        std::fs::remove_file(&path).ok();

        assert_eq!(lines[0], "SNP\tBETA\tSE");
        assert_eq!(lines[1], "rs1\t0.1\t0.05");
        assert_eq!(lines[2], "rs2\t-0.2\t0.06");
        assert_eq!(lines.len(), 3);
    }
}
