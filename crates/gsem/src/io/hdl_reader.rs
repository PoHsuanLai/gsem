//! Reader for HDL LD reference panel pieces.
//!
//! Parses the text-format HDL LD reference directory containing:
//! - `pieces.tsv` or `pieces.txt`: index of chromosome/piece pairs
//! - Per-piece SNP files: `chr{chr}.{piece}.snps.tsv` or `piece.{piece}.snps.txt`

use std::path::Path;

use anyhow::{Context, Result, bail};
use gsem_ldsc::hdl::LdPiece;

/// Load HDL LD pieces from a text-format reference directory.
///
/// Expects a `pieces.tsv` (or `pieces.txt`) index file listing chromosome/piece pairs,
/// and per-piece SNP files with columns: SNP, A1, A2, LD_score.
pub fn load_hdl_pieces(ld_dir: &Path) -> Result<Vec<LdPiece>> {
    let pieces_file = ld_dir.join("pieces.tsv");
    let pieces_alt = ld_dir.join("pieces.txt");
    let pieces_path = if pieces_file.exists() {
        pieces_file
    } else if pieces_alt.exists() {
        pieces_alt
    } else {
        bail!(
            "HDL LD reference directory missing pieces.tsv at {}",
            ld_dir.display()
        );
    };

    let pieces_content = std::fs::read_to_string(&pieces_path)
        .with_context(|| format!("failed to read pieces file: {}", pieces_path.display()))?;

    let mut ld_pieces = Vec::new();
    for line in pieces_content.lines() {
        let line = line.trim();
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("piece")
            || line.starts_with("chr")
        {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }

        let chr_val = fields[0];
        let piece_val = fields[1];
        let snp_file_tsv = ld_dir.join(format!("chr{chr_val}.{piece_val}.snps.tsv"));
        let snp_file_txt = ld_dir.join(format!("piece.{piece_val}.snps.txt"));
        let snp_file = if snp_file_tsv.exists() {
            snp_file_tsv
        } else if snp_file_txt.exists() {
            snp_file_txt
        } else {
            continue;
        };

        let snp_content = match std::fs::read_to_string(&snp_file) {
            Ok(c) => c,
            Err(_) => continue,
        };

        let mut snps = Vec::new();
        let mut a1 = Vec::new();
        let mut a2 = Vec::new();
        let mut ld_scores = Vec::new();

        for sline in snp_content.lines() {
            let sline = sline.trim();
            if sline.is_empty() || sline.starts_with('#') || sline.starts_with("SNP") {
                continue;
            }
            let sf: Vec<&str> = sline.split('\t').collect();
            if sf.len() < 4 {
                continue;
            }
            snps.push(sf[0].to_string());
            a1.push(sf[1].to_string());
            a2.push(sf[2].to_string());
            let ld_val = sf[3].parse::<f64>().map_err(|_| {
                anyhow::anyhow!("Non-numeric LD score '{}' in HDL reference", sf[3])
            })?;
            ld_scores.push(ld_val);
        }

        let m = snps.len();
        if m > 0 {
            // Eigen-decomposition of the block LD matrix (lam / V), required by
            // the HDL likelihood. Stored alongside the SNP file.
            let eigen_file = ld_dir.join(format!("chr{chr_val}.{piece_val}.eigen.tsv"));
            let (eigenvalues, eigenvectors) = gsem_ldsc::hdl::read_eigen_file(&eigen_file)
                .with_context(|| {
                    format!(
                        "HDL piece chr{chr_val}.{piece_val} is missing its eigen file {}",
                        eigen_file.display()
                    )
                })?;
            ld_pieces.push(LdPiece {
                snps,
                a1,
                a2,
                ld_scores,
                eigenvalues,
                eigenvectors,
                m,
            });
        }
    }

    if ld_pieces.is_empty() {
        bail!("no valid LD pieces loaded from {}", ld_dir.display());
    }

    Ok(ld_pieces)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Build a minimal one-piece HDL reference directory: a pieces index, a
    /// 2-SNP per-piece file, and its eigen file (eigenvalues row + vector rows).
    fn write_hdl_dir() -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("gsem_hdl_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();

        let mut pieces = std::fs::File::create(dir.join("pieces.tsv")).unwrap();
        writeln!(pieces, "chr\tpiece\tn").unwrap();
        writeln!(pieces, "1\t1\t2").unwrap();

        let mut snps = std::fs::File::create(dir.join("chr1.1.snps.tsv")).unwrap();
        writeln!(snps, "SNP\tA1\tA2\tLD_score").unwrap();
        writeln!(snps, "rs1\tA\tG\t1.20").unwrap();
        writeln!(snps, "rs2\tC\tT\t0.80").unwrap();

        // Eigen of a 2x2 block: row 0 = eigenvalues, then 2 eigenvector rows.
        let mut eig = std::fs::File::create(dir.join("chr1.1.eigen.tsv")).unwrap();
        writeln!(eig, "1.5\t0.5").unwrap();
        writeln!(eig, "0.7071067811865476\t-0.7071067811865476").unwrap();
        writeln!(eig, "0.7071067811865476\t0.7071067811865476").unwrap();

        dir
    }

    #[test]
    fn test_load_hdl_pieces() {
        let dir = write_hdl_dir();
        let pieces = load_hdl_pieces(&dir).unwrap();
        std::fs::remove_dir_all(&dir).ok();

        assert_eq!(pieces.len(), 1);
        let p = &pieces[0];
        assert_eq!(p.m, 2);
        assert_eq!(p.snps, vec!["rs1", "rs2"]);
        assert_eq!(p.a1, vec!["A", "C"]);
        assert_eq!(p.a2, vec!["G", "T"]);
        assert!((p.ld_scores[0] - 1.20).abs() < 1e-10);
        assert_eq!(p.eigenvalues.len(), 2);
        assert!((p.eigenvalues[0] - 1.5).abs() < 1e-10);
        assert_eq!(p.eigenvectors.nrows(), 2);
        assert_eq!(p.eigenvectors.ncols(), 2);
    }

    #[test]
    fn test_load_hdl_pieces_missing_index_errors() {
        let dir = std::env::temp_dir().join(format!("gsem_hdl_empty_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let res = load_hdl_pieces(&dir);
        std::fs::remove_dir_all(&dir).ok();
        assert!(res.is_err(), "missing pieces.tsv should error");
    }
}
