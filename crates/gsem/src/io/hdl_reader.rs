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
            ld_scores.push(sf[3].parse::<f64>().unwrap_or(0.0));
        }

        let m = snps.len();
        if m > 0 {
            ld_pieces.push(LdPiece {
                snps,
                a1,
                a2,
                ld_scores,
                m,
            });
        }
    }

    if ld_pieces.is_empty() {
        bail!("no valid LD pieces loaded from {}", ld_dir.display());
    }

    Ok(ld_pieces)
}
