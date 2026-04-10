//! Reader for merged summary statistics files (output of sumstats() in R).
//!
//! Expected format: tab-delimited with columns:
//!   SNP, A1, A2, MAF, beta.Trait1, se.Trait1, beta.Trait2, se.Trait2, ...
//!
//! The beta.* and se.* columns are matched by prefix.
//!
//! # Memory layout
//!
//! We deliberately store SNPs as a Structure-of-Arrays rather than
//! `Vec<MergedSnp>`. At ~5M SNPs × 3 traits, the AoS layout cost
//! ~264 bytes/SNP (verified by `size_of` + allocator rounding).
//! Switching to parallel `Vec<String>` / `Vec<u8>` / `Vec<f64>` and
//! a single flat `beta_flat: Vec<f64>` of length `n*k` drops that to
//! ~98 bytes/SNP — measured 1.31 GB → ~0.48 GB on the ANX/OCD/PTSD
//! bench.

use std::io::BufRead;
use std::path::Path;

use anyhow::{Context, Result};

use super::gwas_reader::open_file_reader;
use super::{col_idx, parse_beta_se_columns, parse_beta_se_values};

/// Merged sumstats stored column-wise for cheap load + cheap SNP-major
/// iteration.
#[derive(Debug)]
pub struct MergedSumstats {
    /// Trait names parsed from the `beta.*` column headers.
    pub trait_names: Vec<String>,
    /// Number of traits. `beta_flat.len() == snp.len() * k` always.
    pub k: usize,

    /// SNP rsIDs, one per SNP.
    pub snp: Vec<String>,
    /// Effect allele, one ASCII byte per SNP (A/C/G/T).
    pub a1: Vec<u8>,
    /// Non-effect allele, one ASCII byte per SNP.
    pub a2: Vec<u8>,
    /// Minor allele frequency per SNP.
    pub maf: Vec<f64>,

    /// Row-major beta matrix, length `n_snps * k`.
    pub beta_flat: Vec<f64>,
    /// Row-major SE matrix, length `n_snps * k`.
    pub se_flat: Vec<f64>,

    /// Optional CHR column. Only populated when the merged file had a
    /// `CHR` header. Indexed the same way as `snp`.
    pub chr: Option<Vec<u8>>,
    /// Optional BP column. Only populated when the merged file had a
    /// `BP`/`POS` header.
    pub bp: Option<Vec<u64>>,
}

impl MergedSumstats {
    /// Number of SNPs.
    #[inline]
    pub fn len(&self) -> usize {
        self.snp.len()
    }

    /// `true` when there are no SNPs.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.snp.is_empty()
    }

    /// Number of traits.
    #[inline]
    pub fn k(&self) -> usize {
        self.k
    }

    /// Beta vector for SNP `i`: `&beta_flat[i*k..(i+1)*k]`.
    #[inline]
    pub fn beta_row(&self, i: usize) -> &[f64] {
        let k = self.k;
        &self.beta_flat[i * k..(i + 1) * k]
    }

    /// SE vector for SNP `i`.
    #[inline]
    pub fn se_row(&self, i: usize) -> &[f64] {
        let k = self.k;
        &self.se_flat[i * k..(i + 1) * k]
    }

    /// Build `Vec<&[f64]>` views over the **entire** merged file.
    ///
    /// Identity-path consumers only — if you're filtering, build the
    /// ref vec directly from the filtered index list so you don't
    /// allocate `n_snps` pointers and then throw most away.
    pub fn beta_rows(&self) -> Vec<&[f64]> {
        (0..self.len()).map(|i| self.beta_row(i)).collect()
    }

    /// `Vec<&[f64]>` views for every SE row — same caveat as `beta_rows`.
    pub fn se_rows(&self) -> Vec<&[f64]> {
        (0..self.len()).map(|i| self.se_row(i)).collect()
    }

    /// `2 * maf * (1 - maf)` per SNP — SNP variance used throughout
    /// the GWAS paths.
    pub fn var_snp(&self) -> Vec<f64> {
        self.maf.iter().map(|&m| 2.0 * m * (1.0 - m)).collect()
    }

    /// A1 as a one-character owned `String`. Bytes are range-checked
    /// to `A|C|G|T` at read time, so `self.a1[i] as char` is always a
    /// valid single-char ASCII string — no `unsafe`, no UTF-8
    /// validation.
    #[inline]
    pub fn a1_string(&self, i: usize) -> String {
        String::from(self.a1[i] as char)
    }

    /// A2 as a one-character owned `String`.
    #[inline]
    pub fn a2_string(&self, i: usize) -> String {
        String::from(self.a2[i] as char)
    }
}

/// Read a merged summary statistics file.
///
/// Expected columns (case-insensitive header):
/// - `SNP` — SNP identifier
/// - `A1`, `A2` — single-nucleotide alleles
/// - `MAF` — minor allele frequency
/// - `beta.*` — effect sizes per trait
/// - `se.*`   — standard errors per trait
/// - optional `CHR` + `BP`/`POS` — for Manhattan plots
///
/// Rows with multi-character alleles or non-ACGT bytes are dropped
/// (see the module-level docs). This is consistent with
/// `SumstatsConfig::keep_indel = false` in the upstream merge
/// pipeline, which is the default.
pub fn read_merged_sumstats(path: &Path) -> Result<MergedSumstats> {
    let mut reader = open_file_reader(path)
        .with_context(|| format!("failed to open merged sumstats {}", path.display()))?;

    // Parse header.
    let mut header_line = String::new();
    reader
        .read_line(&mut header_line)
        .with_context(|| format!("failed to read header from {}", path.display()))?;
    let header_trimmed = header_line.trim_end_matches(['\n', '\r']);
    let headers: Vec<String> = header_trimmed
        .split('\t')
        .map(|s| s.trim().to_string())
        .collect();

    let snp_idx = col_idx(&headers, "SNP", "merged sumstats")?;
    let a1_idx = col_idx(&headers, "A1", "merged sumstats")?;
    let a2_idx = col_idx(&headers, "A2", "merged sumstats")?;
    let maf_idx = col_idx(&headers, "MAF", "merged sumstats")?;
    // CHR and BP are optional (needed only for Manhattan plots).
    let chr_idx = headers.iter().position(|h| h.eq_ignore_ascii_case("CHR"));
    let bp_idx = headers
        .iter()
        .position(|h| h.eq_ignore_ascii_case("BP") || h.eq_ignore_ascii_case("POS"));

    let bs = parse_beta_se_columns(&headers, "merged sumstats")?;
    let k = bs.trait_names.len();

    // We need every field we plan to read to exist on every row; use
    // the maximum column index as a quick row-length sanity check.
    let max_idx = *[snp_idx, a1_idx, a2_idx, maf_idx]
        .iter()
        .chain(chr_idx.iter())
        .chain(bp_idx.iter())
        .chain(bs.beta_indices.iter())
        .chain(bs.se_indices.iter())
        .max()
        .expect("at least SNP index is present");

    // Reasonable starting capacity — smaller than the 10M reference
    // reader because merged files are typically narrower. The
    // allocator grows as needed.
    let cap = 5_000_000;
    let mut snp: Vec<String> = Vec::with_capacity(cap);
    let mut a1: Vec<u8> = Vec::with_capacity(cap);
    let mut a2: Vec<u8> = Vec::with_capacity(cap);
    let mut maf: Vec<f64> = Vec::with_capacity(cap);
    let mut beta_flat: Vec<f64> = Vec::with_capacity(cap * k);
    let mut se_flat: Vec<f64> = Vec::with_capacity(cap * k);
    let mut chr_vec: Option<Vec<u8>> = chr_idx.map(|_| Vec::with_capacity(cap));
    let mut bp_vec: Option<Vec<u64>> = bp_idx.map(|_| Vec::with_capacity(cap));

    let mut dropped_indel: usize = 0;

    let mut line = String::new();
    loop {
        line.clear();
        let n = reader
            .read_line(&mut line)
            .with_context(|| format!("failed to read from {}", path.display()))?;
        if n == 0 {
            break;
        }
        let trimmed = line.trim_end_matches(['\n', '\r']);
        if trimmed.is_empty() {
            continue;
        }

        // Merged files are always tab-separated (we produce them
        // that way in `write_merged_sumstats`), so a single split by
        // '\t' suffices. This is cheaper than the whitespace
        // fallback the upstream reader had to keep because the 1000G
        // reference is space-separated.
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() <= max_idx {
            continue;
        }

        // A1/A2 must be single-byte A/C/G/T. Drop indels rather than
        // silently truncating them.
        let a1_bytes = fields[a1_idx].as_bytes();
        let a2_bytes = fields[a2_idx].as_bytes();
        if a1_bytes.len() != 1 || a2_bytes.len() != 1 {
            dropped_indel += 1;
            continue;
        }
        let a1_byte = a1_bytes[0].to_ascii_uppercase();
        let a2_byte = a2_bytes[0].to_ascii_uppercase();
        if !is_acgt(a1_byte) || !is_acgt(a2_byte) {
            dropped_indel += 1;
            continue;
        }

        let maf_val: f64 = match fields[maf_idx].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };

        let Some((beta, se)) = parse_beta_se_values(&fields, &bs.beta_indices, &bs.se_indices)
        else {
            continue;
        };

        // All fields valid — push in parallel.
        snp.push(fields[snp_idx].to_string());
        a1.push(a1_byte);
        a2.push(a2_byte);
        maf.push(maf_val);
        beta_flat.extend_from_slice(&beta);
        se_flat.extend_from_slice(&se);

        if let (Some(i), Some(vec)) = (chr_idx, chr_vec.as_mut()) {
            vec.push(fields.get(i).and_then(|s| s.trim().parse().ok()).unwrap_or(0));
        }
        if let (Some(i), Some(vec)) = (bp_idx, bp_vec.as_mut()) {
            vec.push(fields.get(i).and_then(|s| s.trim().parse().ok()).unwrap_or(0));
        }
    }

    if dropped_indel > 0 {
        log::warn!(
            "read_merged_sumstats: dropped {} row(s) with multi-character or non-ACGT alleles from {}",
            dropped_indel,
            path.display()
        );
    }
    log::info!(
        "Read {} SNPs with {} traits from {}",
        snp.len(),
        k,
        path.display()
    );

    Ok(MergedSumstats {
        trait_names: bs.trait_names,
        k,
        snp,
        a1,
        a2,
        maf,
        beta_flat,
        se_flat,
        chr: chr_vec,
        bp: bp_vec,
    })
}

#[inline]
fn is_acgt(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
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
        assert_eq!(result.len(), 2);
        assert_eq!(result.k(), 2);
        assert_eq!(result.trait_names, vec!["t1", "t2"]);
        assert_eq!(result.snp[0], "rs1");
        assert_eq!(result.a1[0], b'A');
        assert_eq!(result.a2[0], b'G');
        assert_eq!(result.a1_string(0), "A");
        assert_eq!(result.a2_string(0), "G");
        assert!((result.beta_row(0)[0] - 0.05).abs() < 1e-10);
        assert!((result.se_row(1)[1] - 0.03).abs() < 1e-10);
        assert!(result.chr.is_none());
        assert!(result.bp.is_none());

        // var_snp = 2 * 0.3 * 0.7 = 0.42
        let v = result.var_snp();
        assert!((v[0] - 0.42).abs() < 1e-10);

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_merged_sumstats_rejects_indels() {
        let dir = std::env::temp_dir().join("gsem_test_merged_indel");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("test_indel.tsv");

        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "SNP\tA1\tA2\tMAF\tbeta.T1\tse.T1").unwrap();
        writeln!(f, "rs1\tA\tG\t0.3\t0.05\t0.01").unwrap();
        writeln!(f, "rs2\tAC\tG\t0.1\t0.10\t0.02").unwrap(); // indel, dropped
        writeln!(f, "rs3\tC\tT\t0.2\t0.07\t0.015").unwrap();
        writeln!(f, "rs4\tN\tG\t0.15\t0.05\t0.01").unwrap(); // non-ACGT, dropped

        let result = read_merged_sumstats(&path).unwrap();
        assert_eq!(result.len(), 2);
        assert_eq!(result.snp, vec!["rs1".to_string(), "rs3".to_string()]);

        std::fs::remove_dir_all(&dir).ok();
    }
}
