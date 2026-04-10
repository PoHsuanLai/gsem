//! Reader for annotation-specific LD score files used by stratified LDSC.
//!
//! Expects per-chromosome files:
//! - `{dir}/{chr}.l2.ldscore.gz` with columns: CHR, SNP, BP, [annotation columns...]
//! - `{dir}/{chr}.l2.M_5_50` with one value per annotation (whitespace-separated)

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use faer::Mat;
use flate2::read::GzDecoder;

/// Annotation LD score data loaded from disk.
pub struct AnnotLdScores {
    /// SNP identifiers in order.
    pub snps: Vec<String>,
    /// Chromosome per SNP.
    pub chr: Vec<u32>,
    /// Base-pair position per SNP.
    pub bp: Vec<u64>,
    /// Annotation-specific LD scores (n_snps x n_annot).
    pub annot_ld: Mat<f64>,
    /// Weight LD scores per SNP.
    pub w_ld: Vec<f64>,
    /// Names of each annotation.
    pub annotation_names: Vec<String>,
    /// Number of SNPs per annotation (summed across chromosomes).
    pub m_annot: Vec<f64>,
}

/// Read annotation-specific LD score files for the given chromosomes.
///
/// Expects files at `{ld_dir}/{chr}.l2.ldscore.gz` with annotation LD columns,
/// `{ld_dir}/{chr}.l2.M_5_50` with per-annotation M counts, and weight files
/// at `{wld_dir}/{chr}.l2.ldscore.gz`.
pub fn read_annot_ld_scores(
    ld_dir: &Path,
    wld_dir: &Path,
    chromosomes: &[usize],
) -> Result<AnnotLdScores> {
    let mut all_snps = Vec::new();
    let mut all_chr = Vec::new();
    let mut all_bp = Vec::new();
    let mut all_annot_rows: Vec<Vec<f64>> = Vec::new();
    let mut all_w_ld = Vec::new();
    let mut annotation_names: Option<Vec<String>> = None;
    let mut m_annot_total: Option<Vec<f64>> = None;

    for &chr in chromosomes {
        // Read annotation LD scores
        let ld_path = ld_dir.join(format!("{chr}.l2.ldscore.gz"));
        let (snps, chr_vals, bp_vals, annot_data, names) = read_annot_ld_file(&ld_path)?;

        if let Some(ref existing) = annotation_names {
            if *existing != names {
                anyhow::bail!("annotation names differ between chromosomes");
            }
        } else {
            annotation_names = Some(names);
        }

        // Read weight LD scores (just need the L2 column)
        let wld_path = wld_dir.join(format!("{chr}.l2.ldscore.gz"));
        let w_records = read_weight_ld_file(&wld_path)?;

        // Build SNP->weight lookup for this chromosome
        let w_map: HashMap<&str, f64> = w_records
            .iter()
            .map(|(snp, l2)| (snp.as_str(), *l2))
            .collect();

        for (i, snp) in snps.iter().enumerate() {
            let w = w_map.get(snp.as_str()).copied().unwrap_or(1.0);
            all_w_ld.push(w);
            all_snps.push(snp.clone());
            all_chr.push(chr_vals[i]);
            all_bp.push(bp_vals[i]);
            all_annot_rows.push(annot_data[i].clone());
        }

        // Read M_5_50 file
        let m_path = ld_dir.join(format!("{chr}.l2.M_5_50"));
        let m_vals = read_m_annot_file(&m_path)?;
        match m_annot_total {
            Some(ref mut total) => {
                for (t, v) in total.iter_mut().zip(m_vals.iter()) {
                    *t += v;
                }
            }
            None => {
                m_annot_total = Some(m_vals);
            }
        }
    }

    let names = annotation_names.context("no LD score files found")?;
    let m_annot = m_annot_total.context("no M_5_50 files found")?;
    let n_annot = names.len();
    let n_snps = all_snps.len();

    // Build Mat from rows
    let annot_ld = Mat::from_fn(n_snps, n_annot, |i, j| all_annot_rows[i][j]);

    Ok(AnnotLdScores {
        snps: all_snps,
        chr: all_chr,
        bp: all_bp,
        annot_ld,
        w_ld: all_w_ld,
        annotation_names: names,
        m_annot,
    })
}

/// Parsed annotation LD score file data.
// (snp_names, chr, bp, annotation_data[n_snps][n_annot], annotation_names)
type AnnotLdFileData = (Vec<String>, Vec<u32>, Vec<u64>, Vec<Vec<f64>>, Vec<String>);

/// Read a single annotation LD score file.
/// Returns (snp_names, annotation_data[n_snps][n_annot], annotation_names).
fn read_annot_ld_file(path: &Path) -> Result<AnnotLdFileData> {
    let file =
        std::fs::File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    let reader = BufReader::new(GzDecoder::new(file));
    let mut lines = reader.lines();

    let header = lines.next().context("empty LD score file")??;
    let fields: Vec<&str> = if header.contains('\t') {
        header.split('\t').collect()
    } else {
        header.split_whitespace().collect()
    };

    // Find SNP, CHR, BP columns
    let snp_idx = fields
        .iter()
        .position(|&h| h == "SNP")
        .context("SNP column not found")?;
    let chr_idx = fields.iter().position(|&h| h == "CHR");
    let bp_idx = fields.iter().position(|&h| h == "BP");

    // Annotation columns: everything except CHR, SNP, BP, CM, MAF
    let skip_cols: std::collections::HashSet<&str> =
        ["CHR", "SNP", "BP", "CM", "MAF"].iter().copied().collect();
    let annot_indices: Vec<(usize, String)> = fields
        .iter()
        .enumerate()
        .filter(|&(_, name)| !skip_cols.contains(name))
        .map(|(i, &name)| (i, name.to_string()))
        .collect();

    let annotation_names: Vec<String> = annot_indices.iter().map(|(_, n)| n.clone()).collect();
    let n_annot = annot_indices.len();

    let mut snps = Vec::new();
    let mut chrs = Vec::new();
    let mut bps = Vec::new();
    let mut data: Vec<Vec<f64>> = Vec::new();

    for line_result in lines {
        let line = line_result?;
        if line.is_empty() {
            continue;
        }
        let flds: Vec<&str> = if line.contains('\t') {
            line.split('\t').collect()
        } else {
            line.split_whitespace().collect()
        };

        if flds.len() <= snp_idx {
            continue;
        }

        snps.push(flds[snp_idx].to_string());
        chrs.push(
            chr_idx
                .and_then(|i| flds.get(i)?.parse::<u32>().ok())
                .unwrap_or(0),
        );
        bps.push(
            bp_idx
                .and_then(|i| flds.get(i)?.parse::<u64>().ok())
                .unwrap_or(0),
        );
        let mut row = Vec::with_capacity(n_annot);
        for &(col_idx, ref col_name) in &annot_indices {
            let raw = flds
                .get(col_idx)
                .ok_or_else(|| anyhow::anyhow!("Missing column {} in annotation file", col_name))?;
            let val: f64 = raw.parse().map_err(|_| {
                anyhow::anyhow!(
                    "Non-numeric annotation value '{}' in column {}",
                    raw,
                    col_name
                )
            })?;
            row.push(val);
        }
        data.push(row);
    }

    Ok((snps, chrs, bps, data, annotation_names))
}

/// Read a weight LD score file (just SNP and L2 columns).
fn read_weight_ld_file(path: &Path) -> Result<Vec<(String, f64)>> {
    let file =
        std::fs::File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    let reader = BufReader::new(GzDecoder::new(file));
    let mut lines = reader.lines();

    let header = lines.next().context("empty weight file")??;
    let fields: Vec<&str> = if header.contains('\t') {
        header.split('\t').collect()
    } else {
        header.split_whitespace().collect()
    };

    let snp_idx = fields
        .iter()
        .position(|&h| h == "SNP")
        .context("SNP column not found")?;
    let l2_idx = fields
        .iter()
        .position(|&h| h == "L2")
        .context("L2 column not found")?;

    let mut records = Vec::new();
    for line_result in lines {
        let line = line_result?;
        if line.is_empty() {
            continue;
        }
        let flds: Vec<&str> = if line.contains('\t') {
            line.split('\t').collect()
        } else {
            line.split_whitespace().collect()
        };
        if flds.len() <= snp_idx.max(l2_idx) {
            continue;
        }
        let snp = flds[snp_idx].to_string();
        let l2: f64 = flds[l2_idx].parse().unwrap_or(1.0);
        records.push((snp, l2));
    }

    Ok(records)
}

/// Read M_5_50 file with per-annotation counts.
/// Each row has one or more whitespace-separated values (one per annotation).
fn read_m_annot_file(path: &Path) -> Result<Vec<f64>> {
    let content =
        std::fs::read_to_string(path).with_context(|| format!("cannot open {}", path.display()))?;
    content
        .split_whitespace()
        .map(|s| {
            s.parse::<f64>()
                .with_context(|| format!("invalid M value: {s}"))
        })
        .collect()
}

/// Read PLINK `.frq` files for the given chromosomes.
///
/// Expects files at `{frq_dir}/{chr}.frq` with columns: CHR, SNP, A1, A2, MAF, NCHROBS.
/// Returns a map of SNP → MAF for SNPs with MAF in (0.05, 0.95).
pub fn read_frq_files(frq_dir: &Path, chromosomes: &[usize]) -> Result<HashMap<String, f64>> {
    let mut snp_maf = HashMap::new();

    for &chr in chromosomes {
        let path = frq_dir.join(format!("{chr}.frq"));
        if !path.exists() {
            continue;
        }

        let content = std::fs::read_to_string(&path)
            .with_context(|| format!("cannot open {}", path.display()))?;
        let mut lines = content.lines();

        // Parse header to find SNP and MAF columns
        let header = match lines.next() {
            Some(h) => h,
            None => continue,
        };
        let fields: Vec<&str> = header.split_whitespace().collect();
        let snp_idx = fields
            .iter()
            .position(|&h| h == "SNP")
            .context("SNP column not found in .frq file")?;
        let maf_idx = fields
            .iter()
            .position(|&h| h == "MAF")
            .context("MAF column not found in .frq file")?;

        for line in lines {
            let flds: Vec<&str> = line.split_whitespace().collect();
            if flds.len() <= snp_idx.max(maf_idx) {
                continue;
            }
            let maf: f64 = match flds[maf_idx].parse() {
                Ok(v) => v,
                Err(_) => continue,
            };
            // R filters: MAF > 0.05 & MAF < 0.95
            if maf > 0.05 && maf < 0.95 {
                snp_maf.insert(flds[snp_idx].to_string(), maf);
            }
        }
    }

    Ok(snp_maf)
}

/// Filter annotation LD score data to only SNPs present in the frq map.
pub fn filter_annot_by_frq(annot: &mut AnnotLdScores, frq_snps: &HashMap<String, f64>) {
    let keep: Vec<bool> = annot
        .snps
        .iter()
        .map(|s| frq_snps.contains_key(s))
        .collect();
    let n_keep = keep.iter().filter(|&&k| k).count();
    if n_keep == annot.snps.len() {
        return; // nothing to filter
    }

    let n_annot = annot.annot_ld.ncols();
    let mut new_snps = Vec::with_capacity(n_keep);
    let mut new_w_ld = Vec::with_capacity(n_keep);
    let mut new_rows: Vec<Vec<f64>> = Vec::with_capacity(n_keep);

    for (i, &kept) in keep.iter().enumerate() {
        if kept {
            new_snps.push(annot.snps[i].clone());
            new_w_ld.push(annot.w_ld[i]);
            new_rows.push((0..n_annot).map(|j| annot.annot_ld[(i, j)]).collect());
        }
    }

    annot.snps = new_snps;
    annot.w_ld = new_w_ld;
    let n_new = new_rows.len();
    annot.annot_ld = Mat::from_fn(n_new, n_annot, |i, j| new_rows[i][j]);
}
