//! Reader for FUSION TWAS association output (`.dat` files).
//!
//! Port of R GenomicSEM's `read_fusion()`. FUSION's `assoc_test.R` writes one
//! whitespace-delimited `.dat` file per trait with (among others) the columns:
//!   `FILE` (panel/weight path), `ID` (gene), `TWAS.Z`, `HSQ`,
//!   and `PERM.PV` / `PERM.N` when permutation testing was run.
//!
//! For each trait this converts the TWAS Z-statistic into a standardized
//! gene-expression effect and SE, then merges the traits (inner join on
//! Gene + Panel) into the **merged TWAS format** that [`super::twas_reader`]
//! consumes: `HSQ, Panel, Gene, beta.<trait>, se.<trait>, ...`. This is the
//! on-ramp from raw FUSION output into `multiGene` / `userGWAS(TWAS=TRUE)`.

use std::collections::HashMap;
use std::f64::consts::PI;
use std::path::Path;

use anyhow::{Context, Result, bail};

use super::gwas_reader::open_file_reader;
use std::io::BufRead;

/// One gene/panel row, merged across all traits.
#[derive(Debug, Clone)]
pub struct FusionGene {
    pub gene: String,
    pub panel: String,
    pub hsq: f64,
    /// Standardized effect per trait (length k).
    pub beta: Vec<f64>,
    /// SE per trait (length k).
    pub se: Vec<f64>,
}

/// Merged FUSION result, in the same shape as [`super::twas_reader::TwasSumstats`].
#[derive(Debug)]
pub struct FusionSumstats {
    pub genes: Vec<FusionGene>,
    pub trait_names: Vec<String>,
}

/// Per-trait parsed FUSION row before merging.
struct RawRow {
    gene: String,
    panel: String,
    hsq: f64,
    effect: f64,
    se: f64,
}

/// Read and merge FUSION `.dat` association files into the merged TWAS format.
///
/// * `files`       — one FUSION `.dat` path per trait, in LDSC trait order.
/// * `trait_names` — names for the `beta.*` / `se.*` columns (defaults to
///   `1..=k` when `None`, matching R).
/// * `binary`      — per-trait binary/continuous flag. `None` ⇒ all binary
///   (matching R's default, which also prints a warning).
/// * `n`           — per-trait sample size (used in the standardization
///   denominator). Required for every trait.
/// * `perm`        — when true, derive the Z-statistic from the permutation
///   p-value (`PERM.PV`, `PERM.N`) instead of `TWAS.Z`.
pub fn read_fusion(
    files: &[impl AsRef<Path>],
    trait_names: Option<&[String]>,
    binary: Option<&[bool]>,
    n: &[f64],
    perm: bool,
) -> Result<FusionSumstats> {
    let k = files.len();
    if k == 0 {
        bail!("read_fusion: no files provided");
    }
    if n.len() != k {
        bail!("read_fusion: N has length {} but {k} files given", n.len());
    }

    let names: Vec<String> = match trait_names {
        Some(t) => {
            if t.len() != k {
                bail!(
                    "read_fusion: trait_names has length {} but {k} files",
                    t.len()
                );
            }
            t.to_vec()
        }
        None => (1..=k).map(|i| i.to_string()).collect(),
    };
    // R default: all binary when not specified.
    let binary: Vec<bool> = match binary {
        Some(b) => {
            if b.len() != k {
                bail!("read_fusion: binary has length {} but {k} files", b.len());
            }
            b.to_vec()
        }
        None => vec![true; k],
    };

    // Parse each trait file into per-gene rows keyed by (Gene, Panel).
    let per_trait: Vec<HashMap<(String, String), RawRow>> = files
        .iter()
        .enumerate()
        .map(|(i, f)| parse_fusion_file(f.as_ref(), binary[i], n[i], perm))
        .collect::<Result<Vec<_>>>()?;

    // Inner-join across traits on (Gene, Panel), preserving the first trait's
    // gene order (R merges sequentially with all.x=F, all.y=F).
    let first = &per_trait[0];
    let mut genes: Vec<FusionGene> = Vec::new();
    let mut seen: std::collections::HashSet<(String, String)> = std::collections::HashSet::new();

    for (key, row0) in iter_in_order(files[0].as_ref(), first)? {
        if per_trait[1..].iter().all(|m| m.contains_key(&key)) {
            if !seen.insert(key.clone()) {
                continue; // unique() — drop duplicate Gene/Panel pairs
            }
            let mut beta = Vec::with_capacity(k);
            let mut se = Vec::with_capacity(k);
            beta.push(row0.effect);
            se.push(row0.se);
            for m in &per_trait[1..] {
                let r = &m[&key];
                beta.push(r.effect);
                se.push(r.se);
            }
            genes.push(FusionGene {
                gene: key.0,
                panel: key.1,
                hsq: row0.hsq, // HSQ from the first trait, as in R
                beta,
                se,
            });
        }
    }

    Ok(FusionSumstats {
        genes,
        trait_names: names,
    })
}

/// Iterate the first trait's rows in original file order (R keeps trait-1
/// ordering through the sequential merges). Re-reads the file once.
fn iter_in_order(
    path: &Path,
    map: &HashMap<(String, String), RawRow>,
) -> Result<Vec<((String, String), RawRow)>> {
    // Cheap: reconstruct order from the file's gene/panel sequence.
    let reader = open_file_reader(path)?;
    let mut lines = reader.lines();
    let header = lines.next().context("empty FUSION file")??;
    let cols = FusionCols::from_header(&header)?;

    let mut out = Vec::new();
    let mut emitted: std::collections::HashSet<(String, String)> = std::collections::HashSet::new();
    for line in lines {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        let Some(gene) = fields.get(cols.id).map(|s| s.to_string()) else {
            continue;
        };
        let Some(panel) = fields.get(cols.file).map(|s| extract_panel(s)) else {
            continue;
        };
        let key = (gene, panel);
        if let Some(row) = map.get(&key)
            && emitted.insert(key.clone())
        {
            out.push((
                key,
                RawRow {
                    gene: row.gene.clone(),
                    panel: row.panel.clone(),
                    hsq: row.hsq,
                    effect: row.effect,
                    se: row.se,
                },
            ));
        }
    }
    Ok(out)
}

/// Column indices in a FUSION `.dat` header.
struct FusionCols {
    file: usize,
    id: usize,
    twas_z: usize,
    hsq: usize,
    perm_pv: Option<usize>,
    perm_n: Option<usize>,
}

impl FusionCols {
    fn from_header(header: &str) -> Result<Self> {
        let cols: Vec<String> = header.split_whitespace().map(|s| s.to_string()).collect();
        let find = |name: &str| -> Option<usize> {
            cols.iter().position(|c| c.eq_ignore_ascii_case(name))
        };
        Ok(FusionCols {
            file: find("FILE").context("FUSION file: FILE column not found")?,
            id: find("ID").context("FUSION file: ID column not found")?,
            twas_z: find("TWAS.Z").context("FUSION file: TWAS.Z column not found")?,
            hsq: find("HSQ").context("FUSION file: HSQ column not found")?,
            perm_pv: find("PERM.PV"),
            perm_n: find("PERM.N"),
        })
    }
}

/// Extract the panel id from a FUSION `FILE` path: the last two path
/// components (`.../<panel>/<weight>` → `<panel>/<weight>`), matching R's
/// `sub(".*//([^/]+/[^/]+)$|.*/([^/]+/[^/]+)$", ...)`.
fn extract_panel(file: &str) -> String {
    let parts: Vec<&str> = file.split('/').filter(|s| !s.is_empty()).collect();
    match parts.len() {
        0 => file.to_string(),
        1 => parts[0].to_string(),
        _ => format!("{}/{}", parts[parts.len() - 2], parts[parts.len() - 1]),
    }
}

fn parse_fusion_file(
    path: &Path,
    binary: bool,
    n: f64,
    perm: bool,
) -> Result<HashMap<(String, String), RawRow>> {
    let reader = open_file_reader(path)?;
    let mut lines = reader.lines();
    let header = lines
        .next()
        .with_context(|| format!("empty FUSION file: {}", path.display()))??;
    let cols = FusionCols::from_header(&header)?;

    if perm && (cols.perm_pv.is_none() || cols.perm_n.is_none()) {
        bail!(
            "read_fusion(perm=true) needs PERM.PV and PERM.N columns in {}",
            path.display()
        );
    }

    let mut out = HashMap::new();
    for line in lines {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        let max_idx = cols.twas_z.max(cols.hsq).max(cols.id).max(cols.file);
        if fields.len() <= max_idx {
            continue;
        }

        // R reads "." / "NA" / "" as NA and na.omit()s the row.
        let Some(hsq) = parse_num(fields[cols.hsq]) else {
            continue;
        };
        let gene = fields[cols.id].to_string();
        let panel = extract_panel(fields[cols.file]);

        // Z statistic: TWAS.Z, or the permutation-derived Z when perm=true.
        let z = if perm {
            let pv_raw = parse_num(fields[cols.perm_pv.unwrap()]);
            let pn = parse_num(fields[cols.perm_n.unwrap()]);
            let twas_z = parse_num(fields[cols.twas_z]);
            match (pv_raw, pn, twas_z) {
                (Some(pv), Some(pn), Some(tz)) if pn > 0.0 => {
                    // Clamp PERM.PV away from {0,1} exactly as R does.
                    let pv = if pv == 1.0 { (pn - 1.0) / pn } else { pv };
                    let pv = if pv == 0.0 { 1.0 / pn } else { pv };
                    Some(tz.signum() * chisq1_upper_quantile(pv).sqrt())
                }
                _ => None,
            }
        } else {
            parse_num(fields[cols.twas_z])
        };
        let Some(z) = z else { continue };

        // Standardize: binary traits go through the liability conversion,
        // continuous traits use the observed-scale effect directly.
        let (effect, se) = if binary {
            let denom = ((n / 4.0) * hsq).sqrt();
            if denom == 0.0 || !denom.is_finite() {
                continue;
            }
            let effect = z / denom;
            let se = 1.0 / denom;
            let scale = (effect * effect * hsq + PI * PI / 3.0).sqrt();
            if scale == 0.0 || !scale.is_finite() {
                continue;
            }
            (effect / scale, se / scale)
        } else {
            let denom = (n * hsq).sqrt();
            if denom == 0.0 || !denom.is_finite() {
                continue;
            }
            let effect = z / denom;
            // R: continuous SE = |effect / Z| = 1/denom (Z cancels), but keep
            // the explicit form so a perm-derived Z is handled identically.
            let se = (effect / z).abs();
            (effect, se)
        };

        if !effect.is_finite() || !se.is_finite() {
            continue;
        }
        // Inner join is on (Gene, Panel); last write wins on dup keys, matching
        // R's eventual unique().
        out.insert(
            (gene.clone(), panel.clone()),
            RawRow {
                gene,
                panel,
                hsq,
                effect,
                se,
            },
        );
    }
    Ok(out)
}

/// Parse a numeric field, treating R's NA tokens (".", "NA", "") as missing.
fn parse_num(s: &str) -> Option<f64> {
    let t = s.trim();
    if t.is_empty() || t == "." || t.eq_ignore_ascii_case("NA") {
        return None;
    }
    t.parse::<f64>().ok().filter(|v| v.is_finite())
}

/// Upper-tail quantile of the chi-square(df=1) distribution: the `q` such that
/// P(X > q) = `p`. Equivalent to R's `qchisq(p, 1, lower=FALSE)`. For df=1 this
/// is `(Phi^-1(p/2))^2`.
fn chisq1_upper_quantile(p: f64) -> f64 {
    use statrs::distribution::{ContinuousCDF, Normal};
    let z = Normal::standard().inverse_cdf(p / 2.0);
    z * z
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_dat(dir: &Path, name: &str, body: &str) -> std::path::PathBuf {
        let path = dir.join(name);
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(body.as_bytes()).unwrap();
        path
    }

    #[test]
    fn test_extract_panel() {
        assert_eq!(
            extract_panel("/a/b/GTEx_Brain/ENSG1.wgt"),
            "GTEx_Brain/ENSG1.wgt"
        );
        assert_eq!(
            extract_panel("GTEx_Brain/ENSG1.wgt"),
            "GTEx_Brain/ENSG1.wgt"
        );
        assert_eq!(extract_panel("solo"), "solo");
    }

    #[test]
    fn test_read_fusion_continuous_single_trait() {
        let dir = std::env::temp_dir().join("gsem_fusion_cont");
        std::fs::create_dir_all(&dir).unwrap();
        // Continuous: effect = Z / sqrt(N*HSQ), se = 1/sqrt(N*HSQ).
        let p = write_dat(
            &dir,
            "t1.dat",
            "FILE ID TWAS.Z HSQ\n\
             panelA/g1.wgt G1 2.0 0.25\n\
             panelA/g2.wgt G2 -1.0 0.16\n",
        );
        let res = read_fusion(
            &[p],
            Some(&["X".to_string()]),
            Some(&[false]),
            &[1000.0],
            false,
        )
        .unwrap();
        assert_eq!(res.genes.len(), 2);
        let g1 = &res.genes[0];
        let denom = (1000.0_f64 * 0.25).sqrt();
        assert!(
            (g1.beta[0] - 2.0 / denom).abs() < 1e-12,
            "beta {}",
            g1.beta[0]
        );
        assert!((g1.se[0] - 1.0 / denom).abs() < 1e-12, "se {}", g1.se[0]);
        assert_eq!(g1.gene, "G1");
        assert_eq!(g1.panel, "panelA/g1.wgt");

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_read_fusion_binary_liability() {
        let dir = std::env::temp_dir().join("gsem_fusion_bin");
        std::fs::create_dir_all(&dir).unwrap();
        let p = write_dat(
            &dir,
            "t1.dat",
            "FILE ID TWAS.Z HSQ\npanelA/g1.wgt G1 3.0 0.40\n",
        );
        let n = 8000.0_f64;
        let hsq = 0.40_f64;
        let res = read_fusion(&[p], None, None, &[n], false).unwrap(); // None binary => binary
        let denom = ((n / 4.0) * hsq).sqrt();
        let effect = 3.0 / denom;
        let se = 1.0 / denom;
        let scale = (effect * effect * hsq + PI * PI / 3.0).sqrt();
        assert!((res.genes[0].beta[0] - effect / scale).abs() < 1e-12);
        assert!((res.genes[0].se[0] - se / scale).abs() < 1e-12);

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_read_fusion_inner_join_and_na_omit() {
        let dir = std::env::temp_dir().join("gsem_fusion_join");
        std::fs::create_dir_all(&dir).unwrap();
        // Trait 1 has G1,G2,G3 (G3 HSQ is NA -> dropped). Trait 2 has G1,G2 only.
        // Inner join keeps G1,G2.
        let p1 = write_dat(
            &dir,
            "t1.dat",
            "FILE ID TWAS.Z HSQ\n\
             pA/g1 G1 2.0 0.25\n\
             pA/g2 G2 1.0 0.25\n\
             pA/g3 G3 1.5 .\n",
        );
        let p2 = write_dat(
            &dir,
            "t2.dat",
            "FILE ID TWAS.Z HSQ\n\
             pA/g1 G1 1.0 0.25\n\
             pA/g2 G2 2.0 0.25\n",
        );
        let res = read_fusion(
            &[p1, p2],
            Some(&["A".to_string(), "B".to_string()]),
            Some(&[false, false]),
            &[1000.0, 1000.0],
            false,
        )
        .unwrap();
        let names: Vec<&str> = res.genes.iter().map(|g| g.gene.as_str()).collect();
        assert_eq!(names, vec!["G1", "G2"], "inner join + na.omit");
        assert_eq!(res.genes[0].beta.len(), 2, "two traits merged");

        std::fs::remove_dir_all(&dir).ok();
    }
}
