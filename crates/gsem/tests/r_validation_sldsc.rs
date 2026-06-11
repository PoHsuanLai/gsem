//! Stratified LDSC (s_ldsc) R-equivalence with OVERLAPPING annotations.
//!
//! Reference produced by `tests/generate_sldsc_reference.R` running R
//! GenomicSEM on synthetic baselineLD-style inputs (a `base` annotation
//! covering all SNPs + two overlapping binary categories `A`, `B`).
//!
//! R returns per-annotation `S` (overlap-weighted partitioned heritability,
//! `overlap · cats`) and `S_Tau` (raw `tau · M`). gsem must match both.

use std::path::PathBuf;

use faer::Mat;
use serde_json::Value;

fn fixtures_dir() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.pop();
    p.pop();
    p.push("tests");
    p.push("fixtures");
    p
}

fn load_fixture(name: &str) -> Value {
    let path = fixtures_dir().join(format!("{name}.json"));
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("failed to read fixture {}: {e}", path.display()));
    serde_json::from_str(&data).unwrap_or_else(|e| panic!("invalid JSON in {name}.json: {e}"))
}

fn json_to_mat(val: &Value) -> Mat<f64> {
    let rows: Vec<Vec<f64>> = val
        .as_array()
        .unwrap()
        .iter()
        .map(|row| {
            row.as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_f64().unwrap())
                .collect()
        })
        .collect();
    Mat::from_fn(rows.len(), rows[0].len(), |i, j| rows[i][j])
}

fn json_to_strs(val: &Value) -> Vec<String> {
    val.as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap().to_string())
        .collect()
}

/// Per-annotation list of k×k matrices.
fn json_to_mats(val: &Value) -> Vec<Mat<f64>> {
    val.as_array().unwrap().iter().map(json_to_mat).collect()
}

fn assert_close(rust: f64, r: f64, abs_tol: f64, rel_tol: f64, msg: &str) {
    let diff = (rust - r).abs();
    let scale = r.abs().max(1.0);
    assert!(
        diff < abs_tol || diff / scale < rel_tol,
        "{msg}: Rust={rust:.8e} R={r:.8e} absΔ={diff:.3e} relΔ={:.3e}",
        diff / scale
    );
}

fn run_gsem_sldsc() -> (gsem_ldsc::stratified::StratifiedLdscResult, Value) {
    let fix = load_fixture("sldsc_synth");
    let dir = fixtures_dir();
    let munged: Vec<PathBuf> = json_to_strs(&fix["munged_files"])
        .iter()
        .map(|p| dir.join(p))
        .collect();
    let ld_dir = dir.join(fix["ld_dir"].as_str().unwrap());
    let wld_dir = dir.join(fix["wld_dir"].as_str().unwrap());
    let chr = fix["chr"].as_u64().unwrap() as usize;
    let n_blocks = fix["n_blocks"].as_u64().unwrap() as usize;

    let trait_data = gsem::io::gwas_reader::load_trait_data(&munged).unwrap();
    let chromosomes: Vec<usize> = (1..=chr).collect();
    let frq_dir = dir.join(fix["frq_dir"].as_str().unwrap());
    let annot =
        gsem_ldsc::annot_reader::read_annot_ld_scores(&ld_dir, &wld_dir, &chromosomes).unwrap();
    let cross = gsem_ldsc::annot_reader::read_annot_cross(&ld_dir, &frq_dir, &chromosomes).unwrap();

    let k = trait_data.len();
    let cfg = gsem_ldsc::stratified::StratifiedLdscConfig {
        n_blocks,
        rm_flank: false,
        flank_kb: 500,
    };
    let res = gsem_ldsc::stratified::s_ldsc(
        &trait_data,
        &vec![None; k],
        &vec![None; k],
        &annot.annot_ld,
        &annot.w_ld,
        &annot.snps,
        &annot.annotation_names,
        &annot.m_annot,
        &cfg,
        Some(&annot.chr),
        Some(&annot.bp),
        Some(&cross),
    )
    .unwrap();
    (res, fix)
}

// ── Raw per-annotation tau (S_Tau / V_Tau): validates the partitioned
//    regression itself, independent of the overlap weighting. ─────────────────

#[test]
fn test_sldsc_tau_matches_r() {
    let (res, fix) = run_gsem_sldsc();
    let r_s_tau = json_to_mats(&fix["s_tau"]);
    assert_eq!(res.s_tau.len(), r_s_tau.len(), "annotation count (S_Tau)");

    for (a, (rust, r)) in res.s_tau.iter().zip(r_s_tau.iter()).enumerate() {
        for i in 0..r.nrows() {
            for j in 0..r.ncols() {
                assert_close(
                    rust[(i, j)],
                    r[(i, j)],
                    1e-9,
                    1e-4,
                    &format!("s_tau annot {a} [{i},{j}]"),
                );
            }
        }
    }
}

// ── Overlap-weighted partitioned heritability S and its sampling cov V ───────

#[test]
fn test_sldsc_overlap_s_v_match_r() {
    let (res, fix) = run_gsem_sldsc();
    let r_s = json_to_mats(&fix["s"]);
    let r_v = json_to_mats(&fix["v"]);
    assert_eq!(res.s_annot.len(), r_s.len(), "annotation count (S)");

    for (a, (rs, rr)) in res.s_annot.iter().zip(r_s.iter()).enumerate() {
        // S: overlap-weighted partitioned (co)heritability.
        for i in 0..rr.nrows() {
            for j in 0..rr.ncols() {
                assert_close(
                    rs[(i, j)],
                    rr[(i, j)],
                    1e-9,
                    1e-4,
                    &format!("S annot {a} [{i},{j}]"),
                );
            }
        }
    }
    for (a, (rv, rrv)) in res.v_annot.iter().zip(r_v.iter()).enumerate() {
        // V: sampling covariance of S (jackknife; delta-method block transform).
        for i in 0..rrv.nrows() {
            for j in 0..rrv.ncols() {
                assert_close(
                    rv[(i, j)],
                    rrv[(i, j)],
                    1e-10,
                    1e-3,
                    &format!("V annot {a} [{i},{j}]"),
                );
            }
        }
    }
}

// ── Intercept matrix (same as standard LDSC) ────────────────────────────────

#[test]
fn test_sldsc_intercept_matches_r() {
    let (res, fix) = run_gsem_sldsc();
    let r_i = json_to_mat(&fix["i"]);
    for i in 0..r_i.nrows() {
        for j in 0..r_i.ncols() {
            assert_close(
                res.i_mat[(i, j)],
                r_i[(i, j)],
                1e-4,
                1e-3,
                &format!("intercept [{i},{j}]"),
            );
        }
    }
}
