//! R-equivalence tests for the deterministic `S_Full` / `V_Full` construction
//! performed by R GenomicSEM's `multiSNP()`.
//!
//! R's `multiSNP` returns the augmented observed-covariance matrix `S_Full`
//! (k traits + f SNPs) and its sampling covariance `V_Full`; the SEM fit is a
//! separate step (`usermodel` / `userGWAS`). gsem's `build_multi_snp_sv` is the
//! comparable unit. The reference is produced by
//! `tests/generate_multisnp_reference.R` and committed as `multisnp_synth.json`.
//!
//! NOTE on the cross-SNP cross-trait block: R's `multiSNP` weights every such
//! cell by a constant LD value (a documented indexing bug). The fixture uses
//! equal off-diagonal LD so R's constant coincides with the correct per-pair
//! LD that gsem uses; the full V_Full therefore matches to numerical precision.

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

fn json_rows(val: &Value) -> Vec<Vec<f64>> {
    val.as_array()
        .unwrap()
        .iter()
        .map(|row| {
            row.as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_f64().unwrap())
                .collect()
        })
        .collect()
}

fn assert_close(rust: f64, r: f64, tol: f64, msg: &str) {
    let diff = (rust - r).abs();
    assert!(
        diff < tol,
        "{msg}: Rust={rust:.10e} R={r:.10e} diff={diff:.3e} (tol={tol:.0e})"
    );
}

#[test]
fn test_multisnp_s_full_and_v_full_match_r() {
    let fix = load_fixture("multisnp_synth");

    let s_ld = json_to_mat(&fix["s_ld"]);
    let v_ld = json_to_mat(&fix["v_ld"]);
    let i_ld = json_to_mat(&fix["i_ld"]);
    let ld = json_to_mat(&fix["ld"]);
    let beta_rows = json_rows(&fix["beta"]);
    let se_rows = json_rows(&fix["se"]);
    let maf: Vec<f64> = fix["maf"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect();
    let n_snps = maf.len();
    let var_snp: Vec<f64> = maf.iter().map(|m| 2.0 * m * (1.0 - m)).collect();

    let r_s_full = json_to_mat(&fix["s_full"]);
    let r_v_full = json_to_mat(&fix["v_full"]);

    let beta_refs: Vec<&[f64]> = beta_rows.iter().map(Vec::as_slice).collect();
    let se_refs: Vec<&[f64]> = se_rows.iter().map(Vec::as_slice).collect();

    // R's multiSNP default SNPSE corresponds to a raw SE of 0.0005.
    let (s_full, v_full) = gsem::gwas::multi_snp::build_multi_snp_sv(
        &s_ld, &v_ld, &i_ld, &beta_refs, &se_refs, &var_snp, &ld, n_snps, 0.0005,
    );

    // S_Full: exact (deterministic construction).
    assert_eq!(s_full.nrows(), r_s_full.nrows(), "S_Full dim");
    for i in 0..r_s_full.nrows() {
        for j in 0..r_s_full.ncols() {
            assert_close(
                s_full[(i, j)],
                r_s_full[(i, j)],
                1e-10,
                &format!("S_Full[{i},{j}]"),
            );
        }
    }

    // V_Full: full sampling-covariance matrix.
    assert_eq!(v_full.nrows(), r_v_full.nrows(), "V_Full dim");
    for i in 0..r_v_full.nrows() {
        for j in 0..r_v_full.ncols() {
            assert_close(
                v_full[(i, j)],
                r_v_full[(i, j)],
                1e-12,
                &format!("V_Full[{i},{j}]"),
            );
        }
    }
}
