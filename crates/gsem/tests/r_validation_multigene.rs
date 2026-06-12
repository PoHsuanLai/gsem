//! R-equivalence test for the deterministic `S_Full` / `V_Full` construction of
//! R GenomicSEM's `multiGene()`.
//!
//! `multiGene` is algorithmically identical to `multiSNP` except that gene
//! heritabilities (`Genes$HSQ`) play the role of the "variances" and the fixed
//! sampling-SE floor is `1e-8` (vs `5e-4` for SNPs). gsem therefore reuses
//! `build_multi_snp_sv` with gene parameters.
//!
//! NOTE: stock R `multiGene` aborts for k >= 2 traits (`object 'V_SNP' not
//! found`, a typo for `V_Gene`). The reference (`multigene_synth.json`) is
//! produced by a minimally-patched `multiGene` and documents R's intended
//! algorithm; gsem matches it to numerical precision.

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
fn test_multigene_s_full_and_v_full_match_r() {
    let fix = load_fixture("multigene_synth");

    let s_ld = json_to_mat(&fix["s_ld"]);
    let v_ld = json_to_mat(&fix["v_ld"]);
    let i_ld = json_to_mat(&fix["i_ld"]);
    let ld = json_to_mat(&fix["ld"]);
    let beta_rows = json_rows(&fix["beta"]);
    let se_rows = json_rows(&fix["se"]);
    // Gene "variance" is heritability (HSQ), not 2*MAF*(1-MAF).
    let var_gene: Vec<f64> = fix["hsq"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect();
    let n_genes = var_gene.len();

    let r_s_full = json_to_mat(&fix["s_full"]);
    let r_v_full = json_to_mat(&fix["v_full"]);

    let beta_refs: Vec<&[f64]> = beta_rows.iter().map(Vec::as_slice).collect();
    let se_refs: Vec<&[f64]> = se_rows.iter().map(Vec::as_slice).collect();

    // multiGene's default GeneSE floor is 1e-8 (vs 5e-4 for multiSNP).
    let (s_full, v_full) = gsem::gwas::multi_snp::build_multi_snp_sv(
        &s_ld, &v_ld, &i_ld, &beta_refs, &se_refs, &var_gene, &ld, n_genes, 1e-8,
    );

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

    assert_eq!(v_full.nrows(), r_v_full.nrows(), "V_Full dim");
    for i in 0..r_v_full.nrows() {
        for j in 0..r_v_full.ncols() {
            assert_close(
                v_full[(i, j)],
                r_v_full[(i, j)],
                1e-14,
                &format!("V_Full[{i},{j}]"),
            );
        }
    }
}
