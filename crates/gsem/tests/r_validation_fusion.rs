//! R-equivalence for the FUSION/TWAS `.dat` reader (`gsem::io::fusion_reader`),
//! the port of R GenomicSEM's `read_fusion()`.
//!
//! The fixture is produced by `tests/generate_fusion_reference.R`, which writes
//! a pair of FUSION `.dat` files (committed under `tests/fixtures/fusion/`) and
//! the merged output of the real `GenomicSEM::read_fusion` on them. This test
//! replays the same files through the Rust port and asserts the standardized
//! beta/SE match R exactly.
//!
//! Run with: `cargo test -p gsem --test r_validation_fusion`

use std::path::PathBuf;

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

fn json_strs(v: &Value) -> Vec<String> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_str().unwrap().to_string())
        .collect()
}

fn json_f64s(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect()
}

fn json_rows(v: &Value) -> Vec<Vec<f64>> {
    v.as_array().unwrap().iter().map(json_f64s).collect()
}

#[test]
fn test_read_fusion_matches_r() {
    let fix = load_fixture("fusion_read");
    let dir = fixtures_dir();

    let files: Vec<PathBuf> = json_strs(&fix["files"])
        .iter()
        .map(|p| dir.join(p))
        .collect();
    let trait_names = json_strs(&fix["trait_names"]);
    let binary: Vec<bool> = fix["binary"]
        .as_array()
        .unwrap()
        .iter()
        .map(|b| b.as_bool().unwrap())
        .collect();
    let n = json_f64s(&fix["n"]);

    let res =
        gsem::io::fusion_reader::read_fusion(&files, Some(&trait_names), Some(&binary), &n, false)
            .expect("read_fusion should succeed");

    let r_gene = json_strs(&fix["gene"]);
    let r_panel = json_strs(&fix["panel"]);
    let r_hsq = json_f64s(&fix["hsq"]);
    let r_beta = json_rows(&fix["beta"]);
    let r_se = json_rows(&fix["se"]);

    assert_eq!(res.trait_names, trait_names, "trait names");
    assert_eq!(
        res.genes.len(),
        r_gene.len(),
        "gene count: Rust={} R={}",
        res.genes.len(),
        r_gene.len()
    );

    for (i, g) in res.genes.iter().enumerate() {
        assert_eq!(g.gene, r_gene[i], "gene[{i}]");
        assert_eq!(g.panel, r_panel[i], "panel[{i}]");
        assert!((g.hsq - r_hsq[i]).abs() < 1e-12, "HSQ[{i}]");
        for t in 0..trait_names.len() {
            assert!(
                (g.beta[t] - r_beta[i][t]).abs() < 1e-9,
                "beta[{i}][{t}]: Rust={} R={}",
                g.beta[t],
                r_beta[i][t]
            );
            assert!(
                (g.se[t] - r_se[i][t]).abs() < 1e-9,
                "se[{i}][{t}]: Rust={} R={}",
                g.se[t],
                r_se[i][t]
            );
        }
    }
}
