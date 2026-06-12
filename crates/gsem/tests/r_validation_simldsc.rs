//! R-equivalence test for simLDSC's deterministic per-SNP Z covariance (Sigma).
//!
//! simLDSC's random Z draw (R's `MASS::mvrnorm` on R's RNG stream) cannot be
//! reproduced bit-for-bit cross-platform, but the matrix it draws from is fully
//! deterministic. The reference (`simldsc_synth.json`) runs R GenomicSEM
//! simLDSC's own construction algebra and dumps Sigma for several LD-score
//! values; gsem's `per_snp_z_cov` must reproduce each to numerical precision.

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

fn json_to_vec(val: &Value) -> Vec<f64> {
    val.as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect()
}

#[test]
fn test_simldsc_per_snp_sigma_matches_r() {
    let fix = load_fixture("simldsc_synth");
    let s = json_to_mat(&fix["s"]);
    let n_diag = json_to_vec(&fix["n_diag"]);
    let n_overlap = fix["n_overlap"].as_f64().unwrap();
    let int_vec = json_to_vec(&fix["intercepts"]);
    let r_pheno = json_to_mat(&fix["r_pheno"]);
    let m = fix["m"].as_f64().unwrap();
    let ld_values = json_to_vec(&fix["ld_values"]);

    let k = s.nrows();
    // R's `int` is a per-trait vector -> diagonal intercept matrix.
    let intercepts = Mat::from_fn(k, k, |i, j| if i == j { int_vec[i] } else { 0.0 });

    let r_sigmas = fix["sigma"].as_array().unwrap();
    assert_eq!(r_sigmas.len(), ld_values.len(), "sigma count");

    for (li, &ld) in ld_values.iter().enumerate() {
        let r_sigma = json_to_mat(&r_sigmas[li]);
        let sigma = gsem::stats::simulation::per_snp_z_cov(
            &s,
            &n_diag,
            ld,
            m,
            &intercepts,
            Some(&r_pheno),
            n_overlap,
        );
        for i in 0..k {
            for j in 0..k {
                let diff = (sigma[(i, j)] - r_sigma[(i, j)]).abs();
                assert!(
                    diff < 1e-9,
                    "Sigma(ld={ld})[{i},{j}]: Rust={:.10e} R={:.10e} diff={:.3e}",
                    sigma[(i, j)],
                    r_sigma[(i, j)],
                    diff
                );
            }
        }
    }
}
