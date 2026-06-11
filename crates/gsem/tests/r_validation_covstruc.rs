//! R-equivalence tests for the covstruc-derived functions that operate on
//! an LDSC S/V (or a GLS design): summaryGLS and paLDSC.
//!
//! References are produced by `tests/generate_covstruc_reference.R` running
//! R GenomicSEM, and committed as `{summary_gls,paldsc}.json`.

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

fn assert_close(rust: f64, r: f64, tol: f64, msg: &str) {
    let diff = (rust - r).abs();
    assert!(
        diff < tol,
        "{msg}: Rust={rust:.10e} R={r:.10e} diff={diff:.3e} (tol={tol:.0e})"
    );
}

// ── summaryGLS ──────────────────────────────────────────────────────────────

#[test]
fn test_summary_gls_matches_r() {
    let fix = load_fixture("summary_gls");
    let x = json_to_mat(&fix["x"]);
    let y = json_to_vec(&fix["y"]);
    let v = json_to_mat(&fix["v"]);

    let r_betas = json_to_vec(&fix["betas"]);
    let r_se = json_to_vec(&fix["se"]);
    let r_z = json_to_vec(&fix["z"]);
    let r_p = json_to_vec(&fix["pvals"]);

    let res = gsem::stats::gls::summary_gls(&x, &y, &v).expect("GLS should solve");
    assert_eq!(res.beta.len(), r_betas.len(), "GLS coefficient count");

    for i in 0..r_betas.len() {
        assert_close(res.beta[i], r_betas[i], 1e-9, &format!("GLS beta[{i}]"));
        assert_close(res.se[i], r_se[i], 1e-9, &format!("GLS se[{i}]"));
        assert_close(res.z[i], r_z[i], 1e-9, &format!("GLS z[{i}]"));
        assert_close(res.p[i], r_p[i], 1e-9, &format!("GLS p[{i}]"));
    }
}

// ── rgmodel: genetic correlation matrix R and its sampling cov V_R ──────────

#[test]
fn test_rgmodel_matches_r() {
    let fix = load_fixture("rgmodel");
    let s = json_to_mat(&fix["s"]);
    let v = json_to_mat(&fix["v"]);
    let r_r = json_to_mat(&fix["r"]);
    let r_vr = json_to_mat(&fix["v_r"]);

    let res = gsem_sem::rgmodel::run_rgmodel(&s, &v, gsem_sem::EstimationMethod::Dwls)
        .expect("rgmodel should run");

    // Genetic correlation matrix R == cov2cor(S) (deterministic).
    assert_eq!(res.r.nrows(), r_r.nrows(), "rgmodel R dim");
    for i in 0..r_r.nrows() {
        for j in 0..r_r.ncols() {
            assert_close(
                res.r[(i, j)],
                r_r[(i, j)],
                1e-6,
                &format!("rgmodel R[{i},{j}]"),
            );
        }
    }

    // V_R: sampling covariance of the off-diagonal correlations. gsem uses a
    // numerical-Jacobian delta method vs R's sandwich on the standardized fit,
    // so allow a slightly looser tolerance than the (exact) R matrix.
    assert_eq!(res.v_r.nrows(), r_vr.nrows(), "rgmodel V_R dim");
    for i in 0..r_vr.nrows() {
        for j in 0..r_vr.ncols() {
            assert_close(
                res.v_r[(i, j)],
                r_vr[(i, j)],
                1e-5,
                &format!("rgmodel V_R[{i},{j}]"),
            );
        }
    }
}

// ── write.model: factor -> indicator assignment ─────────────────────────────

#[test]
fn test_write_model_structure_matches_r() {
    let fix = load_fixture("write_model");
    let loadings = json_to_mat(&fix["loadings"]);
    let names: Vec<String> = fix["names"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap().to_string())
        .collect();
    let cutoff = fix["cutoff"].as_f64().unwrap();

    // R defaults: fix_resid=TRUE, bifactor=FALSE, mustload=FALSE, common=FALSE.
    let model =
        gsem_sem::write_model::write_model(&loadings, &names, cutoff, true, false, false, false);

    // Parse "F.. =~ NA*Va + Vb" lines into factor -> indicators (gsem's format
    // differs from R's — NA* markers, fixed factor variance — but the
    // factor->indicator assignment is the deterministic, comparable part).
    let mut fac: Vec<String> = Vec::new();
    let mut inds: Vec<Vec<String>> = Vec::new();
    for line in model.lines() {
        let line = line.trim();
        if let Some((lhs, rhs)) = line.split_once("=~") {
            // Skip a bifactor/common "Common_F" line if present (not here).
            fac.push(lhs.trim().to_string());
            let list: Vec<String> = rhs
                .split('+')
                .map(|t| t.trim().replace("NA*", "").trim().to_string())
                .collect();
            inds.push(list);
        }
    }

    let r_fac: Vec<String> = fix["factors"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap().to_string())
        .collect();
    let r_inds: Vec<Vec<String>> = fix["indicators"]
        .as_array()
        .unwrap()
        .iter()
        .map(|row| {
            row.as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_str().unwrap().to_string())
                .collect()
        })
        .collect();

    assert_eq!(fac, r_fac, "write.model factor names");
    assert_eq!(inds, r_inds, "write.model factor->indicator assignment");
}

// ── paLDSC: observed eigenvalue spectrum ────────────────────────────────────

#[test]
fn test_paldsc_observed_eigenvalues_match_r() {
    let fix = load_fixture("paldsc");
    let s = json_to_mat(&fix["s"]);
    let v = json_to_mat(&fix["v"]);
    let r_obs = json_to_vec(&fix["observed_eig"]);

    // The observed eigenvalues are RNG-independent; a tiny n_sim keeps the
    // (unused-here) simulation cheap.
    let res =
        gsem::stats::parallel_analysis::parallel_analysis(&s, &v, 10, 0.95, false, Some(1), None);

    assert_eq!(
        res.observed.len(),
        r_obs.len(),
        "paLDSC observed eigenvalue count: Rust={} R={}",
        res.observed.len(),
        r_obs.len()
    );
    // Compare the sorted spectra (both descending).
    for (i, (&obs, &r)) in res.observed.iter().zip(r_obs.iter()).enumerate() {
        assert_close(obs, r, 1e-9, &format!("paLDSC observed eig[{i}]"));
    }
}
