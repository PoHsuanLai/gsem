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

// ── enrich: model-based functional enrichment ───────────────────────────────

#[test]
fn test_enrich_matches_r() {
    let fix = load_fixture("enrich_synth");
    let s_list: Vec<Mat<f64>> = fix["s"]
        .as_array()
        .unwrap()
        .iter()
        .map(json_to_mat)
        .collect();
    let v_list: Vec<Mat<f64>> = fix["v"]
        .as_array()
        .unwrap()
        .iter()
        .map(json_to_mat)
        .collect();
    let prop = json_to_vec(&fix["prop"]);
    let obs_names: Vec<String> = fix["traits"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap().to_string())
        .collect();
    let annot_names: Vec<String> = fix["annot_names"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap().to_string())
        .collect();
    let model = fix["model"].as_str().unwrap();
    let params: Vec<String> = fix["params"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap().to_string())
        .collect();

    let r_enrich = json_to_vec(&fix["enrichment"]);
    let r_se = json_to_vec(&fix["enrichment_se"]);
    let r_p = json_to_vec(&fix["enrichment_p"]);

    let res = gsem_sem::enrich_model::model_enrichment(
        &s_list,
        &v_list,
        &prop,
        &annot_names,
        &obs_names,
        model,
        &params,
        gsem_sem::enrich_model::FixMode::Regressions,
        gsem_sem::EstimationMethod::Dwls,
    )
    .expect("enrich should run");

    // Single target param (F1~~F1); compare per-annotation enrichment/SE/p.
    let enr = &res.enrichment[0];
    let se = &res.se[0];
    let p = &res.p[0];
    for a in 0..annot_names.len() {
        assert_close(enr[a], r_enrich[a], 1e-4, &format!("enrichment annot {a}"));
        assert_close(se[a], r_se[a], 1e-5, &format!("enrichment_se annot {a}"));
        assert_close(p[a], r_p[a], 1e-4, &format!("enrichment_p annot {a}"));
    }
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

// ── subSV: vech-position subsetting of S and V ──────────────────────────────

#[test]
fn test_subsv_matches_r() {
    use gsem_matrix::vech::{SubsetType, subset_sv};

    let fix = load_fixture("subsv");
    let s = json_to_mat(&fix["s"]);
    let v = json_to_mat(&fix["v"]);
    let r_corr = json_to_mat(&fix["r_corr"]);
    let v_r = json_to_mat(&fix["v_r"]);

    let to_idx =
        |val: &Value| -> Vec<usize> { json_to_vec(val).iter().map(|&x| x as usize).collect() };
    let index_s = to_idx(&fix["index_s"]);
    let index_r = to_idx(&fix["index_r"]);

    let r_sub_s = json_to_vec(&fix["sub_s"]);
    let r_sub_v = json_to_mat(&fix["sub_v"]);
    let r_sub_s_r = json_to_vec(&fix["sub_s_r"]);
    let r_sub_v_r = json_to_mat(&fix["sub_v_r"]);

    // TYPE = "S": full lower triangle incl. diagonal.
    let out = subset_sv(&s, &v, &index_s, SubsetType::WithDiagonal).unwrap();
    assert_eq!(out.sub_s.len(), r_sub_s.len(), "subSV[S] subS length");
    for (i, (&rust, &r)) in out.sub_s.iter().zip(r_sub_s.iter()).enumerate() {
        assert_close(rust, r, 1e-12, &format!("subSV[S] subS[{i}]"));
    }
    for i in 0..r_sub_v.nrows() {
        for j in 0..r_sub_v.ncols() {
            assert_close(
                out.sub_v[(i, j)],
                r_sub_v[(i, j)],
                1e-12,
                &format!("subSV[S] subV[{i},{j}]"),
            );
        }
    }

    // TYPE = "R": strict lower triangle (off-diagonal numbering) on the
    // correlation matrix.
    let out_r = subset_sv(&r_corr, &v_r, &index_r, SubsetType::OffDiagonal).unwrap();
    assert_eq!(out_r.sub_s.len(), r_sub_s_r.len(), "subSV[R] subS length");
    for (i, (&rust, &r)) in out_r.sub_s.iter().zip(r_sub_s_r.iter()).enumerate() {
        assert_close(rust, r, 1e-12, &format!("subSV[R] subS[{i}]"));
    }
    for i in 0..r_sub_v_r.nrows() {
        for j in 0..r_sub_v_r.ncols() {
            assert_close(
                out_r.sub_v[(i, j)],
                r_sub_v_r[(i, j)],
                1e-12,
                &format!("subSV[R] subV[{i},{j}]"),
            );
        }
    }
}

// ── summaryGLSbands: GLS confidence-band data ───────────────────────────────

#[test]
fn test_summary_gls_bands_matches_r() {
    let fix = load_fixture("gls_bands");
    let predictors = json_to_vec(&fix["predictors"]);
    let y = json_to_vec(&fix["y"]);
    let v_y = json_to_mat(&fix["v_y"]);
    let intervals = fix["intervals"].as_u64().unwrap() as usize;
    let band_size = fix["band_size"].as_f64().unwrap();

    let r_betas = json_to_vec(&fix["betas"]);
    let r_grid = json_to_vec(&fix["grid"]);
    let r_line = json_to_vec(&fix["line"]);
    let r_band_se = json_to_vec(&fix["band_se"]);
    let r_upper = json_to_vec(&fix["upper"]);
    let r_lower = json_to_vec(&fix["lower"]);

    let bands = gsem::stats::gls::summary_gls_bands(
        &predictors,
        &y,
        &v_y,
        true,  // INTERCEPT
        false, // QUAD
        &[],   // no CONTROLVARS
        intervals,
        band_size,
    )
    .expect("summary_gls_bands should solve");

    for (i, (&rust, &r)) in bands.fit.beta.iter().zip(r_betas.iter()).enumerate() {
        assert_close(rust, r, 1e-9, &format!("GLSbands beta[{i}]"));
    }
    for (i, (&rust, &r)) in bands.grid.iter().zip(r_grid.iter()).enumerate() {
        assert_close(rust, r, 1e-9, &format!("GLSbands grid[{i}]"));
    }
    for (i, (&rust, &r)) in bands.line.iter().zip(r_line.iter()).enumerate() {
        assert_close(rust, r, 1e-9, &format!("GLSbands line[{i}]"));
    }
    for (i, (&rust, &r)) in bands.band_se.iter().zip(r_band_se.iter()).enumerate() {
        assert_close(rust, r, 1e-9, &format!("GLSbands band_se[{i}]"));
    }
    for (i, (&rust, &r)) in bands.upper.iter().zip(r_upper.iter()).enumerate() {
        assert_close(rust, r, 1e-9, &format!("GLSbands upper[{i}]"));
    }
    for (i, (&rust, &r)) in bands.lower.iter().zip(r_lower.iter()).enumerate() {
        assert_close(rust, r, 1e-9, &format!("GLSbands lower[{i}]"));
    }
}
