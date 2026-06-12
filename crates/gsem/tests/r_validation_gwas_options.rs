//! Per-option R-equivalence for the GWAS path (Phase 1 of the coverage plan).
//!
//! The baseline DWLS / GC=standard / fix_measurement=TRUE case is covered by
//! `r_validation.rs::test_user_gwas_per_snp_match_r`. This file varies ONE
//! behaviour-changing argument at a time off that baseline and checks the Rust
//! port matches R GenomicSEM's `userGWAS` output for that option:
//!
//!   * estimation = "ML"          → `EstimationMethod::Ml`
//!   * GC = "conserv"             → `GcMode::Conservative`
//!   * GC = "none"                → `GcMode::None`
//!   * Q_SNP = TRUE               → `q_snp::compute_q_snp` (the heterogeneity stat)
//!   * std.lv = TRUE              → `parse_model(.., std_lv=true)` (factor scaled
//!     by fixing its variance, first loading freed)
//!
//! (fix_measurement=FALSE is intentionally absent: R's free-measurement fit is
//! computationally singular on this subset, so there is no reference to match.)
//!
//! Reference values are produced by the option-matrix block in
//! `tests/generate_gwas_fixture.R` (real `GenomicSEM::userGWAS` calls on the
//! bench PGC subset) and committed as `gwas_options.json`.
//!
//! Run with: `cargo test -p gsem --test r_validation_gwas_options`

use faer::Mat;
use serde_json::Value;
use std::path::PathBuf;

use gsem::gwas::gc_correction::GcMode;
use gsem::gwas::user_gwas::{SnpParamResult, UserGwasConfig, VariantLabel, run_user_gwas};
use gsem_sem::EstimationMethod;
use gsem_sem::syntax::{Op, parse_model};

// ── Helpers (each integration-test file re-declares its own; see r_validation.rs) ──

fn fixtures_dir() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.pop(); // crates/
    p.pop(); // repo root
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
    let nrows = rows.len();
    let ncols = rows[0].len();
    Mat::from_fn(nrows, ncols, |i, j| rows[i][j])
}

fn json_to_vec(val: &Value) -> Vec<f64> {
    val.as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect()
}

/// Shared S/V/I + SNP inputs carried in gwas_options.json.
struct Inputs {
    s: Mat<f64>,
    v: Mat<f64>,
    i_mat: Mat<f64>,
    model: String,
    std_lv_model: String,
    beta_snp: Vec<Vec<f64>>,
    se_snp: Vec<Vec<f64>>,
    var_snp: Vec<f64>,
    snp_ids: Vec<String>,
}

fn load_inputs(fix: &Value) -> Inputs {
    let s = json_to_mat(&fix["s"]);
    let v = json_to_mat(&fix["v"]);
    let i_mat = json_to_mat(&fix["i_mat"]);
    let model = fix["model"].as_str().unwrap().to_string();
    let std_lv_model = fix["std_lv_model"].as_str().unwrap().to_string();

    let snps_arr = fix["snps"].as_array().unwrap();
    let mut beta_snp = Vec::new();
    let mut se_snp = Vec::new();
    let mut var_snp = Vec::new();
    let mut snp_ids = Vec::new();
    for snp in snps_arr {
        snp_ids.push(snp["SNP"].as_str().unwrap().to_string());
        let maf = snp["MAF"].as_f64().unwrap();
        var_snp.push(2.0 * maf * (1.0 - maf));
        beta_snp.push(json_to_vec(&snp["beta"]));
        se_snp.push(json_to_vec(&snp["se"]));
    }
    Inputs {
        s,
        v,
        i_mat,
        model,
        std_lv_model,
        beta_snp,
        se_snp,
        var_snp,
        snp_ids,
    }
}

/// Build a UserGwasConfig off the baseline, applying one mutation closure.
fn config_for(model_str: &str, mutate: impl FnOnce(&mut UserGwasConfig)) -> UserGwasConfig {
    config_for_with(model_str, false, mutate)
}

/// Like `config_for`, but parses the model with the given `std_lv` flag so the
/// std.lv variant (free first loading, factor variance identified by the flag)
/// is built exactly as R's `userGWAS(std.lv=TRUE)` does.
fn config_for_with(
    model_str: &str,
    std_lv: bool,
    mutate: impl FnOnce(&mut UserGwasConfig),
) -> UserGwasConfig {
    let pt = parse_model(model_str, std_lv).unwrap();
    let mut cfg = UserGwasConfig {
        model: pt,
        estimation: EstimationMethod::Dwls,
        gc: GcMode::Standard,
        max_iter: 500,
        smooth_check: false,
        snp_se: None,
        variant_label: VariantLabel::Snp,
        q_snp: false,
        fix_measurement: true,
        num_threads: None,
    };
    mutate(&mut cfg);
    cfg
}

/// Run userGWAS for a given config and return per-SNP results.
fn run(inp: &Inputs, cfg: &UserGwasConfig) -> Vec<gsem::gwas::user_gwas::SnpResult> {
    let beta_refs: Vec<&[f64]> = inp.beta_snp.iter().map(Vec::as_slice).collect();
    let se_refs: Vec<&[f64]> = inp.se_snp.iter().map(Vec::as_slice).collect();
    run_user_gwas(
        cfg,
        &inp.s,
        &inp.v,
        &inp.i_mat,
        &beta_refs,
        &se_refs,
        &inp.var_snp,
        None,
    )
}

/// Find the F1~SNP regression estimate in a SNP result.
fn snp_effect(res: &gsem::gwas::user_gwas::SnpResult) -> Option<&SnpParamResult> {
    res.params
        .iter()
        .find(|p| p.op == Op::Regression && p.lhs == "F1" && p.rhs == "SNP")
}

/// Shared assertion: for an option variant, every converged SNP's |est| must
/// match R's |est| within `tol`. Common-factor orientation is unidentified, so
/// we compare magnitudes. Returns the number of SNPs compared.
fn assert_est_parity(
    inp: &Inputs,
    r_rows: &[Value],
    cfg: &UserGwasConfig,
    tol: f64,
    label: &str,
) -> usize {
    let results = run(inp, cfg);
    assert_eq!(
        results.len(),
        inp.snp_ids.len(),
        "{label}: SNP count mismatch"
    );
    let mut n = 0;
    for (idx, res) in results.iter().enumerate() {
        if !res.converged {
            continue;
        }
        let r_row = &r_rows[idx];
        let r_snp = r_row["SNP"].as_str().unwrap();
        assert_eq!(
            inp.snp_ids[idx], r_snp,
            "{label}: SNP order mismatch at {idx}"
        );
        let r_est = r_row["est"].as_f64().unwrap();
        let p = snp_effect(res).unwrap_or_else(|| panic!("{label}: no F1~SNP param for {r_snp}"));
        let diff = (p.est.abs() - r_est.abs()).abs();
        assert!(
            diff < tol,
            "{label} |est| for {r_snp}: Rust={:.6} R={r_est:.6} diff={diff:.2e} (tol={tol:.0e})",
            p.est
        );
        n += 1;
    }
    assert!(n >= 15, "{label}: too few SNPs converged ({n}/20)");
    n
}

// ── estimation = "ML" ───────────────────────────────────────────────────────

#[test]
fn test_user_gwas_ml_matches_r() {
    let fix = load_fixture("gwas_options");
    let inp = load_inputs(&fix);
    let r_rows = fix["ml"].as_array().unwrap();
    let cfg = config_for(&inp.model, |c| c.estimation = EstimationMethod::Ml);
    // ML and DWLS differ numerically; tolerance is looser than the DWLS baseline.
    assert_est_parity(&inp, r_rows, &cfg, 5e-3, "userGWAS[ML]");
}

// ── GC modes ────────────────────────────────────────────────────────────────

#[test]
fn test_user_gwas_gc_conservative_matches_r() {
    let fix = load_fixture("gwas_options");
    let inp = load_inputs(&fix);
    let r_rows = fix["gc_conservative"].as_array().unwrap();
    let cfg = config_for(&inp.model, |c| c.gc = GcMode::Conservative);
    assert_est_parity(&inp, r_rows, &cfg, 1e-3, "userGWAS[GC=conservative]");
}

#[test]
fn test_user_gwas_gc_none_matches_r() {
    let fix = load_fixture("gwas_options");
    let inp = load_inputs(&fix);
    let r_rows = fix["gc_none"].as_array().unwrap();
    let cfg = config_for(&inp.model, |c| c.gc = GcMode::None);
    assert_est_parity(&inp, r_rows, &cfg, 1e-3, "userGWAS[GC=none]");
}

// ── Q_SNP heterogeneity statistic (covers q_snp.rs, previously 0%) ──────────

#[test]
fn test_user_gwas_q_snp_matches_r() {
    let fix = load_fixture("gwas_options");
    let inp = load_inputs(&fix);
    let r_rows = fix["q_snp"].as_array().unwrap();
    let cfg = config_for(&inp.model, |c| c.q_snp = true);

    // First confirm the SNP effect still matches with Q_SNP on.
    assert_est_parity(&inp, r_rows, &cfg, 1e-3, "userGWAS[Q_SNP].est");

    // Now the actual Q_SNP statistic — this is the path that was 0% covered.
    let results = run(&inp, &cfg);
    let measure = std::env::var("MEASURE").is_ok();

    // Collect ALL (rust, r) pairs first, then assert — so MEASURE shows the full
    // distribution (a systematic same-sign bias would mean a Rust bug; scatter
    // about zero means optimiser noise propagating through the quadratic form).
    let mut pairs: Vec<(String, f64, f64)> = Vec::new();
    for (idx, res) in results.iter().enumerate() {
        if !res.converged {
            continue;
        }
        let r_q = match r_rows[idx]["q_snp"].as_f64() {
            Some(q) if q.is_finite() => q,
            _ => continue,
        };
        let rust_q = res
            .q_snp
            .unwrap_or_else(|| panic!("Q_SNP not computed for {}", inp.snp_ids[idx]));
        pairs.push((inp.snp_ids[idx].clone(), rust_q, r_q));
        // df is an integer; must match exactly.
        if let Some(r_df) = r_rows[idx]["q_snp_df"].as_f64() {
            assert_eq!(
                res.q_snp_df.map(|d| d as f64),
                Some(r_df),
                "Q_SNP df mismatch for {}",
                inp.snp_ids[idx]
            );
        }
    }

    if measure {
        let mut signed_sum = 0.0;
        for (snp, rust_q, r_q) in &pairs {
            let d = rust_q - r_q;
            signed_sum += d;
            eprintln!(
                "MEASURE Q_SNP {snp}: Rust={rust_q:.6} R={r_q:.6} signed={d:+.2e} rel={:.2e}",
                d.abs() / r_q.abs()
            );
        }
        eprintln!(
            "MEASURE Q_SNP: n={} mean_signed={:+.3e} (near 0 ⇒ noise; one-sided ⇒ bias)",
            pairs.len(),
            signed_sum / pairs.len() as f64
        );
    }

    let n_q = pairs.len();
    for (snp, rust_q, r_q) in &pairs {
        let diff = (rust_q - r_q).abs();
        let rel = diff / r_q.abs();
        // With the LDSC intercept floored at 1.0 (matching R's
        // userGWAS.R:156), Q_SNP agrees with R to optimiser precision — the
        // only residual difference is the per-SNP DWLS fit. 1e-4 relative is a
        // tight parity check, not a fudge factor.
        assert!(
            rel < 1e-4,
            "Q_SNP for {snp}: Rust={rust_q:.6} R={r_q:.6} diff={diff:.2e} rel={rel:.2e}"
        );
    }
    assert!(
        n_q >= 10,
        "Too few SNPs had a comparable Q_SNP value ({n_q}); the statistic may not be wired"
    );
}

// ── std.lv = TRUE ───────────────────────────────────────────────────────────

#[test]
fn test_user_gwas_std_lv_matches_r() {
    let fix = load_fixture("gwas_options");
    let inp = load_inputs(&fix);
    let r_rows = fix["std_lv"].as_array().unwrap();
    // The std.lv model frees the first loading and identifies the factor via
    // std_lv=true (no fixed loading, no "F1 ~~ 1*F1"). The F1~SNP effect is on a
    // different scale than the baseline model, so we compare against R's own
    // std.lv reference rather than the baseline.
    let cfg = config_for_with(&inp.std_lv_model, true, |_| {});
    assert_est_parity(&inp, r_rows, &cfg, 1e-3, "userGWAS[std.lv]");
}
