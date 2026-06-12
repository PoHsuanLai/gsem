//! Integration tests validating the LDSC half of the pipeline (munge,
//! sumstats, ldsc) against R GenomicSEM reference output.
//!
//! Inputs are small synthetic GWAS + LD-score files committed under
//! `tests/fixtures/synth/`, simulated from the standard LDSC generative
//! model. Reference outputs are produced by `tests/generate_synthetic_reference.R`
//! running R GenomicSEM on those exact files, and committed as
//! `{ldsc,munge,sumstats}_synth.json`.
//!
//! Run with: `cargo test -p gsem --test r_validation_ldsc`

use std::collections::HashMap;
use std::path::PathBuf;

use faer::Mat;
use serde_json::Value;

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
    Mat::from_fn(rows.len(), rows[0].len(), |i, j| rows[i][j])
}

fn json_to_vec(val: &Value) -> Vec<f64> {
    val.as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect()
}

fn json_to_strs(val: &Value) -> Vec<String> {
    val.as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap().to_string())
        .collect()
}

/// Max-relative-or-absolute closeness check, robust to the small absolute
/// scale of LDSC genetic-covariance entries.
fn assert_close(rust: f64, r: f64, abs_tol: f64, rel_tol: f64, msg: &str) {
    let diff = (rust - r).abs();
    let scale = r.abs().max(1.0);
    assert!(
        diff < abs_tol || diff / scale < rel_tol,
        "{msg}: Rust={rust:.8e} R={r:.8e} absΔ={diff:.3e} relΔ={:.3e}",
        diff / scale
    );
}

// ── ldsc: S / V / I matrices ────────────────────────────────────────────────

#[test]
fn test_ldsc_synth_matches_r() {
    let fix = load_fixture("ldsc_synth");
    let dir = fixtures_dir();

    let munged: Vec<PathBuf> = json_to_strs(&fix["munged_files"])
        .iter()
        .map(|p| dir.join(p))
        .collect();
    let ld_dir = dir.join(fix["ld_dir"].as_str().unwrap());
    let chr = fix["chr"].as_u64().unwrap() as usize;
    let n_blocks = fix["n_blocks"].as_u64().unwrap() as usize;

    let r_s = json_to_mat(&fix["s"]);
    let r_v = json_to_mat(&fix["v"]);
    let r_i = json_to_mat(&fix["i"]);

    // Mirror the binding's orchestration (bindings/r/src/rust/src/lib.rs).
    let trait_data = gsem::io::gwas_reader::load_trait_data(&munged).unwrap();
    let chromosomes: Vec<usize> = (1..=chr).collect();
    let ld_data = gsem::io::ld_reader::read_ld_scores(&ld_dir, &ld_dir, &chromosomes).unwrap();
    let ld_snps: Vec<String> = ld_data.records.iter().map(|r| r.snp.clone()).collect();
    let ld_scores: Vec<f64> = ld_data.records.iter().map(|r| r.l2).collect();

    let k = trait_data.len();
    let config = gsem_ldsc::LdscConfig {
        n_blocks,
        chisq_max: None,
        num_threads: Some(1),
    };
    let result = gsem_ldsc::ldsc(
        &trait_data,
        &vec![None; k],
        &vec![None; k],
        &ld_scores,
        &ld_data.w_ld,
        &ld_snps,
        ld_data.total_m,
        &config,
        None,
    )
    .unwrap();

    // S (genetic covariance): entries are O(1e-4) here, so check relative.
    for i in 0..k {
        for j in 0..k {
            assert_close(
                result.s[(i, j)],
                r_s[(i, j)],
                1e-9,
                1e-4,
                &format!("ldsc S[{i},{j}]"),
            );
        }
    }
    // I (intercepts): O(1) values, tight absolute check.
    for i in 0..k {
        for j in 0..k {
            assert_close(
                result.i_mat[(i, j)],
                r_i[(i, j)],
                1e-4,
                1e-3,
                &format!("ldsc I[{i},{j}]"),
            );
        }
    }
    // V (sampling covariance from the block jackknife): same block
    // partition and data on both sides, so this should match tightly in
    // relative terms despite the tiny absolute magnitudes.
    let kstar = k * (k + 1) / 2;
    for i in 0..kstar {
        for j in 0..kstar {
            assert_close(
                result.v[(i, j)],
                r_v[(i, j)],
                1e-12,
                1e-3,
                &format!("ldsc V[{i},{j}]"),
            );
        }
    }
}

/// Liability-scale ldsc: binary traits with sample/population prevalences, so R
/// applies the observed->liability conversion. Pins `apply_liability_scale`
/// (the prevalence branch the continuous fixture leaves untested) against R.
#[test]
fn test_ldsc_liability_matches_r() {
    let fix = load_fixture("ldsc_liability");
    let dir = fixtures_dir();

    let munged: Vec<PathBuf> = json_to_strs(&fix["munged_files"])
        .iter()
        .map(|p| dir.join(p))
        .collect();
    let ld_dir = dir.join(fix["ld_dir"].as_str().unwrap());
    let chr = fix["chr"].as_u64().unwrap() as usize;
    let n_blocks = fix["n_blocks"].as_u64().unwrap() as usize;

    let r_s = json_to_mat(&fix["s"]);
    let r_v = json_to_mat(&fix["v"]);
    let r_i = json_to_mat(&fix["i"]);

    let sample_prev: Vec<Option<f64>> = json_to_vec(&fix["sample_prev"])
        .into_iter()
        .map(Some)
        .collect();
    let pop_prev: Vec<Option<f64>> = json_to_vec(&fix["population_prev"])
        .into_iter()
        .map(Some)
        .collect();

    let trait_data = gsem::io::gwas_reader::load_trait_data(&munged).unwrap();
    let chromosomes: Vec<usize> = (1..=chr).collect();
    let ld_data = gsem::io::ld_reader::read_ld_scores(&ld_dir, &ld_dir, &chromosomes).unwrap();
    let ld_snps: Vec<String> = ld_data.records.iter().map(|r| r.snp.clone()).collect();
    let ld_scores: Vec<f64> = ld_data.records.iter().map(|r| r.l2).collect();

    let k = trait_data.len();
    let config = gsem_ldsc::LdscConfig {
        n_blocks,
        chisq_max: None,
        num_threads: Some(1),
    };
    let result = gsem_ldsc::ldsc(
        &trait_data,
        &sample_prev,
        &pop_prev,
        &ld_scores,
        &ld_data.w_ld,
        &ld_snps,
        ld_data.total_m,
        &config,
        None,
    )
    .unwrap();

    // S and V are scaled by the liability conversion; I (intercepts) is not.
    for i in 0..k {
        for j in 0..k {
            assert_close(
                result.s[(i, j)],
                r_s[(i, j)],
                1e-9,
                1e-4,
                &format!("ldsc-liab S[{i},{j}]"),
            );
            assert_close(
                result.i_mat[(i, j)],
                r_i[(i, j)],
                1e-4,
                1e-3,
                &format!("ldsc-liab I[{i},{j}]"),
            );
        }
    }
    let kstar = k * (k + 1) / 2;
    for i in 0..kstar {
        for j in 0..kstar {
            assert_close(
                result.v[(i, j)],
                r_v[(i, j)],
                1e-12,
                1e-3,
                &format!("ldsc-liab V[{i},{j}]"),
            );
        }
    }
}

/// chisq.max filter: a low explicit cutoff drops the top hits before the LDSC
/// regression. Pins the `chisq_max = Some(_)` branch that the default (auto)
/// path leaves untested.
#[test]
fn test_ldsc_chisqmax_matches_r() {
    let fix = load_fixture("ldsc_chisqmax");
    let dir = fixtures_dir();

    let munged: Vec<PathBuf> = json_to_strs(&fix["munged_files"])
        .iter()
        .map(|p| dir.join(p))
        .collect();
    let ld_dir = dir.join(fix["ld_dir"].as_str().unwrap());
    let chr = fix["chr"].as_u64().unwrap() as usize;
    let n_blocks = fix["n_blocks"].as_u64().unwrap() as usize;
    let chisq_max = fix["chisq_max"].as_f64().unwrap();

    let r_s = json_to_mat(&fix["s"]);
    let r_v = json_to_mat(&fix["v"]);
    let r_i = json_to_mat(&fix["i"]);

    let trait_data = gsem::io::gwas_reader::load_trait_data(&munged).unwrap();
    let chromosomes: Vec<usize> = (1..=chr).collect();
    let ld_data = gsem::io::ld_reader::read_ld_scores(&ld_dir, &ld_dir, &chromosomes).unwrap();
    let ld_snps: Vec<String> = ld_data.records.iter().map(|r| r.snp.clone()).collect();
    let ld_scores: Vec<f64> = ld_data.records.iter().map(|r| r.l2).collect();

    let k = trait_data.len();
    let config = gsem_ldsc::LdscConfig {
        n_blocks,
        chisq_max: Some(chisq_max),
        num_threads: Some(1),
    };
    let result = gsem_ldsc::ldsc(
        &trait_data,
        &vec![None; k],
        &vec![None; k],
        &ld_scores,
        &ld_data.w_ld,
        &ld_snps,
        ld_data.total_m,
        &config,
        None,
    )
    .unwrap();

    for i in 0..k {
        for j in 0..k {
            assert_close(
                result.s[(i, j)],
                r_s[(i, j)],
                1e-9,
                1e-4,
                &format!("ldsc-chisqmax S[{i},{j}]"),
            );
            assert_close(
                result.i_mat[(i, j)],
                r_i[(i, j)],
                1e-4,
                1e-3,
                &format!("ldsc-chisqmax I[{i},{j}]"),
            );
        }
    }
    let kstar = k * (k + 1) / 2;
    for i in 0..kstar {
        for j in 0..kstar {
            assert_close(
                result.v[(i, j)],
                r_v[(i, j)],
                1e-12,
                1e-3,
                &format!("ldsc-chisqmax V[{i},{j}]"),
            );
        }
    }
}

// ── munge: Z / N per SNP ────────────────────────────────────────────────────

#[test]
fn test_munge_synth_matches_r() {
    let fix = load_fixture("munge_synth");
    let dir = fixtures_dir();

    let raw = dir.join(fix["raw_file"].as_str().unwrap());
    let hm3 = dir.join(fix["hm3_file"].as_str().unwrap());
    let n = fix["n"].as_f64().unwrap();

    let reference = gsem::munge::read_reference(&hm3).unwrap();
    let config = gsem::munge::MungeConfig {
        info_filter: 0.0,
        maf_filter: 0.01,
        n_override: Some(n),
        column_overrides: None,
    };
    let records = gsem::munge::munge_file(&raw, &reference, &config).unwrap();
    let by_snp: HashMap<&str, &gsem::io::gwas_reader::MungedRecord> =
        records.iter().map(|r| (r.snp.as_str(), r)).collect();

    let r_snp = json_to_strs(&fix["snp"]);
    let r_z = json_to_vec(&fix["z"]);
    let r_n = json_to_vec(&fix["n_out"]);
    let r_a1 = json_to_strs(&fix["a1"]);
    let r_a2 = json_to_strs(&fix["a2"]);

    assert_eq!(
        records.len(),
        r_snp.len(),
        "munge SNP count: Rust={} R={}",
        records.len(),
        r_snp.len()
    );

    for idx in 0..r_snp.len() {
        let rec = by_snp
            .get(r_snp[idx].as_str())
            .unwrap_or_else(|| panic!("munge: SNP {} missing from Rust output", r_snp[idx]));
        assert_close(
            rec.z,
            r_z[idx],
            1e-5,
            1e-6,
            &format!("munge Z {}", r_snp[idx]),
        );
        assert_close(
            rec.n,
            r_n[idx],
            1e-9,
            0.0,
            &format!("munge N {}", r_snp[idx]),
        );
        assert_eq!(rec.a1, r_a1[idx], "munge A1 for {}", r_snp[idx]);
        assert_eq!(rec.a2, r_a2[idx], "munge A2 for {}", r_snp[idx]);
    }
}

// ── sumstats: merged per-SNP betas / SEs ────────────────────────────────────

#[test]
fn test_sumstats_synth_matches_r() {
    let fix = load_fixture("sumstats_synth");
    let dir = fixtures_dir();

    let raw: Vec<PathBuf> = json_to_strs(&fix["raw_files"])
        .iter()
        .map(|p| dir.join(p))
        .collect();
    let raw_refs: Vec<&std::path::Path> = raw.iter().map(|p| p.as_path()).collect();
    let ref_file = dir.join(fix["ref_file"].as_str().unwrap());
    let trait_names = json_to_strs(&fix["trait_names"]);
    let n = json_to_vec(&fix["n"]);
    let k = trait_names.len();

    let config = gsem::sumstats::SumstatsConfig {
        info_filter: 0.0,
        maf_filter: 0.01,
        n_overrides: n.iter().map(|&x| Some(x)).collect(),
        se_logit: vec![false; k],
        ols: vec![true; k],
        linprob: vec![false; k],
        keep_indel: false,
        // R GenomicSEM's `sumstats` only removes strand-ambiguous SNPs when
        // ambig=TRUE; its default (ambig=FALSE) KEEPS them. Mirror that here
        // so the numeric merge is compared on the same SNP set R produced.
        // (NB: the gsemr/genomicsem *bindings* invert this default — they map
        // keep_ambig = ambig, so their default drops ambiguous SNPs. See the
        // ambiguous-default regression test below.)
        keep_ambig: true,
        beta_overrides: vec![None; k],
        direct_filter: false,
        num_threads: Some(1),
    };

    let tmp = std::env::temp_dir().join(format!("gsem_sumstats_synth_{}.tsv", std::process::id()));
    gsem::sumstats::merge_sumstats(&raw_refs, &ref_file, &trait_names, &config, &tmp).unwrap();
    let merged = gsem::io::sumstats_reader::read_merged_sumstats(&tmp).unwrap();
    let _ = std::fs::remove_file(&tmp);

    let mut by_snp: HashMap<String, usize> = HashMap::new();
    for i in 0..merged.len() {
        by_snp.insert(merged.snp[i].clone(), i);
    }

    let r_snp = json_to_strs(&fix["snp"]);
    let r_a1 = json_to_strs(&fix["a1"]);
    let r_a2 = json_to_strs(&fix["a2"]);
    let r_beta = json_to_mat(&fix["beta"]);
    let r_se = json_to_mat(&fix["se"]);

    assert_eq!(
        merged.len(),
        r_snp.len(),
        "sumstats SNP count: Rust={} R={}",
        merged.len(),
        r_snp.len()
    );

    // Full equivalence to R for the default OLS mode: merged SNP set, the
    // reference-aligned alleles, and the standardized betas/SEs
    // (beta = Z/sqrt(N*varSNP), se = 1/sqrt(N*varSNP)). The sign of beta can
    // flip relative to R when A1/A2 orientation differs, so compare |beta|.
    let mut compared = 0;
    for idx in 0..r_snp.len() {
        let Some(&ri) = by_snp.get(&r_snp[idx]) else {
            panic!("sumstats: SNP {} missing from Rust output", r_snp[idx]);
        };
        assert_eq!(
            merged.a1_string(ri),
            r_a1[idx],
            "sumstats A1 {}",
            r_snp[idx]
        );
        assert_eq!(
            merged.a2_string(ri),
            r_a2[idx],
            "sumstats A2 {}",
            r_snp[idx]
        );
        let beta = merged.beta_row(ri);
        let se = merged.se_row(ri);
        for t in 0..k {
            assert_close(
                beta[t].abs(),
                r_beta[(idx, t)].abs(),
                1e-7,
                1e-5,
                &format!("sumstats OLS |beta| {} trait {t}", r_snp[idx]),
            );
            assert_close(
                se[t],
                r_se[(idx, t)],
                1e-7,
                1e-5,
                &format!("sumstats OLS se {} trait {t}", r_snp[idx]),
            );
        }
        compared += 1;
    }
    assert!(compared > 100, "too few SNPs compared: {compared}");
}

/// Regression test for the `ambig` semantics. R GenomicSEM's `sumstats`
/// removes strand-ambiguous SNPs (A/T, C/G) only when ambig=TRUE; its
/// default (ambig=FALSE) KEEPS them. The gsemr/genomicsem bindings map
/// keep_ambig = !ambig, so their default must reproduce R's full SNP set.
/// Here we pin the gsem-core contract those bindings rely on:
///   * keep_ambig=true  -> keeps ambiguous -> matches R's default count
///   * keep_ambig=false -> drops exactly the strand-ambiguous SNPs
#[test]
fn test_sumstats_ambiguous_semantics() {
    let fix = load_fixture("sumstats_synth");
    let dir = fixtures_dir();
    let raw: Vec<PathBuf> = json_to_strs(&fix["raw_files"])
        .iter()
        .map(|p| dir.join(p))
        .collect();
    let raw_refs: Vec<&std::path::Path> = raw.iter().map(|p| p.as_path()).collect();
    let ref_file = dir.join(fix["ref_file"].as_str().unwrap());
    let trait_names = json_to_strs(&fix["trait_names"]);
    let n = json_to_vec(&fix["n"]);
    let k = trait_names.len();

    let run = |keep_ambig: bool| -> usize {
        let config = gsem::sumstats::SumstatsConfig {
            info_filter: 0.0,
            maf_filter: 0.01,
            n_overrides: n.iter().map(|&x| Some(x)).collect(),
            se_logit: vec![false; k],
            ols: vec![true; k],
            linprob: vec![false; k],
            keep_indel: false,
            keep_ambig,
            beta_overrides: vec![None; k],
            direct_filter: false,
            num_threads: Some(1),
        };
        let tmp = std::env::temp_dir().join(format!(
            "gsem_ambig_{}_{}.tsv",
            keep_ambig,
            std::process::id()
        ));
        gsem::sumstats::merge_sumstats(&raw_refs, &ref_file, &trait_names, &config, &tmp).unwrap();
        let m = gsem::io::sumstats_reader::read_merged_sumstats(&tmp).unwrap();
        let _ = std::fs::remove_file(&tmp);
        m.len()
    };

    let r_snp = json_to_strs(&fix["snp"]);
    let r_a1 = json_to_strs(&fix["a1"]);
    let r_a2 = json_to_strs(&fix["a2"]);
    let n_ambiguous = r_a1
        .iter()
        .zip(r_a2.iter())
        .filter(|(a1, a2)| {
            matches!(
                (a1.as_str(), a2.as_str()),
                ("A", "T") | ("T", "A") | ("C", "G") | ("G", "C")
            )
        })
        .count();
    assert!(n_ambiguous > 0, "fixture should contain ambiguous SNPs");

    let kept = run(true);
    let dropped = run(false);

    // keep_ambig=true (R's default behaviour) reproduces R's SNP count.
    assert_eq!(
        kept,
        r_snp.len(),
        "keep_ambig=true should match R's default SNP count ({} vs {})",
        kept,
        r_snp.len()
    );
    // keep_ambig=false drops exactly the strand-ambiguous SNPs.
    assert_eq!(
        kept - dropped,
        n_ambiguous,
        "keep_ambig=false should drop exactly the {n_ambiguous} ambiguous SNPs (dropped {})",
        kept - dropped
    );
}

/// Validate every one of R GenomicSEM's `sumstats` standardization modes
/// (OLS, linprob, se.logit, and the default "none") against R, formula for
/// formula, on the same synthetic raw files. The reference for each mode is
/// produced by running R `sumstats` with the corresponding flags.
#[test]
fn test_sumstats_all_modes_match_r() {
    let fix = load_fixture("sumstats_synth");
    let dir = fixtures_dir();
    let raw: Vec<PathBuf> = json_to_strs(&fix["raw_files"])
        .iter()
        .map(|p| dir.join(p))
        .collect();
    let raw_refs: Vec<&std::path::Path> = raw.iter().map(|p| p.as_path()).collect();
    let ref_file = dir.join(fix["ref_file"].as_str().unwrap());
    let trait_names = json_to_strs(&fix["trait_names"]);
    let n = json_to_vec(&fix["n"]);
    let k = trait_names.len();
    let r_snp = json_to_strs(&fix["snp"]);

    // (mode name, ols, linprob, se_logit)
    let modes = [
        ("ols", true, false, false),
        ("linprob", false, true, false),
        ("se_logit", false, false, true),
        ("none", false, false, false),
    ];

    for (name, ols, linprob, se_logit) in modes {
        let config = gsem::sumstats::SumstatsConfig {
            info_filter: 0.0,
            maf_filter: 0.01,
            n_overrides: n.iter().map(|&x| Some(x)).collect(),
            se_logit: vec![se_logit; k],
            ols: vec![ols; k],
            linprob: vec![linprob; k],
            keep_indel: false,
            keep_ambig: true,
            beta_overrides: vec![None; k],
            direct_filter: false,
            num_threads: Some(1),
        };
        let tmp = std::env::temp_dir().join(format!("gsem_ss_{name}_{}.tsv", std::process::id()));
        gsem::sumstats::merge_sumstats(&raw_refs, &ref_file, &trait_names, &config, &tmp).unwrap();
        let merged = gsem::io::sumstats_reader::read_merged_sumstats(&tmp).unwrap();
        let _ = std::fs::remove_file(&tmp);

        let mut by_snp: HashMap<String, usize> = HashMap::new();
        for i in 0..merged.len() {
            by_snp.insert(merged.snp[i].clone(), i);
        }
        let r_beta = json_to_mat(&fix["modes"][name]["beta"]);
        let r_se = json_to_mat(&fix["modes"][name]["se"]);

        assert_eq!(
            merged.len(),
            r_snp.len(),
            "[{name}] SNP count: Rust={} R={}",
            merged.len(),
            r_snp.len()
        );

        for idx in 0..r_snp.len() {
            let ri = by_snp[&r_snp[idx]];
            let beta = merged.beta_row(ri);
            let se = merged.se_row(ri);
            for t in 0..k {
                // Sign may flip with A1/A2 orientation; compare |beta|.
                assert_close(
                    beta[t].abs(),
                    r_beta[(idx, t)].abs(),
                    1e-7,
                    1e-5,
                    &format!("[{name}] |beta| {} trait {t}", r_snp[idx]),
                );
                assert_close(
                    se[t],
                    r_se[(idx, t)],
                    1e-7,
                    1e-5,
                    &format!("[{name}] se {} trait {t}", r_snp[idx]),
                );
            }
        }
    }
}
