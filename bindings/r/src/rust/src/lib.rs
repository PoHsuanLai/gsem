use extendr_api::prelude::*;

mod conversions;

/// Run LDSC pipeline.
///
/// @param trait_files Character vector of .sumstats.gz file paths
/// @param sample_prev Numeric vector of sample prevalences (NA for continuous)
/// @param pop_prev Numeric vector of population prevalences (NA for continuous)
/// @param ld_dir Path to LD score directory
/// @param wld_dir Path to weight LD score directory (or same as ld_dir)
/// @param n_blocks Number of jackknife blocks (default 200)
/// @return JSON string with S, V, I matrices, N vector, and m
#[extendr]
fn ldsc_rust(
    trait_files: Vec<String>,
    sample_prev: Vec<Rfloat>,
    pop_prev: Vec<Rfloat>,
    ld_dir: &str,
    wld_dir: &str,
    n_blocks: i32,
) -> String {
    let n_blocks = n_blocks as usize;

    let sp: Vec<Option<f64>> = sample_prev
        .iter()
        .map(|v| if v.is_na() { None } else { Some(v.inner()) })
        .collect();
    let pp: Vec<Option<f64>> = pop_prev
        .iter()
        .map(|v| if v.is_na() { None } else { Some(v.inner()) })
        .collect();

    let mut trait_data = Vec::new();
    for path in &trait_files {
        let path = std::path::Path::new(path);
        match gsem::io::gwas_reader::read_sumstats(path) {
            Ok(records) => {
                trait_data.push(gsem_ldsc::TraitSumstats {
                    snp: records.iter().map(|r| r.snp.clone()).collect(),
                    z: records.iter().map(|r| r.z).collect(),
                    n: records.iter().map(|r| r.n).collect(),
                    a1: records.iter().map(|r| r.a1.clone()).collect(),
                    a2: records.iter().map(|r| r.a2.clone()).collect(),
                });
            }
            Err(e) => return format!("{{\"error\": \"{e}\"}}"),
        }
    }

    let ld_path = std::path::Path::new(ld_dir);
    let wld_path = std::path::Path::new(wld_dir);
    let chromosomes: Vec<usize> = (1..=22).collect();
    let ld_data = match gsem::io::ld_reader::read_ld_scores(ld_path, wld_path, &chromosomes) {
        Ok(d) => d,
        Err(e) => return format!("{{\"error\": \"{e}\"}}"),
    };

    let ld_snps: Vec<String> = ld_data.records.iter().map(|r| r.snp.clone()).collect();
    let ld_scores: Vec<f64> = ld_data.records.iter().map(|r| r.l2).collect();

    let config = gsem_ldsc::LdscConfig {
        n_blocks,
        chisq_max: None,
    };

    match gsem_ldsc::ldsc(
        &trait_data,
        &sp,
        &pp,
        &ld_scores,
        &ld_data.w_ld,
        &ld_snps,
        ld_data.total_m,
        &config,
    ) {
        Ok(result) => conversions::ldsc_result_to_json(&result),
        Err(e) => format!("{{\"error\": \"{e}\"}}"),
    }
}

/// Fit a user-specified SEM model.
///
/// @param covstruc_json LDSC result as JSON string
/// @param model lavaan-style model syntax
/// @param estimation Estimation method: "DWLS" or "ML"
/// @return JSON string with converged, objective, and parameter estimates
#[extendr]
fn usermodel_rust(covstruc_json: &str, model: &str, estimation: &str) -> String {
    let ldsc_result = match conversions::json_to_ldsc_result(covstruc_json) {
        Some(r) => r,
        None => return "{\"error\": \"failed to parse covstruc JSON\"}".to_string(),
    };

    let pt = match gsem_sem::syntax::parse_model(model, false) {
        Ok(pt) => pt,
        Err(e) => return format!("{{\"error\": \"model parse error: {e}\"}}"),
    };

    let k = ldsc_result.s.nrows();
    let obs_names: Vec<String> = (0..k).map(|i| format!("V{}", i + 1)).collect();
    let mut sem_model = gsem_sem::model::Model::from_partable(&pt, &obs_names);

    let kstar = k * (k + 1) / 2;
    let v_diag: Vec<f64> = (0..kstar).map(|i| ldsc_result.v[(i, i)]).collect();

    let fit = if estimation.to_uppercase() == "ML" {
        gsem_sem::estimator::fit_ml(&mut sem_model, &ldsc_result.s, 1000)
    } else {
        gsem_sem::estimator::fit_dwls(&mut sem_model, &ldsc_result.s, &v_diag, 1000)
    };

    let params_json: Vec<String> = pt
        .rows
        .iter()
        .enumerate()
        .filter(|(_, row)| row.free > 0)
        .map(|(i, row)| {
            let est = fit.params.get(i).copied().unwrap_or(0.0);
            format!(
                "{{\"lhs\":\"{}\",\"op\":\"{}\",\"rhs\":\"{}\",\"est\":{:.6}}}",
                row.lhs, row.op, row.rhs, est
            )
        })
        .collect();

    format!(
        "{{\"converged\":{},\"objective\":{:.6},\"parameters\":[{}]}}",
        fit.converged,
        fit.objective,
        params_json.join(",")
    )
}

/// Munge GWAS summary statistics files.
///
/// @param files Character vector of GWAS file paths
/// @param hm3 Path to HapMap3 SNP list
/// @param trait_names Character vector of trait names
/// @param info_filter INFO score filter threshold
/// @param maf_filter MAF filter threshold
/// @param out_dir Output directory
/// @return Character vector of output file paths
#[extendr]
fn munge_rust(
    files: Vec<String>,
    hm3: &str,
    trait_names: Vec<String>,
    info_filter: f64,
    maf_filter: f64,
    out_dir: &str,
) -> Vec<String> {
    let hm3_path = std::path::Path::new(hm3);
    let reference = match gsem::munge::read_reference(hm3_path) {
        Ok(r) => r,
        Err(e) => {
            log::error!("Failed to read reference: {e}");
            return vec![];
        }
    };

    let config = gsem::munge::MungeConfig {
        info_filter,
        maf_filter,
        n_override: None,
    };

    let mut output_paths = Vec::new();
    for (i, file) in files.iter().enumerate() {
        let name = trait_names.get(i).map(|s| s.as_str()).unwrap_or("trait");
        let out_path = std::path::Path::new(out_dir).join(format!("{name}.sumstats.gz"));

        match gsem::munge::munge_and_write(
            std::path::Path::new(file),
            &reference,
            &config,
            &out_path,
        ) {
            Ok(()) => output_paths.push(out_path.to_string_lossy().to_string()),
            Err(e) => log::error!("Failed to munge {file}: {e}"),
        }
    }
    output_paths
}

extendr_module! {
    mod gsemr;
    fn ldsc_rust;
    fn usermodel_rust;
    fn munge_rust;
}
