use faer::Mat;
use gsem::gwas::user_gwas::SnpResult;
use gsem::io::sumstats_reader::MergedSnp;
use gsem_ldsc::LdscResult;

pub fn ldsc_result_to_json(result: &LdscResult) -> String {
    result.to_json_string().unwrap_or_default()
}

pub fn json_to_ldsc_result(json: &str) -> Option<LdscResult> {
    LdscResult::from_json_string(json).ok()
}

/// Parse a JSON 2D array (e.g. `[[1,2],[3,4]]`) into a faer Mat.
pub fn json_to_mat(json: &str) -> Option<Mat<f64>> {
    let rows: Vec<Vec<f64>> = serde_json::from_str(json).ok()?;
    if rows.is_empty() {
        return Some(Mat::zeros(0, 0));
    }
    let nrows = rows.len();
    let ncols = rows[0].len();
    Some(Mat::from_fn(nrows, ncols, |i, j| rows[i][j]))
}

pub fn snp_results_to_json(results: &[SnpResult], snps: &[MergedSnp]) -> String {
    let entries: Vec<String> = results
        .iter()
        .map(|r| {
            let snp_name = &snps[r.snp_idx].snp;
            let params: Vec<String> = r
                .params
                .iter()
                .map(|p| {
                    format!(
                        "{{\"lhs\":\"{}\",\"op\":\"{}\",\"rhs\":\"{}\",\"est\":{:.6},\"se\":{:.6},\"z\":{:.4},\"p\":{:.6e}}}",
                        p.lhs, p.op, p.rhs, p.est, p.se, p.z_stat, p.p_value
                    )
                })
                .collect();
            format!(
                "{{\"SNP\":\"{}\",\"chisq\":{:.4},\"df\":{},\"converged\":{},\"params\":[{}]}}",
                snp_name, r.chisq, r.chisq_df, r.converged, params.join(",")
            )
        })
        .collect();
    format!("[{}]", entries.join(","))
}
