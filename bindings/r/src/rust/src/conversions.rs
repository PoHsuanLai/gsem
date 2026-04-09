use faer::Mat;
use gsem::gwas::user_gwas::SnpResult;
use gsem::io::sumstats_reader::MergedSnp;
use gsem_ldsc::LdscResult;

pub fn ldsc_result_to_json(result: &LdscResult) -> String {
    result.to_json_string().unwrap_or_default()
}

/// JSON with additional S_Stand and V_Stand fields (when stand=TRUE).
pub fn ldsc_result_to_json_stand(result: &LdscResult) -> String {
    let s_stand = gsem_matrix::smooth::cov_to_cor(&result.s);

    // Rescale V: V_Stand = diag(scale) * cov2cor(V) * diag(scale)
    // where scale = sqrt(diag(V)) * (vech(S_Stand) / vech(S))
    let k = result.s.nrows();
    let kstar = k * (k + 1) / 2;
    let s_vec = gsem_matrix::vech::vech(&result.s);
    let s_stand_vec = gsem_matrix::vech::vech(&s_stand);
    let scale: Vec<f64> = s_stand_vec
        .iter()
        .zip(s_vec.iter())
        .enumerate()
        .map(|(i, (&st, &orig))| {
            let ratio = if orig.abs() > 1e-30 { st / orig } else { 0.0 };
            result.v[(i, i)].sqrt() * ratio
        })
        .collect();
    let v_cor = gsem_matrix::smooth::cov_to_cor(&result.v);
    let v_stand = Mat::from_fn(kstar, kstar, |i, j| scale[i] * v_cor[(i, j)] * scale[j]);

    // Build JSON with extra fields
    let base = result.to_json_string().unwrap_or_default();
    // Insert S_Stand and V_Stand before the closing }
    if let Some(pos) = base.rfind('}') {
        let s_stand_json = mat_to_json_2d(&s_stand);
        let v_stand_json = mat_to_json_2d(&v_stand);
        format!(
            "{},\"s_stand\":{},\"v_stand\":{}}}",
            &base[..pos],
            s_stand_json,
            v_stand_json
        )
    } else {
        base
    }
}

fn mat_to_json_2d(mat: &Mat<f64>) -> String {
    let rows: Vec<String> = (0..mat.nrows())
        .map(|i| {
            let vals: Vec<String> = (0..mat.ncols())
                .map(|j| format!("{:.15e}", mat[(i, j)]))
                .collect();
            format!("[{}]", vals.join(","))
        })
        .collect();
    format!("[{}]", rows.join(","))
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
