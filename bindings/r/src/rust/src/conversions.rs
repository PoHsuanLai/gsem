use gsem_ldsc::LdscResult;

pub fn ldsc_result_to_json(result: &LdscResult) -> String {
    result.to_json_string().unwrap_or_default()
}

pub fn json_to_ldsc_result(json: &str) -> Option<LdscResult> {
    LdscResult::from_json_string(json).ok()
}
