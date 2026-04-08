use faer::Mat;

use super::gc_correction::GcMode;
use super::user_gwas::{self, SnpResult, UserGwasConfig};

/// Run common factor GWAS.
///
/// This is a simplified version of userGWAS that auto-generates a 1-factor model:
///   F1 =~ NA*V1 + V2 + ... + Vk
///   F1 ~ SNP
///   F1 ~~ 1*F1
///   SNP ~~ SNP
pub fn run_common_factor_gwas(
    trait_names: &[String],
    s_ld: &Mat<f64>,
    v_ld: &Mat<f64>,
    i_ld: &Mat<f64>,
    beta_snp: &[Vec<f64>],
    se_snp: &[Vec<f64>],
    var_snp: &[f64],
    gc: GcMode,
) -> Vec<SnpResult> {
    let _k = trait_names.len();

    // Auto-generate model
    let mut model_lines = Vec::new();

    // Factor loadings: F1 =~ NA*V1 + V2 + ... + Vk
    let mut loading_terms = Vec::new();
    for (i, name) in trait_names.iter().enumerate() {
        if i == 0 {
            loading_terms.push(format!("NA*{name}"));
        } else {
            loading_terms.push(name.clone());
        }
    }
    model_lines.push(format!("F1 =~ {}", loading_terms.join(" + ")));

    // SNP regression
    model_lines.push("F1 ~ SNP".to_string());

    // Fix factor variance
    model_lines.push("F1 ~~ 1*F1".to_string());

    // SNP variance
    model_lines.push("SNP ~~ SNP".to_string());

    let model = model_lines.join("\n");

    let config = UserGwasConfig {
        model,
        estimation: "DWLS".to_string(),
        gc,
        max_iter: 500,
        std_lv: false,
    };

    user_gwas::run_user_gwas(&config, s_ld, v_ld, i_ld, beta_snp, se_snp, var_snp)
}
