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
/// Configuration for common factor GWAS.
pub struct CommonFactorGwasConfig {
    pub estimation: String,
    pub gc: GcMode,
    pub snp_se: Option<f64>,
    pub smooth_check: bool,
}

impl Default for CommonFactorGwasConfig {
    fn default() -> Self {
        Self {
            estimation: "DWLS".to_string(),
            gc: GcMode::Standard,
            snp_se: None,
            smooth_check: false,
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub fn run_common_factor_gwas(
    trait_names: &[String],
    s_ld: &Mat<f64>,
    v_ld: &Mat<f64>,
    i_ld: &Mat<f64>,
    beta_snp: &[Vec<f64>],
    se_snp: &[Vec<f64>],
    var_snp: &[f64],
    cfg: &CommonFactorGwasConfig,
) -> Vec<SnpResult> {
    // Auto-generate model
    let loading = std::iter::once(format!("NA*{}", trait_names[0]))
        .chain(trait_names[1..].iter().cloned())
        .collect::<Vec<_>>()
        .join(" + ");
    let model = format!("F1 =~ {loading}\nF1 ~ SNP\nF1 ~~ 1*F1\nSNP ~~ SNP");

    let config = UserGwasConfig {
        model,
        estimation: cfg.estimation.clone(),
        gc: cfg.gc,
        max_iter: 500,
        std_lv: false,
        smooth_check: cfg.smooth_check,
        snp_se: cfg.snp_se,
        snp_label: "SNP".to_string(),
        q_snp: false,
        fix_measurement: false,
    };

    user_gwas::run_user_gwas(&config, s_ld, v_ld, i_ld, beta_snp, se_snp, var_snp)
}
