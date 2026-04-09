use faer::Mat;

use gsem_sem::EstimationMethod;
use gsem_sem::syntax;

use super::gc_correction::GcMode;
use super::user_gwas::{self, SnpResult, UserGwasConfig};

/// Configuration for common factor GWAS.
///
/// This is a simplified version of userGWAS that auto-generates a 1-factor model:
///   F1 =~ NA*V1 + V2 + ... + Vk
///   F1 ~ SNP
///   F1 ~~ 1*F1
///   SNP ~~ SNP
pub struct CommonFactorGwasConfig {
    /// Estimation method
    pub estimation: EstimationMethod,
    /// Genomic control correction mode
    pub gc: GcMode,
    /// Override for SNP SE (default: 0.0005)
    pub snp_se: Option<f64>,
    /// Log warnings when covariance matrix requires smoothing
    pub smooth_check: bool,
}

impl Default for CommonFactorGwasConfig {
    fn default() -> Self {
        Self {
            estimation: EstimationMethod::Dwls,
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
    // Auto-generate and parse model
    let loading = std::iter::once(format!("NA*{}", trait_names[0]))
        .chain(trait_names[1..].iter().cloned())
        .collect::<Vec<_>>()
        .join(" + ");
    let model_str = format!("F1 =~ {loading}\nF1 ~ SNP\nF1 ~~ 1*F1\nSNP ~~ SNP");
    // This parse cannot fail since we generate the syntax ourselves
    let model =
        syntax::parse_model(&model_str, false).expect("auto-generated model syntax is invalid");

    let config = UserGwasConfig {
        model,
        estimation: cfg.estimation,
        gc: cfg.gc,
        max_iter: 500,
        smooth_check: cfg.smooth_check,
        snp_se: cfg.snp_se,
        variant_label: user_gwas::VariantLabel::Snp,
        q_snp: false,
        fix_measurement: false,
    };

    user_gwas::run_user_gwas(&config, s_ld, v_ld, i_ld, beta_snp, se_snp, var_snp)
}
