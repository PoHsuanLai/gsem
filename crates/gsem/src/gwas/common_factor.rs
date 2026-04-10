use faer::Mat;

use gsem_sem::EstimationMethod;
use gsem_sem::syntax;

use super::gc_correction::GcMode;
use super::user_gwas::{self, SnpResult, UserGwasConfig};

/// Identification strategy for the common factor model.
///
/// In SEM, a latent factor's scale must be identified by either:
///
/// - **FixedVariance**: Fix the factor variance to 1 and free all loadings
///   (`F1 =~ NA*V1 + V2 + V3; F1 ~~ 1*F1`). This is the default and is
///   numerically more stable for our optimizer. Sign of the factor (and thus
///   the sign of `F1 ~ SNP`) is arbitrary.
///
/// - **MarkerIndicator**: Fix the first loading to 1 and free the factor
///   variance (`F1 =~ V1 + V2 + V3`). This is lavaan's default and what R
///   GenomicSEM's `commonfactorGWAS` uses. Sign of `F1 ~ SNP` is tied to the
///   first indicator's direction. Use this if you need exact numerical parity
///   with R GenomicSEM.
///
/// Both parameterizations are mathematically equivalent for well-identified
/// models — the implied covariance matrices are identical and `|F1 ~ SNP|` is
/// the same. They can differ numerically on degenerate or Heywood-prone data.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Identification {
    /// Factor variance fixed to 1, all loadings free (default).
    #[default]
    FixedVariance,
    /// First loading fixed to 1, factor variance free (R GenomicSEM style).
    MarkerIndicator,
}

/// Configuration for common factor GWAS.
pub struct CommonFactorGwasConfig {
    /// Estimation method
    pub estimation: EstimationMethod,
    /// Genomic control correction mode
    pub gc: GcMode,
    /// Override for SNP SE (default: 0.0005)
    pub snp_se: Option<f64>,
    /// Log warnings when covariance matrix requires smoothing
    pub smooth_check: bool,
    /// Factor identification strategy. See [`Identification`] for details.
    pub identification: Identification,
}

impl Default for CommonFactorGwasConfig {
    fn default() -> Self {
        Self {
            estimation: EstimationMethod::Dwls,
            gc: GcMode::Standard,
            snp_se: None,
            smooth_check: false,
            identification: Identification::default(),
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
    on_snp_done: Option<&(dyn Fn() + Sync)>,
) -> Vec<SnpResult> {
    // Auto-generate the common factor model based on the chosen identification
    // strategy. See `Identification` for the tradeoffs between the two.
    // Residual variances are auto-added by the parser.
    let model_str = match cfg.identification {
        Identification::FixedVariance => {
            let loading = std::iter::once(format!("NA*{}", trait_names[0]))
                .chain(trait_names[1..].iter().cloned())
                .collect::<Vec<_>>()
                .join(" + ");
            format!("F1 =~ {loading}\nF1 ~ SNP\nF1 ~~ 1*F1")
        }
        Identification::MarkerIndicator => {
            let loading = trait_names.join(" + ");
            format!("F1 =~ {loading}\nF1 ~ SNP")
        }
    };
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
        fix_measurement: true,
        num_threads: None,
    };

    user_gwas::run_user_gwas(&config, s_ld, v_ld, i_ld, beta_snp, se_snp, var_snp, on_snp_done)
}
