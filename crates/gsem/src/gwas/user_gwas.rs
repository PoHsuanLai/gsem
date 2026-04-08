use faer::Mat;
use gsem_matrix::smooth;
use gsem_sem::estimator;
use gsem_sem::model::Model;
use gsem_sem::sandwich;
use gsem_sem::syntax::{self, ParTable};
use rayon::prelude::*;
use statrs::distribution::ContinuousCDF;

use super::add_snps;
use super::gc_correction::GcMode;

/// Result for a single SNP from userGWAS.
#[derive(Debug, Clone)]
pub struct SnpResult {
    pub snp_idx: usize,
    pub params: Vec<SnpParamResult>,
    pub chisq: f64,
    pub chisq_df: usize,
    pub converged: bool,
    pub warning: Option<String>,
}

/// A single parameter result for one SNP.
#[derive(Debug, Clone)]
pub struct SnpParamResult {
    pub lhs: String,
    pub op: String,
    pub rhs: String,
    pub est: f64,
    pub se: f64,
    pub z_stat: f64,
    pub p_value: f64,
}

/// Configuration for userGWAS.
#[derive(Debug, Clone)]
pub struct UserGwasConfig {
    pub model: String,
    pub estimation: String,
    pub gc: GcMode,
    pub max_iter: usize,
    pub std_lv: bool,
    /// Override for SNP SE (default: 0.0005, treating MAF as fixed).
    /// Matches R's `SNPSE` parameter.
    pub snp_se: Option<f64>,
}

impl Default for UserGwasConfig {
    fn default() -> Self {
        Self {
            model: String::new(),
            estimation: "DWLS".to_string(),
            gc: GcMode::Standard,
            max_iter: 500,
            std_lv: false,
            snp_se: None,
        }
    }
}

/// Run user-specified model GWAS across all SNPs.
///
/// Parallelized via rayon.
pub fn run_user_gwas(
    config: &UserGwasConfig,
    s_ld: &Mat<f64>,
    v_ld: &Mat<f64>,
    i_ld: &Mat<f64>,
    beta_snp: &[Vec<f64>], // n_snps × k
    se_snp: &[Vec<f64>],   // n_snps × k
    var_snp: &[f64],       // n_snps
) -> Vec<SnpResult> {
    let n_snps = var_snp.len();
    let k = s_ld.nrows();

    // Parse model once
    let pt = match syntax::parse_model(&config.model, config.std_lv) {
        Ok(pt) => pt,
        Err(e) => {
            return vec![SnpResult {
                snp_idx: 0,
                params: vec![],
                chisq: f64::NAN,
                chisq_df: 0,
                converged: false,
                warning: Some(format!("model parse error: {e}")),
            }];
        }
    };

    // Process SNPs in parallel
    (0..n_snps)
        .into_par_iter()
        .map(|i| {
            process_single_snp(
                i,
                &pt,
                config,
                s_ld,
                v_ld,
                i_ld,
                &beta_snp[i],
                &se_snp[i],
                var_snp[i],
                k,
            )
        })
        .collect()
}

/// Process a single SNP.
fn process_single_snp(
    snp_idx: usize,
    pt: &ParTable,
    config: &UserGwasConfig,
    s_ld: &Mat<f64>,
    v_ld: &Mat<f64>,
    i_ld: &Mat<f64>,
    beta_snp: &[f64],
    se_snp: &[f64],
    var_snp: f64,
    k: usize,
) -> SnpResult {
    // Build S_Full and V_Full
    let mut s_full = add_snps::build_s_full(s_ld, beta_snp, var_snp, k);
    // R: varSNPSE2 <- (.0005)^2 — small constant, treating MAF as fixed
    let snp_se = config.snp_se.unwrap_or(0.0005);
    let var_snp_se2 = snp_se.powi(2);
    let mut v_full = add_snps::build_v_full(v_ld, se_snp, var_snp, var_snp_se2, i_ld, config.gc, k);

    // Smooth if needed
    smooth::smooth_if_needed(&mut s_full);
    smooth::smooth_if_needed(&mut v_full);

    // Build model
    let mut obs_names = vec!["SNP".to_string()];
    obs_names.extend((0..k).map(|i| format!("V{}", i + 1)));

    let mut model = Model::from_partable(pt, &obs_names);

    // Construct diagonal weight matrix
    let kstar_full = (k + 1) * (k + 2) / 2;
    let v_diag: Vec<f64> = (0..kstar_full).map(|i| v_full[(i, i)]).collect();

    // Fit model
    let fit = if config.estimation.eq_ignore_ascii_case("ML") {
        estimator::fit_ml(&mut model, &s_full, config.max_iter)
    } else {
        estimator::fit_dwls(&mut model, &s_full, &v_diag, config.max_iter)
    };

    // Compute sandwich SEs
    let w_diag = Mat::from_fn(kstar_full, kstar_full, |i, j| {
        if i == j && v_diag[i] > 1e-30 {
            1.0 / v_diag[i]
        } else {
            0.0
        }
    });

    let (se, _ohtt) = sandwich::sandwich_se(&mut model, &w_diag, &v_full);

    // Build parameter results
    let params: Vec<SnpParamResult> = pt
        .rows
        .iter()
        .zip(fit.params.iter().chain(std::iter::repeat(&0.0)))
        .enumerate()
        .filter(|(_, (row, _))| row.free > 0)
        .map(|(i, (row, &est))| {
            let se_val = se.get(i).copied().unwrap_or(0.0);
            let z = if se_val > 0.0 { est / se_val } else { 0.0 };
            let p = if z.is_finite() {
                2.0 * statrs::distribution::Normal::standard().cdf(-z.abs())
            } else {
                1.0
            };
            SnpParamResult {
                lhs: row.lhs.clone(),
                op: row.op.to_string(),
                rhs: row.rhs.clone(),
                est,
                se: se_val,
                z_stat: z,
                p_value: p,
            }
        })
        .collect();

    SnpResult {
        snp_idx,
        params,
        chisq: fit.objective,
        chisq_df: model.df(),
        converged: fit.converged,
        warning: None,
    }
}
