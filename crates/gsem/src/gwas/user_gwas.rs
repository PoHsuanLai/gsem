use faer::Mat;
use gsem_matrix::smooth;
use gsem_sem::estimator;
use gsem_sem::model::Model;
use gsem_sem::sandwich;
use gsem_sem::syntax::{self, ParTable};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use statrs::distribution::ContinuousCDF;

use super::add_snps;
use super::gc_correction::GcMode;

/// Result for a single SNP from userGWAS.
#[derive(Debug, Clone)]
pub struct SnpResult {
    /// Index of this SNP in the input arrays
    pub snp_idx: usize,
    /// Parameter estimates for this SNP
    pub params: Vec<SnpParamResult>,
    /// Model chi-square statistic
    pub chisq: f64,
    /// Model degrees of freedom
    pub chisq_df: usize,
    /// Whether the optimizer converged
    pub converged: bool,
    /// Warning message, if any
    pub warning: Option<String>,
    /// Q_SNP heterogeneity statistic (if computed).
    pub q_snp: Option<f64>,
    /// Degrees of freedom for Q_SNP test.
    pub q_snp_df: Option<usize>,
    /// P-value for Q_SNP test.
    pub q_snp_p: Option<f64>,
}

/// A single parameter result for one SNP.
#[derive(Debug, Clone)]
pub struct SnpParamResult {
    /// Left-hand side variable name
    pub lhs: String,
    /// Operator (e.g., "=~", "~~", "~")
    pub op: String,
    /// Right-hand side variable name
    pub rhs: String,
    /// Point estimate
    pub est: f64,
    /// Standard error (sandwich-corrected)
    pub se: f64,
    /// Z-statistic (est / se)
    pub z_stat: f64,
    /// P-value (two-tailed)
    pub p_value: f64,
}

/// Configuration for userGWAS.
#[derive(Debug, Clone)]
pub struct UserGwasConfig {
    /// Lavaan-style model syntax string
    pub model: String,
    /// Estimation method: "DWLS" or "ML"
    pub estimation: String,
    /// Genomic control correction mode
    pub gc: GcMode,
    /// Maximum optimizer iterations per SNP
    pub max_iter: usize,
    /// Standardize latent variables (free all loadings, fix factor variances to 1)
    pub std_lv: bool,
    /// Log warnings when covariance matrix requires smoothing.
    pub smooth_check: bool,
    /// Override for SNP SE (default: 0.0005, treating MAF as fixed).
    /// Matches R's `SNPSE` parameter.
    pub snp_se: Option<f64>,
    /// Label for the first observed variable in the model.
    /// "SNP" for standard GWAS, "Gene" for TWAS mode.
    pub snp_label: String,
    /// Compute Q_SNP heterogeneity statistic.
    pub q_snp: bool,
    /// Fix measurement parameters at baseline estimates (fit without SNP first).
    /// Dramatically speeds up per-SNP fitting.
    pub fix_measurement: bool,
}

impl Default for UserGwasConfig {
    fn default() -> Self {
        Self {
            model: String::new(),
            estimation: "DWLS".to_string(),
            gc: GcMode::Standard,
            max_iter: 500,
            std_lv: false,
            smooth_check: false,
            snp_se: None,
            snp_label: "SNP".to_string(),
            q_snp: false,
            fix_measurement: false,
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
    let mut pt = match syntax::parse_model(&config.model, config.std_lv) {
        Ok(pt) => pt,
        Err(e) => {
            return vec![SnpResult {
                snp_idx: 0,
                params: vec![],
                chisq: f64::NAN,
                chisq_df: 0,
                converged: false,
                warning: Some(format!("model parse error: {e}")),
                q_snp: None,
                q_snp_df: None,
                q_snp_p: None,
            }];
        }
    };

    // fix_measurement: fit baseline model on S_LD (without SNP),
    // then fix all non-SNP parameters at baseline estimates
    if config.fix_measurement {
        // Build a baseline model using only the trait-trait portion of the model
        // (parameters that don't involve the SNP label)
        let snp_label = &config.snp_label;
        let obs_names: Vec<String> = (0..k).map(|i| format!("V{}", i + 1)).collect();
        // Extract only non-SNP rows from the partable
        let baseline_rows: Vec<_> = pt
            .rows
            .iter()
            .filter(|r| r.lhs != *snp_label && r.rhs != *snp_label)
            .cloned()
            .collect();
        if !baseline_rows.is_empty() {
            let baseline_pt = ParTable { rows: baseline_rows };
            let mut baseline_model = Model::from_partable(&baseline_pt, &obs_names);
            let kstar_ld = k * (k + 1) / 2;
            let v_diag: Vec<f64> = (0..kstar_ld).map(|i| v_ld[(i, i)]).collect();
            let baseline_fit = if config.estimation.to_uppercase() == "ML" {
                estimator::fit_ml(&mut baseline_model, s_ld, config.max_iter, None)
            } else {
                estimator::fit_dwls(&mut baseline_model, s_ld, &v_diag, config.max_iter, None)
            };
            if baseline_fit.converged {
                // Fix baseline parameters in the full model
                let mut free_idx = 0;
                for row in &baseline_pt.rows {
                    if row.free > 0 {
                        if let Some(&est) = baseline_fit.params.get(free_idx) {
                            // Find matching row in full pt and fix it
                            for full_row in &mut pt.rows {
                                if full_row.lhs == row.lhs
                                    && full_row.op == row.op
                                    && full_row.rhs == row.rhs
                                    && full_row.free > 0
                                {
                                    full_row.free = 0;
                                    full_row.value = est;
                                }
                            }
                        }
                        free_idx += 1;
                    }
                }
                log::info!(
                    "fix_measurement: fixed {} baseline params, {} remain free",
                    free_idx,
                    pt.rows.iter().filter(|r| r.free > 0).count()
                );
            } else {
                log::warn!("fix_measurement: baseline model did not converge, using unfixed params");
            }
        }
    }

    // Process SNPs in parallel
    let pb = ProgressBar::new(n_snps as u64);
    pb.set_style(
        ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} SNPs ({eta})",
        )
        .unwrap(),
    );

    let results = (0..n_snps)
        .into_par_iter()
        .map(|i| {
            let result = process_single_snp(
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
            );
            pb.inc(1);
            result
        })
        .collect();

    pb.finish_with_message("complete");
    results
}

/// Process a single SNP.
#[allow(clippy::too_many_arguments)]
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
    let s_smoothed = smooth::smooth_if_needed(&mut s_full);
    let v_smoothed = smooth::smooth_if_needed(&mut v_full);
    if config.smooth_check && (s_smoothed || v_smoothed) {
        log::warn!(
            "SNP {}: covariance matrix required smoothing (S={}, V={})",
            snp_idx,
            s_smoothed,
            v_smoothed,
        );
    }

    // Build model
    let mut obs_names = vec![config.snp_label.clone()];
    obs_names.extend((0..k).map(|i| format!("V{}", i + 1)));

    let mut model = Model::from_partable(pt, &obs_names);

    // Construct diagonal weight matrix
    let kstar_full = (k + 1) * (k + 2) / 2;
    let v_diag: Vec<f64> = (0..kstar_full).map(|i| v_full[(i, i)]).collect();

    // Fit model
    let fit = if config.estimation.eq_ignore_ascii_case("ML") {
        estimator::fit_ml(&mut model, &s_full, config.max_iter, None)
    } else {
        estimator::fit_dwls(&mut model, &s_full, &v_diag, config.max_iter, None)
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
            let se_val = se.get(i).copied().unwrap_or(f64::NAN);
            let z = est / se_val;
            let p = if z.is_finite() {
                2.0 * statrs::distribution::Normal::standard().cdf(-z.abs())
            } else {
                f64::NAN
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

    // Compute Q_SNP if requested
    let (q_snp_val, q_snp_df_val, q_snp_p_val) = if config.q_snp {
        let sigma_hat = model.implied_cov();
        let (q, df, p) = super::q_snp::compute_q_snp(&s_full, &sigma_hat, &v_full)
            .expect("q_snp: matrices must be square");
        (Some(q), Some(df), Some(p))
    } else {
        (None, None, None)
    };

    SnpResult {
        snp_idx,
        params,
        chisq: fit.objective,
        chisq_df: model.df(),
        converged: fit.converged,
        warning: None,
        q_snp: q_snp_val,
        q_snp_df: q_snp_df_val,
        q_snp_p: q_snp_p_val,
    }
}
