//! Model-based functional enrichment — port of R GenomicSEM's `enrich`.
//!
//! Fits a SEM `model` to a baseline annotation's S/V, fixes the
//! regressions/loadings (or covariances/variances, per `fix`) to their
//! baseline estimates, then re-fits the model to each annotation's S/V
//! freeing the remaining parameters. For each target parameter it reports
//!
//! ```text
//!   enrichment    = (est_annot / est_baseline) / Prop      (null = 1)
//!   enrichment_se = (se_annot  / |est_baseline|) / Prop
//!   enrichment_p  = 1 - Phi((enrichment - 1) / enrichment_se)   (1-sided)
//! ```

use std::collections::HashMap;

use anyhow::Result;
use faer::Mat;
use statrs::distribution::{ContinuousCDF, Normal};

use crate::fit::fit_partable;
use crate::syntax::{Op, parse_model};

/// Which parameter classes are fixed to baseline (the others, plus the
/// target `params`, are freed per annotation). Mirrors R's `fix` argument.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FixMode {
    /// Fix regressions + loadings; free variances/covariances (`~~`).
    Regressions,
    /// Fix covariances; free regressions, loadings, and (residual) variances.
    Covariances,
    /// Fix (residual) variances; free regressions, loadings, and covariances.
    Variances,
}

/// Per-parameter, per-annotation enrichment results.
#[derive(Debug, Clone)]
pub struct ModelEnrichResult {
    pub annotations: Vec<String>,
    /// Target parameter keys (e.g. `"F1~~F1"`), whitespace-stripped.
    pub params: Vec<String>,
    /// `[param][annotation]` enrichment estimate (null = 1).
    pub enrichment: Vec<Vec<f64>>,
    /// `[param][annotation]` enrichment standard error.
    pub se: Vec<Vec<f64>>,
    /// `[param][annotation]` 1-sided enrichment p-value.
    pub p: Vec<Vec<f64>>,
}

fn key(lhs: &str, op: Op, rhs: &str) -> String {
    format!("{lhs}{op}{rhs}") // matches R's paste0(lhs, op, rhs), e.g. "F1~~F1"
}

/// Run model-based functional enrichment. `s_list[0]`/`v_list[0]` is the
/// baseline annotation (`prop[0]` should be 1).
#[allow(clippy::too_many_arguments)]
pub fn model_enrichment(
    s_list: &[Mat<f64>],
    v_list: &[Mat<f64>],
    prop: &[f64],
    annot_names: &[String],
    obs_names: &[String],
    model_str: &str,
    params: &[String],
    fix: FixMode,
    estimation: crate::EstimationMethod,
) -> Result<ModelEnrichResult> {
    let params_norm: Vec<String> = params.iter().map(|p| p.replace(' ', "")).collect();
    let pt0 = parse_model(model_str, false).map_err(|e| anyhow::anyhow!("{e}"))?;

    // 1. Fit the free model to the baseline annotation.
    let base = fit_partable(&pt0, obs_names, &s_list[0], &v_list[0], estimation);
    let mut base_est: HashMap<String, f64> = HashMap::new();
    for p in &base.parameters {
        base_est.insert(key(&p.lhs, p.op, &p.rhs), p.est);
    }
    for row in &pt0.rows {
        base_est
            .entry(key(&row.lhs, row.op, &row.rhs))
            .or_insert(row.value); // originally-fixed rows (e.g. marker loading = 1)
    }

    // 2. Build the partable with non-target classes fixed to baseline.
    let mut pt_fixed = pt0.clone();
    let mut free_counter = 0usize;
    for row in &mut pt_fixed.rows {
        let k = key(&row.lhs, row.op, &row.rhs);
        let is_param = params_norm.contains(&k);
        let free = match fix {
            FixMode::Regressions => is_param || row.op == Op::Covariance,
            FixMode::Covariances => {
                is_param
                    || row.op == Op::Regression
                    || row.op == Op::Loading
                    || (row.op == Op::Covariance && row.lhs == row.rhs)
            }
            FixMode::Variances => {
                is_param
                    || row.op == Op::Regression
                    || row.op == Op::Loading
                    || (row.op == Op::Covariance && row.lhs != row.rhs)
            }
        };
        if free {
            free_counter += 1;
            row.free = free_counter;
        } else {
            row.free = 0;
            row.value = *base_est.get(&k).unwrap_or(&row.value);
        }
    }

    // 3-4. Fit the fixed model per annotation and compute enrichment.
    let n_annot = s_list.len();
    let normal = Normal::new(0.0, 1.0).expect("standard normal");
    let mut enrichment = vec![vec![0.0; n_annot]; params_norm.len()];
    let mut se_out = vec![vec![0.0; n_annot]; params_norm.len()];
    let mut p_out = vec![vec![0.0; n_annot]; params_norm.len()];

    for a in 0..n_annot {
        let res_a = fit_partable(&pt_fixed, obs_names, &s_list[a], &v_list[a], estimation);
        let est_a: HashMap<String, (f64, f64)> = res_a
            .parameters
            .iter()
            .map(|p| (key(&p.lhs, p.op, &p.rhs), (p.est, p.se)))
            .collect();

        for (pi, pk) in params_norm.iter().enumerate() {
            let est_base = *base_est.get(pk).unwrap_or(&1.0);
            let (ea, sa) = *est_a.get(pk).unwrap_or(&(0.0, 0.0));
            let enr = (ea / est_base) / prop[a];
            let enr_se = (sa / est_base.abs()) / prop[a];
            enrichment[pi][a] = enr;
            se_out[pi][a] = enr_se;
            p_out[pi][a] = if enr_se > 0.0 {
                1.0 - normal.cdf((enr - 1.0) / enr_se)
            } else {
                0.5
            };
        }
    }

    Ok(ModelEnrichResult {
        annotations: annot_names.to_vec(),
        params: params_norm,
        enrichment,
        se: se_out,
        p: p_out,
    })
}
