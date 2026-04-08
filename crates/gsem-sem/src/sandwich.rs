use faer::Mat;
use faer::linalg::solvers::DenseSolveCore;

use crate::delta;
use crate::model::Model;

/// Compute sandwich-corrected standard errors.
///
/// Algorithm (from commonfactor.R):
///   bread = (Delta' * W * Delta)^{-1}
///   lettuce = W * Delta
///   Ohtt = bread * lettuce' * V * lettuce * bread
///   SE = sqrt(diag(Ohtt))
///
/// Returns (standard_errors, Ohtt_matrix).
pub fn sandwich_se(model: &mut Model, w: &Mat<f64>, v: &Mat<f64>) -> (Vec<f64>, Mat<f64>) {
    let eps = 1e-7;
    let delta_mat = delta::compute_delta(model, eps);
    let n_free = model.n_free();

    // bread = (Delta' * W * Delta)^{-1}
    let dt_w = delta_mat.transpose() * w;
    let dt_w_d = &dt_w * &delta_mat;

    let bread = match dt_w_d.partial_piv_lu().inverse() {
        inv => inv,
    };

    // lettuce = W * Delta
    let lettuce = w * &delta_mat;

    // Ohtt = bread * lettuce' * V * lettuce * bread
    let lt_v = lettuce.transpose() * v;
    let lt_v_l = &lt_v * &lettuce;
    let ohtt = &bread * &lt_v_l * &bread;

    // SE = sqrt(diag(Ohtt))
    let se: Vec<f64> = (0..n_free)
        .map(|i| {
            let val = ohtt[(i, i)];
            if val > 0.0 { val.sqrt() } else { 0.0 }
        })
        .collect();

    (se, ohtt)
}
