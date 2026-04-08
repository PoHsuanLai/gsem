use faer::Mat;

use crate::regression::RegressionResult;

/// Result of delete-one-block jackknife.
#[derive(Debug, Clone)]
pub struct JackknifeResult {
    /// Jackknife standard errors for each coefficient
    pub se: [f64; 2],
    /// Pseudo-values for the slope coefficient (one per block), used for V matrix
    pub pseudo_slope: Vec<f64>,
}

/// Perform delete-one-block jackknife on a regression result.
///
/// For each block b:
///   xtx_loo = xtx_total - xtx_block[b]
///   xty_loo = xty_total - xty_block[b]
///   coef_loo = solve(xtx_loo, xty_loo)
///   pseudo[b] = n_blocks * coef_total - (n_blocks-1) * coef_loo
///
/// SE = sqrt(var(pseudo) / n_blocks)
pub fn jackknife(result: &RegressionResult) -> JackknifeResult {
    let n_blocks = result.xtx_blocks.len();
    let nb = n_blocks as f64;

    let mut pseudo_values = vec![[0.0f64; 2]; n_blocks];

    for b in 0..n_blocks {
        // Leave-one-out XtX and Xty
        let mut xtx_loo = result.xtx_total;
        let mut xty_loo = result.xty_total;
        for i in 0..4 {
            xtx_loo[i] -= result.xtx_blocks[b][i];
        }
        for i in 0..2 {
            xty_loo[i] -= result.xty_blocks[b][i];
        }

        // Solve leave-one-out regression
        let coef_loo = solve_2x2(&xtx_loo, &xty_loo);

        // Pseudo-values: n_blocks * coef_total - (n_blocks-1) * coef_loo
        for i in 0..2 {
            pseudo_values[b][i] = nb * result.coef[i] - (nb - 1.0) * coef_loo[i];
        }
    }

    // Compute variance of pseudo-values
    let mut mean = [0.0f64; 2];
    for pv in &pseudo_values {
        mean[0] += pv[0];
        mean[1] += pv[1];
    }
    mean[0] /= nb;
    mean[1] /= nb;

    let mut var = [0.0f64; 2];
    for pv in &pseudo_values {
        var[0] += (pv[0] - mean[0]).powi(2);
        var[1] += (pv[1] - mean[1]).powi(2);
    }
    var[0] /= nb * (nb - 1.0);
    var[1] /= nb * (nb - 1.0);

    let se = [var[0].sqrt(), var[1].sqrt()];
    let pseudo_slope: Vec<f64> = pseudo_values.iter().map(|pv| pv[0]).collect();

    JackknifeResult { se, pseudo_slope }
}

/// Construct the V matrix (sampling covariance of vech(S)) from jackknife pseudo-values.
///
/// V[i,j] = cov(pseudo_i, pseudo_j) / n_blocks
pub fn construct_v_matrix(all_pseudos: &[Vec<f64>], n_blocks: usize) -> Mat<f64> {
    let kstar = all_pseudos.len();
    let nb = n_blocks as f64;

    // Compute means
    let means: Vec<f64> = all_pseudos
        .iter()
        .map(|p| p.iter().sum::<f64>() / nb)
        .collect();

    // Compute covariance matrix
    let mut v = Mat::zeros(kstar, kstar);
    for i in 0..kstar {
        for j in i..kstar {
            let mut cov = 0.0;
            let n = all_pseudos[i].len().min(all_pseudos[j].len());
            for b in 0..n {
                cov += (all_pseudos[i][b] - means[i]) * (all_pseudos[j][b] - means[j]);
            }
            cov /= nb * (nb - 1.0);
            v[(i, j)] = cov;
            v[(j, i)] = cov;
        }
    }

    v
}

/// Solve a 2x2 linear system (duplicated from regression for independence).
fn solve_2x2(xtx: &[f64; 4], xty: &[f64; 2]) -> [f64; 2] {
    let a = xtx[0];
    let b = xtx[1];
    let c = xtx[2];
    let d = xtx[3];
    let det = a * d - b * c;
    if det.abs() < 1e-30 {
        return [0.0, 0.0];
    }
    let x = (d * xty[0] - b * xty[1]) / det;
    let y = (a * xty[1] - c * xty[0]) / det;
    [x, y]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_v_matrix_symmetric() {
        let pseudos = vec![
            vec![1.0, 2.0, 3.0, 4.0, 5.0],
            vec![2.0, 3.0, 4.0, 5.0, 6.0],
            vec![1.5, 2.5, 3.5, 4.5, 5.5],
        ];
        let v = construct_v_matrix(&pseudos, 5);
        assert_eq!(v.nrows(), 3);
        assert_eq!(v.ncols(), 3);
        for i in 0..3 {
            for j in 0..3 {
                assert!((v[(i, j)] - v[(j, i)]).abs() < 1e-15, "V not symmetric at ({i},{j})");
            }
        }
    }
}
