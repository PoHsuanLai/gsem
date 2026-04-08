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

    let pseudo_values: Vec<[f64; 2]> = result
        .xtx_blocks
        .iter()
        .zip(result.xty_blocks.iter())
        .map(|(xtx_block, xty_block)| {
            // Leave-one-out
            let mut xtx_loo = result.xtx_total;
            let mut xty_loo = result.xty_total;
            for (loo, blk) in xtx_loo.iter_mut().zip(xtx_block.iter()) {
                *loo -= blk;
            }
            for (loo, blk) in xty_loo.iter_mut().zip(xty_block.iter()) {
                *loo -= blk;
            }

            let coef_loo = crate::regression::solve_2x2(&xtx_loo, &xty_loo);
            [
                nb * result.coef[0] - (nb - 1.0) * coef_loo[0],
                nb * result.coef[1] - (nb - 1.0) * coef_loo[1],
            ]
        })
        .collect();

    // Compute mean and variance of pseudo-values
    let mean: [f64; 2] = [
        pseudo_values.iter().map(|pv| pv[0]).sum::<f64>() / nb,
        pseudo_values.iter().map(|pv| pv[1]).sum::<f64>() / nb,
    ];

    let denom = nb * (nb - 1.0);
    let se = [
        (pseudo_values
            .iter()
            .map(|pv| (pv[0] - mean[0]).powi(2))
            .sum::<f64>()
            / denom)
            .sqrt(),
        (pseudo_values
            .iter()
            .map(|pv| (pv[1] - mean[1]).powi(2))
            .sum::<f64>()
            / denom)
            .sqrt(),
    ];

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
    let denom = nb * (nb - 1.0);
    let mut v = Mat::zeros(kstar, kstar);
    for i in 0..kstar {
        for j in i..kstar {
            let cov: f64 = all_pseudos[i]
                .iter()
                .zip(all_pseudos[j].iter())
                .map(|(&a, &b)| (a - means[i]) * (b - means[j]))
                .sum::<f64>()
                / denom;
            v[(i, j)] = cov;
            v[(j, i)] = cov;
        }
    }
    v
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
                assert!(
                    (v[(i, j)] - v[(j, i)]).abs() < 1e-15,
                    "V not symmetric at ({i},{j})"
                );
            }
        }
    }
}
