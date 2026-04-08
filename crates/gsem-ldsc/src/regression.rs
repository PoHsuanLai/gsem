/// Result of block-partitioned weighted least squares regression.
#[derive(Debug, Clone)]
pub struct RegressionResult {
    /// Regression coefficients [slope, intercept]
    pub coef: [f64; 2],
    /// Per-block XtX accumulations (n_blocks × 2 × 2 flattened)
    pub xtx_blocks: Vec<[f64; 4]>,
    /// Per-block Xty accumulations (n_blocks × 2)
    pub xty_blocks: Vec<[f64; 2]>,
    /// Total XtX
    pub xtx_total: [f64; 4],
    /// Total Xty
    pub xty_total: [f64; 2],
}

/// Run block-partitioned weighted least squares regression.
///
/// Regresses y on x (LD scores), with a column of ones for the intercept.
/// Weights are applied to both X and y.
///
/// Returns per-block XtX/Xty for jackknife.
pub fn block_regression(
    ld: &[f64],
    y: &[f64],
    weights: &[f64],
    n_blocks: usize,
) -> RegressionResult {
    let n = ld.len();
    assert_eq!(n, y.len());
    assert_eq!(n, weights.len());

    // Compute block boundaries
    let block_bounds = compute_block_bounds(n, n_blocks);
    let actual_blocks = block_bounds.len() - 1;

    // Initialize per-block accumulators
    let mut xtx_blocks = vec![[0.0f64; 4]; actual_blocks];
    let mut xty_blocks = vec![[0.0f64; 2]; actual_blocks];

    // Accumulate per-block XtX and Xty
    for b in 0..actual_blocks {
        let start = block_bounds[b];
        let end = block_bounds[b + 1];

        for i in start..end {
            let w = weights[i];
            let x0 = ld[i] * w; // weighted LD score
            let x1 = w; // weighted intercept (1 * w)
            let yi = y[i] * w; // weighted response

            // XtX accumulation (2x2 symmetric)
            xtx_blocks[b][0] += x0 * x0; // (0,0)
            xtx_blocks[b][1] += x0 * x1; // (0,1)
            xtx_blocks[b][2] += x1 * x0; // (1,0)
            xtx_blocks[b][3] += x1 * x1; // (1,1)

            // Xty accumulation
            xty_blocks[b][0] += x0 * yi;
            xty_blocks[b][1] += x1 * yi;
        }
    }

    // Sum to get totals
    let mut xtx_total = [0.0f64; 4];
    let mut xty_total = [0.0f64; 2];
    for block in &xtx_blocks {
        for (total, val) in xtx_total.iter_mut().zip(block.iter()) {
            *total += val;
        }
    }
    for block in &xty_blocks {
        for (total, val) in xty_total.iter_mut().zip(block.iter()) {
            *total += val;
        }
    }

    // Solve: coef = solve(XtX, Xty)
    let coef = solve_2x2(&xtx_total, &xty_total);

    RegressionResult {
        coef,
        xtx_blocks,
        xty_blocks,
        xtx_total,
        xty_total,
    }
}

/// Solve a 2x2 linear system: [a b; c d] * [x; y] = [e; f]
pub(crate) fn solve_2x2(xtx: &[f64; 4], xty: &[f64; 2]) -> [f64; 2] {
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

/// Compute block boundaries for n items split into n_blocks blocks.
/// Returns Vec of length n_blocks+1 with boundary indices.
fn compute_block_bounds(n: usize, n_blocks: usize) -> Vec<usize> {
    let mut bounds = Vec::with_capacity(n_blocks + 1);
    for i in 0..=n_blocks {
        bounds.push((i * n) / n_blocks);
    }
    bounds
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_regression() {
        // y = 2*x + 1
        let ld = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![3.0, 5.0, 7.0, 9.0, 11.0];
        let w = vec![1.0; 5];
        let result = block_regression(&ld, &y, &w, 5);
        assert!(
            (result.coef[0] - 2.0).abs() < 1e-10,
            "slope should be 2.0, got {}",
            result.coef[0]
        );
        assert!(
            (result.coef[1] - 1.0).abs() < 1e-10,
            "intercept should be 1.0, got {}",
            result.coef[1]
        );
    }

    #[test]
    fn test_block_bounds() {
        let bounds = compute_block_bounds(10, 3);
        assert_eq!(bounds, vec![0, 3, 6, 10]);
    }

    #[test]
    fn test_solve_2x2() {
        // [2 1; 1 3] * [x; y] = [5; 7] => x=8/5, y=9/5
        let xtx = [2.0, 1.0, 1.0, 3.0];
        let xty = [5.0, 7.0];
        let sol = solve_2x2(&xtx, &xty);
        assert!((sol[0] - 1.6).abs() < 1e-10);
        assert!((sol[1] - 1.8).abs() < 1e-10);
    }
}
