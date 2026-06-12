use faer::Mat;
use faer::linalg::solvers::DenseSolveCore;

/// Result of Generalized Least Squares regression.
#[derive(Debug, Clone)]
pub struct GlsResult {
    /// Coefficient estimates
    pub beta: Vec<f64>,
    /// Standard errors
    pub se: Vec<f64>,
    /// Z-statistics
    pub z: Vec<f64>,
    /// P-values
    pub p: Vec<f64>,
}

/// Run GLS regression on genetic parameters.
///
/// Port of GenomicSEM's `summaryGLS()`.
///
/// beta_gls = (X' V^{-1} X)^{-1} X' V^{-1} y
/// SE = sqrt(diag((X' V^{-1} X)^{-1}))
pub fn summary_gls(x: &Mat<f64>, y: &[f64], v: &Mat<f64>) -> Option<GlsResult> {
    let n = x.nrows();
    let p = x.ncols();
    if n != y.len() || n != v.nrows() {
        return None;
    }

    // V^{-1}
    let v_inv = v.partial_piv_lu().inverse();

    // X' V^{-1}
    let xt_vinv = x.transpose() * &v_inv;

    // X' V^{-1} X
    let xt_vinv_x = &xt_vinv * x;

    // (X' V^{-1} X)^{-1}
    let bread = xt_vinv_x.partial_piv_lu().inverse();

    // X' V^{-1} y
    let y_mat = Mat::from_fn(n, 1, |i, _| y[i]);
    let xt_vinv_y = &xt_vinv * &y_mat;

    // beta = bread * X' V^{-1} y
    let beta_mat = &bread * &xt_vinv_y;

    let beta: Vec<f64> = (0..p).map(|i| beta_mat[(i, 0)]).collect();
    let se: Vec<f64> = (0..p)
        .map(|i| {
            let v = bread[(i, i)];
            if v > 0.0 { v.sqrt() } else { 0.0 }
        })
        .collect();
    let z: Vec<f64> = beta
        .iter()
        .zip(se.iter())
        .map(|(b, s)| if *s > 0.0 { b / s } else { 0.0 })
        .collect();
    let p_vals: Vec<f64> = z
        .iter()
        .map(|&zi| {
            use statrs::distribution::{ContinuousCDF, Normal};
            let n = Normal::standard();
            2.0 * n.cdf(-zi.abs())
        })
        .collect();

    Some(GlsResult {
        beta,
        se,
        z,
        p: p_vals,
    })
}

/// Confidence-band data from [`summary_gls_bands`] — the numeric core of R
/// GenomicSEM's `summaryGLSbands` (the ggplot rendering is intentionally not
/// ported; this returns the data a caller would plot or report).
#[derive(Debug, Clone)]
pub struct GlsBands {
    /// The full GLS fit on the (intercept +) predictor (+ control) design.
    pub fit: GlsResult,
    /// Evenly spaced predictor grid (length `intervals`), from min to max of
    /// the predictor values.
    pub grid: Vec<f64>,
    /// Fitted regression line evaluated at `grid`.
    pub line: Vec<f64>,
    /// Upper band: `line[i] + band_size * se[i]`.
    pub upper: Vec<f64>,
    /// Lower band: `line[i] - band_size * se[i]`.
    pub lower: Vec<f64>,
    /// Per-grid-point SE of the fitted value (the re-centered-intercept SE).
    pub band_se: Vec<f64>,
}

/// Compute GLS regression confidence-band data (port of `summaryGLSbands`'s
/// numeric core).
///
/// Fits `y ~ predictors` by GLS with response covariance `v_y`, then evaluates
/// the fitted line and a ± `band_size`·SE envelope at `intervals` evenly spaced
/// predictor values. The per-point SE is R's approach: re-center the predictor
/// to each grid value and take the SE of the GLS intercept there.
///
/// * `predictors` — one predictor value per response element (length n).
/// * `y`          — response vector (length n).
/// * `v_y`        — n×n response covariance.
/// * `intercept`  — include an intercept column (R's `INTERCEPT`, default true).
/// * `quad`       — add a quadratic predictor² term (R's `QUAD`).
/// * `controls`   — optional control columns (each length n), R's `CONTROLVARS`.
/// * `intervals`  — number of grid points / bands (R's `INTERVALS`, default 20).
/// * `band_size`  — band half-width in SE units (R's `BAND_SIZE`, default 1).
#[allow(clippy::too_many_arguments)]
pub fn summary_gls_bands(
    predictors: &[f64],
    y: &[f64],
    v_y: &Mat<f64>,
    intercept: bool,
    quad: bool,
    controls: &[Vec<f64>],
    intervals: usize,
    band_size: f64,
) -> Option<GlsBands> {
    let n = predictors.len();
    if n == 0 || y.len() != n || v_y.nrows() != n || intervals == 0 {
        return None;
    }
    if controls.iter().any(|c| c.len() != n) {
        return None;
    }

    // Build the main design: [intercept?, predictor, predictor^2?, controls...].
    let build_design = |pred_cols: &[Vec<f64>]| -> Mat<f64> {
        let mut cols: Vec<Vec<f64>> = Vec::new();
        if intercept {
            cols.push(vec![1.0; n]);
        }
        for c in pred_cols {
            cols.push(c.clone());
        }
        for c in controls {
            cols.push(c.clone());
        }
        let p = cols.len();
        Mat::from_fn(n, p, |i, j| cols[j][i])
    };

    let mut main_cols: Vec<Vec<f64>> = vec![predictors.to_vec()];
    if quad {
        main_cols.push(predictors.iter().map(|x| x * x).collect());
    }
    let x = build_design(&main_cols);
    let fit = summary_gls(&x, y, v_y)?;

    // Predictor range for the grid and the re-centering offsets.
    let pmin = predictors.iter().cloned().fold(f64::INFINITY, f64::min);
    let pmax = predictors.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let span = pmax - pmin;

    // Per-interval band SE: re-center the predictor to each grid origin and
    // take the SE of the GLS intercept of the re-centered fit.
    let mut band_se = Vec::with_capacity(intervals);
    for i in 0..intervals {
        let origin = pmin + (i as f64) * span / intervals as f64;
        let centered: Vec<f64> = predictors.iter().map(|x| x - origin).collect();
        let mut cols: Vec<Vec<f64>> = vec![centered.clone()];
        if quad {
            cols.push(centered.iter().map(|x| x * x).collect());
        }
        let xx = build_design(&cols);
        let f = summary_gls(&xx, y, v_y)?;
        // R reads SE[1] — the intercept SE when an intercept is present,
        // otherwise the first coefficient's SE.
        band_se.push(f.se[0]);
    }

    // Evenly spaced grid and the fitted line on it.
    let grid: Vec<f64> = (0..intervals)
        .map(|i| {
            if intervals == 1 {
                pmin
            } else {
                pmin + (i as f64) * span / (intervals as f64 - 1.0)
            }
        })
        .collect();

    // Fitted line: b0 + b1*g (+ b2*g^2). With an intercept the coefficients are
    // [b0, b1, (b2)], otherwise [b1, (b2)] and there is no constant term.
    let (b0, b1, b2) = if intercept {
        let b0 = fit.beta[0];
        let b1 = fit.beta.get(1).copied().unwrap_or(0.0);
        let b2 = if quad {
            fit.beta.get(2).copied().unwrap_or(0.0)
        } else {
            0.0
        };
        (b0, b1, b2)
    } else {
        let b1 = fit.beta[0];
        let b2 = if quad {
            fit.beta.get(1).copied().unwrap_or(0.0)
        } else {
            0.0
        };
        (0.0, b1, b2)
    };
    let line: Vec<f64> = grid.iter().map(|&g| b0 + b1 * g + b2 * g * g).collect();
    let upper: Vec<f64> = line
        .iter()
        .zip(band_se.iter())
        .map(|(&l, &s)| l + band_size * s)
        .collect();
    let lower: Vec<f64> = line
        .iter()
        .zip(band_se.iter())
        .map(|(&l, &s)| l - band_size * s)
        .collect();

    Some(GlsBands {
        fit,
        grid,
        line,
        upper,
        lower,
        band_se,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gls_identity_weights() {
        // With V = I, GLS = OLS
        let x = faer::mat![[1.0, 1.0], [1.0, 2.0], [1.0, 3.0], [1.0, 4.0],];
        let y = vec![2.0, 4.0, 6.0, 8.0]; // y = 2x
        let v = Mat::<f64>::identity(4, 4);
        let result = summary_gls(&x, &y, &v).unwrap();
        assert_eq!(result.beta.len(), 2);
        // Slope should be ~2
        assert!((result.beta[1] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_gls_bands_shapes_and_line() {
        // Linear data y = 0.5 + 0.3x with V = small*I: the fitted line passes
        // through the band centres and every output vector has `intervals` len.
        let pred: Vec<f64> = (0..8).map(|i| -1.5 + 3.0 * (i as f64) / 7.0).collect();
        let y: Vec<f64> = pred.iter().map(|p| 0.5 + 0.3 * p).collect();
        let v = Mat::from_fn(8, 8, |i, j| if i == j { 0.01 } else { 0.0 });
        let b = summary_gls_bands(&pred, &y, &v, true, false, &[], 5, 1.0).unwrap();
        assert_eq!(b.grid.len(), 5);
        assert_eq!(b.line.len(), 5);
        assert_eq!(b.band_se.len(), 5);
        // line == 0.5 + 0.3*grid
        for (g, l) in b.grid.iter().zip(b.line.iter()) {
            assert!((l - (0.5 + 0.3 * g)).abs() < 1e-9);
        }
        // upper/lower straddle the line by band_size*se.
        for i in 0..5 {
            assert!((b.upper[i] - (b.line[i] + b.band_se[i])).abs() < 1e-12);
            assert!((b.lower[i] - (b.line[i] - b.band_se[i])).abs() < 1e-12);
        }
    }

    #[test]
    fn test_gls_bands_quad_branch() {
        // Quadratic data exercises the QUAD predictor^2 design path.
        let pred: Vec<f64> = (0..9).map(|i| -2.0 + 4.0 * (i as f64) / 8.0).collect();
        let y: Vec<f64> = pred.iter().map(|p| 0.2 + 0.1 * p + 0.5 * p * p).collect();
        let v = Mat::from_fn(9, 9, |i, j| if i == j { 0.01 } else { 0.0 });
        let b = summary_gls_bands(&pred, &y, &v, true, true, &[], 6, 1.0).unwrap();
        assert_eq!(b.fit.beta.len(), 3, "intercept + linear + quadratic");
        // Fitted quadratic recovers the generating coefficients.
        assert!((b.fit.beta[0] - 0.2).abs() < 1e-6);
        assert!((b.fit.beta[1] - 0.1).abs() < 1e-6);
        assert!((b.fit.beta[2] - 0.5).abs() < 1e-6);
        assert_eq!(b.line.len(), 6);
    }
}
