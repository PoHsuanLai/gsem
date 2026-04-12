use faer::Mat;
use faer::linalg::solvers::DenseSolveCore;

use crate::syntax::{Op, ParTable};

/// Internal SEM model representation using LISREL all-y parameterization.
///
/// Model-implied covariance: Sigma = Lambda * (I-Beta)^{-1} * Psi * ((I-Beta)^{-1})' * Lambda' + Theta
#[derive(Debug, Clone)]
pub struct Model {
    /// Factor loading matrix (p_observed × m_latent)
    pub lambda: Mat<f64>,
    /// Latent covariance matrix (m × m)
    pub psi: Mat<f64>,
    /// Residual covariance matrix (p × p)
    pub theta: Mat<f64>,
    /// Regression between latents (m × m), 0 for CFA
    pub beta: Mat<f64>,
    /// Observed variable names (order matters)
    pub obs_names: Vec<String>,
    /// Latent variable names
    pub lat_names: Vec<String>,
    /// Map: (matrix_id, row, col) for each free parameter
    pub free_params: Vec<ParamLocation>,
    /// Lower bounds for free parameters (None = unbounded)
    pub lower_bounds: Vec<Option<f64>>,
    /// Equality-constrained cells that mirror a primary free parameter
    aliases: Vec<ParamAlias>,
}

/// Which matrix and position a free parameter belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatrixId {
    Lambda,
    Psi,
    Theta,
    Beta,
}

/// Location of a free parameter within a LISREL matrix.
#[derive(Debug, Clone, Copy)]
pub struct ParamLocation {
    /// Which LISREL matrix this parameter belongs to
    pub matrix: MatrixId,
    /// Row index within the matrix
    pub row: usize,
    /// Column index within the matrix
    pub col: usize,
}

/// An alias: a matrix cell that shares a free parameter with another cell
/// (equality constraint via label, e.g. `a*V1 + a*V2`).
#[derive(Debug, Clone, Copy)]
struct ParamAlias {
    /// Index into `free_params` that this cell mirrors
    source: usize,
    /// The matrix cell to keep in sync
    loc: ParamLocation,
}

impl Model {
    /// Build a Model from a parsed parameter table.
    pub fn from_partable(pt: &ParTable, obs_order: &[String]) -> Self {
        let mut lat_names = pt.latent_vars();
        let obs_names: Vec<String> = if obs_order.is_empty() {
            pt.observed_vars()
        } else {
            obs_order.to_vec()
        };

        // Detect observed variables used as regression predictors (e.g., F1 ~ SNP).
        // These need phantom latents so we can represent the path in the Beta matrix.
        // Phantom latent gets: Lambda[obs, phantom] = 1 (fixed), Theta[obs,obs] = 0 (fixed).
        let mut phantom_map: std::collections::HashMap<String, usize> =
            std::collections::HashMap::new();
        for row in &pt.rows {
            if row.op == Op::Regression {
                let rhs_is_lat = lat_names.iter().any(|n| n == &row.rhs);
                let rhs_is_obs = obs_names.iter().any(|n| n == &row.rhs);
                if !rhs_is_lat && rhs_is_obs && !phantom_map.contains_key(&row.rhs) {
                    let phantom_idx = lat_names.len();
                    phantom_map.insert(row.rhs.clone(), phantom_idx);
                    lat_names.push(row.rhs.clone());
                }
                // Also check lhs
                let lhs_is_lat = lat_names.iter().any(|n| n == &row.lhs);
                let lhs_is_obs = obs_names.iter().any(|n| n == &row.lhs);
                if !lhs_is_lat && lhs_is_obs && !phantom_map.contains_key(&row.lhs) {
                    let phantom_idx = lat_names.len();
                    phantom_map.insert(row.lhs.clone(), phantom_idx);
                    lat_names.push(row.lhs.clone());
                }
            }
        }

        let p = obs_names.len();
        let m = lat_names.len();

        let mut lambda = Mat::zeros(p, m);
        let mut psi = Mat::zeros(m, m);
        let mut theta = Mat::zeros(p, p);
        let mut beta = Mat::zeros(m, m);
        let mut free_params: Vec<ParamLocation> = Vec::new();
        let mut lower_bounds: Vec<Option<f64>> = Vec::new();
        let mut aliases: Vec<ParamAlias> = Vec::new();
        // label -> index into free_params (for equality constraints)
        let mut label_map: std::collections::HashMap<String, usize> =
            std::collections::HashMap::new();

        // Set up phantom latents: Lambda[obs, phantom] = 1.0 (fixed)
        for (obs_name, &phantom_lat_idx) in &phantom_map {
            if let Some(obs_idx) = obs_names.iter().position(|n| n == obs_name) {
                lambda[(obs_idx, phantom_lat_idx)] = 1.0;
            }
        }

        /// Register a free parameter, respecting equality constraints via labels.
        /// If the label has been seen before, add an alias instead of a new free
        /// parameter so both cells share the same optimizer slot.
        #[inline]
        fn add_free(
            loc: ParamLocation,
            label: &Option<String>,
            lower_bound: Option<f64>,
            free_params: &mut Vec<ParamLocation>,
            lower_bounds: &mut Vec<Option<f64>>,
            aliases: &mut Vec<ParamAlias>,
            label_map: &mut std::collections::HashMap<String, usize>,
        ) {
            if let Some(lbl) = label {
                if let Some(&existing_idx) = label_map.get(lbl) {
                    // Equality constraint: alias this cell to the existing free param
                    aliases.push(ParamAlias {
                        source: existing_idx,
                        loc,
                    });
                    return;
                }
                // First occurrence of this label
                label_map.insert(lbl.clone(), free_params.len());
            }
            free_params.push(loc);
            lower_bounds.push(lower_bound);
        }

        for row in &pt.rows {
            match row.op {
                Op::Loading => {
                    let lat_idx = lat_names.iter().position(|n| n == &row.lhs);
                    let obs_idx = obs_names.iter().position(|n| n == &row.rhs);
                    if let (Some(li), Some(oi)) = (lat_idx, obs_idx) {
                        if row.free > 0 {
                            let loc = ParamLocation {
                                matrix: MatrixId::Lambda,
                                row: oi,
                                col: li,
                            };
                            add_free(
                                loc,
                                &row.label,
                                row.lower_bound,
                                &mut free_params,
                                &mut lower_bounds,
                                &mut aliases,
                                &mut label_map,
                            );
                            // Start value
                            lambda[(oi, li)] = if row.value != 0.0 { row.value } else { 0.5 };
                        } else {
                            lambda[(oi, li)] = row.value;
                        }
                    }
                }
                Op::Covariance => {
                    let lhs_lat = lat_names.iter().position(|n| n == &row.lhs);
                    let rhs_lat = lat_names.iter().position(|n| n == &row.rhs);
                    let lhs_obs = obs_names.iter().position(|n| n == &row.lhs);
                    let rhs_obs = obs_names.iter().position(|n| n == &row.rhs);

                    if let (Some(li), Some(ri)) = (lhs_lat, rhs_lat) {
                        // Latent covariance/variance (includes phantom latents)
                        if row.free > 0 {
                            let loc = ParamLocation {
                                matrix: MatrixId::Psi,
                                row: li,
                                col: ri,
                            };
                            add_free(
                                loc,
                                &row.label,
                                row.lower_bound,
                                &mut free_params,
                                &mut lower_bounds,
                                &mut aliases,
                                &mut label_map,
                            );
                            psi[(li, ri)] = if row.value != 0.0 { row.value } else { 0.5 };
                            if li != ri {
                                psi[(ri, li)] = psi[(li, ri)];
                            }
                        } else {
                            psi[(li, ri)] = row.value;
                            if li != ri {
                                psi[(ri, li)] = row.value;
                            }
                        }
                    } else if let (Some(li), Some(ri)) = (lhs_obs, rhs_obs) {
                        // Residual covariance/variance
                        if row.free > 0 {
                            let loc = ParamLocation {
                                matrix: MatrixId::Theta,
                                row: li,
                                col: ri,
                            };
                            add_free(
                                loc,
                                &row.label,
                                row.lower_bound,
                                &mut free_params,
                                &mut lower_bounds,
                                &mut aliases,
                                &mut label_map,
                            );
                            theta[(li, ri)] = if row.value != 0.0 { row.value } else { 0.5 };
                            if li != ri {
                                theta[(ri, li)] = theta[(li, ri)];
                            }
                        } else {
                            theta[(li, ri)] = row.value;
                            if li != ri {
                                theta[(ri, li)] = row.value;
                            }
                        }
                    }
                }
                Op::Regression => {
                    // Resolve both sides: use phantom latent index if the variable is observed
                    let lhs_idx = lat_names.iter().position(|n| n == &row.lhs);
                    let rhs_idx = lat_names.iter().position(|n| n == &row.rhs);

                    if let (Some(li), Some(ri)) = (lhs_idx, rhs_idx) {
                        if row.free > 0 {
                            let loc = ParamLocation {
                                matrix: MatrixId::Beta,
                                row: li,
                                col: ri,
                            };
                            add_free(
                                loc,
                                &row.label,
                                row.lower_bound,
                                &mut free_params,
                                &mut lower_bounds,
                                &mut aliases,
                                &mut label_map,
                            );
                            beta[(li, ri)] = row.value;
                        } else {
                            beta[(li, ri)] = row.value;
                        }
                    }
                }
                Op::Defined => {
                    // Handled separately
                }
            }
        }

        Model {
            lambda,
            psi,
            theta,
            beta,
            obs_names,
            lat_names,
            free_params,
            lower_bounds,
            aliases,
        }
    }

    /// Get current free parameter values as a vector.
    pub fn get_param_vec(&self) -> Vec<f64> {
        self.free_params
            .iter()
            .map(|loc| match loc.matrix {
                MatrixId::Lambda => self.lambda[(loc.row, loc.col)],
                MatrixId::Psi => self.psi[(loc.row, loc.col)],
                MatrixId::Theta => self.theta[(loc.row, loc.col)],
                MatrixId::Beta => self.beta[(loc.row, loc.col)],
            })
            .collect()
    }

    /// Set free parameter values from a vector.
    pub fn set_param_vec(&mut self, x: &[f64]) {
        debug_assert_eq!(x.len(), self.free_params.len());
        // Copy free_params to avoid borrow conflict
        let locs: Vec<ParamLocation> = self.free_params.clone();
        for (i, loc) in locs.iter().enumerate() {
            self.write_cell(loc, x[i]);
        }
        // Propagate equality constraints: aliased cells get the same value
        // as their source free parameter.
        let aliases = self.aliases.clone();
        for alias in &aliases {
            let val = x[alias.source];
            self.write_cell(&alias.loc, val);
        }
    }

    /// Write a value to the matrix cell indicated by `loc`, handling
    /// symmetric off-diagonal entries in Psi and Theta.
    fn write_cell(&mut self, loc: &ParamLocation, val: f64) {
        match loc.matrix {
            MatrixId::Lambda => self.lambda[(loc.row, loc.col)] = val,
            MatrixId::Psi => {
                self.psi[(loc.row, loc.col)] = val;
                if loc.row != loc.col {
                    self.psi[(loc.col, loc.row)] = val;
                }
            }
            MatrixId::Theta => {
                self.theta[(loc.row, loc.col)] = val;
                if loc.row != loc.col {
                    self.theta[(loc.col, loc.row)] = val;
                }
            }
            MatrixId::Beta => self.beta[(loc.row, loc.col)] = val,
        }
    }

    /// Compute model-implied covariance matrix.
    ///
    /// Sigma = Lambda * (I - Beta)^{-1} * Psi * ((I - Beta)^{-1})' * Lambda' + Theta
    pub fn implied_cov(&self) -> Mat<f64> {
        let m = self.lat_names.len();

        // (I - Beta)^{-1}
        let i_minus_beta = Mat::from_fn(m, m, |i, j| {
            (if i == j { 1.0 } else { 0.0 }) - self.beta[(i, j)]
        });

        // For simple CFA (Beta=0), this is just identity, so skip inverse
        let has_beta = self.beta.col_iter().any(|c| c.iter().any(|&x| x != 0.0));

        if has_beta {
            // Need to invert (I - Beta)
            // Use LU decomposition
            let lu = i_minus_beta.partial_piv_lu();
            let inv = lu.inverse();
            // Sigma = Lambda * inv * Psi * inv' * Lambda' + Theta
            let lambda_inv = &self.lambda * &inv;
            let inner = &lambda_inv * &self.psi * lambda_inv.transpose();
            &inner + &self.theta
        } else {
            // Simple CFA: Sigma = Lambda * Psi * Lambda' + Theta
            let lam_psi = &self.lambda * &self.psi;
            let inner = &lam_psi * self.lambda.transpose();
            &inner + &self.theta
        }
    }

    /// Number of free parameters.
    pub fn n_free(&self) -> usize {
        self.free_params.len()
    }

    /// Returns alias (source_index, location) pairs for the Jacobian.
    pub fn alias_pairs(&self) -> Vec<(usize, ParamLocation)> {
        self.aliases.iter().map(|a| (a.source, a.loc)).collect()
    }

    /// Check for negative variance estimates (Heywood cases).
    ///
    /// Returns a list of (variable_name, variance_value) for any diagonal
    /// entries in Theta (residual) or Psi (factor) that are negative.
    pub fn negative_variances(&self) -> Vec<(String, f64)> {
        let mut bad = Vec::new();
        for i in 0..self.obs_names.len() {
            let v = self.theta[(i, i)];
            if v < 0.0 {
                bad.push((self.obs_names[i].clone(), v));
            }
        }
        for i in 0..self.lat_names.len() {
            let v = self.psi[(i, i)];
            if v < 0.0 {
                bad.push((self.lat_names[i].clone(), v));
            }
        }
        bad
    }

    /// Degrees of freedom.
    pub fn df(&self) -> usize {
        let p = self.obs_names.len();
        let n_moments = p * (p + 1) / 2;
        n_moments.saturating_sub(self.n_free())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::syntax::parse_model;

    #[test]
    fn test_model_from_cfa() {
        let pt = parse_model("F1 =~ NA*V1 + V2 + V3\nF1 ~~ 1*F1", false).unwrap();
        let model = Model::from_partable(&pt, &[]);
        assert_eq!(model.obs_names.len(), 3);
        assert_eq!(model.lat_names.len(), 1);
        assert_eq!(model.lambda.nrows(), 3);
        assert_eq!(model.lambda.ncols(), 1);
        // 3 free loadings (NA* frees first) + 3 auto-added V1/V2/V3 residual
        // variances + 0 F1 variance (fixed to 1) = 6 free params
        assert_eq!(model.n_free(), 6);
    }

    #[test]
    fn test_implied_cov_identity() {
        // 1-factor model: F =~ 1*V1 + V2, F ~~ 1*F.
        // Parser auto-adds free residual variances V1~~V1, V2~~V2 (start 0.5).
        let pt = parse_model("F1 =~ V1 + V2\nF1 ~~ 1*F1", false).unwrap();
        let model = Model::from_partable(&pt, &[]);
        let sigma = model.implied_cov();
        // Lambda = [1, 0.5]', Psi = [[1]], Theta = diag(0.5, 0.5)
        // Sigma = Lambda*Psi*Lambda' + Theta
        //       = [[1, 0.5], [0.5, 0.25]] + diag(0.5, 0.5)
        //       = [[1.5, 0.5], [0.5, 0.75]]
        assert_eq!(sigma.nrows(), 2);
        assert!((sigma[(0, 0)] - 1.5).abs() < 1e-10, "got {}", sigma[(0, 0)]);
        assert!((sigma[(0, 1)] - 0.5).abs() < 1e-10, "got {}", sigma[(0, 1)]);
        assert!(
            (sigma[(1, 1)] - 0.75).abs() < 1e-10,
            "got {}",
            sigma[(1, 1)]
        );
    }

    #[test]
    fn test_param_vec_roundtrip() {
        let pt = parse_model("F1 =~ NA*V1 + V2 + V3\nF1 ~~ 1*F1", false).unwrap();
        let mut model = Model::from_partable(&pt, &[]);
        let original = model.get_param_vec();
        let modified: Vec<f64> = original.iter().map(|x| x + 0.1).collect();
        model.set_param_vec(&modified);
        let recovered = model.get_param_vec();
        for (a, b) in modified.iter().zip(recovered.iter()) {
            assert!((a - b).abs() < 1e-15);
        }
    }

    #[test]
    fn test_equality_constraint_reduces_free_params() {
        // Without equality: F1 =~ NA*V1 + V2 → 2 free loadings
        let pt1 = parse_model("F1 =~ NA*V1 + V2\nF1 ~~ 1*F1", false).unwrap();
        let m1 = Model::from_partable(&pt1, &[]);

        // With equality: F1 =~ a*V1 + a*V2 → 1 free loading (shared)
        let pt2 = parse_model("F1 =~ a*V1 + a*V2\nF1 ~~ 1*F1", false).unwrap();
        let m2 = Model::from_partable(&pt2, &[]);

        // Both have auto-added residual variances (V1~~V1, V2~~V2)
        assert_eq!(m1.n_free(), 4, "2 free loadings + 2 residuals");
        assert_eq!(m2.n_free(), 3, "1 shared loading + 2 residuals");
    }

    #[test]
    fn test_equality_constraint_propagates_values() {
        let pt = parse_model("F1 =~ a*V1 + a*V2\nF1 ~~ 1*F1\nV1 ~~ V1\nV2 ~~ V2", false).unwrap();
        let mut model = Model::from_partable(&pt, &["V1".into(), "V2".into()]);

        // Set the shared loading to 0.7, residuals to 0.3 and 0.4
        model.set_param_vec(&[0.7, 0.3, 0.4]);

        // Both Lambda entries should be 0.7
        assert!(
            (model.lambda[(0, 0)] - 0.7).abs() < 1e-15,
            "Lambda[V1,F1] should be 0.7, got {}",
            model.lambda[(0, 0)]
        );
        assert!(
            (model.lambda[(1, 0)] - 0.7).abs() < 1e-15,
            "Lambda[V2,F1] should be 0.7 (aliased), got {}",
            model.lambda[(1, 0)]
        );
    }

    #[test]
    fn test_equality_constraint_implied_cov() {
        // F1 =~ a*V1 + a*V2, F1~~1*F1, V1~~V1, V2~~V2
        // Lambda = [a, a]', Psi = 1, Theta = diag(t1, t2)
        // Sigma = a^2 * [[1,1],[1,1]] + diag(t1, t2)
        let pt = parse_model("F1 =~ a*V1 + a*V2\nF1 ~~ 1*F1\nV1 ~~ V1\nV2 ~~ V2", false).unwrap();
        let mut model = Model::from_partable(&pt, &["V1".into(), "V2".into()]);
        model.set_param_vec(&[0.6, 0.2, 0.3]);

        let sigma = model.implied_cov();
        let a = 0.6;
        assert!((sigma[(0, 0)] - (a * a + 0.2)).abs() < 1e-10);
        assert!((sigma[(1, 1)] - (a * a + 0.3)).abs() < 1e-10);
        assert!((sigma[(0, 1)] - a * a).abs() < 1e-10);
        assert!((sigma[(1, 0)] - a * a).abs() < 1e-10);
    }

    #[test]
    fn test_negative_variances_detected() {
        let pt = parse_model("F1 =~ NA*V1 + V2\nF1 ~~ 1*F1\nV1 ~~ V1\nV2 ~~ V2", false).unwrap();
        let mut model = Model::from_partable(&pt, &["V1".into(), "V2".into()]);
        // Set V2 residual variance negative
        model.set_param_vec(&[0.5, 0.5, 0.3, -0.1]);
        let bad = model.negative_variances();
        assert_eq!(bad.len(), 1);
        assert_eq!(bad[0].0, "V2");
        assert!((bad[0].1 - (-0.1)).abs() < 1e-15);
    }
}
