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
}

/// Which matrix and position a free parameter belongs to.
#[derive(Debug, Clone, Copy)]
pub enum MatrixId {
    Lambda,
    Psi,
    Theta,
    Beta,
}

#[derive(Debug, Clone, Copy)]
pub struct ParamLocation {
    pub matrix: MatrixId,
    pub row: usize,
    pub col: usize,
}

impl Model {
    /// Build a Model from a parsed parameter table.
    pub fn from_partable(pt: &ParTable, obs_order: &[String]) -> Self {
        let lat_names = pt.latent_vars();
        let obs_names: Vec<String> = if obs_order.is_empty() {
            pt.observed_vars()
        } else {
            obs_order.to_vec()
        };

        let p = obs_names.len();
        let m = lat_names.len();

        let mut lambda = Mat::zeros(p, m);
        let mut psi = Mat::zeros(m, m);
        let mut theta = Mat::zeros(p, p);
        let mut beta = Mat::zeros(m, m);
        let mut free_params = Vec::new();
        let mut lower_bounds = Vec::new();

        for row in &pt.rows {
            match row.op {
                Op::Loading => {
                    let lat_idx = lat_names.iter().position(|n| n == &row.lhs);
                    let obs_idx = obs_names.iter().position(|n| n == &row.rhs);
                    if let (Some(li), Some(oi)) = (lat_idx, obs_idx) {
                        if row.free > 0 {
                            free_params.push(ParamLocation {
                                matrix: MatrixId::Lambda,
                                row: oi,
                                col: li,
                            });
                            lower_bounds.push(row.lower_bound);
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
                        // Latent covariance/variance
                        if row.free > 0 {
                            free_params.push(ParamLocation {
                                matrix: MatrixId::Psi,
                                row: li,
                                col: ri,
                            });
                            lower_bounds.push(row.lower_bound);
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
                            free_params.push(ParamLocation {
                                matrix: MatrixId::Theta,
                                row: li,
                                col: ri,
                            });
                            lower_bounds.push(row.lower_bound);
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
                    let lhs_lat = lat_names.iter().position(|n| n == &row.lhs);
                    let rhs_lat = lat_names.iter().position(|n| n == &row.rhs);

                    if let (Some(li), Some(ri)) = (lhs_lat, rhs_lat) {
                        if row.free > 0 {
                            free_params.push(ParamLocation {
                                matrix: MatrixId::Beta,
                                row: li,
                                col: ri,
                            });
                            lower_bounds.push(row.lower_bound);
                            beta[(li, ri)] = 0.1;
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
        assert_eq!(x.len(), self.free_params.len());
        for (i, loc) in self.free_params.iter().enumerate() {
            match loc.matrix {
                MatrixId::Lambda => self.lambda[(loc.row, loc.col)] = x[i],
                MatrixId::Psi => {
                    self.psi[(loc.row, loc.col)] = x[i];
                    if loc.row != loc.col {
                        self.psi[(loc.col, loc.row)] = x[i];
                    }
                }
                MatrixId::Theta => {
                    self.theta[(loc.row, loc.col)] = x[i];
                    if loc.row != loc.col {
                        self.theta[(loc.col, loc.row)] = x[i];
                    }
                }
                MatrixId::Beta => self.beta[(loc.row, loc.col)] = x[i],
            }
        }
    }

    /// Compute model-implied covariance matrix.
    ///
    /// Sigma = Lambda * (I - Beta)^{-1} * Psi * ((I - Beta)^{-1})' * Lambda' + Theta
    pub fn implied_cov(&self) -> Mat<f64> {
        let m = self.lat_names.len();

        // (I - Beta)^{-1}
        let i_minus_beta = {
            let mut ib = Mat::<f64>::identity(m, m);
            for i in 0..m {
                for j in 0..m {
                    ib[(i, j)] -= self.beta[(i, j)];
                }
            }
            ib
        };

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
        // 3 free loadings (NA* frees first) + 0 psi (fixed to 1) = 3 free params
        assert_eq!(model.n_free(), 3);
    }

    #[test]
    fn test_implied_cov_identity() {
        // 1-factor model: F =~ 1*V1 + 1*V2, F ~~ 1*F, no theta
        let pt = parse_model("F1 =~ V1 + V2\nF1 ~~ 1*F1", false).unwrap();
        let model = Model::from_partable(&pt, &[]);
        let sigma = model.implied_cov();
        // Lambda = [1, free]', Psi = [1], Theta = 0
        // With first loading fixed to 1 and second free (starts at 0.5):
        // Sigma[0,0] = 1*1*1 + 0 = 1
        // Sigma[0,1] = 1*1*0.5 = 0.5
        assert_eq!(sigma.nrows(), 2);
        assert!((sigma[(0, 0)] - 1.0).abs() < 1e-10);
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
}
