//! Analytical Jacobian: d(vech(Sigma)) / d(theta).
//!
//! Computes the Delta matrix analytically rather than via finite differences.
//! For a model with p observed variables, m latent variables, and n_free free
//! parameters, the Jacobian is a kstar × n_free matrix where kstar = p*(p+1)/2.

use faer::Mat;
use faer::linalg::solvers::DenseSolveCore;

use crate::model::{MatrixId, Model};

/// Compute the analytical Jacobian matrix Delta[kstar × n_free].
///
/// Delta[k, j] = d(vech(Sigma))_k / d(theta_j)
///
/// where vech iterates column-major lower triangle: for j in 0..p, for i in j..p.
pub fn analytical_jacobian(model: &Model) -> Mat<f64> {
    let p = model.obs_names.len();
    let m = model.lat_names.len();
    let n_free = model.free_params.len();
    let kstar = p * (p + 1) / 2;

    let mut delta = Mat::zeros(kstar, n_free);

    // Use the structural-regression path if the model has ANY Beta parameter
    // (free or fixed non-zero). A free regression slope starts at 0.0
    // (matching lavaan) but still needs the Beta Jacobian path — we can't
    // rely on checking matrix entries alone.
    let has_beta = model.free_params.iter().any(|p| p.matrix == MatrixId::Beta)
        || model.beta.col_iter().any(|c| c.iter().any(|&x| x != 0.0));

    if has_beta {
        fill_delta_with_beta(model, &mut delta, p, m);
    } else {
        fill_delta_cfa(model, &mut delta, p);
    }

    delta
}

/// Build vech index mapping: for a p×p symmetric matrix, returns the vech index
/// for element (i, j) where i >= j. Returns None if i < j.
#[inline]
fn vech_index(p: usize, i: usize, j: usize) -> usize {
    // vech iterates: for col j in 0..p, for row i in j..p
    // Index = sum_{c=0}^{j-1} (p - c) + (i - j)
    //       = j*p - j*(j-1)/2 + (i - j)
    debug_assert!(i >= j && i < p && j < p);
    j * p - j * (j.wrapping_sub(1)) / 2 + (i - j)
}

/// CFA case: Beta = 0, Sigma = Lambda * Psi * Lambda' + Theta
fn fill_delta_cfa(model: &Model, delta: &mut Mat<f64>, p: usize) {
    // Pre-compute PsiLt = Psi * Lambda' (m × p)
    let psi_lt = &model.psi * model.lambda.transpose();

    for (param_j, loc) in model.free_params.iter().enumerate() {
        match loc.matrix {
            MatrixId::Lambda => {
                let r = loc.row; // observed index
                let c = loc.col; // latent index
                // dSigma_{i,j}/dLambda_{r,c} = delta_{i,r} * (PsiLt)_{c,j} + (PsiLt)_{c,i} * delta_{j,r}
                // When j == r and i >= r:
                for i in r..p {
                    let k = vech_index(p, i, r);
                    delta[(k, param_j)] += psi_lt[(c, i)];
                }
                // When i == r and j <= r:
                for j in 0..=r {
                    let k = vech_index(p, r, j);
                    delta[(k, param_j)] += psi_lt[(c, j)];
                }
                // Note: (r, r) was counted twice, which is correct:
                // dSigma_{r,r} = 2 * PsiLt_{c,r}
            }
            MatrixId::Psi => {
                let r = loc.row; // latent index
                let c = loc.col; // latent index
                if r == c {
                    for col_j in 0..p {
                        let lj_r = model.lambda[(col_j, r)];
                        if lj_r == 0.0 {
                            continue;
                        }
                        for row_i in col_j..p {
                            let li_r = model.lambda[(row_i, r)];
                            if li_r == 0.0 {
                                continue;
                            }
                            let k = vech_index(p, row_i, col_j);
                            delta[(k, param_j)] += li_r * lj_r;
                        }
                    }
                } else {
                    for col_j in 0..p {
                        let lj_r = model.lambda[(col_j, r)];
                        let lj_c = model.lambda[(col_j, c)];
                        if lj_r == 0.0 && lj_c == 0.0 {
                            continue;
                        }
                        for row_i in col_j..p {
                            let li_r = model.lambda[(row_i, r)];
                            let li_c = model.lambda[(row_i, c)];
                            let val = li_r * lj_c + li_c * lj_r;
                            if val != 0.0 {
                                let k = vech_index(p, row_i, col_j);
                                delta[(k, param_j)] += val;
                            }
                        }
                    }
                }
            }
            MatrixId::Theta => {
                fill_theta_derivative(delta, param_j, loc.row, loc.col, p);
            }
            MatrixId::Beta => {
                unreachable!("Beta parameter in CFA model (has_beta=false)");
            }
        }
    }
}

/// Structural regression case: Beta != 0
/// Sigma = LA * Psi * (LA)' + Theta where LA = Lambda * A, A = (I-Beta)^{-1}
fn fill_delta_with_beta(model: &Model, delta: &mut Mat<f64>, p: usize, m: usize) {
    // Compute A = (I - Beta)^{-1}
    let i_minus_beta = Mat::from_fn(m, m, |i, j| {
        (if i == j { 1.0 } else { 0.0 }) - model.beta[(i, j)]
    });
    let a_mat = i_minus_beta.partial_piv_lu().inverse();

    // LA = Lambda * A (p × m)
    let la = &model.lambda * &a_mat;

    // PsiAt = Psi * A' (m × m)
    let psi_at = &model.psi * a_mat.transpose();

    // LAPsiAt = LA * PsiAt (p × m)
    let la_psi_at = &la * &psi_at;

    // APsiAtLt = A * PsiAt * Lambda' (m × p)
    let a_psi_at_lt = &a_mat * &psi_at * model.lambda.transpose();

    // For Beta derivatives:
    // M = A * Psi * A' (m × m)
    let m_mat = &a_mat * &psi_at;
    // MLt = M * Lambda' (m × p)
    let m_lt = &m_mat * model.lambda.transpose();
    // LM = Lambda * M (p × m)
    let lm = &model.lambda * &m_mat;
    // AtLt = A' * Lambda' (m × p)
    let at_lt = a_mat.transpose() * model.lambda.transpose();

    for (param_j, loc) in model.free_params.iter().enumerate() {
        match loc.matrix {
            MatrixId::Lambda => {
                let r = loc.row;
                let c = loc.col;
                for i in r..p {
                    let k = vech_index(p, i, r);
                    delta[(k, param_j)] += la_psi_at[(i, c)];
                }
                for j in 0..=r {
                    let k = vech_index(p, r, j);
                    delta[(k, param_j)] += a_psi_at_lt[(c, j)];
                }
            }
            MatrixId::Psi => {
                let r = loc.row;
                let c = loc.col;
                if r == c {
                    for col_j in 0..p {
                        let laj_r = la[(col_j, r)];
                        if laj_r == 0.0 {
                            continue;
                        }
                        for row_i in col_j..p {
                            let lai_r = la[(row_i, r)];
                            if lai_r == 0.0 {
                                continue;
                            }
                            let k = vech_index(p, row_i, col_j);
                            delta[(k, param_j)] += lai_r * laj_r;
                        }
                    }
                } else {
                    for col_j in 0..p {
                        let laj_r = la[(col_j, r)];
                        let laj_c = la[(col_j, c)];
                        if laj_r == 0.0 && laj_c == 0.0 {
                            continue;
                        }
                        for row_i in col_j..p {
                            let lai_r = la[(row_i, r)];
                            let lai_c = la[(row_i, c)];
                            let val = lai_r * laj_c + lai_c * laj_r;
                            if val != 0.0 {
                                let k = vech_index(p, row_i, col_j);
                                delta[(k, param_j)] += val;
                            }
                        }
                    }
                }
            }
            MatrixId::Beta => {
                let r = loc.row;
                let c = loc.col;
                for col_j in 0..p {
                    let mlt_c_j = m_lt[(c, col_j)];
                    let atlt_r_j = at_lt[(r, col_j)];
                    for row_i in col_j..p {
                        let val = la[(row_i, r)] * mlt_c_j + lm[(row_i, c)] * atlt_r_j;
                        if val != 0.0 {
                            let k = vech_index(p, row_i, col_j);
                            delta[(k, param_j)] += val;
                        }
                    }
                }
            }
            MatrixId::Theta => {
                fill_theta_derivative(delta, param_j, loc.row, loc.col, p);
            }
        }
    }
}

/// Fill Theta derivative entries in delta.
fn fill_theta_derivative(delta: &mut Mat<f64>, param_j: usize, r: usize, c: usize, p: usize) {
    if r == c {
        let k = vech_index(p, r, r);
        delta[(k, param_j)] = 1.0;
    } else {
        let (big, small) = if r > c { (r, c) } else { (c, r) };
        let k = vech_index(p, big, small);
        delta[(k, param_j)] = 1.0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Model;
    use crate::syntax::parse_model;
    use gsem_matrix::vech;

    /// Compare analytical Jacobian against numerical finite differences.
    fn numerical_jacobian(model: &mut Model, eps: f64) -> Mat<f64> {
        let params = model.get_param_vec();
        let n_free = params.len();
        let p = model.obs_names.len();
        let kstar = p * (p + 1) / 2;
        let mut delta = Mat::zeros(kstar, n_free);

        for j in 0..n_free {
            let mut p_plus = params.clone();
            let mut p_minus = params.clone();
            p_plus[j] += eps;
            p_minus[j] -= eps;

            model.set_param_vec(&p_plus);
            let sigma_plus = vech::vech(&model.implied_cov()).unwrap();

            model.set_param_vec(&p_minus);
            let sigma_minus = vech::vech(&model.implied_cov()).unwrap();

            for i in 0..kstar {
                delta[(i, j)] = (sigma_plus[i] - sigma_minus[i]) / (2.0 * eps);
            }
        }
        model.set_param_vec(&params);
        delta
    }

    #[test]
    fn test_analytical_vs_numerical_cfa() {
        let model_str = "F1 =~ NA*V1 + V2 + V3\nF1 ~~ 1*F1\nV1 ~~ V1\nV2 ~~ V2\nV3 ~~ V3";
        let pt = parse_model(model_str, false).unwrap();
        let obs: Vec<String> = vec!["V1", "V2", "V3"]
            .into_iter()
            .map(String::from)
            .collect();
        let mut model = Model::from_partable(&pt, &obs);

        model.set_param_vec(&[0.7, 0.8, 0.6, 0.3, 0.4, 0.2]);

        let analytical = analytical_jacobian(&model);
        let numerical = numerical_jacobian(&mut model, 1e-7);

        assert_eq!(analytical.nrows(), numerical.nrows());
        assert_eq!(analytical.ncols(), numerical.ncols());

        for i in 0..analytical.nrows() {
            for j in 0..analytical.ncols() {
                let diff = (analytical[(i, j)] - numerical[(i, j)]).abs();
                assert!(
                    diff < 1e-5,
                    "Mismatch at ({}, {}): analytical={}, numerical={}, diff={}",
                    i,
                    j,
                    analytical[(i, j)],
                    numerical[(i, j)],
                    diff
                );
            }
        }
    }

    #[test]
    fn test_analytical_vs_numerical_2factor() {
        let model_str = "F1 =~ NA*V1 + V2 + V3\nF2 =~ NA*V4 + V5 + V6\n\
                          F1 ~~ 1*F1\nF2 ~~ 1*F2\nF1 ~~ F2\n\
                          V1 ~~ V1\nV2 ~~ V2\nV3 ~~ V3\nV4 ~~ V4\nV5 ~~ V5\nV6 ~~ V6";
        let pt = parse_model(model_str, false).unwrap();
        let obs: Vec<String> = (1..=6).map(|i| format!("V{i}")).collect();
        let mut model = Model::from_partable(&pt, &obs);

        let analytical = analytical_jacobian(&model);
        let numerical = numerical_jacobian(&mut model, 1e-7);

        for i in 0..analytical.nrows() {
            for j in 0..analytical.ncols() {
                let diff = (analytical[(i, j)] - numerical[(i, j)]).abs();
                assert!(
                    diff < 1e-5,
                    "Mismatch at ({}, {}): analytical={}, numerical={}, diff={}",
                    i,
                    j,
                    analytical[(i, j)],
                    numerical[(i, j)],
                    diff
                );
            }
        }
    }

    #[test]
    fn test_analytical_vs_numerical_with_psi_offdiag() {
        let model_str = "F1 =~ NA*V1 + V2\nF2 =~ NA*V3 + V4\n\
                          F1 ~~ 1*F1\nF2 ~~ 1*F2\nF1 ~~ F2\n\
                          V1 ~~ V1\nV2 ~~ V2\nV3 ~~ V3\nV4 ~~ V4";
        let pt = parse_model(model_str, false).unwrap();
        let obs: Vec<String> = (1..=4).map(|i| format!("V{i}")).collect();
        let mut model = Model::from_partable(&pt, &obs);

        let params = model.get_param_vec();
        let modified: Vec<f64> = params
            .iter()
            .enumerate()
            .map(|(i, _)| 0.3 + 0.1 * i as f64)
            .collect();
        model.set_param_vec(&modified);

        let analytical = analytical_jacobian(&model);
        let numerical = numerical_jacobian(&mut model, 1e-7);

        for i in 0..analytical.nrows() {
            for j in 0..analytical.ncols() {
                let diff = (analytical[(i, j)] - numerical[(i, j)]).abs();
                assert!(
                    diff < 1e-5,
                    "Mismatch at ({}, {}): analytical={}, numerical={}, diff={}",
                    i,
                    j,
                    analytical[(i, j)],
                    numerical[(i, j)],
                    diff
                );
            }
        }
    }

    #[test]
    fn test_analytical_vs_numerical_with_beta() {
        let model_str = "F1 =~ NA*V1 + V2 + V3\nF2 =~ NA*V4 + V5 + V6\n\
                          F2 ~ F1\n\
                          F1 ~~ 1*F1\nF2 ~~ F2\n\
                          V1 ~~ V1\nV2 ~~ V2\nV3 ~~ V3\nV4 ~~ V4\nV5 ~~ V5\nV6 ~~ V6";
        let pt = parse_model(model_str, false).unwrap();
        let obs: Vec<String> = (1..=6).map(|i| format!("V{i}")).collect();
        let mut model = Model::from_partable(&pt, &obs);

        let analytical = analytical_jacobian(&model);
        let numerical = numerical_jacobian(&mut model, 1e-7);

        for i in 0..analytical.nrows() {
            for j in 0..analytical.ncols() {
                let diff = (analytical[(i, j)] - numerical[(i, j)]).abs();
                assert!(
                    diff < 1e-5,
                    "Mismatch at ({}, {}): analytical={}, numerical={}, diff={}",
                    i,
                    j,
                    analytical[(i, j)],
                    numerical[(i, j)],
                    diff
                );
            }
        }
    }
}
