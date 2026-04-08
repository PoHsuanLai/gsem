use faer::Mat;

/// Reorder V matrix from user variable ordering to model variable ordering.
///
/// Port of `.rearrange()` from GenomicSEM's utils.R.
///
/// The V matrix is the sampling covariance of vech(S). Its rows/columns correspond
/// to elements of the lower triangle of S. When the SEM model reorders variables
/// internally, V must be permuted accordingly.
pub fn reorder_v(v: &Mat<f64>, user_order: &[String], model_order: &[String]) -> Mat<f64> {
    let k = user_order.len();
    assert_eq!(k, model_order.len());
    let kstar = k * (k + 1) / 2;
    assert_eq!(v.nrows(), kstar);
    assert_eq!(v.ncols(), kstar);

    // Build permutation: where does each user-ordered variable appear in model order?
    let perm: Vec<usize> = model_order
        .iter()
        .map(|name| {
            user_order
                .iter()
                .position(|n| n == name)
                .unwrap_or_else(|| panic!("variable {name} not found in user ordering"))
        })
        .collect();

    // Build vech index mapping
    // user_vech_idx maps (i,j) in user order -> vech position
    // We need to map model vech positions to user vech positions
    let user_vech = vech_index_map(k);
    let model_to_user: Vec<usize> = {
        let mut mapping = Vec::with_capacity(kstar);
        for j in 0..k {
            for i in j..k {
                // (i, j) in model order corresponds to (perm[i], perm[j]) in user order
                let ui = perm[i];
                let uj = perm[j];
                let (ri, rj) = if ui >= uj { (ui, uj) } else { (uj, ui) };
                let user_idx = user_vech[ri * (ri + 1) / 2 + rj];
                mapping.push(user_idx);
            }
        }
        mapping
    };

    // Apply permutation to V
    let mut v_reorder = Mat::zeros(kstar, kstar);
    for i in 0..kstar {
        for j in 0..kstar {
            v_reorder[(i, j)] = v[(model_to_user[i], model_to_user[j])];
        }
    }

    v_reorder
}

/// Build a mapping: for a k×k matrix, vech position of element (i,j) where i>=j.
/// Returns a flat array indexed by i*(i+1)/2 + j.
fn vech_index_map(k: usize) -> Vec<usize> {
    let mut map = vec![0usize; k * (k + 1) / 2 + k]; // extra space
    let mut idx = 0;
    for j in 0..k {
        for i in j..k {
            map[i * (i + 1) / 2 + j] = idx;
            idx += 1;
        }
    }
    map
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reorder_identity() {
        let v = faer::mat![[1.0, 0.5, 0.3], [0.5, 2.0, 0.4], [0.3, 0.4, 3.0],];
        let order = vec!["A".to_string(), "B".to_string()];
        let result = reorder_v(&v, &order, &order);
        for i in 0..3 {
            for j in 0..3 {
                assert!((result[(i, j)] - v[(i, j)]).abs() < 1e-15);
            }
        }
    }

    #[test]
    fn test_reorder_swap() {
        // 2 variables: vech has 3 elements: (0,0), (1,0), (1,1)
        let v = faer::mat![[1.0, 0.2, 0.3], [0.2, 2.0, 0.4], [0.3, 0.4, 3.0],];
        let user = vec!["A".to_string(), "B".to_string()];
        let model = vec!["B".to_string(), "A".to_string()];
        let result = reorder_v(&v, &user, &model);
        // After swapping: vech order should be (B,B), (A,B), (A,A)
        // which maps to original (1,1), (1,0), (0,0) = indices 2, 1, 0
        assert!((result[(0, 0)] - v[(2, 2)]).abs() < 1e-15);
        assert!((result[(2, 2)] - v[(0, 0)]).abs() < 1e-15);
    }
}
