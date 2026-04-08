use faer::Mat;
use gsem_matrix::vech;

use super::gc_correction::{self, GcMode};

/// Build the full S matrix (k+1 × k+1) for a single SNP.
///
/// Port of `.get_S_Full()` from GenomicSEM's utils.R.
///
/// Structure:
///   S_Full[0,0] = var_snp
///   S_Full[0, 1..=k] = var_snp * beta_snp (SNP-phenotype covariances)
///   S_Full[1..=k, 1..=k] = S_LD
pub fn build_s_full(s_ld: &Mat<f64>, beta_snp: &[f64], var_snp: f64, k: usize) -> Mat<f64> {
    let n = k + 1;
    let mut s_full = Mat::zeros(n, n);

    // SNP variance
    s_full[(0, 0)] = var_snp;

    // SNP-phenotype covariances
    for p in 0..k {
        let cov = var_snp * beta_snp[p];
        s_full[(0, p + 1)] = cov;
        s_full[(p + 1, 0)] = cov;
    }

    // LD block
    for i in 0..k {
        for j in 0..k {
            s_full[(i + 1, j + 1)] = s_ld[(i, j)];
        }
    }

    s_full
}

/// Build the full V matrix ((k+1)*(k+2)/2 × (k+1)*(k+2)/2) for a single SNP.
///
/// Port of `.get_V_full()` from GenomicSEM's utils.R.
///
/// Structure:
///   V_Full[0,0] = var_snp_se2
///   V_Full[1..=k, 1..=k] = V_SNP (from GC correction)
///   V_Full[(k+1).., (k+1)..] = V_LD
pub fn build_v_full(
    v_ld: &Mat<f64>,
    se_snp: &[f64],
    var_snp: f64,
    var_snp_se2: f64,
    i_ld: &Mat<f64>,
    gc: GcMode,
    k: usize,
) -> Mat<f64> {
    let kstar_full = (k + 1) * (k + 2) / 2;
    let kstar_ld = k * (k + 1) / 2;
    let mut v_full = Mat::zeros(kstar_full, kstar_full);

    // V_SNP block
    let v_snp = gc_correction::build_v_snp(se_snp, i_ld, var_snp, gc, k);
    let v_snp_vech = vech::vech(&v_snp);

    // Position 0: variance of SNP variance estimate
    v_full[(0, 0)] = var_snp_se2;

    // Positions 1..=k: V_SNP (as vech of k×k)
    for i in 0..v_snp_vech.len().min(k) {
        for j in 0..v_snp_vech.len().min(k) {
            if (1 + i) < kstar_full && (1 + j) < kstar_full {
                v_full[(1 + i, 1 + j)] = v_snp[(i.min(k - 1), j.min(k - 1))];
            }
        }
    }

    // Simple approach: put V_SNP diagonal elements
    for i in 0..k {
        v_full[(1 + i, 1 + i)] = v_snp[(i, i)];
        for j in (i + 1)..k {
            v_full[(1 + i, 1 + j)] = v_snp[(i, j)];
            v_full[(1 + j, 1 + i)] = v_snp[(j, i)];
        }
    }

    // V_LD block
    let ld_offset = k + 1;
    for i in 0..kstar_ld {
        for j in 0..kstar_ld {
            if (ld_offset + i) < kstar_full && (ld_offset + j) < kstar_full {
                v_full[(ld_offset + i, ld_offset + j)] = v_ld[(i, j)];
            }
        }
    }

    v_full
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_s_full_structure() {
        let s_ld = faer::mat![[0.3, 0.1], [0.1, 0.4],];
        let beta_snp = vec![0.05, -0.03];
        let var_snp = 0.25;
        let s = build_s_full(&s_ld, &beta_snp, var_snp, 2);

        assert_eq!(s.nrows(), 3);
        assert_eq!(s.ncols(), 3);
        assert!((s[(0, 0)] - 0.25).abs() < 1e-10);
        assert!((s[(0, 1)] - 0.25 * 0.05).abs() < 1e-10);
        assert!((s[(1, 0)] - 0.25 * 0.05).abs() < 1e-10);
        assert!((s[(1, 1)] - 0.3).abs() < 1e-10);
        // Symmetric
        assert!((s[(0, 2)] - s[(2, 0)]).abs() < 1e-15);
    }

    #[test]
    fn test_build_v_full_dimensions() {
        let v_ld = faer::mat![
            [0.01, 0.001, 0.002],
            [0.001, 0.02, 0.003],
            [0.002, 0.003, 0.03],
        ];
        let se = vec![0.1, 0.2];
        let i_ld = faer::mat![[1.05, 0.02], [0.02, 1.03]];
        let v = build_v_full(&v_ld, &se, 0.25, 0.001, &i_ld, GcMode::Standard, 2);
        // k=2, kstar_full = 3*4/2 = 6
        assert_eq!(v.nrows(), 6);
        assert_eq!(v.ncols(), 6);
    }
}
