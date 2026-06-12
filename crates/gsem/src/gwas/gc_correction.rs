use std::str::FromStr;

use faer::Mat;

/// Genomic control correction mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GcMode {
    Conservative,
    Standard,
    None,
}

impl GcMode {
    /// GC scaling factor for a diagonal element of the LD intercept matrix.
    fn scale(self, i_ld_diag: f64) -> f64 {
        match self {
            GcMode::Conservative => i_ld_diag,
            GcMode::Standard => i_ld_diag.sqrt(),
            GcMode::None => 1.0,
        }
    }
}

impl FromStr for GcMode {
    type Err = std::convert::Infallible;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s.to_lowercase().as_str() {
            "conserv" | "conservative" => GcMode::Conservative,
            "none" | "off" => GcMode::None,
            _ => GcMode::Standard,
        })
    }
}

/// Floor the diagonal of the LD-score intercept matrix at 1.0.
///
/// Port of `diag(I_LD) <- ifelse(diag(I_LD) <= 1, 1, diag(I_LD))` from R
/// GenomicSEM (`userGWAS.R:156`, `commonfactorGWAS.R:86`). The LDSC intercept
/// can fall slightly below 1 from sampling noise; R clamps it before building
/// the per-SNP sampling-variance matrix so SNP variances (and Q_SNP) are not
/// deflated. Off-diagonal elements are left unchanged.
pub fn floor_intercept_diag(i_ld: &Mat<f64>) -> Mat<f64> {
    let n = i_ld.nrows();
    Mat::from_fn(n, n, |r, c| {
        if r == c {
            i_ld[(r, c)].max(1.0)
        } else {
            i_ld[(r, c)]
        }
    })
}

/// Construct V_SNP matrix (k×k) for a single SNP with genomic control correction.
///
/// Port of `.get_V_SNP()` from GenomicSEM's utils.R.
pub fn build_v_snp(
    se_snp: &[f64],
    i_ld: &Mat<f64>,
    var_snp: f64,
    gc: GcMode,
    k: usize,
) -> Mat<f64> {
    let mut v_snp = Mat::zeros(k, k);

    for x in 0..k {
        for y in x..k {
            let val = if x == y {
                (se_snp[x] * var_snp * gc.scale(i_ld[(x, x)])).powi(2)
            } else {
                let base = se_snp[x] * se_snp[y] * var_snp.powi(2);
                base * i_ld[(x, y)] * gc.scale(i_ld[(x, x)]) * gc.scale(i_ld[(y, y)])
            };
            v_snp[(x, y)] = val;
            v_snp[(y, x)] = val;
        }
    }

    v_snp
}

/// Compute GC-adjusted Z-scores for a single SNP.
///
/// Port of `.get_Z_pre()` from GenomicSEM's utils.R.
pub fn gc_adjusted_z(beta: &[f64], se: &[f64], i_ld: &Mat<f64>, gc: GcMode, k: usize) -> Vec<f64> {
    (0..k)
        .map(|x| {
            let denom = se[x] * gc.scale(i_ld[(x, x)]);
            if denom.abs() > 1e-30 {
                beta[x] / denom
            } else {
                0.0
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gc_none_diagonal() {
        let se = vec![0.1, 0.2];
        let i_ld = faer::mat![[1.05, 0.02], [0.02, 1.03]];
        let v = build_v_snp(&se, &i_ld, 0.5, GcMode::None, 2);
        // Diagonal: (se * var_snp)^2 = (0.1 * 0.5)^2 = 0.0025
        assert!((v[(0, 0)] - 0.0025).abs() < 1e-10);
    }

    #[test]
    fn test_floor_intercept_diag() {
        // Diagonal entries <= 1 are clamped to 1; entries > 1 and all
        // off-diagonals are untouched (matches R userGWAS.R:156).
        let i_ld = faer::mat![[0.993, -0.01], [-0.01, 1.0024]];
        let floored = floor_intercept_diag(&i_ld);
        assert!(
            (floored[(0, 0)] - 1.0).abs() < 1e-15,
            "0.993 should floor to 1"
        );
        assert!(
            (floored[(1, 1)] - 1.0024).abs() < 1e-15,
            "1.0024 should be unchanged"
        );
        assert!(
            (floored[(0, 1)] - (-0.01)).abs() < 1e-15,
            "off-diagonal unchanged"
        );
        assert!((floored[(1, 0)] - (-0.01)).abs() < 1e-15);
    }

    #[test]
    fn test_gc_conservative_inflates() {
        let se = vec![0.1];
        let i_ld = faer::mat![[1.1]];
        let v_none = build_v_snp(&se, &i_ld, 0.5, GcMode::None, 1);
        let v_conserv = build_v_snp(&se, &i_ld, 0.5, GcMode::Conservative, 1);
        // Conservative should give larger variance
        assert!(v_conserv[(0, 0)] > v_none[(0, 0)]);
    }

    #[test]
    fn test_gc_adjusted_z() {
        let beta = vec![0.5];
        let se = vec![0.1];
        let i_ld = faer::mat![[1.0]];
        let z = gc_adjusted_z(&beta, &se, &i_ld, GcMode::None, 1);
        assert!((z[0] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_parse_gc_mode() {
        assert_eq!(
            "conservative".parse::<GcMode>().unwrap(),
            GcMode::Conservative
        );
        assert_eq!("conserv".parse::<GcMode>().unwrap(), GcMode::Conservative);
        assert_eq!("standard".parse::<GcMode>().unwrap(), GcMode::Standard);
        assert_eq!("none".parse::<GcMode>().unwrap(), GcMode::None);
        assert_eq!("unknown".parse::<GcMode>().unwrap(), GcMode::Standard);
    }
}
