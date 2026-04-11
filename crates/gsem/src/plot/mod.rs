//! SVG plot generation for GWAS, QQ, and parallel-analysis scree plots.
//!
//! Uses the `plotters` crate with the SVG backend (no system dependencies).
//! Each plot type exposes a free function that writes an SVG file from
//! already-prepared data.  The CLI calls these from its `--save-plot-*`
//! flags on the corresponding analysis subcommands.

pub mod manhattan;
pub mod qq;
pub mod scree;

/// Point in a Manhattan plot: (chromosome, base-pair, -log10(p)).
#[derive(Debug, Clone, Copy)]
pub struct ManhattanPoint {
    pub chr: u8,
    pub bp: u64,
    pub neg_log10_p: f64,
}

/// Expected / observed pair for a QQ plot (both on -log10 scale).
#[derive(Debug, Clone, Copy)]
pub struct QqPoint {
    pub expected: f64,
    pub observed: f64,
}

/// Compute -log10(p) with p clamped to avoid -inf on exactly-zero p values.
pub fn neg_log10_p(p: f64) -> f64 {
    let clamped = p.max(1e-300);
    -clamped.log10()
}

/// Genomic inflation factor lambda = median(chi^2) / 0.4549364.
///
/// Converts two-sided p-values to chi^2(1) statistics via the identity
/// `qchisq(p, 1, lower=F) == qnorm(1 - p/2)^2`, which uses
/// `Normal::inverse_cdf` and is numerically stable across the full
/// range. The `ChiSquared::inverse_cdf` path used earlier hit a
/// statrs 0.18 NaN regression for `p ≳ 0.97`, silently dropping
/// ~3% of SNPs and biasing lambda_gc downward.
pub fn lambda_gc(p_values: &[f64]) -> f64 {
    use statrs::distribution::{ContinuousCDF, Normal};
    let normal = Normal::standard();
    let mut chisq: Vec<f64> = p_values
        .iter()
        .filter(|p| p.is_finite() && **p > 0.0 && **p <= 1.0)
        .map(|p| {
            let z = normal.inverse_cdf(1.0 - *p / 2.0);
            z * z
        })
        .filter(|q| q.is_finite())
        .collect();
    if chisq.is_empty() {
        return f64::NAN;
    }
    chisq.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = chisq.len() / 2;
    let median = if chisq.len().is_multiple_of(2) {
        0.5 * (chisq[mid - 1] + chisq[mid])
    } else {
        chisq[mid]
    };
    median / 0.4549364
}

/// Build QQ plot points from observed p-values.
///
/// Expected p-values are uniform quantiles `(i + 0.5) / n`.  Output is
/// sorted by observed ascending (smallest observed = largest -log10 at the
/// right of the plot in standard orientation).
pub fn build_qq_points(p_values: &[f64]) -> Vec<QqPoint> {
    let mut ps: Vec<f64> = p_values
        .iter()
        .copied()
        .filter(|p| p.is_finite() && *p > 0.0 && *p <= 1.0)
        .collect();
    if ps.is_empty() {
        return Vec::new();
    }
    ps.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = ps.len() as f64;
    ps.iter()
        .enumerate()
        .map(|(i, &p)| {
            let expected = (i as f64 + 0.5) / n;
            QqPoint {
                expected: neg_log10_p(expected),
                observed: neg_log10_p(p),
            }
        })
        .collect()
}

/// Layout info for Manhattan plot: per-chromosome (offset, midpoint).
///
/// Chromosomes are laid out left-to-right in numeric order, each
/// starting at the previous chromosome's end.  Returned vector is
/// indexed by chromosome in the order they appear in `points`.
#[derive(Debug, Clone)]
pub struct ChromLayout {
    /// (chr, cumulative_offset, midpoint_x) per chromosome present
    pub chroms: Vec<(u8, f64, f64)>,
    /// Total x-axis extent
    pub total_x: f64,
}

/// Compute chromosome offsets for Manhattan layout.
///
/// Points must be non-empty.  Returns cumulative offsets so chromosome 1
/// starts at 0, chromosome 2 starts at max(bp on chr1), etc.
pub fn compute_chrom_layout(points: &[ManhattanPoint]) -> ChromLayout {
    let mut chroms: std::collections::BTreeMap<u8, u64> = std::collections::BTreeMap::new();
    for p in points {
        let e = chroms.entry(p.chr).or_insert(0);
        if p.bp > *e {
            *e = p.bp;
        }
    }
    let mut out = Vec::with_capacity(chroms.len());
    let mut offset: f64 = 0.0;
    for (chr, max_bp) in chroms {
        let span = max_bp as f64;
        let mid = offset + span * 0.5;
        out.push((chr, offset, mid));
        offset += span;
    }
    ChromLayout {
        chroms: out,
        total_x: offset,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_neg_log10_p_zero() {
        // Clamps to avoid -inf
        assert!(neg_log10_p(0.0).is_finite());
        assert!(neg_log10_p(0.0) > 200.0);
    }

    #[test]
    fn test_lambda_gc_null() {
        // Uniform p-values should give lambda near 1 (within ~10% at this sample size)
        let ps: Vec<f64> = (1..=999).map(|i| i as f64 / 1000.0).collect();
        let lam = lambda_gc(&ps);
        assert!((lam - 1.0).abs() < 0.1, "lambda = {lam}");
    }

    #[test]
    fn test_build_qq_points_sorted() {
        let ps = vec![0.5, 0.01, 0.1, 0.001];
        let qq = build_qq_points(&ps);
        assert_eq!(qq.len(), 4);
        // Observed should be monotonically non-decreasing (smallest p first → largest -log10)
        // Actually we sorted p ascending, so observed neg_log10 should be descending
        for i in 1..qq.len() {
            assert!(qq[i].observed <= qq[i - 1].observed);
        }
    }

    #[test]
    fn test_chrom_layout_offsets() {
        let pts = vec![
            ManhattanPoint {
                chr: 1,
                bp: 100,
                neg_log10_p: 1.0,
            },
            ManhattanPoint {
                chr: 1,
                bp: 500,
                neg_log10_p: 2.0,
            },
            ManhattanPoint {
                chr: 2,
                bp: 200,
                neg_log10_p: 3.0,
            },
        ];
        let layout = compute_chrom_layout(&pts);
        assert_eq!(layout.chroms.len(), 2);
        assert_eq!(layout.chroms[0].0, 1);
        assert_eq!(layout.chroms[1].0, 2);
        // chr1 starts at 0
        assert_eq!(layout.chroms[0].1, 0.0);
        // chr2 starts at chr1's max bp = 500
        assert_eq!(layout.chroms[1].1, 500.0);
        // total extent = 500 + 200
        assert_eq!(layout.total_x, 700.0);
    }
}
