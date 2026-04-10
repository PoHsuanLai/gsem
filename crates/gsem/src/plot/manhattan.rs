//! Manhattan plot for GWAS results.

use anyhow::Result;
use plotters::prelude::*;
use std::path::Path;

use super::{compute_chrom_layout, ManhattanPoint};

/// Genome-wide significance threshold (5e-8) on -log10 scale.
pub const GWS_NEG_LOG10: f64 = 7.30103;
/// Suggestive significance threshold (1e-5) on -log10 scale.
pub const SUGGESTIVE_NEG_LOG10: f64 = 5.0;

/// Render a Manhattan plot to SVG.
///
/// Points are laid out left-to-right by cumulative genomic position.
/// Alternating chromosome colors improve readability; genome-wide
/// significance (5e-8) and suggestive (1e-5) thresholds are drawn
/// as horizontal reference lines.
///
/// `max_points` thins the plot if the raw input exceeds this many
/// points, keeping all significant SNPs (p < 1e-4) and randomly
/// sampling the rest.  Passing `None` disables thinning.
pub fn manhattan_plot(
    points: &[ManhattanPoint],
    title: &str,
    out: &Path,
    max_points: Option<usize>,
) -> Result<()> {
    if points.is_empty() {
        anyhow::bail!("manhattan_plot: no data points");
    }

    // Optional thinning: keep everything above suggestive, thin below
    let thinned_owned;
    let points_ref: &[ManhattanPoint] = if let Some(cap) = max_points {
        if points.len() > cap {
            thinned_owned = thin_points(points, cap);
            &thinned_owned
        } else {
            points
        }
    } else {
        points
    };

    let layout = compute_chrom_layout(points_ref);
    if layout.total_x <= 0.0 {
        anyhow::bail!("manhattan_plot: zero x-axis extent");
    }
    // Build a fast lookup for chromosome offsets
    let chrom_offsets: std::collections::HashMap<u8, f64> =
        layout.chroms.iter().map(|(c, off, _)| (*c, *off)).collect();

    let y_max = points_ref
        .iter()
        .map(|p| p.neg_log10_p)
        .fold(0.0f64, f64::max)
        .max(GWS_NEG_LOG10 + 1.0)
        * 1.05;

    let root = SVGBackend::new(out, (1200, 600)).into_drawing_area();
    root.fill(&WHITE)
        .map_err(|e| anyhow::anyhow!("fill: {e}"))?;

    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 24).into_font())
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(0.0f64..layout.total_x, 0.0f64..y_max)
        .map_err(|e| anyhow::anyhow!("build chart: {e}"))?;

    chart
        .configure_mesh()
        .y_desc("-log10(p)")
        .x_desc("Chromosome")
        .x_labels(0) // we draw custom labels below
        .disable_x_mesh()
        .draw()
        .map_err(|e| anyhow::anyhow!("mesh: {e}"))?;

    // Draw points colored by chromosome parity
    let color_even = RGBColor(50, 100, 180);
    let color_odd = RGBColor(220, 120, 50);
    chart
        .draw_series(points_ref.iter().map(|p| {
            let offset = chrom_offsets.get(&p.chr).copied().unwrap_or(0.0);
            let x = offset + p.bp as f64;
            let color = if p.chr % 2 == 0 { color_even } else { color_odd };
            Circle::new((x, p.neg_log10_p), 1, color.filled())
        }))
        .map_err(|e| anyhow::anyhow!("points: {e}"))?;

    // Suggestive line (1e-5)
    chart
        .draw_series(std::iter::once(PathElement::new(
            vec![(0.0, SUGGESTIVE_NEG_LOG10), (layout.total_x, SUGGESTIVE_NEG_LOG10)],
            BLUE.mix(0.6).stroke_width(1),
        )))
        .map_err(|e| anyhow::anyhow!("suggestive: {e}"))?;

    // Genome-wide significance line (5e-8)
    chart
        .draw_series(std::iter::once(PathElement::new(
            vec![(0.0, GWS_NEG_LOG10), (layout.total_x, GWS_NEG_LOG10)],
            RED.mix(0.7).stroke_width(1),
        )))
        .map_err(|e| anyhow::anyhow!("gws line: {e}"))?;

    // Custom x-axis labels at chromosome midpoints
    chart
        .draw_series(layout.chroms.iter().map(|(chr, _off, mid)| {
            Text::new(
                chr.to_string(),
                (*mid, 0.0),
                ("sans-serif", 12).into_font().color(&BLACK),
            )
        }))
        .map_err(|e| anyhow::anyhow!("chr labels: {e}"))?;

    root.present()
        .map_err(|e| anyhow::anyhow!("present: {e}"))?;

    log::info!(
        "Manhattan plot written to {} ({} points across {} chromosomes)",
        out.display(),
        points_ref.len(),
        layout.chroms.len()
    );
    Ok(())
}

/// Thin points to `cap`, keeping all with -log10(p) >= 4 (i.e. p <= 1e-4)
/// and uniformly sampling the rest.
fn thin_points(points: &[ManhattanPoint], cap: usize) -> Vec<ManhattanPoint> {
    let (significant, rest): (Vec<_>, Vec<_>) = points
        .iter()
        .copied()
        .partition(|p| p.neg_log10_p >= 4.0);
    if significant.len() >= cap {
        return significant;
    }
    let remaining = cap - significant.len();
    if rest.len() <= remaining {
        let mut out = significant;
        out.extend(rest);
        return out;
    }
    // Deterministic stride sampling so the plot is reproducible
    let stride = rest.len() as f64 / remaining as f64;
    let mut out = significant;
    let mut i = 0.0;
    while (i as usize) < rest.len() && out.len() < cap {
        out.push(rest[i as usize]);
        i += stride;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_thin_points_keeps_significant() {
        let mut points = Vec::new();
        // 10 significant points (p < 1e-4)
        for i in 0..10 {
            points.push(ManhattanPoint { chr: 1, bp: i * 1000, neg_log10_p: 8.0 });
        }
        // 100 null points
        for i in 0..100 {
            points.push(ManhattanPoint { chr: 1, bp: 100000 + i * 1000, neg_log10_p: 0.5 });
        }
        let thinned = thin_points(&points, 30);
        assert!(thinned.len() <= 30);
        // All significant points should remain
        let sig_count = thinned.iter().filter(|p| p.neg_log10_p >= 4.0).count();
        assert_eq!(sig_count, 10);
    }

    #[test]
    fn test_manhattan_plot_smoke() {
        let mut points = Vec::new();
        for chr in 1..=3 {
            for i in 0..20 {
                points.push(ManhattanPoint {
                    chr,
                    bp: (i as u64) * 1_000_000,
                    neg_log10_p: (i as f64) * 0.3,
                });
            }
        }
        // Add a GW-significant hit
        points.push(ManhattanPoint { chr: 2, bp: 15_000_000, neg_log10_p: 9.0 });

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("svg");
        manhattan_plot(&points, "Test Manhattan", &path, None).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("<svg") || content.starts_with("<?xml"));
        assert!(content.contains("-log10"));
    }
}
