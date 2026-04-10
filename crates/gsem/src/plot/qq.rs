//! QQ plot for GWAS p-values.

use anyhow::Result;
use plotters::prelude::*;
use std::path::Path;

use super::{build_qq_points, lambda_gc};

/// Render a QQ plot of observed vs expected -log10(p) values.
///
/// The genomic inflation factor lambda (median chi^2 / 0.4549364) is
/// printed in the plot title.
pub fn qq_plot(p_values: &[f64], title: &str, out: &Path) -> Result<()> {
    let points = build_qq_points(p_values);
    if points.is_empty() {
        anyhow::bail!("qq_plot: no valid p-values");
    }
    let lambda = lambda_gc(p_values);

    let root = SVGBackend::new(out, (800, 800)).into_drawing_area();
    root.fill(&WHITE)
        .map_err(|e| anyhow::anyhow!("fill: {e}"))?;

    let x_max = points.iter().map(|p| p.expected).fold(0.0f64, f64::max) * 1.05;
    let y_max = points
        .iter()
        .map(|p| p.observed)
        .fold(0.0f64, f64::max)
        .max(x_max)
        * 1.05;

    let caption = if lambda.is_finite() {
        format!("{title}  (λ = {lambda:.3})")
    } else {
        title.to_string()
    };

    let mut chart = ChartBuilder::on(&root)
        .caption(caption, ("sans-serif", 24).into_font())
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(0.0f64..x_max, 0.0f64..y_max)
        .map_err(|e| anyhow::anyhow!("build chart: {e}"))?;

    chart
        .configure_mesh()
        .x_desc("Expected -log10(p)")
        .y_desc("Observed -log10(p)")
        .draw()
        .map_err(|e| anyhow::anyhow!("mesh: {e}"))?;

    // y = x reference line
    let diag_max = x_max.min(y_max);
    chart
        .draw_series(LineSeries::new(
            [(0.0, 0.0), (diag_max, diag_max)],
            RED.mix(0.6).stroke_width(1),
        ))
        .map_err(|e| anyhow::anyhow!("diagonal: {e}"))?;

    // QQ points
    chart
        .draw_series(
            points
                .iter()
                .map(|p| Circle::new((p.expected, p.observed), 2, BLUE.filled())),
        )
        .map_err(|e| anyhow::anyhow!("points: {e}"))?;

    root.present()
        .map_err(|e| anyhow::anyhow!("present: {e}"))?;

    log::info!(
        "QQ plot written to {} ({} points, lambda = {:.3})",
        out.display(),
        points.len(),
        lambda
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_qq_plot_smoke() {
        // Mix of small and uniform p-values
        let mut ps: Vec<f64> = (1..=500).map(|i| i as f64 / 500.0).collect();
        ps.push(1e-8);
        ps.push(1e-6);
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("svg");
        qq_plot(&ps, "Test QQ", &path).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("<svg") || content.starts_with("<?xml"));
        assert!(content.contains("Expected"));
    }
}
