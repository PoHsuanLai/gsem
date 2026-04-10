//! Scree plot for parallel-analysis eigenvalues.

use anyhow::{Context, Result};
use plotters::prelude::*;
use std::path::Path;

/// One row of the parallel-analysis table.
#[derive(Debug, Clone, Copy)]
pub struct ScreePoint {
    pub factor: usize,
    pub observed: f64,
    pub simulated_95: f64,
}

/// Render a scree plot to SVG comparing observed vs 95th-percentile
/// simulated eigenvalues.  The suggested number of factors (where the
/// observed first falls below the simulated) is drawn as a vertical
/// reference line.
pub fn scree_plot(points: &[ScreePoint], title: &str, out: &Path) -> Result<()> {
    if points.is_empty() {
        anyhow::bail!("scree_plot: no data points");
    }

    let root = SVGBackend::new(out, (800, 600)).into_drawing_area();
    root.fill(&WHITE)
        .map_err(|e| anyhow::anyhow!("fill: {e}"))?;

    let x_max = points.len() as f64 + 0.5;
    let y_max = points
        .iter()
        .flat_map(|p| [p.observed, p.simulated_95])
        .fold(f64::NEG_INFINITY, f64::max)
        .max(1.0)
        * 1.1;
    let y_min = points
        .iter()
        .flat_map(|p| [p.observed, p.simulated_95])
        .fold(f64::INFINITY, f64::min)
        .min(0.0);

    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 24).into_font())
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(0.5f64..x_max, y_min..y_max)
        .map_err(|e| anyhow::anyhow!("build chart: {e}"))?;

    chart
        .configure_mesh()
        .x_desc("Factor")
        .y_desc("Eigenvalue")
        .x_labels(points.len().min(20))
        .disable_x_mesh()
        .draw()
        .map_err(|e| anyhow::anyhow!("mesh: {e}"))?;

    // Observed line (solid blue)
    chart
        .draw_series(LineSeries::new(
            points
                .iter()
                .map(|p| (p.factor as f64, p.observed)),
            BLUE.stroke_width(2),
        ))
        .map_err(|e| anyhow::anyhow!("observed line: {e}"))?
        .label("Observed")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE.stroke_width(2)));

    // Point markers for observed
    chart
        .draw_series(
            points
                .iter()
                .map(|p| Circle::new((p.factor as f64, p.observed), 4, BLUE.filled())),
        )
        .map_err(|e| anyhow::anyhow!("observed markers: {e}"))?;

    // Simulated 95th percentile line (dashed red)
    chart
        .draw_series(LineSeries::new(
            points
                .iter()
                .map(|p| (p.factor as f64, p.simulated_95)),
            RED.stroke_width(2),
        ))
        .map_err(|e| anyhow::anyhow!("simulated line: {e}"))?
        .label("Simulated (95%)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED.stroke_width(2)));

    chart
        .draw_series(
            points
                .iter()
                .map(|p| TriangleMarker::new((p.factor as f64, p.simulated_95), 5, RED.filled())),
        )
        .map_err(|e| anyhow::anyhow!("simulated markers: {e}"))?;

    // Suggested factor count: first k where observed < simulated
    if let Some(n_factors) = points
        .iter()
        .position(|p| p.observed < p.simulated_95)
    {
        // Use the factor *before* the crossover as the suggested count
        let suggest_x = n_factors as f64 + 0.5;
        if suggest_x >= 0.5 && suggest_x <= x_max {
            chart
                .draw_series(std::iter::once(PathElement::new(
                    vec![(suggest_x, y_min), (suggest_x, y_max)],
                    BLACK.mix(0.4).stroke_width(1),
                )))
                .map_err(|e| anyhow::anyhow!("suggest line: {e}"))?;
        }
    }

    chart
        .configure_series_labels()
        .position(SeriesLabelPosition::UpperRight)
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()
        .map_err(|e| anyhow::anyhow!("legend: {e}"))?;

    root.present()
        .map_err(|e| anyhow::anyhow!("present: {e}"))?;

    log::info!("scree plot written to {}", out.display());
    Ok(())
}

/// Parse a parallel-analysis TSV file.  The CLI's paLDSC subcommand
/// writes three columns: `factor`, `observed`, `simulated_95` (plus a
/// trailing comment line starting with `#`).
pub fn read_scree_tsv(path: &Path) -> Result<Vec<ScreePoint>> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read {}", path.display()))?;
    let mut out = Vec::new();
    for (i, line) in content.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') || i == 0 {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        let factor: usize = match fields[0].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let observed: f64 = match fields[1].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let simulated_95: f64 = match fields[2].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        out.push(ScreePoint { factor, observed, simulated_95 });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scree_plot_smoke() {
        let points = vec![
            ScreePoint { factor: 1, observed: 3.5, simulated_95: 1.2 },
            ScreePoint { factor: 2, observed: 2.0, simulated_95: 1.1 },
            ScreePoint { factor: 3, observed: 0.9, simulated_95: 1.0 },
            ScreePoint { factor: 4, observed: 0.6, simulated_95: 0.95 },
        ];
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("svg");
        scree_plot(&points, "Test Scree", &path).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.starts_with("<?xml") || content.starts_with("<svg"));
        assert!(content.contains("Observed"));
        assert!(content.contains("Simulated"));
    }
}
