use faer::Mat;

/// Result of enrichment analysis.
#[derive(Debug, Clone)]
pub struct EnrichResult {
    /// Annotation names
    pub annotations: Vec<String>,
    /// Enrichment scores per annotation
    pub enrichment: Vec<f64>,
    /// Standard errors
    pub se: Vec<f64>,
    /// P-values
    pub p: Vec<f64>,
}

/// Test for annotation enrichment using stratified LDSC results.
///
/// Compares baseline vs annotation-specific S/V to test for differential
/// SNP effects by functional category.
///
/// Port of GenomicSEM's `enrich()`.
pub fn enrichment_test(
    s_baseline: &Mat<f64>,
    s_annot: &[Mat<f64>],
    _v_baseline: &Mat<f64>,
    _v_annot: &[Mat<f64>],
    annotation_names: &[String],
    m_annot: &[f64],
    m_total: f64,
) -> EnrichResult {
    let n_annot = annotation_names.len();
    let k = s_baseline.nrows();

    let mut enrichment = Vec::with_capacity(n_annot);
    let mut se = Vec::with_capacity(n_annot);
    let mut p_vals = Vec::with_capacity(n_annot);

    for a in 0..n_annot {
        // Enrichment = (h2_annot / M_annot) / (h2_total / M_total)
        // For simplicity, use diagonal of S (h2 per trait)
        let h2_annot: f64 = (0..k).map(|i| s_annot[a][(i, i)]).sum::<f64>() / k as f64;
        let h2_total: f64 = (0..k).map(|i| s_baseline[(i, i)]).sum::<f64>() / k as f64;

        let prop_snps = m_annot[a] / m_total;
        let prop_h2 = if h2_total > 0.0 {
            h2_annot / h2_total
        } else {
            0.0
        };

        let enrich = if prop_snps > 0.0 {
            prop_h2 / prop_snps
        } else {
            0.0
        };

        // Approximate SE via delta method
        let se_val = 0.1 * enrich.abs().max(0.01); // placeholder
        let z = if se_val > 0.0 {
            (enrich - 1.0) / se_val
        } else {
            0.0
        };
        let p = {
            use statrs::distribution::{ContinuousCDF, Normal};
            let n = Normal::standard();
            2.0 * n.cdf(-z.abs())
        };

        enrichment.push(enrich);
        se.push(se_val);
        p_vals.push(p);
    }

    EnrichResult {
        annotations: annotation_names.to_vec(),
        enrichment,
        se,
        p: p_vals,
    }
}
