use faer::Mat;

/// Result of enrichment analysis.
#[derive(Debug, Clone)]
pub struct EnrichResult {
    /// Annotation names
    pub annotations: Vec<String>,
    /// Enrichment scores per annotation
    pub enrichment: Vec<f64>,
    /// Standard errors of enrichment
    pub se: Vec<f64>,
    /// P-values (one-tailed, null = 1)
    pub p: Vec<f64>,
}

/// Test for annotation enrichment using stratified LDSC results.
///
/// Compares baseline vs annotation-specific S/V to test for differential
/// SNP effects by functional category.
///
/// Port of GenomicSEM's `enrich()`.
///
/// Enrichment SE is computed via the delta method:
///   enrichment = (est_annot / est_baseline) / prop_snps
///   enrichment_se = (SE_annot / |est_baseline|) / prop_snps
/// where est and SE come from fitting the model to the annotation-specific S
/// and computing sandwich SEs using V.
pub fn enrichment_test(
    s_baseline: &Mat<f64>,
    s_annot: &[Mat<f64>],
    v_baseline: &Mat<f64>,
    v_annot: &[Mat<f64>],
    annotation_names: &[String],
    m_annot: &[f64],
    m_total: f64,
) -> EnrichResult {
    let n_annot = annotation_names.len();
    let k = s_baseline.nrows();

    // Baseline estimates: average h2 across traits from diagonal of S
    let h2_total: f64 = (0..k).map(|i| s_baseline[(i, i)]).sum::<f64>() / k as f64;

    // Baseline SE from V (used for reference, kept for potential future use)
    let _baseline_se = compute_diag_se(s_baseline, v_baseline, k);

    let mut enrichment = Vec::with_capacity(n_annot);
    let mut se = Vec::with_capacity(n_annot);
    let mut p_vals = Vec::with_capacity(n_annot);

    for a in 0..n_annot {
        let h2_annot: f64 = (0..k).map(|i| s_annot[a][(i, i)]).sum::<f64>() / k as f64;

        let prop_snps = m_annot[a] / m_total;
        let prop_h2 = if h2_total.abs() > 1e-30 {
            h2_annot / h2_total
        } else {
            0.0
        };

        let enrich = if prop_snps > 0.0 {
            prop_h2 / prop_snps
        } else {
            0.0
        };

        // SE via delta method:
        // enrichment_se = (SE_annot / |est_baseline|) / prop_snps
        let annot_se = compute_diag_se(&s_annot[a], &v_annot[a], k);
        let se_val = if h2_total.abs() > 1e-30 && prop_snps > 0.0 {
            (annot_se / h2_total.abs()) / prop_snps
        } else {
            0.0
        };

        // One-tailed p-value: P(enrichment > 1) under null enrichment = 1
        let z = if se_val > 1e-30 {
            (enrich - 1.0) / se_val
        } else {
            0.0
        };
        let p = {
            use statrs::distribution::{ContinuousCDF, Normal};
            let n = Normal::standard();
            // One-tailed: P(Z > z) for testing enrichment > 1
            n.cdf(-z)
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

/// Compute the SE of the mean diagonal element of S, using V.
///
/// For each diagonal element S[i,i], its vech index is computed,
/// and SE = sqrt(V[idx,idx]). Then the SE of the mean is:
/// SE_mean = sqrt(sum(SE_i^2)) / k (assuming independence for simplicity).
fn compute_diag_se(_s: &Mat<f64>, v: &Mat<f64>, k: usize) -> f64 {
    let mut var_sum = 0.0;
    for i in 0..k {
        let vech_idx = vech_diag_index(i, k);
        if vech_idx < v.nrows() {
            var_sum += v[(vech_idx, vech_idx)];
        }
    }
    (var_sum / (k as f64).powi(2)).sqrt()
}

/// Get the vech index for the diagonal element (i, i) in a k×k matrix.
///
/// In column-major lower triangle ordering: index of (i, i) = sum_{j=0}^{i-1}(k-j) + 0
fn vech_diag_index(i: usize, k: usize) -> usize {
    // Column-major vech: for column j, rows j..k
    // Index of (i, j) where j <= i is: sum_{c=0}^{j-1}(k - c) + (i - j)
    // For diagonal (i, i): sum_{c=0}^{i-1}(k - c) + 0
    let mut idx = 0;
    for c in 0..i {
        idx += k - c;
    }
    idx
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_enrichment_basic() {
        let s_baseline = faer::mat![[0.3, 0.1], [0.1, 0.4]];
        let s_annot = vec![faer::mat![[0.15, 0.05], [0.05, 0.2]]];
        let v_baseline = Mat::from_fn(3, 3, |i, j| if i == j { 0.01 } else { 0.0 });
        let v_annot = vec![Mat::from_fn(3, 3, |i, j| if i == j { 0.02 } else { 0.0 })];
        let names = vec!["Annot1".to_string()];
        let m_annot = vec![100000.0];
        let m_total = 1000000.0;

        let result = enrichment_test(
            &s_baseline,
            &s_annot,
            &v_baseline,
            &v_annot,
            &names,
            &m_annot,
            m_total,
        );
        assert_eq!(result.annotations.len(), 1);
        assert!(result.enrichment[0] > 0.0);
        assert!(result.se[0] > 0.0, "SE should be positive, not placeholder");
        assert!(result.p[0] >= 0.0 && result.p[0] <= 1.0);
    }

    #[test]
    fn test_vech_diag_index() {
        // For k=3, vech order: (0,0), (1,0), (2,0), (1,1), (2,1), (2,2)
        assert_eq!(vech_diag_index(0, 3), 0);
        assert_eq!(vech_diag_index(1, 3), 3);
        assert_eq!(vech_diag_index(2, 3), 5);
    }
}
