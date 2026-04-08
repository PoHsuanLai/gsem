use faer::Mat;

/// Auto-generate lavaan model syntax from a factor loading matrix.
///
/// Port of GenomicSEM's `write.model()`.
///
/// For each factor, selects phenotypes with |loading| > cutoff.
/// Generates syntax like: F1 =~ NA*V1 + V2 + V3
pub fn write_model(
    loadings: &Mat<f64>,
    names: &[String],
    cutoff: f64,
    fix_resid: bool,
    bifactor: bool,
) -> String {
    let n_pheno = loadings.nrows();
    let n_factors = loadings.ncols();
    let mut lines = Vec::new();

    // Factor definitions
    for f in 0..n_factors {
        let factor_name = format!("F{}", f + 1);
        let mut indicators = Vec::new();

        for p in 0..n_pheno {
            if loadings[(p, f)].abs() > cutoff {
                indicators.push(names[p].clone());
            }
        }

        if indicators.is_empty() {
            continue;
        }

        // First indicator gets NA* prefix (free loading), factor variance fixed to 1
        let mut terms: Vec<String> = Vec::new();
        terms.push(format!("NA*{}", indicators[0]));
        for ind in &indicators[1..] {
            terms.push(ind.clone());
        }

        lines.push(format!("{} =~ {}", factor_name, terms.join(" + ")));
        lines.push(format!("{} ~~ 1*{}", factor_name, factor_name));
    }

    // Bifactor: add common factor and zero cross-factor covariances
    if bifactor && n_factors > 1 {
        let mut all_indicators = Vec::new();
        for p in 0..n_pheno {
            if (0..n_factors).any(|f| loadings[(p, f)].abs() > cutoff) {
                all_indicators.push(names[p].clone());
            }
        }

        if !all_indicators.is_empty() {
            let mut terms: Vec<String> = Vec::new();
            terms.push(format!("NA*{}", all_indicators[0]));
            for ind in &all_indicators[1..] {
                terms.push(ind.clone());
            }
            lines.push(format!("Common_F =~ {}", terms.join(" + ")));
            lines.push("Common_F ~~ 1*Common_F".to_string());

            // Zero covariances between common and specific factors
            for f in 0..n_factors {
                lines.push(format!("Common_F ~~ 0*F{}", f + 1));
            }
        }

        // Zero covariances between specific factors
        for f1 in 0..n_factors {
            for f2 in (f1 + 1)..n_factors {
                lines.push(format!("F{} ~~ 0*F{}", f1 + 1, f2 + 1));
            }
        }
    }

    // Residual variance constraints
    if fix_resid {
        let mut label_idx = 0;
        for p in 0..n_pheno {
            let has_loading = (0..n_factors).any(|f| loadings[(p, f)].abs() > cutoff);
            if has_loading {
                let label = generate_label(label_idx);
                lines.push(format!("{} ~~ {}*{}", names[p], label, names[p]));
                lines.push(format!("{} > 0.0001", label));
                label_idx += 1;
            }
        }
    }

    lines.join("\n")
}

/// Generate a unique 4-letter label from an index.
fn generate_label(idx: usize) -> String {
    let chars = b"abcdefghijklmnopqrstuvwxyz";
    let n = chars.len();
    let a = chars[idx % n] as char;
    let b = chars[(idx / n) % n] as char;
    let c = chars[(idx / (n * n)) % n] as char;
    let d = chars[(idx / (n * n * n)) % n] as char;
    format!("{d}{c}{b}{a}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_model_single_factor() {
        let loadings = faer::mat![
            [0.8],
            [0.7],
            [0.6],
        ];
        let names = vec!["V1".to_string(), "V2".to_string(), "V3".to_string()];
        let model = write_model(&loadings, &names, 0.3, false, false);
        assert!(model.contains("F1 =~"));
        assert!(model.contains("NA*V1"));
        assert!(model.contains("V2"));
        assert!(model.contains("V3"));
    }

    #[test]
    fn test_write_model_cutoff() {
        let loadings = faer::mat![
            [0.8],
            [0.1],
            [0.6],
        ];
        let names = vec!["V1".to_string(), "V2".to_string(), "V3".to_string()];
        let model = write_model(&loadings, &names, 0.3, false, false);
        // V2 should NOT be included (loading 0.1 < cutoff 0.3)
        assert!(!model.contains("+ V2"));
    }

    #[test]
    fn test_write_model_with_residuals() {
        let loadings = faer::mat![
            [0.8],
            [0.7],
        ];
        let names = vec!["V1".to_string(), "V2".to_string()];
        let model = write_model(&loadings, &names, 0.3, true, false);
        assert!(model.contains("> 0.0001"));
    }

    #[test]
    fn test_generate_label() {
        assert_eq!(generate_label(0), "aaaa");
        assert_eq!(generate_label(1), "aaab");
        assert_eq!(generate_label(26), "aaba");
    }
}
