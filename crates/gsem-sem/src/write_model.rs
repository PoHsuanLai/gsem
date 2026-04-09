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
    mustload: bool,
    common: bool,
) -> String {
    let n_pheno = loadings.nrows();
    let n_factors = loadings.ncols();
    let mut lines = Vec::new();

    // Build per-factor indicator lists
    let mut factor_indicators: Vec<Vec<usize>> = vec![Vec::new(); n_factors];
    for f in 0..n_factors {
        for p in 0..n_pheno {
            if loadings[(p, f)].abs() > cutoff {
                factor_indicators[f].push(p);
            }
        }
    }

    // mustload: ensure every phenotype loads on at least one factor
    if mustload {
        for p in 0..n_pheno {
            let has_any = factor_indicators.iter().any(|inds| inds.contains(&p));
            if !has_any {
                // Find factor with highest absolute loading for this phenotype
                let best_f = (0..n_factors)
                    .max_by(|&a, &b| {
                        loadings[(p, a)]
                            .abs()
                            .partial_cmp(&loadings[(p, b)].abs())
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .unwrap_or(0);
                factor_indicators[best_f].push(p);
            }
        }
    }

    // Factor definitions
    for (f, indicators) in factor_indicators.iter().enumerate() {
        if indicators.is_empty() {
            continue;
        }
        let factor_name = format!("F{}", f + 1);
        let ind_names: Vec<&str> = indicators.iter().map(|&p| names[p].as_str()).collect();

        let terms: Vec<String> = std::iter::once(format!("NA*{}", ind_names[0]))
            .chain(ind_names[1..].iter().map(|s| s.to_string()))
            .collect();

        lines.push(format!("{} =~ {}", factor_name, terms.join(" + ")));
        lines.push(format!("{} ~~ 1*{}", factor_name, factor_name));
    }

    // Common factor (non-bifactor): all phenotypes load, orthogonal to specifics
    if common && !bifactor && n_factors > 0 {
        let all_phenos: Vec<&str> = names.iter().map(|s| s.as_str()).collect();
        if !all_phenos.is_empty() {
            let terms: Vec<String> = std::iter::once(format!("NA*{}", all_phenos[0]))
                .chain(all_phenos[1..].iter().map(|s| s.to_string()))
                .collect();
            lines.push(format!("Common_F =~ {}", terms.join(" + ")));
            lines.push("Common_F ~~ 1*Common_F".to_string());

            // Common factor orthogonal to specific factors
            for f in 0..n_factors {
                if !factor_indicators[f].is_empty() {
                    lines.push(format!("Common_F ~~ 0*F{}", f + 1));
                }
            }
        }
    }

    // Bifactor: add common factor and zero ALL cross-factor covariances
    if bifactor && n_factors > 1 {
        let all_loaded: Vec<&str> = (0..n_pheno)
            .filter(|&p| factor_indicators.iter().any(|inds| inds.contains(&p)))
            .map(|p| names[p].as_str())
            .collect();

        if !all_loaded.is_empty() {
            let terms: Vec<String> = std::iter::once(format!("NA*{}", all_loaded[0]))
                .chain(all_loaded[1..].iter().map(|s| s.to_string()))
                .collect();
            lines.push(format!("Common_F =~ {}", terms.join(" + ")));
            lines.push("Common_F ~~ 1*Common_F".to_string());

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
            let has_loading = factor_indicators.iter().any(|inds| inds.contains(&p));
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
        let loadings = faer::mat![[0.8], [0.7], [0.6],];
        let names = vec!["V1".to_string(), "V2".to_string(), "V3".to_string()];
        let model = write_model(&loadings, &names, 0.3, false, false, false, false);
        assert!(model.contains("F1 =~"));
        assert!(model.contains("NA*V1"));
        assert!(model.contains("V2"));
        assert!(model.contains("V3"));
    }

    #[test]
    fn test_write_model_cutoff() {
        let loadings = faer::mat![[0.8], [0.1], [0.6],];
        let names = vec!["V1".to_string(), "V2".to_string(), "V3".to_string()];
        let model = write_model(&loadings, &names, 0.3, false, false, false, false);
        // V2 should NOT be included (loading 0.1 < cutoff 0.3)
        assert!(!model.contains("+ V2"));
    }

    #[test]
    fn test_write_model_with_residuals() {
        let loadings = faer::mat![[0.8], [0.7],];
        let names = vec!["V1".to_string(), "V2".to_string()];
        let model = write_model(&loadings, &names, 0.3, true, false, false, false);
        assert!(model.contains("> 0.0001"));
    }

    #[test]
    fn test_generate_label() {
        assert_eq!(generate_label(0), "aaaa");
        assert_eq!(generate_label(1), "aaab");
        assert_eq!(generate_label(26), "aaba");
    }
}
