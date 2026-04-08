use crate::error::SemError;

/// Operators in lavaan model syntax.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    /// =~ factor loading
    Loading,
    /// ~~ covariance/variance
    Covariance,
    /// ~ regression
    Regression,
    /// := defined parameter
    Defined,
}

impl std::fmt::Display for Op {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Op::Loading => write!(f, "=~"),
            Op::Covariance => write!(f, "~~"),
            Op::Regression => write!(f, "~"),
            Op::Defined => write!(f, ":="),
        }
    }
}

/// A single row in the parameter table.
#[derive(Debug, Clone)]
pub struct ParRow {
    pub lhs: String,
    pub op: Op,
    pub rhs: String,
    /// 0 = fixed, 1..n = free parameter index
    pub free: usize,
    /// Fixed value or start value
    pub value: f64,
    /// Optional label
    pub label: Option<String>,
    /// Lower bound constraint
    pub lower_bound: Option<f64>,
    /// Expression for defined params (e.g., "a*b")
    pub expression: Option<String>,
}

/// Parsed parameter table from lavaan syntax.
#[derive(Debug, Clone)]
pub struct ParTable {
    pub rows: Vec<ParRow>,
}

/// Parse lavaan model syntax into a parameter table.
///
/// Supported syntax:
/// - `F1 =~ V1 + V2 + V3` (factor loadings)
/// - `F1 ~~ 1*F1` (fix variance)
/// - `V1 ~~ V2` (covariance)
/// - `F1 ~ SNP` (regression)
/// - `ind := a*b` (defined parameter)
/// - `NA*V1` frees first loading
/// - `label*V1` assigns label
/// - `0.5*V1` fixes/starts at value
/// - `label > 0.0001` (constraint, parsed as lower bound)
pub fn parse_model(model: &str, std_lv: bool) -> Result<ParTable, SemError> {
    let mut rows = Vec::new();
    let mut free_counter = 0usize;
    let mut constraints: Vec<(String, f64)> = Vec::new(); // label -> lower bound

    // Split on newlines and semicolons
    let lines: Vec<&str> = model
        .split(['\n', ';'])
        .map(|s| s.trim())
        .filter(|s| !s.is_empty() && !s.starts_with('#'))
        .collect();

    for line in &lines {
        // Check for constraint: label > value
        if let Some(idx) = line.find('>') {
            let label = line[..idx].trim().to_string();
            let val: f64 = line[idx + 1..]
                .trim()
                .parse()
                .map_err(|_| SemError::SyntaxError(format!("invalid constraint: {line}")))?;
            constraints.push((label, val));
            continue;
        }

        // Detect operator
        let (lhs, op, rhs_str) = if let Some(idx) = line.find(":=") {
            let l = line[..idx].trim();
            let r = line[idx + 2..].trim();
            (l, Op::Defined, r)
        } else if let Some(idx) = line.find("=~") {
            let l = line[..idx].trim();
            let r = line[idx + 2..].trim();
            (l, Op::Loading, r)
        } else if let Some(idx) = line.find("~~") {
            let l = line[..idx].trim();
            let r = line[idx + 2..].trim();
            (l, Op::Covariance, r)
        } else if let Some(idx) = line.find('~') {
            let l = line[..idx].trim();
            let r = line[idx + 1..].trim();
            (l, Op::Regression, r)
        } else {
            return Err(SemError::SyntaxError(format!(
                "no operator found in: {line}"
            )));
        };

        let lhs_str = lhs.to_string();

        if op == Op::Defined {
            // Defined parameter: ind := a*b
            free_counter += 1;
            rows.push(ParRow {
                lhs: lhs_str,
                op,
                rhs: String::new(),
                free: free_counter,
                value: 0.0,
                label: None,
                expression: Some(rhs_str.to_string()),
                lower_bound: None,
            });
            continue;
        }

        // Parse RHS terms separated by +
        let terms: Vec<&str> = rhs_str.split('+').map(|s| s.trim()).collect();
        let mut is_first_loading = true;

        for term in &terms {
            let (var_name, fixed_val, label, is_na) = parse_term(term);

            let is_free;
            let value;

            if let Some(fv) = fixed_val {
                // Explicit fixed value: e.g., 1*V1
                is_free = false;
                value = fv;
            } else if is_na {
                // NA* prefix: free this parameter
                is_free = true;
                value = 0.0;
            } else if op == Op::Loading && is_first_loading && !std_lv {
                // Default: first loading per factor is fixed to 1 (unless std.lv)
                is_free = false;
                value = 1.0;
            } else {
                // Free parameter
                is_free = true;
                value = 0.0;
            }

            let free_idx = if is_free {
                free_counter += 1;
                free_counter
            } else {
                0
            };

            rows.push(ParRow {
                lhs: lhs_str.clone(),
                op,
                rhs: var_name.to_string(),
                free: free_idx,
                value,
                label: label.map(String::from),
                expression: None,
                lower_bound: None,
            });

            if op == Op::Loading {
                is_first_loading = false;
            }
        }
    }

    // Apply constraints: find rows with matching labels and set lower bounds
    for (label, lb) in &constraints {
        for row in &mut rows {
            if row.label.as_deref() == Some(label) {
                row.lower_bound = Some(*lb);
            }
        }
    }

    Ok(ParTable { rows })
}

/// Parse a single term like "a*V1", "1*V1", "NA*V1", or just "V1".
/// Returns (variable_name, fixed_value, label, is_na).
fn parse_term(term: &str) -> (&str, Option<f64>, Option<&str>, bool) {
    if let Some(idx) = term.find('*') {
        let prefix = term[..idx].trim();
        let var = term[idx + 1..].trim();

        if prefix.eq_ignore_ascii_case("NA") {
            return (var, None, None, true);
        }

        // Try parsing as number
        if let Ok(val) = prefix.parse::<f64>() {
            return (var, Some(val), None, false);
        }

        // It's a label
        (var, None, Some(prefix), false)
    } else {
        (term.trim(), None, None, false)
    }
}

impl ParTable {
    /// Get all unique observed variable names (appear in rhs of =~ or lhs of ~~ without being latent).
    pub fn observed_vars(&self) -> Vec<String> {
        let latent = self.latent_vars();
        let mut observed = Vec::new();
        let mut seen = std::collections::HashSet::new();

        for row in &self.rows {
            if row.op == Op::Loading && !latent.contains(&row.rhs) && seen.insert(row.rhs.clone()) {
                observed.push(row.rhs.clone());
            }
        }

        // Also include regression targets that aren't latent
        for row in &self.rows {
            if row.op == Op::Regression
                && !latent.contains(&row.rhs)
                && seen.insert(row.rhs.clone())
            {
                observed.push(row.rhs.clone());
            }
        }

        observed
    }

    /// Get all latent variable names (appear on lhs of =~).
    pub fn latent_vars(&self) -> Vec<String> {
        let mut latent = Vec::new();
        let mut seen = std::collections::HashSet::new();
        for row in &self.rows {
            if row.op == Op::Loading && seen.insert(row.lhs.clone()) {
                latent.push(row.lhs.clone());
            }
        }
        latent
    }

    /// Count free parameters.
    pub fn n_free(&self) -> usize {
        self.rows.iter().filter(|r| r.free > 0).count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_cfa() {
        let model = "F1 =~ V1 + V2 + V3\nF1 ~~ 1*F1";
        let pt = parse_model(model, false).unwrap();
        // 3 loadings + 1 variance = 4 rows
        assert_eq!(pt.rows.len(), 4);
        // First loading fixed to 1
        assert_eq!(pt.rows[0].free, 0);
        assert_eq!(pt.rows[0].value, 1.0);
        // Second and third free
        assert!(pt.rows[1].free > 0);
        assert!(pt.rows[2].free > 0);
        // Variance fixed to 1
        assert_eq!(pt.rows[3].free, 0);
        assert_eq!(pt.rows[3].value, 1.0);
    }

    #[test]
    fn test_parse_na_prefix() {
        let model = "F1 =~ NA*V1 + V2 + V3\nF1 ~~ 1*F1";
        let pt = parse_model(model, false).unwrap();
        // First loading should be free (NA* prefix)
        assert!(pt.rows[0].free > 0);
    }

    #[test]
    fn test_parse_regression() {
        let model = "F1 =~ V1 + V2\nF1 ~ SNP\nF1 ~~ 1*F1";
        let pt = parse_model(model, false).unwrap();
        let reg_rows: Vec<_> = pt.rows.iter().filter(|r| r.op == Op::Regression).collect();
        assert_eq!(reg_rows.len(), 1);
        assert_eq!(reg_rows[0].lhs, "F1");
        assert_eq!(reg_rows[0].rhs, "SNP");
    }

    #[test]
    fn test_parse_labels() {
        let model = "F1 =~ a*V1 + b*V2";
        let pt = parse_model(model, false).unwrap();
        assert_eq!(pt.rows[0].label.as_deref(), Some("a"));
        assert_eq!(pt.rows[1].label.as_deref(), Some("b"));
    }

    #[test]
    fn test_parse_defined() {
        let model = "F1 =~ a*V1 + V2\nind := a*b";
        let pt = parse_model(model, false).unwrap();
        let defined: Vec<_> = pt.rows.iter().filter(|r| r.op == Op::Defined).collect();
        assert_eq!(defined.len(), 1);
        assert_eq!(defined[0].expression.as_deref(), Some("a*b"));
    }

    #[test]
    fn test_parse_constraint() {
        let model = "F1 =~ V1 + V2\nV1 ~~ res1*V1\nres1 > 0.0001";
        let pt = parse_model(model, false).unwrap();
        let constrained: Vec<_> = pt.rows.iter().filter(|r| r.lower_bound.is_some()).collect();
        assert_eq!(constrained.len(), 1);
        assert_eq!(constrained[0].lower_bound, Some(0.0001));
    }

    #[test]
    fn test_latent_and_observed() {
        let model = "F1 =~ V1 + V2 + V3\nF2 =~ V3 + V4 + V5";
        let pt = parse_model(model, false).unwrap();
        assert_eq!(pt.latent_vars(), vec!["F1", "F2"]);
        assert_eq!(pt.observed_vars(), vec!["V1", "V2", "V3", "V4", "V5"]);
    }

    #[test]
    fn test_std_lv() {
        let model = "F1 =~ V1 + V2 + V3";
        let pt = parse_model(model, true).unwrap();
        // With std.lv, first loading should be free
        assert!(pt.rows[0].free > 0);
    }
}
