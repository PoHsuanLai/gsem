pub mod column_detect;
pub mod gwas_reader;
pub mod hdl_reader;
pub mod ld_reader;
pub mod sumstats_reader;
pub mod twas_reader;
pub mod writer;

use anyhow::{Context, Result};

/// Case-insensitive column index lookup.
pub(crate) fn col_idx(headers: &[String], name: &str, context: &str) -> Result<usize> {
    let upper = name.to_uppercase();
    headers
        .iter()
        .position(|h| h.to_uppercase() == upper)
        .with_context(|| format!("{name} column not found in {context}"))
}

/// Parsed beta/se column info from a header row.
pub(crate) struct BetaSeColumns {
    pub trait_names: Vec<String>,
    pub beta_indices: Vec<usize>,
    pub se_indices: Vec<usize>,
}

/// Find beta.* and se.* columns in headers (case-insensitive prefix match).
pub(crate) fn parse_beta_se_columns(headers: &[String], context: &str) -> Result<BetaSeColumns> {
    let mut beta_cols: Vec<(usize, String)> = Vec::new();
    let mut se_cols: Vec<(usize, String)> = Vec::new();

    for (i, h) in headers.iter().enumerate() {
        let lower = h.to_lowercase();
        if let Some(name) = lower.strip_prefix("beta.") {
            beta_cols.push((i, name.to_string()));
        } else if let Some(name) = lower.strip_prefix("se.") {
            se_cols.push((i, name.to_string()));
        }
    }

    if beta_cols.is_empty() {
        anyhow::bail!("no beta.* columns found in {context}");
    }
    if beta_cols.len() != se_cols.len() {
        anyhow::bail!(
            "mismatched beta ({}) and se ({}) column counts in {context}",
            beta_cols.len(),
            se_cols.len()
        );
    }

    let trait_names: Vec<String> = beta_cols.iter().map(|(_, name)| name.clone()).collect();
    let beta_indices: Vec<usize> = beta_cols.iter().map(|(i, _)| *i).collect();
    let se_indices: Vec<usize> = se_cols.iter().map(|(i, _)| *i).collect();

    Ok(BetaSeColumns {
        trait_names,
        beta_indices,
        se_indices,
    })
}

/// Parse beta and se values from a row's fields. Returns None if any value is
/// missing or non-finite.
pub(crate) fn parse_beta_se_values(
    fields: &[&str],
    beta_indices: &[usize],
    se_indices: &[usize],
) -> Option<(Vec<f64>, Vec<f64>)> {
    let k = beta_indices.len();
    let mut beta = Vec::with_capacity(k);
    let mut se = Vec::with_capacity(k);

    for t in 0..k {
        match (
            fields[beta_indices[t]].parse::<f64>(),
            fields[se_indices[t]].parse::<f64>(),
        ) {
            (Ok(b), Ok(s)) if b.is_finite() && s.is_finite() => {
                beta.push(b);
                se.push(s);
            }
            _ => return None,
        }
    }

    Some((beta, se))
}
