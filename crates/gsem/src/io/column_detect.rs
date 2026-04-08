use std::collections::HashMap;

/// Canonical column names and their known aliases.
/// Ported from GenomicSEM R's `.get_renamed_colnames()` in utils.R.
fn build_alias_map() -> HashMap<&'static str, Vec<&'static str>> {
    let mut m = HashMap::new();
    m.insert(
        "SNP",
        vec![
            "SNP",
            "SNPID",
            "RSID",
            "RS_NUMBER",
            "RS_NUMBERS",
            "MARKERNAME",
            "ID",
            "PREDICTOR",
            "SNP_ID",
            "VARIANTID",
            "VARIANT_ID",
            "RSIDS",
            "RS_ID",
        ],
    );
    m.insert(
        "A1",
        vec![
            "A1",
            "ALLELE1",
            "EFFECT_ALLELE",
            "INC_ALLELE",
            "REFERENCE_ALLELE",
            "EA",
            "REF",
        ],
    );
    m.insert(
        "A2",
        vec![
            "A2",
            "ALLELE2",
            "ALLELE0",
            "OTHER_ALLELE",
            "NON_EFFECT_ALLELE",
            "DEC_ALLELE",
            "OA",
            "NEA",
            "ALT",
            "A0",
        ],
    );
    m.insert(
        "effect",
        vec![
            "OR",
            "B",
            "BETA",
            "LOG_ODDS",
            "EFFECTS",
            "EFFECT",
            "SIGNED_SUMSTAT",
            "EST",
            "BETA1",
            "LOGOR",
        ],
    );
    m.insert("INFO", vec!["INFO", "IMPINFO"]);
    m.insert(
        "P",
        vec![
            "P",
            "PVALUE",
            "PVAL",
            "P_VALUE",
            "P-VALUE",
            "P.VALUE",
            "P_VAL",
            "GC_PVALUE",
            "WALD_P",
        ],
    );
    m.insert(
        "N",
        vec![
            "N",
            "WEIGHT",
            "NCOMPLETESAMPLES",
            "TOTALSAMPLESIZE",
            "TOTALN",
            "TOTAL_N",
            "N_COMPLETE_SAMPLES",
            "SAMPLESIZE",
            "NEFF",
            "N_EFF",
            "N_EFFECTIVE",
            "SUMNEFF",
        ],
    );
    m.insert(
        "MAF",
        vec![
            "MAF",
            "CEUAF",
            "FREQ1",
            "EAF",
            "FREQ1.HAPMAP",
            "FREQALLELE1HAPMAPCEU",
            "FREQ.ALLELE1.HAPMAPCEU",
            "EFFECT_ALLELE_FREQ",
            "FREQ.A1",
            "A1FREQ",
            "ALLELEFREQ",
            "EFFECT_ALLELE_FREQUENCY",
        ],
    );
    m.insert(
        "Z",
        vec![
            "Z",
            "ZSCORE",
            "Z-SCORE",
            "ZSTATISTIC",
            "ZSTAT",
            "Z-STATISTIC",
        ],
    );
    m.insert(
        "SE",
        vec![
            "STDERR",
            "SE",
            "STDERRLOGOR",
            "SEBETA",
            "STANDARDERROR",
            "STANDARD_ERROR",
        ],
    );
    m.insert("DIRECTION", vec!["DIRECTION", "DIREC", "DIRE", "SIGN"]);
    // Special N columns that need doubling
    m.insert("NEFFDIV2", vec!["NEFFDIV2"]);
    m.insert("NEFF_HALF", vec!["NEFF_HALF"]);
    m
}

/// Result of column detection: maps canonical names to column indices.
#[derive(Debug, Clone)]
pub struct DetectedColumns {
    /// Maps canonical name (e.g., "SNP", "P", "effect") to column index
    pub columns: HashMap<String, usize>,
    /// The original header names
    pub headers: Vec<String>,
}

impl DetectedColumns {
    /// Get column index for a canonical name, if detected.
    pub fn get(&self, canonical: &str) -> Option<usize> {
        self.columns.get(canonical).copied()
    }

    /// Check if a canonical column was detected.
    pub fn has(&self, canonical: &str) -> bool {
        self.columns.contains_key(canonical)
    }
}

/// Detect canonical column names from file headers.
///
/// Case-insensitive matching against known aliases. Returns a mapping
/// of canonical names to column indices.
pub fn detect_columns(headers: &[String]) -> DetectedColumns {
    let alias_map = build_alias_map();
    let mut columns = HashMap::new();

    // Uppercase all headers for case-insensitive matching
    let upper_headers: Vec<String> = headers.iter().map(|h| h.to_uppercase()).collect();

    for (canonical, aliases) in &alias_map {
        for (idx, header) in upper_headers.iter().enumerate() {
            if aliases.contains(&header.as_str()) {
                columns.insert(canonical.to_string(), idx);
                break;
            }
        }
    }

    DetectedColumns {
        columns,
        headers: headers.to_vec(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_detection() {
        let headers: Vec<String> = vec!["SNP", "A1", "A2", "BETA", "P", "N"]
            .into_iter()
            .map(String::from)
            .collect();
        let detected = detect_columns(&headers);
        assert_eq!(detected.get("SNP"), Some(0));
        assert_eq!(detected.get("A1"), Some(1));
        assert_eq!(detected.get("A2"), Some(2));
        assert_eq!(detected.get("effect"), Some(3));
        assert_eq!(detected.get("P"), Some(4));
        assert_eq!(detected.get("N"), Some(5));
    }

    #[test]
    fn test_case_insensitive() {
        let headers: Vec<String> =
            vec!["rsid", "allele1", "allele2", "beta", "pvalue", "samplesize"]
                .into_iter()
                .map(String::from)
                .collect();
        let detected = detect_columns(&headers);
        assert_eq!(detected.get("SNP"), Some(0));
        assert_eq!(detected.get("A1"), Some(1));
        assert_eq!(detected.get("A2"), Some(2));
        assert_eq!(detected.get("effect"), Some(3));
        assert_eq!(detected.get("P"), Some(4));
        assert_eq!(detected.get("N"), Some(5));
    }

    #[test]
    fn test_alternative_names() {
        let headers: Vec<String> = vec!["VARIANT_ID", "EA", "NEA", "OR", "PVAL", "NEFF", "EAF"]
            .into_iter()
            .map(String::from)
            .collect();
        let detected = detect_columns(&headers);
        assert_eq!(detected.get("SNP"), Some(0));
        assert_eq!(detected.get("A1"), Some(1));
        assert_eq!(detected.get("A2"), Some(2));
        assert_eq!(detected.get("effect"), Some(3));
        assert_eq!(detected.get("P"), Some(4));
        assert_eq!(detected.get("N"), Some(5));
        assert_eq!(detected.get("MAF"), Some(6));
    }

    #[test]
    fn test_missing_columns() {
        let headers: Vec<String> = vec!["SNP", "P"].into_iter().map(String::from).collect();
        let detected = detect_columns(&headers);
        assert!(detected.has("SNP"));
        assert!(detected.has("P"));
        assert!(!detected.has("A1"));
        assert!(!detected.has("effect"));
    }
}
