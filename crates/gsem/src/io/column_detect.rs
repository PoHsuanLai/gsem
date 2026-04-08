use std::collections::HashMap;
use std::sync::LazyLock;

/// Inverted alias map: uppercase alias → canonical name.
/// Built once, shared across all calls to `detect_columns`.
static ALIAS_TO_CANONICAL: LazyLock<HashMap<&'static str, &'static str>> = LazyLock::new(|| {
    let entries: &[(&str, &[&str])] = &[
        (
            "SNP",
            &[
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
        ),
        (
            "A1",
            &[
                "A1",
                "ALLELE1",
                "EFFECT_ALLELE",
                "INC_ALLELE",
                "REFERENCE_ALLELE",
                "EA",
                "REF",
            ],
        ),
        (
            "A2",
            &[
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
        ),
        (
            "effect",
            &[
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
        ),
        ("INFO", &["INFO", "IMPINFO"]),
        (
            "P",
            &[
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
        ),
        (
            "N",
            &[
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
        ),
        (
            "MAF",
            &[
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
        ),
        (
            "Z",
            &[
                "Z",
                "ZSCORE",
                "Z-SCORE",
                "ZSTATISTIC",
                "ZSTAT",
                "Z-STATISTIC",
            ],
        ),
        (
            "SE",
            &[
                "STDERR",
                "SE",
                "STDERRLOGOR",
                "SEBETA",
                "STANDARDERROR",
                "STANDARD_ERROR",
            ],
        ),
        ("DIRECTION", &["DIRECTION", "DIREC", "DIRE", "SIGN"]),
        ("NEFFDIV2", &["NEFFDIV2"]),
        ("NEFF_HALF", &["NEFF_HALF"]),
    ];

    let mut map = HashMap::new();
    for &(canonical, aliases) in entries {
        for &alias in aliases {
            map.insert(alias, canonical);
        }
    }
    map
});

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
    let mut columns = HashMap::new();

    for (idx, header) in headers.iter().enumerate() {
        let upper = header.to_uppercase();
        if let Some(&canonical) = ALIAS_TO_CANONICAL.get(upper.as_str()) {
            columns.entry(canonical.to_string()).or_insert(idx);
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
