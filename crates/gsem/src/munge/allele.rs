/// Result of matching file alleles against reference alleles.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlleleMatch {
    /// Alleles match directly (A1==A1_ref, A2==A2_ref)
    Match,
    /// Alleles are flipped (A1==A2_ref, A2==A1_ref) -- effect sign should be negated
    Flipped,
    /// Alleles don't match at all -- SNP should be removed
    NoMatch,
}

/// Check if a pair of alleles is strand-ambiguous (A/T or C/G).
pub fn is_strand_ambiguous(a1: &str, a2: &str) -> bool {
    matches!(
        (a1, a2),
        ("A", "T") | ("T", "A") | ("C", "G") | ("G", "C")
    )
}

/// Get the complement of an allele (strand flip).
pub fn complement(allele: &str) -> &'static str {
    match allele {
        "A" => "T",
        "T" => "A",
        "C" => "G",
        "G" => "C",
        _ => "",
    }
}

/// Match file alleles against reference alleles.
///
/// Returns whether alleles match, are flipped, or don't match.
pub fn alleles_match(file_a1: &str, file_a2: &str, ref_a1: &str, ref_a2: &str) -> AlleleMatch {
    let fa1 = file_a1.to_uppercase();
    let fa2 = file_a2.to_uppercase();
    let ra1 = ref_a1.to_uppercase();
    let ra2 = ref_a2.to_uppercase();

    // Direct match
    if fa1 == ra1 && fa2 == ra2 {
        return AlleleMatch::Match;
    }

    // Flipped
    if fa1 == ra2 && fa2 == ra1 {
        return AlleleMatch::Flipped;
    }

    // Try strand complement
    let ca1 = complement(&fa1);
    let ca2 = complement(&fa2);

    if ca1 == ra1 && ca2 == ra2 {
        return AlleleMatch::Match;
    }

    if ca1 == ra2 && ca2 == ra1 {
        return AlleleMatch::Flipped;
    }

    AlleleMatch::NoMatch
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_direct_match() {
        assert_eq!(alleles_match("A", "G", "A", "G"), AlleleMatch::Match);
    }

    #[test]
    fn test_flipped() {
        assert_eq!(alleles_match("G", "A", "A", "G"), AlleleMatch::Flipped);
    }

    #[test]
    fn test_complement_match() {
        assert_eq!(alleles_match("T", "C", "A", "G"), AlleleMatch::Match);
    }

    #[test]
    fn test_complement_flipped() {
        assert_eq!(alleles_match("C", "T", "A", "G"), AlleleMatch::Flipped);
    }

    #[test]
    fn test_no_match() {
        assert_eq!(alleles_match("A", "C", "A", "G"), AlleleMatch::NoMatch);
    }

    #[test]
    fn test_strand_ambiguous() {
        assert!(is_strand_ambiguous("A", "T"));
        assert!(is_strand_ambiguous("C", "G"));
        assert!(!is_strand_ambiguous("A", "C"));
        assert!(!is_strand_ambiguous("G", "T"));
    }

    #[test]
    fn test_case_insensitive() {
        assert_eq!(alleles_match("a", "g", "A", "G"), AlleleMatch::Match);
    }
}
