//! A BLOSUM62-derived "is this substitution scored positively?" table.
//!
//! Gblocks uses a similarity criterion (the program ships a Gonnet 120 matrix)
//! when deciding whether a *protein* column is conserved enough to anchor a
//! block flank: residues that are merely *similar* to the column's most common
//! residue are counted alongside identical ones.
//!
//! **This deliberately duplicates data that also lives in `tolviewer-align`'s
//! substitution matrices.** `tolviewer-clean` and `tolviewer-align` are
//! siblings in the dependency graph (see `docs/ARCHITECTURE.md`); both depend
//! only on `tolviewer-core`. Making cleaning depend on an alignment engine
//! just to ask a 20x20 boolean question would couple the two for no benefit
//! and would drag a large crate into the GUI's per-keystroke path. The table
//! below is a tiny, self-contained slice of that information, so it is kept
//! local on purpose.
//!
//! We substitute BLOSUM62 for Gonnet 120 because BLOSUM62 is the matrix used
//! everywhere else in TOLViewer; the two agree on the great majority of pairs
//! (both encode the classic exchange groups: I/L/M/V, F/W/Y, D/E/N/Q/K/R,
//! S/T/A, with C, G and P similar only to themselves).

/// The 20 standard amino acids, in the canonical BLOSUM row/column order.
pub(crate) const AMINO_ACIDS: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";

/// Bit `j` of `POSITIVE[i]` is set when BLOSUM62 scores
/// `AMINO_ACIDS[i]` against `AMINO_ACIDS[j]` **greater than zero**.
///
/// The table is symmetric, which [`tests::table_is_symmetric`] checks, and
/// every diagonal bit is set (BLOSUM62 has no non-positive diagonal entry).
const POSITIVE: [u32; 20] = [
    0x0000_8001, // A: A S
    0x0000_0822, // R: R Q K
    0x0000_810C, // N: N D H S
    0x0000_004C, // D: N D E
    0x0000_0010, // C: C
    0x0000_0862, // Q: R Q E K
    0x0000_0868, // E: D Q E K
    0x0000_0080, // G: G
    0x0004_0104, // H: N H Y
    0x0008_1600, // I: I L M V
    0x0008_1600, // L: I L M V
    0x0000_0862, // K: R Q E K
    0x0008_1600, // M: I L M V
    0x0006_2000, // F: F W Y
    0x0000_4000, // P: P
    0x0001_8005, // S: A N S T
    0x0001_8000, // T: S T
    0x0006_2000, // W: F W Y
    0x0006_2100, // Y: H F W Y
    0x0008_1600, // V: I L M V
];

/// Index of `c` (case-insensitive) among [`AMINO_ACIDS`], or `None` for gaps,
/// ambiguity codes and anything else outside the standard twenty.
#[inline]
pub(crate) fn amino_acid_index(c: u8) -> Option<usize> {
    let c = c.to_ascii_uppercase();
    AMINO_ACIDS.iter().position(|&a| a == c)
}

/// The set of amino acids similar to `index`, as a bit mask over
/// [`AMINO_ACIDS`] positions. Includes `index` itself.
#[inline]
pub(crate) fn similar_mask(index: usize) -> u32 {
    POSITIVE[index]
}

/// True when BLOSUM62 gives `a` against `b` a positive score, i.e. the two
/// residues belong to the same exchange group. Case-insensitive; false if
/// either residue is not one of the standard twenty.
///
/// The hot path in [`crate::gblocks`] uses [`similar_mask`] instead, which
/// answers the same question for a whole column at once; this pairwise form
/// exists so the table can be checked residue by residue.
#[cfg(test)]
#[inline]
pub(crate) fn is_similar(a: u8, b: u8) -> bool {
    match (amino_acid_index(a), amino_acid_index(b)) {
        (Some(i), Some(j)) => POSITIVE[i] & (1 << j) != 0,
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_is_symmetric() {
        for i in 0..20 {
            for j in 0..20 {
                assert_eq!(
                    POSITIVE[i] & (1 << j) != 0,
                    POSITIVE[j] & (1 << i) != 0,
                    "asymmetry at {}/{}",
                    AMINO_ACIDS[i] as char,
                    AMINO_ACIDS[j] as char
                );
            }
        }
    }

    #[test]
    fn every_residue_is_similar_to_itself() {
        for (i, &aa) in AMINO_ACIDS.iter().enumerate() {
            assert!(POSITIVE[i] & (1 << i) != 0, "{} not self-similar", aa as char);
        }
    }

    #[test]
    fn classic_exchange_groups() {
        // Hydrophobic aliphatics form a clique.
        for &a in b"ILMV" {
            for &b in b"ILMV" {
                assert!(is_similar(a, b), "{}/{} should be similar", a as char, b as char);
            }
        }
        // Aromatics.
        assert!(is_similar(b'F', b'Y') && is_similar(b'W', b'Y') && is_similar(b'F', b'W'));
        // Acidic / amide / basic neighbourhoods.
        assert!(is_similar(b'D', b'E') && is_similar(b'K', b'R') && is_similar(b'Q', b'E'));
        // Small polar.
        assert!(is_similar(b'S', b'T') && is_similar(b'A', b'S'));
    }

    #[test]
    fn loners_are_similar_only_to_themselves() {
        for &lone in b"CGP" {
            for &other in AMINO_ACIDS.iter() {
                assert_eq!(
                    is_similar(lone, other),
                    lone == other,
                    "{}/{}",
                    lone as char,
                    other as char
                );
            }
        }
    }

    #[test]
    fn known_negatives() {
        assert!(!is_similar(b'A', b'T')); // BLOSUM62 A/T is 0, not positive
        assert!(!is_similar(b'F', b'I'));
        assert!(!is_similar(b'R', b'E'));
        assert!(!is_similar(b'G', b'S'));
    }

    #[test]
    fn case_insensitive_and_unknown_residues() {
        assert!(is_similar(b'i', b'L'));
        assert!(!is_similar(b'X', b'X'));
        assert!(!is_similar(b'-', b'-'));
        assert_eq!(amino_acid_index(b'a'), Some(0));
        assert_eq!(amino_acid_index(b'?'), None);
    }
}
