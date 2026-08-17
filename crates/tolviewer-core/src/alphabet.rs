//! Residue alphabets and the canonical gap / missing characters.

/// The single canonical gap character. Everything in [`GAP_CHARS`] is
/// normalised to this on load.
pub const GAP: u8 = b'-';
/// Characters treated as gaps on input.
pub const GAP_CHARS: &[u8] = b"-.~";
/// Ambiguity character used when a residue is unknown.
pub const MISSING: u8 = b'?';

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Alphabet {
    Dna,
    Rna,
    Protein,
}

impl Alphabet {
    pub fn is_nucleotide(self) -> bool {
        matches!(self, Alphabet::Dna | Alphabet::Rna)
    }

    pub fn name(self) -> &'static str {
        match self {
            Alphabet::Dna => "DNA",
            Alphabet::Rna => "RNA",
            Alphabet::Protein => "protein",
        }
    }

    /// NEXUS `datatype=` token.
    pub fn nexus_datatype(self) -> &'static str {
        match self {
            Alphabet::Dna => "dna",
            Alphabet::Rna => "rna",
            Alphabet::Protein => "protein",
        }
    }

    /// Residue characters that are meaningful in this alphabet, uppercase,
    /// excluding gap and missing.
    pub fn symbols(self) -> &'static [u8] {
        match self {
            Alphabet::Dna => b"ACGTRYSWKMBDHVN",
            Alphabet::Rna => b"ACGURYSWKMBDHVN",
            Alphabet::Protein => b"ACDEFGHIKLMNPQRSTVWYBZXJUO",
        }
    }

    /// The unambiguous "core" residues, used for consensus and composition.
    pub fn core_symbols(self) -> &'static [u8] {
        match self {
            Alphabet::Dna => b"ACGT",
            Alphabet::Rna => b"ACGU",
            Alphabet::Protein => b"ACDEFGHIKLMNPQRSTVWY",
        }
    }

    pub fn is_valid(self, c: u8) -> bool {
        let c = c.to_ascii_uppercase();
        c == GAP || c == MISSING || self.symbols().contains(&c)
    }

    /// True for residues that carry no information (gap, `?`, `N`/`X`).
    pub fn is_ambiguous(self, c: u8) -> bool {
        let c = c.to_ascii_uppercase();
        if c == GAP || c == MISSING {
            return true;
        }
        match self {
            Alphabet::Dna | Alphabet::Rna => !b"ACGTU".contains(&c),
            Alphabet::Protein => !self.core_symbols().contains(&c),
        }
    }

    /// Guess the alphabet from residue content.
    ///
    /// Uses the usual heuristic: if at least 85% of the non-gap, non-missing
    /// residues are A/C/G/T/U/N, call it nucleotide, and pick RNA over DNA when
    /// U outnumbers T.
    pub fn guess(residues: impl IntoIterator<Item = u8>) -> Alphabet {
        let mut nuc = 0usize;
        let mut total = 0usize;
        let mut t = 0usize;
        let mut u = 0usize;
        for c in residues {
            let c = c.to_ascii_uppercase();
            if c == GAP || c == MISSING || c.is_ascii_whitespace() || c.is_ascii_digit() {
                continue;
            }
            total += 1;
            match c {
                b'T' => {
                    t += 1;
                    nuc += 1;
                }
                b'U' => {
                    u += 1;
                    nuc += 1;
                }
                b'A' | b'C' | b'G' | b'N' => nuc += 1,
                _ => {}
            }
        }
        if total == 0 || (nuc as f64) / (total as f64) >= 0.85 {
            if u > t {
                Alphabet::Rna
            } else {
                Alphabet::Dna
            }
        } else {
            Alphabet::Protein
        }
    }

    /// Complement a nucleotide, preserving case and IUPAC ambiguity.
    /// Returns the input unchanged for protein.
    pub fn complement(self, c: u8) -> u8 {
        if !self.is_nucleotide() {
            return c;
        }
        let upper = c.to_ascii_uppercase();
        let comp = match upper {
            b'A' => {
                if self == Alphabet::Rna {
                    b'U'
                } else {
                    b'T'
                }
            }
            b'T' | b'U' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'R' => b'Y',
            b'Y' => b'R',
            b'S' => b'S',
            b'W' => b'W',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            other => other,
        };
        if c.is_ascii_lowercase() {
            comp.to_ascii_lowercase()
        } else {
            comp
        }
    }
}

/// Normalise a raw input byte: map alternative gap characters to [`GAP`],
/// leave everything else alone.
#[inline]
pub fn normalize_residue(c: u8) -> u8 {
    if GAP_CHARS.contains(&c) {
        GAP
    } else {
        c
    }
}

#[inline]
pub fn is_gap(c: u8) -> bool {
    c == GAP || GAP_CHARS.contains(&c)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn guesses_dna() {
        assert_eq!(Alphabet::guess(*b"ACGTACGTNNAC"), Alphabet::Dna);
    }

    #[test]
    fn guesses_rna() {
        assert_eq!(Alphabet::guess(*b"ACGUACGUUUAC"), Alphabet::Rna);
    }

    #[test]
    fn guesses_protein() {
        assert_eq!(Alphabet::guess(*b"MKVLWIPQRSTY"), Alphabet::Protein);
    }

    #[test]
    fn gaps_do_not_sway_the_guess() {
        assert_eq!(Alphabet::guess(*b"---ACGT---ACGT"), Alphabet::Dna);
    }

    #[test]
    fn complement_preserves_case() {
        assert_eq!(Alphabet::Dna.complement(b'a'), b't');
        assert_eq!(Alphabet::Dna.complement(b'A'), b'T');
        assert_eq!(Alphabet::Rna.complement(b'A'), b'U');
        assert_eq!(Alphabet::Protein.complement(b'K'), b'K');
    }
}
