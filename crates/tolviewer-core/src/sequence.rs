//! A single named sequence, aligned or not.

use crate::alphabet::{is_gap, normalize_residue, Alphabet, GAP};

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Sequence {
    /// Short identifier, e.g. the first word of a FASTA header.
    pub id: String,
    /// Free-text remainder of the header, if any.
    pub description: String,
    /// Residues, gaps included. ASCII, case preserved.
    pub residues: Vec<u8>,
    /// Per-residue Phred scores from FASTQ / trace files, if known. When
    /// present this is the same length as `residues`, with 0 at gaps.
    pub quality: Option<Vec<u8>>,
    /// User flag: excluded rows are drawn dimmed and skipped by exports that
    /// ask for it.
    pub hidden: bool,
}

impl Sequence {
    pub fn new(id: impl Into<String>, residues: impl Into<Vec<u8>>) -> Self {
        let residues: Vec<u8> = residues.into();
        Sequence {
            id: id.into(),
            description: String::new(),
            residues: residues.into_iter().map(normalize_residue).collect(),
            quality: None,
            hidden: false,
        }
    }

    /// Length including gaps (i.e. the number of alignment columns spanned).
    pub fn len(&self) -> usize {
        self.residues.len()
    }

    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }

    /// Length excluding gaps.
    pub fn ungapped_len(&self) -> usize {
        self.residues.iter().filter(|&&c| !is_gap(c)).count()
    }

    /// The header line as it would appear in FASTA, without the `>`.
    pub fn header(&self) -> String {
        if self.description.is_empty() {
            self.id.clone()
        } else {
            format!("{} {}", self.id, self.description)
        }
    }

    pub fn set_header(&mut self, header: &str) {
        let header = header.trim();
        match header.split_once(char::is_whitespace) {
            Some((id, rest)) => {
                self.id = id.to_string();
                self.description = rest.trim().to_string();
            }
            None => {
                self.id = header.to_string();
                self.description.clear();
            }
        }
    }

    /// Residues with all gaps removed.
    pub fn ungapped(&self) -> Vec<u8> {
        self.residues.iter().copied().filter(|&c| !is_gap(c)).collect()
    }

    /// Map an alignment column to the 0-based residue index, or `None` if the
    /// column holds a gap. Linear scan; the GUI caches this per row.
    pub fn residue_index_at(&self, column: usize) -> Option<usize> {
        if column >= self.residues.len() || is_gap(self.residues[column]) {
            return None;
        }
        Some(self.residues[..column].iter().filter(|&&c| !is_gap(c)).count())
    }

    /// Reverse complement in place (no-op for protein). Quality is reversed too.
    pub fn reverse_complement(&mut self, alphabet: Alphabet) {
        self.residues.reverse();
        if alphabet.is_nucleotide() {
            for c in &mut self.residues {
                *c = alphabet.complement(*c);
            }
        }
        if let Some(q) = &mut self.quality {
            q.reverse();
        }
    }

    /// Pad with gaps on the right so the row spans `width` columns.
    pub fn pad_to(&mut self, width: usize) {
        if self.residues.len() < width {
            let extra = width - self.residues.len();
            self.residues.resize(width, GAP);
            if let Some(q) = &mut self.quality {
                q.resize(q.len() + extra, 0);
            }
        }
    }

    /// Mean quality over non-gap positions, if quality is known.
    pub fn mean_quality(&self) -> Option<f32> {
        let q = self.quality.as_ref()?;
        let mut sum = 0u64;
        let mut n = 0u64;
        for (i, &score) in q.iter().enumerate() {
            if self.residues.get(i).is_some_and(|&c| !is_gap(c)) {
                sum += score as u64;
                n += 1;
            }
        }
        (n > 0).then(|| sum as f32 / n as f32)
    }

    /// Fraction of non-gap residues that are ambiguous (N, X, ?, IUPAC codes).
    pub fn ambiguity_fraction(&self, alphabet: Alphabet) -> f32 {
        let mut amb = 0usize;
        let mut total = 0usize;
        for &c in &self.residues {
            if is_gap(c) {
                continue;
            }
            total += 1;
            if alphabet.is_ambiguous(c) {
                amb += 1;
            }
        }
        if total == 0 {
            0.0
        } else {
            amb as f32 / total as f32
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_round_trip() {
        let mut s = Sequence::new("x", *b"ACGT");
        s.set_header("Homo_sapiens COI partial cds");
        assert_eq!(s.id, "Homo_sapiens");
        assert_eq!(s.description, "COI partial cds");
        assert_eq!(s.header(), "Homo_sapiens COI partial cds");
    }

    #[test]
    fn alternative_gap_chars_are_normalised() {
        let s = Sequence::new("x", *b"AC.GT~A");
        assert_eq!(s.residues, b"AC-GT-A");
    }

    #[test]
    fn residue_index_skips_gaps() {
        let s = Sequence::new("x", *b"A-CG");
        assert_eq!(s.residue_index_at(0), Some(0));
        assert_eq!(s.residue_index_at(1), None);
        assert_eq!(s.residue_index_at(2), Some(1));
        assert_eq!(s.residue_index_at(3), Some(2));
    }

    #[test]
    fn revcomp_of_gapped_dna() {
        let mut s = Sequence::new("x", *b"AACG-T");
        s.reverse_complement(Alphabet::Dna);
        assert_eq!(s.residues, b"A-CGTT");
    }
}
