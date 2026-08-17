//! Per-column and per-alignment statistics used by the viewer's quality track.

use crate::alignment::Alignment;
use crate::alphabet::{is_gap, Alphabet, GAP};

#[derive(Debug, Clone, Copy, Default)]
pub struct ColumnStats {
    /// Number of rows with a non-gap residue here.
    pub occupancy: u32,
    /// Number of rows in total.
    pub rows: u32,
    /// Most frequent non-gap, unambiguous residue (uppercase), if any.
    pub majority: Option<u8>,
    /// Count of the majority residue.
    pub majority_count: u32,
    /// Number of distinct non-gap, unambiguous residues.
    pub distinct: u8,
    /// Number of rows with an ambiguous residue (N, X, ?, IUPAC).
    pub ambiguous: u32,
}

impl ColumnStats {
    /// Fraction of rows with a residue here, 0..=1.
    pub fn occupancy_fraction(&self) -> f32 {
        if self.rows == 0 {
            0.0
        } else {
            self.occupancy as f32 / self.rows as f32
        }
    }

    /// Fraction of occupied rows carrying the majority residue, 0..=1.
    /// Columns with no residues score 0.
    pub fn identity(&self) -> f32 {
        if self.occupancy == 0 {
            0.0
        } else {
            self.majority_count as f32 / self.occupancy as f32
        }
    }

    /// True when every occupied row has the same unambiguous residue.
    pub fn is_conserved(&self) -> bool {
        self.occupancy > 0 && self.distinct == 1 && self.ambiguous == 0
    }

    /// A 0..=1 quality score combining occupancy and identity, used to draw the
    /// summary bar above the alignment.
    pub fn quality(&self) -> f32 {
        self.occupancy_fraction() * self.identity()
    }
}

/// Consensus sequence plus the per-column statistics it was derived from.
#[derive(Debug, Clone, Default)]
pub struct Consensus {
    pub residues: Vec<u8>,
    pub columns: Vec<ColumnStats>,
}

impl Consensus {
    /// Compute the consensus and column stats.
    ///
    /// A column becomes the majority residue when that residue reaches
    /// `threshold` of the occupied rows (0.5 gives a simple majority rule);
    /// otherwise it becomes `N`/`X`. Columns below `min_occupancy` become gaps.
    pub fn compute(
        aln: &Alignment,
        alphabet: Alphabet,
        threshold: f32,
        min_occupancy: f32,
    ) -> Consensus {
        let width = aln.width();
        let rows = aln.len() as u32;
        let mut columns = Vec::with_capacity(width);
        let mut residues = Vec::with_capacity(width);
        let unknown = if alphabet.is_nucleotide() { b'N' } else { b'X' };
        let core = alphabet.core_symbols();

        // Counts indexed by residue byte; 256 entries is smaller than any
        // hashing scheme and keeps the inner loop branch-free.
        let mut counts = [0u32; 256];
        for col in 0..width {
            counts.fill(0);
            let mut occupancy = 0u32;
            let mut ambiguous = 0u32;
            for c in aln.column(col) {
                if is_gap(c) {
                    continue;
                }
                occupancy += 1;
                let up = c.to_ascii_uppercase();
                if core.contains(&up) {
                    counts[up as usize] += 1;
                } else {
                    ambiguous += 1;
                }
            }
            let mut majority = None;
            let mut majority_count = 0u32;
            let mut distinct = 0u8;
            for &sym in core {
                let n = counts[sym as usize];
                if n > 0 {
                    distinct += 1;
                    if n > majority_count {
                        majority_count = n;
                        majority = Some(sym);
                    }
                }
            }
            let stats =
                ColumnStats { occupancy, rows, majority, majority_count, distinct, ambiguous };
            let occ = stats.occupancy_fraction();
            residues.push(if occ < min_occupancy {
                GAP
            } else if stats.identity() >= threshold {
                majority.unwrap_or(unknown)
            } else {
                unknown
            });
            columns.push(stats);
        }
        Consensus { residues, columns }
    }

    /// Fraction of columns that are fully conserved.
    pub fn conserved_fraction(&self) -> f32 {
        if self.columns.is_empty() {
            return 0.0;
        }
        let n = self.columns.iter().filter(|c| c.is_conserved()).count();
        n as f32 / self.columns.len() as f32
    }
}

/// Percent identity between two rows over the columns where both have a
/// residue. Returns `None` when they never overlap.
pub fn pairwise_identity(a: &[u8], b: &[u8]) -> Option<f32> {
    let mut same = 0usize;
    let mut compared = 0usize;
    for (&x, &y) in a.iter().zip(b) {
        if is_gap(x) || is_gap(y) {
            continue;
        }
        compared += 1;
        if x.eq_ignore_ascii_case(&y) {
            same += 1;
        }
    }
    (compared > 0).then(|| same as f32 / compared as f32)
}

/// Mean pairwise identity across all row pairs. On large alignments this is
/// O(n^2 * width); callers should sample instead when `aln.len()` is big.
pub fn mean_pairwise_identity(aln: &Alignment) -> Option<f32> {
    let n = aln.len();
    if n < 2 {
        return None;
    }
    let mut sum = 0.0f64;
    let mut pairs = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            if let Some(id) =
                pairwise_identity(&aln.sequences[i].residues, &aln.sequences[j].residues)
            {
                sum += id as f64;
                pairs += 1;
            }
        }
    }
    (pairs > 0).then(|| (sum / pairs as f64) as f32)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::Sequence;

    fn aln() -> Alignment {
        Alignment::new(
            "t",
            vec![
                Sequence::new("a", *b"ACGT-A"),
                Sequence::new("b", *b"ACGT-C"),
                Sequence::new("c", *b"AC-T-G"),
            ],
        )
    }

    #[test]
    fn consensus_of_a_conserved_column() {
        let c = Consensus::compute(&aln(), Alphabet::Dna, 0.5, 0.5);
        assert_eq!(c.residues[0], b'A');
        assert!(c.columns[0].is_conserved());
        assert_eq!(c.columns[0].occupancy, 3);
    }

    #[test]
    fn column_with_a_gap_still_reaches_consensus() {
        let c = Consensus::compute(&aln(), Alphabet::Dna, 0.5, 0.5);
        // column 2 is G, G, gap -> occupancy 2/3, identity 1.0
        assert_eq!(c.residues[2], b'G');
        assert!((c.columns[2].occupancy_fraction() - 2.0 / 3.0).abs() < 1e-6);
    }

    #[test]
    fn all_gap_column_is_a_gap_in_consensus() {
        let c = Consensus::compute(&aln(), Alphabet::Dna, 0.5, 0.5);
        assert_eq!(c.residues[4], GAP);
        assert_eq!(c.columns[4].occupancy, 0);
        assert_eq!(c.columns[4].quality(), 0.0);
    }

    #[test]
    fn no_majority_becomes_n() {
        let c = Consensus::compute(&aln(), Alphabet::Dna, 0.9, 0.5);
        // column 5 is A, C, G -> no residue reaches 90%
        assert_eq!(c.residues[5], b'N');
        assert_eq!(c.columns[5].distinct, 3);
    }

    #[test]
    fn identity_ignores_gap_columns() {
        assert_eq!(pairwise_identity(b"ACGT", b"ACGA"), Some(0.75));
        assert_eq!(pairwise_identity(b"AC--", b"AC-T"), Some(1.0));
        assert_eq!(pairwise_identity(b"----", b"ACGT"), None);
    }

    #[test]
    fn mean_identity_over_three_rows() {
        let m = mean_pairwise_identity(&aln()).unwrap();
        assert!(m > 0.5 && m < 1.0, "unexpected mean identity {m}");
    }
}
