//! Simple complementary filters the GUI offers next to Gblocks.
//!
//! None of these mutate: each returns a mask or a range that the caller can
//! preview, combine with a Gblocks mask, or hand to
//! [`Alignment::keep_columns`](tolviewer_core::Alignment::keep_columns).

use std::ops::Range;

use tolviewer_core::alphabet::is_gap;
use tolviewer_core::{Alignment, MISSING};

/// True for characters that carry no residue: gaps and the missing marker.
#[inline]
fn is_blank(c: u8) -> bool {
    is_gap(c) || c == MISSING
}

/// A keep-mask over columns: `true` where the column's gap fraction is at most
/// `max_gap_fraction`.
///
/// The mask has one entry per alignment column, so it can be passed straight
/// to [`Alignment::keep_columns`](tolviewer_core::Alignment::keep_columns).
/// Rows shorter than the alignment width count as gapped in the columns they
/// do not reach. `max_gap_fraction` is clamped to 0..=1, so 0.0 keeps only
/// gap-free columns and 1.0 keeps everything.
pub fn remove_gappy_columns(aln: &Alignment, max_gap_fraction: f32) -> Vec<bool> {
    let width = aln.width();
    let rows = aln.len();
    let limit = max_gap_fraction.clamp(0.0, 1.0);
    if rows == 0 {
        return vec![true; width];
    }
    let mut mask = vec![true; width];
    for (col, keep) in mask.iter_mut().enumerate() {
        let gaps = aln
            .sequences
            .iter()
            .filter(|s| s.residues.get(col).copied().is_none_or(is_blank))
            .count();
        *keep = gaps as f32 / rows as f32 <= limit;
    }
    mask
}

/// A keep-mask over rows: `true` where the sequence's gap fraction, measured
/// over the full alignment width, is at most `max_gap_fraction`.
///
/// The mask has one entry per sequence, in row order. `max_gap_fraction` is
/// clamped to 0..=1.
pub fn remove_gappy_sequences(aln: &Alignment, max_gap_fraction: f32) -> Vec<bool> {
    let width = aln.width();
    let limit = max_gap_fraction.clamp(0.0, 1.0);
    aln.sequences
        .iter()
        .map(|s| {
            if width == 0 {
                return true;
            }
            let residues = s.residues.iter().filter(|&&c| !is_blank(c)).count();
            let gaps = width - residues;
            gaps as f32 / width as f32 <= limit
        })
        .collect()
}

/// Trim ragged 5'/3' ends: the range from the first to just past the last
/// column whose occupancy (fraction of rows carrying a residue) is at least
/// `min_occupancy`.
///
/// Returns an empty range (`0..0`) when no column qualifies. Columns *inside*
/// the range are not examined — this only squares off the ends, which is what
/// distinguishes it from [`remove_gappy_columns`].
pub fn trim_ends(aln: &Alignment, min_occupancy: f32) -> Range<usize> {
    let width = aln.width();
    let rows = aln.len();
    if width == 0 || rows == 0 {
        return 0..0;
    }
    let threshold = min_occupancy.clamp(0.0, 1.0);
    let occupied = |col: usize| {
        let n = aln
            .sequences
            .iter()
            .filter(|s| s.residues.get(col).copied().is_some_and(|c| !is_blank(c)))
            .count();
        n as f32 / rows as f32 >= threshold
    };
    match (0..width).find(|&c| occupied(c)) {
        None => 0..0,
        Some(start) => {
            // `start` qualifies, so `rfind` cannot fail.
            let end = (start..width).rev().find(|&c| occupied(c)).unwrap_or(start);
            start..end + 1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tolviewer_core::Sequence;

    fn aln(rows: &[&str]) -> Alignment {
        Alignment::new(
            "test",
            rows.iter()
                .enumerate()
                .map(|(i, r)| Sequence::new(format!("s{i}"), r.as_bytes().to_vec()))
                .collect(),
        )
    }

    #[test]
    fn gappy_columns_by_fraction() {
        // gap counts per column: 0, 1, 2, 4
        let a = aln(&["ACG-", "ACG-", "AC--", "A---"]);
        assert_eq!(remove_gappy_columns(&a, 0.0), vec![true, false, false, false]);
        assert_eq!(remove_gappy_columns(&a, 0.25), vec![true, true, false, false]);
        assert_eq!(remove_gappy_columns(&a, 0.5), vec![true, true, true, false]);
        assert_eq!(remove_gappy_columns(&a, 1.0), vec![true; 4]);
    }

    #[test]
    fn gappy_columns_treats_short_rows_as_gapped() {
        let a = aln(&["ACGT", "AC"]);
        assert_eq!(remove_gappy_columns(&a, 0.0), vec![true, true, false, false]);
    }

    #[test]
    fn gappy_columns_does_not_mutate_and_clamps() {
        let a = aln(&["ACGT", "A-GT"]);
        let before = a.clone();
        assert_eq!(remove_gappy_columns(&a, -5.0), vec![true, false, true, true]);
        assert_eq!(remove_gappy_columns(&a, 9.0), vec![true; 4]);
        assert_eq!(a, before);
    }

    #[test]
    fn gappy_sequences_by_fraction() {
        // gap fractions: 0, 0.25, 0.5, 1.0
        let a = aln(&["ACGT", "ACG-", "AC--", "----"]);
        assert_eq!(remove_gappy_sequences(&a, 0.0), vec![true, false, false, false]);
        assert_eq!(remove_gappy_sequences(&a, 0.25), vec![true, true, false, false]);
        assert_eq!(remove_gappy_sequences(&a, 0.5), vec![true, true, true, false]);
        assert_eq!(remove_gappy_sequences(&a, 1.0), vec![true; 4]);
    }

    #[test]
    fn gappy_sequences_measures_against_the_full_width() {
        // The short row is 100% gap once padded to width 4.
        let a = aln(&["ACGT", ""]);
        assert_eq!(remove_gappy_sequences(&a, 0.9), vec![true, false]);
    }

    #[test]
    fn trim_ends_squares_off_ragged_ends() {
        // occupancy per column: 1/4, 2/4, 4/4, 4/4, 4/4, 2/4, 1/4
        let a = aln(&["ACGTACG", "-CGTAC-", "--GTA--", "--GTA--"]);
        assert_eq!(trim_ends(&a, 1.0), 2..5);
        assert_eq!(trim_ends(&a, 0.5), 1..6);
        assert_eq!(trim_ends(&a, 0.25), 0..7);
        assert_eq!(trim_ends(&a, 0.0), 0..7);
    }

    #[test]
    fn trim_ends_keeps_interior_gaps() {
        let a = aln(&["-AAAAA-", "-A---A-"]);
        // Column 3 is empty but sits inside the trimmed range.
        assert_eq!(trim_ends(&a, 1.0), 1..6);
    }

    #[test]
    fn trim_ends_on_hopeless_or_empty_input() {
        let a = aln(&["----", "----"]);
        assert_eq!(trim_ends(&a, 0.5), 0..0);
        assert_eq!(trim_ends(&aln(&[]), 0.5), 0..0);
        assert_eq!(trim_ends(&aln(&["", ""]), 0.5), 0..0);
        assert!(remove_gappy_columns(&aln(&[]), 0.5).is_empty());
        assert!(remove_gappy_sequences(&aln(&[]), 0.5).is_empty());
    }

    #[test]
    fn missing_marker_counts_as_a_gap_everywhere() {
        let a = aln(&["A?A", "AAA"]);
        assert_eq!(remove_gappy_columns(&a, 0.0), vec![true, false, true]);
        assert_eq!(remove_gappy_sequences(&a, 0.0), vec![false, true]);
    }
}
