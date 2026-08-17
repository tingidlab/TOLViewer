//! The Gblocks column-selection algorithm.
//!
//! Reimplemented from the published description, not from the original C
//! source:
//!
//! * Castresana J. (2000) "Selection of conserved blocks from multiple
//!   alignments for their use in phylogenetic analysis."
//!   *Mol Biol Evol* 17(4):540-552.
//! * The Gblocks 0.91b/1.0 documentation ("Gblocks Documentation", J.
//!   Castresana), which is what the shipped program actually implements.
//! * Talavera G. & Castresana J. (2007) *Syst Biol* 56(4):564-577, for the
//!   relaxed parameter set.

use std::ops::Range;

use tolviewer_core::alphabet::is_gap;
use tolviewer_core::{Alignment, Alphabet, Error, Result, GAP, MISSING};

use crate::similarity::{amino_acid_index, similar_mask, AMINO_ACIDS};

/// How gap-bearing columns are treated (Gblocks' `b5`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GapPolicy {
    /// No gap positions are allowed in the final selection: a column holding a
    /// single gap already disqualifies it.
    None,
    /// A column is a gap position only when half or more of the sequences have
    /// a gap there.
    Half,
    /// Gaps never disqualify a column; gap positions are not treated
    /// differently from any other position.
    All,
}

impl GapPolicy {
    /// The Gblocks command-line letter for this policy (`n`, `h`, `a`).
    pub fn letter(self) -> char {
        match self {
            GapPolicy::None => 'n',
            GapPolicy::Half => 'h',
            GapPolicy::All => 'a',
        }
    }

    /// Human-readable name for menus.
    pub fn name(self) -> &'static str {
        match self {
            GapPolicy::None => "none",
            GapPolicy::Half => "half",
            GapPolicy::All => "all",
        }
    }

    /// All policies, in increasing permissiveness, for menu construction.
    pub fn all() -> &'static [GapPolicy] {
        &[GapPolicy::None, GapPolicy::Half, GapPolicy::All]
    }

    /// True when a column with `gaps` gaps out of `rows` sequences is still
    /// eligible for selection under this policy.
    #[inline]
    fn allows(self, gaps: usize, rows: usize) -> bool {
        match self {
            GapPolicy::None => gaps == 0,
            // "Only positions where 50% or more of the sequences have a gap
            // are treated as a gap position" (Gblocks documentation).
            GapPolicy::Half => gaps * 2 < rows,
            GapPolicy::All => true,
        }
    }
}

/// The five Gblocks parameters, plus the protein similarity switch.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GblocksParams {
    /// b1: minimum number of sequences for a conserved position. Gblocks
    /// requires > n/2; default is `n / 2 + 1`.
    pub min_seqs_conserved: usize,
    /// b2: minimum number of sequences for a flank position. Must be >= b1;
    /// default is `ceil(n * 0.85)`.
    pub min_seqs_flank: usize,
    /// b3: maximum number of contiguous non-conserved positions. Default 8.
    pub max_contiguous_nonconserved: usize,
    /// b4: minimum length of a block. Default 10.
    pub min_block_length: usize,
    /// b5: allowed gap positions. Default [`GapPolicy::None`].
    pub gaps: GapPolicy,
    /// Treat similar residues (positive substitution score) as conserved for
    /// the flank test, as Gblocks does for protein.
    ///
    /// Ignored for nucleotide alignments, where there is no meaningful
    /// similarity relation between bases: [`gblocks`] falls back to identity.
    pub use_similarity: bool,
}

impl GblocksParams {
    /// Gblocks' own defaults for `n` sequences: b1 = `n / 2 + 1`,
    /// b2 = `ceil(n * 0.85)` (never below b1), b3 = 8, b4 = 10,
    /// b5 = [`GapPolicy::None`], similarity on.
    pub fn defaults(n_seqs: usize) -> Self {
        let b1 = n_seqs / 2 + 1;
        // ceil(n * 0.85) == ceil(17n / 20), in integer arithmetic.
        let b2 = (n_seqs * 17).div_ceil(20);
        GblocksParams {
            min_seqs_conserved: b1,
            min_seqs_flank: b2.max(b1),
            max_contiguous_nonconserved: 8,
            min_block_length: 10,
            gaps: GapPolicy::None,
            use_similarity: true,
        }
    }

    /// The relaxed settings of Talavera & Castresana (2007): b2 = b1, b3 = 8,
    /// b4 = 5, b5 = [`GapPolicy::Half`].
    pub fn relaxed(n_seqs: usize) -> Self {
        let b1 = n_seqs / 2 + 1;
        GblocksParams {
            min_seqs_conserved: b1,
            min_seqs_flank: b1,
            max_contiguous_nonconserved: 8,
            min_block_length: 5,
            gaps: GapPolicy::Half,
            use_similarity: true,
        }
    }

    /// Reject impossible combinations for an alignment of `n_seqs` sequences.
    ///
    /// Each message names the offending parameter and the range that would be
    /// valid, so the GUI can show it next to the offending slider.
    pub fn validate(&self, n_seqs: usize) -> Result<()> {
        if n_seqs == 0 {
            return Err(Error::algorithm(
                "Gblocks needs at least one sequence; this alignment has none",
            ));
        }
        let half = n_seqs / 2;
        if self.min_seqs_conserved <= half {
            return Err(Error::algorithm(format!(
                "b1 (minimum sequences for a conserved position) is {}, but it must be more than \
                 half of the {n_seqs} sequences: use {} to {n_seqs}",
                self.min_seqs_conserved,
                half + 1
            )));
        }
        if self.min_seqs_conserved > n_seqs {
            return Err(Error::algorithm(format!(
                "b1 (minimum sequences for a conserved position) is {}, which is more than the \
                 {n_seqs} sequences available: use {} to {n_seqs}",
                self.min_seqs_conserved,
                half + 1
            )));
        }
        if self.min_seqs_flank < self.min_seqs_conserved {
            return Err(Error::algorithm(format!(
                "b2 (minimum sequences for a flank position) is {}, but it must be at least b1 = \
                 {}: use {} to {n_seqs}",
                self.min_seqs_flank, self.min_seqs_conserved, self.min_seqs_conserved
            )));
        }
        if self.min_seqs_flank > n_seqs {
            return Err(Error::algorithm(format!(
                "b2 (minimum sequences for a flank position) is {}, which is more than the \
                 {n_seqs} sequences available: use {} to {n_seqs}",
                self.min_seqs_flank, self.min_seqs_conserved
            )));
        }
        if self.max_contiguous_nonconserved < 1 {
            return Err(Error::algorithm(
                "b3 (maximum contiguous non-conserved positions) is 0, but it must be at least 1",
            ));
        }
        if self.min_block_length < 2 {
            return Err(Error::algorithm(format!(
                "b4 (minimum block length) is {}, but it must be at least 2",
                self.min_block_length
            )));
        }
        Ok(())
    }
}

impl Default for GblocksParams {
    /// Defaults for a hypothetical 100-sequence alignment. Prefer
    /// [`GblocksParams::defaults`], which knows the real row count.
    fn default() -> Self {
        GblocksParams::defaults(100)
    }
}

/// Per-column classification, for the GUI's cleaning track.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ColumnFlag {
    /// At least b1 sequences share the most common residue.
    Conserved,
    /// At least b2 sequences share (or are similar to) the most common
    /// residue, so the column can anchor a block flank.
    HighlyConserved,
    /// Too variable to be part of a block on its own.
    NonConserved,
    /// Disqualified by the gap policy b5. Counts as non-conserved everywhere
    /// the algorithm asks "is this position conserved?".
    GapRich,
}

impl ColumnFlag {
    /// True for [`ColumnFlag::Conserved`] and [`ColumnFlag::HighlyConserved`].
    pub fn is_conserved(self) -> bool {
        matches!(self, ColumnFlag::Conserved | ColumnFlag::HighlyConserved)
    }

    /// True when this column may sit at the edge of a block.
    pub fn is_flank(self) -> bool {
        matches!(self, ColumnFlag::HighlyConserved)
    }
}

/// The outcome of running [`gblocks`] on an alignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GblocksResult {
    /// One entry per alignment column: true = keep.
    pub mask: Vec<bool>,
    /// Contiguous kept ranges, in ascending order and never touching.
    pub blocks: Vec<Range<usize>>,
    /// Conservation class of every column, whether kept or not.
    pub flags: Vec<ColumnFlag>,
    /// Number of kept columns (the width of [`GblocksResult::apply`]'s output).
    pub kept: usize,
    /// Number of columns examined.
    pub total: usize,
}

impl GblocksResult {
    /// Fraction of columns kept, 0..=1. An empty alignment scores 0.
    pub fn kept_fraction(&self) -> f32 {
        if self.total == 0 {
            0.0
        } else {
            self.kept as f32 / self.total as f32
        }
    }

    /// A new alignment with only the kept columns. Row order, names,
    /// descriptions and per-row quality are preserved.
    pub fn apply(&self, aln: &Alignment) -> Result<Alignment> {
        if self.mask.len() != aln.width() {
            return Err(Error::out_of_range(format!(
                "this Gblocks result covers {} columns but the alignment is {} wide",
                self.mask.len(),
                aln.width()
            )));
        }
        let mut out = aln.clone();
        out.pad_to_width();
        out.keep_columns(&self.mask)?;
        Ok(out)
    }

    /// The Gblocks-style mask line (`"  ####  ..."`) for display and export:
    /// one character per column, `#` where the column was selected and a space
    /// where it was rejected.
    pub fn mask_line(&self) -> String {
        self.mask.iter().map(|&keep| if keep { '#' } else { ' ' }).collect()
    }
}

/// Per-column tallies, filled once per column by the single O(rows * columns)
/// pass and then consulted by the classification.
#[derive(Clone, Copy, Default)]
struct ColumnCounts {
    /// Sequences sharing the most common residue.
    identities: usize,
    /// Sequences whose residue is identical *or similar* to the most common
    /// one. Equal to `identities` unless the similarity rule applies.
    similarities: usize,
    /// Sequences with a gap (or `?`) here.
    gaps: usize,
}

/// True for characters that carry no residue: the canonical gap and its
/// aliases, plus the missing-data marker `?`.
#[inline]
fn is_blank(c: u8) -> bool {
    is_gap(c) || c == MISSING
}

/// Select the conserved blocks of `aln`.
///
/// Returns [`Error::NotAligned`] for ragged input and [`Error::Algorithm`] for
/// parameters that cannot apply to this many sequences (see
/// [`GblocksParams::validate`]).
///
/// # Step order
///
/// The published paper (Castresana 2000) and the shipped program differ
/// slightly, and this implementation follows **the program**, because that is
/// what users compare their results against:
///
/// 1. classify every column (this is where the gap policy b5 acts: a
///    disqualified column is non-conserved by definition — the paper's own
///    wording is "nonconserved: fewer than IS identical residues *or there is
///    a gap*");
/// 2. reject stretches of more than b3 contiguous non-conserved columns;
/// 3. trim the flanks of every surviving block inwards until both ends are
///    highly conserved;
/// 4. gap cleaning: drop the gap columns that are still inside blocks, along
///    with the non-conserved columns adjacent to them, up to the next
///    conserved column;
/// 5. **finally** drop blocks shorter than b4.
///
/// The paper describes *two* length filters — `BL1` (default 15) applied right
/// after the flank step and `BL2` (default 10) applied after gap cleaning —
/// but the program exposes only one, `b4`, and its documentation is explicit
/// that it acts last: "Minimum Length Of A Block: blocks smaller than this
/// value **after gap cleaning** are rejected". So b4 is applied once, at the
/// end, to the final blocks. This also makes the parameter monotone (raising
/// b4 can only remove columns), which matters for the GUI slider.
///
/// Step 3 *trims* rather than extends. The documentation says "flanks are
/// examined and positions are **removed** until blocks are surrounded by
/// highly conserved positions at both flanks"; growing a block outwards would
/// add exactly the ambiguous columns the flank rule exists to exclude, and
/// would leave the block's own edges merely conserved rather than highly
/// conserved, contradicting the stated invariant.
pub fn gblocks(aln: &Alignment, params: &GblocksParams) -> Result<GblocksResult> {
    aln.require_aligned()?;
    let rows = aln.len();
    let width = aln.width();
    params.validate(rows)?;
    if width == 0 {
        return Ok(GblocksResult {
            mask: Vec::new(),
            blocks: Vec::new(),
            flags: Vec::new(),
            kept: 0,
            total: 0,
        });
    }

    // Gblocks' similarity rule only means anything for amino acids.
    let alphabet = aln.alphabet_hint().unwrap_or_else(|| {
        Alphabet::guess(
            aln.sequences.iter().take(100).flat_map(|s| s.residues.iter().take(5000).copied()),
        )
    });
    let similarity = params.use_similarity && !alphabet.is_nucleotide();

    let counts = tally_columns(aln, width, similarity);
    let flags = classify(&counts, rows, params);

    // Step 2: reject long non-conserved stretches.
    let mut mask = vec![true; width];
    for run in nonconserved_runs(&flags) {
        if run.len() > params.max_contiguous_nonconserved {
            mask[run].fill(false);
        }
    }

    // Step 3: trim each surviving block until both flanks are highly conserved.
    for block in selected_runs(&mask) {
        let mut start = block.start;
        let mut end = block.end;
        while start < end && !flags[start].is_flank() {
            start += 1;
        }
        while end > start && !flags[end - 1].is_flank() {
            end -= 1;
        }
        mask[block.start..start].fill(false);
        mask[end..block.end].fill(false);
    }

    // Step 4: gap cleaning. Every gap column still selected goes, and with it
    // the non-conserved columns beside it, "until a conserved position is
    // reached". A maximal non-conserved run either contains a gap column, in
    // which case the whole run is condemned, or it does not and survives, so
    // one pass over the runs does the whole job.
    if params.gaps != GapPolicy::All {
        for run in nonconserved_runs(&flags) {
            if flags[run.clone()].contains(&ColumnFlag::GapRich) {
                mask[run].fill(false);
            }
        }
    }

    // Step 5: b4, applied last, to the final blocks. Adjacent survivors have
    // already merged because `selected_runs` walks maximal runs of the mask.
    let mut blocks = Vec::new();
    for block in selected_runs(&mask) {
        if block.len() < params.min_block_length {
            mask[block].fill(false);
        } else {
            blocks.push(block);
        }
    }

    let kept = blocks.iter().map(|b| b.len()).sum();
    Ok(GblocksResult { mask, blocks, flags, kept, total: width })
}

/// One O(rows * columns) pass: for every column, how many sequences share the
/// most common residue, how many are similar to it, and how many are gaps.
///
/// The 256-entry count array is allocated once and reused, so the column loop
/// allocates nothing.
fn tally_columns(aln: &Alignment, width: usize, similarity: bool) -> Vec<ColumnCounts> {
    let mut out = vec![ColumnCounts::default(); width];
    let mut counts = [0u32; 256];
    for (col, slot) in out.iter_mut().enumerate() {
        counts.fill(0);
        let mut gaps = 0usize;
        for seq in &aln.sequences {
            let c = seq.residues.get(col).copied().unwrap_or(GAP);
            if is_blank(c) {
                gaps += 1;
            } else {
                counts[c.to_ascii_uppercase() as usize] += 1;
            }
        }
        // Most common residue. Scanning the fixed 256-entry table keeps this
        // independent of the alphabet, so unusual characters (IUPAC codes,
        // selenocysteine, ...) are counted like any other residue.
        let mut best = 0u32;
        let mut best_byte = 0u8;
        for (byte, &n) in counts.iter().enumerate() {
            if n > best {
                best = n;
                best_byte = byte as u8;
            }
        }
        slot.identities = best as usize;
        slot.gaps = gaps;
        slot.similarities = best as usize;
        if similarity && best > 0 {
            if let Some(index) = amino_acid_index(best_byte) {
                let mask = similar_mask(index);
                let mut total = 0u32;
                for (j, &aa) in AMINO_ACIDS.iter().enumerate() {
                    if mask & (1 << j) != 0 {
                        total += counts[aa as usize];
                    }
                }
                slot.similarities = total as usize;
            }
        }
    }
    out
}

/// Turn the per-column tallies into conservation classes.
///
/// The gap policy is applied *first*: a column the policy disqualifies is
/// `GapRich` and is treated as non-conserved from then on, so it can neither
/// anchor a flank nor break a long non-conserved stretch.
fn classify(counts: &[ColumnCounts], rows: usize, params: &GblocksParams) -> Vec<ColumnFlag> {
    counts
        .iter()
        .map(|c| {
            if !params.gaps.allows(c.gaps, rows) {
                return ColumnFlag::GapRich;
            }
            // The b1 test always counts identities. The b2 (flank) test may
            // also count similar residues, but only once the column already
            // has a majority of identities: "a position needs to have a number
            // of identities bigger than half the number of sequences to start
            // adding more values from similar amino acids" (Gblocks
            // documentation). Without that guard a column of twenty different
            // hydrophobics would qualify as a flank.
            let flank_count = if c.identities * 2 > rows { c.similarities } else { c.identities };
            if flank_count >= params.min_seqs_flank {
                ColumnFlag::HighlyConserved
            } else if c.identities >= params.min_seqs_conserved {
                ColumnFlag::Conserved
            } else {
                ColumnFlag::NonConserved
            }
        })
        .collect()
}

/// Maximal runs of columns that are not conserved (`NonConserved` or
/// `GapRich`).
fn nonconserved_runs(flags: &[ColumnFlag]) -> Vec<Range<usize>> {
    runs(flags.len(), |i| !flags[i].is_conserved())
}

/// Maximal runs of selected columns.
fn selected_runs(mask: &[bool]) -> Vec<Range<usize>> {
    runs(mask.len(), |i| mask[i])
}

/// Maximal runs of indices in `0..len` for which `pred` holds.
fn runs(len: usize, pred: impl Fn(usize) -> bool) -> Vec<Range<usize>> {
    let mut out: Vec<Range<usize>> = Vec::new();
    let mut start = None;
    for i in 0..len {
        match (pred(i), start) {
            (true, None) => start = Some(i),
            (false, Some(s)) => {
                out.push(s..i);
                start = None;
            }
            _ => {}
        }
    }
    if let Some(s) = start {
        out.push(s..len);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use tolviewer_core::Sequence;

    /// Build an alignment from rows given as strings, named `s0`, `s1`, ...
    fn aln(rows: &[&str]) -> Alignment {
        Alignment::new(
            "test",
            rows.iter()
                .enumerate()
                .map(|(i, r)| Sequence::new(format!("s{i}"), r.as_bytes().to_vec()))
                .collect(),
        )
    }

    /// A DNA alignment: similarity is meaningless there, so the tests below
    /// exercise pure identity counting unless they say otherwise.
    fn dna(rows: &[&str]) -> Alignment {
        let mut a = aln(rows);
        a.set_alphabet(Alphabet::Dna);
        a
    }

    fn run(a: &Alignment, p: &GblocksParams) -> GblocksResult {
        gblocks(a, p).expect("gblocks should succeed")
    }

    // ---- classification -------------------------------------------------

    #[test]
    fn fully_conserved_block_is_kept_whole() {
        let a = dna(&["ACGTACGTACGT", "ACGTACGTACGT", "ACGTACGTACGT", "ACGTACGTACGT"]);
        let r = run(&a, &GblocksParams::defaults(4));
        assert_eq!(r.total, 12);
        assert_eq!(r.kept, 12);
        assert_eq!(r.blocks, vec![0..12]);
        assert_eq!(r.mask_line(), "############");
        assert!(r.flags.iter().all(|f| *f == ColumnFlag::HighlyConserved));
        assert_eq!(r.kept_fraction(), 1.0);
    }

    #[test]
    fn classification_uses_b1_and_b2() {
        // 4 rows. Column 0 is A,A,A,A; column 1 is A,A,A,C; column 2 is
        // A,C,G,T. With b1 = 3 and b2 = 4 that is one of each class.
        let a = dna(&["AAA", "AAC", "AAG", "ACT"]);
        let mut p = GblocksParams::defaults(4);
        p.min_seqs_conserved = 3;
        p.min_seqs_flank = 4;
        let r = run(&a, &p);
        assert_eq!(
            r.flags,
            vec![ColumnFlag::HighlyConserved, ColumnFlag::Conserved, ColumnFlag::NonConserved]
        );
    }

    #[test]
    fn eight_nonconserved_columns_are_tolerated_but_nine_are_not() {
        // 4 identical flanks of 10 conserved columns, with a variable middle.
        let flank = "ACGTACGTAC";
        let make = |gap_len: usize| {
            let mut rows: Vec<String> = Vec::new();
            for i in 0..4 {
                // Each row gets a different residue in the middle, so no
                // residue reaches b1 = 3 there.
                let middle: String = std::iter::repeat_n(b"ACGT"[i] as char, gap_len).collect();
                rows.push(format!("{flank}{middle}{flank}"));
            }
            rows
        };
        let eight = make(8);
        let refs: Vec<&str> = eight.iter().map(|s| s.as_str()).collect();
        let r = run(&dna(&refs), &GblocksParams::defaults(4));
        assert_eq!(r.blocks, vec![0..28], "a run of 8 non-conserved columns must stay inside");
        assert_eq!(r.kept, 28);

        let nine = make(9);
        let refs: Vec<&str> = nine.iter().map(|s| s.as_str()).collect();
        let r = run(&dna(&refs), &GblocksParams::defaults(4));
        assert_eq!(r.blocks, vec![0..10, 19..29], "a run of 9 must split the block");
        assert_eq!(r.kept, 20);
        assert_eq!(&r.mask_line()[10..19], "         ");
    }

    #[test]
    fn block_of_nine_is_dropped_and_ten_is_kept_under_default_b4() {
        // A conserved island surrounded by non-conserved stretches longer than
        // b3, so the island is the only candidate block.
        let make = |island: usize| {
            let mut rows = Vec::new();
            for i in 0..4 {
                let n: String = std::iter::repeat_n(b"ACGT"[i] as char, 10).collect();
                rows.push(format!("{n}{}{n}", "ACGTACGTAC".get(..island).unwrap()));
            }
            rows
        };
        for (island, expect_kept) in [(9usize, 0usize), (10, 10)] {
            let rows = make(island);
            let refs: Vec<&str> = rows.iter().map(|s| s.as_str()).collect();
            let r = run(&dna(&refs), &GblocksParams::defaults(4));
            assert_eq!(r.kept, expect_kept, "island of {island} columns");
        }
    }

    #[test]
    fn flank_extension_trims_to_a_highly_conserved_column() {
        // 4 rows, b1 = 3, b2 = 4. The first two and last two columns are only
        // 3/4 conserved, so they are `Conserved` but cannot anchor a flank.
        let a = dna(&["AAAAAAAAAAAAAA", "AAAAAAAAAAAAAA", "AAAAAAAAAAAAAA", "CCAAAAAAAAAACC"]);
        let p = GblocksParams::defaults(4);
        assert_eq!((p.min_seqs_conserved, p.min_seqs_flank), (3, 4));
        let r = run(&a, &p);
        assert_eq!(r.flags[0], ColumnFlag::Conserved);
        assert_eq!(r.flags[2], ColumnFlag::HighlyConserved);
        assert_eq!(r.blocks, vec![2..12], "flanks must be trimmed in to a b2 column");
        assert_eq!(r.mask_line(), "  ##########  ");
    }

    #[test]
    fn a_block_with_no_highly_conserved_column_disappears() {
        let a = dna(&["AAAAAAAAAAAA", "AAAAAAAAAAAA", "AAAAAAAAAAAA", "CCCCCCCCCCCC"]);
        let mut p = GblocksParams::defaults(4);
        p.min_seqs_flank = 4;
        let r = run(&a, &p);
        assert_eq!(r.kept, 0);
        assert!(r.blocks.is_empty());
        assert_eq!(r.mask_line(), "            ");
    }

    // ---- gap policies ---------------------------------------------------

    /// One isolated gap in the middle of an otherwise perfect block, in 1 of 4
    /// rows (25% gaps, so `Half` still allows it).
    fn gappy() -> Alignment {
        dna(&[
            "AAAAAAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAAAAAAAA",
            "AAAAAAAAAA-AAAAAAAAA",
        ])
    }

    #[test]
    fn gap_policy_none_splits_the_block_and_b4_then_kills_the_short_side() {
        let mut p = GblocksParams::defaults(4);
        p.gaps = GapPolicy::None;
        let r = run(&gappy(), &p);
        assert_eq!(r.flags[10], ColumnFlag::GapRich);
        // Column 10 goes; the two sides are 10 and 9 long, so b4 = 10 keeps
        // only the left one.
        assert_eq!(r.blocks, vec![0..10]);
        assert_eq!(r.kept, 10);
        assert_eq!(r.mask_line(), "##########          ");
    }

    #[test]
    fn gap_policy_half_keeps_a_column_gapped_in_a_quarter_of_rows() {
        let mut p = GblocksParams::defaults(4);
        p.gaps = GapPolicy::Half;
        let r = run(&gappy(), &p);
        assert_eq!(r.flags[10], ColumnFlag::Conserved, "3 of 4 rows have A, so b1 but not b2");
        assert_eq!(r.blocks, vec![0..20]);
        assert_eq!(r.kept, 20);
    }

    #[test]
    fn gap_policy_all_keeps_everything() {
        let mut p = GblocksParams::defaults(4);
        p.gaps = GapPolicy::All;
        let r = run(&gappy(), &p);
        assert_eq!(r.blocks, vec![0..20]);
        assert_eq!(r.kept, 20);
        assert!(!r.flags.contains(&ColumnFlag::GapRich));
    }

    #[test]
    fn half_policy_rejects_a_column_gapped_in_exactly_half_the_rows() {
        let a = dna(&["AAAA", "AAAA", "AA-A", "AA-A"]);
        let mut p = GblocksParams::defaults(4);
        p.gaps = GapPolicy::Half;
        p.min_block_length = 2;
        let r = run(&a, &p);
        assert_eq!(r.flags[2], ColumnFlag::GapRich, "50% gaps is already too many");
    }

    #[test]
    fn gap_cleaning_also_removes_adjacent_nonconserved_columns() {
        // Column 10 is a gap column, column 11 is non-conserved. Both must go,
        // and the cleaning must stop at column 12, which is conserved.
        let a = dna(&[
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            "AAAAAAAAAAACAAAAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAGAAAAAAAAAAAAAAAAAA",
            "AAAAAAAAAA-TAAAAAAAAAAAAAAAAAA",
        ]);
        let mut p = GblocksParams::defaults(4);
        p.gaps = GapPolicy::None;
        p.min_block_length = 2;
        let r = run(&a, &p);
        assert_eq!(r.flags[10], ColumnFlag::GapRich);
        assert_eq!(r.flags[11], ColumnFlag::NonConserved);
        assert_eq!(r.blocks, vec![0..10, 12..30]);
    }

    // ---- similarity -----------------------------------------------------

    #[test]
    fn similarity_promotes_a_protein_flank() {
        // 4 rows; column 5 is I,I,I,V: 3 identities (a majority, so the
        // similarity top-up applies) and V is similar to I, giving 4.
        let mut a =
            aln(&["KKKKKIKKKKKKKKKK", "KKKKKIKKKKKKKKKK", "KKKKKIKKKKKKKKKK", "KKKKKVKKKKKKKKKK"]);
        a.set_alphabet(Alphabet::Protein);
        let mut p = GblocksParams::defaults(4);
        assert_eq!((p.min_seqs_conserved, p.min_seqs_flank), (3, 4));

        p.use_similarity = false;
        let r = gblocks(&a, &p).unwrap();
        assert_eq!(r.flags[5], ColumnFlag::Conserved);

        p.use_similarity = true;
        let r = gblocks(&a, &p).unwrap();
        assert_eq!(r.flags[5], ColumnFlag::HighlyConserved);
    }

    #[test]
    fn similarity_needs_a_majority_of_identities_first() {
        // I,L,M,V are pairwise similar but no residue is a majority, so the
        // column stays non-conserved however similar its residues are.
        let mut a = aln(&["ILMV", "LMVI", "MVIL", "VILM"]);
        a.set_alphabet(Alphabet::Protein);
        let mut p = GblocksParams::defaults(4);
        p.use_similarity = true;
        p.min_block_length = 2;
        let r = gblocks(&a, &p).unwrap();
        assert!(r.flags.iter().all(|f| *f == ColumnFlag::NonConserved));
        assert_eq!(r.kept, 0);
    }

    #[test]
    fn similarity_is_ignored_for_nucleotides() {
        let a = dna(&["ACGT", "ACGT", "ACGT", "ACGT"]);
        let mut p = GblocksParams::defaults(4);
        p.use_similarity = true;
        p.min_block_length = 2;
        // Identity alone already classifies these; the point is that no
        // amino-acid table is consulted, so A/G are not "similar".
        let r = run(&a, &p);
        assert!(r.flags.iter().all(|f| *f == ColumnFlag::HighlyConserved));
    }

    // ---- parameters -----------------------------------------------------

    #[test]
    fn defaults_match_gblocks() {
        let p = GblocksParams::defaults(20);
        assert_eq!(p.min_seqs_conserved, 11); // 20/2 + 1
        assert_eq!(p.min_seqs_flank, 17); // ceil(20 * 0.85)
        assert_eq!(p.max_contiguous_nonconserved, 8);
        assert_eq!(p.min_block_length, 10);
        assert_eq!(p.gaps, GapPolicy::None);
        assert!(p.use_similarity);
        assert_eq!(GblocksParams::defaults(7).min_seqs_flank, 6); // ceil(5.95)
        assert_eq!(GblocksParams::defaults(3).min_seqs_conserved, 2);
        // b2 is never allowed below b1.
        let small = GblocksParams::defaults(2);
        assert!(small.min_seqs_flank >= small.min_seqs_conserved);
    }

    #[test]
    fn relaxed_matches_talavera_and_castresana() {
        let p = GblocksParams::relaxed(20);
        assert_eq!(p.min_seqs_conserved, 11);
        assert_eq!(p.min_seqs_flank, 11);
        assert_eq!(p.max_contiguous_nonconserved, 8);
        assert_eq!(p.min_block_length, 5);
        assert_eq!(p.gaps, GapPolicy::Half);
    }

    #[test]
    fn defaults_and_relaxed_validate() {
        for n in 1..64 {
            GblocksParams::defaults(n).validate(n).unwrap_or_else(|e| panic!("n = {n}: {e}"));
            GblocksParams::relaxed(n).validate(n).unwrap_or_else(|e| panic!("n = {n}: {e}"));
        }
    }

    fn rejection_message(p: &GblocksParams, n: usize) -> String {
        match p.validate(n) {
            Err(Error::Algorithm(m)) => m,
            other => panic!("expected Error::Algorithm, got {other:?}"),
        }
    }

    #[test]
    fn validate_names_the_offending_parameter() {
        let base = GblocksParams::defaults(10);

        // n == 0
        let m = rejection_message(&base, 0);
        assert!(m.contains("at least one sequence"), "{m}");

        // b1 <= n/2
        let mut p = base.clone();
        p.min_seqs_conserved = 5;
        p.min_seqs_flank = 9;
        let m = rejection_message(&p, 10);
        assert!(m.contains("b1") && m.contains("half") && m.contains("6 to 10"), "{m}");

        // b1 > n
        let mut p = base.clone();
        p.min_seqs_conserved = 11;
        p.min_seqs_flank = 11;
        let m = rejection_message(&p, 10);
        assert!(m.contains("b1") && m.contains("more than the 10 sequences"), "{m}");

        // b2 < b1
        let mut p = base.clone();
        p.min_seqs_flank = 5;
        let m = rejection_message(&p, 10);
        assert!(m.contains("b2") && m.contains("at least b1 = 6"), "{m}");

        // b2 > n
        let mut p = base.clone();
        p.min_seqs_flank = 11;
        let m = rejection_message(&p, 10);
        assert!(m.contains("b2") && m.contains("more than the 10 sequences"), "{m}");

        // b3 < 1
        let mut p = base.clone();
        p.max_contiguous_nonconserved = 0;
        let m = rejection_message(&p, 10);
        assert!(m.contains("b3") && m.contains("at least 1"), "{m}");

        // b4 < 2
        let mut p = base.clone();
        p.min_block_length = 1;
        let m = rejection_message(&p, 10);
        assert!(m.contains("b4") && m.contains("at least 2"), "{m}");
    }

    #[test]
    fn gblocks_propagates_validation_errors() {
        let a = dna(&["ACGT", "ACGT"]);
        let mut p = GblocksParams::defaults(2);
        p.min_block_length = 1;
        assert!(matches!(gblocks(&a, &p), Err(Error::Algorithm(_))));
    }

    #[test]
    fn ragged_input_is_rejected() {
        let a = aln(&["ACGT", "AC"]);
        assert!(matches!(gblocks(&a, &GblocksParams::defaults(2)), Err(Error::NotAligned)));
    }

    #[test]
    fn empty_alignment_is_rejected_but_zero_width_is_not() {
        let empty = Alignment::new("e", Vec::new());
        assert!(matches!(gblocks(&empty, &GblocksParams::defaults(1)), Err(Error::Algorithm(_))));

        let zero_width = dna(&["", ""]);
        let r = run(&zero_width, &GblocksParams::defaults(2));
        assert_eq!((r.total, r.kept), (0, 0));
        assert_eq!(r.kept_fraction(), 0.0);
        assert_eq!(r.mask_line(), "");
    }

    // ---- result invariants ----------------------------------------------

    #[test]
    fn apply_agrees_with_mask_blocks_and_mask_line() {
        let a = gappy();
        let r = run(&a, &GblocksParams::defaults(4));
        let out = r.apply(&a).unwrap();

        assert_eq!(out.width(), r.kept);
        assert_eq!(out.len(), a.len());
        let before: Vec<&str> = a.sequences.iter().map(|s| s.id.as_str()).collect();
        let after: Vec<&str> = out.sequences.iter().map(|s| s.id.as_str()).collect();
        assert_eq!(before, after);

        // mask <-> blocks
        let from_blocks: Vec<bool> = {
            let mut m = vec![false; r.total];
            for b in &r.blocks {
                m[b.clone()].fill(true);
            }
            m
        };
        assert_eq!(from_blocks, r.mask);
        assert_eq!(r.kept, r.mask.iter().filter(|&&k| k).count());

        // mask <-> mask_line
        let line = r.mask_line();
        assert_eq!(line.chars().count(), r.total);
        for (c, &keep) in line.chars().zip(&r.mask) {
            assert_eq!(c, if keep { '#' } else { ' ' });
        }

        // blocks are ascending, non-empty and never adjacent.
        for w in r.blocks.windows(2) {
            assert!(w[0].end < w[1].start, "blocks {:?} and {:?} should have merged", w[0], w[1]);
        }
        assert!(r.blocks.iter().all(|b| !b.is_empty()));
    }

    #[test]
    fn apply_rejects_a_mismatched_alignment() {
        let a = gappy();
        let r = run(&a, &GblocksParams::defaults(4));
        let other = dna(&["ACGT", "ACGT", "ACGT", "ACGT"]);
        assert!(matches!(r.apply(&other), Err(Error::OutOfRange(_))));
    }

    #[test]
    fn deterministic_and_idempotent() {
        let a = gappy();
        let p = GblocksParams::defaults(4);
        let first = run(&a, &p);
        let second = run(&a, &p);
        assert_eq!(first, second, "gblocks must be deterministic");

        // Re-running on the cleaned alignment keeps everything.
        let cleaned = first.apply(&a).unwrap();
        let again = run(&cleaned, &p);
        assert_eq!(again.kept, again.total);
        assert_eq!(again.blocks, vec![0..cleaned.width()]);
        assert_eq!(again.apply(&cleaned).unwrap(), cleaned);
    }

    #[test]
    fn gap_policy_helpers() {
        assert_eq!(GapPolicy::all().len(), 3);
        assert_eq!(GapPolicy::None.letter(), 'n');
        assert_eq!(GapPolicy::Half.name(), "half");
        assert!(ColumnFlag::HighlyConserved.is_conserved());
        assert!(ColumnFlag::HighlyConserved.is_flank());
        assert!(ColumnFlag::Conserved.is_conserved() && !ColumnFlag::Conserved.is_flank());
        assert!(!ColumnFlag::GapRich.is_conserved());
        assert_eq!(GblocksParams::default(), GblocksParams::defaults(100));
    }

    #[test]
    fn missing_data_marker_counts_as_a_gap() {
        let a = dna(&["AAAA", "AAAA", "AAAA", "AA?A"]);
        let mut p = GblocksParams::defaults(4);
        p.min_block_length = 2;
        let r = run(&a, &p);
        assert_eq!(r.flags[2], ColumnFlag::GapRich);
    }

    #[test]
    fn residues_are_compared_case_insensitively() {
        let a = dna(&["acgtacgtacgt", "ACGTACGTACGT", "acGTacGTacGT", "ACgtACgtACgt"]);
        let r = run(&a, &GblocksParams::defaults(4));
        assert_eq!(r.kept, 12);
        // ... and the kept residues keep their original case.
        let out = r.apply(&a).unwrap();
        assert_eq!(out.sequences[0].residues, b"acgtacgtacgt");
    }

    // ---- performance ----------------------------------------------------

    /// 500 sequences x 5000 columns. Run with
    /// `cargo test -p tolviewer-clean --release -- --ignored --nocapture`.
    #[test]
    #[ignore = "benchmark; run explicitly in release mode"]
    fn benchmark_500x5000() {
        let rows = 500;
        let cols = 5000;
        let letters = b"ACDEFGHIKLMNPQRSTVWY";
        let mut state = 0x2545_F491_4F6C_DD1Du64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let consensus: Vec<u8> = (0..cols).map(|_| letters[(next() % 20) as usize]).collect();
        let mut sequences = Vec::with_capacity(rows);
        for i in 0..rows {
            let mut residues = Vec::with_capacity(cols);
            for &c in &consensus {
                residues.push(match next() % 100 {
                    0..=79 => c,
                    80..=94 => letters[(next() % 20) as usize],
                    _ => GAP,
                });
            }
            sequences.push(Sequence::new(format!("s{i}"), residues));
        }
        let mut a = Alignment::new("bench", sequences);
        a.set_alphabet(Alphabet::Protein);

        // Both parameter sets, because they exercise different branches: the
        // defaults reject almost everything (b5 = None on 5% gappy data), the
        // relaxed settings keep almost everything.
        for (label, p) in
            [("default", GblocksParams::defaults(rows)), ("relaxed", GblocksParams::relaxed(rows))]
        {
            // Warm up, then time ten runs (the GUI re-runs this per keystroke).
            let mut kept = gblocks(&a, &p).unwrap().kept;
            let start = std::time::Instant::now();
            for _ in 0..10 {
                kept = gblocks(&a, &p).unwrap().kept;
            }
            let each = start.elapsed() / 10;
            println!("gblocks {label} on {rows}x{cols}: {each:?} per run, kept {kept} columns");
            assert!(each < std::time::Duration::from_secs(2), "too slow: {each:?}");
        }
    }
}
