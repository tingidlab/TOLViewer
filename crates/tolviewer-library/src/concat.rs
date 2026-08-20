//! Concatenating per-locus alignments into one supermatrix.
//!
//! A multi-locus phylogeny is built from one alignment per gene, joined
//! end-to-end so that each row is one specimen across every gene. The awkward
//! part is never the joining; it is deciding that `TL-2213_18S_F` and
//! `TL_2213_28S` are the same animal. [`crate::naming`] does that, and this
//! module turns the resulting groups into a matrix.
//!
//! Specimens missing from a locus are filled with gaps for its whole width,
//! which is what every phylogenetic program expects and treats as missing data.
//! Which specimens were filled in is reported rather than hidden, because a row
//! that is 80% gaps is usually a name that failed to match rather than a real
//! gap in the sampling.

use std::ops::Range;

use tolviewer_core::{Alignment, Error, Result, Sequence, GAP};

use crate::naming::{sample_key, MatchOptions};

/// How to join alignments.
#[derive(Debug, Clone, PartialEq)]
pub struct ConcatOptions {
    /// How names are reduced to samples before matching.
    pub matching: MatchOptions,
    /// Keep samples that are absent from some of the alignments, padding them
    /// with gaps. With this off, only samples present in every alignment are
    /// kept, which is the strict "complete matrix" behaviour.
    pub include_partial: bool,
}

impl Default for ConcatOptions {
    fn default() -> Self {
        ConcatOptions { matching: MatchOptions::default(), include_partial: true }
    }
}

/// One input alignment's span in the concatenated matrix.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Partition {
    pub name: String,
    /// Columns, 0-based and half-open, as [`Alignment`] indexes them.
    pub range: Range<usize>,
}

impl Partition {
    /// The 1-based inclusive span NEXUS and RAxML write.
    pub fn as_charset(&self) -> String {
        format!("{}-{}", self.range.start + 1, self.range.end)
    }
}

/// A sample that was not in every alignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MissingSample {
    /// The name the sample appears under in the result.
    pub sample: String,
    /// The alignments it was absent from.
    pub absent_from: Vec<String>,
}

/// The supermatrix and an account of how it was built.
#[derive(Debug, Clone)]
pub struct ConcatResult {
    pub alignment: Alignment,
    pub partitions: Vec<Partition>,
    /// Samples padded with gaps for at least one locus, in output row order.
    pub missing: Vec<MissingSample>,
    /// Samples found in every input.
    pub complete: usize,
    /// Samples dropped because `include_partial` was off.
    pub dropped: Vec<String>,
}

impl ConcatResult {
    /// The partitions as a NEXUS `sets` block, ready to paste into a file for
    /// a partitioned analysis.
    pub fn nexus_charsets(&self) -> String {
        let mut out = String::from("begin sets;\n");
        for p in &self.partitions {
            out.push_str(&format!("    charset {} = {};\n", sanitize(&p.name), p.as_charset()));
        }
        out.push_str("end;\n");
        out
    }

    /// The partitions in RAxML's format, one line per locus.
    pub fn raxml_partitions(&self, model: &str) -> String {
        self.partitions
            .iter()
            .map(|p| format!("{model}, {} = {}\n", sanitize(&p.name), p.as_charset()))
            .collect()
    }
}

/// NEXUS and RAxML both choke on whitespace and punctuation in a set name.
fn sanitize(name: &str) -> String {
    let cleaned: String =
        name.chars().map(|c| if c.is_ascii_alphanumeric() || c == '_' { c } else { '_' }).collect();
    let trimmed = cleaned.trim_matches('_');
    if trimmed.is_empty() {
        "locus".to_string()
    } else if trimmed.starts_with(|c: char| c.is_ascii_digit()) {
        // A leading digit is not a legal identifier in either format.
        format!("p_{trimmed}")
    } else {
        trimmed.to_string()
    }
}

/// Join `parts` end-to-end, matching rows by sample.
///
/// Every input must already be aligned; concatenating ragged rows would
/// silently shift every locus after the first. Row order follows the first
/// alignment, then each later one contributes its new samples in its own order,
/// so a matrix built from the same files is byte-identical every time.
pub fn concatenate(parts: &[&Alignment], opts: &ConcatOptions) -> Result<ConcatResult> {
    if parts.len() < 2 {
        return Err(Error::format("concatenating needs at least two alignments"));
    }
    for aln in parts {
        aln.require_aligned().map_err(|_| {
            Error::format(format!(
                "'{}' has rows of different lengths, so it is not an alignment yet; \
                 align it before concatenating",
                aln.name
            ))
        })?;
    }

    // Index each alignment by sample key, refusing an ambiguous one: two rows
    // of one locus claiming the same specimen is a data problem the user has to
    // resolve, and picking one silently would put the wrong gene in the matrix.
    let mut indexes: Vec<Vec<(String, usize)>> = Vec::with_capacity(parts.len());
    for aln in parts {
        let mut index: Vec<(String, usize)> = Vec::with_capacity(aln.len());
        for (row, seq) in aln.sequences.iter().enumerate() {
            let key = sample_key(&seq.id, &opts.matching);
            if let Some((_, other)) = index.iter().find(|(k, _)| *k == key) {
                return Err(Error::format(format!(
                    "in '{}', '{}' and '{}' both look like sample '{}'; \
                     rename one or turn off name matching",
                    aln.name, aln.sequences[*other].id, seq.id, key
                )));
            }
            index.push((key, row));
        }
        indexes.push(index);
    }

    // Output row order: first alignment's order, then newcomers in the order
    // each later alignment introduces them.
    let mut order: Vec<(String, String)> = Vec::new(); // (key, display name)
    for (aln, index) in parts.iter().zip(&indexes) {
        for (key, row) in index {
            if !order.iter().any(|(k, _)| k == key) {
                order.push((key.clone(), aln.sequences[*row].id.clone()));
            }
        }
    }

    let widths: Vec<usize> = parts.iter().map(|a| a.width()).collect();
    let mut partitions = Vec::with_capacity(parts.len());
    let mut at = 0;
    for (aln, &width) in parts.iter().zip(&widths) {
        partitions.push(Partition { name: aln.name.clone(), range: at..at + width });
        at += width;
    }
    let total = at;

    let mut sequences = Vec::new();
    let mut missing = Vec::new();
    let mut dropped = Vec::new();
    let mut complete = 0;

    for (key, display) in &order {
        let mut absent: Vec<String> = Vec::new();
        let mut residues = Vec::with_capacity(total);
        for ((aln, index), &width) in parts.iter().zip(&indexes).zip(&widths) {
            match index.iter().find(|(k, _)| k == key) {
                Some((_, row)) => {
                    let seq = &aln.sequences[*row];
                    residues.extend_from_slice(&seq.residues);
                    // `require_aligned` passed, so this only pads the degenerate
                    // zero-row case; keeping it makes the width exact regardless.
                    residues.resize(residues.len() + width.saturating_sub(seq.len()), GAP);
                }
                None => {
                    absent.push(aln.name.clone());
                    residues.resize(residues.len() + width, GAP);
                }
            }
        }
        if absent.is_empty() {
            complete += 1;
        } else if !opts.include_partial {
            dropped.push(display.clone());
            continue;
        } else {
            missing.push(MissingSample { sample: display.clone(), absent_from: absent });
        }
        sequences.push(Sequence::new(display.clone(), residues));
    }

    let name = parts.iter().map(|a| a.name.as_str()).collect::<Vec<_>>().join("+");
    let mut alignment = Alignment::new(name, sequences);
    if let Some(first) = parts.first().and_then(|a| a.alphabet_hint()) {
        alignment.set_alphabet(first);
    }
    Ok(ConcatResult { alignment, partitions, missing, complete, dropped })
}

/// What concatenating *would* do, without building the matrix.
///
/// The GUI shows this before committing, because a name-matching mistake is
/// invisible in the result but obvious in a table of which samples were found
/// in which locus.
pub fn preview(parts: &[&Alignment], opts: &ConcatOptions) -> Vec<SamplePreview> {
    let mut previews: Vec<SamplePreview> = Vec::new();
    for (i, aln) in parts.iter().enumerate() {
        for seq in &aln.sequences {
            let key = sample_key(&seq.id, &opts.matching);
            match previews.iter_mut().find(|p| p.key == key) {
                Some(p) => p.found_in.push((i, seq.id.clone())),
                None => previews.push(SamplePreview {
                    key,
                    display: seq.id.clone(),
                    found_in: vec![(i, seq.id.clone())],
                }),
            }
        }
    }
    previews
}

/// One sample as the matcher sees it, for the confirmation table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SamplePreview {
    /// The key the names were reduced to.
    pub key: String,
    /// The name the output row would carry.
    pub display: String,
    /// `(alignment index, the name it appears under there)`.
    pub found_in: Vec<(usize, String)>,
}

impl SamplePreview {
    /// Is this sample in every one of `total` alignments?
    pub fn is_complete(&self, total: usize) -> bool {
        let mut seen: Vec<usize> = self.found_in.iter().map(|(i, _)| *i).collect();
        seen.sort_unstable();
        seen.dedup();
        seen.len() == total
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aln(name: &str, rows: &[(&str, &str)]) -> Alignment {
        Alignment::new(
            name,
            rows.iter().map(|(id, r)| Sequence::new(*id, r.as_bytes().to_vec())).collect(),
        )
    }

    fn ssu() -> Alignment {
        aln("18S", &[("TL-2213_18S_F", "ACGTACGT"), ("TL-2214_18S_F", "ACGAACGT")])
    }

    fn lsu() -> Alignment {
        aln("28S", &[("TL_2214_28S", "TTTTGGGG"), ("TL_2213_28S", "TTTTCCCC")])
    }

    #[test]
    fn matching_names_are_joined_across_loci() {
        let r = concatenate(&[&ssu(), &lsu()], &ConcatOptions::default()).unwrap();
        assert_eq!(r.alignment.len(), 2);
        assert_eq!(r.alignment.width(), 16);
        assert!(r.alignment.is_aligned());
        assert_eq!(r.complete, 2);
        assert!(r.missing.is_empty());
        // Row order follows the first alignment, not the second.
        assert_eq!(r.alignment.sequences[0].id, "TL-2213_18S_F");
        assert_eq!(r.alignment.sequences[0].residues, b"ACGTACGTTTTTCCCC");
        assert_eq!(r.alignment.sequences[1].residues, b"ACGAACGTTTTTGGGG");
    }

    #[test]
    fn partitions_span_the_matrix_exactly() {
        let r = concatenate(&[&ssu(), &lsu()], &ConcatOptions::default()).unwrap();
        assert_eq!(r.partitions[0], Partition { name: "18S".into(), range: 0..8 });
        assert_eq!(r.partitions[1], Partition { name: "28S".into(), range: 8..16 });
        assert_eq!(r.partitions.last().unwrap().range.end, r.alignment.width());
        assert!(r.nexus_charsets().contains("charset p_18S = 1-8;"), "{}", r.nexus_charsets());
        assert!(r.nexus_charsets().contains("charset p_28S = 9-16;"));
        assert!(r.raxml_partitions("GTR").contains("GTR, p_18S = 1-8"));
    }

    #[test]
    fn a_sample_missing_from_a_locus_is_gapped_and_reported() {
        let extra = aln("COI", &[("TL-2213_COI", "GGGG"), ("TL-9999_COI", "CCCC")]);
        let r = concatenate(&[&ssu(), &lsu(), &extra], &ConcatOptions::default()).unwrap();
        assert_eq!(r.alignment.len(), 3);
        assert_eq!(r.complete, 1, "only TL-2213 is in all three");

        let row = r.alignment.find_by_id("TL-2214_18S_F").unwrap();
        assert_eq!(&r.alignment.sequences[row].residues[16..], b"----");

        let ninety_nine = r.missing.iter().find(|m| m.sample == "TL-9999_COI").unwrap();
        assert_eq!(ninety_nine.absent_from, vec!["18S".to_string(), "28S".to_string()]);
        assert!(r.dropped.is_empty());
    }

    #[test]
    fn a_strict_matrix_drops_incomplete_samples() {
        let extra = aln("COI", &[("TL-2213_COI", "GGGG")]);
        let opts = ConcatOptions { include_partial: false, ..Default::default() };
        let r = concatenate(&[&ssu(), &lsu(), &extra], &opts).unwrap();
        assert_eq!(r.alignment.len(), 1);
        assert_eq!(r.alignment.sequences[0].id, "TL-2213_18S_F");
        assert_eq!(r.dropped, vec!["TL-2214_18S_F".to_string()]);
        assert!(r.missing.is_empty());
    }

    #[test]
    fn two_rows_claiming_one_sample_is_an_error_naming_both() {
        let clashing = aln("18S", &[("TL-2213_18S_F", "ACGT"), ("TL-2213_18S_R", "ACGA")]);
        let e = concatenate(&[&clashing, &lsu()], &ConcatOptions::default()).unwrap_err();
        let text = e.to_string();
        assert!(text.contains("TL-2213_18S_F") && text.contains("TL-2213_18S_R"), "{text}");
        assert!(matches!(e, Error::Format(_)));
    }

    #[test]
    fn exact_matching_keeps_differently_written_names_apart() {
        let opts = ConcatOptions { matching: MatchOptions::exact(), ..Default::default() };
        let r = concatenate(&[&ssu(), &lsu()], &opts).unwrap();
        assert_eq!(r.alignment.len(), 4, "nothing matches, so every row is its own sample");
        assert_eq!(r.complete, 0);
    }

    #[test]
    fn ragged_input_is_refused_by_name() {
        let ragged = aln("18S", &[("a_18S", "ACGT"), ("b_18S", "ACG")]);
        let e = concatenate(&[&ragged, &lsu()], &ConcatOptions::default()).unwrap_err();
        assert!(e.to_string().contains("18S"), "{e}");
        assert!(e.to_string().contains("align it"), "{e}");
    }

    #[test]
    fn fewer_than_two_alignments_is_not_a_concatenation() {
        assert!(concatenate(&[], &ConcatOptions::default()).is_err());
        assert!(concatenate(&[&ssu()], &ConcatOptions::default()).is_err());
    }

    #[test]
    fn the_preview_shows_what_matched_before_anything_is_built() {
        let extra = aln("COI", &[("TL-2213_COI", "GGGG")]);
        let p = preview(&[&ssu(), &lsu(), &extra], &ConcatOptions::default());
        assert_eq!(p.len(), 2);
        let first = p.iter().find(|s| s.key == "tl_2213").unwrap();
        assert!(first.is_complete(3));
        assert_eq!(first.found_in.len(), 3);
        assert_eq!(first.found_in[2].1, "TL-2213_COI");
        let second = p.iter().find(|s| s.key == "tl_2214").unwrap();
        assert!(!second.is_complete(3));
    }

    #[test]
    fn set_names_are_made_legal_without_colliding_by_accident() {
        assert_eq!(sanitize("18S rRNA"), "p_18S_rRNA");
        assert_eq!(sanitize("wingless"), "wingless");
        assert_eq!(sanitize("__"), "locus");
        assert_eq!(sanitize(""), "locus");
        assert_eq!(sanitize("COI (Folmer)"), "COI__Folmer");
    }

    #[test]
    fn concatenation_is_deterministic() {
        let a = concatenate(&[&ssu(), &lsu()], &ConcatOptions::default()).unwrap();
        let b = concatenate(&[&ssu(), &lsu()], &ConcatOptions::default()).unwrap();
        assert_eq!(a.alignment.sequences, b.alignment.sequences);
    }
}
