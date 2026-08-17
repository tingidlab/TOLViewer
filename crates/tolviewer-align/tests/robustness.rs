//! Degenerate inputs, cancellation, and the invariants `realign_region` and
//! `add_to_alignment` must never break.
//!
//! None of these are performance tests; they exist because every one of them
//! is a way to make an aligner panic, hang, or silently corrupt a user's data.

mod common;

use std::sync::atomic::{AtomicUsize, Ordering};

use tolviewer_align::{
    add_to_alignment, align, realign_region, AlignParams, Engine, NoProgress, Progress,
};
use tolviewer_core::alphabet::is_gap;
use tolviewer_core::{Alignment, Error, Sequence};

fn aln(rows: &[(&str, &str)]) -> Alignment {
    Alignment::new(
        "t",
        rows.iter().map(|(id, r)| Sequence::new(*id, r.as_bytes().to_vec())).collect(),
    )
}

fn all_params() -> Vec<AlignParams> {
    Engine::all().iter().map(|&e| AlignParams::for_engine(e)).collect()
}

/// Every row must keep exactly its input residues, in order.
fn assert_residues_preserved(input: &Alignment, output: &Alignment) {
    assert_eq!(input.len(), output.len(), "row count changed");
    for (i, (a, b)) in input.sequences.iter().zip(output.sequences.iter()).enumerate() {
        assert_eq!(a.id, b.id, "row {i} was renamed");
        assert_eq!(a.ungapped(), b.ungapped(), "row {i} lost or gained residues");
    }
}

#[test]
fn zero_sequences() {
    let a = Alignment::new("empty", Vec::new());
    for p in all_params() {
        let out = align(&a, &p, &NoProgress).expect("no sequences is not an error");
        assert_eq!(out.len(), 0);
        assert!(out.is_aligned());
    }
}

#[test]
fn one_sequence() {
    let a = aln(&[("only", "ACGT--ACGT")]);
    for p in all_params() {
        let out = align(&a, &p, &NoProgress).expect("one sequence is not an error");
        assert_eq!(out.len(), 1);
        // Gaps are stripped first, so a lone sequence comes back ungapped.
        assert_eq!(out.sequences[0].residues, b"ACGTACGT".to_vec());
    }
}

#[test]
fn duplicate_sequences() {
    let a = aln(&[("a", "ACGTACGT"), ("b", "ACGTACGT"), ("c", "ACGTACGT"), ("d", "ACGTACGT")]);
    for p in all_params() {
        let out = align(&a, &p, &NoProgress).expect("duplicates are not an error");
        assert!(out.is_aligned());
        assert_eq!(out.width(), 8, "{:?} invented columns", p.engine);
        assert_residues_preserved(&a, &out);
    }
}

#[test]
fn zero_length_sequences() {
    let a = aln(&[("a", ""), ("b", ""), ("c", "")]);
    for p in all_params() {
        let out = align(&a, &p, &NoProgress).expect("empty rows are not an error");
        assert_eq!(out.len(), 3);
        assert!(out.is_aligned());
        assert_eq!(out.width(), 0);
    }
}

#[test]
fn rows_that_are_only_gaps() {
    let a = aln(&[("a", "ACGTACGT"), ("gapsonly", "--------"), ("c", "ACGTACGT")]);
    for p in all_params() {
        let out = align(&a, &p, &NoProgress).expect("all-gap rows are not an error");
        assert!(out.is_aligned());
        assert_residues_preserved(&a, &out);
        assert!(out.sequences[1].residues.iter().all(|&c| is_gap(c)));
    }
}

#[test]
fn mixture_of_empty_and_non_empty() {
    let a = aln(&[("a", "ACGTACGTAC"), ("b", ""), ("c", "ACGTACGTAC"), ("d", "--")]);
    for p in all_params() {
        let out = align(&a, &p, &NoProgress).expect("mixed lengths are not an error");
        assert!(out.is_aligned());
        assert_residues_preserved(&a, &out);
    }
}

#[test]
fn wildly_unequal_lengths() {
    let a = aln(&[
        ("short", "ACGT"),
        ("long", "TTTTTTTTTTACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"),
        ("mid", "TTTTTACGTTTTT"),
    ]);
    for p in all_params() {
        let out = align(&a, &p, &NoProgress).expect("ragged input is not an error");
        assert!(out.is_aligned());
        assert_residues_preserved(&a, &out);
    }
}

#[test]
fn lowercase_and_ambiguity_codes_survive() {
    let a = aln(&[("a", "acgtNNRYacgt"), ("b", "acgtNNRYacgt"), ("c", "acgtRYacgt")]);
    for p in all_params() {
        let out = align(&a, &p, &NoProgress).expect("masked residues are not an error");
        assert!(out.is_aligned());
        assert_residues_preserved(&a, &out);
        assert_eq!(out.sequences[0].residues.iter().filter(|c| c.is_ascii_lowercase()).count(), 8);
    }
}

/// A `Progress` that cancels after `limit` ticks.
struct CancelAfter {
    seen: AtomicUsize,
    limit: usize,
}

impl CancelAfter {
    fn new(limit: usize) -> Self {
        CancelAfter { seen: AtomicUsize::new(0), limit }
    }
}

impl Progress for CancelAfter {
    fn tick(&self, _fraction: f32, _message: &str) -> bool {
        self.seen.fetch_add(1, Ordering::SeqCst) < self.limit
    }
}

#[test]
fn cancellation_returns_cancelled_and_does_not_hang() {
    // Enough sequences that every engine ticks more than a handful of times.
    let seqs: Vec<Sequence> = (0..40)
        .map(|i| {
            let mut r = b"ACGTACGTTTGACCATTGACCAGGTTACGATCGAT".to_vec();
            let pos = i % r.len();
            r[pos] = b'A';
            Sequence::new(format!("s{i}"), r)
        })
        .collect();
    let a = Alignment::new("cancel", seqs);
    for p in all_params() {
        for limit in [0usize, 1, 3] {
            let prog = CancelAfter::new(limit);
            let r = align(&a, &p, &prog);
            assert!(
                matches!(r, Err(Error::Cancelled)),
                "{:?} with limit {limit} returned {:?}",
                p.engine,
                r.map(|x| x.width())
            );
        }
    }
}

#[test]
fn cancellation_of_realign_and_add() {
    let a = aln(&[
        ("a", "ACGTACGTACGTACGTACGT"),
        ("b", "ACGTACGTTTACGTACGTAC"),
        ("c", "ACGTACGTACGTACGTACGT"),
        ("d", "ACGTACGAACGTACGTACGT"),
    ]);
    let p = AlignParams::for_engine(Engine::Muscle);
    let r = realign_region(&a, 4..16, &p, &CancelAfter::new(0));
    assert!(matches!(r, Err(Error::Cancelled)));
    let q = vec![Sequence::new("q1", b"ACGTACGTACGTACGTACGT".to_vec())];
    let r = add_to_alignment(&a, &q, &p, &CancelAfter::new(0));
    assert!(matches!(r, Err(Error::Cancelled)));
}

#[test]
fn realign_region_keeps_the_rest_byte_identical() {
    let a = aln(&[
        ("a", "AAAA-CGTACGTAC-TTTT"),
        ("b", "AAAACCGTTACGTACTTTT"),
        ("c", "AAAA--GTACGGTACTTTT"),
        ("d", "AAAACCGTACGTACCTTTT"),
    ]);
    let width = a.width();
    for p in all_params() {
        for (start, end) in [(4usize, 15usize), (0, 6), (10, width), (7, 7), (0, width)] {
            let out = realign_region(&a, start..end, &p, &NoProgress)
                .unwrap_or_else(|e| panic!("realign {start}..{end} failed: {e}"));
            assert!(out.is_aligned(), "ragged result for {start}..{end}");
            assert_residues_preserved(&a, &out);

            // Everything left of the selection is untouched.
            for (i, s) in out.sequences.iter().enumerate() {
                assert_eq!(
                    &s.residues[..start],
                    &a.sequences[i].residues[..start],
                    "row {i}: columns before {start} changed"
                );
                // ... and so is everything right of it, allowing for the block
                // having changed width.
                let tail = width - end;
                let out_tail = &s.residues[s.residues.len() - tail..];
                assert_eq!(
                    out_tail,
                    &a.sequences[i].residues[end..],
                    "row {i}: columns after {end} changed"
                );
            }
        }
    }
}

#[test]
fn realign_region_handles_silly_ranges() {
    let a = aln(&[("a", "ACGT"), ("b", "ACGT")]);
    let p = AlignParams::default();
    for (s, e) in [(0usize, 0usize), (4, 4), (2, 2), (0, 99), (99, 99)] {
        let out = realign_region(&a, s..e, &p, &NoProgress).expect("clamped, not rejected");
        assert!(out.is_aligned());
        assert_residues_preserved(&a, &out);
    }
}

#[test]
fn realign_region_on_an_empty_alignment() {
    let a = Alignment::new("empty", Vec::new());
    let out = realign_region(&a, 0..10, &AlignParams::default(), &NoProgress).unwrap();
    assert_eq!(out.len(), 0);
}

#[test]
fn add_to_alignment_keeps_the_profile_columns() {
    let profile =
        aln(&[("p1", "ACGTACGTACGTACGT"), ("p2", "ACGT--GTACGTACGT"), ("p3", "ACGTACGTACGTACGT")]);
    let query = vec![
        Sequence::new("q1", b"ACGTACGTTTTACGTACGT".to_vec()),
        Sequence::new("q2", b"ACGTACGTACGTACGT".to_vec()),
    ];
    for p in all_params() {
        let out = add_to_alignment(&profile, &query, &p, &NoProgress).expect("adds");
        assert!(out.is_aligned());
        assert_eq!(out.len(), 5);
        assert_eq!(out.sequences[3].id, "q1");
        assert_eq!(out.sequences[4].id, "q2");
        for (i, s) in out.sequences.iter().take(3).enumerate() {
            assert_eq!(s.id, profile.sequences[i].id);
            assert_eq!(s.ungapped(), profile.sequences[i].ungapped());
        }
        for (i, s) in out.sequences.iter().skip(3).enumerate() {
            assert_eq!(s.ungapped(), query[i].residues);
        }

        // Deleting the columns that are gaps in every profile row must restore
        // the profile exactly: the query can add columns, never re-gap the
        // profile's own alignment.
        let recovered: Vec<Vec<u8>> = (0..3)
            .map(|r| {
                (0..out.width())
                    .filter(|&c| (0..3).any(|q| !is_gap(out.sequences[q].residues[c])))
                    .map(|c| out.sequences[r].residues[c])
                    .collect()
            })
            .collect();
        for (r, row) in recovered.iter().enumerate() {
            assert_eq!(
                row, &profile.sequences[r].residues,
                "{:?}: profile row {r} was re-gapped",
                p.engine
            );
        }
    }
}

#[test]
fn add_to_alignment_degenerate_inputs() {
    let profile = aln(&[("p1", "ACGTACGT"), ("p2", "ACGTACGT")]);
    let p = AlignParams::default();

    // No queries at all.
    let out = add_to_alignment(&profile, &[], &p, &NoProgress).unwrap();
    assert_eq!(out.len(), 2);

    // Empty profile.
    let empty = Alignment::new("e", Vec::new());
    let q =
        vec![Sequence::new("q1", b"ACGTACGT".to_vec()), Sequence::new("q2", b"ACGTTACGT".to_vec())];
    let out = add_to_alignment(&empty, &q, &p, &NoProgress).unwrap();
    assert_eq!(out.len(), 2);
    assert!(out.is_aligned());

    // Empty query sequence.
    let out =
        add_to_alignment(&profile, &[Sequence::new("q", Vec::new())], &p, &NoProgress).unwrap();
    assert_eq!(out.len(), 3);
    assert!(out.is_aligned());
    assert_eq!(out.sequences[2].ungapped_len(), 0);
}

#[test]
fn explicit_thread_count_is_honoured_without_touching_the_global_pool() {
    let a = common::degapped(&common::simulate(16, 200, common::Rates::low(), b"ACGT", 3, 0xFEED));
    let mut p = AlignParams::for_engine(Engine::Muscle);
    p.threads = 2;
    let two = align(&a, &p, &NoProgress).expect("aligns on a private pool");
    p.threads = 0;
    let all = align(&a, &p, &NoProgress).expect("aligns on the global pool");
    // The result must not depend on how many threads did the work.
    assert_eq!(two.width(), all.width());
    for (x, y) in two.sequences.iter().zip(all.sequences.iter()) {
        assert_eq!(x.residues, y.residues);
    }
}

#[test]
fn engines_are_deterministic() {
    let a = common::degapped(&common::simulate(
        12,
        200,
        common::Rates::moderate(),
        b"ACDEFGHIKLMNPQRSTVWY",
        3,
        0x1234,
    ));
    for p in all_params() {
        let first = align(&a, &p, &NoProgress).unwrap();
        let second = align(&a, &p, &NoProgress).unwrap();
        for (x, y) in first.sequences.iter().zip(second.sequences.iter()) {
            assert_eq!(x.residues, y.residues, "{:?} is not deterministic", p.engine);
        }
    }
}

#[test]
fn unusual_residue_bytes_do_not_panic() {
    let a = aln(&[
        ("a", "ACGT*?XNBZJUO"),
        ("b", "ACGT*?XNBZJUO"),
        ("c", "ACGTXNBZ"),
        ("d", "!@#$%^&*()"),
    ]);
    for p in all_params() {
        let out = align(&a, &p, &NoProgress).expect("odd bytes are not an error");
        assert!(out.is_aligned());
        assert_residues_preserved(&a, &out);
    }
}
