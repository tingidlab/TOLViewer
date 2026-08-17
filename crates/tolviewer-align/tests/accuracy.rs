//! Accuracy of the three engines against known reference alignments.
//!
//! The method is the standard one: take an alignment whose columns are known
//! to be correct, throw the gaps away, re-align, and measure how much of the
//! reference came back. Two measures are reported:
//!
//! * **SP** - the fraction of residue *pairs* the reference aligns that the
//!   test alignment also aligns. Forgiving, and the usual headline number.
//! * **TC** - the fraction of reference *columns* reproduced exactly. Harsh: a
//!   single misplaced residue loses the whole column.
//!
//! The thresholds below are the measured numbers rounded down with a little
//! headroom, not aspirations. Every one of them was read off a run of this
//! file; see the crate summary for the values.

mod common;

use common::{column_accuracy, degapped, mean_identity, reference, simulate, sp_accuracy, Rates};
use tolviewer_align::{align, AlignParams, Engine, NoProgress};
use tolviewer_core::Alignment;

/// Re-align a reference and return `(SP, TC)`.
fn score(reference: &Alignment, engine: Engine) -> (f32, f32) {
    let input = degapped(reference);
    let params = AlignParams::for_engine(engine);
    let out = align(&input, &params, &NoProgress).expect("alignment succeeds");
    assert!(out.is_aligned(), "{} returned a ragged alignment", engine.name());
    assert_eq!(out.len(), reference.len());
    for (i, s) in out.sequences.iter().enumerate() {
        assert_eq!(
            s.ungapped(),
            input.sequences[i].residues,
            "{} changed the residues of row {i}",
            engine.name()
        );
    }
    (sp_accuracy(reference, &out), column_accuracy(reference, &out))
}

fn report(name: &str, reference: &Alignment) -> Vec<(Engine, f32, f32)> {
    let mut out = Vec::new();
    for &e in Engine::all() {
        let (sp, tc) = score(reference, e);
        println!(
            "{name:<28} {:<8} identity {:.3}  SP {:.4}  TC {:.4}",
            e.name(),
            mean_identity(reference),
            sp,
            tc
        );
        out.push((e, sp, tc));
    }
    out
}

/// Build one reference row from a base sequence by replacing the given ranges
/// with gaps and applying point substitutions. Deriving the rows this way -
/// rather than typing them out - guarantees the reference really is the
/// alignment the edits imply, so an aligner that disagrees with it is wrong.
fn derived(base: &str, gaps: &[(usize, usize)], subs: &[(usize, u8)]) -> String {
    let mut row: Vec<u8> = base.as_bytes().to_vec();
    for &(pos, c) in subs {
        row[pos] = c;
    }
    for &(start, len) in gaps {
        for c in row.iter_mut().skip(start).take(len) {
            *c = b'-';
        }
    }
    String::from_utf8(row).expect("ASCII")
}

const DNA_BASE: &str = "ATGCGTACGTTAGCCCGATTACAGGTATCGCTAGCTAGCATCGATCGTTA";
const PROTEIN_BASE: &str = "MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLRWVQT";

/// A hand-built five-sequence DNA case with two unambiguous indels.
fn dna_case() -> Alignment {
    reference(&[
        ("alpha", &derived(DNA_BASE, &[], &[])),
        ("beta", &derived(DNA_BASE, &[], &[(21, b'T')])),
        ("gamma", &derived(DNA_BASE, &[(13, 3)], &[])),
        ("delta", &derived(DNA_BASE, &[(30, 6)], &[])),
        ("eps", &derived(DNA_BASE, &[], &[(6, b'C')])),
    ])
}

/// A protein case: a conserved core with two deleted loops.
fn protein_case() -> Alignment {
    reference(&[
        ("p1", &derived(PROTEIN_BASE, &[], &[])),
        ("p2", &derived(PROTEIN_BASE, &[], &[(1, b'R')])),
        ("p3", &derived(PROTEIN_BASE, &[(26, 6)], &[])),
        ("p4", &derived(PROTEIN_BASE, &[(49, 6)], &[])),
        ("p5", &derived(PROTEIN_BASE, &[], &[(39, b'A')])),
        ("p6", &derived(PROTEIN_BASE, &[(26, 6)], &[(12, b'M')])),
    ])
}

#[test]
fn hand_built_dna_case() {
    let r = dna_case();
    for (engine, sp, tc) in report("hand-built DNA", &r) {
        assert!(sp >= 0.98, "{} SP {sp:.4} below 0.98", engine.name());
        assert!(tc >= 0.95, "{} TC {tc:.4} below 0.95", engine.name());
    }
}

#[test]
fn hand_built_protein_case() {
    let r = protein_case();
    for (engine, sp, tc) in report("hand-built protein", &r) {
        assert!(sp >= 0.98, "{} SP {sp:.4} below 0.98", engine.name());
        assert!(tc >= 0.95, "{} TC {tc:.4} below 0.95", engine.name());
    }
}

#[test]
fn simulated_dna_low_divergence() {
    let r = simulate(12, 300, Rates::low(), b"ACGT", 3, 0xA11CE);
    // Measured: SP 0.943-0.969, TC 0.808-0.864.
    for (engine, sp, tc) in report("simulated DNA, low", &r) {
        assert!(sp >= 0.90, "{} SP {sp:.4} below 0.90", engine.name());
        assert!(tc >= 0.75, "{} TC {tc:.4} below 0.75", engine.name());
    }
}

#[test]
fn simulated_protein_low_divergence() {
    let r = simulate(12, 300, Rates::low(), b"ACDEFGHIKLMNPQRSTVWY", 3, 0xB0B);
    // Measured: SP 0.970-0.974, TC 0.880-0.887.
    for (engine, sp, tc) in report("simulated protein, low", &r) {
        assert!(sp >= 0.94, "{} SP {sp:.4} below 0.94", engine.name());
        assert!(tc >= 0.82, "{} TC {tc:.4} below 0.82", engine.name());
    }
}

#[test]
fn simulated_dna_moderate_divergence() {
    let r = simulate(12, 300, Rates::moderate(), b"ACGT", 3, 0xC0FFEE);
    // 69 % mean identity with frequent indels; measured SP 0.648-0.681.
    for (engine, sp, _tc) in report("simulated DNA, moderate", &r) {
        assert!(sp >= 0.60, "{} SP {sp:.4} below 0.60", engine.name());
    }
}

#[test]
fn simulated_protein_moderate_divergence() {
    let r = simulate(12, 300, Rates::moderate(), b"ACDEFGHIKLMNPQRSTVWY", 3, 0xD15EA5E);
    // 62 % mean identity with frequent indels; measured SP 0.707 (Clustal,
    // a single progressive pass) to 0.770 (MUSCLE).
    for (engine, sp, _tc) in report("simulated protein, moderate", &r) {
        assert!(sp >= 0.68, "{} SP {sp:.4} below 0.68", engine.name());
    }
}

/// Averaged over a spread of simulated cases, the two refining engines must be
/// at least as good as the single-pass baseline. This is the property the
/// architecture document asks for, and it is what iteration is *for*.
#[test]
fn refining_engines_are_no_worse_than_clustal_on_average() {
    let cases: Vec<Alignment> = vec![
        simulate(10, 250, Rates::low(), b"ACGT", 3, 1),
        simulate(10, 250, Rates::moderate(), b"ACGT", 3, 2),
        simulate(14, 200, Rates::low(), b"ACDEFGHIKLMNPQRSTVWY", 3, 3),
        simulate(14, 200, Rates::moderate(), b"ACDEFGHIKLMNPQRSTVWY", 3, 4),
        simulate(8, 400, Rates::moderate(), b"ACGT", 3, 5),
        simulate(8, 400, Rates::moderate(), b"ACDEFGHIKLMNPQRSTVWY", 3, 6),
    ];
    let mut totals = [0.0f32; 3];
    for (n, case) in cases.iter().enumerate() {
        for (k, &engine) in Engine::all().iter().enumerate() {
            let mut p = AlignParams::for_engine(engine);
            p.iterations = p.iterations.max(2);
            let input = degapped(case);
            let out = align(&input, &p, &NoProgress).expect("alignment succeeds");
            let sp = sp_accuracy(case, &out);
            println!("case {n} {:<8} SP {sp:.4}", engine.name());
            totals[k] += sp;
        }
    }
    let mean: Vec<f32> = totals.iter().map(|t| t / cases.len() as f32).collect();
    println!("mean SP: Clustal {:.4}  MUSCLE {:.4}  MAFFT {:.4}", mean[0], mean[1], mean[2]);
    assert!(
        mean[1] >= mean[0] - 0.005,
        "MUSCLE mean SP {:.4} worse than Clustal {:.4}",
        mean[1],
        mean[0]
    );
    assert!(
        mean[2] >= mean[0] - 0.005,
        "MAFFT mean SP {:.4} worse than Clustal {:.4}",
        mean[2],
        mean[0]
    );
}

/// Identical inputs must come back identical, with no gaps invented.
#[test]
fn duplicate_sequences_align_to_themselves() {
    let r = reference(&[
        ("a", "ACGTACGTACGTACGT"),
        ("b", "ACGTACGTACGTACGT"),
        ("c", "ACGTACGTACGTACGT"),
    ]);
    for &engine in Engine::all() {
        let out = align(&degapped(&r), &AlignParams::for_engine(engine), &NoProgress).unwrap();
        assert_eq!(out.width(), 16, "{} invented columns", engine.name());
    }
}

/// A 200 x ~1000 alignment must finish. Marked `#[ignore]` because it is slow
/// in a debug build; run with `cargo test --release -- --ignored`.
#[test]
#[ignore = "performance check; run with --release --ignored"]
fn two_hundred_sequences_of_a_thousand_columns() {
    let r = simulate(200, 1000, Rates::low(), b"ACGT", 8, 0x5EED);
    let input = degapped(&r);
    assert_eq!(input.len(), 200);
    for &engine in Engine::all() {
        let start = std::time::Instant::now();
        let out = align(&input, &AlignParams::for_engine(engine), &NoProgress)
            .expect("large alignment succeeds");
        let elapsed = start.elapsed();
        assert!(out.is_aligned());
        assert_eq!(out.len(), 200);
        println!(
            "{:<8} 200 x {} columns in {:.2} s (SP {:.4})",
            engine.name(),
            out.width(),
            elapsed.as_secs_f64(),
            sp_accuracy(&r, &out)
        );
    }
}

/// Long sequences: the case the linear-space pairwise path and the FFT band
/// exist for. Two 20 kb "organelle genomes" would need a 400-million-cell
/// traceback if aligned the naive way.
#[test]
#[ignore = "performance check; run with --release --ignored"]
fn long_sequence_pair() {
    let r = simulate(2, 20_000, Rates::low(), b"ACGT", 1, 0x0DDBA11);
    let input = degapped(&r);
    for &engine in Engine::all() {
        let start = std::time::Instant::now();
        let out =
            align(&input, &AlignParams::for_engine(engine), &NoProgress).expect("long pair aligns");
        assert!(out.is_aligned());
        println!(
            "{:<8} 2 x ~20 kb in {:.2} s (width {}, SP {:.4})",
            engine.name(),
            start.elapsed().as_secs_f64(),
            out.width(),
            sp_accuracy(&r, &out)
        );
    }
}

/// A handful of long sequences: this is where the FFT band actually engages,
/// because the profiles are wide enough to exceed the banding threshold.
#[test]
#[ignore = "performance check; run with --release --ignored"]
fn several_long_sequences() {
    let r = simulate(6, 12_000, Rates::low(), b"ACGT", 3, 0x0C7013);
    let input = degapped(&r);
    for &engine in Engine::all() {
        let start = std::time::Instant::now();
        let out = align(&input, &AlignParams::for_engine(engine), &NoProgress)
            .expect("long sequences align");
        assert!(out.is_aligned());
        println!(
            "{:<8} 6 x ~12 kb in {:.2} s (width {}, SP {:.4})",
            engine.name(),
            start.elapsed().as_secs_f64(),
            out.width(),
            sp_accuracy(&r, &out)
        );
    }
}
