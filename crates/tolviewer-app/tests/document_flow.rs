//! End-to-end flows through the document model, without opening a window.
//!
//! These exercise the paths a user actually takes — open a file, edit it,
//! align it, clean it, export it, read it back — across all four library
//! crates, which is where integration bugs live.

use std::path::PathBuf;

use tolviewer_align::{AlignParams, Engine, NoProgress};
use tolviewer_app::Document;
use tolviewer_clean::GblocksParams;
use tolviewer_core::{Alphabet, EditOp};
use tolviewer_io::{Format, WriteOptions};

fn example(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../testdata/examples").join(name)
}

fn open(name: &str) -> Document {
    let path = example(name);
    let alignment = tolviewer_io::read_file(&path)
        .unwrap_or_else(|e| panic!("could not read {}: {e}", path.display()));
    let format = Format::from_path(&path).expect("example files have known extensions");
    Document::new(alignment, Some(path), format)
}

#[test]
fn opens_the_dna_example_and_detects_dna() {
    let doc = open("tingidae_COI.fasta");
    assert_eq!(doc.rows(), 20);
    assert_eq!(doc.alphabet(), Alphabet::Dna);
    assert!(doc.width() > 600);
    assert_eq!(doc.alignment.sequences[0].id, "Tingis_cardui");
    assert_eq!(doc.alignment.sequences[0].description, "COI partial cds");
}

#[test]
fn opens_the_protein_example_and_detects_protein() {
    let doc = open("wingless_protein.fasta");
    assert_eq!(doc.rows(), 12);
    assert_eq!(doc.alphabet(), Alphabet::Protein);
}

#[test]
fn edit_then_undo_restores_the_file_byte_for_byte() {
    let mut doc = open("tingidae_COI.fasta");
    let before =
        tolviewer_io::write_string(&doc.alignment, Format::Fasta, &WriteOptions::default())
            .expect("FASTA always writes");

    doc.apply(EditOp::SetResidue { row: 3, col: 40, residue: b'N' }).unwrap();
    doc.apply(EditOp::RemoveSequence { row: 0 }).unwrap();
    doc.apply(EditOp::InsertGap { row: 1, col: 10 }).unwrap();
    assert_eq!(doc.rows(), 19);

    while doc.undo.can_undo() {
        doc.undo().unwrap();
    }
    let after = tolviewer_io::write_string(&doc.alignment, Format::Fasta, &WriteOptions::default())
        .expect("FASTA always writes");
    assert_eq!(before, after, "undoing every edit must restore the original file");
}

#[test]
fn aligning_makes_the_rows_rectangular_and_keeps_every_residue() {
    let mut doc = open("tingidae_COI.fasta");
    // The example has indels, so the raw file is not an alignment.
    assert!(!doc.alignment.is_aligned());

    let ungapped_before: Vec<Vec<u8>> =
        doc.alignment.sequences.iter().map(|s| s.ungapped()).collect();
    let names_before: Vec<String> = doc.alignment.sequences.iter().map(|s| s.id.clone()).collect();

    let params = AlignParams::for_engine(Engine::Clustal);
    let aligned = tolviewer_align::align(&doc.alignment, &params, &NoProgress).expect("alignment");
    doc.replace("align", aligned).unwrap();

    assert!(doc.alignment.is_aligned(), "the aligner must return rectangular rows");
    let names_after: Vec<String> = doc.alignment.sequences.iter().map(|s| s.id.clone()).collect();
    assert_eq!(names_before, names_after, "row order and names must survive alignment");
    for (i, before) in ungapped_before.iter().enumerate() {
        assert_eq!(
            &doc.alignment.sequences[i].ungapped(),
            before,
            "alignment changed the residues of row {i}"
        );
    }
}

#[test]
fn every_engine_produces_a_usable_alignment() {
    let doc = open("tingidae_COI.fasta");
    for engine in Engine::all() {
        let params = AlignParams::for_engine(*engine);
        let aligned = tolviewer_align::align(&doc.alignment, &params, &NoProgress)
            .unwrap_or_else(|e| panic!("{} failed: {e}", engine.name()));
        assert!(aligned.is_aligned(), "{} returned ragged rows", engine.name());
        assert_eq!(aligned.len(), doc.rows(), "{} lost sequences", engine.name());
        assert!(
            aligned.width()
                >= doc.alignment.sequences.iter().map(|s| s.ungapped_len()).max().unwrap(),
            "{} produced an alignment shorter than its longest sequence",
            engine.name()
        );
    }
}

#[test]
fn cleaning_an_alignment_only_ever_removes_columns() {
    let doc = open("tingidae_COI.fasta");
    let params = AlignParams::for_engine(Engine::Clustal);
    let aligned = tolviewer_align::align(&doc.alignment, &params, &NoProgress).expect("alignment");

    let gb = GblocksParams::defaults(aligned.len());
    let result = tolviewer_clean::gblocks(&aligned, &gb).expect("gblocks");
    assert_eq!(result.total, aligned.width());
    assert_eq!(result.mask.len(), aligned.width());
    assert_eq!(result.kept, result.mask.iter().filter(|&&k| k).count());

    let cleaned = result.apply(&aligned).expect("applying the mask");
    assert_eq!(cleaned.len(), aligned.len(), "cleaning must not drop sequences");
    assert_eq!(cleaned.width(), result.kept);
    assert!(cleaned.width() <= aligned.width());

    // The relaxed settings must never keep less than the strict ones.
    let relaxed = tolviewer_clean::gblocks(&aligned, &GblocksParams::relaxed(aligned.len()))
        .expect("gblocks relaxed");
    assert!(
        relaxed.kept >= result.kept,
        "relaxed settings kept {} columns, strict kept {}",
        relaxed.kept,
        result.kept
    );
}

#[test]
fn exports_round_trip_through_every_writable_format() {
    let doc = open("tingidae_COI.fasta");
    let params = AlignParams::for_engine(Engine::Clustal);
    let aligned = tolviewer_align::align(&doc.alignment, &params, &NoProgress).expect("alignment");

    for format in Format::all().iter().filter(|f| f.can_write()) {
        let text = tolviewer_io::write_string(&aligned, *format, &WriteOptions::default())
            .unwrap_or_else(|e| panic!("writing {} failed: {e}", format.name()));
        let parsed = tolviewer_io::parse(text.as_bytes(), *format, "round-trip")
            .unwrap_or_else(|e| panic!("re-reading {} failed: {e}", format.name()));

        assert_eq!(parsed.len(), aligned.len(), "{} lost sequences", format.name());
        assert_eq!(parsed.width(), aligned.width(), "{} changed the width", format.name());
        for (i, (a, b)) in aligned.sequences.iter().zip(&parsed.sequences).enumerate() {
            assert_eq!(
                a.residues.to_ascii_uppercase(),
                b.residues.to_ascii_uppercase(),
                "{} corrupted row {i} ({})",
                format.name(),
                a.id
            );
            // Strict PHYLIP is allowed to shorten names; nothing else is.
            if *format != Format::Phylip {
                assert_eq!(a.id, b.id, "{} changed a name", format.name());
            }
        }
    }
}

#[test]
fn deleting_selected_columns_shrinks_every_row_equally() {
    let mut doc = open("tingidae_COI.fasta");
    doc.alignment.pad_to_width();
    let width = doc.width();
    doc.apply(EditOp::DeleteColumns { start: 10, end: 40 }).unwrap();
    assert_eq!(doc.width(), width - 30);
    for seq in &doc.alignment.sequences {
        assert_eq!(seq.len(), width - 30);
    }
    doc.undo().unwrap();
    assert_eq!(doc.width(), width);
}

#[test]
fn writing_to_disk_is_atomic_and_readable() {
    let doc = open("tingidae_COI.fasta");
    let dir = std::env::temp_dir().join(format!("tolviewer-test-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("temp dir");
    let path = dir.join("out.fasta");

    tolviewer_io::write_file(&doc.alignment, &path, Format::Fasta, &WriteOptions::default())
        .expect("write");
    let reread = tolviewer_io::read_file(&path).expect("read back");
    assert_eq!(reread.len(), doc.rows());

    // No temporary files should be left behind.
    let leftovers: Vec<_> = std::fs::read_dir(&dir)
        .expect("listing")
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| n != "out.fasta")
        .collect();
    assert!(leftovers.is_empty(), "write left temporary files behind: {leftovers:?}");

    std::fs::remove_dir_all(&dir).ok();
}
