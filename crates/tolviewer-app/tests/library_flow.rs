//! End-to-end flows through the library, without opening a window.
//!
//! The library's own tests cover the tree and the save policy. What is tested
//! here is the glue the app adds on top: that a document opened from an entry
//! remembers where it came from, and that the chromatogram keeps pointing at
//! the right peaks as that document is edited, undone and redone — which is
//! the one piece of state that is derived rather than stored, and so the one
//! most likely to drift.

use std::path::{Path, PathBuf};

use tolviewer_app::chromatogram::{Link, TraceView};
use tolviewer_app::Document;
use tolviewer_core::EditOp;
use tolviewer_io::{Format, WriteOptions};
use tolviewer_library::{Library, SaveChoice, SaveTarget};

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("tolviewer-app-flow-{name}"));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn trace_file(dir: &Path) -> PathBuf {
    let source =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../testdata/ab1/tingidae_COI_F.ab1");
    let path = dir.join("TL-2213_COI_F.ab1");
    std::fs::copy(&source, &path).unwrap();
    path
}

/// Open a library entry the way the app does.
fn open(library: &Library, id: tolviewer_library::NodeId) -> Document {
    let entry = library.entry(id).expect("entry should exist");
    let alignment = entry.load().expect("entry should load");
    let trace = entry.load_trace().ok().map(|t| TraceView::new(t, 0));
    Document::new(alignment, Some(entry.source_path().to_path_buf()), entry.format)
        .from_library(id, trace)
}

#[test]
fn a_document_opened_from_a_trace_carries_its_chromatogram() {
    let dir = scratch("open-trace");
    let mut library = Library::new("p");
    let id = library.add_file(None, &trace_file(&dir)).unwrap();

    let mut doc = open(&library, id);
    assert_eq!(doc.origin, Some(id));
    assert_eq!(doc.rows(), 1);
    assert_eq!(doc.format, Format::Ab1);

    let view = doc.trace_view().expect("a trace document has a chromatogram");
    assert!(matches!(view.link(), Link::At { offset: 0, identity: 1.0 }));
    assert!(view.samples_are_present());
}

#[test]
fn the_chromatogram_follows_an_edit_and_comes_back_on_undo() {
    let dir = scratch("edit-undo");
    let mut library = Library::new("p");
    let id = library.add_file(None, &trace_file(&dir)).unwrap();
    let mut doc = open(&library, id);

    // What base 100 sits on before anything is touched.
    let peak_100 = doc.trace_view().unwrap().sample_for_residue(100);
    assert!(peak_100.is_some());

    // Retype one call: the link holds, and the panel can tell that the row and
    // the instrument now disagree there.
    doc.apply(EditOp::SetResidue { row: 0, col: 100, residue: b'A' }).unwrap();
    let view = doc.trace_view().unwrap();
    assert!(matches!(view.link(), Link::At { offset: 0, .. }));
    assert_eq!(view.sample_for_residue(100), peak_100, "one edit must not move the peaks");
    assert_eq!(doc.alignment.sequences[0].residues[100], b'A');

    // Trim the first 40 bases. Every remaining base keeps its own peak.
    let mut trimmed = doc.alignment.clone();
    trimmed.sequences[0].residues = trimmed.sequences[0].residues[40..].to_vec();
    trimmed.sequences[0].quality = None;
    doc.replace("trim", trimmed).unwrap();
    let view = doc.trace_view().unwrap();
    assert!(matches!(view.link(), Link::At { offset: 40, .. }), "{:?}", view.link());
    assert_eq!(view.sample_for_residue(60), peak_100, "base 100 is base 60 after a 40-base trim");

    // Undo, and the link goes back on its own.
    doc.undo().unwrap();
    let view = doc.trace_view().unwrap();
    assert!(matches!(view.link(), Link::At { offset: 0, .. }));
    assert_eq!(view.sample_for_residue(100), peak_100);

    // Redo, and it follows again.
    doc.redo().unwrap();
    assert!(matches!(doc.trace_view().unwrap().link(), Link::At { offset: 40, .. }));
}

#[test]
fn aligning_a_trace_row_does_not_break_the_chromatogram() {
    let dir = scratch("gapped");
    let mut library = Library::new("p");
    let id = library.add_file(None, &trace_file(&dir)).unwrap();
    let mut doc = open(&library, id);
    let peak_50 = doc.trace_view().unwrap().sample_for_residue(50);

    // Gaps are what aligning a read against others puts in. They occupy a
    // column but not a call, so residue 50 is still call 50.
    doc.apply(EditOp::InsertGap { row: 0, col: 10 }).unwrap();
    doc.apply(EditOp::InsertGap { row: 0, col: 30 }).unwrap();
    let view = doc.trace_view().unwrap();
    assert!(matches!(view.link(), Link::At { offset: 0, identity: 1.0 }));
    assert_eq!(view.sample_for_residue(50), peak_50);
}

#[test]
fn a_trace_document_saves_through_the_library_and_never_over_the_ab1() {
    let dir = scratch("save-trace");
    let path = trace_file(&dir);
    let untouched = std::fs::read(&path).unwrap();
    let mut library = Library::new("p");
    let id = library.add_file(None, &path).unwrap();
    let mut doc = open(&library, id);

    doc.apply(EditOp::SetResidue { row: 0, col: 30, residue: b'A' }).unwrap();
    assert!(doc.is_dirty());

    // An .ab1 cannot be written, so a copy is the only option — the app must
    // not be able to offer "replace the original" at all.
    let target = library.save_target(id).unwrap();
    assert!(target.needs_confirmation());
    assert!(!target.can_overwrite());
    assert!(matches!(target, SaveTarget::MustCopy(..)));
    assert_eq!(target.path().extension().unwrap(), "fasta");

    let copy = target.path().to_path_buf();
    library
        .save_entry(id, &doc.alignment, SaveChoice::NewCopy(copy.clone()), &WriteOptions::default())
        .unwrap();
    doc.mark_saved();

    assert_eq!(std::fs::read(&path).unwrap(), untouched, "the trace file was written to");
    assert!(!doc.is_dirty());
    // Reopening now reads the copy, and the chromatogram still comes from the
    // .ab1, because that is where the signal lives.
    let reopened = open(&library, id);
    assert_eq!(reopened.alignment.sequences[0].residues[30], b'A');
    assert!(reopened.trace.is_none(), "the entry is a FASTA copy now, not a trace");
}

#[test]
fn a_plain_alignment_document_has_no_chromatogram_and_says_so() {
    let dir = scratch("no-trace");
    let path = dir.join("aln.fasta");
    std::fs::write(&path, ">a\nACGT\n>b\nACGA\n").unwrap();
    let mut library = Library::new("p");
    let id = library.add_file(None, &path).unwrap();
    let mut doc = open(&library, id);
    assert!(doc.trace_view().is_none());
    assert_eq!(doc.origin, Some(id));
    // And it can be written back in place, unlike a trace.
    assert!(library.save_target(id).unwrap().can_overwrite());
}
