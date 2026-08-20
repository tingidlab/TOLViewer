//! The whole thing, once: the workflow a lace bug project actually goes
//! through, from a plate of traces to a concatenated matrix.
//!
//! Each step is a real operation on real files in a scratch directory, so this
//! catches the things unit tests do not: that the pieces compose, that the
//! files on disk are what the library says they are, and that nothing the user
//! did not ask for got written.

use std::path::{Path, PathBuf};

use tolviewer_align::{AlignParams, Engine, NoProgress};
use tolviewer_core::{Alignment, Sequence};
use tolviewer_io::{Format, WriteOptions};
use tolviewer_library::{
    concat, primer, store, ConcatOptions, EntryKind, Library, NodeId, Primer, SaveChoice,
    SaveTarget, TrimOptions,
};

const LCO: &str = "GGTCAACAAATCATAAAGATATTGG";
const HCO: &str = "TAAACTTCAGGGTGACCAAAAAATCA";

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("tolviewer-workflow-{name}"));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn rc(s: &str) -> String {
    s.bytes().rev().map(|c| tolviewer_core::Alphabet::Dna.complement(c) as char).collect()
}

fn write_fasta(path: &Path, rows: &[(&str, &str)]) {
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    let text: String = rows.iter().map(|(id, seq)| format!(">{id}\n{seq}\n")).collect();
    std::fs::write(path, text).unwrap();
}

/// A read as it comes off the machine: junk, forward primer, insert, the
/// reverse primer's complement, more junk.
fn amplicon(insert: &str) -> String {
    format!("NNCTGA{LCO}{insert}{}GGNTAC", rc(HCO))
}

/// Four specimens' worth of 18S, each with the primers still on.
fn inserts() -> [(&'static str, String); 4] {
    let core = "ACGTTGGCCATTGGCCAATTGGCCAATTGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAA";
    [
        ("TL-2213", core.to_string()),
        ("TL-2214", core.replacen("ACGT", "ACGA", 1)),
        ("TL-2215", core.replacen("ACGT", "ACGC", 1)),
        ("TL-2216", core.replacen("ACGT", "ACGG", 1)),
    ]
}

#[test]
fn a_lace_bug_project_from_reads_to_a_supermatrix() {
    let dir = scratch("full");

    // ---- the tree the user arranges ------------------------------------
    let mut library = Library::new("Lace bug project");
    let project = library.add_folder(None, "Lace bug project").unwrap();
    let ssu = library.add_folder(Some(project), "18S").unwrap();
    let lsu = library.add_folder(Some(project), "28S").unwrap();
    assert_eq!(library.path_of(ssu), "Lace bug project / 18S");

    // ---- files stay where the facility put them ------------------------
    let reads = dir.join("reads");
    let mut originals: Vec<(PathBuf, Vec<u8>)> = Vec::new();
    let mut read_ids: Vec<NodeId> = Vec::new();
    for (sample, insert) in inserts() {
        let path = reads.join(format!("{sample}_18S_F.fasta"));
        write_fasta(&path, &[(&format!("{sample}_18S_F"), &amplicon(&insert))]);
        originals.push((path.clone(), std::fs::read(&path).unwrap()));
        read_ids.push(library.add_file(Some(ssu), &path).unwrap());
    }
    // One specimen was sequenced from the reverse primer, so its file holds the
    // other strand.
    let reverse_path = reads.join("TL-2216_18S_R.fasta");
    let (_, forward_2216) = inserts()[3].clone();
    write_fasta(&reverse_path, &[("TL-2216_18S_R", &rc(&amplicon(&forward_2216)))]);
    originals.push((reverse_path.clone(), std::fs::read(&reverse_path).unwrap()));
    let reverse_read = library.add_file(Some(ssu), &reverse_path).unwrap();

    // A real trace, with a chromatogram behind it.
    let trace_source =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../testdata/ab1/tingidae_COI_F.ab1");
    let trace_path = reads.join("TL-2213_COI_F.ab1");
    std::fs::copy(&trace_source, &trace_path).unwrap();
    let coi = library.add_folder(Some(project), "COI").unwrap();
    let trace_entry = library.add_file(Some(coi), &trace_path).unwrap();
    assert_eq!(library.entry(trace_entry).unwrap().kind, EntryKind::Trace);
    assert_eq!(library.get(trace_entry).unwrap().name, "TL-2213_COI_F");
    assert!(library.entry(trace_entry).unwrap().load_trace().unwrap().samples() > 1000);

    // ---- reversing is a flag, not a rewrite ----------------------------
    let as_read = library.load(reverse_read).unwrap().sequences[0].residues.clone();
    library.set_reversed(reverse_read, true).unwrap();
    let flipped = library.load(reverse_read).unwrap().sequences[0].residues.clone();
    assert_ne!(as_read, flipped);
    assert_eq!(flipped, amplicon(&forward_2216).into_bytes());

    // ---- primers map, and trimming takes them off ----------------------
    library.primers.push(Primer::new("LCO1490", LCO).unwrap());
    library.primers.push(Primer::new("HCO2198", HCO).unwrap());
    library.touch();

    let read = library.load(read_ids[0]).unwrap();
    let residues = read.sequences[0].ungapped();
    let hits = library.primers.map(&residues, 0.1);
    assert_eq!(hits.len(), 2, "both primers should bind: {hits:?}");

    let plan = primer::plan_trim(&library.primers, &residues, &TrimOptions::default());
    assert!(plan.trims_anything());
    assert_eq!(&residues[plan.range.clone()], inserts()[0].1.as_bytes());

    // Trimming every read gives the inserts, which is what gets aligned.
    let mut trimmed_18s: Vec<Sequence> = Vec::new();
    for &id in read_ids.iter().take(3).chain([&reverse_read]) {
        let loaded = library.load(id).unwrap();
        for seq in &loaded.sequences {
            let residues = seq.ungapped();
            let plan = primer::plan_trim(&library.primers, &residues, &TrimOptions::default());
            assert!(plan.trims_anything(), "{} kept its primers", seq.id);
            trimmed_18s.push(Sequence::new(seq.id.clone(), residues[plan.range].to_vec()));
        }
    }
    assert_eq!(trimmed_18s.len(), 4);

    // ---- select several sequences and align them -----------------------
    let params = AlignParams { engine: Engine::Clustal, ..AlignParams::default() };
    let aligned =
        tolviewer_align::align(&Alignment::new("18S", trimmed_18s), &params, &NoProgress).unwrap();
    assert!(aligned.is_aligned());
    assert_eq!(aligned.len(), 4);

    let ssu_path = dir.join("18S.fasta");
    tolviewer_io::write_file(&aligned, &ssu_path, Format::Fasta, &WriteOptions::default()).unwrap();
    let ssu_entry = library.add_file(Some(ssu), &ssu_path).unwrap();
    assert_eq!(library.entry(ssu_entry).unwrap().kind, EntryKind::Alignment);

    // ---- extract one sequence out of that alignment --------------------
    let extracted = library
        .add_selection(Some(ssu), ssu_entry, vec!["TL-2213_18S_F".to_string()], "TL-2213 18S")
        .unwrap();
    let one = library.load(extracted).unwrap();
    assert_eq!(one.len(), 1);
    assert_eq!(one.sequences[0].id, "TL-2213_18S_F");

    // It cannot be written back over the alignment it is part of.
    let target = library.save_target(extracted).unwrap();
    assert!(!target.can_overwrite());
    assert!(library
        .save_entry(extracted, &one, SaveChoice::Overwrite, &WriteOptions::default())
        .is_err());

    // ---- a second locus, and a supermatrix -----------------------------
    let lsu_path = dir.join("28S.fasta");
    write_fasta(
        &lsu_path,
        &[
            ("TL-2213_28S", "TTTTGGGGAAAACCCC"),
            ("TL-2214_28S", "TTTTGGGGAAAACCCA"),
            ("TL-2215_28S", "TTTTGGGGAAAACCCG"),
            // TL-2216 was not sequenced for 28S.
            ("TL-9999_28S", "TTTTGGGGAAAACCCT"),
        ],
    );
    let lsu_entry = library.add_file(Some(lsu), &lsu_path).unwrap();

    let mut a = library.load(ssu_entry).unwrap();
    let mut b = library.load(lsu_entry).unwrap();
    a.name = "18S".into();
    b.name = "28S".into();
    let result = concat::concatenate(&[&a, &b], &ConcatOptions::default()).unwrap();

    // Four specimens from 18S plus the one only 28S has.
    assert_eq!(result.alignment.len(), 5);
    assert_eq!(result.complete, 3, "TL-2216 has no 28S and TL-9999 has no 18S");
    assert_eq!(result.alignment.width(), a.width() + b.width());
    assert!(result.alignment.is_aligned());
    assert_eq!(result.partitions[1].range, a.width()..a.width() + b.width());
    assert!(result.nexus_charsets().contains("charset p_18S"));

    // The reverse read matched its own specimen despite being named _R.
    let missing_names: Vec<&str> = result.missing.iter().map(|m| m.sample.as_str()).collect();
    assert!(missing_names.contains(&"TL-2216_18S_R"), "{missing_names:?}");
    assert!(missing_names.contains(&"TL-9999_28S"), "{missing_names:?}");

    // ---- save, reopen, and find everything where it was ----------------
    let library_path = dir.join("lace-bugs.tolvlib");
    store::save(&mut library, &library_path).unwrap();
    assert!(!library.is_dirty());

    let reopened = store::load(&library_path).unwrap();
    assert_eq!(reopened.name, "Lace bug project");
    assert_eq!(reopened.walk().len(), library.walk().len());
    assert_eq!(reopened.primers.len(), 2);
    assert!(reopened.broken_entries().is_empty());
    // The orientation flag and the extract's row selection came back.
    let reopened_entries = reopened.entries_under(None);
    assert!(reopened_entries.iter().any(|&id| reopened.entry(id).unwrap().reversed));
    assert!(reopened_entries.iter().any(|&id| reopened.entry(id).unwrap().select.is_some()));

    // ---- and nothing wrote to the lab's files --------------------------
    for (path, before) in &originals {
        assert_eq!(&std::fs::read(path).unwrap(), before, "{} was modified", path.display());
    }
    assert_eq!(
        std::fs::read(&trace_path).unwrap(),
        std::fs::read(&trace_source).unwrap(),
        "the trace was modified"
    );
}

#[test]
fn editing_a_read_asks_once_and_then_keeps_one_copy() {
    let dir = scratch("save-policy");
    let original = dir.join("TL-2213_18S_F.fasta");
    write_fasta(&original, &[("TL-2213_18S_F", &amplicon(&inserts()[0].1))]);
    let untouched = std::fs::read(&original).unwrap();

    let mut library = Library::new("p");
    let id = library.add_file(None, &original).unwrap();

    // The first save would land on the lab's file, so it has to be asked about.
    let target = library.save_target(id).unwrap();
    assert!(target.needs_confirmation());
    assert_eq!(target, SaveTarget::Original(original.clone()));

    // The user says "save a copy".
    let copy = library.entry(id).unwrap().suggested_copy();
    let mut edited = library.load(id).unwrap();
    edited.sequences[0].residues.truncate(80);
    let written = library
        .save_entry(id, &edited, SaveChoice::NewCopy(copy.clone()), &WriteOptions::default())
        .unwrap();
    assert_eq!(written, copy);
    assert_eq!(std::fs::read(&original).unwrap(), untouched);

    // Every later edit goes to that copy without asking, and without making
    // another one.
    for length in [70usize, 60, 50] {
        assert!(!library.save_target(id).unwrap().needs_confirmation());
        let mut again = library.load(id).unwrap();
        again.sequences[0].residues.truncate(length);
        let path = library
            .save_entry(id, &again, SaveChoice::Overwrite, &WriteOptions::default())
            .unwrap();
        assert_eq!(path, copy);
    }
    let files: Vec<String> = std::fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .collect();
    assert_eq!(files.len(), 2, "four edits should leave one extra file: {files:?}");
    assert_eq!(library.load(id).unwrap().sequences[0].len(), 50);
    assert_eq!(std::fs::read(&original).unwrap(), untouched);
}

#[test]
fn the_user_may_still_choose_to_replace_the_original() {
    let dir = scratch("overwrite");
    let original = dir.join("read.fasta");
    write_fasta(&original, &[("a", "ACGTACGT")]);
    let mut library = Library::new("p");
    let id = library.add_file(None, &original).unwrap();

    let mut edited = library.load(id).unwrap();
    edited.sequences[0].residues = b"TTTT".to_vec();
    library.save_entry(id, &edited, SaveChoice::Overwrite, &WriteOptions::default()).unwrap();

    assert!(std::fs::read_to_string(&original).unwrap().contains("TTTT"));
    assert_eq!(std::fs::read_dir(&dir).unwrap().count(), 1, "no copy should have been made");
    // Still the original, so the next save asks again rather than assuming.
    assert!(library.save_target(id).unwrap().needs_confirmation());
}

#[test]
fn a_library_moved_to_another_machine_still_finds_its_files() {
    let here = scratch("portable-here");
    let there = scratch("portable-there");
    write_fasta(&here.join("reads/a.fasta"), &[("a", "ACGT")]);
    write_fasta(&here.join("reads/b.fasta"), &[("b", "ACGA")]);

    let mut library = Library::new("p");
    let folder = library.add_folder(None, "18S").unwrap();
    library.add_file(Some(folder), &here.join("reads/a.fasta")).unwrap();
    library.add_file(Some(folder), &here.join("reads/b.fasta")).unwrap();
    let path = here.join("p.tolvlib");
    store::save(&mut library, &path).unwrap();

    // Copy the project directory wholesale, as a shared drive or a zip would.
    std::fs::create_dir_all(there.join("reads")).unwrap();
    for name in ["reads/a.fasta", "reads/b.fasta", "p.tolvlib"] {
        std::fs::copy(here.join(name), there.join(name)).unwrap();
    }
    std::fs::remove_dir_all(&here).unwrap();

    let moved = store::load(&there.join("p.tolvlib")).unwrap();
    assert!(moved.broken_entries().is_empty(), "the moved project did not resolve");
    let (all, failed) = moved.gather(moved.roots(), "everything");
    assert!(failed.is_empty());
    assert_eq!(all.len(), 2);
}
