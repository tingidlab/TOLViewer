//! Every file in `testdata/io` must be detected and read correctly.

use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};

use tolviewer_io::{read_file, Format};

fn fixture_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../testdata/io")
}

/// file name, format `sniff` must report, rows, columns.
const EXPECTED: &[(&str, Format, usize, usize)] = &[
    ("simple.fasta", Format::Fasta, 3, 24),
    ("reads.fastq", Format::Fastq, 3, 12),
    ("interleaved.phy", Format::Phylip, 4, 20),
    ("strict_names.phy", Format::Phylip, 2, 30),
    ("relaxed.phy", Format::PhylipRelaxed, 3, 24),
    ("simple.nex", Format::Nexus, 3, 16),
    ("interleaved.nex", Format::Nexus, 3, 20),
    ("conserved.aln", Format::Clustal, 3, 28),
    ("pfam.sto", Format::Stockholm, 2, 20),
    ("align.msf", Format::Msf, 3, 16),
    ("records.gb", Format::Genbank, 2, 24),
];

#[test]
fn every_fixture_is_listed_in_the_table() {
    let listed: HashSet<&str> = EXPECTED.iter().map(|e| e.0).collect();
    let mut found = 0;
    for entry in fs::read_dir(fixture_dir()).expect("testdata/io must exist") {
        let path = entry.unwrap().path();
        if !path.is_file() {
            continue;
        }
        let name = path.file_name().unwrap().to_string_lossy().into_owned();
        assert!(listed.contains(name.as_str()), "{name} is not in EXPECTED");
        found += 1;
    }
    assert_eq!(found, EXPECTED.len(), "a fixture in EXPECTED is missing");
}

#[test]
fn fixtures_sniff_and_read_with_the_expected_dimensions() {
    for &(name, format, rows, cols) in EXPECTED {
        let path = fixture_dir().join(name);
        let bytes = fs::read(&path).unwrap();
        assert_eq!(Format::sniff(&bytes), Some(format), "sniffing {name}");
        let aln = read_file(&path).unwrap_or_else(|e| panic!("reading {name}: {e}"));
        assert_eq!(aln.len(), rows, "rows of {name}");
        assert_eq!(aln.width(), cols, "columns of {name}");
        assert_eq!(aln.name, Path::new(name).file_stem().unwrap().to_str().unwrap());
    }
}

#[test]
fn fasta_keeps_descriptions_and_case() {
    let aln = read_file(&fixture_dir().join("simple.fasta")).unwrap();
    assert_eq!(aln.sequences[0].id, "alpha");
    assert_eq!(aln.sequences[0].description, "Homo sapiens COI partial cds");
    // lowercase masking is preserved
    assert!(aln.sequences[2].residues.starts_with(b"acgt"));
}

#[test]
fn fastq_quality_is_normalised_to_phred_scores() {
    let aln = read_file(&fixture_dir().join("reads.fastq")).unwrap();
    let q = aln.sequences[0].quality.as_ref().unwrap();
    assert_eq!(q, &vec![40u8; 12]);
    assert_eq!(aln.sequences[2].quality.as_ref().unwrap()[0], 0);
}

#[test]
fn strict_phylip_names_may_contain_spaces() {
    let aln = read_file(&fixture_dir().join("strict_names.phy")).unwrap();
    assert_eq!(aln.sequences[0].id, "Homo sapi");
    assert_eq!(aln.sequences[1].id, "Pan trogl");
}

#[test]
fn interleaved_phylip_blocks_are_concatenated_in_order() {
    let aln = read_file(&fixture_dir().join("interleaved.phy")).unwrap();
    assert_eq!(aln.sequences[2].id, "Seq_C");
    assert_eq!(aln.sequences[2].residues, b"ACGT--GTACGTACGTACGT");
    assert_eq!(aln.sequences[3].residues, b"ACGTACGTACGTACGT--GT");
}

#[test]
fn nexus_matchchar_and_quoted_names() {
    let mut aln = read_file(&fixture_dir().join("simple.nex")).unwrap();
    assert_eq!(aln.sequences[2].id, "Homo sapiens");
    // `.` copies from the first row, `-` stays a gap.
    assert_eq!(aln.sequences[1].residues, b"ACGTACGT----ACGT");
    assert_eq!(aln.alphabet(), tolviewer_core::Alphabet::Dna);
}

#[test]
fn nexus_datatype_sets_the_alphabet() {
    let mut aln = read_file(&fixture_dir().join("interleaved.nex")).unwrap();
    assert_eq!(aln.alphabet(), tolviewer_core::Alphabet::Protein);
    assert_eq!(aln.sequences[2].residues, b"MKVL--PQRSTYACDEFGHI");
}

#[test]
fn clustal_conservation_lines_are_not_sequences() {
    let aln = read_file(&fixture_dir().join("conserved.aln")).unwrap();
    let ids: Vec<&str> = aln.sequences.iter().map(|s| s.id.as_str()).collect();
    assert_eq!(ids, vec!["alpha", "beta", "gamma"]);
    assert_eq!(aln.sequences[1].residues, b"ACGTACGTACGTACGT--GTGGGGCCCC");
}

#[test]
fn stockholm_annotation_is_skipped() {
    let aln = read_file(&fixture_dir().join("pfam.sto")).unwrap();
    let ids: Vec<&str> = aln.sequences.iter().map(|s| s.id.as_str()).collect();
    assert_eq!(ids, vec!["alpha/1-20", "beta/1-20"]);
    assert_eq!(aln.sequences[1].residues, b"ACGUACGU--GUAAAAUUUU");
}

#[test]
fn msf_dots_are_gaps() {
    let aln = read_file(&fixture_dir().join("align.msf")).unwrap();
    assert_eq!(aln.sequences[1].residues, b"MKVLWIPQRSTY--DE");
}

#[test]
fn genbank_records_become_sequences() {
    let aln = read_file(&fixture_dir().join("records.gb")).unwrap();
    assert_eq!(aln.sequences[0].id, "AB000001");
    assert_eq!(
        aln.sequences[0].description,
        "Example species mitochondrial cytochrome oxidase subunit I (COI) gene, partial cds"
    );
    assert_eq!(aln.sequences[1].residues, b"ggccggccggcc");
}
