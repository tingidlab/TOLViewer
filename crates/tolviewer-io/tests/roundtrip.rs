//! Write then re-read: names and residues must survive every writable format.

use std::fs;
use std::path::PathBuf;

use tolviewer_core::{Alignment, Sequence};
use tolviewer_io::{parse, read_file, write_file, write_string, Format, LineEnding, WriteOptions};

/// Names are chosen to survive strict PHYLIP's 10-character field unchanged,
/// so the same expectations hold for every format.
fn fixture() -> Alignment {
    let mut aln = Alignment::new(
        "fixture",
        vec![
            Sequence::new("alpha", *b"ACGTACGTACGTACGTACGT"),
            Sequence::new("beta", *b"ACGTACGT----ACGTACGT"),
            Sequence::new("gamma", *b"acgtacgtacgtacgtacgt"),
            Sequence::new("delta_9999", *b"ACGT--GTACGTACGTACGT"),
        ],
    );
    aln.sequences[0].description = "a description".to_string();
    aln
}

fn writable() -> Vec<Format> {
    Format::all().iter().copied().filter(|f| f.can_write()).collect()
}

fn scratch(test: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("tolviewer-io-{test}"));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn write_string_round_trips_names_and_residues() {
    let aln = fixture();
    for format in writable() {
        for opts in [
            WriteOptions::default(),
            WriteOptions { interleaved: true, block_width: 7, ..Default::default() },
            WriteOptions { line_width: 0, line_ending: LineEnding::Crlf, ..Default::default() },
        ] {
            let text = write_string(&aln, format, &opts)
                .unwrap_or_else(|e| panic!("writing {}: {e}", format.name()));
            let back = parse(text.as_bytes(), format, "fixture")
                .unwrap_or_else(|e| panic!("re-reading {}:\n{text}\n{e}", format.name()));
            assert_eq!(back.len(), aln.len(), "{}", format.name());
            for (want, got) in aln.sequences.iter().zip(&back.sequences) {
                assert_eq!(got.id, want.id, "{} ids", format.name());
                assert_eq!(
                    got.residues,
                    want.residues,
                    "{} residues of {}",
                    format.name(),
                    want.id
                );
            }
        }
    }
}

#[test]
fn write_file_then_read_file_detects_the_format_it_wrote() {
    let aln = fixture();
    let dir = scratch("roundtrip-files");
    for format in writable() {
        let path = dir.join(format!("out.{}", format.extensions()[0]));
        write_file(&aln, &path, format, &WriteOptions::default()).unwrap();
        let back = read_file(&path).unwrap_or_else(|e| panic!("{}: {e}", format.name()));
        assert_eq!(back.len(), aln.len(), "{}", format.name());
        for (want, got) in aln.sequences.iter().zip(&back.sequences) {
            assert_eq!(got.id, want.id, "{}", format.name());
            assert_eq!(got.residues, want.residues, "{}", format.name());
        }
    }
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn descriptions_survive_only_in_fasta_and_fastq() {
    let aln = fixture();
    for format in writable() {
        let text = write_string(&aln, format, &WriteOptions::default()).unwrap();
        let back = parse(text.as_bytes(), format, "fixture").unwrap();
        let kept = back.sequences[0].description == "a description";
        assert_eq!(kept, matches!(format, Format::Fasta | Format::Fastq), "{}", format.name());
    }
}

#[test]
fn uppercase_option_applies_everywhere() {
    let aln = fixture();
    let opts = WriteOptions { uppercase: true, ..Default::default() };
    for format in writable() {
        let text = write_string(&aln, format, &opts).unwrap();
        let back = parse(text.as_bytes(), format, "fixture").unwrap();
        assert_eq!(back.sequences[2].residues, b"ACGTACGTACGTACGTACGT", "{}", format.name());
    }
}

#[test]
fn a_fixture_file_survives_conversion_to_every_writable_format() {
    let src = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../testdata/io/interleaved.phy");
    let aln = read_file(&src).unwrap();
    for format in writable() {
        let text = write_string(&aln, format, &WriteOptions::default()).unwrap();
        let back = parse(text.as_bytes(), format, "x").unwrap();
        let ids: Vec<&str> = back.sequences.iter().map(|s| s.id.as_str()).collect();
        assert_eq!(ids, vec!["Seq_A", "Seq_B", "Seq_C", "Seq_D"], "{}", format.name());
        assert_eq!(back.sequences[2].residues, aln.sequences[2].residues, "{}", format.name());
    }
}

#[test]
fn strict_phylip_truncates_long_names_but_keeps_them_unique() {
    let aln = Alignment::new(
        "t",
        vec![
            Sequence::new("Drosophila_melanogaster", *b"ACGT"),
            Sequence::new("Drosophila_simulans", *b"ACGT"),
            Sequence::new("Drosophila_yakuba", *b"ACGT"),
        ],
    );
    let text = write_string(&aln, Format::Phylip, &WriteOptions::default()).unwrap();
    let back = parse(text.as_bytes(), Format::Phylip, "t").unwrap();
    let ids: Vec<&str> = back.sequences.iter().map(|s| s.id.as_str()).collect();
    assert_eq!(ids, vec!["Drosophila", "Drosophil2", "Drosophil3"]);

    // Relaxed PHYLIP keeps them whole.
    let text = write_string(&aln, Format::PhylipRelaxed, &WriteOptions::default()).unwrap();
    let back = parse(text.as_bytes(), Format::PhylipRelaxed, "t").unwrap();
    assert_eq!(back.sequences[0].id, "Drosophila_melanogaster");
}

#[test]
fn names_with_spaces_are_sanitized_by_default_and_quoted_when_not() {
    let aln = Alignment::new(
        "t",
        vec![Sequence::new("Homo sapiens", *b"ACGT"), Sequence::new("Pan troglodytes", *b"ACGT")],
    );
    let text = write_string(&aln, Format::Nexus, &WriteOptions::default()).unwrap();
    let back = parse(text.as_bytes(), Format::Nexus, "t").unwrap();
    assert_eq!(back.sequences[0].id, "Homo_sapiens");

    let opts = WriteOptions { sanitize_names: false, ..Default::default() };
    let text = write_string(&aln, Format::Nexus, &opts).unwrap();
    let back = parse(text.as_bytes(), Format::Nexus, "t").unwrap();
    assert_eq!(back.sequences[0].id, "Homo sapiens");
}

#[test]
fn fastq_quality_round_trips() {
    let mut aln = Alignment::new("t", vec![Sequence::new("read1", *b"ACGT")]);
    aln.sequences[0].quality = Some(vec![40, 30, 20, 2]);
    let text = write_string(&aln, Format::Fastq, &WriteOptions::default()).unwrap();
    let back = parse(text.as_bytes(), Format::Fastq, "t").unwrap();
    assert_eq!(back.sequences[0].quality.as_ref().unwrap(), &vec![40, 30, 20, 2]);
}

#[test]
fn matrix_formats_refuse_ragged_input() {
    let ragged =
        Alignment::new("t", vec![Sequence::new("a", *b"ACGTACGT"), Sequence::new("b", *b"ACGT")]);
    for format in writable() {
        let result = write_string(&ragged, format, &WriteOptions::default());
        match format {
            Format::Fasta | Format::Fastq => assert!(result.is_ok(), "{}", format.name()),
            _ => assert!(
                matches!(result, Err(tolviewer_core::Error::Format(_))),
                "{} should refuse ragged input",
                format.name()
            ),
        }
    }
}

#[test]
fn hidden_rows_are_excluded_by_default() {
    let mut aln = fixture();
    aln.sequences[1].hidden = true;
    for format in writable() {
        let text = write_string(&aln, format, &WriteOptions::default()).unwrap();
        let back = parse(text.as_bytes(), format, "t").unwrap();
        assert_eq!(back.len(), 3, "{}", format.name());
        let opts = WriteOptions { include_hidden: true, ..Default::default() };
        let text = write_string(&aln, format, &opts).unwrap();
        let back = parse(text.as_bytes(), format, "t").unwrap();
        assert_eq!(back.len(), 4, "{}", format.name());
    }
}
