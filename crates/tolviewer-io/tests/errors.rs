//! Malformed input must produce a helpful error, never a panic.

use std::fs;

use tolviewer_core::Error;
use tolviewer_io::{parse, read_file, Format};

fn parse_error(bytes: &[u8], format: Format) -> (Option<usize>, String) {
    match parse(bytes, format, "t") {
        Err(Error::Parse { line, message, .. }) => (line, message),
        Err(other) => panic!("expected a parse error, got: {other}"),
        Ok(_) => panic!("expected a parse error, but the file was accepted"),
    }
}

#[test]
fn phylip_header_disagreeing_with_the_content() {
    let (line, message) =
        parse_error(b"3 20\nSeq_A     ACGTACGTAC\nSeq_B     ACGTACGTAC\n", Format::Phylip);
    assert!(line.is_some(), "the error should name a line");
    assert!(message.contains("3") || message.contains("20"), "unhelpful message: {message}");

    // Wrong nchar, right ntax.
    let (_, message) =
        parse_error(b"2 40\nSeq_A     ACGTACGTAC\nSeq_B     ACGTACGTAC\n", Format::Phylip);
    assert!(message.contains("40"), "unhelpful message: {message}");
}

#[test]
fn fastq_sequence_and_quality_lengths_must_match() {
    let (line, message) = parse_error(b"@r1\nACGT\n+\nII\n@r2\nAC\n+\n!!\n", Format::Fastq);
    assert_eq!(line, Some(1));
    assert!(message.contains("quality"), "unhelpful message: {message}");
}

#[test]
fn nexus_with_an_unterminated_matrix() {
    let src = b"#NEXUS\nbegin data;\n  dimensions ntax=2 nchar=4;\n  format datatype=dna;\n  matrix\n  a ACGT\n  b ACGT\n";
    let (line, message) = parse_error(src, Format::Nexus);
    assert_eq!(line, Some(5));
    assert!(message.contains("matrix"), "unhelpful message: {message}");
}

#[test]
fn nexus_with_an_unterminated_comment() {
    let (line, message) =
        parse_error(b"#NEXUS\n[ where does this end?\nbegin data;\n", Format::Nexus);
    assert_eq!(line, Some(2));
    assert!(message.contains("comment"), "unhelpful message: {message}");
}

#[test]
fn a_file_that_is_no_known_format() {
    let dir = std::env::temp_dir().join("tolviewer-io-errors");
    fs::create_dir_all(&dir).unwrap();
    let path = dir.join("notes.dat");
    fs::write(&path, b"Lab notebook, 14 March.\nNothing to see here.\n").unwrap();
    match read_file(&path) {
        Err(Error::Parse { message, .. }) => {
            assert!(message.contains("notes.dat"), "unhelpful message: {message}");
        }
        Err(other) => panic!("expected a parse error, got: {other}"),
        Ok(_) => panic!("a text file was accepted as an alignment"),
    }
    let _ = fs::remove_file(&path);
}

#[test]
fn empty_and_truncated_files_do_not_panic() {
    for format in Format::all() {
        for bytes in [b"".as_slice(), b"\n\n\n", b">", b"@", b"#NEXUS\n", b"2 4\n"] {
            // Any outcome is fine as long as it is not a panic.
            let _ = parse(bytes, *format, "t");
        }
    }
}

#[test]
fn a_fasta_file_read_as_phylip_fails_clearly() {
    let (_, message) = parse_error(b">a\nACGT\n>b\nACGT\n", Format::Phylip);
    assert!(message.contains("header line"), "unhelpful message: {message}");
}
