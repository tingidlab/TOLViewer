//! FASTA reader and writer.
//!
//! The reader is deliberately permissive: CRLF, blank lines anywhere, old-style
//! `;` comment lines, and whitespace or digits inside sequence lines are all
//! accepted, and case is preserved. An empty sequence body is allowed.

use tolviewer_core::{Alignment, Error, Result, Sequence};

use crate::options::WriteOptions;
use crate::util::{chunks, decode, lines, push_residues, residue_str, rows, sequence_from_header, Out};

const FORMAT: &str = "FASTA";

/// Parse FASTA bytes into an alignment named `name`.
pub(crate) fn parse(bytes: &[u8], name: &str) -> Result<Alignment> {
    let text = decode(bytes);
    let mut seqs: Vec<Sequence> = Vec::new();
    let mut header: Option<String> = None;
    let mut residues: Vec<u8> = Vec::new();

    for (n, line) in lines(&text) {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix('>') {
            if let Some(h) = header.take() {
                seqs.push(sequence_from_header(&h, std::mem::take(&mut residues)));
            }
            header = Some(rest.trim().to_string());
            continue;
        }
        if trimmed.starts_with(';') {
            continue; // old-style comment
        }
        if header.is_none() {
            return Err(Error::parse(
                FORMAT,
                Some(n),
                "sequence data before the first '>' header line",
            ));
        }
        push_residues(&mut residues, line);
    }
    if let Some(h) = header.take() {
        seqs.push(sequence_from_header(&h, residues));
    }
    if seqs.is_empty() {
        return Err(Error::parse(FORMAT, None, "no '>' header line found"));
    }
    Ok(Alignment::new(name, seqs))
}

/// Render an alignment as FASTA.
pub(crate) fn write(aln: &Alignment, opts: &WriteOptions) -> Result<String> {
    let mut out = Out::new(opts.line_ending);
    for row in rows(aln, opts) {
        out.line(format!(">{}", row.header()));
        if row.residues.is_empty() {
            continue;
        }
        for chunk in chunks(&row.residues, opts.line_width) {
            out.line(residue_str(chunk));
        }
    }
    Ok(out.finish())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::options::LineEnding;

    #[test]
    fn reads_basic_records() {
        let a = parse(b">a one\nACGT\n>b\nAC\nGT\n", "t").unwrap();
        assert_eq!(a.len(), 2);
        assert_eq!(a.sequences[0].id, "a");
        assert_eq!(a.sequences[0].description, "one");
        assert_eq!(a.sequences[0].residues, b"ACGT");
        assert_eq!(a.sequences[1].residues, b"ACGT");
    }

    #[test]
    fn tolerates_crlf_blank_lines_comments_and_junk() {
        let src = b"\r\n; a comment\r\n>a  spaced  desc \r\n AC GT 60\r\n\r\nac.gt\r\n";
        let a = parse(src, "t").unwrap();
        assert_eq!(a.sequences[0].id, "a");
        assert_eq!(a.sequences[0].description, "spaced  desc");
        // digits and whitespace stripped, case preserved, '.' normalised to '-'
        assert_eq!(a.sequences[0].residues, b"ACGTac-gt");
    }

    #[test]
    fn empty_body_is_allowed() {
        let a = parse(b">empty\n>b\nAC\n", "t").unwrap();
        assert_eq!(a.len(), 2);
        assert!(a.sequences[0].residues.is_empty());
    }

    #[test]
    fn data_before_header_is_an_error_with_a_line_number() {
        let e = parse(b"ACGT\n>a\nACGT\n", "t").unwrap_err();
        match e {
            Error::Parse { line, .. } => assert_eq!(line, Some(1)),
            other => panic!("wrong error: {other}"),
        }
    }

    #[test]
    fn writes_wrapped_and_unwrapped() {
        let aln = Alignment::new("t", vec![Sequence::new("a", *b"ACGTACGT")]);
        let mut o = WriteOptions {
            line_width: 4,
            ..Default::default()
        };
        assert_eq!(write(&aln, &o).unwrap(), ">a\nACGT\nACGT\n");
        o.line_width = 0;
        assert_eq!(write(&aln, &o).unwrap(), ">a\nACGTACGT\n");
        o.line_ending = LineEnding::Crlf;
        assert_eq!(write(&aln, &o).unwrap(), ">a\r\nACGTACGT\r\n");
    }

    #[test]
    fn hidden_rows_are_skipped_unless_asked_for() {
        let mut aln = Alignment::new(
            "t",
            vec![Sequence::new("a", *b"AC"), Sequence::new("b", *b"GT")],
        );
        aln.sequences[1].hidden = true;
        let mut o = WriteOptions::default();
        assert_eq!(write(&aln, &o).unwrap(), ">a\nAC\n");
        o.include_hidden = true;
        assert_eq!(write(&aln, &o).unwrap(), ">a\nAC\n>b\nGT\n");
    }
}
