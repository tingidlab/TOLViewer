//! FASTQ reader and writer.
//!
//! Both the classic four-line layout and the multi-line variant are read.
//! Quality is stored in [`tolviewer_core::Sequence::quality`] as Phred scores:
//! the file's offset is detected once for the whole file (Phred+33 unless some
//! byte is above `J` and none is below `;`, which can only be Phred+64).
//! Output is always Phred+33.

use tolviewer_core::{Alignment, Error, Result};

use crate::options::WriteOptions;
use crate::util::{decode, lines, residue_str, rows, sequence_from_header, Out};

const FORMAT: &str = "FASTQ";

/// Highest raw byte that is legal in Sanger (Phred+33) data: `J` = Q41.
const SANGER_MAX: u8 = b'J';
/// Lowest raw byte that Phred+64 data can contain: `;` = Q-5 in Solexa.
const SOLEXA_MIN: u8 = b';';

struct Record {
    header: String,
    residues: Vec<u8>,
    quality: Vec<u8>,
}

/// Parse FASTQ bytes into an alignment named `name`.
pub(crate) fn parse(bytes: &[u8], name: &str) -> Result<Alignment> {
    let text = decode(bytes);
    let all: Vec<(usize, &str)> = lines(&text).collect();
    let mut i = 0usize;
    let mut records: Vec<Record> = Vec::new();

    while i < all.len() {
        let (n, line) = all[i];
        if line.trim().is_empty() {
            i += 1;
            continue;
        }
        let header = match line.trim_end().strip_prefix('@') {
            Some(h) => h.trim().to_string(),
            None => {
                return Err(Error::parse(
                    FORMAT,
                    Some(n),
                    format!("expected a record starting with '@', found {:?}", truncate(line)),
                ))
            }
        };
        i += 1;

        // Sequence lines, up to the '+' separator.
        let mut residues = Vec::new();
        let mut saw_plus = false;
        while i < all.len() {
            let (_, l) = all[i];
            if l.starts_with('+') {
                saw_plus = true;
                i += 1;
                break;
            }
            residues.extend(l.trim().bytes().filter(|c| !c.is_ascii_whitespace()));
            i += 1;
        }
        if !saw_plus {
            return Err(Error::parse(
                FORMAT,
                Some(n),
                "truncated record: no '+' separator line before end of file",
            ));
        }

        // Quality lines, until we have as many characters as residues.
        let mut quality = Vec::new();
        while i < all.len() && quality.len() < residues.len() {
            let (_, l) = all[i];
            quality.extend(l.trim_end().bytes().filter(|c| !c.is_ascii_whitespace()));
            i += 1;
        }
        if quality.len() != residues.len() {
            return Err(Error::parse(
                FORMAT,
                Some(n),
                format!(
                    "record '{}' has {} residues but {} quality values",
                    header,
                    residues.len(),
                    quality.len()
                ),
            ));
        }
        records.push(Record { header, residues, quality });
    }

    if records.is_empty() {
        return Err(Error::parse(FORMAT, None, "no FASTQ records found"));
    }
    let offset = detect_offset(&records);
    let seqs = records
        .into_iter()
        .map(|r| {
            let mut s = sequence_from_header(&r.header, r.residues);
            s.quality = Some(r.quality.iter().map(|&q| q.saturating_sub(offset)).collect());
            s
        })
        .collect();
    Ok(Alignment::new(name, seqs))
}

/// 33 or 64: the ASCII offset the file's quality strings use.
fn detect_offset(records: &[Record]) -> u8 {
    let mut min = u8::MAX;
    let mut max = 0u8;
    for r in records {
        for &q in &r.quality {
            min = min.min(q);
            max = max.max(q);
        }
    }
    if max > SANGER_MAX && min >= SOLEXA_MIN {
        64
    } else {
        33
    }
}

fn truncate(s: &str) -> String {
    s.chars().take(30).collect()
}

/// Render an alignment as FASTQ (four lines per record, Phred+33).
///
/// Rows without quality are written as all-zero scores (`!`).
pub(crate) fn write(aln: &Alignment, opts: &WriteOptions) -> Result<String> {
    let mut out = Out::new(opts.line_ending);
    for row in rows(aln, opts) {
        out.line(format!("@{}", row.header()));
        out.line(residue_str(&row.residues));
        out.line("+");
        let mut q = Vec::with_capacity(row.residues.len());
        for i in 0..row.residues.len() {
            let score = row.quality.as_ref().and_then(|v| v.get(i).copied()).unwrap_or(0);
            q.push(33u8.saturating_add(score.min(93)));
        }
        out.line(residue_str(&q));
    }
    Ok(out.finish())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_four_line_records() {
        let a = parse(b"@r1 desc\nACGT\n+\nIIII\n@r2\nAC\n+r2\n!!\n", "t").unwrap();
        assert_eq!(a.len(), 2);
        assert_eq!(a.sequences[0].id, "r1");
        assert_eq!(a.sequences[0].description, "desc");
        assert_eq!(a.sequences[0].quality.as_ref().unwrap(), &vec![40; 4]);
        assert_eq!(a.sequences[1].quality.as_ref().unwrap(), &vec![0, 0]);
    }

    #[test]
    fn reads_multi_line_records() {
        let a = parse(b"@r1\nACGT\nACGT\n+\nIIII\nIIII\n", "t").unwrap();
        assert_eq!(a.sequences[0].residues, b"ACGTACGT");
        assert_eq!(a.sequences[0].quality.as_ref().unwrap().len(), 8);
    }

    #[test]
    fn detects_phred64() {
        // All bytes in 'a'..'h' (97..104): above 'J' and above ';'.
        let a = parse(b"@r1\nACGTACGT\n+\nabcdefgh\n", "t").unwrap();
        let q = a.sequences[0].quality.as_ref().unwrap();
        assert_eq!(q[0], 97 - 64);
        assert_eq!(q[7], 104 - 64);
    }

    #[test]
    fn keeps_phred33_when_a_low_byte_is_present() {
        // '!' (33) is below ';' so this cannot be Phred+64 even though 'h' is high.
        let a = parse(b"@r1\nAC\n+\n!h\n", "t").unwrap();
        let q = a.sequences[0].quality.as_ref().unwrap();
        assert_eq!(q[0], 0);
        assert_eq!(q[1], 104 - 33);
    }

    #[test]
    fn length_mismatch_is_an_error_with_a_line_number() {
        let e = parse(b"@r1\nACGT\n+\nII\n", "t").unwrap_err();
        match e {
            Error::Parse { line, message, .. } => {
                assert_eq!(line, Some(1));
                assert!(message.contains("4 residues"), "{message}");
            }
            other => panic!("wrong error: {other}"),
        }
    }

    #[test]
    fn missing_plus_line_is_an_error() {
        assert!(parse(b"@r1\nACGT\n", "t").is_err());
    }

    #[test]
    fn crlf_is_tolerated() {
        let a = parse(b"@r1\r\nACGT\r\n+\r\nIIII\r\n", "t").unwrap();
        assert_eq!(a.sequences[0].residues, b"ACGT");
        assert_eq!(a.sequences[0].quality.as_ref().unwrap().len(), 4);
    }

    #[test]
    fn writes_bang_for_rows_without_quality() {
        let aln = Alignment::new("t", vec![tolviewer_core::Sequence::new("a", *b"ACGT")]);
        assert_eq!(write(&aln, &WriteOptions::default()).unwrap(), "@a\nACGT\n+\n!!!!\n");
    }
}
