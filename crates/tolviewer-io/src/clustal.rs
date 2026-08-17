//! Clustal (`.aln`) reader and writer.
//!
//! Blocks are `name  residues [count]`, optionally followed by a conservation
//! line made only of ` `, `.`, `:` and `*`, which must never be mistaken for a
//! sequence. Blocks are concatenated per name, in first-seen order.

use std::collections::HashMap;

use tolviewer_core::{Alignment, Error, Result, Sequence};

use crate::options::WriteOptions;
use crate::util::{
    decode, is_conservation_line, lines, push_residues, require_rectangular, residue_str, rows, Out,
};

const FORMAT: &str = "Clustal";

/// Parse Clustal bytes into an alignment named `name`.
pub(crate) fn parse(bytes: &[u8], name: &str) -> Result<Alignment> {
    let text = decode(bytes);
    let mut order: Vec<String> = Vec::new();
    let mut index: HashMap<String, usize> = HashMap::new();
    let mut residues: Vec<Vec<u8>> = Vec::new();
    let mut seen_header = false;

    for (n, line) in lines(&text) {
        if line.trim().is_empty() {
            continue;
        }
        if !seen_header {
            seen_header = true;
            // Any first line mentioning CLUSTAL (or the MUSCLE/Kalign variants
            // of the same header) is a header, not data.
            if line.to_ascii_uppercase().contains("CLUSTAL")
                || line.to_ascii_uppercase().contains("MULTIPLE SEQUENCE ALIGNMENT")
            {
                continue;
            }
        }
        if is_conservation_line(line) {
            continue;
        }
        let mut tokens = line.split_whitespace();
        let id = match tokens.next() {
            Some(t) => t.to_string(),
            None => continue,
        };
        let rest: Vec<&str> = tokens.collect();
        // A trailing running total is dropped (push_residues would too).
        let rest = match rest.split_last() {
            Some((last, head)) if last.bytes().all(|c| c.is_ascii_digit()) => head,
            _ => &rest[..],
        };
        if rest.is_empty() {
            return Err(Error::parse(
                FORMAT,
                Some(n),
                format!("line for '{id}' has a name but no residues"),
            ));
        }
        let row = *index.entry(id.clone()).or_insert_with(|| {
            order.push(id.clone());
            residues.push(Vec::new());
            residues.len() - 1
        });
        for part in rest {
            push_residues(&mut residues[row], part);
        }
    }

    if order.is_empty() {
        return Err(Error::parse(FORMAT, None, "no sequence blocks found"));
    }
    let seqs = order
        .into_iter()
        .zip(residues)
        .map(|(id, r)| Sequence::new(id, r))
        .collect();
    Ok(Alignment::new(name, seqs))
}

/// Render an alignment as Clustal.
pub(crate) fn write(aln: &Alignment, opts: &WriteOptions) -> Result<String> {
    let rows = rows(aln, opts);
    let ncols = require_rectangular(&rows, "Clustal")?;
    let width = opts.effective_block_width();
    let name_width = rows
        .iter()
        .map(|r| r.id.chars().count())
        .max()
        .unwrap_or(0)
        .max(10)
        + 6;

    let mut out = Out::new(opts.line_ending);
    out.line("CLUSTAL W (1.81) multiple sequence alignment");
    out.blank();
    out.blank();

    let blocks = ncols.div_ceil(width).max(1);
    for b in 0..blocks {
        let start = b * width;
        let end = (start + width).min(ncols);
        if b > 0 {
            out.blank();
        }
        for row in &rows {
            let chunk = residue_str(&row.residues[start..end]);
            out.line(format!("{:<name_width$}{}", row.id, chunk));
        }
        out.line(format!(
            "{:<name_width$}{}",
            "",
            conservation(&rows, start, end)
        ));
    }
    Ok(out.finish())
}

/// `*` under fully conserved, non-gap columns; a space elsewhere.
fn conservation(rows: &[crate::util::Row], start: usize, end: usize) -> String {
    (start..end)
        .map(|c| {
            let mut it = rows.iter().map(|r| r.residues[c].to_ascii_uppercase());
            match it.next() {
                Some(first) if !tolviewer_core::alphabet::is_gap(first) && it.all(|x| x == first) => {
                    '*'
                }
                _ => ' ',
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    const SRC: &[u8] = b"CLUSTAL W (1.81) multiple sequence alignment\n\n\
Seq1            ACGTACGTAC      10\n\
Seq2            ACGTTTGTAC      10\n\
                ****  ****\n\n\
Seq1            GGGG\n\
Seq2            GGGG\n\
                ****\n";

    #[test]
    fn reads_blocks_and_skips_conservation_lines() {
        let a = parse(SRC, "t").unwrap();
        assert_eq!(a.len(), 2);
        assert_eq!(a.sequences[0].id, "Seq1");
        assert_eq!(a.sequences[0].residues, b"ACGTACGTACGGGG");
        assert_eq!(a.sequences[1].residues, b"ACGTTTGTACGGGG");
    }

    #[test]
    fn accepts_other_headers_and_crlf() {
        let src = b"MUSCLE (3.8) multiple sequence alignment\r\n\r\nA  ACGT\r\nB  ACGT\r\n";
        let a = parse(src, "t").unwrap();
        assert_eq!(a.len(), 2);
        assert_eq!(a.sequences[0].residues, b"ACGT");
    }

    #[test]
    fn a_dot_only_line_is_conservation_not_sequence() {
        let src = b"CLUSTAL\n\nA  ACGT\nB  ACGT\n  ..::\n";
        let a = parse(src, "t").unwrap();
        assert_eq!(a.len(), 2);
    }

    #[test]
    fn round_trips() {
        let aln = Alignment::new(
            "t",
            vec![
                Sequence::new("alpha", *b"ACGTACGTAC"),
                Sequence::new("beta", *b"ACGTTTGTAC"),
            ],
        );
        let opts = WriteOptions {
            block_width: 4,
            ..Default::default()
        };
        let text = write(&aln, &opts).unwrap();
        let back = parse(text.as_bytes(), "t").unwrap();
        assert_eq!(back.sequences[0].id, "alpha");
        assert_eq!(back.sequences[0].residues, b"ACGTACGTAC");
        assert_eq!(back.sequences[1].residues, b"ACGTTTGTAC");
    }

    #[test]
    fn refuses_ragged_input() {
        let aln = Alignment::new(
            "t",
            vec![Sequence::new("a", *b"ACGT"), Sequence::new("b", *b"AC")],
        );
        assert!(matches!(
            write(&aln, &WriteOptions::default()),
            Err(Error::Format(_))
        ));
    }
}
