//! Stockholm (Pfam/Rfam) reader and writer.
//!
//! `#=GF`/`#=GS`/`#=GC`/`#=GR` annotation lines are skipped rather than
//! rejected; they are *not* preserved on write. Multi-block files are
//! concatenated per name; parsing stops at the `//` terminator.

use std::collections::HashMap;

use tolviewer_core::{Alignment, Error, Result, Sequence};

use crate::options::WriteOptions;
use crate::util::{decode, lines, push_residues, require_rectangular, residue_str, rows, Out};

const FORMAT: &str = "Stockholm";

/// Parse Stockholm bytes into an alignment named `name`.
pub(crate) fn parse(bytes: &[u8], name: &str) -> Result<Alignment> {
    let text = decode(bytes);
    let mut order: Vec<String> = Vec::new();
    let mut index: HashMap<String, usize> = HashMap::new();
    let mut residues: Vec<Vec<u8>> = Vec::new();

    for (n, line) in lines(&text) {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue; // header and #=GF/#=GS/#=GC/#=GR annotation
        }
        if trimmed == "//" {
            break;
        }
        let mut tokens = trimmed.split_whitespace();
        let id = match tokens.next() {
            Some(t) => t.to_string(),
            None => continue,
        };
        let rest: Vec<&str> = tokens.collect();
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
        return Err(Error::parse(FORMAT, None, "no sequence lines found"));
    }
    let seqs = order.into_iter().zip(residues).map(|(id, r)| Sequence::new(id, r)).collect();
    Ok(Alignment::new(name, seqs))
}

/// Render an alignment as Stockholm 1.0.
pub(crate) fn write(aln: &Alignment, opts: &WriteOptions) -> Result<String> {
    let rows = rows(aln, opts);
    let ncols = require_rectangular(&rows, "Stockholm")?;
    let name_width = rows.iter().map(|r| r.id.chars().count()).max().unwrap_or(0) + 2;

    let mut out = Out::new(opts.line_ending);
    out.line("# STOCKHOLM 1.0");
    out.blank();
    if opts.interleaved {
        let width = opts.effective_block_width();
        let blocks = ncols.div_ceil(width).max(1);
        for b in 0..blocks {
            if b > 0 {
                out.blank();
            }
            for row in &rows {
                let start = b * width;
                let end = (start + width).min(ncols);
                let chunk = residue_str(&row.residues[start..end]);
                out.line(format!("{:<name_width$}{}", row.id, chunk));
            }
        }
    } else {
        for row in &rows {
            out.line(format!("{:<name_width$}{}", row.id, residue_str(&row.residues)));
        }
    }
    out.line("//");
    Ok(out.finish())
}

#[cfg(test)]
mod tests {
    use super::*;

    const SRC: &[u8] = b"# STOCKHOLM 1.0\n#=GF ID test\n\n\
seq1 ACGT\n\
seq2 ACGA\n\
#=GC SS_cons ....\n\n\
seq1 TTTT\n\
seq2 TTTA\n\
//\n";

    #[test]
    fn reads_blocks_and_skips_annotation() {
        let a = parse(SRC, "t").unwrap();
        assert_eq!(a.len(), 2);
        assert_eq!(a.sequences[0].residues, b"ACGTTTTT");
        assert_eq!(a.sequences[1].residues, b"ACGATTTA");
    }

    #[test]
    fn stops_at_the_terminator() {
        let a = parse(b"# STOCKHOLM 1.0\nx ACGT\n//\ny GGGG\n", "t").unwrap();
        assert_eq!(a.len(), 1);
    }

    #[test]
    fn round_trips() {
        let aln = Alignment::new(
            "t",
            vec![Sequence::new("alpha", *b"ACGTACGT"), Sequence::new("beta", *b"ACGT--GT")],
        );
        let text = write(&aln, &WriteOptions::default()).unwrap();
        assert!(text.starts_with("# STOCKHOLM 1.0\n"));
        assert!(text.trim_end().ends_with("//"));
        let back = parse(text.as_bytes(), "t").unwrap();
        assert_eq!(back.sequences[1].id, "beta");
        assert_eq!(back.sequences[1].residues, b"ACGT--GT");
    }

    #[test]
    fn interleaved_round_trips() {
        let aln = Alignment::new("t", vec![Sequence::new("a", *b"ACGTACGT")]);
        let opts = WriteOptions { interleaved: true, block_width: 3, ..Default::default() };
        let text = write(&aln, &opts).unwrap();
        let back = parse(text.as_bytes(), "t").unwrap();
        assert_eq!(back.sequences[0].residues, b"ACGTACGT");
    }
}
