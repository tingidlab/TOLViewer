//! PHYLIP reader and writer, strict (10-column names) and relaxed.
//!
//! Both the sequential and the interleaved layouts are read; the reader picks
//! one by looking at whether the first `ntax` data lines already carry `nchar`
//! residues, and falls back to the other layout if its first choice does not
//! validate against the header.

use std::collections::HashSet;

use tolviewer_core::{Alignment, Error, Result, Sequence};

use crate::format::phylip_header;
use crate::options::WriteOptions;
use crate::util::{decode, lines, push_residues, require_rectangular, residue_str, rows, Out};

const FORMAT: &str = "PHYLIP";

/// The strict name field width.
const NAME_WIDTH: usize = 10;

type Rows = Vec<(String, Vec<u8>)>;

/// Parse PHYLIP bytes. `strict` selects the 10-column name field.
pub(crate) fn parse(bytes: &[u8], name: &str, strict: bool) -> Result<Alignment> {
    let text = decode(bytes);
    let all: Vec<(usize, &str)> = lines(&text).collect();

    let header_pos = all
        .iter()
        .position(|(_, l)| !l.trim().is_empty())
        .ok_or_else(|| Error::parse(FORMAT, None, "file is empty"))?;
    let (hline, htext) = all[header_pos];
    let (ntax, nchar) = phylip_header(htext.trim()).ok_or_else(|| {
        Error::parse(FORMAT, Some(hline), "expected a header line of two numbers, 'ntax nchar'")
    })?;

    let data: Vec<(usize, &str)> =
        all[header_pos + 1..].iter().copied().filter(|(_, l)| !l.trim().is_empty()).collect();

    if ntax == 0 {
        return Ok(Alignment::new(name, Vec::new()));
    }

    // If the first ntax lines are already full-length rows the file is
    // sequential (or a single interleaved block, which parses identically).
    let first_block_complete =
        data.len() >= ntax && data[..ntax].iter().all(|(_, l)| residue_count(l, strict) == nchar);

    let rows = if first_block_complete {
        sequential(&data, ntax, nchar, strict)
    } else {
        interleaved(&data, ntax, nchar, strict)
            .or_else(|first| sequential(&data, ntax, nchar, strict).map_err(|_| first))
    }?;

    let seqs: Vec<Sequence> =
        rows.into_iter().map(|(id, residues)| Sequence::new(id, residues)).collect();
    Ok(Alignment::new(name, seqs))
}

/// Residues on a data line once its name field is removed.
fn residue_count(line: &str, strict: bool) -> usize {
    let (_, rest) = split_name(line, strict);
    let mut v = Vec::new();
    push_residues(&mut v, rest);
    v.len()
}

/// Split a data line into its name and the rest of the line.
fn split_name(line: &str, strict: bool) -> (String, &str) {
    if strict {
        let (name, rest) = split_at_chars(line, NAME_WIDTH);
        (name.trim().to_string(), rest)
    } else {
        match line.find(char::is_whitespace) {
            Some(i) => (line[..i].to_string(), &line[i..]),
            None => (line.to_string(), ""),
        }
    }
}

fn split_at_chars(line: &str, n: usize) -> (&str, &str) {
    match line.char_indices().nth(n) {
        Some((i, _)) => line.split_at(i),
        None => (line, ""),
    }
}

/// One sequence after another, each optionally wrapped over several lines.
fn sequential(data: &[(usize, &str)], ntax: usize, nchar: usize, strict: bool) -> Result<Rows> {
    let mut out: Rows = Vec::with_capacity(ntax);
    let mut i = 0usize;
    while out.len() < ntax {
        let (n, line) = *data.get(i).ok_or_else(|| {
            Error::parse(
                FORMAT,
                data.last().map(|(n, _)| *n),
                format!("header declares {ntax} sequences but the file only holds {}", out.len()),
            )
        })?;
        i += 1;
        let (id, rest) = split_name(line, strict);
        let mut residues = Vec::with_capacity(nchar);
        push_residues(&mut residues, rest);
        while residues.len() < nchar && i < data.len() {
            push_residues(&mut residues, data[i].1);
            i += 1;
        }
        if residues.len() != nchar {
            return Err(count_error(n, &id, residues.len(), nchar));
        }
        out.push((id, residues));
    }
    if i < data.len() {
        return Err(Error::parse(
            FORMAT,
            Some(data[i].0),
            format!("unexpected extra data after the {ntax} declared sequences"),
        ));
    }
    Ok(out)
}

/// Blocks of `ntax` lines; only the first block carries names.
fn interleaved(data: &[(usize, &str)], ntax: usize, nchar: usize, strict: bool) -> Result<Rows> {
    if data.len() < ntax {
        return Err(Error::parse(
            FORMAT,
            data.last().map(|(n, _)| *n),
            format!(
                "header declares {ntax} sequences but the first block only holds {}",
                data.len()
            ),
        ));
    }
    let mut names = Vec::with_capacity(ntax);
    let mut residues: Vec<Vec<u8>> = Vec::with_capacity(ntax);
    for (_, line) in &data[..ntax] {
        let (id, rest) = split_name(line, strict);
        let mut v = Vec::with_capacity(nchar);
        push_residues(&mut v, rest);
        names.push(id);
        residues.push(v);
    }

    let mut i = ntax;
    while i < data.len() {
        for r in 0..ntax {
            if i >= data.len() {
                break;
            }
            let (_, line) = data[i];
            i += 1;
            // Some writers repeat the names in every block; drop them.
            let line = strip_repeated_name(line, &names[r], strict);
            push_residues(&mut residues[r], line);
        }
    }

    for (r, v) in residues.iter().enumerate() {
        if v.len() != nchar {
            return Err(count_error(data[r].0, &names[r], v.len(), nchar));
        }
    }
    Ok(names.into_iter().zip(residues).collect())
}

fn strip_repeated_name<'a>(line: &'a str, name: &str, strict: bool) -> &'a str {
    if name.is_empty() {
        return line;
    }
    if strict {
        let (head, rest) = split_at_chars(line, NAME_WIDTH);
        if head.trim() == name && !rest.is_empty() {
            return rest;
        }
        line
    } else {
        match line.strip_prefix(name) {
            Some(rest) if rest.starts_with(char::is_whitespace) => rest,
            _ => line,
        }
    }
}

fn count_error(line: usize, id: &str, got: usize, want: usize) -> Error {
    Error::parse(
        FORMAT,
        Some(line),
        format!("sequence '{id}' has {got} residues but the header declares {want}"),
    )
}

/// Render an alignment as PHYLIP. `strict` writes 10-column names.
pub(crate) fn write(aln: &Alignment, opts: &WriteOptions, strict: bool) -> Result<String> {
    let strict = strict || opts.strict_phylip_names;
    let rows = rows(aln, opts);
    let nchar = require_rectangular(&rows, "PHYLIP")?;
    let ids: Vec<String> = rows.iter().map(|r| r.id.clone()).collect();
    let labels: Vec<String> = if strict {
        strict_names(&ids)?.iter().map(|s| format!("{s:<NAME_WIDTH$}")).collect()
    } else {
        let width = ids.iter().map(|s| s.chars().count()).max().unwrap_or(0);
        ids.iter().map(|s| format!("{s:<width$} ")).collect()
    };
    let blank: String = " ".repeat(labels.first().map_or(0, |l| l.chars().count()));

    let mut out = Out::new(opts.line_ending);
    out.line(format!("{} {}", rows.len(), nchar));
    let width = opts.effective_block_width();

    if opts.interleaved {
        let blocks = nchar.div_ceil(width).max(1);
        for b in 0..blocks {
            if b > 0 {
                out.blank();
            }
            for (row, label) in rows.iter().zip(&labels) {
                let start = b * width;
                let end = (start + width).min(nchar);
                let chunk = residue_str(&row.residues[start..end]).into_owned();
                if b == 0 {
                    out.line(format!("{label}{chunk}"));
                } else {
                    out.line(chunk);
                }
            }
        }
    } else {
        for (row, label) in rows.iter().zip(&labels) {
            let mut first = true;
            if row.residues.is_empty() {
                out.line(label.trim_end());
            }
            for chunk in row.residues.chunks(width) {
                let chunk = residue_str(chunk);
                if first {
                    out.line(format!("{label}{chunk}"));
                    first = false;
                } else {
                    out.line(format!("{blank}{chunk}"));
                }
            }
        }
    }
    Ok(out.finish())
}

/// Truncate names to 10 characters, keeping them unique by replacing the tail
/// with a counter. Fails only when no unique 10-character name exists.
fn strict_names(ids: &[String]) -> Result<Vec<String>> {
    let mut used: HashSet<String> = HashSet::new();
    let mut out = Vec::with_capacity(ids.len());
    for id in ids {
        let base: String = id.chars().take(NAME_WIDTH).collect();
        if used.insert(base.clone()) {
            out.push(base);
            continue;
        }
        let mut chosen = None;
        // A free candidate must appear within one try per existing name.
        for counter in 2..(ids.len() as u64 + 12) {
            let suffix = counter.to_string();
            if suffix.len() >= NAME_WIDTH {
                break;
            }
            let stem: String = id.chars().take(NAME_WIDTH - suffix.len()).collect();
            let cand = format!("{stem}{suffix}");
            if used.insert(cand.clone()) {
                chosen = Some(cand);
                break;
            }
        }
        match chosen {
            Some(c) => out.push(c),
            None => {
                return Err(Error::format(format!(
                    "cannot make the name '{id}' unique within the 10 characters \
                     strict PHYLIP allows; rename the sequences or use relaxed PHYLIP"
                )))
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const SEQUENTIAL: &[u8] = b"2 8\nSeq1      ACGTACGT\nSeq2      ACGTTTTT\n";
    const INTERLEAVED: &[u8] = b"2 8\nSeq1      ACGT\nSeq2      ACGT\n\nACGT\nTTTT\n";

    #[test]
    fn reads_sequential() {
        let a = parse(SEQUENTIAL, "t", true).unwrap();
        assert_eq!(a.len(), 2);
        assert_eq!(a.sequences[0].id, "Seq1");
        assert_eq!(a.sequences[1].residues, b"ACGTTTTT");
    }

    #[test]
    fn reads_interleaved() {
        let a = parse(INTERLEAVED, "t", true).unwrap();
        assert_eq!(a.sequences[0].residues, b"ACGTACGT");
        assert_eq!(a.sequences[1].residues, b"ACGTTTTT");
    }

    #[test]
    fn reads_wrapped_sequential() {
        let src = b"2 8\nSeq1      ACGT\nACGT\nSeq2      ACGT\nTTTT\n";
        let a = parse(src, "t", true).unwrap();
        assert_eq!(a.sequences[0].residues, b"ACGTACGT");
        assert_eq!(a.sequences[1].id, "Seq2");
        assert_eq!(a.sequences[1].residues, b"ACGTTTTT");
    }

    #[test]
    fn reads_interleaved_with_repeated_names() {
        let src = b"2 8\nSeq1      ACGT\nSeq2      ACGT\n\nSeq1      ACGT\nSeq2      TTTT\n";
        let a = parse(src, "t", true).unwrap();
        assert_eq!(a.sequences[0].residues, b"ACGTACGT");
        assert_eq!(a.sequences[1].residues, b"ACGTTTTT");
    }

    #[test]
    fn strict_names_may_contain_spaces() {
        let src = b"2 4\nHomo sap  ACGT\nPan trog  ACGT\n";
        let a = parse(src, "t", true).unwrap();
        assert_eq!(a.sequences[0].id, "Homo sap");
        assert_eq!(a.sequences[1].id, "Pan trog");
        assert_eq!(a.sequences[0].residues, b"ACGT");
    }

    #[test]
    fn relaxed_names_may_be_long() {
        let src = b"2 4\nAVeryLongTaxonName ACGT\nShort              ACGT\n";
        let a = parse(src, "t", false).unwrap();
        assert_eq!(a.sequences[0].id, "AVeryLongTaxonName");
        assert_eq!(a.sequences[1].id, "Short");
    }

    #[test]
    fn header_disagreement_is_an_error() {
        let e = parse(b"2 10\nSeq1      ACGT\nSeq2      ACGT\n", "t", true).unwrap_err();
        match e {
            Error::Parse { message, line, .. } => {
                assert!(message.contains("10"), "{message}");
                assert!(line.is_some());
            }
            other => panic!("wrong error: {other}"),
        }
        assert!(parse(b"3 4\nSeq1      ACGT\n", "t", true).is_err());
        assert!(parse(b"not a header\nSeq1 ACGT\n", "t", true).is_err());
    }

    #[test]
    fn digits_and_crlf_in_sequence_lines_are_ignored() {
        let src = b"2 8\r\nSeq1      ACGT ACGT 8\r\nSeq2      ACGT TTTT 8\r\n";
        let a = parse(src, "t", true).unwrap();
        assert_eq!(a.sequences[0].residues, b"ACGTACGT");
    }

    #[test]
    fn strict_name_truncation_stays_unique() {
        let ids: Vec<String> = vec![
            "Drosophila_melanogaster".into(),
            "Drosophila_simulans".into(),
            "Drosophila_yakuba".into(),
        ];
        let names = strict_names(&ids).unwrap();
        assert_eq!(names[0], "Drosophila");
        assert_eq!(names[1], "Drosophil2");
        assert_eq!(names[2], "Drosophil3");
        let unique: HashSet<&String> = names.iter().collect();
        assert_eq!(unique.len(), 3);
    }

    #[test]
    fn writes_and_re_reads_interleaved() {
        let aln = Alignment::new(
            "t",
            vec![Sequence::new("alpha", *b"ACGTACGTAC"), Sequence::new("beta", *b"ACGTTTTTAC")],
        );
        let opts = WriteOptions { interleaved: true, block_width: 4, ..Default::default() };
        let text = write(&aln, &opts, true).unwrap();
        assert!(text.starts_with("2 10\n"), "{text}");
        let back = parse(text.as_bytes(), "t", true).unwrap();
        assert_eq!(back.sequences[1].id, "beta");
        assert_eq!(back.sequences[1].residues, b"ACGTTTTTAC");
    }

    #[test]
    fn refuses_ragged_input() {
        let aln =
            Alignment::new("t", vec![Sequence::new("a", *b"ACGT"), Sequence::new("b", *b"AC")]);
        let e = write(&aln, &WriteOptions::default(), true).unwrap_err();
        assert!(matches!(e, Error::Format(_)));
    }
}
