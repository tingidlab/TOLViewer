//! GenBank flat file reader. Writing GenBank is not supported.
//!
//! Only the parts that map onto [`tolviewer_core::Sequence`] are kept: the
//! LOCUS name becomes the id, DEFINITION the description, and the ORIGIN block
//! the residues. Features and references are ignored. Several records in one
//! file become several sequences.

use tolviewer_core::{Alignment, Error, Result, Sequence};

use crate::util::{decode, lines, push_residues};

const FORMAT: &str = "GenBank";

#[derive(Default)]
struct Record {
    id: String,
    definition: String,
    residues: Vec<u8>,
}

/// Finish the record in progress, if any.
fn flush(cur: &mut Option<Record>, seqs: &mut Vec<Sequence>) {
    if let Some(r) = cur.take() {
        let mut s = Sequence::new(r.id, r.residues);
        s.description = r.definition.trim_end_matches('.').trim().to_string();
        seqs.push(s);
    }
}

/// Parse GenBank bytes into an alignment named `name`.
pub(crate) fn parse(bytes: &[u8], name: &str) -> Result<Alignment> {
    let text = decode(bytes);
    let mut seqs: Vec<Sequence> = Vec::new();
    let mut cur: Option<Record> = None;
    let mut in_origin = false;
    let mut in_definition = false;

    for (_, line) in lines(&text) {
        let trimmed = line.trim_end();
        if let Some(rest) = trimmed.strip_prefix("LOCUS") {
            flush(&mut cur, &mut seqs);
            in_origin = false;
            in_definition = false;
            let id = rest.split_whitespace().next().unwrap_or("").to_string();
            cur = Some(Record { id, ..Default::default() });
            continue;
        }
        if trimmed.starts_with("//") {
            flush(&mut cur, &mut seqs);
            in_origin = false;
            in_definition = false;
            continue;
        }
        let record = match cur.as_mut() {
            Some(r) => r,
            None => continue, // preamble before the first LOCUS
        };
        if let Some(rest) = trimmed.strip_prefix("DEFINITION") {
            record.definition = rest.trim().to_string();
            in_definition = true;
            in_origin = false;
            continue;
        }
        if trimmed.starts_with("ORIGIN") {
            in_origin = true;
            in_definition = false;
            continue;
        }
        if in_definition {
            // Continuation lines of DEFINITION are indented.
            if trimmed.starts_with(' ') && !trimmed.trim().is_empty() {
                record.definition.push(' ');
                record.definition.push_str(trimmed.trim());
                continue;
            }
            in_definition = false;
        }
        if in_origin {
            push_residues(&mut record.residues, trimmed);
        } else if !trimmed.starts_with(' ') && !trimmed.trim().is_empty() {
            // A new top-level keyword ends any pending section.
            in_definition = false;
        }
    }
    flush(&mut cur, &mut seqs);

    if seqs.is_empty() {
        return Err(Error::parse(
            FORMAT,
            None,
            "no LOCUS line found: this does not look like a GenBank file",
        ));
    }
    Ok(Alignment::new(name, seqs))
}

#[cfg(test)]
mod tests {
    use super::*;

    // Written without escapes so the significant leading whitespace survives.
    const SRC: &str = "LOCUS       AB000001  12 bp  DNA  linear  INV 01-JAN-2000
DEFINITION  Example species mitochondrial COI gene,
            partial cds.
ACCESSION   AB000001
FEATURES             Location/Qualifiers
     source          1..12
ORIGIN
        1 acgtacg tacgt
//
LOCUS       AB000002  4 bp  DNA  linear  INV 01-JAN-2000
DEFINITION  Second record.
ORIGIN
        1 ggcc
//
";

    #[test]
    fn reads_locus_definition_and_origin() {
        let a = parse(SRC.as_bytes(), "t").unwrap();
        assert_eq!(a.len(), 2);
        assert_eq!(a.sequences[0].id, "AB000001");
        assert_eq!(
            a.sequences[0].description,
            "Example species mitochondrial COI gene, partial cds"
        );
        // Case is preserved and the coordinate numbers are dropped.
        assert_eq!(a.sequences[0].residues, b"acgtacgtacgt");
        assert_eq!(a.sequences[1].id, "AB000002");
        assert_eq!(a.sequences[1].residues, b"ggcc");
    }

    #[test]
    fn features_are_not_mistaken_for_sequence() {
        let a = parse(SRC.as_bytes(), "t").unwrap();
        assert_eq!(a.sequences[0].residues.len(), 12);
    }

    #[test]
    fn a_file_without_locus_is_an_error() {
        assert!(parse(b"ORIGIN\n 1 acgt\n//\n", "t").is_err());
    }
}
