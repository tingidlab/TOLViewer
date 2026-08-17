//! GCG MSF / PileUp reader. Writing MSF is not supported.
//!
//! The header runs up to the `//` separator and declares one
//! `Name: X Len: N Check: N Weight: N` line per sequence; the blocks that
//! follow repeat the names. `.` is the MSF gap character and is normalised to
//! [`tolviewer_core::GAP`] by `Sequence::new`.

use std::collections::HashMap;

use tolviewer_core::{Alignment, Error, Result, Sequence};

use crate::util::{decode, lines, push_residues};

const FORMAT: &str = "MSF";

/// Parse MSF bytes into an alignment named `name`.
pub(crate) fn parse(bytes: &[u8], name: &str) -> Result<Alignment> {
    let text = decode(bytes);
    let mut order: Vec<String> = Vec::new();
    let mut index: HashMap<String, usize> = HashMap::new();
    let mut residues: Vec<Vec<u8>> = Vec::new();
    let mut in_body = false;

    for (_, line) in lines(&text) {
        let trimmed = line.trim();
        if !in_body {
            if trimmed == "//" {
                in_body = true;
                continue;
            }
            // "Name: seq1  Len: 60  Check: 1234  Weight: 1.00"
            if let Some(rest) = find_name_field(trimmed) {
                add(&mut order, &mut index, &mut residues, rest);
            }
            continue;
        }
        if trimmed.is_empty() {
            continue;
        }
        let mut tokens = trimmed.split_whitespace();
        let id = match tokens.next() {
            Some(t) => t,
            None => continue,
        };
        // Coordinate ruler lines hold only numbers.
        if id.bytes().all(|c| c.is_ascii_digit()) {
            continue;
        }
        let row = add(&mut order, &mut index, &mut residues, id);
        for part in tokens {
            push_residues(&mut residues[row], part);
        }
    }

    if !in_body {
        return Err(Error::parse(
            FORMAT,
            None,
            "no '//' separator: this does not look like an MSF file",
        ));
    }
    if order.is_empty() {
        return Err(Error::parse(FORMAT, None, "no sequences found"));
    }
    let seqs = order.into_iter().zip(residues).map(|(id, r)| Sequence::new(id, r)).collect();
    Ok(Alignment::new(name, seqs))
}

/// Index of `id`, adding an empty row for it if this is the first sighting.
fn add(
    order: &mut Vec<String>,
    index: &mut HashMap<String, usize>,
    residues: &mut Vec<Vec<u8>>,
    id: &str,
) -> usize {
    *index.entry(id.to_string()).or_insert_with(|| {
        order.push(id.to_string());
        residues.push(Vec::new());
        residues.len() - 1
    })
}

/// The token after `Name:` on a header line, if any.
fn find_name_field(line: &str) -> Option<&str> {
    let pos = line.find("Name:").or_else(|| line.find("NAME:"))?;
    line[pos + 5..].split_whitespace().next()
}

#[cfg(test)]
mod tests {
    use super::*;

    const SRC: &[u8] = b"!!AA_MULTIPLE_ALIGNMENT 1.0\n\n\
 test.msf  MSF: 8  Type: P  Check: 0  ..\n\n\
 Name: alpha  Len: 8  Check: 1234  Weight: 1.00\n\
 Name: beta   Len: 8  Check: 5678  Weight: 1.00\n\n//\n\n\
           1                                        50\n\
alpha      MKVL ..LA\n\
beta       MKVL WWLA\n\n\
alpha      \n";

    #[test]
    fn reads_names_and_blocks() {
        let a = parse(SRC, "t").unwrap();
        assert_eq!(a.len(), 2);
        assert_eq!(a.sequences[0].id, "alpha");
        // '.' is the MSF gap character.
        assert_eq!(a.sequences[0].residues, b"MKVL--LA");
        assert_eq!(a.sequences[1].residues, b"MKVLWWLA");
    }

    #[test]
    fn header_order_wins_over_block_order() {
        let src =
            b"!!NA_MULTIPLE_ALIGNMENT 1.0\n Name: a Len: 4\n Name: b Len: 4\n//\nb ACGT\na TTTT\n";
        let a = parse(src, "t").unwrap();
        assert_eq!(a.sequences[0].id, "a");
        assert_eq!(a.sequences[0].residues, b"TTTT");
    }

    #[test]
    fn missing_separator_is_an_error() {
        assert!(parse(b"!!AA_MULTIPLE_ALIGNMENT\nName: a Len: 4\n", "t").is_err());
    }
}
