//! Small helpers shared by the format readers and writers.
//!
//! Nothing in here is public: the crate's contract is the handful of items
//! re-exported from `lib.rs`.

use std::borrow::Cow;

use tolviewer_core::{Alignment, Error, Result, Sequence};

use crate::options::{LineEnding, WriteOptions};

/// Decode input bytes as text, dropping a UTF-8 BOM and replacing invalid
/// sequences with U+FFFD. Sequence files are ASCII in practice, and a lossy
/// decode keeps a stray high byte from failing a whole file.
pub(crate) fn decode(bytes: &[u8]) -> Cow<'_, str> {
    let bytes = bytes.strip_prefix(&[0xEF, 0xBB, 0xBF]).unwrap_or(bytes);
    String::from_utf8_lossy(bytes)
}

/// Iterate `(1-based line number, line)` pairs, tolerating CRLF.
pub(crate) fn lines(text: &str) -> impl Iterator<Item = (usize, &str)> {
    text.split('\n').enumerate().map(|(i, l)| (i + 1, l.strip_suffix('\r').unwrap_or(l)))
}

/// Append the residues of `s` to `out`, discarding whitespace and digits
/// (sequence lines in the wild carry column rulers and running counts).
pub(crate) fn push_residues(out: &mut Vec<u8>, s: &str) {
    out.extend(s.bytes().filter(|c| !c.is_ascii_whitespace() && !c.is_ascii_digit()));
}

/// Residues as text. Residues are ASCII by construction.
pub(crate) fn residue_str(residues: &[u8]) -> Cow<'_, str> {
    String::from_utf8_lossy(residues)
}

/// True when every byte of `s` is one of the conservation symbols Clustal and
/// friends draw under a block (` `, `.`, `:`, `*`), and at least one is not a
/// space. Such lines must never be read as sequence data.
pub(crate) fn is_conservation_line(s: &str) -> bool {
    let mut any = false;
    for b in s.bytes() {
        match b {
            b' ' | b'\t' => {}
            b'.' | b':' | b'*' => any = true,
            _ => return false,
        }
    }
    any
}

/// Replace characters that are illegal or ambiguous in NEXUS/PHYLIP/Newick
/// names with `_`.
pub(crate) fn sanitize_name(name: &str) -> String {
    name.chars()
        .map(|c| {
            if c.is_whitespace()
                || matches!(c, '\'' | '"' | '(' | ')' | '[' | ']' | ':' | ';' | ',')
            {
                '_'
            } else {
                c
            }
        })
        .collect()
}

/// One row prepared for output: hidden rows dropped, names sanitized and
/// residues upper-cased according to `opts`.
pub(crate) struct Row {
    pub id: String,
    pub description: String,
    pub residues: Vec<u8>,
    pub quality: Option<Vec<u8>>,
}

impl Row {
    /// The FASTA-style header, `id` plus description when there is one.
    pub fn header(&self) -> String {
        if self.description.is_empty() {
            self.id.clone()
        } else {
            format!("{} {}", self.id, self.description)
        }
    }
}

/// Prepare the rows of `aln` for writing.
pub(crate) fn rows(aln: &Alignment, opts: &WriteOptions) -> Vec<Row> {
    aln.sequences
        .iter()
        .filter(|s| opts.include_hidden || !s.hidden)
        .map(|s| Row {
            id: if opts.sanitize_names { sanitize_name(&s.id) } else { s.id.clone() },
            description: s.description.clone(),
            residues: if opts.uppercase {
                s.residues.iter().map(|c| c.to_ascii_uppercase()).collect()
            } else {
                s.residues.clone()
            },
            quality: s.quality.clone(),
        })
        .collect()
}

/// The common column count of `rows`, or `Error::Format` when they are ragged.
///
/// PHYLIP, NEXUS, Clustal and Stockholm are all matrix formats: writing a
/// ragged set would silently produce a file nothing can read back.
pub(crate) fn require_rectangular(rows: &[Row], format: &str) -> Result<usize> {
    let mut widths = rows.iter().map(|r| r.residues.len());
    let first = widths.next().unwrap_or(0);
    if widths.all(|w| w == first) {
        Ok(first)
    } else {
        Err(Error::format(format!(
            "{format} needs an alignment: the sequences have different lengths \
             (pad or align them first)"
        )))
    }
}

/// Chunk residues for wrapped output. A width of 0 means "one line".
pub(crate) fn chunks(residues: &[u8], width: usize) -> std::slice::Chunks<'_, u8> {
    if width == 0 || residues.is_empty() {
        residues.chunks(residues.len().max(1))
    } else {
        residues.chunks(width)
    }
}

/// A text buffer that emits the line ending chosen in [`WriteOptions`].
pub(crate) struct Out {
    buf: String,
    nl: &'static str,
}

impl Out {
    pub fn new(ending: LineEnding) -> Self {
        Out {
            buf: String::new(),
            nl: match ending {
                LineEnding::Lf => "\n",
                LineEnding::Crlf => "\r\n",
            },
        }
    }

    /// Append `s` followed by a line ending.
    pub fn line(&mut self, s: impl AsRef<str>) {
        self.buf.push_str(s.as_ref());
        self.buf.push_str(self.nl);
    }

    /// Append an empty line.
    pub fn blank(&mut self) {
        self.buf.push_str(self.nl);
    }

    pub fn finish(self) -> String {
        self.buf
    }
}

/// Build a sequence, splitting `header` into id and description.
pub(crate) fn sequence_from_header(header: &str, residues: Vec<u8>) -> Sequence {
    let mut s = Sequence::new("", residues);
    s.set_header(header);
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lines_handle_crlf_and_number_from_one() {
        let got: Vec<_> = lines("a\r\nb\n").collect();
        assert_eq!(got[0], (1, "a"));
        assert_eq!(got[1], (2, "b"));
    }

    #[test]
    fn bom_is_stripped() {
        assert_eq!(decode(b"\xef\xbb\xbf>x"), ">x");
    }

    #[test]
    fn conservation_lines_detected() {
        assert!(is_conservation_line("  ***  ::..*"));
        assert!(!is_conservation_line("   "));
        assert!(!is_conservation_line("ACGT"));
        assert!(!is_conservation_line("----"));
    }

    #[test]
    fn residues_drop_whitespace_and_digits() {
        let mut v = Vec::new();
        push_residues(&mut v, "ACG T-A 60");
        assert_eq!(v, b"ACGT-A");
    }

    #[test]
    fn sanitize_replaces_illegal_characters() {
        assert_eq!(sanitize_name("Homo sapiens (COI):1"), "Homo_sapiens__COI__1");
    }
}
