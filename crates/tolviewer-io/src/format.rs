//! The set of supported file formats, plus extension and content sniffing.

use std::path::Path;

use crate::util::{decode, lines};

/// A sequence or alignment file format.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Format {
    /// FASTA, `>id description` followed by residues.
    Fasta,
    /// FASTQ, four-line records with Phred quality.
    Fastq,
    /// Strict (10-character name field) PHYLIP.
    Phylip,
    /// Relaxed PHYLIP: names of any length, separated by whitespace.
    PhylipRelaxed,
    /// NEXUS `data`/`characters` blocks.
    Nexus,
    /// Clustal `.aln`.
    Clustal,
    /// Stockholm (Pfam/Rfam).
    Stockholm,
    /// GCG MSF / PileUp. Read only.
    Msf,
    /// GenBank flat file. Read only.
    Genbank,
}

const ALL: &[Format] = &[
    Format::Fasta,
    Format::Fastq,
    Format::Phylip,
    Format::PhylipRelaxed,
    Format::Nexus,
    Format::Clustal,
    Format::Stockholm,
    Format::Msf,
    Format::Genbank,
];

impl Format {
    /// Human name for menus, e.g. "FASTA".
    pub fn name(self) -> &'static str {
        match self {
            Format::Fasta => "FASTA",
            Format::Fastq => "FASTQ",
            Format::Phylip => "PHYLIP (strict)",
            Format::PhylipRelaxed => "PHYLIP (relaxed)",
            Format::Nexus => "NEXUS",
            Format::Clustal => "Clustal",
            Format::Stockholm => "Stockholm",
            Format::Msf => "MSF",
            Format::Genbank => "GenBank",
        }
    }

    /// Lowercase extensions without the dot, first is the default for saving.
    pub fn extensions(self) -> &'static [&'static str] {
        match self {
            Format::Fasta => &["fasta", "fa", "fas", "fna", "faa", "fst", "mpfa", "seq"],
            Format::Fastq => &["fastq", "fq"],
            Format::Phylip => &["phy", "phylip", "ph"],
            Format::PhylipRelaxed => &["phy", "phylip"],
            Format::Nexus => &["nex", "nexus", "nxs"],
            Format::Clustal => &["aln", "clustal", "clw"],
            Format::Stockholm => &["sto", "stk", "stockholm", "sth"],
            Format::Msf => &["msf"],
            Format::Genbank => &["gb", "gbk", "genbank", "gbff"],
        }
    }

    /// Every format can be read.
    pub fn can_read(self) -> bool {
        true
    }

    /// MSF and GenBank are read-only; everything else round-trips.
    pub fn can_write(self) -> bool {
        !matches!(self, Format::Msf | Format::Genbank)
    }

    /// All formats, for menu construction.
    pub fn all() -> &'static [Format] {
        ALL
    }

    /// Guess from a file extension.
    ///
    /// PHYLIP's extensions are shared between the strict and relaxed variants;
    /// this returns [`Format::Phylip`] for them, and reading falls back to the
    /// other variant if the file does not fit. Prefer [`Format::sniff`].
    pub fn from_path(path: &Path) -> Option<Format> {
        let ext = path.extension()?.to_str()?.to_ascii_lowercase();
        // `.gz` and friends are not handled, but a doubled extension like
        // `aln.txt` is common enough to be worth ignoring `txt`.
        if ext == "txt" {
            let stem = Path::new(path.file_stem()?);
            return Format::from_path(stem);
        }
        ALL.iter().copied().find(|f| f.extensions().contains(&ext.as_str()))
    }

    /// Guess from the first bytes of the file. Prefer this over the extension
    /// when they disagree and the content is unambiguous.
    pub fn sniff(bytes: &[u8]) -> Option<Format> {
        // Only the head matters, and callers may hand us a whole genome.
        let head = &bytes[..bytes.len().min(64 * 1024)];
        let text = decode(head);
        let mut data: Vec<&str> = Vec::new();
        for (_, line) in lines(&text) {
            if line.trim().is_empty() && data.is_empty() {
                continue; // leading blank lines
            }
            data.push(line);
            if data.len() >= 12 {
                break;
            }
        }
        let first = data.first()?.trim_start();
        let upper = first.to_ascii_uppercase();

        if upper.starts_with("#NEXUS") {
            return Some(Format::Nexus);
        }
        if upper.starts_with("# STOCKHOLM") || upper.starts_with("#STOCKHOLM") {
            return Some(Format::Stockholm);
        }
        if upper.contains("CLUSTAL") {
            return Some(Format::Clustal);
        }
        if upper.starts_with("LOCUS ") || upper.starts_with("LOCUS\t") {
            return Some(Format::Genbank);
        }
        if first.starts_with("!!") || upper.contains("PILEUP") || upper.contains("MSF:") {
            return Some(Format::Msf);
        }
        if first.starts_with('>') {
            return Some(Format::Fasta);
        }
        if first.starts_with(';') && data.iter().any(|l| l.starts_with('>')) {
            // Old-style FASTA with leading comments.
            return Some(Format::Fasta);
        }
        if first.starts_with('@') && data.iter().skip(1).take(6).any(|l| l.starts_with('+')) {
            return Some(Format::Fastq);
        }
        if phylip_header(first).is_some() {
            let second = data.get(1).copied().unwrap_or("");
            return Some(if strict_phylip_line(second) {
                Format::Phylip
            } else {
                Format::PhylipRelaxed
            });
        }
        None
    }
}

/// `ntax nchar` and nothing else.
pub(crate) fn phylip_header(line: &str) -> Option<(usize, usize)> {
    let mut it = line.split_whitespace();
    let ntax: usize = it.next()?.parse().ok()?;
    let nchar: usize = it.next()?.parse().ok()?;
    if it.next().is_some() {
        return None;
    }
    Some((ntax, nchar))
}

/// Does this data line look like strict PHYLIP, i.e. a 10-column name field?
///
/// Relaxed if the first token is longer than 10 characters (it would be cut in
/// half), or if the boundary at column 10 falls inside a run of non-blanks.
pub(crate) fn strict_phylip_line(line: &str) -> bool {
    let b = line.as_bytes();
    if b.len() <= 10 {
        return false;
    }
    let token = line.split_whitespace().next().unwrap_or("");
    if token.len() > 10 {
        return false;
    }
    !(!b[9].is_ascii_whitespace() && !b[10].is_ascii_whitespace())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extensions_map_back_to_formats() {
        for f in Format::all() {
            let ext = f.extensions()[0];
            let path = format!("x.{ext}");
            let got = Format::from_path(Path::new(&path)).unwrap();
            // PHYLIP's extensions are shared; both variants resolve to strict.
            if *f == Format::PhylipRelaxed {
                assert_eq!(got, Format::Phylip);
            } else {
                assert_eq!(got, *f);
            }
        }
    }

    #[test]
    fn double_extension_is_unwrapped() {
        assert_eq!(Format::from_path(Path::new("x.aln.txt")), Some(Format::Clustal));
    }

    #[test]
    fn sniffs_the_obvious_ones() {
        assert_eq!(Format::sniff(b">a\nACGT\n"), Some(Format::Fasta));
        assert_eq!(Format::sniff(b"@r1\nACGT\n+\nIIII\n"), Some(Format::Fastq));
        assert_eq!(Format::sniff(b"#NEXUS\nbegin data;\n"), Some(Format::Nexus));
        assert_eq!(
            Format::sniff(b"CLUSTAL W (1.81) multiple sequence alignment\n"),
            Some(Format::Clustal)
        );
        assert_eq!(Format::sniff(b"# STOCKHOLM 1.0\nseq ACGT\n//\n"), Some(Format::Stockholm));
        assert_eq!(Format::sniff(b"LOCUS       X 12 bp DNA linear\n"), Some(Format::Genbank));
        assert_eq!(Format::sniff(b"!!AA_MULTIPLE_ALIGNMENT 1.0\n"), Some(Format::Msf));
        assert_eq!(Format::sniff(b"not a sequence file\n"), None);
    }

    #[test]
    fn sniff_skips_leading_blank_lines_and_bom() {
        assert_eq!(Format::sniff(b"\n\n\n>a\nACGT\n"), Some(Format::Fasta));
        assert_eq!(Format::sniff(b"\xef\xbb\xbf#NEXUS\n"), Some(Format::Nexus));
        assert_eq!(Format::sniff(b"\r\n>a\r\nACGT\r\n"), Some(Format::Fasta));
    }

    #[test]
    fn sniff_separates_strict_and_relaxed_phylip() {
        assert_eq!(
            Format::sniff(b"2 8\nSeq1      ACGTACGT\nSeq2      ACGTACGT\n"),
            Some(Format::Phylip)
        );
        assert_eq!(
            Format::sniff(b"2 8\nSeq1 ACGTACGT\nSeq2 ACGTACGT\n"),
            Some(Format::PhylipRelaxed)
        );
        assert_eq!(
            Format::sniff(b"2 8\nAVeryLongTaxonName ACGTACGT\nB ACGTACGT\n"),
            Some(Format::PhylipRelaxed)
        );
    }

    #[test]
    fn phylip_header_rejects_extra_tokens() {
        assert_eq!(phylip_header(" 4 60 "), Some((4, 60)));
        assert_eq!(phylip_header("4 60 I"), None);
        assert_eq!(phylip_header("4"), None);
        assert_eq!(phylip_header("ACGT ACGT"), None);
    }

    #[test]
    fn read_write_capability_flags() {
        assert!(Format::Msf.can_read() && !Format::Msf.can_write());
        assert!(Format::Genbank.can_read() && !Format::Genbank.can_write());
        assert!(Format::Fasta.can_write());
    }
}
