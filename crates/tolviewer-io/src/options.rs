//! Options controlling how alignments are written.

/// Line ending to emit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LineEnding {
    /// Unix `\n`.
    Lf,
    /// DOS `\r\n`.
    Crlf,
}

/// Knobs for the writers. Readers need no options: they accept everything.
#[derive(Debug, Clone)]
pub struct WriteOptions {
    /// Residues per line for sequential formats (FASTA); 0 = one line per seq.
    pub line_width: usize,
    /// Write interleaved blocks (PHYLIP/NEXUS/Clustal).
    pub interleaved: bool,
    /// Residues per interleaved block.
    pub block_width: usize,
    /// Uppercase all residues on output.
    pub uppercase: bool,
    /// Truncate/pad names to 10 chars for strict PHYLIP.
    ///
    /// [`crate::Format::Phylip`] always does this; the flag additionally
    /// applies it to [`crate::Format::PhylipRelaxed`].
    pub strict_phylip_names: bool,
    /// Include rows flagged `hidden`.
    pub include_hidden: bool,
    /// Replace characters illegal in the target format (whitespace, quotes,
    /// parentheses, brackets, colons, semicolons, commas) in names with `_`.
    pub sanitize_names: bool,
    /// Line ending to emit.
    pub line_ending: LineEnding,
}

impl Default for WriteOptions {
    fn default() -> Self {
        WriteOptions {
            line_width: 60,
            interleaved: false,
            block_width: 60,
            uppercase: false,
            strict_phylip_names: false,
            include_hidden: false,
            sanitize_names: true,
            line_ending: LineEnding::Lf,
        }
    }
}

impl WriteOptions {
    /// The block width to actually use, substituting the default when the
    /// caller asked for 0 in a format that cannot put everything on one line.
    pub(crate) fn effective_block_width(&self) -> usize {
        if self.block_width == 0 {
            60
        } else {
            self.block_width
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_match_the_documented_contract() {
        let o = WriteOptions::default();
        assert_eq!(o.line_width, 60);
        assert!(!o.interleaved);
        assert_eq!(o.block_width, 60);
        assert!(!o.uppercase);
        assert!(!o.strict_phylip_names);
        assert!(!o.include_hidden);
        assert!(o.sanitize_names);
        assert_eq!(o.line_ending, LineEnding::Lf);
    }
}
