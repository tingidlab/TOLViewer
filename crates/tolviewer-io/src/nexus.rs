//! NEXUS reader and writer.
//!
//! The reader understands `data` and `characters` blocks: `dimensions`,
//! `format` (datatype, interleave, gap, missing, matchchar) and `matrix`, in
//! both the interleaved and the sequential layouts. `[square bracket comments]`
//! nest, `'single quoted names'` may contain spaces (with `''` for a literal
//! quote), and blocks that are not understood (trees, assumptions, sets,
//! mrbayes, ...) are skipped rather than rejected.
//!
//! The writer emits a minimal `#NEXUS` + `begin data;` file.

use std::collections::HashMap;

use tolviewer_core::{Alignment, Alphabet, Error, Result, Sequence, GAP, MISSING};

use crate::options::WriteOptions;
use crate::util::{
    decode, push_residues, require_rectangular, residue_str, rows as prepare_rows, Out,
};

const FORMAT: &str = "NEXUS";

// ---------------------------------------------------------------- tokenizer

#[derive(Debug, Clone)]
struct Token {
    text: String,
    line: usize,
    quoted: bool,
}

impl Token {
    fn is(&self, word: &str) -> bool {
        !self.quoted && self.text.eq_ignore_ascii_case(word)
    }
}

/// Split NEXUS text into words, `;` and `=`, dropping nested comments and
/// unwrapping single-quoted names.
fn tokenize(text: &str) -> Result<Vec<Token>> {
    let chars: Vec<char> = text.chars().collect();
    let mut toks: Vec<Token> = Vec::new();
    let mut cur = String::new();
    let mut cur_line = 1usize;
    let mut line = 1usize;
    let mut i = 0usize;

    while i < chars.len() {
        let c = chars[i];
        match c {
            '[' => {
                flush(&mut toks, &mut cur, cur_line);
                let start = line;
                let mut depth = 1usize;
                i += 1;
                while i < chars.len() && depth > 0 {
                    match chars[i] {
                        '[' => depth += 1,
                        ']' => depth -= 1,
                        '\n' => line += 1,
                        _ => {}
                    }
                    i += 1;
                }
                if depth > 0 {
                    return Err(Error::parse(FORMAT, Some(start), "unterminated [comment"));
                }
                continue;
            }
            '\'' => {
                flush(&mut toks, &mut cur, cur_line);
                let start = line;
                let mut value = String::new();
                i += 1;
                loop {
                    match chars.get(i) {
                        None => {
                            return Err(Error::parse(
                                FORMAT,
                                Some(start),
                                "unterminated quoted name",
                            ))
                        }
                        Some('\'') => {
                            if chars.get(i + 1) == Some(&'\'') {
                                value.push('\'');
                                i += 2;
                            } else {
                                i += 1;
                                break;
                            }
                        }
                        Some(&ch) => {
                            if ch == '\n' {
                                line += 1;
                            }
                            value.push(ch);
                            i += 1;
                        }
                    }
                }
                toks.push(Token { text: value, line: start, quoted: true });
                continue;
            }
            ';' | '=' => {
                flush(&mut toks, &mut cur, cur_line);
                toks.push(Token { text: c.to_string(), line, quoted: false });
            }
            '\n' => {
                flush(&mut toks, &mut cur, cur_line);
                line += 1;
            }
            c if c.is_whitespace() => flush(&mut toks, &mut cur, cur_line),
            c => {
                if cur.is_empty() {
                    cur_line = line;
                }
                cur.push(c);
            }
        }
        i += 1;
    }
    flush(&mut toks, &mut cur, cur_line);
    Ok(toks)
}

fn flush(toks: &mut Vec<Token>, cur: &mut String, line: usize) {
    if !cur.is_empty() {
        toks.push(Token { text: std::mem::take(cur), line, quoted: false });
    }
}

// ------------------------------------------------------------------- reader

/// State gathered from one `data`/`characters` block.
#[derive(Default)]
struct Block {
    ntax: Option<usize>,
    nchar: Option<usize>,
    datatype: Option<Alphabet>,
    interleave: Option<bool>,
    gap: Option<u8>,
    missing: Option<u8>,
    matchchar: Option<u8>,
}

/// Parse NEXUS bytes into an alignment named `name`.
pub(crate) fn parse(bytes: &[u8], name: &str) -> Result<Alignment> {
    let text = decode(bytes);
    let toks = tokenize(&text)?;
    if !toks.first().is_some_and(|t| t.is("#nexus")) {
        return Err(Error::parse(FORMAT, Some(1), "missing the '#NEXUS' first line"));
    }

    let mut i = 1usize;
    while i < toks.len() {
        if !toks[i].is("begin") {
            i += 1;
            continue;
        }
        let block_name = toks.get(i + 1).map(|t| t.text.to_ascii_lowercase()).unwrap_or_default();
        let mut body_start = (i + 2).min(toks.len());
        if toks.get(body_start).is_some_and(|t| t.text == ";") {
            body_start += 1;
        }
        let end = find_block_end(&toks, body_start);
        let body = &toks[body_start.min(end)..end];
        if block_name == "data" || block_name.starts_with("characters") {
            return parse_data_block(body, name);
        }
        // Not a block we understand (taxa, trees, assumptions, sets, mrbayes,
        // codons, ...): skip it wholesale.
        i = (end + 2).max(i + 1);
    }
    Err(Error::parse(FORMAT, None, "no 'begin data;' or 'begin characters;' block found"))
}

/// Index of the `end`/`endblock` token closing the block that starts at `from`.
fn find_block_end(toks: &[Token], from: usize) -> usize {
    let mut i = from;
    while i < toks.len() {
        if toks[i].is("end") || toks[i].is("endblock") {
            return i;
        }
        i += 1;
    }
    toks.len()
}

fn parse_data_block(body: &[Token], name: &str) -> Result<Alignment> {
    let mut block = Block::default();
    let mut matrix: Option<&[Token]> = None;

    // Commands are separated by ';'.
    let mut start = 0usize;
    let mut i = 0usize;
    while i <= body.len() {
        let at_end = i == body.len();
        if !at_end && !(body[i].text == ";" && !body[i].quoted) {
            i += 1;
            continue;
        }
        let cmd = &body[start..i];
        if let Some(head) = cmd.first() {
            let keyword = head.text.to_ascii_lowercase();
            if at_end && keyword == "matrix" {
                // The command ran to the end of the block without its ';'.
                return Err(Error::parse(
                    FORMAT,
                    Some(head.line),
                    "unterminated matrix: no ';' closing the matrix command",
                ));
            }
            match keyword.as_str() {
                "dimensions" => read_dimensions(cmd, &mut block),
                "format" => read_format(cmd, &mut block),
                "matrix" => matrix = Some(&cmd[1..]),
                _ => {}
            }
        }
        i += 1;
        start = i;
    }

    let matrix = matrix
        .ok_or_else(|| Error::parse(FORMAT, None, "the data block has no 'matrix' command"))?;
    let rows = read_matrix(matrix, &block)?;

    if let Some(ntax) = block.ntax {
        if rows.len() != ntax {
            return Err(Error::parse(
                FORMAT,
                matrix.first().map(|t| t.line),
                format!("dimensions declares ntax={ntax} but the matrix holds {} rows", rows.len()),
            ));
        }
    }

    let mut seqs = Vec::with_capacity(rows.len());
    for row in &rows {
        if let Some(nchar) = block.nchar {
            if row.residues.len() != nchar {
                return Err(Error::parse(
                    FORMAT,
                    Some(row.line),
                    format!(
                        "dimensions declares nchar={nchar} but '{}' has {} characters",
                        row.id,
                        row.residues.len()
                    ),
                ));
            }
        }
        seqs.push(Sequence::new(row.id.clone(), row.residues.clone()));
    }
    let mut aln = Alignment::new(name, seqs);
    if let Some(a) = block.datatype {
        aln.set_alphabet(a);
    }
    Ok(aln)
}

fn read_dimensions(cmd: &[Token], block: &mut Block) {
    for (key, value) in pairs(cmd) {
        match key.as_str() {
            "ntax" => block.ntax = value.and_then(|v| v.parse().ok()),
            "nchar" => block.nchar = value.and_then(|v| v.parse().ok()),
            _ => {}
        }
    }
}

fn read_format(cmd: &[Token], block: &mut Block) {
    for (key, value) in pairs(cmd) {
        match key.as_str() {
            "datatype" => {
                block.datatype = match value.unwrap_or_default().to_ascii_lowercase().as_str() {
                    "dna" | "nucleotide" => Some(Alphabet::Dna),
                    "rna" => Some(Alphabet::Rna),
                    "protein" | "aminoacid" => Some(Alphabet::Protein),
                    _ => None,
                }
            }
            "interleave" => {
                block.interleave = Some(!matches!(
                    value.unwrap_or("yes").to_ascii_lowercase().as_str(),
                    "no" | "false"
                ))
            }
            "gap" => block.gap = first_byte(value),
            "missing" => block.missing = first_byte(value),
            "matchchar" => block.matchchar = first_byte(value),
            _ => {}
        }
    }
}

fn first_byte(value: Option<&str>) -> Option<u8> {
    value.and_then(|v| v.bytes().next())
}

/// `key`, `key=value` and `key = value` pairs following a command keyword.
fn pairs(cmd: &[Token]) -> Vec<(String, Option<&str>)> {
    let mut out = Vec::new();
    let mut i = 1usize; // skip the command keyword
    while i < cmd.len() {
        let key = cmd[i].text.to_ascii_lowercase();
        if cmd.get(i + 1).is_some_and(|t| t.text == "=" && !t.quoted) {
            let value = cmd.get(i + 2).map(|t| t.text.as_str());
            out.push((key, value));
            i += 3;
        } else {
            out.push((key, None));
            i += 1;
        }
    }
    out
}

struct Row {
    id: String,
    residues: Vec<u8>,
    line: usize,
}

/// Turn the matrix tokens into rows, honouring interleaving, `matchchar` and
/// the declared gap/missing characters.
fn read_matrix(matrix: &[Token], block: &Block) -> Result<Vec<Row>> {
    // Group tokens by source line: one line is one matrix row (or one row of
    // one interleaved block).
    let mut groups: Vec<Vec<&Token>> = Vec::new();
    for t in matrix {
        match groups.last_mut() {
            Some(g) if g[0].line == t.line => g.push(t),
            _ => groups.push(vec![t]),
        }
    }

    // Interleaved unless told otherwise; when the format command is silent,
    // a repeated leading name gives it away.
    let interleaved = block.interleave.unwrap_or_else(|| {
        let mut seen: HashMap<&str, ()> = HashMap::new();
        groups.iter().any(|g| seen.insert(g[0].text.as_str(), ()).is_some())
    });

    let mut order: Vec<String> = Vec::new();
    let mut index: HashMap<String, usize> = HashMap::new();
    let mut rows: Vec<Row> = Vec::new();
    let mut cur: Option<usize> = None;

    for g in &groups {
        let continuing = !interleaved
            && cur.is_some_and(|c| block.nchar.is_some_and(|n| rows[c].residues.len() < n));
        if continuing {
            let c = cur.unwrap_or(0);
            for t in g {
                push_residues(&mut rows[c].residues, &t.text);
            }
            continue;
        }
        let id = g[0].text.clone();
        let row = match index.get(&id) {
            Some(&r) => r,
            None => {
                order.push(id.clone());
                rows.push(Row { id: id.clone(), residues: Vec::new(), line: g[0].line });
                index.insert(id, rows.len() - 1);
                rows.len() - 1
            }
        };
        for t in &g[1..] {
            push_residues(&mut rows[row].residues, &t.text);
        }
        cur = Some(row);
    }

    // `matchchar` copies from the first row; do this before gap normalisation
    // because '.' is a popular matchchar.
    if let Some(mc) = block.matchchar {
        if !rows.is_empty() {
            let first = rows[0].residues.clone();
            for row in rows.iter_mut().skip(1) {
                for (c, r) in row.residues.iter_mut().enumerate() {
                    if *r == mc {
                        *r = first.get(c).copied().unwrap_or(GAP);
                    }
                }
            }
        }
    }
    for row in &mut rows {
        for r in &mut row.residues {
            if Some(*r) == block.gap {
                *r = GAP;
            } else if Some(*r) == block.missing {
                *r = MISSING;
            }
        }
    }
    Ok(rows)
}

// ------------------------------------------------------------------- writer

/// Render an alignment as a minimal NEXUS `data` block.
pub(crate) fn write(aln: &Alignment, opts: &WriteOptions) -> Result<String> {
    let rows = prepare_rows(aln, opts);
    let nchar = require_rectangular(&rows, "NEXUS")?;
    let alphabet = aln
        .alphabet_hint()
        .unwrap_or_else(|| Alphabet::guess(rows.iter().flat_map(|r| r.residues.iter().copied())));
    let names: Vec<String> = rows.iter().map(|r| quote_name(&r.id)).collect();
    let name_width = names.iter().map(|n| n.chars().count()).max().unwrap_or(0) + 2;

    let mut out = Out::new(opts.line_ending);
    out.line("#NEXUS");
    out.blank();
    out.line("begin data;");
    out.line(format!("    dimensions ntax={} nchar={};", rows.len(), nchar));
    out.line(format!(
        "    format datatype={} missing={} gap={}{};",
        alphabet.nexus_datatype(),
        MISSING as char,
        GAP as char,
        if opts.interleaved { " interleave=yes" } else { "" }
    ));
    out.line("    matrix");

    let width = opts.effective_block_width();
    let blocks = if opts.interleaved { nchar.div_ceil(width).max(1) } else { 1 };
    for b in 0..blocks {
        if b > 0 {
            out.blank();
        }
        for (row, name) in rows.iter().zip(&names) {
            let (start, end) = if opts.interleaved {
                (b * width, ((b + 1) * width).min(nchar))
            } else {
                (0, nchar)
            };
            out.line(format!("{:<name_width$}{}", name, residue_str(&row.residues[start..end])));
        }
    }
    out.line("    ;");
    out.line("end;");
    Ok(out.finish())
}

/// Quote a name if NEXUS would not read it back as a single token.
fn quote_name(name: &str) -> String {
    let needs = name.is_empty()
        || name.chars().any(|c| {
            c.is_whitespace()
                || matches!(
                    c,
                    '(' | ')'
                        | '['
                        | ']'
                        | '{'
                        | '}'
                        | '/'
                        | '\\'
                        | ','
                        | ';'
                        | ':'
                        | '='
                        | '*'
                        | '\''
                        | '"'
                        | '`'
                        | '<'
                        | '>'
                        | '-'
                )
        });
    if needs {
        format!("'{}'", name.replace('\'', "''"))
    } else {
        name.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SEQUENTIAL: &[u8] = b"#NEXUS\n\
[ a comment [ nested ] still a comment ]\n\
begin taxa;\n  dimensions ntax=2;\n  taxlabels a b;\nend;\n\
begin data;\n\
  dimensions ntax=2 nchar=8;\n\
  format datatype=dna missing=? gap=-;\n\
  matrix\n\
    alpha ACGTACGT\n\
    'Homo sapiens' ACGT--GT\n\
  ;\n\
end;\n\
begin trees;\n  tree x = ((a,b),c);\nend;\n";

    const INTERLEAVED: &[u8] = b"#NEXUS\nbegin data;\n\
dimensions ntax=2 nchar=8;\n\
format datatype=protein interleave=yes gap=- missing=?;\n\
matrix\n\
alpha MKVL\n\
beta  MKVL\n\n\
alpha WIPQ\n\
beta  WIPA\n\
;\nend;\n";

    #[test]
    fn reads_sequential_with_comments_and_quoted_names() {
        let mut a = parse(SEQUENTIAL, "t").unwrap();
        assert_eq!(a.len(), 2);
        assert_eq!(a.sequences[0].id, "alpha");
        assert_eq!(a.sequences[1].id, "Homo sapiens");
        assert_eq!(a.sequences[1].residues, b"ACGT--GT");
        assert_eq!(a.alphabet(), Alphabet::Dna);
    }

    #[test]
    fn reads_interleaved() {
        let mut a = parse(INTERLEAVED, "t").unwrap();
        assert_eq!(a.sequences[0].residues, b"MKVLWIPQ");
        assert_eq!(a.sequences[1].residues, b"MKVLWIPA");
        assert_eq!(a.alphabet(), Alphabet::Protein);
    }

    #[test]
    fn detects_interleaving_without_the_flag() {
        let src = b"#NEXUS\nbegin data;\ndimensions ntax=2 nchar=8;\nformat datatype=dna;\nmatrix\na ACGT\nb ACGT\na TTTT\nb GGGG\n;\nend;\n";
        let a = parse(src, "t").unwrap();
        assert_eq!(a.sequences[0].residues, b"ACGTTTTT");
        assert_eq!(a.sequences[1].residues, b"ACGTGGGG");
    }

    #[test]
    fn reads_wrapped_sequential_rows() {
        let src = b"#NEXUS\nbegin data;\ndimensions ntax=2 nchar=8;\nformat datatype=dna;\nmatrix\na ACGT\nTTTT\nb ACGT\nGGGG\n;\nend;\n";
        let a = parse(src, "t").unwrap();
        assert_eq!(a.sequences[0].residues, b"ACGTTTTT");
        assert_eq!(a.sequences[1].residues, b"ACGTGGGG");
    }

    #[test]
    fn expands_matchchar() {
        let src = b"#NEXUS\nbegin data;\ndimensions ntax=2 nchar=4;\nformat datatype=dna matchchar=.;\nmatrix\na ACGT\nb .C.T\n;\nend;\n";
        let a = parse(src, "t").unwrap();
        assert_eq!(a.sequences[1].residues, b"ACGT");
    }

    #[test]
    fn normalises_declared_gap_and_missing() {
        let src = b"#NEXUS\nbegin data;\ndimensions ntax=1 nchar=4;\nformat datatype=dna gap=_ missing=N;\nmatrix\na A_CN\n;\nend;\n";
        let a = parse(src, "t").unwrap();
        assert_eq!(a.sequences[0].residues, b"A-C?");
    }

    #[test]
    fn ignores_blocks_it_does_not_understand() {
        let src = b"#NEXUS\nbegin mrbayes;\n lset nst=6;\n mcmc ngen=100;\nend;\nbegin sets;\n charset gene1 = 1-4;\nend;\nbegin data;\ndimensions ntax=1 nchar=4;\nformat datatype=dna;\nmatrix\na ACGT\n;\nend;\n";
        let a = parse(src, "t").unwrap();
        assert_eq!(a.len(), 1);
        assert_eq!(a.sequences[0].residues, b"ACGT");
    }

    #[test]
    fn unterminated_matrix_is_an_error() {
        let src = b"#NEXUS\nbegin data;\ndimensions ntax=2 nchar=4;\nformat datatype=dna;\nmatrix\na ACGT\nb ACGT\n";
        let e = parse(src, "t").unwrap_err();
        match e {
            Error::Parse { message, line, .. } => {
                assert!(message.contains("matrix"), "{message}");
                assert_eq!(line, Some(5));
            }
            other => panic!("wrong error: {other}"),
        }
    }

    #[test]
    fn unterminated_comment_is_an_error() {
        assert!(parse(b"#NEXUS\n[ oops\nbegin data;\n", "t").is_err());
    }

    #[test]
    fn dimension_mismatch_is_an_error() {
        let src = b"#NEXUS\nbegin data;\ndimensions ntax=2 nchar=3;\nformat datatype=dna;\nmatrix\na ACGT\nb ACGT\n;\nend;\n";
        let e = parse(src, "t").unwrap_err();
        assert!(format!("{e}").contains("nchar=3"), "{e}");

        let src = b"#NEXUS\nbegin data;\ndimensions ntax=3 nchar=4;\nformat datatype=dna;\nmatrix\na ACGT\nb ACGT\n;\nend;\n";
        let e = parse(src, "t").unwrap_err();
        assert!(format!("{e}").contains("ntax=3"), "{e}");
    }

    #[test]
    fn missing_nexus_header_is_an_error() {
        assert!(parse(b"begin data;\nend;\n", "t").is_err());
    }

    #[test]
    fn round_trips_quoting_names_that_need_it() {
        let aln = Alignment::new(
            "t",
            vec![Sequence::new("alpha", *b"ACGTACGT"), Sequence::new("Homo sapiens", *b"ACGT--GT")],
        );
        let opts = WriteOptions { sanitize_names: false, ..Default::default() };
        let text = write(&aln, &opts).unwrap();
        assert!(text.contains("'Homo sapiens'"), "{text}");
        let back = parse(text.as_bytes(), "t").unwrap();
        assert_eq!(back.sequences[1].id, "Homo sapiens");
        assert_eq!(back.sequences[1].residues, b"ACGT--GT");
    }

    #[test]
    fn interleaved_write_round_trips() {
        let aln = Alignment::new(
            "t",
            vec![Sequence::new("alpha", *b"ACGTACGTAC"), Sequence::new("beta", *b"ACGTTTGTAC")],
        );
        let opts = WriteOptions { interleaved: true, block_width: 4, ..Default::default() };
        let text = write(&aln, &opts).unwrap();
        let back = parse(text.as_bytes(), "t").unwrap();
        assert_eq!(back.sequences[0].residues, b"ACGTACGTAC");
        assert_eq!(back.sequences[1].residues, b"ACGTTTGTAC");
    }

    #[test]
    fn refuses_ragged_input() {
        let aln =
            Alignment::new("t", vec![Sequence::new("a", *b"ACGT"), Sequence::new("b", *b"AC")]);
        assert!(matches!(write(&aln, &WriteOptions::default()), Err(Error::Format(_))));
    }
}
