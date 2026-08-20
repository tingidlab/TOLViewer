//! Reading and writing the `.tolvlib` library file.
//!
//! The file is plain text, tab-indented to mirror the tree, and versioned by
//! its first line. It records *where things are*, never what is in them:
//!
//! ```text
//! TOLViewer-library 1
//! name "Lace bug project"
//! folder "Lace bug project"
//! \tfolder "18S"
//! \t\tentry "TL-2213"
//! \t\t\tfile "reads/TL-2213_18S_F.ab1"
//! \t\t\tformat ab1
//! \t\t\tkind trace
//! \t\t\treversed
//! primer "LCO1490" "GGTCAACAAATCATAAAGATATTGG"
//! ```
//!
//! Paths under the library file's own directory are stored relative to it, so a
//! project folder can be copied to another machine or a shared drive and still
//! open. Anything outside stays absolute, because rewriting it relative would
//! be a guess about a directory layout that is not ours.
//!
//! Text is quoted and backslash-escaped, so a folder called `18S "run 2"` or a
//! path with a tab in it survives the round trip. Unknown keys inside a known
//! block are skipped rather than rejected, so a library written by a later
//! version still opens here with whatever this version understands.

use std::fmt::Write as _;
use std::path::{Component, Path, PathBuf};

use tolviewer_core::{Error, Result};
use tolviewer_io::Format;

use crate::library::{Entry, EntryKind, Folder, Library, Node, NodeKind};
use crate::primer::{Primer, PrimerSet};

/// The format's name in error messages.
const FORMAT: &str = "library";
/// Bumped only when an old reader could misread a new file. Adding keys does
/// not need it, because unknown keys are skipped.
const VERSION: u32 = 1;
const MAGIC: &str = "TOLViewer-library";

/// The extension a library file gets.
pub const EXTENSION: &str = "tolvlib";

/// Render a library as the text of a `.tolvlib` file.
///
/// `at` is where the file will live, and decides which paths can be stored
/// relative.
pub fn write_string(library: &Library, at: &Path) -> String {
    let base = at.parent().unwrap_or_else(|| Path::new("."));
    let mut out = String::new();
    let _ = writeln!(out, "{MAGIC} {VERSION}");
    let _ = writeln!(out, "name {}", quote(&library.name));

    for (id, depth) in library.walk() {
        let Some(node) = library.get(id) else { continue };
        let indent = "\t".repeat(depth);
        match &node.kind {
            NodeKind::Folder(folder) => {
                let _ = writeln!(out, "{indent}folder {}", quote(&node.name));
                if !folder.expanded {
                    let _ = writeln!(out, "{indent}\tcollapsed");
                }
                if !folder.note.is_empty() {
                    let _ = writeln!(out, "{indent}\tnote {}", quote(&folder.note));
                }
            }
            NodeKind::Entry(entry) => {
                let _ = writeln!(out, "{indent}entry {}", quote(&node.name));
                let inner = format!("{indent}\t");
                let _ = writeln!(out, "{inner}file {}", quote(&relative(&entry.origin, base)));
                let _ = writeln!(out, "{inner}format {}", format_token(entry.format));
                let _ = writeln!(out, "{inner}kind {}", kind_token(entry.kind));
                if let Some(working) = &entry.working {
                    let _ = writeln!(out, "{inner}working {}", quote(&relative(working, base)));
                }
                if let Some(select) = &entry.select {
                    for id in select {
                        let _ = writeln!(out, "{inner}select {}", quote(id));
                    }
                }
                if entry.reversed {
                    let _ = writeln!(out, "{inner}reversed");
                }
                if !entry.note.is_empty() {
                    let _ = writeln!(out, "{inner}note {}", quote(&entry.note));
                }
            }
        }
    }

    for primer in library.primers.primers() {
        let sequence = String::from_utf8_lossy(&primer.sequence).into_owned();
        let _ = writeln!(out, "primer {} {}", quote(&primer.name), quote(&sequence));
    }
    out
}

/// Write a library to `path`, atomically, and mark it saved.
pub fn save(library: &mut Library, path: &Path) -> Result<()> {
    let text = write_string(library, path);
    let tmp = temp_path(path);
    match std::fs::write(&tmp, text.as_bytes()).and_then(|()| std::fs::rename(&tmp, path)) {
        Ok(()) => {
            library.path = Some(path.to_path_buf());
            library.set_saved_revision();
            Ok(())
        }
        Err(e) => {
            let _ = std::fs::remove_file(&tmp);
            Err(Error::Io(e))
        }
    }
}

/// Read a library file.
pub fn load(path: &Path) -> Result<Library> {
    let text = std::fs::read_to_string(path)?;
    let mut library = parse(&text, path)?;
    library.path = Some(path.to_path_buf());
    library.set_saved_revision();
    Ok(library)
}

/// `dir/.name.tolviewer-tmp`, so the rename onto the target is atomic.
fn temp_path(path: &Path) -> PathBuf {
    let file = path
        .file_name()
        .map(|f| f.to_string_lossy().into_owned())
        .unwrap_or_else(|| "library".to_string());
    let dir = path.parent().unwrap_or_else(|| Path::new("."));
    dir.join(format!(".{file}.tolviewer-tmp"))
}

/// Parse library text. `at` is the path it was read from, used to resolve
/// relative entry paths.
pub fn parse(text: &str, at: &Path) -> Result<Library> {
    let base = at.parent().unwrap_or_else(|| Path::new(".")).to_path_buf();
    let mut lines = text.lines().enumerate();

    let (line_no, header) =
        lines.next().ok_or_else(|| Error::parse(FORMAT, Some(1), "the file is empty"))?;
    let version = header
        .strip_prefix(MAGIC)
        .and_then(|rest| rest.trim().parse::<u32>().ok())
        .ok_or_else(|| {
            Error::parse(
                FORMAT,
                Some(line_no + 1),
                format!("this is not a TOLViewer library file (expected a '{MAGIC}' header)"),
            )
        })?;
    if version > VERSION {
        return Err(Error::parse(
            FORMAT,
            Some(line_no + 1),
            format!(
                "this library was written by a newer version of TOLViewer \
                 (format {version}, this build understands {VERSION})"
            ),
        ));
    }

    let mut library = Library::new("Untitled library");
    let mut primers = PrimerSet::default();
    // The folder each depth is currently inside, so an indented line knows its
    // parent without the file having to name it.
    let mut stack: Vec<crate::library::NodeId> = Vec::new();
    // The entry currently being filled in, and the depth its keys sit at.
    let mut pending: Option<(usize, String, EntryDraft)> = None;

    for (index, raw) in lines {
        let line_no = index + 1;
        let depth = raw.chars().take_while(|&c| c == '\t').count();
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let (key, rest) = split_key(line);

        // An entry's keys are indented one deeper than the entry itself.
        // Anything at or above the entry's own depth ends it.
        if let Some((entry_depth, _, _)) = &pending {
            if depth <= *entry_depth {
                finish_entry(&mut library, &mut pending, &mut stack, line_no)?;
            }
        }

        match key {
            "name" => library.name = unquote(rest, line_no)?,
            "folder" => {
                let name = unquote(rest, line_no)?;
                let node = Node {
                    name,
                    parent: None,
                    kind: NodeKind::Folder(Folder { expanded: true, ..Folder::default() }),
                };
                library.push_at_depth(depth, node, &mut stack);
            }
            "entry" => {
                let name = unquote(rest, line_no)?;
                pending = Some((depth, name, EntryDraft::default()));
            }
            "primer" => {
                let (name, sequence) = two_strings(rest, line_no)?;
                primers.push(
                    Primer::new(name, &sequence)
                        .map_err(|e| Error::parse(FORMAT, Some(line_no), e.to_string()))?,
                );
            }
            // Folder-level keys.
            "collapsed" => {
                if let Some(folder) =
                    stack.last().and_then(|&id| library.get_mut(id)).and_then(|n| n.folder_mut())
                {
                    folder.expanded = false;
                }
            }
            // Entry-level keys, and `note`, which both kinds have.
            other => {
                if let Some((_, _, draft)) = &mut pending {
                    draft.set(other, rest, &base, line_no)?;
                } else if other == "note" {
                    let note = unquote(rest, line_no)?;
                    if let Some(folder) = stack
                        .last()
                        .and_then(|&id| library.get_mut(id))
                        .and_then(|n| n.folder_mut())
                    {
                        folder.note = note;
                    }
                }
                // Anything else is a key from a newer version: skip it.
            }
        }
    }
    finish_entry(&mut library, &mut pending, &mut stack, text.lines().count())?;
    library.primers = primers;
    Ok(library)
}

/// Attach the entry that was being read, if there is one.
fn finish_entry(
    library: &mut Library,
    pending: &mut Option<(usize, String, EntryDraft)>,
    stack: &mut Vec<crate::library::NodeId>,
    line_no: usize,
) -> Result<()> {
    let Some((depth, name, draft)) = pending.take() else { return Ok(()) };
    let entry = draft.build(&name, line_no)?;
    let node = Node { name, parent: None, kind: NodeKind::Entry(Box::new(entry)) };
    library.push_at_depth(depth, node, stack);
    Ok(())
}

/// An entry's keys as they are read, before it is known whether they are
/// complete.
#[derive(Default)]
struct EntryDraft {
    origin: Option<PathBuf>,
    format: Option<Format>,
    working: Option<PathBuf>,
    select: Vec<String>,
    reversed: bool,
    kind: Option<EntryKind>,
    note: String,
}

impl EntryDraft {
    fn set(&mut self, key: &str, rest: &str, base: &Path, line_no: usize) -> Result<()> {
        match key {
            "file" => self.origin = Some(absolute(&unquote(rest, line_no)?, base)),
            "working" => self.working = Some(absolute(&unquote(rest, line_no)?, base)),
            "select" => self.select.push(unquote(rest, line_no)?),
            "reversed" => self.reversed = true,
            "note" => self.note = unquote(rest, line_no)?,
            "format" => {
                self.format = Some(format_from_token(rest.trim()).ok_or_else(|| {
                    Error::parse(FORMAT, Some(line_no), format!("unknown format '{}'", rest.trim()))
                })?)
            }
            "kind" => self.kind = kind_from_token(rest.trim()),
            // A key this version does not know; a newer TOLViewer wrote it.
            _ => {}
        }
        Ok(())
    }

    fn build(self, name: &str, line_no: usize) -> Result<Entry> {
        let origin = self.origin.ok_or_else(|| {
            Error::parse(
                FORMAT,
                Some(line_no),
                format!("entry '{name}' does not say which file it is"),
            )
        })?;
        // The format is only a hint about a file we do not own; if the line is
        // missing, fall back to the extension rather than refusing to open the
        // whole library over one entry.
        let format = self.format.or_else(|| Format::from_path(&origin)).unwrap_or(Format::Fasta);
        Ok(Entry {
            origin,
            format,
            working: self.working,
            select: (!self.select.is_empty()).then_some(self.select),
            reversed: self.reversed,
            kind: self.kind.unwrap_or(EntryKind::Sequences),
            note: self.note,
        })
    }
}

/// The first word and the remainder.
fn split_key(line: &str) -> (&str, &str) {
    match line.split_once(char::is_whitespace) {
        Some((key, rest)) => (key, rest.trim_start()),
        None => (line, ""),
    }
}

/// Two quoted strings on one line, as `primer` uses.
fn two_strings(rest: &str, line_no: usize) -> Result<(String, String)> {
    let (first, remainder) = take_quoted(rest, line_no)?;
    let (second, _) = take_quoted(remainder.trim_start(), line_no)?;
    Ok((first, second))
}

fn unquote(rest: &str, line_no: usize) -> Result<String> {
    Ok(take_quoted(rest, line_no)?.0)
}

/// Read one `"..."` string, returning it and whatever followed.
///
/// An unquoted value is accepted as the rest of the line, so a hand-edited file
/// with `name Lace bug project` still opens.
fn take_quoted(rest: &str, line_no: usize) -> Result<(String, &str)> {
    let rest = rest.trim_start();
    let Some(body) = rest.strip_prefix('"') else {
        return Ok((rest.trim().to_string(), ""));
    };
    let mut out = String::with_capacity(body.len());
    let mut chars = body.char_indices();
    while let Some((i, c)) = chars.next() {
        match c {
            '"' => return Ok((out, &body[i + 1..])),
            '\\' => match chars.next() {
                Some((_, 'n')) => out.push('\n'),
                Some((_, 't')) => out.push('\t'),
                Some((_, 'r')) => out.push('\r'),
                Some((_, other)) => out.push(other),
                None => break,
            },
            other => out.push(other),
        }
    }
    Err(Error::parse(FORMAT, Some(line_no), "a quoted value is missing its closing quote"))
}

fn quote(value: &str) -> String {
    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for c in value.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\t' => out.push_str("\\t"),
            '\r' => out.push_str("\\r"),
            other => out.push(other),
        }
    }
    out.push('"');
    out
}

/// `path` relative to `base` when it sits under it, and absolute otherwise.
///
/// Only a plain prefix match is used: walking up out of `base` with `..` would
/// produce something longer and more fragile than the absolute path it replaced.
fn relative(path: &Path, base: &Path) -> String {
    let text = match path.strip_prefix(base) {
        Ok(rest) if rest.as_os_str().is_empty() => path.to_string_lossy().into_owned(),
        // Always store `/` so a library written on Windows opens on Linux.
        Ok(rest) => {
            rest.components().map(|c| c.as_os_str().to_string_lossy()).collect::<Vec<_>>().join("/")
        }
        Err(_) => path.to_string_lossy().into_owned(),
    };
    text
}

/// Resolve a stored path against the library's directory.
fn absolute(stored: &str, base: &Path) -> PathBuf {
    let path = PathBuf::from(stored.replace('/', std::path::MAIN_SEPARATOR_STR));
    if path.components().next() == Some(Component::CurDir) || path.is_relative() {
        base.join(path)
    } else {
        path
    }
}

/// A stable on-disk name for a format.
///
/// Deliberately not [`Format::name`], which is a menu label and may be
/// reworded. The match is exhaustive, so adding a format is a compile error
/// here rather than a file that cannot be reopened.
fn format_token(format: Format) -> &'static str {
    match format {
        Format::Fasta => "fasta",
        Format::Fastq => "fastq",
        Format::Phylip => "phylip",
        Format::PhylipRelaxed => "phylip-relaxed",
        Format::Nexus => "nexus",
        Format::Clustal => "clustal",
        Format::Stockholm => "stockholm",
        Format::Msf => "msf",
        Format::Genbank => "genbank",
        Format::Ab1 => "ab1",
    }
}

fn format_from_token(token: &str) -> Option<Format> {
    Format::all().iter().copied().find(|&f| format_token(f) == token)
}

fn kind_token(kind: EntryKind) -> &'static str {
    match kind {
        EntryKind::Sequences => "sequences",
        EntryKind::Alignment => "alignment",
        EntryKind::Trace => "trace",
    }
}

fn kind_from_token(token: &str) -> Option<EntryKind> {
    match token {
        "sequences" => Some(EntryKind::Sequences),
        "alignment" => Some(EntryKind::Alignment),
        "trace" => Some(EntryKind::Trace),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::library::NodeId;

    fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("tolviewer-store-{name}"));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn fasta(dir: &Path, name: &str) -> PathBuf {
        let path = dir.join(name);
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        std::fs::write(&path, ">a\nACGT\n").unwrap();
        path
    }

    /// The library from the request, with a file under each locus.
    fn project(dir: &Path) -> (Library, PathBuf) {
        let mut lib = Library::new("Lace bug project");
        let project = lib.add_folder(None, "Lace bug project").unwrap();
        let ssu = lib.add_folder(Some(project), "18S").unwrap();
        let lsu = lib.add_folder(Some(project), "28S").unwrap();
        lib.add_file(Some(ssu), &fasta(dir, "reads/TL-2213_18S.fasta")).unwrap();
        lib.add_file(Some(lsu), &fasta(dir, "reads/TL-2213_28S.fasta")).unwrap();
        lib.get_mut(lsu).unwrap().folder_mut().unwrap().expanded = false;
        lib.primers.push(Primer::new("LCO1490", "GGTCAACAAATCATAAAGATATTGG").unwrap());
        (lib, dir.join("project.tolvlib"))
    }

    /// Compare two libraries by everything the file is supposed to carry.
    fn same_shape(a: &Library, b: &Library) {
        assert_eq!(a.name, b.name);
        assert_eq!(a.primers, b.primers);
        let (aw, bw) = (a.walk(), b.walk());
        assert_eq!(aw.len(), bw.len(), "different number of nodes");
        for ((ai, ad), (bi, bd)) in aw.iter().zip(&bw) {
            assert_eq!(ad, bd, "depth differs");
            let (an, bn) = (a.get(*ai).unwrap(), b.get(*bi).unwrap());
            assert_eq!(an.name, bn.name);
            assert_eq!(an.entry(), bn.entry());
            assert_eq!(
                an.folder().map(|f| (f.expanded, f.note.clone())),
                bn.folder().map(|f| (f.expanded, f.note.clone()))
            );
        }
    }

    #[test]
    fn a_project_round_trips_through_the_file() {
        let dir = scratch("roundtrip");
        let (lib, at) = project(&dir);
        let text = write_string(&lib, &at);
        let back = parse(&text, &at).unwrap();
        same_shape(&lib, &back);
        assert_eq!(back.walk().len(), 5);
    }

    #[test]
    fn saving_and_loading_goes_through_disk_intact() {
        let dir = scratch("disk");
        let (mut lib, at) = project(&dir);
        save(&mut lib, &at).unwrap();
        assert!(!lib.is_dirty());
        assert_eq!(lib.path.as_deref(), Some(at.as_path()));
        assert!(!dir.join(".project.tolvlib.tolviewer-tmp").exists(), "no temp file left");

        let back = load(&at).unwrap();
        same_shape(&lib, &back);
        assert!(!back.is_dirty(), "a library just loaded has nothing to save");
        // The entries still resolve to real files.
        assert!(back.broken_entries().is_empty());
        let entry = back.entries_under(None)[0];
        assert_eq!(back.load(entry).unwrap().len(), 1);
    }

    #[test]
    fn paths_inside_the_project_are_stored_relative_so_it_can_move() {
        let dir = scratch("relative");
        let (mut lib, at) = project(&dir);
        save(&mut lib, &at).unwrap();
        let text = std::fs::read_to_string(&at).unwrap();
        assert!(text.contains("\"reads/TL-2213_18S.fasta\""), "{text}");
        assert!(!text.contains(dir.to_str().unwrap()), "an absolute path leaked in");

        // Move the whole project directory and it still opens.
        let moved = scratch("relative-moved");
        for name in ["reads/TL-2213_18S.fasta", "reads/TL-2213_28S.fasta"] {
            let to = moved.join(name);
            std::fs::create_dir_all(to.parent().unwrap()).unwrap();
            std::fs::copy(dir.join(name), to).unwrap();
        }
        let moved_at = moved.join("project.tolvlib");
        std::fs::copy(&at, &moved_at).unwrap();
        let back = load(&moved_at).unwrap();
        assert!(back.broken_entries().is_empty(), "the moved project did not resolve");
        assert_eq!(
            back.entry(back.entries_under(None)[0]).unwrap().origin,
            moved.join("reads/TL-2213_18S.fasta")
        );
    }

    #[test]
    fn a_file_outside_the_project_keeps_its_absolute_path() {
        let dir = scratch("outside");
        let elsewhere = scratch("outside-archive");
        let far = fasta(&elsewhere, "archive.fasta");
        let mut lib = Library::new("l");
        lib.add_file(None, &far).unwrap();
        let at = dir.join("p.tolvlib");
        let text = write_string(&lib, &at);
        assert!(text.contains(far.to_str().unwrap()), "{text}");
        let back = parse(&text, &at).unwrap();
        assert_eq!(back.entry(back.entries_under(None)[0]).unwrap().origin, far);
    }

    #[test]
    fn an_entrys_working_copy_and_selection_survive() {
        let dir = scratch("working");
        let mut lib = Library::new("l");
        let source = lib.add_file(None, &fasta(&dir, "msa.fasta")).unwrap();
        let one = lib.add_selection(None, source, vec!["a".into()], "just a").unwrap();
        lib.set_reversed(one, true).unwrap();
        lib.get_mut(one).unwrap().entry_mut().unwrap().working = Some(dir.join("a.edited.fasta"));
        lib.get_mut(one).unwrap().entry_mut().unwrap().note = "checked by eye".into();

        let at = dir.join("p.tolvlib");
        let back = parse(&write_string(&lib, &at), &at).unwrap();
        let entry = back.entry(back.entries_under(None)[1]).unwrap();
        assert_eq!(entry.select.as_deref(), Some(&["a".to_string()][..]));
        assert_eq!(entry.working, Some(dir.join("a.edited.fasta")));
        assert!(entry.reversed);
        assert_eq!(entry.note, "checked by eye");
    }

    #[test]
    fn awkward_names_survive_quoting() {
        let dir = scratch("quoting");
        let mut lib = Library::new("18S \"run 2\"\tand\\3");
        let folder = lib.add_folder(None, "a\nb").unwrap();
        lib.get_mut(folder).unwrap().folder_mut().unwrap().note = "tab\there".into();
        let at = dir.join("p.tolvlib");
        let back = parse(&write_string(&lib, &at), &at).unwrap();
        assert_eq!(back.name, "18S \"run 2\"\tand\\3");
        assert_eq!(back.get(back.roots()[0]).unwrap().name, "a\nb");
        assert_eq!(back.get(back.roots()[0]).unwrap().folder().unwrap().note, "tab\there");
    }

    #[test]
    fn collapsed_folders_reopen_collapsed() {
        let dir = scratch("collapsed");
        let (lib, at) = project(&dir);
        let back = parse(&write_string(&lib, &at), &at).unwrap();
        let states: Vec<bool> = back
            .walk()
            .iter()
            .filter_map(|(id, _)| back.get(*id).unwrap().folder().map(|f| f.expanded))
            .collect();
        assert_eq!(states, vec![true, true, false]);
    }

    #[test]
    fn every_format_has_a_token_that_maps_back() {
        for &f in Format::all() {
            let token = format_token(f);
            assert_eq!(format_from_token(token), Some(f), "{token}");
            assert!(!token.contains(' '), "tokens must be single words: {token}");
        }
        assert_eq!(format_from_token("nonsense"), None);
        for k in [EntryKind::Sequences, EntryKind::Alignment, EntryKind::Trace] {
            assert_eq!(kind_from_token(kind_token(k)), Some(k));
        }
    }

    #[test]
    fn a_file_that_is_not_a_library_is_refused_by_name() {
        let at = Path::new("/tmp/p.tolvlib");
        let e = parse(">a\nACGT\n", at).unwrap_err();
        assert!(e.to_string().contains("not a TOLViewer library"), "{e}");
        assert!(parse("", at).is_err());
    }

    #[test]
    fn a_newer_format_version_is_refused_with_an_explanation() {
        let at = Path::new("/tmp/p.tolvlib");
        let e = parse(&format!("{MAGIC} 99\nname \"x\"\n"), at).unwrap_err();
        assert!(e.to_string().contains("newer version"), "{e}");
    }

    #[test]
    fn unknown_keys_from_a_future_version_are_skipped_not_fatal() {
        let at = Path::new("/tmp/p.tolvlib");
        let text = format!(
            "{MAGIC} 1\nname \"p\"\ncolour \"blue\"\nfolder \"f\"\n\tstarred\n\tentry \"e\"\n\
             \t\tfile \"/tmp/x.fasta\"\n\t\tformat fasta\n\t\tconfidence 0.9\n"
        );
        let lib = parse(&text, at).unwrap();
        assert_eq!(lib.name, "p");
        assert_eq!(lib.walk().len(), 2);
        assert_eq!(
            lib.entry(lib.entries_under(None)[0]).unwrap().origin,
            PathBuf::from("/tmp/x.fasta")
        );
    }

    #[test]
    fn an_entry_with_no_file_is_a_clear_error() {
        let at = Path::new("/tmp/p.tolvlib");
        let text = format!("{MAGIC} 1\nentry \"nameless\"\n\tformat fasta\n");
        let e = parse(&text, at).unwrap_err();
        assert!(e.to_string().contains("nameless"), "{e}");
        assert!(e.to_string().contains("which file"), "{e}");
    }

    #[test]
    fn an_unterminated_quote_is_reported_with_its_line() {
        let at = Path::new("/tmp/p.tolvlib");
        let e = parse(&format!("{MAGIC} 1\nname \"unfinished\n"), at).unwrap_err();
        match e {
            Error::Parse { line, message, .. } => {
                assert_eq!(line, Some(2));
                assert!(message.contains("closing quote"), "{message}");
            }
            other => panic!("wrong error: {other}"),
        }
    }

    #[test]
    fn a_hand_written_file_without_quotes_still_opens() {
        let at = Path::new("/tmp/p.tolvlib");
        let text = format!("{MAGIC} 1\nname Lace bug project\nfolder 18S\n");
        let lib = parse(&text, at).unwrap();
        assert_eq!(lib.name, "Lace bug project");
        assert_eq!(lib.get(lib.roots()[0]).unwrap().name, "18S");
    }

    #[test]
    fn blank_lines_and_comments_are_ignored() {
        let at = Path::new("/tmp/p.tolvlib");
        let text = format!("{MAGIC} 1\n\n# written by hand\nname \"p\"\n\nfolder \"f\"\n");
        let lib = parse(&text, at).unwrap();
        assert_eq!(lib.walk().len(), 1);
    }

    #[test]
    fn primers_round_trip_and_a_bad_one_names_its_line() {
        let dir = scratch("primers");
        let (lib, at) = project(&dir);
        let back = parse(&write_string(&lib, &at), &at).unwrap();
        assert_eq!(back.primers.len(), 1);
        assert_eq!(back.primers.get(0).unwrap().name, "LCO1490");

        let bad = format!("{MAGIC} 1\nprimer \"p\" \"GGZZ\"\n");
        let e = parse(&bad, &at).unwrap_err();
        match e {
            Error::Parse { line, .. } => assert_eq!(line, Some(2)),
            other => panic!("wrong error: {other}"),
        }
    }

    #[test]
    fn deep_nesting_keeps_its_shape() {
        let dir = scratch("deep");
        let mut lib = Library::new("l");
        let mut parent: Option<NodeId> = None;
        for level in 0..6 {
            parent = Some(lib.add_folder(parent, format!("level {level}")).unwrap());
        }
        lib.add_file(parent, &fasta(&dir, "deep.fasta")).unwrap();
        let at = dir.join("p.tolvlib");
        let back = parse(&write_string(&lib, &at), &at).unwrap();
        same_shape(&lib, &back);
        assert_eq!(back.walk().last().unwrap().1, 6);
    }

    #[test]
    fn an_empty_library_round_trips() {
        let at = Path::new("/tmp/p.tolvlib");
        let lib = Library::new("nothing yet");
        let back = parse(&write_string(&lib, at), at).unwrap();
        assert_eq!(back.name, "nothing yet");
        assert!(back.is_empty());
    }
}
