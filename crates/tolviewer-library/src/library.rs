//! The project library: folders, subfolders and the files they point at.
//!
//! A library is a tree the user arranges — "Lace bug project" with "18S" and
//! "28S" under it — whose leaves are *references* to sequence files. The files
//! themselves are never moved, copied or rewritten by adding them to a library:
//! a library is an index over data that stays where the sequencing facility put
//! it. That is what makes it safe to point one at a read-only archive directory
//! or a shared drive.
//!
//! Because the files are the lab's originals, writing to one is treated as a
//! decision rather than a default. See [`SaveTarget`] and
//! [`Library::save_entry`].
//!
//! The tree is an arena: nodes live in one `Vec` and refer to each other by
//! [`NodeId`]. Ids are never reused, so a stale id from the GUI resolves to
//! `None` instead of silently addressing whatever was allocated in its place.

use std::path::{Path, PathBuf};

use tolviewer_core::{Alignment, Error, Result, Sequence};
use tolviewer_io::{ab1, Format, WriteOptions};

/// A handle to a folder or entry. Stable for the life of the library.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NodeId(u32);

impl NodeId {
    /// The raw number, for persisting selection state. Only meaningful within
    /// one library.
    pub fn raw(self) -> u32 {
        self.0
    }
}

/// What an entry's file holds, which decides what opening it does.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EntryKind {
    /// One or more sequences that are not aligned to each other.
    Sequences,
    /// Rows of equal length: an alignment.
    Alignment,
    /// A Sanger chromatogram, which opens in the trace viewer.
    Trace,
}

impl EntryKind {
    pub fn name(self) -> &'static str {
        match self {
            EntryKind::Sequences => "sequences",
            EntryKind::Alignment => "alignment",
            EntryKind::Trace => "trace",
        }
    }

    /// Work out what a file holds, from its format and its contents.
    pub fn of(format: Format, alignment: &Alignment) -> EntryKind {
        if format == Format::Ab1 {
            EntryKind::Trace
        } else if alignment.len() > 1 && alignment.is_aligned() {
            EntryKind::Alignment
        } else {
            EntryKind::Sequences
        }
    }
}

/// A leaf: one named reference to sequence data on disk.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Entry {
    /// The file this entry came from. The library treats it as the lab's data
    /// and never writes here without an explicit [`SaveChoice::Overwrite`].
    pub origin: PathBuf,
    pub format: Format,
    /// Set the first time the user diverts edits to a copy. Every later save
    /// goes straight here, so a sequence edited ten times leaves one extra
    /// file behind rather than ten.
    pub working: Option<PathBuf>,
    /// Which sequences of the file this entry stands for, by id. `None` means
    /// the whole file, which is the usual case; `Some` is what extracting a row
    /// out of an alignment produces.
    pub select: Option<Vec<String>>,
    /// Show the sequences reverse complemented. Stored as a flag because
    /// flipping a read is cheap, lossless and frequently undone — rewriting the
    /// file for it would be both slower and more destructive.
    pub reversed: bool,
    pub kind: EntryKind,
    /// The operator's note, shown on hover.
    pub note: String,
}

impl Entry {
    /// The file a read would come from: the working copy once there is one.
    pub fn source_path(&self) -> &Path {
        self.working.as_deref().unwrap_or(&self.origin)
    }

    /// Read this entry's sequences, applying its selection and orientation.
    ///
    /// Reading always goes through [`Entry::source_path`], so once edits have
    /// been diverted to a copy the original is never touched again.
    pub fn load(&self) -> Result<Alignment> {
        let path = self.source_path();
        let mut alignment = tolviewer_io::read_file_as(path, self.effective_format())?;
        if let Some(wanted) = &self.select {
            let rows: Vec<usize> =
                wanted.iter().filter_map(|id| alignment.find_by_id(id)).collect();
            if rows.is_empty() {
                return Err(Error::format(format!(
                    "none of the sequences this entry refers to are still in {}; \
                     it may have been rewritten by something else",
                    path.display()
                )));
            }
            let width = alignment.width();
            alignment = alignment.subset(&rows, 0..width);
        }
        if self.reversed {
            let alphabet = alignment.alphabet();
            if !alphabet.is_nucleotide() {
                return Err(Error::format(format!(
                    "{} is marked reversed but reads as {}; \
                     only nucleotide sequences have a reverse complement",
                    path.display(),
                    alphabet.name()
                )));
            }
            for seq in &mut alignment.sequences {
                seq.reverse_complement(alphabet);
            }
        }
        Ok(alignment)
    }

    /// Read the chromatogram, oriented the way the entry is.
    ///
    /// Fails for anything that is not a trace: there is no signal to show.
    pub fn load_trace(&self) -> Result<ab1::Trace> {
        if self.kind != EntryKind::Trace || self.effective_format() != Format::Ab1 {
            return Err(Error::format(format!(
                "{} is not a trace file, so it has no chromatogram",
                self.source_path().display()
            )));
        }
        let mut trace = ab1::read_file(self.source_path())?;
        if self.reversed {
            trace.reverse_complement();
        }
        Ok(trace)
    }

    /// The format to read the current source with.
    ///
    /// A working copy is always written as FASTA, because the formats an entry
    /// can arrive in include several that cannot be written back.
    fn effective_format(&self) -> Format {
        match &self.working {
            Some(path) => Format::from_path(path).unwrap_or(Format::Fasta),
            None => self.format,
        }
    }

    /// Where a save would go, and what it would cost.
    pub fn save_target(&self) -> SaveTarget {
        if let Some(path) = &self.working {
            return SaveTarget::WorkingCopy(path.clone());
        }
        // A slice of a file cannot be written back over the whole file without
        // destroying its siblings, and a read-only format cannot be written at
        // all. Either way only a copy will do.
        if self.select.is_some() {
            return SaveTarget::MustCopy(self.suggested_copy(), CopyReason::PartOfAFile);
        }
        if !self.format.can_write() {
            return SaveTarget::MustCopy(self.suggested_copy(), CopyReason::ReadOnlyFormat);
        }
        SaveTarget::Original(self.origin.clone())
    }

    /// A path for the working copy: the origin's name with `.edited` before the
    /// extension, made unique against what is already on disk.
    ///
    /// Traces and other read-only formats become FASTA, since their format
    /// cannot be written.
    pub fn suggested_copy(&self) -> PathBuf {
        let dir = self.origin.parent().unwrap_or_else(|| Path::new("."));
        let stem = self
            .origin
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "sequences".to_string());
        let ext = if self.format.can_write() {
            self.origin
                .extension()
                .map(|e| e.to_string_lossy().into_owned())
                .unwrap_or_else(|| self.format.extensions()[0].to_string())
        } else {
            Format::Fasta.extensions()[0].to_string()
        };
        let mut candidate = dir.join(format!("{stem}.edited.{ext}"));
        let mut n = 2;
        // Never propose a name that would clobber something. The loop is
        // bounded because each try is a different name.
        while candidate.exists() {
            candidate = dir.join(format!("{stem}.edited-{n}.{ext}"));
            n += 1;
            if n > 1000 {
                break;
            }
        }
        candidate
    }
}

/// Why an entry cannot be written back over its origin.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CopyReason {
    /// The entry is some of a file's sequences, not all of them.
    PartOfAFile,
    /// The origin is in a format this program can read but not write.
    ReadOnlyFormat,
}

impl CopyReason {
    /// A sentence for the dialog, explaining why there is no "overwrite"
    /// button.
    pub fn explain(self) -> &'static str {
        match self {
            CopyReason::PartOfAFile => {
                "This entry is part of a larger file. Writing it back would \
                 discard the other sequences in that file, so it has to be \
                 saved separately."
            }
            CopyReason::ReadOnlyFormat => {
                "TOLViewer can read this file's format but not write it, so the \
                 edits have to go somewhere else."
            }
        }
    }
}

/// Where saving an entry would put the data.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SaveTarget {
    /// A copy the library already made. Nothing original is at risk, so this
    /// needs no confirmation.
    WorkingCopy(PathBuf),
    /// The file the entry was imported from. Overwriting it replaces the lab's
    /// data, so the user is asked first.
    Original(PathBuf),
    /// Only a copy is possible, at the suggested path.
    MustCopy(PathBuf, CopyReason),
}

impl SaveTarget {
    /// Does the user have to be asked before this write happens?
    pub fn needs_confirmation(&self) -> bool {
        !matches!(self, SaveTarget::WorkingCopy(_))
    }

    /// May the user choose to overwrite, or is a copy the only option?
    pub fn can_overwrite(&self) -> bool {
        !matches!(self, SaveTarget::MustCopy(..))
    }

    pub fn path(&self) -> &Path {
        match self {
            SaveTarget::WorkingCopy(p) | SaveTarget::Original(p) | SaveTarget::MustCopy(p, _) => p,
        }
    }
}

/// The user's answer to the overwrite question.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SaveChoice {
    /// Replace the origin file. Refused when the target says it cannot be done.
    Overwrite,
    /// Write to a copy and remember it, so later saves go there without asking
    /// again.
    NewCopy(PathBuf),
}

/// A folder in the tree.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Folder {
    pub children: Vec<NodeId>,
    /// Whether the GUI draws it open. Persisted so a library reopens looking
    /// the way it was left.
    pub expanded: bool,
    pub note: String,
}

/// A folder or an entry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NodeKind {
    Folder(Folder),
    Entry(Box<Entry>),
}

/// One node of the tree.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Node {
    pub name: String,
    pub parent: Option<NodeId>,
    pub kind: NodeKind,
}

impl Node {
    pub fn is_folder(&self) -> bool {
        matches!(self.kind, NodeKind::Folder(_))
    }

    pub fn entry(&self) -> Option<&Entry> {
        match &self.kind {
            NodeKind::Entry(e) => Some(e),
            NodeKind::Folder(_) => None,
        }
    }

    pub fn entry_mut(&mut self) -> Option<&mut Entry> {
        match &mut self.kind {
            NodeKind::Entry(e) => Some(e),
            NodeKind::Folder(_) => None,
        }
    }

    pub fn folder(&self) -> Option<&Folder> {
        match &self.kind {
            NodeKind::Folder(f) => Some(f),
            NodeKind::Entry(_) => None,
        }
    }

    pub fn folder_mut(&mut self) -> Option<&mut Folder> {
        match &mut self.kind {
            NodeKind::Folder(f) => Some(f),
            NodeKind::Entry(_) => None,
        }
    }
}

/// A project library.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Library {
    /// What the library is called, shown at the root of the tree.
    pub name: String,
    /// The project's primers, used by the mapping and trimming commands.
    pub primers: crate::primer::PrimerSet,
    /// Where the library file itself lives. Entry paths are stored relative to
    /// this directory when they sit under it, so a project folder can be moved
    /// or shared without breaking.
    pub path: Option<PathBuf>,
    nodes: Vec<Option<Node>>,
    roots: Vec<NodeId>,
    /// Bumped by every change, so the GUI can tell whether the library needs
    /// saving without comparing trees.
    revision: u64,
    saved_revision: u64,
}

impl Default for Library {
    fn default() -> Self {
        Library::new("Untitled library")
    }
}

impl Library {
    pub fn new(name: impl Into<String>) -> Self {
        Library {
            name: name.into(),
            primers: crate::primer::PrimerSet::default(),
            path: None,
            nodes: Vec::new(),
            roots: Vec::new(),
            revision: 0,
            saved_revision: 0,
        }
    }

    // ---- reading the tree ----------------------------------------------

    /// The top-level folders and entries, in order.
    pub fn roots(&self) -> &[NodeId] {
        &self.roots
    }

    pub fn get(&self, id: NodeId) -> Option<&Node> {
        self.nodes.get(id.0 as usize).and_then(|n| n.as_ref())
    }

    pub fn get_mut(&mut self, id: NodeId) -> Option<&mut Node> {
        self.touch();
        self.nodes.get_mut(id.0 as usize).and_then(|n| n.as_mut())
    }

    pub fn entry(&self, id: NodeId) -> Option<&Entry> {
        self.get(id)?.entry()
    }

    /// The children of a folder, or the roots when `parent` is `None`.
    pub fn children(&self, parent: Option<NodeId>) -> &[NodeId] {
        match parent {
            None => &self.roots,
            Some(id) => match self.get(id).and_then(|n| n.folder()) {
                Some(f) => &f.children,
                None => &[],
            },
        }
    }

    /// Every node, depth first, as `(id, depth)`. This is the order the tree is
    /// drawn and saved in.
    pub fn walk(&self) -> Vec<(NodeId, usize)> {
        let mut out = Vec::new();
        for &root in &self.roots {
            self.walk_into(root, 0, &mut out);
        }
        out
    }

    fn walk_into(&self, id: NodeId, depth: usize, out: &mut Vec<(NodeId, usize)>) {
        out.push((id, depth));
        if let Some(folder) = self.get(id).and_then(|n| n.folder()) {
            for &child in &folder.children {
                self.walk_into(child, depth + 1, out);
            }
        }
    }

    /// Every entry under `id` (or the whole library when `None`), depth first.
    pub fn entries_under(&self, id: Option<NodeId>) -> Vec<NodeId> {
        let mut out = Vec::new();
        let start: Vec<NodeId> = match id {
            Some(id) if self.get(id).is_some_and(|n| !n.is_folder()) => vec![id],
            other => self.children(other).to_vec(),
        };
        for node in start {
            self.collect_entries(node, &mut out);
        }
        out
    }

    fn collect_entries(&self, id: NodeId, out: &mut Vec<NodeId>) {
        match self.get(id).map(|n| &n.kind) {
            Some(NodeKind::Entry(_)) => out.push(id),
            Some(NodeKind::Folder(f)) => {
                for &child in &f.children {
                    self.collect_entries(child, out);
                }
            }
            None => {}
        }
    }

    /// `Lace bug project / 18S / TL-2213`, for tooltips and error messages.
    pub fn path_of(&self, id: NodeId) -> String {
        let mut parts = Vec::new();
        let mut at = Some(id);
        while let Some(current) = at {
            match self.get(current) {
                Some(node) => {
                    parts.push(node.name.clone());
                    at = node.parent;
                }
                None => break,
            }
        }
        parts.reverse();
        parts.join(" / ")
    }

    pub fn len(&self) -> usize {
        self.nodes.iter().filter(|n| n.is_some()).count()
    }

    pub fn is_empty(&self) -> bool {
        self.roots.is_empty()
    }

    pub fn revision(&self) -> u64 {
        self.revision
    }

    /// Has anything changed since the library was last saved or loaded?
    pub fn is_dirty(&self) -> bool {
        self.revision != self.saved_revision
    }

    pub fn mark_saved(&mut self) {
        self.saved_revision = self.revision;
    }

    /// Mark the library changed.
    ///
    /// Every method here calls this for itself; it is public because
    /// [`Library::primers`] is a plain field that callers edit directly, and an
    /// added primer has to count as something worth saving.
    pub fn touch(&mut self) {
        self.revision += 1;
    }

    // ---- building the tree ---------------------------------------------

    /// Add a folder under `parent`, or at the top level when `parent` is
    /// `None`.
    pub fn add_folder(
        &mut self,
        parent: Option<NodeId>,
        name: impl Into<String>,
    ) -> Result<NodeId> {
        self.check_folder(parent)?;
        let node = Node {
            name: name.into(),
            parent,
            kind: NodeKind::Folder(Folder { expanded: true, ..Folder::default() }),
        };
        Ok(self.attach(parent, node))
    }

    /// Add an entry for a file already on disk, reading it once to work out
    /// what it holds.
    ///
    /// The file is not copied or modified; the library only remembers where it
    /// is.
    pub fn add_file(&mut self, parent: Option<NodeId>, path: &Path) -> Result<NodeId> {
        self.check_folder(parent)?;
        let format = tolviewer_io::sniff_file(path)?;
        let alignment = tolviewer_io::read_file_as(path, format)?;
        let kind = EntryKind::of(format, &alignment);
        // A trace is named for the sample the operator typed, which is more
        // use than the file name the facility generated.
        let name = match (kind, alignment.sequences.first()) {
            (EntryKind::Trace, Some(seq)) if !seq.id.is_empty() => seq.id.clone(),
            _ => path
                .file_stem()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_else(|| "sequences".to_string()),
        };
        let entry = Entry {
            origin: path.to_path_buf(),
            format,
            working: None,
            select: None,
            reversed: false,
            kind,
            note: String::new(),
        };
        Ok(self.attach(parent, Node { name, parent, kind: NodeKind::Entry(Box::new(entry)) }))
    }

    /// Add an entry standing for named rows of an existing file: what
    /// extracting sequences out of an alignment produces.
    ///
    /// The alignment is left alone. The new entry reads its rows back out of
    /// the same file on demand, so nothing is duplicated on disk and nothing
    /// can be overwritten by accident — saving an edit to it will insist on a
    /// copy.
    pub fn add_selection(
        &mut self,
        parent: Option<NodeId>,
        source: NodeId,
        ids: Vec<String>,
        name: impl Into<String>,
    ) -> Result<NodeId> {
        self.check_folder(parent)?;
        if ids.is_empty() {
            return Err(Error::format("select at least one sequence to extract"));
        }
        let entry = self
            .entry(source)
            .ok_or_else(|| Error::out_of_range("that entry is no longer in the library"))?;
        let extracted = Entry {
            origin: entry.origin.clone(),
            format: entry.format,
            // The extract reads from wherever the source currently reads from,
            // so an extract taken from an edited alignment sees the edits.
            working: entry.working.clone(),
            select: Some(ids),
            reversed: entry.reversed,
            kind: EntryKind::Sequences,
            note: String::new(),
        };
        Ok(self.attach(
            parent,
            Node { name: name.into(), parent, kind: NodeKind::Entry(Box::new(extracted)) },
        ))
    }

    /// Refuse to file something under an entry: only folders have children.
    fn check_folder(&self, parent: Option<NodeId>) -> Result<()> {
        match parent {
            None => Ok(()),
            Some(id) => match self.get(id) {
                Some(node) if node.is_folder() => Ok(()),
                Some(node) => Err(Error::format(format!(
                    "'{}' is a sequence, not a folder, so nothing can go inside it",
                    node.name
                ))),
                None => Err(Error::out_of_range("that folder is no longer in the library")),
            },
        }
    }

    fn attach(&mut self, parent: Option<NodeId>, node: Node) -> NodeId {
        let id = NodeId(self.nodes.len() as u32);
        self.nodes.push(Some(node));
        match parent {
            Some(p) => {
                if let Some(folder) = self.nodes[p.0 as usize].as_mut().and_then(|n| n.folder_mut())
                {
                    folder.children.push(id);
                }
            }
            None => self.roots.push(id),
        }
        self.touch();
        id
    }

    /// Remove a node and everything under it. The files it referred to are not
    /// touched — removing something from a library is filing, not deleting.
    ///
    /// Returns the ids that went away.
    pub fn remove(&mut self, id: NodeId) -> Vec<NodeId> {
        if self.get(id).is_none() {
            return Vec::new();
        }
        let mut doomed = Vec::new();
        self.collect_subtree(id, &mut doomed);
        let parent = self.get(id).and_then(|n| n.parent);
        match parent {
            Some(p) => {
                if let Some(folder) = self.nodes[p.0 as usize].as_mut().and_then(|n| n.folder_mut())
                {
                    folder.children.retain(|&c| c != id);
                }
            }
            None => self.roots.retain(|&c| c != id),
        }
        for &node in &doomed {
            self.nodes[node.0 as usize] = None;
        }
        self.touch();
        doomed
    }

    fn collect_subtree(&self, id: NodeId, out: &mut Vec<NodeId>) {
        out.push(id);
        if let Some(folder) = self.get(id).and_then(|n| n.folder()) {
            for &child in &folder.children {
                self.collect_subtree(child, out);
            }
        }
    }

    pub fn rename(&mut self, id: NodeId, name: impl Into<String>) -> Result<()> {
        let node = self
            .get_mut(id)
            .ok_or_else(|| Error::out_of_range("that item is no longer in the library"))?;
        node.name = name.into();
        Ok(())
    }

    /// Move `id` into `parent` (or to the top level), at `index` among its new
    /// siblings.
    ///
    /// Refuses to put a folder inside itself, which would detach the subtree
    /// from the tree entirely.
    pub fn move_to(&mut self, id: NodeId, parent: Option<NodeId>, index: usize) -> Result<()> {
        self.check_folder(parent)?;
        if self.get(id).is_none() {
            return Err(Error::out_of_range("that item is no longer in the library"));
        }
        if let Some(target) = parent {
            if target == id || self.is_ancestor(id, target) {
                return Err(Error::format("a folder cannot be moved inside itself"));
            }
        }
        let old_parent = self.get(id).and_then(|n| n.parent);
        match old_parent {
            Some(p) => {
                if let Some(folder) = self.nodes[p.0 as usize].as_mut().and_then(|n| n.folder_mut())
                {
                    folder.children.retain(|&c| c != id);
                }
            }
            None => self.roots.retain(|&c| c != id),
        }
        if let Some(node) = self.nodes[id.0 as usize].as_mut() {
            node.parent = parent;
        }
        let siblings = match parent {
            Some(p) => match self.nodes[p.0 as usize].as_mut().and_then(|n| n.folder_mut()) {
                Some(folder) => &mut folder.children,
                None => return Err(Error::format("that folder went away mid-move")),
            },
            None => &mut self.roots,
        };
        siblings.insert(index.min(siblings.len()), id);
        self.touch();
        Ok(())
    }

    /// Is `ancestor` somewhere above `id`?
    pub fn is_ancestor(&self, ancestor: NodeId, id: NodeId) -> bool {
        let mut at = self.get(id).and_then(|n| n.parent);
        while let Some(current) = at {
            if current == ancestor {
                return true;
            }
            at = self.get(current).and_then(|n| n.parent);
        }
        false
    }

    // ---- loading and saving entries ------------------------------------

    /// Load one entry's sequences.
    pub fn load(&self, id: NodeId) -> Result<Alignment> {
        self.require_entry(id)?.load()
    }

    /// Gather the sequences of every entry under `ids` into one set, ready to
    /// be aligned.
    ///
    /// Names are made unique by appending a counter, because two reads of the
    /// same specimen from different files legitimately share a name and an
    /// aligner with duplicate row names produces output nobody can interpret.
    /// Entries that cannot be read are reported alongside what did load, so one
    /// missing file does not lose the rest of the batch.
    pub fn gather(&self, ids: &[NodeId], name: &str) -> (Alignment, Vec<(NodeId, Error)>) {
        let mut sequences: Vec<Sequence> = Vec::new();
        let mut failed = Vec::new();
        let mut seen: Vec<NodeId> = Vec::new();
        for &id in ids {
            for entry_id in self.entries_under(Some(id)) {
                if seen.contains(&entry_id) {
                    continue; // a folder and one of its entries were both selected
                }
                seen.push(entry_id);
                match self.load(entry_id) {
                    Ok(alignment) => sequences.extend(alignment.sequences),
                    Err(e) => failed.push((entry_id, e)),
                }
            }
        }
        let mut alignment = Alignment::new(name, sequences);
        alignment.deduplicate_ids();
        (alignment, failed)
    }

    /// Where saving `id` would write, and whether the user must be asked.
    pub fn save_target(&self, id: NodeId) -> Result<SaveTarget> {
        Ok(self.require_entry(id)?.save_target())
    }

    /// Write `alignment` back for entry `id`, according to `choice`.
    ///
    /// On [`SaveChoice::NewCopy`] the entry remembers the copy, so this is the
    /// last time the question is asked: later saves land on the copy without
    /// prompting, and one edited sequence leaves one extra file behind however
    /// many times it is edited.
    ///
    /// Returns the path written.
    pub fn save_entry(
        &mut self,
        id: NodeId,
        alignment: &Alignment,
        choice: SaveChoice,
        options: &WriteOptions,
    ) -> Result<PathBuf> {
        let entry = self.require_entry(id)?;
        let target = entry.save_target();
        let (path, format, working) = match choice {
            SaveChoice::Overwrite => {
                if !target.can_overwrite() {
                    return Err(Error::format(match &target {
                        SaveTarget::MustCopy(_, why) => why.explain().to_string(),
                        _ => "this entry cannot be written back in place".to_string(),
                    }));
                }
                // Overwriting a working copy keeps it; overwriting the origin
                // means the user chose the original, so any working copy is
                // now stale and is forgotten.
                match &target {
                    SaveTarget::WorkingCopy(p) => {
                        (p.clone(), entry.effective_format(), entry.working.clone())
                    }
                    _ => (entry.origin.clone(), entry.format, None),
                }
            }
            SaveChoice::NewCopy(path) => {
                let format =
                    Format::from_path(&path).filter(|f| f.can_write()).unwrap_or(Format::Fasta);
                (path.clone(), format, Some(path))
            }
        };

        write_padding_if_needed(alignment, &path, format, options)?;

        let entry = self
            .nodes
            .get_mut(id.0 as usize)
            .and_then(|n| n.as_mut())
            .and_then(|n| n.entry_mut())
            .ok_or_else(|| Error::out_of_range("that entry is no longer in the library"))?;
        entry.working = working;
        // `select` and `reversed` describe how to turn the file's contents into
        // what the entry shows. What was just written *is* what the entry
        // shows, so both are now satisfied and must be cleared — leaving
        // `reversed` set would reverse the saved sequence a second time the
        // next time it was read.
        entry.select = None;
        entry.reversed = false;
        entry.kind = EntryKind::of(format, alignment);
        self.touch();
        Ok(path)
    }

    /// Flip an entry's orientation. Nothing is written; the flag is applied
    /// whenever the entry is read.
    pub fn set_reversed(&mut self, id: NodeId, reversed: bool) -> Result<()> {
        let entry = self
            .get_mut(id)
            .and_then(|n| n.entry_mut())
            .ok_or_else(|| Error::format("only sequences can be reversed, not folders"))?;
        entry.reversed = reversed;
        Ok(())
    }

    fn require_entry(&self, id: NodeId) -> Result<&Entry> {
        match self.get(id) {
            Some(node) => node.entry().ok_or_else(|| {
                Error::format(format!("'{}' is a folder, not a sequence", node.name))
            }),
            None => Err(Error::out_of_range("that entry is no longer in the library")),
        }
    }

    /// Entries whose file has gone missing, for the health check the GUI runs
    /// after loading a library.
    pub fn broken_entries(&self) -> Vec<NodeId> {
        self.entries_under(None)
            .into_iter()
            .filter(|&id| self.entry(id).is_some_and(|e| !e.source_path().exists()))
            .collect()
    }

    // ---- used by the on-disk format ------------------------------------

    pub(crate) fn set_saved_revision(&mut self) {
        self.saved_revision = self.revision;
    }

    /// Rebuild a library from a flat depth-first list, as the reader produces.
    pub(crate) fn push_at_depth(
        &mut self,
        depth: usize,
        node: Node,
        stack: &mut Vec<NodeId>,
    ) -> NodeId {
        stack.truncate(depth);
        let parent = stack.last().copied();
        let mut node = node;
        node.parent = parent;
        let is_folder = node.is_folder();
        let id = self.attach(parent, node);
        if is_folder {
            stack.push(id);
        }
        id
    }
}

/// Write `alignment` to `path`, padding ragged rows if the format insists on a
/// rectangle. The in-memory data is left alone, matching what the editor does.
fn write_padding_if_needed(
    alignment: &Alignment,
    path: &Path,
    format: Format,
    options: &WriteOptions,
) -> Result<()> {
    match tolviewer_io::write_file(alignment, path, format, options) {
        Ok(()) => Ok(()),
        Err(Error::Format(msg)) if !alignment.is_aligned() => {
            let mut padded = alignment.clone();
            padded.pad_to_width();
            tolviewer_io::write_file(&padded, path, format, options).map_err(|_| Error::Format(msg))
        }
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("tolviewer-library-{name}"));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn fasta(dir: &Path, name: &str, rows: &[(&str, &str)]) -> PathBuf {
        let path = dir.join(name);
        let text: String = rows.iter().map(|(id, seq)| format!(">{id}\n{seq}\n")).collect();
        std::fs::write(&path, text).unwrap();
        path
    }

    /// The example from the request: a project with two loci under it.
    fn project(dir: &Path) -> (Library, NodeId, NodeId) {
        let mut lib = Library::new("Lace bug project");
        let project = lib.add_folder(None, "Lace bug project").unwrap();
        let ssu = lib.add_folder(Some(project), "18S").unwrap();
        let lsu = lib.add_folder(Some(project), "28S").unwrap();
        let f = fasta(dir, "TL-2213_18S.fasta", &[("TL-2213_18S", "ACGTACGT")]);
        lib.add_file(Some(ssu), &f).unwrap();
        (lib, ssu, lsu)
    }

    #[test]
    fn folders_nest_and_report_their_path() {
        let dir = scratch("nest");
        let (lib, ssu, _) = project(&dir);
        let entry = lib.entries_under(Some(ssu))[0];
        assert_eq!(lib.path_of(entry), "Lace bug project / 18S / TL-2213_18S");
        assert_eq!(lib.walk().len(), 4);
        assert_eq!(lib.walk()[0].1, 0);
        assert_eq!(lib.walk()[1].1, 1);
    }

    #[test]
    fn nothing_can_be_filed_inside_a_sequence() {
        let dir = scratch("inside");
        let (mut lib, ssu, _) = project(&dir);
        let entry = lib.entries_under(Some(ssu))[0];
        let e = lib.add_folder(Some(entry), "nope").unwrap_err();
        assert!(e.to_string().contains("not a folder"), "{e}");
        let f = fasta(&dir, "x.fasta", &[("a", "ACGT")]);
        assert!(lib.add_file(Some(entry), &f).is_err());
    }

    #[test]
    fn adding_a_file_does_not_touch_it() {
        let dir = scratch("insitu");
        let path = fasta(&dir, "reads.fasta", &[("a", "ACGT"), ("b", "ACGA")]);
        let before = std::fs::read(&path).unwrap();
        let mut lib = Library::new("l");
        let id = lib.add_file(None, &path).unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), before);
        assert_eq!(lib.entry(id).unwrap().origin, path);
        assert_eq!(lib.entry(id).unwrap().kind, EntryKind::Alignment);
        assert_eq!(lib.load(id).unwrap().len(), 2);
    }

    #[test]
    fn removing_an_entry_leaves_its_file_alone() {
        let dir = scratch("remove");
        let (mut lib, ssu, _) = project(&dir);
        let entry = lib.entries_under(Some(ssu))[0];
        let path = lib.entry(entry).unwrap().origin.clone();
        let gone = lib.remove(ssu);
        assert_eq!(gone.len(), 2, "the folder and its entry");
        assert!(lib.get(entry).is_none());
        assert!(path.exists(), "removing from a library must not delete data");
    }

    #[test]
    fn a_stale_id_resolves_to_nothing_rather_than_to_the_wrong_node() {
        let dir = scratch("stale");
        let (mut lib, ssu, lsu) = project(&dir);
        let entry = lib.entries_under(Some(ssu))[0];
        lib.remove(entry);
        assert!(lib.get(entry).is_none());
        // A new node must not inherit the dead id.
        let fresh = lib.add_folder(Some(lsu), "new").unwrap();
        assert_ne!(fresh, entry);
        assert!(lib.load(entry).is_err());
    }

    #[test]
    fn reversing_is_a_flag_not_a_rewrite() {
        let dir = scratch("reverse");
        let path = fasta(&dir, "read.fasta", &[("a", "AAAACCCC")]);
        let before = std::fs::read(&path).unwrap();
        let mut lib = Library::new("l");
        let id = lib.add_file(None, &path).unwrap();
        lib.set_reversed(id, true).unwrap();
        assert_eq!(lib.load(id).unwrap().sequences[0].residues, b"GGGGTTTT");
        assert_eq!(std::fs::read(&path).unwrap(), before, "the file must be untouched");
        lib.set_reversed(id, false).unwrap();
        assert_eq!(lib.load(id).unwrap().sequences[0].residues, b"AAAACCCC");
    }

    #[test]
    fn folders_cannot_be_reversed() {
        let dir = scratch("revfolder");
        let (mut lib, ssu, _) = project(&dir);
        assert!(lib.set_reversed(ssu, true).is_err());
    }

    #[test]
    fn saving_over_the_original_is_flagged_but_allowed() {
        let dir = scratch("overwrite");
        let path = fasta(&dir, "reads.fasta", &[("a", "ACGT")]);
        let mut lib = Library::new("l");
        let id = lib.add_file(None, &path).unwrap();

        let target = lib.save_target(id).unwrap();
        assert_eq!(target, SaveTarget::Original(path.clone()));
        assert!(target.needs_confirmation());
        assert!(target.can_overwrite());

        let mut edited = lib.load(id).unwrap();
        edited.sequences[0].residues = b"TTTT".to_vec();
        let written =
            lib.save_entry(id, &edited, SaveChoice::Overwrite, &WriteOptions::default()).unwrap();
        assert_eq!(written, path);
        assert!(std::fs::read_to_string(&path).unwrap().contains("TTTT"));
        // Still the original, so the next save asks again.
        assert!(lib.save_target(id).unwrap().needs_confirmation());
    }

    #[test]
    fn saving_as_a_copy_is_asked_once_and_then_never_again() {
        let dir = scratch("copy");
        let path = fasta(&dir, "reads.fasta", &[("a", "ACGT")]);
        let original = std::fs::read(&path).unwrap();
        let mut lib = Library::new("l");
        let id = lib.add_file(None, &path).unwrap();

        let copy = lib.entry(id).unwrap().suggested_copy();
        assert_eq!(copy.file_name().unwrap(), "reads.edited.fasta");

        let mut edited = lib.load(id).unwrap();
        edited.sequences[0].residues = b"TTTT".to_vec();
        let written = lib
            .save_entry(id, &edited, SaveChoice::NewCopy(copy.clone()), &WriteOptions::default())
            .unwrap();
        assert_eq!(written, copy);
        assert_eq!(std::fs::read(&path).unwrap(), original, "the original is untouched");

        // From now on it is a working copy: no confirmation, same file.
        let target = lib.save_target(id).unwrap();
        assert_eq!(target, SaveTarget::WorkingCopy(copy.clone()));
        assert!(!target.needs_confirmation());
        assert_eq!(lib.load(id).unwrap().sequences[0].residues, b"TTTT");

        edited.sequences[0].residues = b"GGGG".to_vec();
        let again =
            lib.save_entry(id, &edited, SaveChoice::Overwrite, &WriteOptions::default()).unwrap();
        assert_eq!(again, copy, "further edits must land on the same copy, not pile up");
        assert_eq!(std::fs::read_dir(&dir).unwrap().count(), 2);
        assert_eq!(std::fs::read(&path).unwrap(), original);
    }

    #[test]
    fn a_suggested_copy_never_lands_on_an_existing_file() {
        let dir = scratch("unique");
        let path = fasta(&dir, "reads.fasta", &[("a", "ACGT")]);
        fasta(&dir, "reads.edited.fasta", &[("a", "ACGT")]);
        let mut lib = Library::new("l");
        let id = lib.add_file(None, &path).unwrap();
        assert_eq!(
            lib.entry(id).unwrap().suggested_copy().file_name().unwrap(),
            "reads.edited-2.fasta"
        );
    }

    #[test]
    fn extracting_a_row_cannot_overwrite_the_alignment_it_came_from() {
        let dir = scratch("extract");
        let path = fasta(&dir, "msa.fasta", &[("a", "AC-GT"), ("b", "ACGGT"), ("c", "AC-GA")]);
        let mut lib = Library::new("l");
        let msa = lib.add_file(None, &path).unwrap();
        let one = lib.add_selection(None, msa, vec!["b".to_string()], "b (extracted)").unwrap();

        let loaded = lib.load(one).unwrap();
        assert_eq!(loaded.len(), 1);
        assert_eq!(loaded.sequences[0].id, "b");
        assert_eq!(loaded.sequences[0].residues, b"ACGGT");

        let target = lib.save_target(one).unwrap();
        assert!(!target.can_overwrite(), "writing one row back would destroy a and c");
        assert_eq!(
            target,
            SaveTarget::MustCopy(lib.entry(one).unwrap().suggested_copy(), CopyReason::PartOfAFile)
        );
        let e = lib
            .save_entry(one, &loaded, SaveChoice::Overwrite, &WriteOptions::default())
            .unwrap_err();
        assert!(e.to_string().contains("part of a larger file"), "{e}");
        // The alignment is still intact.
        assert_eq!(lib.load(msa).unwrap().len(), 3);
    }

    #[test]
    fn an_extract_saved_to_a_copy_becomes_a_file_in_its_own_right() {
        let dir = scratch("extract-save");
        let path = fasta(&dir, "msa.fasta", &[("a", "ACGT"), ("b", "ACGA")]);
        let mut lib = Library::new("l");
        let msa = lib.add_file(None, &path).unwrap();
        let one = lib.add_selection(None, msa, vec!["b".to_string()], "b").unwrap();
        let loaded = lib.load(one).unwrap();
        let copy = dir.join("b.fasta");
        lib.save_entry(one, &loaded, SaveChoice::NewCopy(copy.clone()), &WriteOptions::default())
            .unwrap();
        let entry = lib.entry(one).unwrap();
        assert!(entry.select.is_none(), "the copy is the whole file now");
        assert_eq!(entry.save_target(), SaveTarget::WorkingCopy(copy));
        assert_eq!(lib.load(one).unwrap().len(), 1);
    }

    #[test]
    fn a_reversed_entry_saved_to_a_copy_is_saved_the_way_it_reads() {
        let dir = scratch("reverse-save");
        let path = fasta(&dir, "read.fasta", &[("a", "AAAACCCC")]);
        let mut lib = Library::new("l");
        let id = lib.add_file(None, &path).unwrap();
        lib.set_reversed(id, true).unwrap();
        let shown = lib.load(id).unwrap();
        assert_eq!(shown.sequences[0].residues, b"GGGGTTTT");

        let copy = dir.join("read.rc.fasta");
        lib.save_entry(id, &shown, SaveChoice::NewCopy(copy), &WriteOptions::default()).unwrap();
        assert!(!lib.entry(id).unwrap().reversed, "the copy is already reversed");
        assert_eq!(
            lib.load(id).unwrap().sequences[0].residues,
            b"GGGGTTTT",
            "reading it back must not flip it a second time"
        );
    }

    #[test]
    fn overwriting_a_reversed_entry_does_not_reverse_it_again() {
        let dir = scratch("reverse-overwrite");
        let path = fasta(&dir, "read.fasta", &[("a", "AAAACCCC")]);
        let mut lib = Library::new("l");
        let id = lib.add_file(None, &path).unwrap();
        lib.set_reversed(id, true).unwrap();

        // What the user sees, and therefore what saving must produce.
        let shown = lib.load(id).unwrap();
        assert_eq!(shown.sequences[0].residues, b"GGGGTTTT");
        lib.save_entry(id, &shown, SaveChoice::Overwrite, &WriteOptions::default()).unwrap();

        assert!(std::fs::read_to_string(&path).unwrap().contains("GGGGTTTT"));
        assert!(!lib.entry(id).unwrap().reversed, "the file is already the right way round");
        assert_eq!(
            lib.load(id).unwrap().sequences[0].residues,
            b"GGGGTTTT",
            "reading it back must not flip it a second time"
        );
    }

    #[test]
    fn a_read_only_format_can_only_be_saved_to_a_copy() {
        let dir = scratch("readonly");
        let ab1 =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../../testdata/ab1/tingidae_COI_F.ab1");
        let path = dir.join("read.ab1");
        std::fs::copy(&ab1, &path).unwrap();
        let mut lib = Library::new("l");
        let id = lib.add_file(None, &path).unwrap();
        assert_eq!(lib.entry(id).unwrap().kind, EntryKind::Trace);
        assert_eq!(lib.get(id).unwrap().name, "TL-2213_COI_F", "traces are named for the sample");

        let target = lib.save_target(id).unwrap();
        assert!(!target.can_overwrite());
        assert_eq!(target.path().extension().unwrap(), "fasta", "an ab1 cannot be written");
        match target {
            SaveTarget::MustCopy(_, why) => assert_eq!(why, CopyReason::ReadOnlyFormat),
            other => panic!("wrong target: {other:?}"),
        }
        // The chromatogram is available, and follows the entry's orientation.
        let forward = lib.entry(id).unwrap().load_trace().unwrap();
        lib.set_reversed(id, true).unwrap();
        let mut reversed = lib.entry(id).unwrap().load_trace().unwrap();
        reversed.reverse_complement();
        assert_eq!(reversed.calls, forward.calls);
    }

    #[test]
    fn gathering_a_folder_collects_every_read_under_it() {
        let dir = scratch("gather");
        let mut lib = Library::new("l");
        let project = lib.add_folder(None, "project").unwrap();
        let ssu = lib.add_folder(Some(project), "18S").unwrap();
        lib.add_file(Some(ssu), &fasta(&dir, "a.fasta", &[("TL-1", "ACGT")])).unwrap();
        lib.add_file(Some(ssu), &fasta(&dir, "b.fasta", &[("TL-2", "ACGA")])).unwrap();
        let lsu = lib.add_folder(Some(project), "28S").unwrap();
        lib.add_file(Some(lsu), &fasta(&dir, "c.fasta", &[("TL-1", "TTTT")])).unwrap();

        let (all, failed) = lib.gather(&[project], "everything");
        assert!(failed.is_empty());
        assert_eq!(all.len(), 3);
        // Two files hold a sequence called TL-1; the names must not collide.
        let ids: Vec<&str> = all.sequences.iter().map(|s| s.id.as_str()).collect();
        assert_eq!(ids.len(), 3);
        let mut unique = ids.clone();
        unique.sort_unstable();
        unique.dedup();
        assert_eq!(unique.len(), 3, "duplicate names survived: {ids:?}");

        let (one_locus, _) = lib.gather(&[ssu], "18S");
        assert_eq!(one_locus.len(), 2);
    }

    #[test]
    fn selecting_a_folder_and_a_file_inside_it_does_not_double_up() {
        let dir = scratch("double");
        let mut lib = Library::new("l");
        let folder = lib.add_folder(None, "f").unwrap();
        let entry = lib.add_file(Some(folder), &fasta(&dir, "a.fasta", &[("a", "ACGT")])).unwrap();
        let (all, _) = lib.gather(&[folder, entry], "x");
        assert_eq!(all.len(), 1);
    }

    #[test]
    fn a_missing_file_is_reported_without_losing_the_rest_of_the_batch() {
        let dir = scratch("missing");
        let mut lib = Library::new("l");
        let good = lib.add_file(None, &fasta(&dir, "a.fasta", &[("a", "ACGT")])).unwrap();
        let doomed = fasta(&dir, "b.fasta", &[("b", "ACGA")]);
        let bad = lib.add_file(None, &doomed).unwrap();
        std::fs::remove_file(&doomed).unwrap();

        assert_eq!(lib.broken_entries(), vec![bad]);
        let (all, failed) = lib.gather(&[good, bad], "x");
        assert_eq!(all.len(), 1, "the readable file must still come through");
        assert_eq!(failed.len(), 1);
        assert_eq!(failed[0].0, bad);
    }

    #[test]
    fn moving_a_folder_into_itself_is_refused() {
        let dir = scratch("move");
        let (mut lib, ssu, lsu) = project(&dir);
        let project = lib.get(ssu).unwrap().parent.unwrap();
        assert!(lib.move_to(project, Some(ssu), 0).is_err());
        assert!(lib.move_to(ssu, Some(ssu), 0).is_err());
        // A legal move keeps the tree walkable.
        lib.move_to(ssu, Some(lsu), 0).unwrap();
        assert_eq!(lib.get(ssu).unwrap().parent, Some(lsu));
        assert!(lib.children(Some(project)).contains(&lsu));
        assert!(!lib.children(Some(project)).contains(&ssu));
        assert_eq!(lib.walk().len(), 4);
    }

    #[test]
    fn moving_to_the_top_level_works_and_keeps_order() {
        let dir = scratch("move-root");
        let (mut lib, ssu, _) = project(&dir);
        lib.move_to(ssu, None, 0).unwrap();
        assert_eq!(lib.roots()[0], ssu);
        assert_eq!(lib.get(ssu).unwrap().parent, None);
        assert_eq!(lib.roots().len(), 2);
    }

    #[test]
    fn editing_the_primer_list_makes_the_library_worth_saving() {
        let mut lib = Library::new("l");
        lib.mark_saved();
        lib.primers.push(crate::primer::Primer::new("p", "ACGT").unwrap());
        lib.touch();
        assert!(lib.is_dirty());
    }

    #[test]
    fn the_dirty_flag_follows_edits_and_saves() {
        let dir = scratch("dirty");
        let mut lib = Library::new("l");
        lib.mark_saved();
        assert!(!lib.is_dirty());
        lib.add_file(None, &fasta(&dir, "a.fasta", &[("a", "ACGT")])).unwrap();
        assert!(lib.is_dirty());
        lib.mark_saved();
        assert!(!lib.is_dirty());
        lib.add_folder(None, "f").unwrap();
        assert!(lib.is_dirty());
    }
}
