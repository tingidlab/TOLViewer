//! The library panel's state: what is selected, what is open, and the
//! questions waiting to be answered.
//!
//! The library itself lives in `tolviewer-library` and knows nothing about the
//! GUI. This is the layer between: it holds the multi-selection the tree
//! widget maintains, the pending confirmations the in-situ save policy raises,
//! and the settings the library-wide commands run with.

use std::path::PathBuf;

use tolviewer_library::{ConcatOptions, Library, NodeId, SaveTarget, TrimOptions};

/// The library and everything the panel needs to draw and drive it.
pub struct LibraryState {
    pub library: Library,
    /// Every selected node, in the order it was clicked. Folders and entries
    /// can be selected together; commands resolve a folder to the entries under
    /// it.
    pub selected: Vec<NodeId>,
    /// The node the next shift-click extends from.
    pub anchor: Option<NodeId>,
    /// A node being renamed inline, and the text so far.
    pub renaming: Option<(NodeId, String)>,
    /// Whether the panel is shown at all.
    pub visible: bool,
    pub concat: ConcatOptions,
    pub trim: TrimOptions,
    /// Documents that came from a library entry, so closing one can be matched
    /// back up. Keyed by entry.
    pub open: Vec<(NodeId, usize)>,
}

impl Default for LibraryState {
    fn default() -> Self {
        LibraryState {
            library: Library::default(),
            selected: Vec::new(),
            anchor: None,
            renaming: None,
            visible: true,
            concat: ConcatOptions::default(),
            trim: TrimOptions::default(),
            open: Vec::new(),
        }
    }
}

impl LibraryState {
    /// The single selected node, when exactly one is selected.
    pub fn only_selected(&self) -> Option<NodeId> {
        match self.selected.as_slice() {
            [one] => Some(*one),
            _ => None,
        }
    }

    /// Where a newly created folder or file should go: inside the selected
    /// folder, or beside the selected entry, or at the top level.
    pub fn insertion_parent(&self) -> Option<NodeId> {
        let id = self.only_selected()?;
        let node = self.library.get(id)?;
        if node.is_folder() {
            Some(id)
        } else {
            node.parent
        }
    }

    pub fn is_selected(&self, id: NodeId) -> bool {
        self.selected.contains(&id)
    }

    /// Replace the selection with one node.
    pub fn select_only(&mut self, id: NodeId) {
        self.selected.clear();
        self.selected.push(id);
        self.anchor = Some(id);
    }

    /// Add or remove one node, as ctrl-click does.
    pub fn toggle(&mut self, id: NodeId) {
        match self.selected.iter().position(|&s| s == id) {
            Some(at) => {
                self.selected.remove(at);
            }
            None => self.selected.push(id),
        }
        self.anchor = Some(id);
    }

    /// Select everything visible between the anchor and `id`, as shift-click
    /// does. Ranges are taken over the tree as drawn, so what is selected is
    /// what the user saw between the two clicks.
    pub fn extend_to(&mut self, id: NodeId) {
        let Some(anchor) = self.anchor else {
            self.select_only(id);
            return;
        };
        let order: Vec<NodeId> = self.visible_order();
        let (Some(from), Some(to)) =
            (order.iter().position(|&n| n == anchor), order.iter().position(|&n| n == id))
        else {
            self.select_only(id);
            return;
        };
        let (lo, hi) = if from <= to { (from, to) } else { (to, from) };
        self.selected = order[lo..=hi].to_vec();
    }

    /// The nodes the tree currently draws, top to bottom: collapsed folders
    /// hide their children.
    pub fn visible_order(&self) -> Vec<NodeId> {
        let mut out = Vec::new();
        for &root in self.library.roots() {
            self.visible_into(root, &mut out);
        }
        out
    }

    fn visible_into(&self, id: NodeId, out: &mut Vec<NodeId>) {
        out.push(id);
        if let Some(folder) = self.library.get(id).and_then(|n| n.folder()) {
            if folder.expanded {
                for &child in &folder.children {
                    self.visible_into(child, out);
                }
            }
        }
    }

    /// Drop selected ids that no longer exist, after a removal.
    pub fn prune(&mut self) {
        self.selected.retain(|&id| self.library.get(id).is_some());
        if self.anchor.is_some_and(|a| self.library.get(a).is_none()) {
            self.anchor = self.selected.last().copied();
        }
        self.open.retain(|(id, _)| self.library.get(*id).is_some());
    }

    /// Every entry the current selection resolves to, folders expanded, in tree
    /// order and without duplicates.
    pub fn selected_entries(&self) -> Vec<NodeId> {
        let mut out: Vec<NodeId> = Vec::new();
        for &id in &self.selected {
            for entry in self.library.entries_under(Some(id)) {
                if !out.contains(&entry) {
                    out.push(entry);
                }
            }
        }
        out
    }

    /// The document index showing `entry`, if one is open.
    pub fn document_for(&self, entry: NodeId) -> Option<usize> {
        self.open.iter().find(|(id, _)| *id == entry).map(|(_, doc)| *doc)
    }

    pub fn note_open(&mut self, entry: NodeId, doc: usize) {
        self.open.retain(|(id, _)| *id != entry);
        self.open.push((entry, doc));
    }

    /// Keep the document indices right after one is closed.
    pub fn document_closed(&mut self, closed: usize) {
        self.open.retain(|(_, doc)| *doc != closed);
        for (_, doc) in &mut self.open {
            if *doc > closed {
                *doc -= 1;
            }
        }
    }
}

/// A save the user has to answer for: writing here would replace data the
/// library did not create.
pub struct PendingSave {
    /// The document being saved.
    pub doc: usize,
    pub entry: NodeId,
    pub target: SaveTarget,
    /// The path a copy would go to, editable in the dialog.
    pub copy_to: PathBuf,
    /// What the entry is called, for the wording of the question.
    pub name: String,
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("tolviewer-app-lib-{name}"));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn fasta(dir: &Path, name: &str) -> PathBuf {
        let path = dir.join(name);
        std::fs::write(&path, ">a\nACGT\n").unwrap();
        path
    }

    /// A project with two loci and one read in each.
    fn state(dir: &Path) -> (LibraryState, Vec<NodeId>) {
        let mut s = LibraryState::default();
        let project = s.library.add_folder(None, "Lace bug project").unwrap();
        let ssu = s.library.add_folder(Some(project), "18S").unwrap();
        let a = s.library.add_file(Some(ssu), &fasta(dir, "a.fasta")).unwrap();
        let lsu = s.library.add_folder(Some(project), "28S").unwrap();
        let b = s.library.add_file(Some(lsu), &fasta(dir, "b.fasta")).unwrap();
        (s, vec![project, ssu, a, lsu, b])
    }

    #[test]
    fn clicking_replaces_the_selection_and_ctrl_click_adds_to_it() {
        let dir = scratch("select");
        let (mut s, ids) = state(&dir);
        s.select_only(ids[2]);
        assert_eq!(s.selected, vec![ids[2]]);
        s.toggle(ids[4]);
        assert_eq!(s.selected, vec![ids[2], ids[4]]);
        s.toggle(ids[2]);
        assert_eq!(s.selected, vec![ids[4]], "ctrl-clicking again deselects");
        s.select_only(ids[0]);
        assert_eq!(s.selected, vec![ids[0]]);
    }

    #[test]
    fn shift_click_takes_the_run_the_user_can_see() {
        let dir = scratch("shift");
        let (mut s, ids) = state(&dir);
        s.select_only(ids[1]); // 18S
        s.extend_to(ids[3]); // 28S
        assert_eq!(s.selected, vec![ids[1], ids[2], ids[3]]);

        // With 18S collapsed, its read is not on screen and not in the range.
        s.library.get_mut(ids[1]).unwrap().folder_mut().unwrap().expanded = false;
        s.select_only(ids[1]);
        s.extend_to(ids[3]);
        assert_eq!(s.selected, vec![ids[1], ids[3]]);
    }

    #[test]
    fn shift_click_works_in_both_directions() {
        let dir = scratch("shift-back");
        let (mut s, ids) = state(&dir);
        s.select_only(ids[4]);
        s.extend_to(ids[1]);
        assert_eq!(s.selected, vec![ids[1], ids[2], ids[3], ids[4]]);
    }

    #[test]
    fn selecting_a_folder_means_every_read_under_it() {
        let dir = scratch("entries");
        let (mut s, ids) = state(&dir);
        s.select_only(ids[0]);
        assert_eq!(s.selected_entries(), vec![ids[2], ids[4]]);
        // A folder and a read inside it must not yield the read twice.
        s.toggle(ids[2]);
        assert_eq!(s.selected_entries(), vec![ids[2], ids[4]]);
    }

    #[test]
    fn a_new_folder_goes_inside_a_selected_folder_and_beside_a_selected_read() {
        let dir = scratch("insertion");
        let (mut s, ids) = state(&dir);
        s.select_only(ids[1]);
        assert_eq!(s.insertion_parent(), Some(ids[1]));
        s.select_only(ids[2]);
        assert_eq!(s.insertion_parent(), Some(ids[1]), "beside the read, not inside it");
        s.selected.clear();
        assert_eq!(s.insertion_parent(), None);
        s.select_only(ids[1]);
        s.toggle(ids[3]);
        assert_eq!(s.insertion_parent(), None, "an ambiguous selection picks nothing");
    }

    #[test]
    fn removing_a_folder_prunes_the_selection_and_the_open_list() {
        let dir = scratch("prune");
        let (mut s, ids) = state(&dir);
        s.select_only(ids[2]);
        s.toggle(ids[4]);
        s.note_open(ids[2], 0);
        s.note_open(ids[4], 1);
        s.library.remove(ids[1]);
        s.prune();
        assert_eq!(s.selected, vec![ids[4]]);
        assert_eq!(s.document_for(ids[2]), None);
        assert_eq!(s.document_for(ids[4]), Some(1));
    }

    #[test]
    fn closing_a_document_renumbers_the_ones_after_it() {
        let dir = scratch("close");
        let (mut s, ids) = state(&dir);
        s.note_open(ids[2], 0);
        s.note_open(ids[4], 1);
        s.document_closed(0);
        assert_eq!(s.document_for(ids[2]), None);
        assert_eq!(s.document_for(ids[4]), Some(0));
    }

    #[test]
    fn reopening_an_entry_replaces_its_old_document() {
        let dir = scratch("reopen");
        let (mut s, ids) = state(&dir);
        s.note_open(ids[2], 0);
        s.note_open(ids[2], 3);
        assert_eq!(s.open.len(), 1);
        assert_eq!(s.document_for(ids[2]), Some(3));
    }

    #[test]
    fn a_collapsed_tree_still_orders_what_it_shows() {
        let dir = scratch("order");
        let (mut s, ids) = state(&dir);
        assert_eq!(s.visible_order(), ids);
        s.library.get_mut(ids[0]).unwrap().folder_mut().unwrap().expanded = false;
        assert_eq!(s.visible_order(), vec![ids[0]]);
    }
}
