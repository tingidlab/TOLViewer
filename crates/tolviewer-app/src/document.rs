//! One open file: the alignment, its edit history, view state and the derived
//! data the canvas needs.

use std::path::PathBuf;

use tolviewer_core::{Alignment, Alphabet, Consensus, EditOp, Result, UndoStack};
use tolviewer_io::Format;

use crate::selection::{Cell, Selection, SelectionMode};

/// Derived data that must be recomputed when the alignment changes. Keyed on
/// the undo stack's revision so a stale cache is impossible to use by accident.
#[derive(Default)]
struct Derived {
    revision: Option<u64>,
    consensus: Consensus,
}

pub struct Document {
    pub alignment: Alignment,
    pub undo: UndoStack,
    /// Where it came from, and where Save writes back to.
    pub path: Option<PathBuf>,
    /// Format it was read as, used as the default for Save.
    pub format: Format,
    pub selection: Selection,
    /// First visible column and row, driven by the scroll area.
    pub scroll_col: f32,
    pub scroll_row: f32,
    /// Cell width in points; height is derived from it.
    pub zoom: f32,
    /// Columns kept by the last cleaning run, for the overlay track.
    pub clean_mask: Option<Vec<bool>>,
    /// Revision the clean mask was computed at; the overlay is hidden once the
    /// alignment moves on.
    clean_mask_revision: u64,
    /// Revision at the last successful save.
    saved_revision: u64,
    derived: Derived,
    alphabet: Alphabet,
}

impl Document {
    pub fn new(mut alignment: Alignment, path: Option<PathBuf>, format: Format) -> Self {
        let alphabet = alignment.alphabet();
        Document {
            alignment,
            undo: UndoStack::new(),
            path,
            format,
            selection: Selection::default(),
            scroll_col: 0.0,
            scroll_row: 0.0,
            zoom: 12.0,
            clean_mask: None,
            clean_mask_revision: 0,
            saved_revision: 0,
            derived: Derived::default(),
            alphabet,
        }
    }

    /// Tab label: file name, or the alignment name for unsaved documents.
    pub fn title(&self) -> String {
        let base = self
            .path
            .as_ref()
            .and_then(|p| p.file_name())
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| {
                if self.alignment.name.is_empty() {
                    "untitled".to_string()
                } else {
                    self.alignment.name.clone()
                }
            });
        if self.is_dirty() {
            format!("{base} *")
        } else {
            base
        }
    }

    pub fn is_dirty(&self) -> bool {
        self.undo.revision() != self.saved_revision
    }

    pub fn mark_saved(&mut self) {
        self.saved_revision = self.undo.revision();
    }

    pub fn alphabet(&self) -> Alphabet {
        self.alphabet
    }

    pub fn set_alphabet(&mut self, alphabet: Alphabet) {
        self.alphabet = alphabet;
        self.alignment.set_alphabet(alphabet);
        self.derived.revision = None;
    }

    pub fn rows(&self) -> usize {
        self.alignment.len()
    }

    pub fn width(&self) -> usize {
        self.alignment.width()
    }

    /// Consensus and per-column statistics, recomputed only when the alignment
    /// has changed since the last call.
    pub fn consensus(&mut self) -> &Consensus {
        let rev = self.undo.revision();
        if self.derived.revision != Some(rev) {
            self.derived.consensus =
                Consensus::compute(&self.alignment, self.alphabet, 0.5, 0.05);
            self.derived.revision = Some(rev);
        }
        &self.derived.consensus
    }

    /// The cleaning mask, but only while it still matches the alignment.
    pub fn live_clean_mask(&self) -> Option<&[bool]> {
        let mask = self.clean_mask.as_deref()?;
        (self.clean_mask_revision == self.undo.revision() && mask.len() == self.width())
            .then_some(mask)
    }

    pub fn set_clean_mask(&mut self, mask: Vec<bool>) {
        self.clean_mask_revision = self.undo.revision();
        self.clean_mask = Some(mask);
    }

    /// Apply an edit through the undo stack. The alphabet is left alone: a
    /// user typing an unexpected letter should see it, not have the whole
    /// document reinterpreted.
    pub fn apply(&mut self, op: EditOp) -> Result<()> {
        self.undo.apply(&mut self.alignment, op)?;
        self.clamp_selection();
        Ok(())
    }

    /// Replace the whole alignment as one undoable step (align, clean, sort).
    pub fn replace(&mut self, label: &str, alignment: Alignment) -> Result<()> {
        self.apply(EditOp::Replace { label: label.to_string(), alignment: Box::new(alignment) })
    }

    pub fn undo(&mut self) -> Result<Option<String>> {
        let r = self.undo.undo(&mut self.alignment)?;
        self.clamp_selection();
        Ok(r)
    }

    pub fn redo(&mut self) -> Result<Option<String>> {
        let r = self.undo.redo(&mut self.alignment)?;
        self.clamp_selection();
        Ok(r)
    }

    /// Keep the caret inside the grid after rows or columns disappear.
    pub fn clamp_selection(&mut self) {
        let rows = self.rows();
        let cols = self.width();
        if rows == 0 || cols == 0 {
            self.selection = Selection::default();
            return;
        }
        for cell in [&mut self.selection.anchor, &mut self.selection.cursor] {
            cell.row = cell.row.min(rows - 1);
            cell.col = cell.col.min(cols - 1);
        }
    }

    /// The rows the user is acting on: the selected rows, or every row when
    /// nothing is selected.
    pub fn target_rows(&self) -> Vec<usize> {
        if self.selection.active {
            self.selection.rows(self.rows()).collect()
        } else {
            (0..self.rows()).collect()
        }
    }

    /// The columns the user is acting on, or the full width when nothing is
    /// selected.
    pub fn target_cols(&self) -> std::ops::Range<usize> {
        if self.selection.active {
            self.selection.cols(self.width())
        } else {
            0..self.width()
        }
    }

    /// Select everything.
    pub fn select_all(&mut self) {
        let (rows, cols) = (self.rows(), self.width());
        if rows == 0 || cols == 0 {
            return;
        }
        self.selection.place(Cell::new(0, 0), SelectionMode::Cells);
        self.selection.extend_to(Cell::new(rows - 1, cols - 1));
    }

    /// Ungapped residue number at the caret, for the status bar. `None` when
    /// the caret sits on a gap.
    pub fn caret_residue_number(&self) -> Option<usize> {
        let seq = self.alignment.sequences.get(self.selection.cursor.row)?;
        seq.residue_index_at(self.selection.cursor.col).map(|i| i + 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tolviewer_core::Sequence;

    fn doc() -> Document {
        Document::new(
            Alignment::new(
                "t",
                vec![
                    Sequence::new("a", *b"ACGT"),
                    Sequence::new("b", *b"ACGT"),
                    Sequence::new("c", *b"ATGT"),
                ],
            ),
            None,
            Format::Fasta,
        )
    }

    #[test]
    fn starts_clean_and_becomes_dirty_on_edit() {
        let mut d = doc();
        assert!(!d.is_dirty());
        d.apply(EditOp::SetResidue { row: 0, col: 0, residue: b'G' }).unwrap();
        assert!(d.is_dirty());
        d.mark_saved();
        assert!(!d.is_dirty());
    }

    #[test]
    fn undo_returns_to_clean() {
        let mut d = doc();
        d.apply(EditOp::SetResidue { row: 0, col: 0, residue: b'G' }).unwrap();
        d.undo().unwrap();
        // The revision advances on undo too, so the document stays dirty; that
        // is deliberate, matching editors that cannot prove the file matches.
        assert!(d.is_dirty());
        assert_eq!(d.alignment.sequences[0].residues, b"ACGT");
    }

    #[test]
    fn consensus_is_cached_until_an_edit() {
        let mut d = doc();
        assert_eq!(d.consensus().residues, b"ACGT");
        d.apply(EditOp::SetResidue { row: 0, col: 1, residue: b'T' }).unwrap();
        d.apply(EditOp::SetResidue { row: 1, col: 1, residue: b'T' }).unwrap();
        assert_eq!(d.consensus().residues, b"ATGT");
    }

    #[test]
    fn clean_mask_goes_stale_after_an_edit() {
        let mut d = doc();
        d.set_clean_mask(vec![true, true, false, true]);
        assert!(d.live_clean_mask().is_some());
        d.apply(EditOp::SetResidue { row: 0, col: 0, residue: b'G' }).unwrap();
        assert!(d.live_clean_mask().is_none());
    }

    #[test]
    fn selection_is_clamped_when_rows_go_away() {
        let mut d = doc();
        d.select_all();
        d.apply(EditOp::RemoveSequence { row: 2 }).unwrap();
        assert!(d.selection.cursor.row < d.rows());
    }

    #[test]
    fn target_rows_defaults_to_everything() {
        let mut d = doc();
        assert_eq!(d.target_rows(), vec![0, 1, 2]);
        d.selection.place(Cell::new(1, 0), SelectionMode::Rows);
        d.selection.extend_to(Cell::new(2, 0));
        assert_eq!(d.target_rows(), vec![1, 2]);
    }

    #[test]
    fn caret_residue_number_is_one_based_and_skips_gaps() {
        let mut d = Document::new(
            Alignment::new("t", vec![Sequence::new("a", *b"A-CG")]),
            None,
            Format::Fasta,
        );
        d.selection.place(Cell::new(0, 2), SelectionMode::Cells);
        assert_eq!(d.caret_residue_number(), Some(2));
        d.selection.place(Cell::new(0, 1), SelectionMode::Cells);
        assert_eq!(d.caret_residue_number(), None);
    }

    #[test]
    fn title_marks_unsaved_changes() {
        let mut d = doc();
        assert_eq!(d.title(), "t");
        d.apply(EditOp::SetResidue { row: 0, col: 0, residue: b'G' }).unwrap();
        assert_eq!(d.title(), "t *");
    }
}
