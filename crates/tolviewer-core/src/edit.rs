//! Undoable edits.
//!
//! Every mutation the GUI performs goes through [`UndoStack::apply`], which
//! records enough information to invert it. Bulk operations (align, clean,
//! sort) are recorded as a whole-alignment snapshot, which is simple and fast
//! enough: alignments that fit in the viewer fit in memory several times over.

use crate::alignment::Alignment;
use crate::alphabet::is_gap;
use crate::error::Result;
use crate::sequence::Sequence;

#[derive(Debug, Clone)]
pub enum EditOp {
    /// Overwrite one residue.
    SetResidue { row: usize, col: usize, residue: u8 },
    /// Overwrite a rectangular block, row-major, one Vec per row.
    SetBlock { row: usize, col: usize, residues: Vec<Vec<u8>> },
    /// Insert `count` all-gap columns before `at`.
    InsertColumns { at: usize, count: usize },
    /// Delete columns `[start, end)`.
    DeleteColumns { start: usize, end: usize },
    /// Insert one gap into one row.
    InsertGap { row: usize, col: usize },
    /// Delete one position from one row.
    DeleteAt { row: usize, col: usize },
    /// Remove a whole sequence.
    RemoveSequence { row: usize },
    /// Insert a sequence at `row`.
    InsertSequence { row: usize, seq: Box<Sequence> },
    /// Reorder one row.
    MoveSequence { from: usize, to: usize },
    /// Rename one row.
    Rename { row: usize, id: String, description: String },
    /// Replace the entire alignment (used for align / clean / sort / degap).
    Replace { label: String, alignment: Box<Alignment> },
}

impl EditOp {
    /// Short human-readable label for the Edit menu ("Undo delete columns").
    pub fn label(&self) -> String {
        match self {
            EditOp::SetResidue { .. } => "edit residue".into(),
            EditOp::SetBlock { .. } => "edit block".into(),
            EditOp::InsertColumns { count, .. } => format!("insert {count} column(s)"),
            EditOp::DeleteColumns { start, end } => format!("delete {} column(s)", end - start),
            EditOp::InsertGap { .. } => "insert gap".into(),
            EditOp::DeleteAt { .. } => "delete residue".into(),
            EditOp::RemoveSequence { .. } => "remove sequence".into(),
            EditOp::InsertSequence { .. } => "add sequence".into(),
            EditOp::MoveSequence { .. } => "move sequence".into(),
            EditOp::Rename { .. } => "rename sequence".into(),
            EditOp::Replace { label, .. } => label.clone(),
        }
    }
}

/// The inverse of an applied [`EditOp`], stored on the undo stack.
#[derive(Debug, Clone)]
struct Inverse {
    label: String,
    op: InverseOp,
}

#[derive(Debug, Clone)]
enum InverseOp {
    SetResidue { row: usize, col: usize, residue: u8 },
    SetBlock { row: usize, col: usize, residues: Vec<Vec<u8>> },
    DeleteColumns { start: usize, end: usize },
    RestoreColumns { at: usize, columns: Vec<Vec<u8>> },
    DeleteAt { row: usize, col: usize },
    InsertGapValue { row: usize, col: usize, residue: u8 },
    InsertSequence { row: usize, seq: Box<Sequence> },
    RemoveSequence { row: usize },
    MoveSequence { from: usize, to: usize },
    Rename { row: usize, id: String, description: String },
    Replace { alignment: Box<Alignment> },
    /// Restore each row to a recorded length. Column operations pad ragged rows
    /// to a common width as a side effect; this undoes that padding. Rows are
    /// only shortened when the excess is entirely gaps.
    SetRowLengths(Vec<usize>),
    /// Apply several inverses in order.
    Compound(Vec<InverseOp>),
}

#[derive(Debug, Default)]
pub struct UndoStack {
    undo: Vec<Inverse>,
    redo: Vec<Inverse>,
    limit: usize,
    /// Incremented on every successful mutation; the GUI compares it against
    /// the value at last save to decide whether the document is dirty, and
    /// against cached derived data (consensus, masks) to invalidate it.
    revision: u64,
}

impl UndoStack {
    pub fn new() -> Self {
        UndoStack { undo: Vec::new(), redo: Vec::new(), limit: 200, revision: 0 }
    }

    pub fn with_limit(limit: usize) -> Self {
        UndoStack { limit, ..UndoStack::new() }
    }

    pub fn revision(&self) -> u64 {
        self.revision
    }

    pub fn can_undo(&self) -> bool {
        !self.undo.is_empty()
    }

    pub fn can_redo(&self) -> bool {
        !self.redo.is_empty()
    }

    pub fn undo_label(&self) -> Option<&str> {
        self.undo.last().map(|i| i.label.as_str())
    }

    pub fn redo_label(&self) -> Option<&str> {
        self.redo.last().map(|i| i.label.as_str())
    }

    pub fn clear(&mut self) {
        self.undo.clear();
        self.redo.clear();
    }

    /// Apply `op` to `aln` and push its inverse. A new edit discards the redo
    /// stack, as in every other editor.
    pub fn apply(&mut self, aln: &mut Alignment, op: EditOp) -> Result<()> {
        let label = op.label();
        // Column operations pad ragged rows to a common width; record the
        // shape first so undo can put it back exactly.
        let shape: Vec<usize> = aln.sequences.iter().map(|s| s.len()).collect();
        let primary = Self::apply_inner(aln, op)?;
        let inverse = if aln.sequences.len() == shape.len()
            && aln.sequences.iter().zip(&shape).any(|(s, &l)| s.len() != l)
        {
            InverseOp::Compound(vec![primary, InverseOp::SetRowLengths(shape)])
        } else {
            primary
        };
        self.redo.clear();
        self.undo.push(Inverse { label, op: inverse });
        if self.undo.len() > self.limit {
            self.undo.remove(0);
        }
        self.revision += 1;
        Ok(())
    }

    pub fn undo(&mut self, aln: &mut Alignment) -> Result<Option<String>> {
        let Some(entry) = self.undo.pop() else { return Ok(None) };
        let label = entry.label.clone();
        let redo = Self::apply_inverse(aln, entry.op)?;
        self.redo.push(Inverse { label: label.clone(), op: redo });
        self.revision += 1;
        Ok(Some(label))
    }

    pub fn redo(&mut self, aln: &mut Alignment) -> Result<Option<String>> {
        let Some(entry) = self.redo.pop() else { return Ok(None) };
        let label = entry.label.clone();
        let undo = Self::apply_inverse(aln, entry.op)?;
        self.undo.push(Inverse { label: label.clone(), op: undo });
        self.revision += 1;
        Ok(Some(label))
    }

    fn apply_inner(aln: &mut Alignment, op: EditOp) -> Result<InverseOp> {
        Ok(match op {
            EditOp::SetResidue { row, col, residue } => {
                let old = aln.set(row, col, residue)?;
                InverseOp::SetResidue { row, col, residue: old }
            }
            EditOp::SetBlock { row, col, residues } => {
                let mut old = Vec::with_capacity(residues.len());
                for (dr, block) in residues.iter().enumerate() {
                    let mut old_row = Vec::with_capacity(block.len());
                    for (dc, &c) in block.iter().enumerate() {
                        old_row.push(aln.set(row + dr, col + dc, c)?);
                    }
                    old.push(old_row);
                }
                InverseOp::SetBlock { row, col, residues: old }
            }
            EditOp::InsertColumns { at, count } => {
                aln.insert_columns(at, count)?;
                InverseOp::DeleteColumns { start: at, end: at + count }
            }
            EditOp::DeleteColumns { start, end } => {
                let removed = aln.delete_columns(start, end)?;
                InverseOp::RestoreColumns { at: start, columns: removed }
            }
            EditOp::InsertGap { row, col } => {
                aln.insert_gap(row, col)?;
                InverseOp::DeleteAt { row, col }
            }
            EditOp::DeleteAt { row, col } => {
                let removed = aln.delete_at(row, col)?;
                InverseOp::InsertGapValue { row, col, residue: removed }
            }
            EditOp::RemoveSequence { row } => {
                let seq = aln.remove_sequence(row)?;
                InverseOp::InsertSequence { row, seq: Box::new(seq) }
            }
            EditOp::InsertSequence { row, seq } => {
                aln.insert_sequence(row, *seq)?;
                InverseOp::RemoveSequence { row }
            }
            EditOp::MoveSequence { from, to } => {
                aln.move_sequence(from, to)?;
                InverseOp::MoveSequence { from: to, to: from }
            }
            EditOp::Rename { row, id, description } => {
                let seq = aln
                    .sequences
                    .get_mut(row)
                    .ok_or_else(|| crate::Error::out_of_range(format!("row {row}")))?;
                let old = InverseOp::Rename {
                    row,
                    id: std::mem::replace(&mut seq.id, id),
                    description: std::mem::replace(&mut seq.description, description),
                };
                old
            }
            EditOp::Replace { alignment, .. } => {
                let old = std::mem::replace(aln, *alignment);
                InverseOp::Replace { alignment: Box::new(old) }
            }
        })
    }

    fn apply_inverse(aln: &mut Alignment, op: InverseOp) -> Result<InverseOp> {
        Ok(match op {
            InverseOp::SetResidue { row, col, residue } => {
                let old = aln.set(row, col, residue)?;
                InverseOp::SetResidue { row, col, residue: old }
            }
            InverseOp::SetBlock { row, col, residues } => {
                let mut old = Vec::with_capacity(residues.len());
                for (dr, block) in residues.iter().enumerate() {
                    let mut old_row = Vec::with_capacity(block.len());
                    for (dc, &c) in block.iter().enumerate() {
                        old_row.push(aln.set(row + dr, col + dc, c)?);
                    }
                    old.push(old_row);
                }
                InverseOp::SetBlock { row, col, residues: old }
            }
            InverseOp::DeleteColumns { start, end } => {
                let removed = aln.delete_columns(start, end)?;
                InverseOp::RestoreColumns { at: start, columns: removed }
            }
            InverseOp::RestoreColumns { at, columns } => {
                let count = columns.first().map_or(0, |c| c.len());
                aln.restore_columns(at, &columns)?;
                InverseOp::DeleteColumns { start: at, end: at + count }
            }
            InverseOp::DeleteAt { row, col } => {
                let removed = aln.delete_at(row, col)?;
                InverseOp::InsertGapValue { row, col, residue: removed }
            }
            InverseOp::InsertGapValue { row, col, residue } => {
                aln.insert_gap(row, col)?;
                aln.set(row, col, residue)?;
                InverseOp::DeleteAt { row, col }
            }
            InverseOp::InsertSequence { row, seq } => {
                aln.insert_sequence(row, *seq)?;
                InverseOp::RemoveSequence { row }
            }
            InverseOp::RemoveSequence { row } => {
                let seq = aln.remove_sequence(row)?;
                InverseOp::InsertSequence { row, seq: Box::new(seq) }
            }
            InverseOp::MoveSequence { from, to } => {
                aln.move_sequence(from, to)?;
                InverseOp::MoveSequence { from: to, to: from }
            }
            InverseOp::Rename { row, id, description } => {
                let seq = aln
                    .sequences
                    .get_mut(row)
                    .ok_or_else(|| crate::Error::out_of_range(format!("row {row}")))?;
                InverseOp::Rename {
                    row,
                    id: std::mem::replace(&mut seq.id, id),
                    description: std::mem::replace(&mut seq.description, description),
                }
            }
            InverseOp::Replace { alignment } => {
                let old = std::mem::replace(aln, *alignment);
                InverseOp::Replace { alignment: Box::new(old) }
            }
            InverseOp::SetRowLengths(lengths) => {
                let current: Vec<usize> = aln.sequences.iter().map(|s| s.len()).collect();
                for (s, &want) in aln.sequences.iter_mut().zip(&lengths) {
                    if s.len() < want {
                        s.pad_to(want);
                    } else if s.len() > want && s.residues[want..].iter().all(|&c| is_gap(c)) {
                        // Only ever drop trailing padding, never real residues.
                        s.residues.truncate(want);
                        if let Some(q) = &mut s.quality {
                            q.truncate(want.min(q.len()));
                        }
                    }
                }
                InverseOp::SetRowLengths(current)
            }
            InverseOp::Compound(ops) => {
                let mut inverses = Vec::with_capacity(ops.len());
                for op in ops {
                    inverses.push(Self::apply_inverse(aln, op)?);
                }
                inverses.reverse();
                InverseOp::Compound(inverses)
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aln() -> Alignment {
        Alignment::new(
            "t",
            vec![Sequence::new("a", *b"ACGT"), Sequence::new("b", *b"AGGT"), Sequence::new("c", *b"ATGT")],
        )
    }

    #[test]
    fn set_residue_undo_redo() {
        let mut a = aln();
        let mut u = UndoStack::new();
        u.apply(&mut a, EditOp::SetResidue { row: 0, col: 1, residue: b'T' }).unwrap();
        assert_eq!(a.sequences[0].residues, b"ATGT");
        u.undo(&mut a).unwrap();
        assert_eq!(a.sequences[0].residues, b"ACGT");
        u.redo(&mut a).unwrap();
        assert_eq!(a.sequences[0].residues, b"ATGT");
    }

    #[test]
    fn delete_columns_undo_restores_exactly() {
        let mut a = aln();
        let before = a.clone();
        let mut u = UndoStack::new();
        u.apply(&mut a, EditOp::DeleteColumns { start: 1, end: 3 }).unwrap();
        assert_eq!(a.width(), 2);
        u.undo(&mut a).unwrap();
        assert_eq!(a, before);
        u.redo(&mut a).unwrap();
        assert_eq!(a.width(), 2);
    }

    #[test]
    fn delete_at_undo_restores_the_residue_not_a_gap() {
        let mut a = aln();
        let mut u = UndoStack::new();
        u.apply(&mut a, EditOp::DeleteAt { row: 1, col: 1 }).unwrap();
        assert_eq!(a.sequences[1].residues, b"AGT");
        u.undo(&mut a).unwrap();
        assert_eq!(a.sequences[1].residues, b"AGGT");
    }

    #[test]
    fn remove_sequence_undo_restores_position() {
        let mut a = aln();
        let mut u = UndoStack::new();
        u.apply(&mut a, EditOp::RemoveSequence { row: 1 }).unwrap();
        assert_eq!(a.len(), 2);
        u.undo(&mut a).unwrap();
        assert_eq!(a.len(), 3);
        assert_eq!(a.sequences[1].id, "b");
    }

    #[test]
    fn replace_round_trips_and_labels() {
        let mut a = aln();
        let before = a.clone();
        let mut u = UndoStack::new();
        let new = Alignment::new("t", vec![Sequence::new("z", *b"AAAA")]);
        u.apply(&mut a, EditOp::Replace { label: "align".into(), alignment: Box::new(new) }).unwrap();
        assert_eq!(u.undo_label(), Some("align"));
        u.undo(&mut a).unwrap();
        assert_eq!(a, before);
    }

    #[test]
    fn new_edit_clears_redo() {
        let mut a = aln();
        let mut u = UndoStack::new();
        u.apply(&mut a, EditOp::SetResidue { row: 0, col: 0, residue: b'G' }).unwrap();
        u.undo(&mut a).unwrap();
        assert!(u.can_redo());
        u.apply(&mut a, EditOp::SetResidue { row: 0, col: 1, residue: b'G' }).unwrap();
        assert!(!u.can_redo());
    }

    #[test]
    fn many_edits_undo_back_to_start() {
        let mut a = aln();
        let before = a.clone();
        let mut u = UndoStack::new();
        u.apply(&mut a, EditOp::InsertGap { row: 0, col: 2 }).unwrap();
        u.apply(&mut a, EditOp::SetResidue { row: 1, col: 0, residue: b'G' }).unwrap();
        u.apply(&mut a, EditOp::InsertColumns { at: 0, count: 2 }).unwrap();
        u.apply(&mut a, EditOp::MoveSequence { from: 0, to: 2 }).unwrap();
        while u.can_undo() {
            u.undo(&mut a).unwrap();
        }
        assert_eq!(a, before);
    }
}
