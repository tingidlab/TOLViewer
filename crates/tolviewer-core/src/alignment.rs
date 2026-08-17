//! A set of sequences, aligned or not, plus the column/row edits the GUI needs.

use crate::alphabet::{is_gap, Alphabet, GAP};
use crate::error::{Error, Result};
use crate::sequence::Sequence;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Alignment {
    /// Display name, usually the file stem.
    pub name: String,
    pub sequences: Vec<Sequence>,
    alphabet: Option<Alphabet>,
}

impl Alignment {
    pub fn new(name: impl Into<String>, sequences: Vec<Sequence>) -> Self {
        Alignment { name: name.into(), sequences, alphabet: None }
    }

    /// Number of sequences (rows).
    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }

    /// Number of alignment columns, i.e. the longest row.
    pub fn width(&self) -> usize {
        self.sequences.iter().map(|s| s.len()).max().unwrap_or(0)
    }

    /// True when every row has the same length, so column operations are safe.
    pub fn is_aligned(&self) -> bool {
        let mut lens = self.sequences.iter().map(|s| s.len());
        match lens.next() {
            None => true,
            Some(first) => lens.all(|l| l == first),
        }
    }

    /// The alphabet, guessed from content on first use and cached.
    pub fn alphabet(&mut self) -> Alphabet {
        if let Some(a) = self.alphabet {
            return a;
        }
        // Sample up to 100 rows and 5000 residues each; plenty for a guess and
        // keeps this cheap on huge files.
        let sample =
            self.sequences.iter().take(100).flat_map(|s| s.residues.iter().take(5000).copied());
        let a = Alphabet::guess(sample);
        self.alphabet = Some(a);
        a
    }

    /// The alphabet if already known, without guessing.
    pub fn alphabet_hint(&self) -> Option<Alphabet> {
        self.alphabet
    }

    pub fn set_alphabet(&mut self, alphabet: Alphabet) {
        self.alphabet = Some(alphabet);
    }

    /// Drop the cached alphabet so it is re-guessed after a bulk change.
    pub fn invalidate_alphabet(&mut self) {
        self.alphabet = None;
    }

    /// Right-pad every row with gaps so all rows span [`Alignment::width`].
    pub fn pad_to_width(&mut self) {
        let w = self.width();
        for s in &mut self.sequences {
            s.pad_to(w);
        }
    }

    pub fn require_aligned(&self) -> Result<()> {
        if self.is_aligned() {
            Ok(())
        } else {
            Err(Error::NotAligned)
        }
    }

    pub fn get(&self, row: usize, col: usize) -> Option<u8> {
        self.sequences.get(row)?.residues.get(col).copied()
    }

    /// Overwrite a single residue. Rows shorter than `col` are padded first.
    pub fn set(&mut self, row: usize, col: usize, residue: u8) -> Result<u8> {
        let seq =
            self.sequences.get_mut(row).ok_or_else(|| Error::out_of_range(format!("row {row}")))?;
        if col >= seq.residues.len() {
            seq.pad_to(col + 1);
        }
        let old = seq.residues[col];
        seq.residues[col] = residue;
        Ok(old)
    }

    /// Insert `count` gap columns before column `at`, across all rows.
    pub fn insert_columns(&mut self, at: usize, count: usize) -> Result<()> {
        if count == 0 {
            return Ok(());
        }
        let w = self.width();
        if at > w {
            return Err(Error::out_of_range(format!("column {at} (width {w})")));
        }
        self.pad_to_width();
        for s in &mut self.sequences {
            s.residues.splice(at..at, std::iter::repeat_n(GAP, count));
            if let Some(q) = &mut s.quality {
                if at <= q.len() {
                    q.splice(at..at, std::iter::repeat_n(0, count));
                }
            }
        }
        Ok(())
    }

    /// Delete columns `[start, end)` across all rows, returning what was removed
    /// so the edit can be undone.
    pub fn delete_columns(&mut self, start: usize, end: usize) -> Result<Vec<Vec<u8>>> {
        let w = self.width();
        if start > end || end > w {
            return Err(Error::out_of_range(format!("columns {start}..{end} (width {w})")));
        }
        if start == end {
            return Ok(vec![Vec::new(); self.sequences.len()]);
        }
        self.pad_to_width();
        let mut removed = Vec::with_capacity(self.sequences.len());
        for s in &mut self.sequences {
            removed.push(s.residues.drain(start..end).collect());
            if let Some(q) = &mut s.quality {
                let e = end.min(q.len());
                let b = start.min(e);
                q.drain(b..e);
            }
        }
        Ok(removed)
    }

    /// Re-insert columns previously removed by [`Alignment::delete_columns`].
    pub fn restore_columns(&mut self, at: usize, columns: &[Vec<u8>]) -> Result<()> {
        if columns.len() != self.sequences.len() {
            return Err(Error::out_of_range("restored column count does not match row count"));
        }
        for (s, block) in self.sequences.iter_mut().zip(columns) {
            if at > s.residues.len() {
                s.pad_to(at);
            }
            s.residues.splice(at..at, block.iter().copied());
            if let Some(q) = &mut s.quality {
                if at <= q.len() {
                    q.splice(at..at, std::iter::repeat_n(0, block.len()));
                }
            }
        }
        Ok(())
    }

    /// Keep only the columns whose mask entry is true. `mask` must cover the
    /// full width. Returns the number of columns removed.
    pub fn keep_columns(&mut self, mask: &[bool]) -> Result<usize> {
        let w = self.width();
        if mask.len() != w {
            return Err(Error::out_of_range(format!(
                "mask has {} entries but alignment is {w} columns wide",
                mask.len()
            )));
        }
        self.pad_to_width();
        let kept = mask.iter().filter(|&&k| k).count();
        for s in &mut self.sequences {
            let mut out = Vec::with_capacity(kept);
            for (i, &c) in s.residues.iter().enumerate() {
                if mask[i] {
                    out.push(c);
                }
            }
            s.residues = out;
            if let Some(q) = &mut s.quality {
                let mut oq = Vec::with_capacity(kept);
                for (i, &score) in q.iter().enumerate() {
                    if mask.get(i).copied().unwrap_or(false) {
                        oq.push(score);
                    }
                }
                *q = oq;
            }
        }
        Ok(w - kept)
    }

    /// Insert a single gap into one row at `col`, pushing the rest right.
    pub fn insert_gap(&mut self, row: usize, col: usize) -> Result<()> {
        let seq =
            self.sequences.get_mut(row).ok_or_else(|| Error::out_of_range(format!("row {row}")))?;
        if col > seq.residues.len() {
            seq.pad_to(col);
        }
        seq.residues.insert(col, GAP);
        if let Some(q) = &mut seq.quality {
            if col <= q.len() {
                q.insert(col, 0);
            }
        }
        Ok(())
    }

    /// Delete one position from one row, pulling the rest left. Returns the
    /// removed residue.
    pub fn delete_at(&mut self, row: usize, col: usize) -> Result<u8> {
        let seq =
            self.sequences.get_mut(row).ok_or_else(|| Error::out_of_range(format!("row {row}")))?;
        if col >= seq.residues.len() {
            return Err(Error::out_of_range(format!("column {col} in row {row}")));
        }
        let removed = seq.residues.remove(col);
        if let Some(q) = &mut seq.quality {
            if col < q.len() {
                q.remove(col);
            }
        }
        Ok(removed)
    }

    /// Column indices that are entirely gaps.
    pub fn all_gap_columns(&self) -> Vec<usize> {
        let w = self.width();
        (0..w)
            .filter(|&c| {
                self.sequences.iter().all(|s| s.residues.get(c).is_none_or(|&x| is_gap(x)))
            })
            .collect()
    }

    /// Remove every column that is a gap in all rows. Returns the count removed.
    pub fn remove_all_gap_columns(&mut self) -> usize {
        let w = self.width();
        if w == 0 {
            return 0;
        }
        let mut mask = vec![false; w];
        for (c, keep) in mask.iter_mut().enumerate() {
            *keep = self.sequences.iter().any(|s| s.residues.get(c).is_some_and(|&x| !is_gap(x)));
        }
        self.keep_columns(&mask).unwrap_or(0)
    }

    /// Strip all gaps from all rows, leaving an unaligned set.
    pub fn degap(&mut self) {
        for s in &mut self.sequences {
            let keep: Vec<bool> = s.residues.iter().map(|&c| !is_gap(c)).collect();
            if let Some(q) = &mut s.quality {
                let mut oq = Vec::new();
                for (i, &score) in q.iter().enumerate() {
                    if keep.get(i).copied().unwrap_or(false) {
                        oq.push(score);
                    }
                }
                *q = oq;
            }
            s.residues = s.residues.iter().copied().filter(|&c| !is_gap(c)).collect();
        }
    }

    pub fn remove_sequence(&mut self, row: usize) -> Result<Sequence> {
        if row >= self.sequences.len() {
            return Err(Error::out_of_range(format!("row {row}")));
        }
        Ok(self.sequences.remove(row))
    }

    pub fn insert_sequence(&mut self, row: usize, seq: Sequence) -> Result<()> {
        if row > self.sequences.len() {
            return Err(Error::out_of_range(format!("row {row}")));
        }
        self.sequences.insert(row, seq);
        Ok(())
    }

    /// Move a row to a new index, shifting the rows in between.
    pub fn move_sequence(&mut self, from: usize, to: usize) -> Result<()> {
        if from >= self.sequences.len() || to >= self.sequences.len() {
            return Err(Error::out_of_range(format!("rows {from} -> {to}")));
        }
        let s = self.sequences.remove(from);
        self.sequences.insert(to, s);
        Ok(())
    }

    /// A new alignment holding the given rows and column range. Used for
    /// "export selection" and for previewing a cleaned block.
    pub fn subset(&self, rows: &[usize], cols: std::ops::Range<usize>) -> Alignment {
        let sequences = rows
            .iter()
            .filter_map(|&r| self.sequences.get(r))
            .map(|s| {
                let end = cols.end.min(s.residues.len());
                let start = cols.start.min(end);
                let mut out = s.clone();
                out.residues = s.residues[start..end].to_vec();
                out.quality = s.quality.as_ref().map(|q| {
                    let e = end.min(q.len());
                    let b = start.min(e);
                    q[b..e].to_vec()
                });
                out
            })
            .collect();
        Alignment { name: self.name.clone(), sequences, alphabet: self.alphabet }
    }

    pub fn find_by_id(&self, id: &str) -> Option<usize> {
        self.sequences.iter().position(|s| s.id == id)
    }

    /// Make every `id` unique by appending `_2`, `_3`, ... to duplicates.
    /// PHYLIP and NEXUS both require this. Returns the number renamed.
    pub fn deduplicate_ids(&mut self) -> usize {
        use std::collections::HashMap;
        let mut seen: HashMap<String, usize> = HashMap::new();
        let mut renamed = 0;
        for s in &mut self.sequences {
            let count = {
                let entry = seen.entry(s.id.clone()).or_insert(0);
                *entry += 1;
                *entry
            };
            if count > 1 {
                let mut n = count;
                let mut candidate = format!("{}_{}", s.id, n);
                while seen.contains_key(&candidate) {
                    n += 1;
                    candidate = format!("{}_{}", s.id, n);
                }
                seen.insert(s.id.clone(), n);
                seen.insert(candidate.clone(), 1);
                s.id = candidate;
                renamed += 1;
            }
        }
        renamed
    }

    /// Iterate the residues of one column, padding short rows with gaps.
    pub fn column(&self, col: usize) -> impl Iterator<Item = u8> + '_ {
        self.sequences.iter().map(move |s| s.residues.get(col).copied().unwrap_or(GAP))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aln() -> Alignment {
        Alignment::new(
            "test",
            vec![
                Sequence::new("a", *b"ACGT-A"),
                Sequence::new("b", *b"ACGT-C"),
                Sequence::new("c", *b"AC-T-G"),
            ],
        )
    }

    #[test]
    fn width_and_aligned() {
        let a = aln();
        assert_eq!(a.width(), 6);
        assert_eq!(a.len(), 3);
        assert!(a.is_aligned());
    }

    #[test]
    fn ragged_is_not_aligned_and_pads() {
        let mut a =
            Alignment::new("t", vec![Sequence::new("a", *b"ACGT"), Sequence::new("b", *b"AC")]);
        assert!(!a.is_aligned());
        a.pad_to_width();
        assert!(a.is_aligned());
        assert_eq!(a.sequences[1].residues, b"AC--");
    }

    #[test]
    fn delete_and_restore_columns_round_trips() {
        let mut a = aln();
        let before = a.clone();
        let removed = a.delete_columns(2, 4).unwrap();
        assert_eq!(a.width(), 4);
        assert_eq!(a.sequences[0].residues, b"AC-A");
        a.restore_columns(2, &removed).unwrap();
        assert_eq!(a, before);
    }

    #[test]
    fn removes_all_gap_columns_only() {
        let mut a = aln();
        // column 4 is all gaps; column 2 is not (two G's).
        assert_eq!(a.all_gap_columns(), vec![4]);
        assert_eq!(a.remove_all_gap_columns(), 1);
        assert_eq!(a.sequences[0].residues, b"ACGTA");
        assert_eq!(a.sequences[2].residues, b"AC-TG");
    }

    #[test]
    fn insert_gap_shifts_only_one_row() {
        let mut a = aln();
        a.insert_gap(1, 0).unwrap();
        assert_eq!(a.sequences[1].residues, b"-ACGT-C");
        assert_eq!(a.sequences[0].residues, b"ACGT-A");
        assert!(!a.is_aligned());
    }

    #[test]
    fn delete_at_pulls_left() {
        let mut a = aln();
        assert_eq!(a.delete_at(0, 0).unwrap(), b'A');
        assert_eq!(a.sequences[0].residues, b"CGT-A");
    }

    #[test]
    fn keep_columns_rejects_wrong_mask_length() {
        let mut a = aln();
        assert!(a.keep_columns(&[true; 3]).is_err());
    }

    #[test]
    fn degap_leaves_ragged_rows() {
        let mut a = aln();
        a.degap();
        assert_eq!(a.sequences[0].residues, b"ACGTA");
        assert_eq!(a.sequences[2].residues, b"ACTG");
        assert!(!a.is_aligned());
    }

    #[test]
    fn deduplicate_ids_renames_only_repeats() {
        let mut a = Alignment::new(
            "t",
            vec![
                Sequence::new("x", *b"A"),
                Sequence::new("x", *b"C"),
                Sequence::new("y", *b"G"),
                Sequence::new("x", *b"T"),
            ],
        );
        assert_eq!(a.deduplicate_ids(), 2);
        let ids: Vec<&str> = a.sequences.iter().map(|s| s.id.as_str()).collect();
        assert_eq!(ids, vec!["x", "x_2", "y", "x_3"]);
    }

    #[test]
    fn subset_clips_to_row_length() {
        let a = aln();
        let s = a.subset(&[0, 2], 1..4);
        assert_eq!(s.len(), 2);
        assert_eq!(s.sequences[0].residues, b"CGT");
        assert_eq!(s.sequences[1].residues, b"C-T");
    }
}
