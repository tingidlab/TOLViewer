//! Cursor and selection over the alignment grid.

use std::ops::Range;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Cell {
    pub row: usize,
    pub col: usize,
}

impl Cell {
    pub fn new(row: usize, col: usize) -> Self {
        Cell { row, col }
    }
}

/// What a drag selects. Column and row modes let the user grab whole
/// columns/rows by dragging in the rulers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SelectionMode {
    Cells,
    Columns,
    Rows,
}

#[derive(Debug, Clone)]
pub struct Selection {
    /// Where the current drag or shift-extension started.
    pub anchor: Cell,
    /// Where the caret is now. Edits happen here.
    pub cursor: Cell,
    pub mode: SelectionMode,
    /// False when only the caret is placed and nothing is highlighted.
    pub active: bool,
}

impl Default for Selection {
    fn default() -> Self {
        Selection {
            anchor: Cell::default(),
            cursor: Cell::default(),
            mode: SelectionMode::Cells,
            active: false,
        }
    }
}

impl Selection {
    /// Place the caret, collapsing any selection.
    pub fn place(&mut self, cell: Cell, mode: SelectionMode) {
        self.anchor = cell;
        self.cursor = cell;
        self.mode = mode;
        self.active = false;
    }

    /// Extend the selection to `cell`, keeping the anchor.
    pub fn extend_to(&mut self, cell: Cell) {
        self.cursor = cell;
        self.active = true;
    }

    pub fn clear(&mut self) {
        self.active = false;
        self.anchor = self.cursor;
    }

    /// Selected rows, inclusive of both ends. In `Columns` mode this covers
    /// every row, so callers must clamp with the alignment's row count.
    pub fn rows(&self, row_count: usize) -> Range<usize> {
        if row_count == 0 {
            return 0..0;
        }
        match self.mode {
            SelectionMode::Columns => 0..row_count,
            _ => {
                let (a, b) = minmax(self.anchor.row, self.cursor.row);
                a.min(row_count - 1)..(b + 1).min(row_count)
            }
        }
    }

    /// Selected columns, half-open. In `Rows` mode this covers the full width.
    pub fn cols(&self, width: usize) -> Range<usize> {
        if width == 0 {
            return 0..0;
        }
        match self.mode {
            SelectionMode::Rows => 0..width,
            _ => {
                let (a, b) = minmax(self.anchor.col, self.cursor.col);
                a.min(width - 1)..(b + 1).min(width)
            }
        }
    }

    /// True when `cell` lies inside the highlighted region.
    pub fn contains(&self, cell: Cell, rows: usize, cols: usize) -> bool {
        if !self.active {
            return false;
        }
        self.rows(rows).contains(&cell.row) && self.cols(cols).contains(&cell.col)
    }

    /// Number of selected cells, for the status bar.
    pub fn cell_count(&self, rows: usize, cols: usize) -> usize {
        if !self.active {
            return 0;
        }
        self.rows(rows).len() * self.cols(cols).len()
    }

    /// Move the caret, clamping to the grid. `extend` keeps the anchor so
    /// shift+arrow grows the selection.
    pub fn move_caret(&mut self, drow: isize, dcol: isize, rows: usize, cols: usize, extend: bool) {
        if rows == 0 || cols == 0 {
            return;
        }
        let row = clamp_step(self.cursor.row, drow, rows);
        let col = clamp_step(self.cursor.col, dcol, cols);
        self.cursor = Cell::new(row, col);
        if extend {
            self.active = true;
        } else {
            self.anchor = self.cursor;
            self.active = false;
        }
    }
}

fn minmax(a: usize, b: usize) -> (usize, usize) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

fn clamp_step(value: usize, delta: isize, limit: usize) -> usize {
    let next = value as isize + delta;
    next.clamp(0, limit as isize - 1) as usize
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn place_collapses_the_selection() {
        let mut s = Selection::default();
        s.extend_to(Cell::new(5, 5));
        assert!(s.active);
        s.place(Cell::new(1, 1), SelectionMode::Cells);
        assert!(!s.active);
        assert_eq!(s.cell_count(10, 10), 0);
    }

    #[test]
    fn selection_normalises_backwards_drags() {
        let mut s = Selection::default();
        s.place(Cell::new(4, 8), SelectionMode::Cells);
        s.extend_to(Cell::new(2, 3));
        assert_eq!(s.rows(10), 2..5);
        assert_eq!(s.cols(20), 3..9);
    }

    #[test]
    fn column_mode_spans_every_row() {
        let mut s = Selection::default();
        s.place(Cell::new(3, 2), SelectionMode::Columns);
        s.extend_to(Cell::new(3, 4));
        assert_eq!(s.rows(7), 0..7);
        assert_eq!(s.cols(20), 2..5);
    }

    #[test]
    fn row_mode_spans_every_column() {
        let mut s = Selection::default();
        s.place(Cell::new(1, 0), SelectionMode::Rows);
        s.extend_to(Cell::new(2, 0));
        assert_eq!(s.rows(7), 1..3);
        assert_eq!(s.cols(20), 0..20);
    }

    #[test]
    fn selection_is_clamped_to_the_grid() {
        let mut s = Selection::default();
        s.place(Cell::new(0, 0), SelectionMode::Cells);
        s.extend_to(Cell::new(99, 99));
        assert_eq!(s.rows(3), 0..3);
        assert_eq!(s.cols(5), 0..5);
    }

    #[test]
    fn caret_stops_at_the_edges() {
        let mut s = Selection::default();
        s.move_caret(-1, -1, 4, 4, false);
        assert_eq!(s.cursor, Cell::new(0, 0));
        s.move_caret(100, 100, 4, 4, false);
        assert_eq!(s.cursor, Cell::new(3, 3));
    }

    #[test]
    fn shift_arrow_extends_from_the_anchor() {
        let mut s = Selection::default();
        s.place(Cell::new(1, 1), SelectionMode::Cells);
        s.move_caret(0, 1, 5, 5, true);
        s.move_caret(0, 1, 5, 5, true);
        assert!(s.active);
        assert_eq!(s.cols(5), 1..4);
        assert_eq!(s.anchor, Cell::new(1, 1));
    }

    #[test]
    fn empty_alignment_yields_empty_ranges() {
        let s = Selection::default();
        assert_eq!(s.rows(0), 0..0);
        assert_eq!(s.cols(0), 0..0);
    }
}
