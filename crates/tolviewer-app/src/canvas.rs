//! The alignment canvas: a virtualised grid of residues with a frozen name
//! gutter on the left and frozen ruler/quality/consensus tracks on top.
//!
//! Nothing off screen is ever touched. Per frame the canvas paints at most
//! (visible rows x visible columns) cells, and adjacent cells sharing a colour
//! are merged into a single rectangle, so a 5000 x 100000 alignment costs the
//! same as a small one.

use egui::{
    Align2, Color32, CursorIcon, FontFamily, FontId, Pos2, Rect, Response, Sense, Stroke,
    StrokeKind, UiBuilder, Vec2,
};
use tolviewer_core::{is_gap, EditOp, GAP};

use crate::document::Document;
use crate::selection::{Cell, SelectionMode};
use crate::theme::{background, ink_for, CellContext, ColorScheme};

/// Below this cell width the letters are unreadable, so the canvas paints
/// colour blocks only. Above it, every cell gets its glyph.
const GLYPH_MIN_WIDTH: f32 = 6.0;
const RULER_H: f32 = 18.0;
const QUALITY_H: f32 = 26.0;
const MASK_H: f32 = 8.0;
pub const MIN_ZOOM: f32 = 1.0;
pub const MAX_ZOOM: f32 = 40.0;

/// View settings shared by every document.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ViewSettings {
    pub scheme: ColorScheme,
    pub show_quality_track: bool,
    pub show_consensus: bool,
    pub show_ruler: bool,
    /// Width of the name gutter in points.
    pub name_width: f32,
    /// Draw residues that match the consensus as dots.
    pub dot_matches: bool,
}

impl Default for ViewSettings {
    fn default() -> Self {
        ViewSettings {
            scheme: ColorScheme::Residue,
            show_quality_track: true,
            show_consensus: true,
            show_ruler: true,
            name_width: 170.0,
            dot_matches: false,
        }
    }
}

/// What the canvas wants the application to do. Edits are returned rather than
/// applied so the app can report errors in one place and keep the undo stack's
/// ownership clear.
#[derive(Debug)]
pub enum CanvasAction {
    Edit(EditOp),
    /// The user double-clicked a name and wants to rename that row.
    BeginRename(usize),
    /// Scroll the caret into view on the next frame.
    ScrollToCaret,
}

pub struct AlignmentCanvas<'a> {
    pub doc: &'a mut Document,
    pub view: &'a ViewSettings,
    pub id_salt: egui::Id,
}

impl AlignmentCanvas<'_> {
    /// Draw the canvas and collect the user's intent.
    pub fn show(self, ui: &mut egui::Ui, actions: &mut Vec<CanvasAction>) -> Response {
        let AlignmentCanvas { doc, view, id_salt } = self;
        let dark = ui.visuals().dark_mode;
        let full = ui.available_rect_before_wrap();

        let cell_w = doc.zoom;
        let cell_h = (doc.zoom * 1.4).max(4.0);
        let rows = doc.rows();
        let cols = doc.width();

        let has_mask = doc.live_clean_mask().is_some();
        let header_h = (if view.show_ruler { RULER_H } else { 0.0 })
            + (if view.show_quality_track { QUALITY_H } else { 0.0 })
            + (if view.show_consensus { cell_h } else { 0.0 })
            + (if has_mask { MASK_H } else { 0.0 });

        // The gutter takes at most half the window, but in a window narrower
        // than twice the minimum gutter that upper bound falls below the lower
        // one, so raise it rather than inverting the range.
        const MIN_GUTTER: f32 = 60.0;
        let max_gutter = (full.width() * 0.5).max(MIN_GUTTER);
        let gutter_w = view.name_width.clamp(MIN_GUTTER, max_gutter).min(full.width().max(0.0));
        let grid_rect = Rect::from_min_max(full.min + Vec2::new(gutter_w, header_h), full.max);
        if grid_rect.width() <= 1.0 || grid_rect.height() <= 1.0 {
            return ui.allocate_rect(full, Sense::hover());
        }

        // --- the scrolling residue grid -------------------------------------
        let mut grid = GridPaint {
            offset: Vec2::ZERO,
            visible_cols: 0..0,
            visible_rows: 0..0,
            origin: grid_rect.min,
        };
        let mut response = ui.allocate_rect(full, Sense::hover());

        let scroll_out = ui
            .scope_builder(UiBuilder::new().max_rect(grid_rect), |ui| {
                egui::ScrollArea::both()
                    .id_salt(id_salt.with("grid"))
                    .auto_shrink([false, false])
                    .show_viewport(ui, |ui, viewport| {
                        let content = Vec2::new(cols as f32 * cell_w, rows as f32 * cell_h);
                        let (rect, resp) = ui.allocate_exact_size(content, Sense::click_and_drag());
                        grid.origin = rect.min;
                        grid.visible_cols =
                            visible_range(viewport.min.x, viewport.max.x, cell_w, cols);
                        grid.visible_rows =
                            visible_range(viewport.min.y, viewport.max.y, cell_h, rows);
                        paint_grid(ui, doc, view, &grid, cell_w, cell_h, dark);
                        resp
                    })
            })
            .inner;
        grid.offset = scroll_out.state.offset;
        let grid_response = scroll_out.inner;

        handle_grid_input(doc, &grid_response, &grid, cell_w, cell_h, actions);

        // --- frozen panes, painted after the grid so they sit on top ---------
        let header_rect = Rect::from_min_max(
            Pos2::new(full.min.x + gutter_w, full.min.y),
            Pos2::new(full.max.x, full.min.y + header_h),
        );
        let names_rect = Rect::from_min_max(
            Pos2::new(full.min.x, full.min.y + header_h),
            Pos2::new(full.min.x + gutter_w, full.max.y),
        );
        let corner_rect =
            Rect::from_min_max(full.min, Pos2::new(full.min.x + gutter_w, full.min.y + header_h));

        let header_response =
            paint_header(ui, doc, view, header_rect, &grid, cell_w, cell_h, dark, id_salt);
        let names_response = paint_names(ui, doc, view, names_rect, &grid, cell_h, dark, id_salt);
        paint_corner(ui, corner_rect, dark);

        handle_header_input(doc, &header_response, header_rect, &grid, cell_w);
        handle_names_input(doc, &names_response, names_rect, &grid, cell_h, actions);

        response |= grid_response;
        response |= header_response;
        response |= names_response;
        response
    }
}

/// Geometry shared between the grid and the frozen panes.
struct GridPaint {
    offset: Vec2,
    visible_cols: std::ops::Range<usize>,
    visible_rows: std::ops::Range<usize>,
    /// Screen position of content cell (0, 0), already shifted by the scroll.
    origin: Pos2,
}

fn visible_range(min: f32, max: f32, step: f32, count: usize) -> std::ops::Range<usize> {
    if count == 0 || step <= 0.0 {
        return 0..0;
    }
    let first = (min / step).floor().max(0.0) as usize;
    // One extra cell each way so partially visible cells are still painted.
    let last = ((max / step).ceil() as usize + 1).min(count);
    first.min(count)..last
}

fn paint_grid(
    ui: &egui::Ui,
    doc: &mut Document,
    view: &ViewSettings,
    grid: &GridPaint,
    cell_w: f32,
    cell_h: f32,
    dark: bool,
) {
    let alphabet = doc.alphabet();
    let scheme = view.scheme;
    let draw_glyphs = cell_w >= GLYPH_MIN_WIDTH;
    let font = FontId::new((cell_w * 1.5).min(cell_h * 0.82).max(4.0), FontFamily::Monospace);

    // Consensus is needed by two schemes and by dot-matching; take a copy of
    // just the visible slice so the borrow on `doc` ends here.
    let needs_consensus = scheme.needs_consensus() || view.dot_matches;
    let (consensus, stats) = if needs_consensus {
        let c = doc.consensus();
        (
            c.residues[grid.visible_cols.start.min(c.residues.len())
                ..grid.visible_cols.end.min(c.residues.len())]
                .to_vec(),
            c.columns[grid.visible_cols.start.min(c.columns.len())
                ..grid.visible_cols.end.min(c.columns.len())]
                .to_vec(),
        )
    } else {
        (Vec::new(), Vec::new())
    };

    let painter = ui.painter();
    let mut shapes = Vec::with_capacity(grid.visible_rows.len() * 8);
    let mut glyphs = Vec::new();

    for row in grid.visible_rows.clone() {
        let Some(seq) = doc.alignment.sequences.get(row) else { continue };
        let y = grid.origin.y + row as f32 * cell_h;
        let dim_row = seq.hidden;

        // Coalesce runs of identical colour into one rectangle.
        let mut run_start: Option<(usize, Option<Color32>)> = None;
        for col in grid.visible_cols.clone() {
            let residue = seq.residues.get(col).copied().unwrap_or(GAP);
            let idx = col - grid.visible_cols.start;
            let cx = CellContext {
                residue,
                alphabet,
                consensus: consensus.get(idx).copied(),
                stats: stats.get(idx),
                quality: seq.quality.as_ref().and_then(|q| q.get(col).copied()),
                dark,
            };
            let mut bg = background(scheme, &cx);
            if dim_row {
                bg = bg.map(|c| c.gamma_multiply(0.45));
            }
            match &run_start {
                Some((_, prev)) if *prev == bg => {}
                _ => {
                    if let Some((start, Some(color))) = run_start {
                        shapes.push(egui::Shape::rect_filled(
                            run_rect(grid.origin.x, start, col, cell_w, y, cell_h),
                            0.0,
                            color,
                        ));
                    }
                    run_start = Some((col, bg));
                }
            }

            if draw_glyphs && !is_gap(residue) {
                let shown = if view.dot_matches
                    && consensus.get(idx).is_some_and(|&c| c.eq_ignore_ascii_case(&residue))
                {
                    b'.'
                } else {
                    residue
                };
                glyphs.push((
                    Pos2::new(grid.origin.x + col as f32 * cell_w + cell_w * 0.5, y + cell_h * 0.5),
                    shown as char,
                    ink_for(bg, dark),
                ));
            }
        }
        if let Some((start, Some(color))) = run_start {
            shapes.push(egui::Shape::rect_filled(
                run_rect(grid.origin.x, start, grid.visible_cols.end, cell_w, y, cell_h),
                0.0,
                color,
            ));
        }
    }
    painter.extend(shapes);

    // Glyphs go in a second pass so they are never covered by a later run.
    let mut buf = [0u8; 4];
    for (pos, ch, color) in glyphs {
        painter.text(pos, Align2::CENTER_CENTER, ch.encode_utf8(&mut buf), font.clone(), color);
    }

    paint_selection(ui, doc, grid, cell_w, cell_h);
}

fn run_rect(origin_x: f32, start: usize, end: usize, cell_w: f32, y: f32, cell_h: f32) -> Rect {
    Rect::from_min_max(
        Pos2::new(origin_x + start as f32 * cell_w, y),
        Pos2::new(origin_x + end as f32 * cell_w, y + cell_h),
    )
}

fn paint_selection(ui: &egui::Ui, doc: &Document, grid: &GridPaint, cell_w: f32, cell_h: f32) {
    let rows = doc.rows();
    let cols = doc.width();
    if rows == 0 || cols == 0 {
        return;
    }
    let painter = ui.painter();
    let accent = ui.visuals().selection.bg_fill;

    if doc.selection.active {
        let r = doc.selection.rows(rows);
        let c = doc.selection.cols(cols);
        let rect = Rect::from_min_max(
            Pos2::new(
                grid.origin.x + c.start as f32 * cell_w,
                grid.origin.y + r.start as f32 * cell_h,
            ),
            Pos2::new(grid.origin.x + c.end as f32 * cell_w, grid.origin.y + r.end as f32 * cell_h),
        );
        painter.rect_filled(rect, 0.0, accent.gamma_multiply(0.28));
        painter.rect_stroke(rect, 0.0, Stroke::new(1.0, accent), StrokeKind::Inside);
    }

    let cur = doc.selection.cursor;
    if cur.row < rows && cur.col < cols {
        let rect = Rect::from_min_size(
            Pos2::new(
                grid.origin.x + cur.col as f32 * cell_w,
                grid.origin.y + cur.row as f32 * cell_h,
            ),
            Vec2::new(cell_w, cell_h),
        );
        painter.rect_stroke(
            rect,
            0.0,
            Stroke::new(1.5, ui.visuals().strong_text_color()),
            StrokeKind::Outside,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn paint_header(
    ui: &mut egui::Ui,
    doc: &mut Document,
    view: &ViewSettings,
    rect: Rect,
    grid: &GridPaint,
    cell_w: f32,
    cell_h: f32,
    dark: bool,
    id_salt: egui::Id,
) -> Response {
    let response = ui.interact(rect, id_salt.with("header"), Sense::click_and_drag());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, ui.visuals().extreme_bg_color);

    let cols = doc.width();
    let x_of = |col: usize| rect.min.x + col as f32 * cell_w - grid.offset.x;
    let mut y = rect.min.y;

    if view.show_ruler {
        let ruler =
            Rect::from_min_max(Pos2::new(rect.min.x, y), Pos2::new(rect.max.x, y + RULER_H));
        let step = tick_step(cell_w);
        let font = FontId::new(10.0, FontFamily::Proportional);
        let ink = ui.visuals().weak_text_color();
        // Start at the first multiple of `step` at or after the first visible
        // column, counting from 1 as biologists number alignment positions.
        let first = grid.visible_cols.start.max(1);
        let mut col = first.div_ceil(step) * step;
        while col <= grid.visible_cols.end && col <= cols {
            let x = x_of(col - 1) + cell_w * 0.5;
            if x >= rect.min.x && x <= rect.max.x {
                painter.line_segment(
                    [Pos2::new(x, ruler.max.y - 4.0), Pos2::new(x, ruler.max.y)],
                    Stroke::new(1.0, ink),
                );
                painter.text(
                    Pos2::new(x, ruler.min.y + 1.0),
                    Align2::CENTER_TOP,
                    col.to_string(),
                    font.clone(),
                    ink,
                );
            }
            col += step;
        }
        y = ruler.max.y;
    }

    if view.show_quality_track {
        let track =
            Rect::from_min_max(Pos2::new(rect.min.x, y), Pos2::new(rect.max.x, y + QUALITY_H));
        let stats: Vec<f32> = {
            let c = doc.consensus();
            grid.visible_cols
                .clone()
                .map(|i| c.columns.get(i).map_or(0.0, |s| s.quality()))
                .collect()
        };
        let mut bars = Vec::with_capacity(stats.len());
        for (i, &q) in stats.iter().enumerate() {
            let col = grid.visible_cols.start + i;
            let h = (track.height() - 2.0) * q.clamp(0.0, 1.0);
            let bar = Rect::from_min_max(
                Pos2::new(x_of(col), track.max.y - 1.0 - h),
                Pos2::new(x_of(col) + cell_w.max(1.0), track.max.y - 1.0),
            );
            bars.push(egui::Shape::rect_filled(bar, 0.0, quality_color(q, dark)));
        }
        painter.extend(bars);
        y = track.max.y;
    }

    if view.show_consensus {
        let track = Rect::from_min_max(Pos2::new(rect.min.x, y), Pos2::new(rect.max.x, y + cell_h));
        let alphabet = doc.alphabet();
        let residues: Vec<u8> = {
            let c = doc.consensus();
            grid.visible_cols.clone().map(|i| c.residues.get(i).copied().unwrap_or(GAP)).collect()
        };
        let font = FontId::new((cell_w * 1.5).min(cell_h * 0.82).max(4.0), FontFamily::Monospace);
        let mut buf = [0u8; 4];
        for (i, &residue) in residues.iter().enumerate() {
            let col = grid.visible_cols.start + i;
            let cx = CellContext {
                residue,
                alphabet,
                consensus: None,
                stats: None,
                quality: None,
                dark,
            };
            let bg = background(view.scheme, &cx);
            let cell =
                Rect::from_min_size(Pos2::new(x_of(col), track.min.y), Vec2::new(cell_w, cell_h));
            if let Some(bg) = bg {
                painter.rect_filled(cell, 0.0, bg.gamma_multiply(0.7));
            }
            if cell_w >= GLYPH_MIN_WIDTH {
                painter.text(
                    cell.center(),
                    Align2::CENTER_CENTER,
                    (residue as char).encode_utf8(&mut buf),
                    font.clone(),
                    ink_for(bg, dark),
                );
            }
        }
        y = track.max.y;
    }

    if let Some(mask) = doc.live_clean_mask() {
        let track = Rect::from_min_max(Pos2::new(rect.min.x, y), Pos2::new(rect.max.x, y + MASK_H));
        let keep = Color32::from_rgb(0x4C, 0xAF, 0x6E);
        let drop = ui.visuals().widgets.inactive.bg_fill;
        let mut shapes = Vec::new();
        for col in grid.visible_cols.clone() {
            let color = if mask.get(col).copied().unwrap_or(false) { keep } else { drop };
            shapes.push(egui::Shape::rect_filled(
                Rect::from_min_max(
                    Pos2::new(x_of(col), track.min.y + 1.0),
                    Pos2::new(x_of(col) + cell_w.max(1.0), track.max.y - 1.0),
                ),
                0.0,
                color,
            ));
        }
        painter.extend(shapes);
    }

    // Highlight the selected columns across every track.
    if doc.selection.active && cols > 0 {
        let c = doc.selection.cols(cols);
        let sel = Rect::from_min_max(
            Pos2::new(x_of(c.start), rect.min.y),
            Pos2::new(x_of(c.end), rect.max.y),
        );
        painter.rect_filled(sel, 0.0, ui.visuals().selection.bg_fill.gamma_multiply(0.18));
    }
    painter.line_segment(
        [Pos2::new(rect.min.x, rect.max.y - 0.5), Pos2::new(rect.max.x, rect.max.y - 0.5)],
        ui.visuals().widgets.noninteractive.bg_stroke,
    );
    if response.hovered() {
        ui.ctx().set_cursor_icon(CursorIcon::ResizeHorizontal);
    }
    response
}

#[allow(clippy::too_many_arguments)]
fn paint_names(
    ui: &mut egui::Ui,
    doc: &Document,
    _view: &ViewSettings,
    rect: Rect,
    grid: &GridPaint,
    cell_h: f32,
    _dark: bool,
    id_salt: egui::Id,
) -> Response {
    let response = ui.interact(rect, id_salt.with("names"), Sense::click_and_drag());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, ui.visuals().extreme_bg_color);

    let rows = doc.rows();
    let selected = doc.selection.active.then(|| doc.selection.rows(rows));
    let font = FontId::new((cell_h * 0.62).clamp(8.0, 14.0), FontFamily::Proportional);

    for row in grid.visible_rows.clone() {
        let Some(seq) = doc.alignment.sequences.get(row) else { continue };
        let y = rect.min.y + row as f32 * cell_h - grid.offset.y;
        let line = Rect::from_min_max(Pos2::new(rect.min.x, y), Pos2::new(rect.max.x, y + cell_h));
        if selected.as_ref().is_some_and(|r| r.contains(&row)) {
            painter.rect_filled(line, 0.0, ui.visuals().selection.bg_fill.gamma_multiply(0.3));
        }
        if row == doc.selection.cursor.row {
            painter.rect_filled(line, 0.0, ui.visuals().selection.bg_fill.gamma_multiply(0.12));
        }
        let ink =
            if seq.hidden { ui.visuals().weak_text_color() } else { ui.visuals().text_color() };
        // Clip so long names cannot spill into the residue grid.
        painter.with_clip_rect(line.intersect(rect)).text(
            Pos2::new(rect.min.x + 6.0, line.center().y),
            Align2::LEFT_CENTER,
            &seq.id,
            font.clone(),
            ink,
        );
    }
    painter.line_segment(
        [Pos2::new(rect.max.x - 0.5, rect.min.y), Pos2::new(rect.max.x - 0.5, rect.max.y)],
        ui.visuals().widgets.noninteractive.bg_stroke,
    );
    response
}

fn paint_corner(ui: &egui::Ui, rect: Rect, _dark: bool) {
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, ui.visuals().extreme_bg_color);
    painter.line_segment(
        [Pos2::new(rect.min.x, rect.max.y - 0.5), Pos2::new(rect.max.x, rect.max.y - 0.5)],
        ui.visuals().widgets.noninteractive.bg_stroke,
    );
    painter.line_segment(
        [Pos2::new(rect.max.x - 0.5, rect.min.y), Pos2::new(rect.max.x - 0.5, rect.max.y)],
        ui.visuals().widgets.noninteractive.bg_stroke,
    );
}

fn quality_color(q: f32, dark: bool) -> Color32 {
    let c = if q > 0.9 {
        Color32::from_rgb(0x5E, 0xA9, 0x6F)
    } else if q > 0.6 {
        Color32::from_rgb(0xC8, 0xA8, 0x4A)
    } else {
        Color32::from_rgb(0xC4, 0x6A, 0x5A)
    };
    if dark {
        c.gamma_multiply(0.85)
    } else {
        c
    }
}

/// Choose a ruler tick interval so labels stay about 60 points apart and land
/// on round numbers.
fn tick_step(cell_w: f32) -> usize {
    let wanted = (60.0 / cell_w.max(0.1)).ceil() as usize;
    for &step in &[1usize, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 5000, 10000] {
        if step >= wanted {
            return step;
        }
    }
    100_000
}

fn cell_at(
    pos: Pos2,
    grid: &GridPaint,
    cell_w: f32,
    cell_h: f32,
    rows: usize,
    cols: usize,
) -> Cell {
    let col = ((pos.x - grid.origin.x) / cell_w).floor().max(0.0) as usize;
    let row = ((pos.y - grid.origin.y) / cell_h).floor().max(0.0) as usize;
    Cell::new(row.min(rows.saturating_sub(1)), col.min(cols.saturating_sub(1)))
}

fn handle_grid_input(
    doc: &mut Document,
    response: &Response,
    grid: &GridPaint,
    cell_w: f32,
    cell_h: f32,
    _actions: &mut [CanvasAction],
) {
    let rows = doc.rows();
    let cols = doc.width();
    if rows == 0 || cols == 0 {
        return;
    }
    if response.hovered() {
        response.ctx.set_cursor_icon(CursorIcon::Text);
    }
    if let Some(pos) = response.interact_pointer_pos() {
        let cell = cell_at(pos, grid, cell_w, cell_h, rows, cols);
        let shift = response.ctx.input(|i| i.modifiers.shift);
        if response.drag_started() || (response.clicked() && !shift) {
            doc.selection.place(cell, SelectionMode::Cells);
        } else if response.dragged() || (response.clicked() && shift) {
            doc.selection.extend_to(cell);
        }
    }
}

fn handle_header_input(
    doc: &mut Document,
    response: &Response,
    rect: Rect,
    grid: &GridPaint,
    cell_w: f32,
) {
    let cols = doc.width();
    let rows = doc.rows();
    if cols == 0 || rows == 0 {
        return;
    }
    if let Some(pos) = response.interact_pointer_pos() {
        let col = (((pos.x - rect.min.x + grid.offset.x) / cell_w).floor().max(0.0) as usize)
            .min(cols - 1);
        let cell = Cell::new(doc.selection.cursor.row.min(rows - 1), col);
        if response.drag_started() || response.clicked() {
            doc.selection.place(cell, SelectionMode::Columns);
            doc.selection.extend_to(cell);
        } else if response.dragged() {
            doc.selection.extend_to(cell);
        }
    }
}

fn handle_names_input(
    doc: &mut Document,
    response: &Response,
    rect: Rect,
    grid: &GridPaint,
    cell_h: f32,
    actions: &mut Vec<CanvasAction>,
) {
    let rows = doc.rows();
    let cols = doc.width();
    if rows == 0 || cols == 0 {
        return;
    }
    if let Some(pos) = response.interact_pointer_pos() {
        let row = (((pos.y - rect.min.y + grid.offset.y) / cell_h).floor().max(0.0) as usize)
            .min(rows - 1);
        let cell = Cell::new(row, doc.selection.cursor.col.min(cols - 1));
        if response.double_clicked() {
            actions.push(CanvasAction::BeginRename(row));
        } else if response.drag_started() || response.clicked() {
            doc.selection.place(cell, SelectionMode::Rows);
            doc.selection.extend_to(cell);
        } else if response.dragged() {
            doc.selection.extend_to(cell);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tick_steps_are_round_and_grow_as_you_zoom_out() {
        assert_eq!(tick_step(60.0), 1);
        assert_eq!(tick_step(12.0), 5);
        assert!(tick_step(1.0) >= 50);
        assert!(tick_step(0.2) > tick_step(2.0));
    }

    #[test]
    fn visible_range_covers_the_viewport_with_a_margin() {
        let r = visible_range(100.0, 200.0, 10.0, 1000);
        assert!(r.start <= 10 && r.end >= 20);
        assert!(r.end <= 1000);
    }

    #[test]
    fn visible_range_is_clamped_to_the_alignment() {
        assert_eq!(visible_range(0.0, 1000.0, 10.0, 5), 0..5);
        assert_eq!(visible_range(0.0, 100.0, 10.0, 0), 0..0);
    }

    #[test]
    fn visible_range_handles_a_scrolled_viewport_past_the_end() {
        let r = visible_range(9990.0, 10090.0, 10.0, 1000);
        assert_eq!(r, 999..1000);
    }
}
