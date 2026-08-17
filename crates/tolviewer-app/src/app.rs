//! The application: documents, menus, keyboard editing and background jobs.

use std::path::{Path, PathBuf};

use egui::{Key, Modifiers};
use tolviewer_align::{AlignParams, Engine, MatrixChoice, TreeMethod};
use tolviewer_clean::{GapPolicy, GblocksParams, GblocksResult};
use tolviewer_core::{Alignment, Alphabet, EditOp, Error, Result, Sequence, GAP};
use tolviewer_io::{Format, WriteOptions};

use crate::canvas::{AlignmentCanvas, CanvasAction, ViewSettings, MAX_ZOOM, MIN_ZOOM};
use crate::document::Document;
use crate::selection::SelectionMode;
use crate::tasks::{self, Task, TaskOutcome};
use crate::theme::ColorScheme;

const MAX_RECENT: usize = 12;

/// A transient message in the status bar.
struct Notice {
    text: String,
    is_error: bool,
    /// Seconds remaining before it fades.
    ttl: f32,
}

#[derive(Default)]
pub(crate) struct Dialogs {
    align: bool,
    clean: bool,
    export: bool,
    about: bool,
    goto: bool,
    rename: Option<(usize, String)>,
    /// Documents the user must decide about before quitting.
    confirm_close: Option<usize>,
}

/// Settings persisted between runs.
#[derive(serde::Serialize, serde::Deserialize)]
#[serde(default)]
struct Persisted {
    view: ViewSettings,
    recent: Vec<PathBuf>,
    export_format: String,
    write: PersistedWriteOptions,
}

impl Default for Persisted {
    fn default() -> Self {
        Persisted {
            view: ViewSettings::default(),
            recent: Vec::new(),
            export_format: "FASTA".to_string(),
            write: PersistedWriteOptions::default(),
        }
    }
}

/// `WriteOptions` is not serialisable (it lives in another crate), so the few
/// fields worth remembering are mirrored here.
#[derive(serde::Serialize, serde::Deserialize)]
#[serde(default)]
struct PersistedWriteOptions {
    line_width: usize,
    interleaved: bool,
    block_width: usize,
    uppercase: bool,
    strict_phylip_names: bool,
    include_hidden: bool,
}

impl Default for PersistedWriteOptions {
    fn default() -> Self {
        PersistedWriteOptions {
            line_width: 60,
            interleaved: false,
            block_width: 60,
            uppercase: false,
            strict_phylip_names: false,
            include_hidden: false,
        }
    }
}

impl PersistedWriteOptions {
    fn to_options(&self) -> WriteOptions {
        WriteOptions {
            line_width: self.line_width,
            interleaved: self.interleaved,
            block_width: self.block_width,
            uppercase: self.uppercase,
            strict_phylip_names: self.strict_phylip_names,
            include_hidden: self.include_hidden,
            ..WriteOptions::default()
        }
    }
}

pub struct TolViewerApp {
    docs: Vec<Document>,
    current: usize,
    view: ViewSettings,
    recent: Vec<PathBuf>,
    tasks: Vec<Task>,
    dialogs: Dialogs,
    notices: Vec<Notice>,
    align_params: AlignParams,
    clean_params: Option<GblocksParams>,
    /// Result of the last cleaning run, waiting for the user to apply it.
    pending_clean: Option<GblocksResult>,
    export_format: Format,
    write_options: PersistedWriteOptions,
    goto_input: String,
    /// Set when the user has confirmed they want to quit despite unsaved work.
    allow_exit: bool,
    /// Set while a quit is waiting on the unsaved-changes dialog.
    quit_pending: bool,
}

impl TolViewerApp {
    pub fn new(cc: &eframe::CreationContext<'_>, paths: Vec<PathBuf>) -> Self {
        let persisted: Persisted =
            cc.storage.and_then(|s| eframe::get_value(s, eframe::APP_KEY)).unwrap_or_default();

        let export_format = Format::all()
            .iter()
            .copied()
            .find(|f| f.name() == persisted.export_format && f.can_write())
            .unwrap_or(Format::Fasta);

        let mut app = TolViewerApp {
            docs: Vec::new(),
            current: 0,
            view: persisted.view,
            recent: persisted.recent,
            tasks: Vec::new(),
            dialogs: Dialogs::default(),
            notices: Vec::new(),
            align_params: AlignParams::default(),
            clean_params: None,
            pending_clean: None,
            export_format,
            write_options: persisted.write,
            goto_input: String::new(),
            allow_exit: false,
            quit_pending: false,
        };
        for path in paths {
            app.open_path(&path);
        }
        app
    }

    // ---- documents -----------------------------------------------------

    fn doc(&self) -> Option<&Document> {
        self.docs.get(self.current)
    }

    fn doc_mut(&mut self) -> Option<&mut Document> {
        self.docs.get_mut(self.current)
    }

    fn open_path(&mut self, path: &Path) {
        match tolviewer_io::read_file(path) {
            Ok(alignment) => {
                let format = Format::from_path(path).unwrap_or(Format::Fasta);
                let rows = alignment.len();
                let width = alignment.width();
                self.docs.push(Document::new(alignment, Some(path.to_path_buf()), format));
                self.current = self.docs.len() - 1;
                self.remember(path);
                self.info(format!(
                    "opened {} ({rows} sequences, {width} columns)",
                    path.file_name().unwrap_or_default().to_string_lossy()
                ));
            }
            Err(e) => self.error(format!("could not open {}: {e}", path.display())),
        }
    }

    fn remember(&mut self, path: &Path) {
        self.recent.retain(|p| p != path);
        self.recent.insert(0, path.to_path_buf());
        self.recent.truncate(MAX_RECENT);
    }

    fn open_dialog(&mut self) {
        let mut dialog = rfd::FileDialog::new().set_title("Open sequences or alignment");
        for format in Format::all().iter().filter(|f| f.can_read()) {
            dialog = dialog.add_filter(format.name(), format.extensions());
        }
        dialog = dialog.add_filter("All files", &["*"]);
        if let Some(paths) = dialog.pick_files() {
            for path in paths {
                self.open_path(&path);
            }
        }
    }

    fn save(&mut self, save_as: bool) {
        let Some(doc) = self.doc() else { return };
        let format = doc.format;
        let existing = doc.path.clone();
        // Plain "Save" on a document that already has a path writes straight
        // back to it; everything else asks.
        let path = if let (false, Some(path)) = (save_as, existing.clone()) {
            path
        } else {
            let mut dialog = rfd::FileDialog::new().set_title("Save alignment");
            if let Some(p) = &existing {
                if let Some(dir) = p.parent() {
                    dialog = dialog.set_directory(dir);
                }
                if let Some(name) = p.file_name() {
                    dialog = dialog.set_file_name(name.to_string_lossy());
                }
            } else {
                dialog = dialog.set_file_name(format!("alignment.{}", format.extensions()[0]));
            }
            for f in Format::all().iter().filter(|f| f.can_write()) {
                dialog = dialog.add_filter(f.name(), f.extensions());
            }
            match dialog.save_file() {
                Some(p) => p,
                None => return,
            }
        };
        // A "Save as" to a different extension should honour the extension.
        let format = Format::from_path(&path).filter(|f| f.can_write()).unwrap_or(format);
        let options = self.write_options.to_options();
        let alignment = self.doc().expect("checked above").alignment.clone();
        match write_padding_if_needed(&alignment, &path, format, &options) {
            Ok(padded) => {
                let doc = self.doc_mut().expect("checked above");
                doc.path = Some(path.clone());
                doc.format = format;
                doc.mark_saved();
                self.remember(&path);
                self.info(format!("saved {}", path.display()));
                if padded {
                    self.info(format!(
                        "{} needs equal-length rows, so short rows were padded with gaps in the \
                         saved file (the document itself is unchanged)",
                        format.name()
                    ));
                }
            }
            Err(e) => self.error(format!("could not save: {e}")),
        }
    }

    fn close_current(&mut self, force: bool) {
        let Some(doc) = self.doc() else { return };
        if doc.is_dirty() && !force {
            self.dialogs.confirm_close = Some(self.current);
            return;
        }
        let closing = self.current;
        self.docs.remove(closing);
        self.tasks.retain(|t| t.doc != closing);
        for task in &mut self.tasks {
            if task.doc > closing {
                task.doc -= 1;
            }
        }
        if self.current >= self.docs.len() {
            self.current = self.docs.len().saturating_sub(1);
        }
    }

    // ---- notices -------------------------------------------------------

    fn info(&mut self, text: impl Into<String>) {
        self.notices.push(Notice { text: text.into(), is_error: false, ttl: 6.0 });
    }

    fn error(&mut self, text: impl Into<String>) {
        self.notices.push(Notice { text: text.into(), is_error: true, ttl: 14.0 });
    }

    /// Run an edit, reporting failure instead of unwinding.
    fn edit(&mut self, op: EditOp) {
        let result = match self.doc_mut() {
            Some(doc) => doc.apply(op),
            None => return,
        };
        if let Err(e) = result {
            self.error(e.to_string());
        }
    }

    fn report(&mut self, result: Result<()>) {
        if let Err(e) = result {
            self.error(e.to_string());
        }
    }

    // ---- editing commands ----------------------------------------------

    fn delete_selection(&mut self) {
        let Some(doc) = self.doc() else { return };
        if !doc.selection.active {
            let cell = doc.selection.cursor;
            if doc.rows() > 0 && doc.width() > 0 {
                self.edit(EditOp::DeleteAt { row: cell.row, col: cell.col });
            }
            return;
        }
        match doc.selection.mode {
            SelectionMode::Columns => {
                let cols = doc.selection.cols(doc.width());
                self.edit(EditOp::DeleteColumns { start: cols.start, end: cols.end });
                if let Some(doc) = self.doc_mut() {
                    doc.selection.clear();
                }
            }
            SelectionMode::Rows => {
                let rows: Vec<usize> = doc.selection.rows(doc.rows()).collect();
                // Remove from the bottom up so earlier indices stay valid.
                for &row in rows.iter().rev() {
                    self.edit(EditOp::RemoveSequence { row });
                }
                if let Some(doc) = self.doc_mut() {
                    doc.selection.clear();
                }
            }
            SelectionMode::Cells => {
                // Blanking to gaps keeps every row the same length, which is
                // what you want inside an alignment.
                let rows = doc.selection.rows(doc.rows());
                let cols = doc.selection.cols(doc.width());
                let block = vec![vec![GAP; cols.len()]; rows.len()];
                self.edit(EditOp::SetBlock { row: rows.start, col: cols.start, residues: block });
            }
        }
    }

    /// Type a residue at the caret and step right.
    fn type_residue(&mut self, residue: u8) {
        let Some(doc) = self.doc() else { return };
        if doc.rows() == 0 || doc.width() == 0 {
            return;
        }
        let cell = doc.selection.cursor;
        let (rows, cols) = (doc.rows(), doc.width());
        self.edit(EditOp::SetResidue { row: cell.row, col: cell.col, residue });
        if let Some(doc) = self.doc_mut() {
            doc.selection.move_caret(0, 1, rows, cols, false);
        }
    }

    fn copy_selection_as_fasta(&mut self, ctx: &egui::Context) {
        let Some(doc) = self.doc() else { return };
        let rows = doc.target_rows();
        let cols = doc.target_cols();
        let subset = doc.alignment.subset(&rows, cols);
        match tolviewer_io::write_string(&subset, Format::Fasta, &WriteOptions::default()) {
            Ok(text) => {
                ctx.copy_text(text);
                self.info(format!("copied {} sequence(s) as FASTA", subset.len()));
            }
            Err(e) => self.error(e.to_string()),
        }
    }

    fn paste_fasta(&mut self, text: &str) {
        match tolviewer_io::parse(text.as_bytes(), Format::Fasta, "pasted") {
            Ok(parsed) if !parsed.is_empty() => {
                let n = parsed.len();
                if self.docs.is_empty() {
                    self.docs.push(Document::new(parsed, None, Format::Fasta));
                    self.current = 0;
                } else {
                    let at = self.doc().map_or(0, |d| d.rows());
                    for (i, seq) in parsed.sequences.into_iter().enumerate() {
                        self.edit(EditOp::InsertSequence { row: at + i, seq: Box::new(seq) });
                    }
                }
                self.info(format!("pasted {n} sequence(s)"));
            }
            Ok(_) => self.error("the clipboard held no sequences"),
            Err(e) => self.error(format!("the clipboard is not FASTA: {e}")),
        }
    }

    fn reverse_complement_selection(&mut self) {
        let Some(doc) = self.doc() else { return };
        let alphabet = doc.alphabet();
        if !alphabet.is_nucleotide() {
            self.error("reverse complement only applies to nucleotide sequences");
            return;
        }
        let rows = doc.target_rows();
        let mut next = doc.alignment.clone();
        for &row in &rows {
            if let Some(seq) = next.sequences.get_mut(row) {
                seq.reverse_complement(alphabet);
            }
        }
        let n = rows.len();
        let result = self.doc_mut().expect("checked above").replace("reverse complement", next);
        self.report(result);
        self.info(format!("reverse complemented {n} sequence(s)"));
    }

    fn sort_sequences(&mut self, by_length: bool) {
        let Some(doc) = self.doc() else { return };
        let mut next = doc.alignment.clone();
        if by_length {
            next.sequences.sort_by_key(|s| std::cmp::Reverse(s.ungapped_len()));
        } else {
            next.sequences.sort_by(|a, b| natural_cmp(&a.id, &b.id));
        }
        let result = self.doc_mut().expect("checked above").replace("sort sequences", next);
        self.report(result);
    }

    fn degap(&mut self) {
        let Some(doc) = self.doc() else { return };
        let mut next = doc.alignment.clone();
        next.degap();
        let result = self.doc_mut().expect("checked above").replace("remove all gaps", next);
        self.report(result);
    }

    fn remove_empty_columns(&mut self) {
        let Some(doc) = self.doc() else { return };
        let mut next = doc.alignment.clone();
        let removed = next.remove_all_gap_columns();
        if removed == 0 {
            self.info("no all-gap columns to remove");
            return;
        }
        let result = self.doc_mut().expect("checked above").replace("remove empty columns", next);
        self.report(result);
        self.info(format!("removed {removed} all-gap column(s)"));
    }

    fn apply_column_mask(&mut self, mask: Vec<bool>, label: &str) {
        let Some(doc) = self.doc() else { return };
        let mut next = doc.alignment.clone();
        match next.keep_columns(&mask) {
            Ok(removed) => {
                let result = self.doc_mut().expect("checked above").replace(label, next);
                self.report(result);
                self.info(format!("{label}: removed {removed} column(s)"));
            }
            Err(e) => self.error(e.to_string()),
        }
    }

    // ---- background jobs -----------------------------------------------

    fn start_align(&mut self, ctx: &egui::Context, selection_only: bool) {
        let Some(doc) = self.doc() else { return };
        if doc.rows() < 2 {
            self.error("aligning needs at least two sequences");
            return;
        }
        if self.tasks.iter().any(|t| t.doc == self.current && !t.is_finished()) {
            self.error("this alignment already has a job running");
            return;
        }
        let alignment = doc.alignment.clone();
        let params = self.align_params.clone();
        let task = if selection_only {
            let cols = doc.selection.cols(doc.width());
            if cols.len() < 2 {
                self.error("select at least two columns to realign");
                return;
            }
            tasks::realign_region(ctx, self.current, alignment, cols, params)
        } else {
            tasks::align(ctx, self.current, alignment, params)
        };
        self.tasks.push(task);
    }

    fn start_clean(&mut self, ctx: &egui::Context) {
        let Some(doc) = self.doc() else { return };
        if !doc.alignment.is_aligned() {
            self.error("cleaning needs an alignment; align the sequences first");
            return;
        }
        let params =
            self.clean_params.clone().unwrap_or_else(|| GblocksParams::defaults(doc.rows()));
        if let Err(e) = params.validate(doc.rows()) {
            self.error(e.to_string());
            return;
        }
        let alignment = doc.alignment.clone();
        self.tasks.push(tasks::clean(ctx, self.current, alignment, params));
    }

    fn poll_tasks(&mut self) {
        let mut finished = Vec::new();
        for (i, task) in self.tasks.iter_mut().enumerate() {
            if let Some(result) = task.poll() {
                finished.push((i, task.doc, result));
            }
        }
        for (_, doc_index, result) in finished {
            match result {
                Ok(TaskOutcome::Alignment { label, alignment }) => {
                    let width = alignment.width();
                    if let Some(doc) = self.docs.get_mut(doc_index) {
                        let r = doc.replace(&label, *alignment);
                        if let Err(e) = r {
                            self.error(e.to_string());
                        } else {
                            self.info(format!("{label}: {width} columns"));
                        }
                    }
                }
                Ok(TaskOutcome::Clean(result)) => {
                    let (kept, total) = (result.kept, result.total);
                    if let Some(doc) = self.docs.get_mut(doc_index) {
                        doc.set_clean_mask(result.mask.clone());
                    }
                    self.pending_clean = Some(*result);
                    self.dialogs.clean = true;
                    self.info(format!(
                        "cleaning would keep {kept} of {total} columns ({:.0}%)",
                        if total == 0 { 0.0 } else { 100.0 * kept as f32 / total as f32 }
                    ));
                }
                Err(Error::Cancelled) => self.info("cancelled"),
                Err(e) => self.error(e.to_string()),
            }
        }
        self.tasks.retain(|t| !t.is_finished());
    }

    // ---- keyboard ------------------------------------------------------

    /// `viewport_height` is the height of the canvas area, used to size a
    /// PageUp/PageDown step.
    fn handle_keys(&mut self, ctx: &egui::Context, viewport_height: f32) {
        if ctx.memory(|m| m.focused().is_some()) {
            // A text field has focus; leave the keys alone.
            return;
        }
        let mut typed: Vec<u8> = Vec::new();
        let events = ctx.input(|i| i.events.clone());
        for event in &events {
            if let egui::Event::Text(text) = event {
                for ch in text.chars() {
                    if ch.is_ascii_alphabetic() || matches!(ch, '-' | '.' | '?' | '*') {
                        typed.push(ch as u8);
                    }
                }
            }
        }

        let cmd = |key: Key| ctx.input_mut(|i| i.consume_key(Modifiers::COMMAND, key));
        let cmd_shift =
            |key: Key| ctx.input_mut(|i| i.consume_key(Modifiers::COMMAND | Modifiers::SHIFT, key));

        if cmd(Key::O) {
            self.open_dialog();
        }
        if cmd(Key::S) {
            self.save(false);
        }
        if cmd_shift(Key::S) {
            self.save(true);
        }
        if cmd(Key::W) {
            self.close_current(false);
        }
        if cmd(Key::A) {
            if let Some(doc) = self.doc_mut() {
                doc.select_all();
            }
        }
        if cmd(Key::C) {
            self.copy_selection_as_fasta(ctx);
        }
        if cmd(Key::V) {
            // egui delivers paste as an event; the clipboard is not readable
            // directly, so pick it up from the event stream.
            for event in &events {
                if let egui::Event::Paste(text) = event {
                    self.paste_fasta(text);
                }
            }
        }
        if cmd(Key::Z) {
            let r = self.doc_mut().map(|d| d.undo());
            match r {
                Some(Ok(Some(label))) => self.info(format!("undid {label}")),
                Some(Ok(None)) => self.info("nothing to undo"),
                Some(Err(e)) => self.error(e.to_string()),
                None => {}
            }
        }
        if cmd_shift(Key::Z) || cmd(Key::Y) {
            let r = self.doc_mut().map(|d| d.redo());
            match r {
                Some(Ok(Some(label))) => self.info(format!("redid {label}")),
                Some(Ok(None)) => self.info("nothing to redo"),
                Some(Err(e)) => self.error(e.to_string()),
                None => {}
            }
        }
        if cmd(Key::G) {
            self.dialogs.goto = true;
        }
        if cmd(Key::Plus) || cmd(Key::Equals) {
            self.zoom_by(1.25);
        }
        if cmd(Key::Minus) {
            self.zoom_by(0.8);
        }

        let Some(doc) = self.doc_mut() else { return };
        let (rows, cols) = (doc.rows(), doc.width());
        if rows == 0 || cols == 0 {
            return;
        }
        let shift = ctx.input(|i| i.modifiers.shift);
        let page = (viewport_height / (doc.zoom * 1.4)).max(1.0) as isize;

        let mut moves: Vec<(isize, isize)> = Vec::new();
        ctx.input(|i| {
            for (key, d) in [
                (Key::ArrowLeft, (0, -1)),
                (Key::ArrowRight, (0, 1)),
                (Key::ArrowUp, (-1, 0)),
                (Key::ArrowDown, (1, 0)),
                (Key::PageUp, (-page, 0)),
                (Key::PageDown, (page, 0)),
            ] {
                if i.key_pressed(key) {
                    moves.push(d);
                }
            }
            if i.key_pressed(Key::Home) {
                moves.push((0, -(cols as isize)));
            }
            if i.key_pressed(Key::End) {
                moves.push((0, cols as isize));
            }
        });
        for (dr, dc) in moves {
            doc.selection.move_caret(dr, dc, rows, cols, shift);
        }
        if ctx.input(|i| i.key_pressed(Key::Escape)) {
            doc.selection.clear();
        }

        let delete = ctx.input(|i| i.key_pressed(Key::Delete) || i.key_pressed(Key::Backspace));
        let backspace_only =
            ctx.input(|i| i.key_pressed(Key::Backspace) && !i.key_pressed(Key::Delete));
        let insert_gap = ctx.input(|i| i.key_pressed(Key::Insert));

        // Decide what to do while the document is borrowed, then act once the
        // borrow has ended, since the edit helpers need `&mut self`.
        let cell = doc.selection.cursor;
        let column_selection = doc.selection.active && doc.selection.mode == SelectionMode::Columns;
        let first_selected_col = doc.selection.cols(cols).start;
        let has_selection = doc.selection.active;
        if delete && backspace_only && !has_selection && cell.col > 0 {
            doc.selection.move_caret(0, -1, rows, cols, false);
        }

        if insert_gap {
            if column_selection {
                self.edit(EditOp::InsertColumns { at: first_selected_col, count: 1 });
            } else {
                self.edit(EditOp::InsertGap { row: cell.row, col: cell.col });
            }
            return;
        }
        if delete {
            self.delete_selection();
            return;
        }
        for residue in typed {
            self.type_residue(residue);
        }
    }

    fn zoom_by(&mut self, factor: f32) {
        if let Some(doc) = self.doc_mut() {
            doc.zoom = (doc.zoom * factor).clamp(MIN_ZOOM, MAX_ZOOM);
        }
    }
}

/// Write `alignment` to `path`, and if the format rejects it for having rows
/// of different lengths, pad a copy with trailing gaps and try once more.
///
/// PHYLIP, NEXUS, Clustal and Stockholm all require a rectangular matrix. A
/// user who has just typed into the last row should get their file, not a
/// refusal — but the document is left alone, and the caller says what happened.
/// Returns whether padding was needed.
fn write_padding_if_needed(
    alignment: &Alignment,
    path: &Path,
    format: Format,
    options: &WriteOptions,
) -> Result<bool> {
    match tolviewer_io::write_file(alignment, path, format, options) {
        Ok(()) => Ok(false),
        Err(e) => {
            // Only a formatting refusal is worth retrying; a disk error is not.
            if !matches!(e, Error::Format(_)) || alignment.is_aligned() {
                return Err(e);
            }
            let mut padded = alignment.clone();
            padded.pad_to_width();
            tolviewer_io::write_file(&padded, path, format, options)?;
            Ok(true)
        }
    }
}

/// Compare names so `seq2` sorts before `seq10`.
fn natural_cmp(a: &str, b: &str) -> std::cmp::Ordering {
    let mut ai = a.char_indices().peekable();
    let mut bi = b.char_indices().peekable();
    loop {
        match (ai.peek().copied(), bi.peek().copied()) {
            (None, None) => return std::cmp::Ordering::Equal,
            (None, Some(_)) => return std::cmp::Ordering::Less,
            (Some(_), None) => return std::cmp::Ordering::Greater,
            (Some((ax, ac)), Some((bx, bc))) => {
                if ac.is_ascii_digit() && bc.is_ascii_digit() {
                    let an = take_number(a, ax, &mut ai);
                    let bn = take_number(b, bx, &mut bi);
                    match an.cmp(&bn) {
                        std::cmp::Ordering::Equal => continue,
                        other => return other,
                    }
                }
                let ord = ac.to_ascii_lowercase().cmp(&bc.to_ascii_lowercase());
                if ord != std::cmp::Ordering::Equal {
                    return ord;
                }
                ai.next();
                bi.next();
            }
        }
    }
}

fn take_number(
    s: &str,
    start: usize,
    it: &mut std::iter::Peekable<std::str::CharIndices<'_>>,
) -> u128 {
    let mut end = start;
    while let Some(&(i, c)) = it.peek() {
        if c.is_ascii_digit() {
            end = i + c.len_utf8();
            it.next();
        } else {
            break;
        }
    }
    // Very long digit runs cannot be a real accession number; fall back to 0
    // so they compare by the surrounding text instead of panicking.
    s[start..end].parse().unwrap_or(0)
}

impl eframe::App for TolViewerApp {
    fn save(&mut self, storage: &mut dyn eframe::Storage) {
        let persisted = Persisted {
            view: self.view.clone(),
            recent: self.recent.clone(),
            export_format: self.export_format.name().to_string(),
            write: PersistedWriteOptions {
                line_width: self.write_options.line_width,
                interleaved: self.write_options.interleaved,
                block_width: self.write_options.block_width,
                uppercase: self.write_options.uppercase,
                strict_phylip_names: self.write_options.strict_phylip_names,
                include_hidden: self.write_options.include_hidden,
            },
        };
        eframe::set_value(storage, eframe::APP_KEY, &persisted);
    }

    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();
        let ctx = &ctx;
        self.poll_tasks();
        self.take_dropped_files(ctx);
        self.handle_keys(ctx, ui.max_rect().height());

        crate::ui::menu_bar(self, ui);
        crate::ui::tab_bar(self, ui);
        crate::ui::status_bar(self, ui);
        crate::ui::side_panel(self, ui);
        crate::ui::dialogs(self, ctx);

        let panel_fill = ui.visuals().panel_fill;
        egui::CentralPanel::no_frame().frame(egui::Frame::NONE.fill(panel_fill)).show(ui, |ui| {
            if self.docs.is_empty() {
                crate::ui::welcome(self, ui);
                return;
            }
            let mut actions = Vec::new();
            let current = self.current;
            let view = self.view.clone();
            if let Some(doc) = self.docs.get_mut(current) {
                AlignmentCanvas { doc, view: &view, id_salt: egui::Id::new(("canvas", current)) }
                    .show(ui, &mut actions);
            }
            for action in actions {
                match action {
                    CanvasAction::Edit(op) => self.edit(op),
                    CanvasAction::BeginRename(row) => {
                        let name = self
                            .doc()
                            .and_then(|d| d.alignment.sequences.get(row))
                            .map(|s| s.header())
                            .unwrap_or_default();
                        self.dialogs.rename = Some((row, name));
                    }
                    CanvasAction::ScrollToCaret => {}
                }
            }
        });

        // Age out notices.
        let dt = ctx.input(|i| i.stable_dt).min(0.1);
        for notice in &mut self.notices {
            notice.ttl -= dt;
        }
        self.notices.retain(|n| n.ttl > 0.0);
        if !self.tasks.is_empty() {
            ctx.request_repaint();
        }

        if ctx.input(|i| i.viewport().close_requested()) && !self.allow_exit {
            if let Some(index) = self.docs.iter().position(|d| d.is_dirty()) {
                ctx.send_viewport_cmd(egui::ViewportCommand::CancelClose);
                self.current = index;
                self.quit_pending = true;
                self.dialogs.confirm_close = Some(index);
            }
        }
    }
}

impl TolViewerApp {
    fn take_dropped_files(&mut self, ctx: &egui::Context) {
        let dropped: Vec<PathBuf> =
            ctx.input(|i| i.raw.dropped_files.iter().map(|f| f.path().to_path_buf()).collect());
        for path in dropped {
            self.open_path(&path);
        }
    }
}

// The UI module needs read/write access to a good deal of the app state; these
// accessors keep the fields private without making `ui.rs` unwieldy.
impl TolViewerApp {
    pub(crate) fn documents(&self) -> &[Document] {
        &self.docs
    }
    pub(crate) fn current_index(&self) -> usize {
        self.current
    }
    pub(crate) fn set_current(&mut self, index: usize) {
        if index < self.docs.len() {
            self.current = index;
        }
    }
    pub(crate) fn current_doc(&self) -> Option<&Document> {
        self.doc()
    }
    pub(crate) fn current_doc_mut(&mut self) -> Option<&mut Document> {
        self.doc_mut()
    }
    pub(crate) fn view_settings(&mut self) -> &mut ViewSettings {
        &mut self.view
    }
    pub(crate) fn recent_files(&self) -> &[PathBuf] {
        &self.recent
    }
    pub(crate) fn clear_recent(&mut self) {
        self.recent.clear();
    }
    pub(crate) fn running_tasks(&self) -> &[Task] {
        &self.tasks
    }
    pub(crate) fn notices_iter(&self) -> impl Iterator<Item = (&str, bool)> {
        self.notices.iter().map(|n| (n.text.as_str(), n.is_error))
    }
    pub(crate) fn dialogs_mut(&mut self) -> &mut Dialogs {
        &mut self.dialogs
    }
}

// Commands invoked from menus, kept here so `ui.rs` stays declarative.
impl TolViewerApp {
    pub(crate) fn cmd_open(&mut self) {
        self.open_dialog();
    }
    pub(crate) fn cmd_open_path(&mut self, path: PathBuf) {
        self.open_path(&path);
    }
    pub(crate) fn cmd_save(&mut self, save_as: bool) {
        self.save(save_as);
    }
    pub(crate) fn cmd_close(&mut self) {
        self.close_current(false);
    }
    pub(crate) fn cmd_force_close(&mut self, index: usize) {
        self.current = index;
        self.close_current(true);
    }
    pub(crate) fn cmd_undo(&mut self) {
        if let Some(doc) = self.doc_mut() {
            let r = doc.undo();
            self.report(r.map(|_| ()));
        }
    }
    pub(crate) fn cmd_redo(&mut self) {
        if let Some(doc) = self.doc_mut() {
            let r = doc.redo();
            self.report(r.map(|_| ()));
        }
    }
    pub(crate) fn cmd_select_all(&mut self) {
        if let Some(doc) = self.doc_mut() {
            doc.select_all();
        }
    }
    pub(crate) fn cmd_copy(&mut self, ctx: &egui::Context) {
        self.copy_selection_as_fasta(ctx);
    }
    pub(crate) fn cmd_delete_selection(&mut self) {
        self.delete_selection();
    }
    pub(crate) fn cmd_reverse_complement(&mut self) {
        self.reverse_complement_selection();
    }
    pub(crate) fn cmd_sort(&mut self, by_length: bool) {
        self.sort_sequences(by_length);
    }
    pub(crate) fn cmd_degap(&mut self) {
        self.degap();
    }
    pub(crate) fn cmd_remove_empty_columns(&mut self) {
        self.remove_empty_columns();
    }
    pub(crate) fn cmd_apply_mask(&mut self, mask: Vec<bool>, label: &str) {
        self.apply_column_mask(mask, label);
    }
    pub(crate) fn cmd_align(&mut self, ctx: &egui::Context, selection_only: bool) {
        self.start_align(ctx, selection_only);
    }
    pub(crate) fn cmd_clean(&mut self, ctx: &egui::Context) {
        self.start_clean(ctx);
    }
    pub(crate) fn cmd_zoom(&mut self, factor: f32) {
        self.zoom_by(factor);
    }
    pub(crate) fn cmd_set_alphabet(&mut self, alphabet: Alphabet) {
        if let Some(doc) = self.doc_mut() {
            doc.set_alphabet(alphabet);
        }
    }
    pub(crate) fn cmd_hide_selected(&mut self, hidden: bool) {
        let Some(doc) = self.doc() else { return };
        let rows = doc.target_rows();
        let mut next = doc.alignment.clone();
        for &row in &rows {
            if let Some(seq) = next.sequences.get_mut(row) {
                seq.hidden = hidden;
            }
        }
        let label = if hidden { "hide sequences" } else { "show sequences" };
        let r = self.doc_mut().expect("checked above").replace(label, next);
        self.report(r);
    }
    pub(crate) fn cmd_rename(&mut self, row: usize, header: &str) {
        let mut seq = Sequence::default();
        seq.set_header(header);
        self.edit(EditOp::Rename { row, id: seq.id, description: seq.description });
    }
    pub(crate) fn cmd_new_from_selection(&mut self) {
        let Some(doc) = self.doc() else { return };
        let rows = doc.target_rows();
        let cols = doc.target_cols();
        let mut subset = doc.alignment.subset(&rows, cols);
        subset.name = format!("{} (subset)", doc.alignment.name);
        let format = doc.format;
        self.docs.push(Document::new(subset, None, format));
        self.current = self.docs.len() - 1;
    }
    pub(crate) fn cmd_export(&mut self) {
        self.dialogs.export = true;
    }
    pub(crate) fn cmd_goto(&mut self, column: usize) {
        let Some(doc) = self.doc_mut() else { return };
        let cols = doc.width();
        if cols == 0 {
            return;
        }
        let col = column.saturating_sub(1).min(cols - 1);
        let row = doc.selection.cursor.row;
        doc.selection.place(crate::selection::Cell::new(row, col), SelectionMode::Cells);
        // Centre the caret; the scroll area picks this up on the next frame.
        doc.scroll_col = col as f32;
    }
    /// Called after the unsaved-changes dialog is answered. If a quit was
    /// waiting on it and nothing is dirty any more, let the quit through;
    /// otherwise move on to the next unsaved document.
    pub(crate) fn cmd_resume_quit(&mut self, ctx: &egui::Context) {
        if !self.quit_pending {
            return;
        }
        match self.docs.iter().position(|d| d.is_dirty()) {
            Some(next) => {
                self.current = next;
                self.dialogs.confirm_close = Some(next);
            }
            None => {
                self.quit_pending = false;
                self.allow_exit = true;
                ctx.send_viewport_cmd(egui::ViewportCommand::Close);
            }
        }
    }

    /// The user cancelled out of the dialog, so abandon the quit too.
    pub(crate) fn cmd_cancel_quit(&mut self) {
        self.quit_pending = false;
    }

    pub(crate) fn align_params_mut(&mut self) -> &mut AlignParams {
        &mut self.align_params
    }
    pub(crate) fn clean_params_mut(&mut self) -> &mut Option<GblocksParams> {
        &mut self.clean_params
    }
    pub(crate) fn pending_clean(&self) -> Option<&GblocksResult> {
        self.pending_clean.as_ref()
    }
    pub(crate) fn take_pending_clean(&mut self) -> Option<GblocksResult> {
        self.pending_clean.take()
    }
    pub(crate) fn export_format_mut(&mut self) -> &mut Format {
        &mut self.export_format
    }
    pub(crate) fn goto_input_mut(&mut self) -> &mut String {
        &mut self.goto_input
    }
    pub(crate) fn export_to(&mut self, path: PathBuf, format: Format, selection_only: bool) {
        let Some(doc) = self.doc() else { return };
        let alignment = if selection_only {
            doc.alignment.subset(&doc.target_rows(), doc.target_cols())
        } else {
            doc.alignment.clone()
        };
        let options = self.write_options.to_options();
        match write_padding_if_needed(&alignment, &path, format, &options) {
            Ok(padded) => {
                self.remember(&path);
                self.info(format!("exported {} to {}", format.name(), path.display()));
                if padded {
                    self.info(format!(
                        "{} needs equal-length rows, so short rows were padded with gaps",
                        format.name()
                    ));
                }
            }
            Err(e) => self.error(format!("export failed: {e}")),
        }
    }
}

/// Field access for `ui.rs`, which owns the dialog widgets.
impl Dialogs {
    pub(crate) fn align(&mut self) -> &mut bool {
        &mut self.align
    }
    pub(crate) fn clean(&mut self) -> &mut bool {
        &mut self.clean
    }
    pub(crate) fn export(&mut self) -> &mut bool {
        &mut self.export
    }
    pub(crate) fn about(&mut self) -> &mut bool {
        &mut self.about
    }
    pub(crate) fn goto(&mut self) -> &mut bool {
        &mut self.goto
    }
    pub(crate) fn rename(&mut self) -> &mut Option<(usize, String)> {
        &mut self.rename
    }
    pub(crate) fn confirm_close(&mut self) -> &mut Option<usize> {
        &mut self.confirm_close
    }
}

/// Default parameter sets offered in the align dialog.
pub(crate) fn engine_defaults(engine: Engine) -> AlignParams {
    AlignParams::for_engine(engine)
}

pub(crate) const MATRIX_CHOICES: &[(MatrixChoice, &str)] = &[
    (MatrixChoice::Auto, "Automatic"),
    (MatrixChoice::Blosum62, "BLOSUM62"),
    (MatrixChoice::Blosum45, "BLOSUM45"),
    (MatrixChoice::Blosum80, "BLOSUM80"),
    (MatrixChoice::Pam250, "PAM250"),
    (MatrixChoice::Iub, "IUB (DNA)"),
    (MatrixChoice::ClustalDna, "Clustal DNA"),
    (MatrixChoice::Identity, "Identity"),
];

pub(crate) const TREE_CHOICES: &[(TreeMethod, &str)] =
    &[(TreeMethod::NeighborJoining, "Neighbour joining"), (TreeMethod::Upgma, "UPGMA")];

pub(crate) const GAP_POLICIES: &[(GapPolicy, &str)] = &[
    (GapPolicy::None, "No gaps allowed"),
    (GapPolicy::Half, "Gaps in up to half"),
    (GapPolicy::All, "Gaps allowed"),
];

pub(crate) const SCHEMES: &[ColorScheme] = ColorScheme::ALL;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn natural_order_puts_seq2_before_seq10() {
        let mut names = vec!["seq10", "seq2", "seq1", "Seq3"];
        names.sort_by(|a, b| natural_cmp(a, b));
        assert_eq!(names, vec!["seq1", "seq2", "Seq3", "seq10"]);
    }

    #[test]
    fn natural_order_is_case_insensitive_but_total() {
        assert_eq!(natural_cmp("abc", "ABC"), std::cmp::Ordering::Equal);
        assert_eq!(natural_cmp("a", "ab"), std::cmp::Ordering::Less);
    }

    #[test]
    fn natural_order_survives_absurd_digit_runs() {
        let long = "x".to_string() + &"9".repeat(60);
        // Must not panic on integer overflow.
        let _ = natural_cmp(&long, "x1");
    }
}
