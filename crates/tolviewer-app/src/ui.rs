//! Menus, panels and dialogs.
//!
//! Everything here reads and writes `TolViewerApp` through its `cmd_*`
//! accessors so the widget code stays declarative and the state transitions
//! all live in `app.rs`.

use std::path::PathBuf;

use egui::{Align, Layout, RichText};
use tolviewer_align::Engine;
use tolviewer_clean::GblocksParams;
use tolviewer_core::Alphabet;
use tolviewer_io::Format;
use tolviewer_library::{EntryKind, NodeId, SaveChoice, SaveTarget};

use crate::app::{TolViewerApp, GAP_POLICIES, MATRIX_CHOICES, SCHEMES, TREE_CHOICES};
use crate::canvas::{MAX_ZOOM, MIN_ZOOM};
use crate::chromatogram::{Chromatogram, Link, TraceAction};

pub fn menu_bar(app: &mut TolViewerApp, ui: &mut egui::Ui) {
    let ctx = ui.ctx().clone();
    let ctx = &ctx;
    egui::Panel::top(egui::Id::new("menu")).show(ui, |ui| {
        egui::MenuBar::new().ui(ui, |ui| {
            let has_doc = app.current_doc().is_some();

            ui.menu_button("File", |ui| {
                if ui.button("Open…").clicked() {
                    app.cmd_open();
                    ui.close();
                }
                ui.menu_button("Open recent", |ui| {
                    let recent: Vec<_> = app.recent_files().to_vec();
                    if recent.is_empty() {
                        ui.label(RichText::new("nothing yet").italics());
                    }
                    for path in recent {
                        let label = path
                            .file_name()
                            .map(|s| s.to_string_lossy().into_owned())
                            .unwrap_or_else(|| path.display().to_string());
                        if ui.button(label).on_hover_text(path.display().to_string()).clicked() {
                            app.cmd_open_path(path);
                            ui.close();
                        }
                    }
                    ui.separator();
                    if ui.button("Clear list").clicked() {
                        app.clear_recent();
                        ui.close();
                    }
                });
                ui.separator();
                if ui.add_enabled(has_doc, egui::Button::new("Save")).clicked() {
                    app.cmd_save(false);
                    ui.close();
                }
                if ui.add_enabled(has_doc, egui::Button::new("Save as…")).clicked() {
                    app.cmd_save(true);
                    ui.close();
                }
                if ui.add_enabled(has_doc, egui::Button::new("Export…")).clicked() {
                    app.cmd_export();
                    ui.close();
                }
                ui.separator();
                if ui.add_enabled(has_doc, egui::Button::new("Close")).clicked() {
                    app.cmd_close();
                    ui.close();
                }
                if ui.button("Quit").clicked() {
                    ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                }
            });

            ui.menu_button("Edit", |ui| {
                let (undo_label, redo_label) = match app.current_doc() {
                    Some(doc) => (
                        doc.undo.undo_label().map(|s| format!("Undo {s}")),
                        doc.undo.redo_label().map(|s| format!("Redo {s}")),
                    ),
                    None => (None, None),
                };
                if ui
                    .add_enabled(
                        undo_label.is_some(),
                        egui::Button::new(undo_label.unwrap_or_else(|| "Undo".into())),
                    )
                    .clicked()
                {
                    app.cmd_undo();
                    ui.close();
                }
                if ui
                    .add_enabled(
                        redo_label.is_some(),
                        egui::Button::new(redo_label.unwrap_or_else(|| "Redo".into())),
                    )
                    .clicked()
                {
                    app.cmd_redo();
                    ui.close();
                }
                ui.separator();
                if ui.add_enabled(has_doc, egui::Button::new("Select all")).clicked() {
                    app.cmd_select_all();
                    ui.close();
                }
                if ui.add_enabled(has_doc, egui::Button::new("Copy as FASTA")).clicked() {
                    app.cmd_copy(ctx);
                    ui.close();
                }
                if ui
                    .add_enabled(has_doc, egui::Button::new("Delete selection"))
                    .on_hover_text(
                        "Columns and rows are removed; a cell selection is blanked to gaps",
                    )
                    .clicked()
                {
                    app.cmd_delete_selection();
                    ui.close();
                }
                ui.separator();
                if ui
                    .add_enabled(has_doc, egui::Button::new("New document from selection"))
                    .clicked()
                {
                    app.cmd_new_from_selection();
                    ui.close();
                }
                if ui.add_enabled(has_doc, egui::Button::new("Go to column…")).clicked() {
                    *app.dialogs_mut().goto() = true;
                    ui.close();
                }
            });

            ui.menu_button("Sequences", |ui| {
                if ui.add_enabled(has_doc, egui::Button::new("Sort by name")).clicked() {
                    app.cmd_sort(false);
                    ui.close();
                }
                if ui.add_enabled(has_doc, egui::Button::new("Sort by length")).clicked() {
                    app.cmd_sort(true);
                    ui.close();
                }
                ui.separator();
                if ui.add_enabled(has_doc, egui::Button::new("Reverse complement")).clicked() {
                    app.cmd_reverse_complement();
                    ui.close();
                }
                if ui.add_enabled(has_doc, egui::Button::new("Remove all gaps")).clicked() {
                    app.cmd_degap();
                    ui.close();
                }
                ui.separator();
                if ui.add_enabled(has_doc, egui::Button::new("Hide selected")).clicked() {
                    app.cmd_hide_selected(true);
                    ui.close();
                }
                if ui.add_enabled(has_doc, egui::Button::new("Show all")).clicked() {
                    app.cmd_hide_selected(false);
                    ui.close();
                }
                ui.separator();
                ui.menu_button("Sequence type", |ui| {
                    let current = app.current_doc().map(|d| d.alphabet());
                    for alphabet in [Alphabet::Dna, Alphabet::Rna, Alphabet::Protein] {
                        if ui.radio(current == Some(alphabet), alphabet.name()).clicked() {
                            app.cmd_set_alphabet(alphabet);
                            ui.close();
                        }
                    }
                });
            });

            ui.menu_button("Library", |ui| {
                let entries = app.library().selected_entries().len();
                if ui.button("New library").clicked() {
                    app.cmd_library_new();
                    ui.close();
                }
                if ui.button("Open library…").clicked() {
                    app.cmd_library_open();
                    ui.close();
                }
                if ui.button("Save library").clicked() {
                    app.cmd_library_save(false);
                    ui.close();
                }
                if ui.button("Save library as…").clicked() {
                    app.cmd_library_save(true);
                    ui.close();
                }
                ui.separator();
                if ui.button("Add files…").clicked() {
                    app.cmd_library_add_files();
                    ui.close();
                }
                if ui
                    .button("Add a folder of files…")
                    .on_hover_text("Everything readable under a directory, as one folder")
                    .clicked()
                {
                    app.cmd_library_add_folder_of_files();
                    ui.close();
                }
                if ui.button("New folder…").clicked() {
                    *app.dialogs_mut().new_folder() = Some(String::new());
                    ui.close();
                }
                ui.separator();
                if ui.add_enabled(entries >= 2, egui::Button::new("Align selected")).clicked() {
                    app.cmd_library_align_selection(ctx);
                    ui.close();
                }
                if ui
                    .add_enabled(entries >= 2, egui::Button::new("Concatenate selected…"))
                    .on_hover_text(
                        "Join per-locus alignments into one matrix, matching samples by name",
                    )
                    .clicked()
                {
                    app.cmd_library_concatenate();
                    ui.close();
                }
                if ui
                    .add_enabled(has_doc, egui::Button::new("Extract selected sequences"))
                    .on_hover_text("Put the rows selected in this alignment into the library")
                    .clicked()
                {
                    app.cmd_library_extract();
                    ui.close();
                }
                ui.separator();
                if ui.button("Primers…").clicked() {
                    *app.dialogs_mut().primers() = true;
                    ui.close();
                }
                if ui.add_enabled(entries >= 1, egui::Button::new("Map primers")).clicked() {
                    app.cmd_library_map_primers();
                    ui.close();
                }
                if ui.add_enabled(entries >= 1, egui::Button::new("Trim primers…")).clicked() {
                    *app.dialogs_mut().trim() = true;
                    ui.close();
                }
                ui.separator();
                let visible = app.library().visible;
                if ui.checkbox(&mut { visible }, "Show the library panel").clicked() {
                    app.library_mut().visible = !visible;
                    ui.close();
                }
            });

            ui.menu_button("Align", |ui| {
                for engine in Engine::all() {
                    if ui
                        .add_enabled(
                            has_doc,
                            egui::Button::new(format!("Align with {}", engine.name())),
                        )
                        .clicked()
                    {
                        app.align_params_mut().engine = *engine;
                        app.cmd_align(ctx, false);
                        ui.close();
                    }
                }
                ui.separator();
                if ui.add_enabled(has_doc, egui::Button::new("Realign selected columns")).clicked()
                {
                    app.cmd_align(ctx, true);
                    ui.close();
                }
                if ui.add_enabled(has_doc, egui::Button::new("Alignment settings…")).clicked() {
                    *app.dialogs_mut().align() = true;
                    ui.close();
                }
            });

            ui.menu_button("Clean", |ui| {
                if ui
                    .add_enabled(has_doc, egui::Button::new("Gblocks…"))
                    .on_hover_text("Select conserved blocks (Castresana 2000)")
                    .clicked()
                {
                    *app.dialogs_mut().clean() = true;
                    ui.close();
                }
                ui.separator();
                if ui.add_enabled(has_doc, egui::Button::new("Remove all-gap columns")).clicked() {
                    app.cmd_remove_empty_columns();
                    ui.close();
                }
                if ui
                    .add_enabled(has_doc, egui::Button::new("Remove columns over 50% gaps"))
                    .clicked()
                {
                    if let Some(doc) = app.current_doc() {
                        let mask = tolviewer_clean::remove_gappy_columns(&doc.alignment, 0.5);
                        app.cmd_apply_mask(mask, "remove gappy columns");
                    }
                    ui.close();
                }
                if ui.add_enabled(has_doc, egui::Button::new("Trim ragged ends")).clicked() {
                    if let Some(doc) = app.current_doc() {
                        let range = tolviewer_clean::trim_ends(&doc.alignment, 0.5);
                        let width = doc.width();
                        let mask: Vec<bool> = (0..width).map(|c| range.contains(&c)).collect();
                        app.cmd_apply_mask(mask, "trim ends");
                    }
                    ui.close();
                }
            });

            ui.menu_button("View", |ui| {
                ui.menu_button("Colour scheme", |ui| {
                    let current = app.view_settings().scheme;
                    for &scheme in SCHEMES {
                        if ui.radio(current == scheme, scheme.name()).clicked() {
                            app.view_settings().scheme = scheme;
                            ui.close();
                        }
                    }
                });
                let view = app.view_settings();
                ui.checkbox(&mut view.show_ruler, "Position ruler");
                ui.checkbox(&mut view.show_quality_track, "Quality track");
                ui.checkbox(&mut view.show_consensus, "Consensus row");
                ui.checkbox(&mut view.dot_matches, "Dots for consensus matches");
                ui.separator();
                if ui.button("Zoom in").clicked() {
                    app.cmd_zoom(1.25);
                }
                if ui.button("Zoom out").clicked() {
                    app.cmd_zoom(0.8);
                }
            });

            ui.menu_button("Help", |ui| {
                if ui.button("About TOLViewer").clicked() {
                    *app.dialogs_mut().about() = true;
                    ui.close();
                }
            });

            ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                if let Some(doc) = app.current_doc_mut() {
                    ui.add(
                        egui::Slider::new(&mut doc.zoom, MIN_ZOOM..=MAX_ZOOM)
                            .show_value(false)
                            .text("zoom"),
                    );
                }
            });
        });
    });
}

pub fn tab_bar(app: &mut TolViewerApp, ui: &mut egui::Ui) {
    if app.documents().len() < 2 {
        return;
    }
    egui::Panel::top(egui::Id::new("tabs")).show(ui, |ui| {
        ui.horizontal_wrapped(|ui| {
            let current = app.current_index();
            let titles: Vec<String> = app.documents().iter().map(|d| d.title()).collect();
            for (i, title) in titles.into_iter().enumerate() {
                if ui.selectable_label(i == current, title).clicked() {
                    app.set_current(i);
                }
            }
        });
    });
}

pub fn status_bar(app: &mut TolViewerApp, ui: &mut egui::Ui) {
    egui::Panel::bottom(egui::Id::new("status")).show(ui, |ui| {
        ui.horizontal(|ui| {
            if let Some(doc) = app.current_doc() {
                let rows = doc.rows();
                let cols = doc.width();
                let cur = doc.selection.cursor;
                ui.label(format!("{rows} seqs x {cols} cols"));
                ui.separator();
                ui.label(format!("row {} col {}", cur.row + 1, cur.col + 1));
                if let Some(n) = doc.caret_residue_number() {
                    ui.label(format!("(residue {n})"));
                }
                let selected = doc.selection.cell_count(rows, cols);
                if selected > 0 {
                    ui.separator();
                    let r = doc.selection.rows(rows).len();
                    let c = doc.selection.cols(cols).len();
                    ui.label(format!("selected {r} x {c}"));
                }
                if !doc.alignment.is_aligned() {
                    ui.separator();
                    ui.label(RichText::new("not aligned").color(ui.visuals().warn_fg_color))
                        .on_hover_text("rows have different lengths; column edits will pad them");
                }
            } else {
                ui.label("no document open");
            }

            ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                if let Some((text, is_error)) = app.notices_iter().last() {
                    let color = if is_error {
                        ui.visuals().error_fg_color
                    } else {
                        ui.visuals().weak_text_color()
                    };
                    ui.label(RichText::new(text).color(color));
                }
            });
        });

        for task in app.running_tasks() {
            ui.horizontal(|ui| {
                ui.add(
                    egui::ProgressBar::new(task.progress.fraction())
                        .desired_width(240.0)
                        .text(format!("{}: {}", task.label, task.progress.message())),
                );
                if ui.button("Cancel").clicked() {
                    task.progress.request_cancel();
                }
            });
        }
    });
}

pub fn side_panel(app: &mut TolViewerApp, ui: &mut egui::Ui) {
    if app.current_doc().is_none() {
        return;
    }
    egui::Panel::right(egui::Id::new("info")).resizable(true).default_size(230.0).show(ui, |ui| {
        ui.heading("Alignment");
        let Some(doc) = app.current_doc_mut() else { return };
        let rows = doc.rows();
        let cols = doc.width();
        let alphabet = doc.alphabet();
        egui::Grid::new("info-grid").num_columns(2).striped(true).show(ui, |ui| {
            ui.label("Sequences");
            ui.label(rows.to_string());
            ui.end_row();
            ui.label("Columns");
            ui.label(cols.to_string());
            ui.end_row();
            ui.label("Type");
            ui.label(alphabet.name());
            ui.end_row();
            let conserved = doc.consensus().conserved_fraction();
            ui.label("Identical cols");
            ui.label(format!("{:.1}%", conserved * 100.0));
            ui.end_row();
            let gappy = doc.alignment.all_gap_columns().len();
            ui.label("All-gap cols");
            ui.label(gappy.to_string());
            ui.end_row();
        });

        ui.separator();
        ui.heading("Sequences");
        let selected_rows = doc.selection.rows(rows);
        let active = doc.selection.active;
        egui::ScrollArea::vertical().show(ui, |ui| {
            for (i, seq) in doc.alignment.sequences.iter().enumerate() {
                let marked = active && selected_rows.contains(&i);
                let ambiguity = seq.ambiguity_fraction(alphabet);
                let mut text = RichText::new(format!("{}  ({} nt)", seq.id, seq.ungapped_len()));
                if seq.hidden {
                    text = text.weak().strikethrough();
                }
                if marked {
                    text = text.strong();
                }
                let response = ui.label(text);
                let mut hover = format!(
                    "{}\nungapped length {}\nambiguous {:.1}%",
                    seq.header(),
                    seq.ungapped_len(),
                    ambiguity * 100.0
                );
                if let Some(q) = seq.mean_quality() {
                    hover.push_str(&format!("\nmean Phred {q:.1}"));
                }
                response.on_hover_text(hover);
            }
        });
    });
}

/// The library tree down the left-hand side.
///
/// Folders and sequences are one list; clicking selects, ctrl-click adds,
/// shift-click takes a run, and a double-click opens. Commands act on the
/// selection, resolving a folder to everything under it, which is what makes
/// "select these six reads and align them" a single gesture.
pub fn library_panel(app: &mut TolViewerApp, ui: &mut egui::Ui) {
    if !app.library().visible {
        return;
    }
    let ctx = ui.ctx().clone();
    egui::Panel::left(egui::Id::new("library")).resizable(true).default_size(260.0).show(
        ui,
        |ui| {
            ui.horizontal(|ui| {
                ui.heading("Library");
                ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                    if ui.small_button("✖").on_hover_text("Hide the library panel").clicked() {
                        app.library_mut().visible = false;
                    }
                });
            });
            let name = app.library().library.name.clone();
            let dirty = app.library().library.is_dirty();
            ui.label(
                RichText::new(if dirty { format!("{name} *") } else { name.clone() }).strong(),
            )
            .on_hover_text(match &app.library().library.path {
                Some(p) => p.display().to_string(),
                None => "not saved yet".to_string(),
            });

            ui.horizontal_wrapped(|ui| {
                if ui.button("Add files…").clicked() {
                    app.cmd_library_add_files();
                }
                if ui.button("New folder").clicked() {
                    *app.dialogs_mut().new_folder() = Some(String::new());
                }
            });
            ui.separator();

            let mut actions: Vec<TreeAction> = Vec::new();
            egui::ScrollArea::both().auto_shrink([false, false]).show(ui, |ui| {
                if app.library().library.is_empty() {
                    ui.add_space(8.0);
                    ui.label(
                        RichText::new(
                            "Empty. Add the files for a project and arrange them into \
                             folders — one per gene, or per sequencing run.",
                        )
                        .weak(),
                    );
                    return;
                }
                let roots = app.library().library.roots().to_vec();
                for root in roots {
                    draw_node(app, ui, root, 0, &mut actions);
                }
            });

            for action in actions {
                match action {
                    TreeAction::Click { id, ctrl, shift } => {
                        let state = app.library_mut();
                        if ctrl {
                            state.toggle(id);
                        } else if shift {
                            state.extend_to(id);
                        } else {
                            state.select_only(id);
                        }
                    }
                    TreeAction::Open(id) => app.cmd_library_open_entry(id),
                    TreeAction::Toggle(id) => {
                        if let Some(folder) =
                            app.library_mut().library.get_mut(id).and_then(|n| n.folder_mut())
                        {
                            folder.expanded = !folder.expanded;
                        }
                    }
                    TreeAction::BeginRename(id) => {
                        let name = app
                            .library()
                            .library
                            .get(id)
                            .map(|n| n.name.clone())
                            .unwrap_or_default();
                        app.library_mut().renaming = Some((id, name));
                    }
                    TreeAction::Rename(id, name) => {
                        app.cmd_library_rename(id, &name);
                        app.library_mut().renaming = None;
                    }
                    TreeAction::Typing(id, text) => {
                        app.library_mut().renaming = Some((id, text));
                    }
                    TreeAction::Remove => app.cmd_library_remove(),
                    TreeAction::Reverse => app.cmd_library_reverse(),
                    TreeAction::Align => app.cmd_library_align_selection(&ctx),
                    TreeAction::OpenSelection => app.cmd_library_open_selection(),
                    TreeAction::Concatenate => app.cmd_library_concatenate(),
                    TreeAction::MapPrimers => app.cmd_library_map_primers(),
                    TreeAction::TrimPrimers => app.cmd_library_trim_primers(),
                }
            }

            ui.separator();
            let entries = app.library().selected_entries().len();
            ui.horizontal_wrapped(|ui| {
                if ui
                    .add_enabled(entries >= 2, egui::Button::new("Align selected"))
                    .on_hover_text("Gather the selected sequences into a new tab and align them")
                    .clicked()
                {
                    app.cmd_library_align_selection(&ctx);
                }
                if ui
                    .add_enabled(entries >= 1, egui::Button::new("Reverse"))
                    .on_hover_text(
                        "Show these sequences reverse complemented. The files are not changed.",
                    )
                    .clicked()
                {
                    app.cmd_library_reverse();
                }
            });
            ui.label(
                RichText::new(if entries == 0 {
                    "nothing selected".to_string()
                } else {
                    format!("{entries} sequence file(s) selected")
                })
                .weak()
                .small(),
            );
        },
    );
}

/// What a click in the tree meant. Collected while the tree is drawn and acted
/// on afterwards, because the handlers need `&mut app` and the tree is holding
/// it.
enum TreeAction {
    Click {
        id: NodeId,
        ctrl: bool,
        shift: bool,
    },
    /// The rename field's text as it stands, kept so it survives the frame.
    Typing(NodeId, String),
    Open(NodeId),
    Toggle(NodeId),
    BeginRename(NodeId),
    Rename(NodeId, String),
    Remove,
    Reverse,
    Align,
    OpenSelection,
    Concatenate,
    MapPrimers,
    TrimPrimers,
}

/// One row of the tree, copied out of the library so the widgets can borrow
/// `app` mutably.
struct Row {
    selected: bool,
    is_folder: bool,
    expanded: bool,
    name: String,
    children: Vec<NodeId>,
    /// The text being typed, when this row is the one being renamed.
    renaming: Option<String>,
    kind: Option<EntryKind>,
    /// The file this entry points at has gone. The row still shows — silently
    /// hiding it would leave the user wondering where their read went.
    missing: bool,
    reversed: bool,
    hover: Option<String>,
}

fn read_row(app: &TolViewerApp, id: NodeId) -> Option<Row> {
    let state = app.library();
    let node = state.library.get(id)?;
    let entry = node.entry();
    let missing = entry.is_some_and(|e| !e.source_path().exists());
    let reversed = entry.is_some_and(|e| e.reversed);
    let hover = entry.map(|e| {
        let mut text = format!(
            "{}\n{}\n{}",
            state.library.path_of(id),
            e.kind.name(),
            e.source_path().display()
        );
        if e.working.is_some() {
            text.push_str("\n\nEdits are being kept in this copy; the original is untouched.");
        }
        if e.select.is_some() {
            text.push_str("\n\nPart of a larger file, so edits must be saved separately.");
        }
        if reversed {
            text.push_str("\n\nShown reverse complemented.");
        }
        if missing {
            text.push_str("\n\nThis file is missing.");
        }
        if !e.note.is_empty() {
            text.push_str(&format!("\n\n{}", e.note));
        }
        text
    });
    Some(Row {
        selected: state.is_selected(id),
        is_folder: node.is_folder(),
        expanded: node.folder().is_some_and(|f| f.expanded),
        name: node.name.clone(),
        children: node.folder().map(|f| f.children.clone()).unwrap_or_default(),
        renaming: state.renaming.clone().filter(|(r, _)| *r == id).map(|(_, text)| text),
        kind: entry.map(|e| e.kind),
        missing,
        reversed,
        hover,
    })
}

fn draw_node(
    app: &mut TolViewerApp,
    ui: &mut egui::Ui,
    id: NodeId,
    depth: usize,
    actions: &mut Vec<TreeAction>,
) {
    // Everything the row needs is copied out first: the widgets below need
    // `&mut app` for the context menu, so nothing may still be borrowing it.
    let Row {
        selected,
        is_folder,
        expanded,
        name,
        children,
        renaming,
        kind,
        missing,
        reversed,
        hover,
    } = match read_row(app, id) {
        Some(row) => row,
        None => return,
    };

    ui.horizontal(|ui| {
        ui.add_space(depth as f32 * 12.0);
        if is_folder {
            let arrow = if expanded { "▾" } else { "▸" };
            if ui.add(egui::Label::new(arrow).sense(egui::Sense::click())).clicked() {
                actions.push(TreeAction::Toggle(id));
            }
        } else {
            ui.add_space(12.0);
        }

        if let Some(mut text) = renaming {
            let response = ui.add(egui::TextEdit::singleline(&mut text).desired_width(150.0));
            response.request_focus();
            if response.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                actions.push(TreeAction::Rename(id, text));
            } else if ui.input(|i| i.key_pressed(egui::Key::Escape)) {
                actions.push(TreeAction::Rename(id, name.clone()));
            } else {
                actions.push(TreeAction::Typing(id, text));
            }
            return;
        }

        let icon = match (is_folder, kind) {
            (true, _) => "📁",
            (_, Some(EntryKind::Trace)) => "📈",
            (_, Some(EntryKind::Alignment)) => "▤",
            _ => "≡",
        };
        let mut label = RichText::new(format!("{icon} {name}{}", if reversed { " ↩" } else { "" }));
        if selected {
            label = label.strong();
        }
        if missing {
            label = label.strikethrough().color(ui.visuals().error_fg_color);
        }
        let response = ui.add(egui::Label::new(label).sense(egui::Sense::click()));
        let response = match &hover {
            Some(text) => response.on_hover_text(text),
            None => response,
        };

        if response.clicked() {
            let (ctrl, shift) = ui.input(|i| (i.modifiers.command, i.modifiers.shift));
            actions.push(TreeAction::Click { id, ctrl, shift });
        }
        if response.double_clicked() {
            if is_folder {
                actions.push(TreeAction::Toggle(id));
            } else {
                actions.push(TreeAction::Click { id, ctrl: false, shift: false });
                actions.push(TreeAction::Open(id));
            }
        }
        response.context_menu(|ui| {
            // Right-clicking something outside the selection acts on it alone,
            // which is what every file manager does.
            if !selected {
                actions.push(TreeAction::Click { id, ctrl: false, shift: false });
            }
            if !is_folder && ui.button("Open").clicked() {
                actions.push(TreeAction::Open(id));
                ui.close();
            }
            if ui.button("Open selected together").clicked() {
                actions.push(TreeAction::OpenSelection);
                ui.close();
            }
            if ui.button("Align selected…").clicked() {
                actions.push(TreeAction::Align);
                ui.close();
            }
            if ui.button("Concatenate selected").clicked() {
                actions.push(TreeAction::Concatenate);
                ui.close();
            }
            ui.separator();
            if ui.button("Reverse complement").clicked() {
                actions.push(TreeAction::Reverse);
                ui.close();
            }
            if ui.button("Map primers").clicked() {
                actions.push(TreeAction::MapPrimers);
                ui.close();
            }
            if ui.button("Trim primers").clicked() {
                actions.push(TreeAction::TrimPrimers);
                ui.close();
            }
            ui.separator();
            if ui.button("Rename").clicked() {
                actions.push(TreeAction::BeginRename(id));
                ui.close();
            }
            if ui
                .button("Remove from library")
                .on_hover_text("Takes it out of the tree. The file itself is not deleted.")
                .clicked()
            {
                actions.push(TreeAction::Remove);
                ui.close();
            }
        });
    });

    if is_folder && expanded {
        for child in children {
            draw_node(app, ui, child, depth + 1, actions);
        }
    }
}

/// The residue at ungapped index `index` of a row, for showing which base is
/// currently called.
fn current_call(residues: &[u8], index: usize) -> Option<u8> {
    residues.iter().copied().filter(|&c| !tolviewer_core::is_gap(c)).nth(index)
}

/// The chromatogram, under the alignment, for documents opened from a trace.
pub fn trace_panel(app: &mut TolViewerApp, ui: &mut egui::Ui) {
    if app.current_doc().is_none() {
        return;
    }
    if app.current_doc().is_none_or(|d| d.trace.is_none()) {
        return;
    }
    let mut actions: Vec<TraceAction> = Vec::new();
    egui::Panel::bottom(egui::Id::new("trace")).resizable(true).default_size(190.0).show(
        ui,
        |ui| {
            let Some(doc) = app.current_doc_mut() else { return };
            let row = doc.trace.as_ref().map(|t| t.row).unwrap_or(0);
            let caret_col = doc.selection.cursor.col;
            let residues = match doc.alignment.sequences.get(row) {
                Some(seq) => seq.residues.clone(),
                None => return,
            };
            let caret =
                doc.alignment.sequences.get(row).and_then(|s| s.residue_index_at(caret_col));
            let Some(view) = doc.trace_view() else { return };
            let calls = view.trace.len();
            let samples = view.trace.samples();
            let spacing = view.trace.mean_peak_spacing();
            let sample_name = view.trace.sample_name.clone();
            let comment = view.trace.comment.clone();
            let link = view.link();

            ui.horizontal(|ui| {
                ui.label(RichText::new(&sample_name).strong());
                ui.label(RichText::new(format!("{calls} calls")).weak());
                if !comment.is_empty() {
                    ui.label(RichText::new(comment).weak());
                }
                match link {
                    Link::At { identity, .. } if identity < 1.0 => {
                        ui.label(
                            RichText::new(format!("{:.0}% of calls unchanged", identity * 100.0))
                                .weak(),
                        );
                    }
                    Link::Lost => {
                        ui.label(
                            RichText::new("the sequence no longer matches the trace")
                                .color(ui.visuals().warn_fg_color),
                        );
                    }
                    _ => {}
                }
                ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                    let Some(doc) = app.current_doc_mut() else { return };
                    let Some(view) = doc.trace.as_mut() else { return };
                    ui.add(
                        egui::Slider::new(&mut view.gain, 0.2..=6.0).logarithmic(true).text("gain"),
                    );
                });
            });

            // Vetting a call: with one selected, overrule the instrument
            // without having to move back to the alignment to type. The signal
            // stays as it was — only the call changes, and undo takes it back.
            if let Some(residue) = caret {
                ui.horizontal(|ui| {
                    let called = current_call(&residues, residue);
                    ui.label(RichText::new(format!("base {}", residue + 1)).weak());
                    for &base in b"ACGTN" {
                        let label = RichText::new((base as char).to_string()).monospace();
                        let label = if called == Some(base) { label.strong() } else { label };
                        if ui
                            .add(egui::Button::new(label).min_size(egui::Vec2::new(22.0, 0.0)))
                            .on_hover_text("Call this base here")
                            .clicked()
                        {
                            actions.push(TraceAction::Recall { residue, base });
                        }
                    }
                    ui.label(
                        RichText::new("or type over it in the alignment above").weak().small(),
                    );
                });
            }

            // Show a window of calls around the caret, so moving through the
            // sequence walks the trace.
            let Some(doc) = app.current_doc_mut() else { return };
            let zoom = doc.zoom;
            let Some(view) = doc.trace_view() else { return };
            let visible_calls = (ui.available_width() / zoom.max(2.0)).max(4.0);
            // Match the canvas's zoom: a call gets as much room here as a cell
            // gets there, so the trace scrolls in step with the sequence.
            let span = (visible_calls * spacing).max(1.0);
            let centre = caret
                .and_then(|i| view.sample_for_residue(i))
                .map(|s| s as f32)
                .unwrap_or(span / 2.0);
            let scroll = (centre - span / 2.0).clamp(0.0, (samples as f32 - span).max(0.0));

            let height = ui.available_height().max(60.0);
            Chromatogram {
                view,
                residues: &residues,
                caret,
                samples_per_view: span,
                scroll,
                height,
            }
            .show(ui, &mut actions);
        },
    );

    for action in actions {
        match action {
            TraceAction::Select(residue) => app.cmd_goto_residue(residue),
            TraceAction::Recall { residue, base } => app.cmd_recall(residue, base),
        }
    }
}

pub fn welcome(app: &mut TolViewerApp, ui: &mut egui::Ui) {
    ui.vertical_centered(|ui| {
        ui.add_space(60.0);
        ui.heading("TOLViewer");
        ui.label("View, edit, align and clean DNA and protein alignments.");
        ui.label(
            RichText::new(
                "The library panel on the left keeps a project's files organised without \
                 moving them.",
            )
            .weak(),
        );
        ui.add_space(18.0);
        ui.horizontal(|ui| {
            // Centred by the enclosing `vertical_centered`, so the pair of
            // buttons sits under the heading rather than against the left edge.
            ui.add_space(ui.available_width() / 2.0 - 110.0);
            if ui.button("Open a file…").clicked() {
                app.cmd_open();
            }
            if ui
                .button("Open a library…")
                .on_hover_text("A saved project: folders of reads, alignments and primers")
                .clicked()
            {
                app.cmd_library_open();
            }
        });
        ui.add_space(8.0);
        ui.label(RichText::new("or drop a FASTA, PHYLIP, NEXUS, Clustal or .ab1 file here").weak());
        let recent: Vec<_> = app.recent_files().to_vec();
        if !recent.is_empty() {
            ui.add_space(24.0);
            ui.label(RichText::new("Recent").strong());
            for path in recent.into_iter().take(6) {
                let label = path
                    .file_name()
                    .map(|s| s.to_string_lossy().into_owned())
                    .unwrap_or_else(|| path.display().to_string());
                if ui.link(label).on_hover_text(path.display().to_string()).clicked() {
                    app.cmd_open_path(path);
                }
            }
        }
    });
}

/// The question the in-situ policy raises: this save would land on a file the
/// library did not create.
///
/// It is the one dialog in the program that must not have a harmless default,
/// so neither button is pre-selected and the wording says which file is at
/// risk by name.
fn library_save_dialog(app: &mut TolViewerApp, ctx: &egui::Context) {
    if app.pending_save().is_none() {
        return;
    }
    let (name, target, mut copy_to) = {
        let pending = app.pending_save().expect("checked above");
        (pending.name.clone(), pending.target.clone(), pending.copy_to.display().to_string())
    };
    let mut answer: Option<Option<SaveChoice>> = None;

    egui::Window::new("Save changes")
        .collapsible(false)
        .resizable(false)
        .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
        .show(ctx, |ui| {
            ui.set_max_width(520.0);
            ui.label(RichText::new(&name).strong());
            ui.add_space(6.0);
            match &target {
                SaveTarget::Original(path) => {
                    ui.label("Saving would replace the original file:");
                    ui.label(RichText::new(path.display().to_string()).monospace());
                    ui.add_space(6.0);
                    ui.label(
                        "That file is the one the library points at, so replacing it \
                         replaces the data every other entry reads from.",
                    );
                }
                SaveTarget::MustCopy(_, why) => {
                    ui.label(why.explain());
                }
                SaveTarget::WorkingCopy(_) => {}
            }

            ui.add_space(10.0);
            ui.separator();
            ui.label("Save a copy instead, and keep every later edit in it:");
            ui.horizontal(|ui| {
                let response = ui.add(
                    egui::TextEdit::singleline(&mut copy_to).desired_width(360.0).hint_text("path"),
                );
                if response.changed() {
                    if let Some(pending) = app.pending_save_mut() {
                        pending.copy_to = PathBuf::from(copy_to.clone());
                    }
                }
                if ui.button("Browse…").clicked() {
                    app.cmd_browse_for_copy();
                }
            });
            ui.label(
                RichText::new(
                    "From then on this sequence saves straight to that copy, without asking \
                     again — so editing it repeatedly leaves one extra file, not a pile of them.",
                )
                .weak()
                .small(),
            );

            ui.add_space(12.0);
            ui.horizontal(|ui| {
                if ui.button("Save as a copy").clicked() {
                    let path = app
                        .pending_save()
                        .map(|p| p.copy_to.clone())
                        .unwrap_or_else(|| PathBuf::from(&copy_to));
                    answer = Some(Some(SaveChoice::NewCopy(path)));
                }
                let overwrite =
                    RichText::new("Replace the original").color(ui.visuals().error_fg_color);
                if target.can_overwrite() && ui.button(overwrite).clicked() {
                    answer = Some(Some(SaveChoice::Overwrite));
                }
                if ui.button("Cancel").clicked() {
                    answer = Some(None);
                }
            });
        });

    match answer {
        Some(Some(choice)) => app.cmd_answer_save(choice),
        Some(None) => app.cmd_cancel_save(),
        None => {}
    }
}

/// Naming a new library folder.
fn new_folder_dialog(app: &mut TolViewerApp, ctx: &egui::Context) {
    let Some(mut name) = app.dialogs_mut().new_folder().clone() else { return };
    let mut done = None;
    egui::Window::new("New folder").collapsible(false).resizable(false).show(ctx, |ui| {
        let parent = app
            .library()
            .insertion_parent()
            .map(|id| app.library().library.path_of(id))
            .unwrap_or_else(|| "the top level".to_string());
        ui.label(format!("Inside {parent}"));
        let response =
            ui.add(egui::TextEdit::singleline(&mut name).hint_text("18S").desired_width(240.0));
        response.request_focus();
        let entered = response.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
        ui.horizontal(|ui| {
            if ui.button("Create").clicked() || entered {
                done = Some(true);
            }
            if ui.button("Cancel").clicked() {
                done = Some(false);
            }
        });
    });
    match done {
        Some(true) => {
            app.cmd_library_add_folder(&name);
            *app.dialogs_mut().new_folder() = None;
        }
        Some(false) => *app.dialogs_mut().new_folder() = None,
        None => *app.dialogs_mut().new_folder() = Some(name),
    }
}

/// The library's primer list.
fn primers_dialog(app: &mut TolViewerApp, ctx: &egui::Context) {
    if !*app.dialogs_mut().primers() {
        return;
    }
    let mut open = true;
    let mut remove: Option<usize> = None;
    let mut add = false;
    egui::Window::new("Primers").open(&mut open).resizable(true).default_width(460.0).show(
        ctx,
        |ui| {
            ui.label(
                RichText::new(
                    "The primers this project amplified with. They are saved in the library \
                     and used by Map primers and Trim primers.",
                )
                .weak(),
            );
            ui.separator();
            let primers: Vec<(String, String)> = app
                .library()
                .library
                .primers
                .primers()
                .iter()
                .map(|p| (p.name.clone(), String::from_utf8_lossy(&p.sequence).into_owned()))
                .collect();
            if primers.is_empty() {
                ui.label(RichText::new("none yet").italics());
            }
            egui::Grid::new("primer-grid").num_columns(4).striped(true).show(ui, |ui| {
                for (i, (name, sequence)) in primers.iter().enumerate() {
                    ui.label(name);
                    ui.label(RichText::new(sequence).monospace());
                    ui.label(RichText::new(format!("{} nt", sequence.len())).weak());
                    if ui.small_button("Remove").clicked() {
                        remove = Some(i);
                    }
                    ui.end_row();
                }
            });
            ui.separator();
            ui.horizontal(|ui| {
                let draft = app.primer_draft_mut();
                ui.add(
                    egui::TextEdit::singleline(&mut draft.0)
                        .hint_text("LCO1490")
                        .desired_width(110.0),
                );
                ui.add(
                    egui::TextEdit::singleline(&mut draft.1)
                        .hint_text("GGTCAACAAATCATAAAGATATTGG")
                        .desired_width(240.0),
                );
                if ui.button("Add").clicked() {
                    add = true;
                }
            });
            ui.label(
                RichText::new("IUPAC codes are allowed, so a degenerate primer works as written.")
                    .weak()
                    .small(),
            );
        },
    );
    if let Some(i) = remove {
        app.cmd_remove_primer(i);
    }
    if add {
        app.cmd_add_primer();
    }
    if !open {
        *app.dialogs_mut().primers() = false;
    }
}

/// How much of a read counts as primer, and whether it comes off.
fn trim_dialog(app: &mut TolViewerApp, ctx: &egui::Context) {
    if !*app.dialogs_mut().trim() {
        return;
    }
    let mut open = true;
    let mut run = false;
    egui::Window::new("Trim primers").open(&mut open).resizable(false).show(ctx, |ui| {
        {
            let trim = &mut app.library_mut().trim;
            ui.add(
                egui::Slider::new(&mut trim.max_mismatch_fraction, 0.0..=0.4)
                    .text("mismatches allowed")
                    .custom_formatter(|v, _| format!("{:.0}%", v * 100.0)),
            )
            .on_hover_text(
                "The start of a Sanger read is where the basecaller is least sure, so a \
                 primer rarely matches perfectly.",
            );
            ui.add(
                egui::Slider::new(&mut trim.search_window, 0..=400).text("search window (bases)"),
            )
            .on_hover_text(
                "Only look for the primers this far in from each end. 0 searches the whole \
                 read, which risks matching a repeat in the middle of it.",
            );
            ui.checkbox(&mut trim.keep_primers, "Keep the primer sequence")
                .on_hover_text("Trim the junk outside the primers but leave the primers in.");
        }
        ui.separator();
        ui.label(
            RichText::new(
                "Trimming opens each sequence and edits it there. Nothing is written until \
                 you save, and you are asked before any original file is replaced.",
            )
            .weak(),
        );
        ui.add_space(8.0);
        ui.horizontal(|ui| {
            let entries = app.library().selected_entries().len();
            if ui
                .add_enabled(entries > 0, egui::Button::new(format!("Trim {entries} selected")))
                .clicked()
            {
                run = true;
            }
            if ui
                .button("Map only")
                .on_hover_text("Report where the primers bind, and change nothing")
                .clicked()
            {
                app.cmd_library_map_primers();
            }
        });
    });
    if run {
        app.cmd_library_trim_primers();
    }
    if !open {
        *app.dialogs_mut().trim() = false;
    }
}

/// What the primer run found.
fn primer_report_dialog(app: &mut TolViewerApp, ctx: &egui::Context) {
    if app.dialogs_mut().primer_report().is_none() {
        return;
    }
    let lines = app.dialogs_mut().primer_report().clone().unwrap_or_default();
    let mut open = true;
    let mut copy = false;
    egui::Window::new("Primer report").open(&mut open).resizable(true).default_width(560.0).show(
        ctx,
        |ui| {
            egui::ScrollArea::vertical().max_height(360.0).show(ui, |ui| {
                for line in &lines {
                    ui.label(RichText::new(line).monospace().small());
                }
                if lines.is_empty() {
                    ui.label(RichText::new("nothing to report").italics());
                }
            });
            ui.separator();
            if ui.button("Copy").clicked() {
                copy = true;
            }
        },
    );
    if copy {
        ctx.copy_text(lines.join("\n"));
    }
    if !open {
        *app.dialogs_mut().primer_report() = None;
    }
}

/// The concatenation: how names were matched, and what the matrix came out as.
fn concat_dialog(app: &mut TolViewerApp, ctx: &egui::Context) {
    if !*app.dialogs_mut().concat() {
        return;
    }
    let mut open = true;
    let mut keep = false;
    let mut discard = false;
    let mut copy_charsets: Option<String> = None;

    egui::Window::new("Concatenate").open(&mut open).resizable(true).default_width(620.0).show(
        ctx,
        |ui| {
            {
                let concat = &mut app.library_mut().concat;
                ui.checkbox(
                    &mut concat.matching.strip_suffixes,
                    "Match samples across loci by name",
                )
                .on_hover_text(
                    "Strips the locus and direction from a name, so TL-2213_18S_F and \
                     TL_2213_28S are recognised as one specimen.",
                );
                ui.checkbox(
                    &mut concat.include_partial,
                    "Keep samples that are missing from some loci",
                )
                .on_hover_text("Missing loci are filled with gaps, which is what phylogenetic programs read as missing data.");
            }
            ui.label(
                RichText::new("Re-run Concatenate after changing these.").weak().small(),
            );
            ui.separator();

            let Some(result) = app.concat_result_view() else {
                ui.label("Select two or more alignments in the library and run Concatenate.");
                return;
            };
            ui.label(format!(
                "{} samples across {} loci, {} columns. {} are complete.",
                result.alignment.len(),
                result.partitions.len(),
                result.alignment.width(),
                result.complete
            ));

            ui.add_space(6.0);
            egui::ScrollArea::vertical().max_height(280.0).show(ui, |ui| {
                egui::Grid::new("partitions").num_columns(2).striped(true).show(ui, |ui| {
                    for p in &result.partitions {
                        ui.label(&p.name);
                        ui.label(RichText::new(p.as_charset()).monospace());
                        ui.end_row();
                    }
                });
                if !result.missing.is_empty() {
                    ui.add_space(8.0);
                    ui.label(
                        RichText::new(format!(
                            "{} sample(s) were padded with gaps for at least one locus. \
                             A sample missing everywhere but one locus is usually a name \
                             that failed to match rather than a gap in the sampling:",
                            result.missing.len()
                        ))
                        .color(ui.visuals().warn_fg_color),
                    );
                    for m in result.missing.iter().take(40) {
                        ui.label(
                            RichText::new(format!(
                                "{} — missing from {}",
                                m.sample,
                                m.absent_from.join(", ")
                            ))
                            .small(),
                        );
                    }
                }
                if !result.dropped.is_empty() {
                    ui.add_space(8.0);
                    ui.label(format!("{} incomplete sample(s) were dropped.", result.dropped.len()));
                }
            });

            ui.separator();
            ui.horizontal(|ui| {
                if ui.button("Open as a new tab").clicked() {
                    keep = true;
                }
                if ui
                    .button("Copy NEXUS charsets")
                    .on_hover_text("For a partitioned analysis downstream")
                    .clicked()
                {
                    copy_charsets = Some(result.nexus_charsets());
                }
                if ui.button("Discard").clicked() {
                    discard = true;
                }
            });
        },
    );

    if let Some(text) = copy_charsets {
        ctx.copy_text(text);
    }
    if keep {
        app.cmd_keep_concatenation();
    }
    if discard || !open {
        app.clear_concatenation();
        *app.dialogs_mut().concat() = false;
    }
}

pub fn dialogs(app: &mut TolViewerApp, ctx: &egui::Context) {
    align_dialog(app, ctx);
    clean_dialog(app, ctx);
    export_dialog(app, ctx);
    goto_dialog(app, ctx);
    rename_dialog(app, ctx);
    close_dialog(app, ctx);
    about_dialog(app, ctx);
    new_folder_dialog(app, ctx);
    primers_dialog(app, ctx);
    trim_dialog(app, ctx);
    primer_report_dialog(app, ctx);
    concat_dialog(app, ctx);
    // Last, so it draws over anything else: it is the only dialog that stands
    // between the user and a file being replaced.
    library_save_dialog(app, ctx);
}

fn align_dialog(app: &mut TolViewerApp, ctx: &egui::Context) {
    let mut open = *app.dialogs_mut().align();
    if !open {
        return;
    }
    let rows = app.current_doc().map_or(0, |d| d.rows());
    let mut start = false;
    egui::Window::new("Alignment settings")
        .open(&mut open)
        .resizable(false)
        .show(ctx, |ui| {
            let mut engine_changed = None;
            {
                let params = app.align_params_mut();
                ui.horizontal(|ui| {
                    ui.label("Engine");
                    for engine in Engine::all() {
                        if ui.radio_value(&mut params.engine, *engine, engine.name()).changed() {
                            engine_changed = Some(*engine);
                        }
                    }
                });
                ui.label(
                    RichText::new(match params.engine {
                        Engine::Clustal => "Progressive alignment on a neighbour-joining guide tree. The fastest of the three, and the most predictable.",
                        Engine::Muscle => "Progressive draft, then iterative refinement. Usually the most accurate on divergent sets, but much the slowest: refinement re-aligns across every tree edge each round.",
                        Engine::Mafft => "FFT-accelerated group-to-group alignment. Its advantage grows with the number of sequences; on small sets it is no faster than Clustal.",
                    })
                    .weak()
                    .small(),
                );
                if let Some(warning) = cost_warning(params.engine, params.iterations, rows) {
                    ui.colored_label(ui.visuals().warn_fg_color, warning);
                }
                ui.separator();
                egui::Grid::new("align-grid").num_columns(2).show(ui, |ui| {
                    ui.label("Substitution matrix");
                    egui::ComboBox::from_id_salt("matrix")
                        .selected_text(
                            MATRIX_CHOICES
                                .iter()
                                .find(|(m, _)| *m == params.matrix)
                                .map(|(_, n)| *n)
                                .unwrap_or("Automatic"),
                        )
                        .show_ui(ui, |ui| {
                            for (choice, name) in MATRIX_CHOICES {
                                ui.selectable_value(&mut params.matrix, *choice, *name);
                            }
                        });
                    ui.end_row();

                    ui.label("Gap open");
                    ui.add(egui::DragValue::new(&mut params.gap_open).speed(0.1).range(0.0..=100.0));
                    ui.end_row();

                    ui.label("Gap extend");
                    ui.add(egui::DragValue::new(&mut params.gap_extend).speed(0.05).range(0.0..=50.0));
                    ui.end_row();

                    ui.label("Terminal gap factor");
                    ui.add(
                        egui::Slider::new(&mut params.terminal_gap_factor, 0.0..=1.0)
                            .text("0 = free ends"),
                    );
                    ui.end_row();

                    ui.label("Guide tree");
                    egui::ComboBox::from_id_salt("tree")
                        .selected_text(
                            TREE_CHOICES
                                .iter()
                                .find(|(t, _)| *t == params.tree)
                                .map(|(_, n)| *n)
                                .unwrap_or("Neighbour joining"),
                        )
                        .show_ui(ui, |ui| {
                            for (choice, name) in TREE_CHOICES {
                                ui.selectable_value(&mut params.tree, *choice, *name);
                            }
                        });
                    ui.end_row();

                    ui.label("Refinement rounds");
                    ui.add(egui::DragValue::new(&mut params.iterations).range(0..=32));
                    ui.end_row();

                    ui.label("Threads");
                    ui.add(
                        egui::DragValue::new(&mut params.threads)
                            .range(0..=256)
                            .custom_formatter(|n, _| {
                                if n == 0.0 { "all cores".into() } else { format!("{n}") }
                            }),
                    );
                    ui.end_row();
                });
            }
            if let Some(engine) = engine_changed {
                // Switching engine should bring its own tuned defaults, not
                // carry over penalties calibrated for a different scoring model.
                *app.align_params_mut() = crate::app::engine_defaults(engine);
            }
            ui.separator();
            ui.horizontal(|ui| {
                if ui.button("Align now").clicked() {
                    start = true;
                }
                if ui.button("Reset to defaults").clicked() {
                    let engine = app.align_params_mut().engine;
                    *app.align_params_mut() = crate::app::engine_defaults(engine);
                }
            });
        });
    if start {
        app.cmd_align(ctx, false);
        open = false;
    }
    *app.dialogs_mut().align() = open;
}

fn clean_dialog(app: &mut TolViewerApp, ctx: &egui::Context) {
    let mut open = *app.dialogs_mut().clean();
    if !open {
        return;
    }
    let rows = app.current_doc().map_or(0, |d| d.rows());
    if app.clean_params_mut().is_none() {
        *app.clean_params_mut() = Some(GblocksParams::defaults(rows));
    }
    let mut run = false;
    let mut apply = false;
    let mut validation: Option<String> = None;

    egui::Window::new("Gblocks — select conserved blocks").open(&mut open).resizable(false).show(
        ctx,
        |ui| {
            {
                let params = app.clean_params_mut().as_mut().expect("set above");
                egui::Grid::new("gblocks-grid").num_columns(2).show(ui, |ui| {
                    ui.label("b1 minimum for a conserved position");
                    ui.add(
                        egui::DragValue::new(&mut params.min_seqs_conserved).range(1..=rows.max(1)),
                    );
                    ui.end_row();
                    ui.label("b2 minimum for a flank position");
                    ui.add(egui::DragValue::new(&mut params.min_seqs_flank).range(1..=rows.max(1)));
                    ui.end_row();
                    ui.label("b3 maximum contiguous non-conserved");
                    ui.add(
                        egui::DragValue::new(&mut params.max_contiguous_nonconserved)
                            .range(1..=1000),
                    );
                    ui.end_row();
                    ui.label("b4 minimum block length");
                    ui.add(egui::DragValue::new(&mut params.min_block_length).range(2..=1000));
                    ui.end_row();
                    ui.label("b5 allowed gap positions");
                    egui::ComboBox::from_id_salt("gaps")
                        .selected_text(
                            GAP_POLICIES
                                .iter()
                                .find(|(p, _)| *p == params.gaps)
                                .map(|(_, n)| *n)
                                .unwrap_or("No gaps allowed"),
                        )
                        .show_ui(ui, |ui| {
                            for (policy, name) in GAP_POLICIES {
                                ui.selectable_value(&mut params.gaps, *policy, *name);
                            }
                        });
                    ui.end_row();
                    ui.label("Count similar residues");
                    ui.checkbox(&mut params.use_similarity, "").on_hover_text(
                        "For protein: treat positively scoring residues as conserved",
                    );
                    ui.end_row();
                });
                if let Err(e) = params.validate(rows) {
                    validation = Some(e.to_string());
                }
            }
            ui.horizontal(|ui| {
                if ui.button("Gblocks defaults").clicked() {
                    *app.clean_params_mut() = Some(GblocksParams::defaults(rows));
                }
                if ui
                    .button("Relaxed (Talavera & Castresana)")
                    .on_hover_text("Less stringent settings recommended for phylogenetics")
                    .clicked()
                {
                    *app.clean_params_mut() = Some(GblocksParams::relaxed(rows));
                }
            });
            if let Some(message) = &validation {
                ui.colored_label(ui.visuals().error_fg_color, message);
            }
            ui.separator();
            if let Some(result) = app.pending_clean() {
                ui.label(format!(
                    "Preview: keeping {} of {} columns ({:.0}%) in {} block(s).",
                    result.kept,
                    result.total,
                    result.kept_fraction() * 100.0,
                    result.blocks.len()
                ));
                ui.label(
                    RichText::new("The green track above the alignment shows the kept columns.")
                        .weak()
                        .small(),
                );
            } else {
                ui.label(RichText::new("Run a preview to see which columns survive.").weak());
            }
            ui.separator();
            ui.horizontal(|ui| {
                if ui.add_enabled(validation.is_none(), egui::Button::new("Preview")).clicked() {
                    run = true;
                }
                if ui
                    .add_enabled(
                        app.pending_clean().is_some(),
                        egui::Button::new("Apply to alignment"),
                    )
                    .clicked()
                {
                    apply = true;
                }
            });
        },
    );

    if run {
        app.cmd_clean(ctx);
    }
    if apply {
        if let Some(result) = app.take_pending_clean() {
            app.cmd_apply_mask(result.mask, "Gblocks");
        }
        open = false;
    }
    *app.dialogs_mut().clean() = open;
}

fn export_dialog(app: &mut TolViewerApp, ctx: &egui::Context) {
    let mut open = *app.dialogs_mut().export();
    if !open {
        return;
    }
    let mut do_export: Option<bool> = None;
    egui::Window::new("Export alignment").open(&mut open).resizable(false).show(ctx, |ui| {
        let format = *app.export_format_mut();
        egui::ComboBox::from_label("Format").selected_text(format.name()).show_ui(ui, |ui| {
            for f in Format::all().iter().filter(|f| f.can_write()) {
                ui.selectable_value(app.export_format_mut(), *f, f.name());
            }
        });
        ui.separator();
        ui.label(
            RichText::new(match format {
                Format::Phylip => "Strict PHYLIP truncates names to 10 characters.",
                Format::PhylipRelaxed => {
                    "Relaxed PHYLIP keeps full names; RAxML and IQ-TREE read this."
                }
                Format::Nexus => "NEXUS carries the sequence type and quotes awkward names.",
                _ => "",
            })
            .weak()
            .small(),
        );
        ui.separator();
        ui.horizontal(|ui| {
            if ui.button("Export everything…").clicked() {
                do_export = Some(false);
            }
            if ui.button("Export selection…").clicked() {
                do_export = Some(true);
            }
        });
    });

    if let Some(selection_only) = do_export {
        let format = *app.export_format_mut();
        let suggested = app
            .current_doc()
            .and_then(|d| d.path.clone())
            .and_then(|p| p.file_stem().map(|s| s.to_string_lossy().into_owned()))
            .unwrap_or_else(|| "alignment".to_string());
        let picked = rfd::FileDialog::new()
            .set_title("Export alignment")
            .set_file_name(format!("{suggested}.{}", format.extensions()[0]))
            .add_filter(format.name(), format.extensions())
            .save_file();
        if let Some(path) = picked {
            app.export_to(path, format, selection_only);
            open = false;
        }
    }
    *app.dialogs_mut().export() = open;
}

fn goto_dialog(app: &mut TolViewerApp, ctx: &egui::Context) {
    let mut open = *app.dialogs_mut().goto();
    if !open {
        return;
    }
    let mut go: Option<usize> = None;
    let mut bad = false;
    egui::Window::new("Go to column").open(&mut open).resizable(false).show(ctx, |ui| {
        let response = ui.add(
            egui::TextEdit::singleline(app.goto_input_mut())
                .hint_text("column number")
                .desired_width(120.0),
        );
        response.request_focus();
        let submit = response.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
        if ui.button("Go").clicked() || submit {
            match app.goto_input_mut().trim().parse::<usize>() {
                Ok(n) if n >= 1 => go = Some(n),
                _ => bad = true,
            }
        }
        if bad {
            ui.colored_label(ui.visuals().error_fg_color, "enter a column number of 1 or more");
        }
    });
    if let Some(column) = go {
        app.cmd_goto(column);
        app.goto_input_mut().clear();
        open = false;
    }
    *app.dialogs_mut().goto() = open;
}

fn rename_dialog(app: &mut TolViewerApp, ctx: &egui::Context) {
    let Some((row, mut name)) = app.dialogs_mut().rename().clone() else { return };
    let mut open = true;
    let mut commit = false;
    egui::Window::new("Rename sequence").open(&mut open).resizable(false).show(ctx, |ui| {
        let response = ui.add(
            egui::TextEdit::singleline(&mut name)
                .desired_width(320.0)
                .hint_text("name description"),
        );
        response.request_focus();
        let submit = response.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
        if ui.button("Rename").clicked() || submit {
            commit = true;
        }
    });
    if commit {
        app.cmd_rename(row, &name);
        *app.dialogs_mut().rename() = None;
    } else if !open {
        *app.dialogs_mut().rename() = None;
    } else {
        *app.dialogs_mut().rename() = Some((row, name));
    }
}

fn close_dialog(app: &mut TolViewerApp, ctx: &egui::Context) {
    let Some(index) = *app.dialogs_mut().confirm_close() else { return };
    let title = app.documents().get(index).map(|d| d.title()).unwrap_or_default();
    let mut decided: Option<bool> = None;
    egui::Window::new("Unsaved changes")
        .collapsible(false)
        .resizable(false)
        .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
        .show(ctx, |ui| {
            ui.label(format!("{title} has unsaved changes."));
            ui.add_space(8.0);
            ui.horizontal(|ui| {
                if ui.button("Save").clicked() {
                    app.set_current(index);
                    app.cmd_save(false);
                    decided = Some(true);
                }
                if ui.button("Discard").clicked() {
                    app.cmd_force_close(index);
                    decided = Some(true);
                }
                if ui.button("Cancel").clicked() {
                    decided = Some(false);
                }
            });
        });
    match decided {
        Some(true) => {
            *app.dialogs_mut().confirm_close() = None;
            // Saving or discarding may have been the last thing standing
            // between the user and the quit they asked for.
            app.cmd_resume_quit(ctx);
        }
        Some(false) => {
            *app.dialogs_mut().confirm_close() = None;
            app.cmd_cancel_quit();
        }
        None => {}
    }
}

fn about_dialog(app: &mut TolViewerApp, ctx: &egui::Context) {
    let mut open = *app.dialogs_mut().about();
    if !open {
        return;
    }
    egui::Window::new("About TOLViewer").open(&mut open).resizable(false).show(ctx, |ui| {
        ui.heading("TOLViewer");
        ui.label(format!("version {}", env!("CARGO_PKG_VERSION")));
        ui.add_space(8.0);
        ui.label(
            "A sequence and alignment viewer/editor with a project library, Sanger trace \
             viewing, and built-in alignment and cleaning.",
        );
        ui.add_space(8.0);
        ui.label(
            RichText::new(
                "Alignment engines reimplement the published algorithms of ClustalW \
                 (Thompson et al. 1994), MUSCLE (Edgar 2004) and MAFFT (Katoh et al. 2002). \
                 Cleaning reimplements Gblocks (Castresana 2000). Results will not be \
                 byte-identical to those programs.",
            )
            .small()
            .weak(),
        );
        ui.add_space(8.0);
        ui.hyperlink_to("github.com/tingidlab/TOLViewer", "https://github.com/tingidlab/TOLViewer");
    });
    *app.dialogs_mut().about() = open;
}

/// Warn before the user commits to a run that will take minutes.
///
/// The thresholds come from the measured benchmarks in
/// `crates/tolviewer-align/tests/accuracy.rs`: at 200 sequences of ~1000
/// columns Clustal takes about 2 s and MAFFT about 5 s, while MUSCLE with two
/// refinement rounds takes about 34 s, and refinement scales far worse than
/// linearly in the number of sequences.
fn cost_warning(engine: Engine, iterations: usize, rows: usize) -> Option<String> {
    if rows < 100 {
        return None;
    }
    match engine {
        Engine::Muscle if iterations > 0 => Some(format!(
            "{rows} sequences with {iterations} refinement round(s) may take several minutes. \
             Set refinement rounds to 0, or use MAFFT, for a faster first look.",
        )),
        _ if rows >= 1000 => {
            Some(format!("{rows} sequences will take a while; the job runs in the background and can be cancelled."))
        }
        _ => None,
    }
}
