//! Menus, panels and dialogs.
//!
//! Everything here reads and writes `TolViewerApp` through its `cmd_*`
//! accessors so the widget code stays declarative and the state transitions
//! all live in `app.rs`.

use egui::{Align, Layout, RichText};
use tolviewer_align::Engine;
use tolviewer_clean::GblocksParams;
use tolviewer_core::Alphabet;
use tolviewer_io::Format;

use crate::app::{TolViewerApp, GAP_POLICIES, MATRIX_CHOICES, SCHEMES, TREE_CHOICES};
use crate::canvas::{MAX_ZOOM, MIN_ZOOM};

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

pub fn welcome(app: &mut TolViewerApp, ui: &mut egui::Ui) {
    ui.vertical_centered(|ui| {
        ui.add_space(60.0);
        ui.heading("TOLViewer");
        ui.label("View, edit, align and clean DNA and protein alignments.");
        ui.add_space(18.0);
        if ui.button("Open a file…").clicked() {
            app.cmd_open();
        }
        ui.add_space(8.0);
        ui.label(RichText::new("or drop a FASTA, PHYLIP, NEXUS or Clustal file here").weak());
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

pub fn dialogs(app: &mut TolViewerApp, ctx: &egui::Context) {
    align_dialog(app, ctx);
    clean_dialog(app, ctx);
    export_dialog(app, ctx);
    goto_dialog(app, ctx);
    rename_dialog(app, ctx);
    close_dialog(app, ctx);
    about_dialog(app, ctx);
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
        ui.label("A sequence and alignment viewer/editor with built-in alignment and cleaning.");
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
