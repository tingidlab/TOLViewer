//! Headless rendering tests for the alignment canvas.
//!
//! egui can run a full pass without a window, so these paint real alignments
//! and inspect the resulting shapes. That catches two classes of bug the unit
//! tests cannot: panics in the painting code (bad rects, out-of-range indices)
//! and regressions in the virtualisation, which is the whole reason the canvas
//! stays fast on large files.

use egui::{Context, Pos2, RawInput, Rect, Vec2};
use tolviewer_app::canvas::{AlignmentCanvas, ViewSettings};
use tolviewer_app::Document;
use tolviewer_core::{Alignment, Sequence};
use tolviewer_io::Format;

/// Paint `doc` once into an off-screen context of the given size and count the
/// primitives produced.
fn paint(doc: &mut Document, view: &ViewSettings, size: Vec2) -> usize {
    let ctx = Context::default();
    // A first pass lets egui settle its layout and scroll state; the second is
    // the one measured, so counts are not distorted by first-frame guesses.
    let mut shapes = 0;
    for _ in 0..2 {
        let input = RawInput {
            screen_rect: Some(Rect::from_min_size(Pos2::ZERO, size)),
            ..Default::default()
        };
        let mut output = ctx.run_ui(input, |ui| {
            let mut actions = Vec::new();
            AlignmentCanvas { doc, view, id_salt: egui::Id::new("test-canvas") }
                .show(ui, &mut actions);
        });
        shapes = output.shapes.iter().map(|clipped| count_shape(&clipped.shape)).sum();
        // A real backend uploads these to the GPU; egui panics if they are
        // dropped unhandled, so discard them explicitly.
        output.textures_delta.clear();
    }
    shapes
}

fn count_shape(shape: &egui::Shape) -> usize {
    match shape {
        egui::Shape::Vec(v) => v.iter().map(count_shape).sum(),
        egui::Shape::Noop => 0,
        _ => 1,
    }
}

fn synthetic(rows: usize, cols: usize) -> Document {
    let bases = b"ACGT";
    let sequences = (0..rows)
        .map(|r| {
            let residues: Vec<u8> = (0..cols)
                .map(|c| {
                    // A deterministic pattern with some variation per row, so
                    // consensus and conservation have something to chew on.
                    if (r + c) % 37 == 0 {
                        b'-'
                    } else {
                        bases[(c + r / 5) % 4]
                    }
                })
                .collect();
            Sequence::new(format!("seq{r:04}"), residues)
        })
        .collect();
    Document::new(Alignment::new("synthetic", sequences), None, Format::Fasta)
}

#[test]
fn paints_a_small_alignment_without_panicking() {
    let mut doc = synthetic(20, 200);
    let shapes = paint(&mut doc, &ViewSettings::default(), Vec2::new(1200.0, 800.0));
    assert!(shapes > 0, "the canvas produced nothing to draw");
}

#[test]
fn a_huge_alignment_costs_no_more_to_paint_than_a_screenful() {
    let view = ViewSettings::default();
    let size = Vec2::new(1200.0, 800.0);

    let mut small = synthetic(40, 400);
    let small_shapes = paint(&mut small, &view, size);

    // Two orders of magnitude more data, same window.
    let mut huge = synthetic(4000, 40_000);
    let huge_shapes = paint(&mut huge, &view, size);

    assert!(
        huge_shapes < small_shapes * 4,
        "virtualisation regressed: {small_shapes} shapes for 40x400 but {huge_shapes} for 4000x40000"
    );
}

#[test]
fn painting_a_huge_alignment_is_fast() {
    let mut doc = synthetic(4000, 40_000);
    let view = ViewSettings::default();
    let size = Vec2::new(1600.0, 1000.0);
    // Warm the caches first so this measures steady-state painting.
    paint(&mut doc, &view, size);

    let start = std::time::Instant::now();
    paint(&mut doc, &view, size);
    let elapsed = start.elapsed();
    // Two passes happen inside `paint`, and debug builds are slow; this is a
    // regression guard, not a frame-rate target.
    assert!(
        elapsed < std::time::Duration::from_secs(2),
        "painting a 4000 x 40000 alignment took {elapsed:?}"
    );
}

#[test]
fn zooming_all_the_way_out_still_paints() {
    let mut doc = synthetic(500, 5000);
    doc.zoom = 1.0;
    let shapes = paint(&mut doc, &ViewSettings::default(), Vec2::new(1200.0, 800.0));
    assert!(shapes > 0, "nothing drawn at minimum zoom");
}

#[test]
fn zooming_all_the_way_in_still_paints() {
    let mut doc = synthetic(50, 500);
    doc.zoom = 40.0;
    let shapes = paint(&mut doc, &ViewSettings::default(), Vec2::new(1200.0, 800.0));
    assert!(shapes > 0, "nothing drawn at maximum zoom");
}

#[test]
fn every_colour_scheme_paints() {
    for &scheme in tolviewer_app::theme::ColorScheme::ALL {
        let mut doc = synthetic(30, 300);
        let view = ViewSettings { scheme, ..ViewSettings::default() };
        let shapes = paint(&mut doc, &view, Vec2::new(1000.0, 700.0));
        assert!(shapes > 0, "{scheme:?} drew nothing");
    }
}

#[test]
fn an_empty_document_paints_without_panicking() {
    let mut doc = Document::new(Alignment::new("empty", Vec::new()), None, Format::Fasta);
    paint(&mut doc, &ViewSettings::default(), Vec2::new(800.0, 600.0));
}

#[test]
fn a_ragged_alignment_paints_without_panicking() {
    let doc_rows = vec![
        Sequence::new("long", vec![b'A'; 500]),
        Sequence::new("short", vec![b'C'; 3]),
        Sequence::new("empty", Vec::new()),
    ];
    let mut doc = Document::new(Alignment::new("ragged", doc_rows), None, Format::Fasta);
    let shapes = paint(&mut doc, &ViewSettings::default(), Vec2::new(900.0, 600.0));
    assert!(shapes > 0);
}

#[test]
fn a_tiny_window_paints_without_panicking() {
    let mut doc = synthetic(100, 1000);
    // Smaller than the name gutter, which must not produce a negative rect.
    paint(&mut doc, &ViewSettings::default(), Vec2::new(80.0, 40.0));
}

#[test]
fn the_clean_mask_track_paints() {
    let mut doc = synthetic(20, 200);
    let mask: Vec<bool> = (0..doc.width()).map(|c| c % 3 != 0).collect();
    doc.set_clean_mask(mask);
    assert!(doc.live_clean_mask().is_some());
    let shapes = paint(&mut doc, &ViewSettings::default(), Vec2::new(1000.0, 700.0));
    assert!(shapes > 0);
}

#[test]
fn a_selection_paints_its_highlight() {
    let mut doc = synthetic(20, 200);
    let bare = paint(&mut doc, &ViewSettings::default(), Vec2::new(1000.0, 700.0));
    doc.select_all();
    let selected = paint(&mut doc, &ViewSettings::default(), Vec2::new(1000.0, 700.0));
    assert!(selected > bare, "selecting everything added no highlight ({bare} -> {selected})");
}
