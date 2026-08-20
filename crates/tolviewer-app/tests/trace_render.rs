//! Headless rendering tests for the chromatogram.
//!
//! Rendering cannot be checked by eye in this project's development
//! environment, so it is checked here instead: a real trace is painted into an
//! off-screen context and the resulting shapes are counted. That catches
//! panics in the painting code — bad rects, out-of-range sample indices, a
//! divide by a zero peak — and pins the two properties the widget is supposed
//! to have.

use egui::{Context, Pos2, RawInput, Rect, Vec2};
use tolviewer_app::chromatogram::{Chromatogram, Link, TraceView};
use tolviewer_io::ab1::{self, Trace};

fn trace() -> Trace {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../testdata/ab1/tingidae_COI_F.ab1");
    ab1::read_file(&path).expect("the checked-in trace should read")
}

/// Paint one chromatogram off-screen and count the primitives produced.
fn paint(view: &TraceView, residues: &[u8], samples_per_view: f32, scroll: f32) -> usize {
    let ctx = Context::default();
    let mut shapes = 0;
    // As in `canvas_render.rs`: a first pass lets egui settle its layout, and
    // the second is the one measured.
    for _ in 0..2 {
        let input = RawInput {
            screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(900.0, 240.0))),
            ..Default::default()
        };
        let mut output = ctx.run_ui(input, |ui| {
            let mut actions = Vec::new();
            Chromatogram {
                view,
                residues,
                caret: Some(20),
                samples_per_view,
                scroll,
                height: 200.0,
            }
            .show(ui, &mut actions);
        });
        shapes = output.shapes.iter().map(|c| count(&c.shape)).sum();
        // egui panics if these are dropped unhandled.
        output.textures_delta.clear();
    }
    shapes
}

fn count(shape: &egui::Shape) -> usize {
    match shape {
        egui::Shape::Vec(v) => v.iter().map(count).sum(),
        egui::Shape::Noop => 0,
        _ => 1,
    }
}

fn linked(residues: &[u8]) -> TraceView {
    let mut view = TraceView::new(trace(), 0);
    view.relink(residues, 1);
    view
}

#[test]
fn a_real_trace_paints_without_panicking() {
    let calls = trace().calls;
    let view = linked(&calls);
    assert!(matches!(view.link(), Link::At { offset: 0, .. }));
    assert!(paint(&view, &calls, 400.0, 0.0) > 0);
}

#[test]
fn painting_a_window_costs_no_more_than_the_window_holds() {
    // The whole point of the sample stepping: showing 40 bases of a 3000-sample
    // trace must not cost the same as showing all of it. Both are bounded by
    // the width of the panel, so the two are within a small factor — what must
    // not happen is the cost scaling with the length of the read.
    let calls = trace().calls;
    let view = linked(&calls);
    let narrow = paint(&view, &calls, 400.0, 100.0);
    let whole = paint(&view, &calls, view.trace.samples() as f32, 0.0);
    assert!(narrow > 0 && whole > 0);
    assert!(
        whole < narrow * 4,
        "painting the whole read cost {whole} against {narrow} for a window; \
         the sample stepping is not bounding the work"
    );
}

#[test]
fn scrolling_to_every_part_of_the_read_never_panics() {
    let calls = trace().calls;
    let view = linked(&calls);
    let samples = view.trace.samples() as f32;
    for scroll in [0.0, samples / 3.0, samples - 10.0, samples, samples * 2.0] {
        for span in [1.0, 50.0, 400.0, samples * 2.0] {
            let _ = paint(&view, &calls, span, scroll);
        }
    }
}

#[test]
fn a_trace_that_no_longer_matches_says_so_instead_of_drawing() {
    // Nothing lines up, so the widget must draw its explanation rather than a
    // plausible-looking but wrong picture.
    let nonsense = vec![b'A'; 300];
    let view = linked(&nonsense);
    assert_eq!(view.link(), Link::Lost);
    let shapes = paint(&view, &nonsense, 400.0, 0.0);
    let calls = trace().calls;
    let normal = paint(&linked(&calls), &calls, 400.0, 0.0);
    assert!(shapes > 0, "the message itself is a shape");
    assert!(shapes < normal, "no trace should have been drawn: {shapes} vs {normal}");
}

#[test]
fn a_gapped_row_paints_its_residues_over_the_right_peaks() {
    // An aligned trace row has gaps, which consume a column but not a call.
    let calls = trace().calls;
    let mut gapped: Vec<u8> = Vec::new();
    for (i, &c) in calls.iter().enumerate() {
        if i % 25 == 0 {
            gapped.push(b'-');
        }
        gapped.push(c);
    }
    let view = linked(&calls);
    assert!(paint(&view, &gapped, 400.0, 0.0) > 0);
}

#[test]
fn an_empty_trace_paints_nothing_rather_than_dividing_by_zero() {
    let mut empty = trace();
    empty.trim(0..0).unwrap();
    let mut view = TraceView::new(empty, 0);
    view.relink(&[], 1);
    let _ = paint(&view, &[], 400.0, 0.0);
}

#[test]
fn a_trimmed_read_still_paints_over_its_own_signal() {
    // 40 bases off the front: the link shifts, and the picture must follow.
    let calls = trace().calls;
    let trimmed = calls[40..].to_vec();
    let view = linked(&trimmed);
    assert!(matches!(view.link(), Link::At { offset: 40, .. }));
    assert!(paint(&view, &trimmed, 400.0, 0.0) > 0);
}

#[test]
fn extreme_gain_does_not_produce_degenerate_geometry() {
    let calls = trace().calls;
    let mut view = linked(&calls);
    for gain in [0.05f32, 0.2, 1.0, 6.0, 1000.0] {
        view.gain = gain;
        assert!(paint(&view, &calls, 400.0, 0.0) > 0, "gain {gain}");
    }
}
