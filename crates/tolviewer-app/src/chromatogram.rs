//! The chromatogram: the four-channel signal a Sanger basecall was made from.
//!
//! A base call is a judgement about a wobbly analogue trace, and the calls at
//! the ends of a read — and at every heterozygous site — are the ones a person
//! has to look at. This draws the trace under the sequence so they can, and
//! lets them retype a call they disagree with.
//!
//! ## Keeping the picture honest after an edit
//!
//! The sequence is edited through the ordinary undo stack, which knows nothing
//! about traces. Trimming twenty bases off the 5' end therefore leaves the
//! row's first residue sitting over the trace's twenty-first peak, and undoing
//! it puts things back. Rather than trying to keep a second history in step,
//! [`TraceView::relink`] re-derives the correspondence from the data itself:
//! it finds the offset at which the row's residues line up against the file's
//! calls. That is right after an edit, after an undo, and after a redo, because
//! it never remembers anything it could get wrong.
//!
//! When no offset fits — the row has been rewritten past recognition — the
//! panel says so instead of drawing a plausible lie.

use egui::{Align2, FontId, Pos2, Rect, Sense, Stroke, StrokeKind, Vec2};
use tolviewer_io::ab1::Trace;

use crate::theme::base_color;

/// How the row's residues correspond to the trace's calls.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Link {
    /// Residue `i` of the row is call `i + offset` of the trace, and
    /// `identity` of them agree.
    At { offset: usize, identity: f32 },
    /// Nothing lines up: the sequence has been changed past the point where
    /// the trace can be matched to it.
    Lost,
}

/// Below this fraction of matching bases the offset is a coincidence, not a
/// correspondence. Sanger reads have real errors, and a person may have
/// retyped several calls, so the bar is not high — but a wrong offset lines up
/// at about a quarter for DNA, so anything over half is unambiguous.
const MIN_IDENTITY: f32 = 0.55;

/// A trace and the document row it belongs to.
pub struct TraceView {
    pub trace: Trace,
    /// Which row of the document the calls were read into.
    pub row: usize,
    link: Link,
    /// The document revision `link` was computed at.
    linked_at: Option<u64>,
    /// Vertical scale, as a multiple of the trace's own tallest peak.
    pub gain: f32,
}

impl TraceView {
    pub fn new(trace: Trace, row: usize) -> Self {
        TraceView { trace, row, link: Link::Lost, linked_at: None, gain: 1.0 }
    }

    pub fn link(&self) -> Link {
        self.link
    }

    /// Recompute the correspondence, unless it was already done at this
    /// revision.
    pub fn relink(&mut self, residues: &[u8], revision: u64) {
        if self.linked_at == Some(revision) {
            return;
        }
        self.link = find_offset(residues, &self.trace.calls);
        self.linked_at = Some(revision);
    }

    /// Does the trace carry any signal? A trimmed-to-nothing trace has calls
    /// but no samples to draw them over.
    pub fn samples_are_present(&self) -> bool {
        self.trace.samples() > 0
    }

    /// The trace sample residue `index` of the row peaks at.
    pub fn sample_for_residue(&self, index: usize) -> Option<u32> {
        match self.link {
            Link::At { offset, .. } => self.trace.peaks.get(index + offset).copied(),
            Link::Lost => None,
        }
    }

    /// The instrument's confidence in residue `index`, if it still corresponds
    /// to a call.
    pub fn quality_for_residue(&self, index: usize) -> Option<u8> {
        let Link::At { offset, .. } = self.link else { return None };
        self.trace.quality.as_ref()?.get(index + offset).copied()
    }

    /// The call the instrument made at residue `index`, which may differ from
    /// what the row says once someone has edited it.
    pub fn call_for_residue(&self, index: usize) -> Option<u8> {
        let Link::At { offset, .. } = self.link else { return None };
        self.trace.calls.get(index + offset).copied()
    }
}

/// Find the offset at which `residues` sit inside `calls`.
///
/// An exhaustive scan: reads are ~1000 bases and this runs once per edit, so
/// the million comparisons in the worst case are far cheaper than the
/// bookkeeping an incremental scheme would need. Ties go to the smallest
/// offset, which is the one an untrimmed read has.
fn find_offset(residues: &[u8], calls: &[u8]) -> Link {
    if residues.is_empty() || calls.is_empty() || residues.len() > calls.len() {
        return Link::Lost;
    }
    let mut best = (0usize, 0usize); // (offset, matches)
    for offset in 0..=calls.len() - residues.len() {
        let matches = residues
            .iter()
            .zip(&calls[offset..])
            .filter(|(a, b)| a.eq_ignore_ascii_case(b))
            .count();
        if matches > best.1 {
            best = (offset, matches);
            if matches == residues.len() {
                break; // cannot do better
            }
        }
    }
    let identity = best.1 as f32 / residues.len() as f32;
    if identity >= MIN_IDENTITY {
        Link::At { offset: best.0, identity }
    } else {
        Link::Lost
    }
}

/// What the user did in the chromatogram.
pub enum TraceAction {
    /// Move the caret to this residue of the row.
    Select(usize),
    /// Replace the residue at this index with this base.
    Recall { residue: usize, base: u8 },
}

/// The chromatogram widget.
pub struct Chromatogram<'a> {
    pub view: &'a TraceView,
    /// The row as it currently reads, gaps and all.
    pub residues: &'a [u8],
    /// The residue index the caret is on, if it is on one.
    pub caret: Option<usize>,
    /// Trace samples across the full width of the panel.
    pub samples_per_view: f32,
    /// The first sample drawn.
    pub scroll: f32,
    pub height: f32,
}

/// Bases whose call is worth a second look get a marked background.
const LOW_QUALITY: u8 = 20;

/// Below this much room per call, the letters are dropped.
///
/// The same rule the alignment canvas uses, and for the same reason: a glyph
/// narrower than this is unreadable, so painting one per call across a whole
/// 1000-base read is a thousand shapes nobody can read. The trace lines stay —
/// they are what the zoomed-out view is for.
const MIN_LABEL_WIDTH: f32 = 6.0;

impl Chromatogram<'_> {
    /// Draw the trace, returning what the user did to it.
    pub fn show(self, ui: &mut egui::Ui, actions: &mut Vec<TraceAction>) {
        let width = ui.available_width().max(1.0);
        let (rect, response) =
            ui.allocate_exact_size(Vec2::new(width, self.height), Sense::click());
        let painter = ui.painter_at(rect);
        let visuals = ui.visuals();
        painter.rect_filled(rect, 0.0, visuals.extreme_bg_color);

        if let Link::Lost = self.view.link() {
            painter.text(
                rect.center(),
                Align2::CENTER_CENTER,
                "The sequence no longer matches this trace, so the peaks cannot be lined up \
                 with it. Undo the edits, or reopen the trace.",
                FontId::proportional(13.0),
                visuals.warn_fg_color,
            );
            return;
        }

        let samples = self.view.trace.samples();
        if samples == 0 {
            return;
        }
        let per_sample = width / self.samples_per_view.max(1.0);
        let first = self.scroll.max(0.0);
        let last = (first + self.samples_per_view).min(samples as f32);
        let x_of = |sample: f32| rect.left() + (sample - first) * per_sample;

        // Room for the letters along the top.
        let label_height = 18.0;
        let plot = Rect::from_min_max(
            Pos2::new(rect.left(), rect.top() + label_height),
            rect.right_bottom(),
        );
        let scale =
            plot.height() / (self.view.trace.peak_signal() as f32 / self.view.gain.max(0.05));

        // The four channels, back to front so the tallest is not buried.
        for channel in 0..4 {
            let base = self.view.trace.channel_bases[channel];
            let color = base_color(base, visuals.dark_mode);
            let signal = &self.view.trace.channels[channel];
            let mut points: Vec<Pos2> = Vec::new();
            // One point per sample is more than a screen can show; step so the
            // cost is bounded by the width of the panel, not the length of the
            // read.
            let step = ((last - first) / width.max(1.0)).max(1.0);
            let mut at = first;
            while at <= last {
                let index = at as usize;
                let Some(&value) = signal.get(index) else { break };
                points.push(Pos2::new(
                    x_of(at),
                    plot.bottom() - (value as f32 * scale).min(plot.height()),
                ));
                at += step;
            }
            if points.len() > 1 {
                painter.add(egui::Shape::line(points, Stroke::new(1.2, color)));
            }
        }

        // The calls: the row's own residues, over the peaks they sit on. At a
        // zoom where they would be illegible they are dropped entirely, so the
        // cost of a frame is bounded by the width of the panel rather than by
        // the length of the read.
        let per_call = per_sample * self.view.trace.mean_peak_spacing();
        if per_call < MIN_LABEL_WIDTH {
            self.handle_click(&response, rect, first, per_sample, actions);
            return;
        }
        let font = FontId::monospace(13.0);
        let mut residue_index = 0usize;
        for &residue in self.residues {
            if tolviewer_core::is_gap(residue) {
                continue;
            }
            let index = residue_index;
            residue_index += 1;
            let Some(sample) = self.view.sample_for_residue(index) else { continue };
            let sample = sample as f32;
            if sample < first || sample > last {
                continue;
            }
            let x = x_of(sample);
            let cell = Rect::from_center_size(
                Pos2::new(x, rect.top() + label_height / 2.0),
                Vec2::new((per_sample * 12.0).clamp(8.0, 26.0), label_height),
            );

            // Flag the two things worth a second look: a call the instrument
            // was unsure of, and one a person has since changed.
            let quality = self.view.quality_for_residue(index);
            let edited = self
                .view
                .call_for_residue(index)
                .is_some_and(|called| !called.eq_ignore_ascii_case(&residue));
            if edited {
                painter.rect_filled(cell, 2.0, visuals.selection.bg_fill.gamma_multiply(0.6));
            } else if quality.is_some_and(|q| q < LOW_QUALITY) {
                painter.rect_filled(cell, 2.0, visuals.warn_fg_color.gamma_multiply(0.25));
            }
            if self.caret == Some(index) {
                painter.rect_stroke(
                    cell,
                    2.0,
                    Stroke::new(1.5, visuals.selection.stroke.color),
                    StrokeKind::Inside,
                );
            }
            painter.text(
                cell.center(),
                Align2::CENTER_CENTER,
                (residue as char).to_string(),
                font.clone(),
                base_color(residue, visuals.dark_mode),
            );
            // A tick from the letter down to its peak, so a crowded region is
            // still readable.
            painter.line_segment(
                [Pos2::new(x, cell.bottom()), Pos2::new(x, plot.top() + 3.0)],
                Stroke::new(1.0, visuals.weak_text_color()),
            );
        }

        self.handle_click(&response, rect, first, per_sample, actions);
    }

    /// Clicking picks the nearest call, which is what a person means when they
    /// click at a peak. Available at every zoom, including the one where the
    /// letters are not drawn.
    fn handle_click(
        &self,
        response: &egui::Response,
        rect: Rect,
        first: f32,
        per_sample: f32,
        actions: &mut Vec<TraceAction>,
    ) {
        if let Some(pos) = response.interact_pointer_pos() {
            let sample = first + (pos.x - rect.left()) / per_sample;
            if let Some(index) = nearest_residue(self.view, self.residues, sample) {
                actions.push(TraceAction::Select(index));
            }
        }
    }
}

/// The residue whose peak is closest to `sample`.
///
/// Distances are compared in `f64`: a pointer position far outside the trace
/// (which a drag off the edge of the panel produces) makes every `f32` distance
/// round to the same value, and the "nearest" peak becomes whichever was tested
/// first rather than the one at the end of the read.
fn nearest_residue(view: &TraceView, residues: &[u8], sample: f32) -> Option<usize> {
    let sample = f64::from(sample);
    let mut best: Option<(usize, f64)> = None;
    let mut index = 0usize;
    for &residue in residues {
        if tolviewer_core::is_gap(residue) {
            continue;
        }
        if let Some(at) = view.sample_for_residue(index) {
            let distance = (f64::from(at) - sample).abs();
            if best.is_none_or(|(_, d)| distance < d) {
                best = Some((index, distance));
            }
        }
        index += 1;
    }
    best.map(|(i, _)| i)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tolviewer_io::ab1;

    fn trace() -> Trace {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../testdata/ab1");
        ab1::read_file(&dir.join("tingidae_COI_F.ab1")).unwrap()
    }

    #[test]
    fn an_untouched_read_links_at_offset_zero() {
        let t = trace();
        let calls = t.calls.clone();
        let mut view = TraceView::new(t, 0);
        view.relink(&calls, 1);
        match view.link() {
            Link::At { offset, identity } => {
                assert_eq!(offset, 0);
                assert_eq!(identity, 1.0);
            }
            Link::Lost => panic!("an unedited read must line up with itself"),
        }
        assert_eq!(view.sample_for_residue(0), Some(view.trace.peaks[0]));
    }

    #[test]
    fn trimming_the_start_shifts_the_link_by_exactly_what_was_cut() {
        let t = trace();
        let trimmed = t.calls[40..].to_vec();
        let mut view = TraceView::new(t, 0);
        view.relink(&trimmed, 1);
        assert_eq!(view.link(), Link::At { offset: 40, identity: 1.0 });
        // Residue 0 of the trimmed row is call 40, and sits on its peak.
        assert_eq!(view.sample_for_residue(0), Some(view.trace.peaks[40]));
        assert_eq!(view.call_for_residue(0), Some(view.trace.calls[40]));
    }

    #[test]
    fn undoing_a_trim_puts_the_link_back() {
        let t = trace();
        let whole = t.calls.clone();
        let trimmed = whole[40..].to_vec();
        let mut view = TraceView::new(t, 0);
        view.relink(&trimmed, 1);
        assert!(matches!(view.link(), Link::At { offset: 40, .. }));
        // Undo restores the row; the link must follow it back without anything
        // having been remembered.
        view.relink(&whole, 2);
        assert!(matches!(view.link(), Link::At { offset: 0, .. }));
    }

    #[test]
    fn a_few_retyped_calls_do_not_break_the_link() {
        let t = trace();
        let mut edited = t.calls.clone();
        for i in [10, 50, 51, 200] {
            edited[i] = b'A';
        }
        let mut view = TraceView::new(t, 0);
        view.relink(&edited, 1);
        match view.link() {
            Link::At { offset, identity } => {
                assert_eq!(offset, 0);
                assert!(identity > 0.98, "{identity}");
            }
            Link::Lost => panic!("four edits must not lose the trace"),
        }
        // The panel can tell which residues no longer match the instrument.
        assert_eq!(view.call_for_residue(10), Some(t2_call(10)));
        fn t2_call(i: usize) -> u8 {
            trace().calls[i]
        }
    }

    #[test]
    fn an_unrelated_sequence_is_reported_lost_rather_than_lined_up_anyway() {
        let t = trace();
        let nonsense = vec![b'A'; 300];
        let mut view = TraceView::new(t, 0);
        view.relink(&nonsense, 1);
        assert_eq!(view.link(), Link::Lost);
        assert_eq!(view.sample_for_residue(0), None);
    }

    #[test]
    fn a_row_longer_than_the_trace_is_lost_not_a_panic() {
        let t = trace();
        let too_long = vec![b'A'; t.calls.len() + 10];
        let mut view = TraceView::new(t, 0);
        view.relink(&too_long, 1);
        assert_eq!(view.link(), Link::Lost);
        let mut view = TraceView::new(trace(), 0);
        view.relink(&[], 1);
        assert_eq!(view.link(), Link::Lost);
    }

    #[test]
    fn relinking_is_skipped_when_nothing_has_changed() {
        let t = trace();
        let calls = t.calls.clone();
        let mut view = TraceView::new(t, 0);
        view.relink(&calls, 7);
        // A different row at the same revision must not be recomputed: the
        // revision is the promise that nothing changed.
        view.relink(&[b'A'; 300], 7);
        assert!(matches!(view.link(), Link::At { offset: 0, .. }));
        view.relink(&[b'A'; 300], 8);
        assert_eq!(view.link(), Link::Lost);
    }

    #[test]
    fn the_reverse_complement_of_a_read_links_to_the_reversed_trace() {
        let mut t = trace();
        t.reverse_complement();
        let calls = t.calls.clone();
        let mut view = TraceView::new(t, 0);
        view.relink(&calls, 1);
        assert!(matches!(view.link(), Link::At { offset: 0, identity: 1.0 }));
    }

    #[test]
    fn clicking_finds_the_call_nearest_the_pointer() {
        let t = trace();
        let calls = t.calls.clone();
        let peaks = t.peaks.clone();
        let mut view = TraceView::new(t, 0);
        view.relink(&calls, 1);
        // Exactly on the peak of residue 30, and slightly off it.
        assert_eq!(nearest_residue(&view, &calls, peaks[30] as f32), Some(30));
        assert_eq!(nearest_residue(&view, &calls, peaks[30] as f32 + 2.0), Some(30));
        assert_eq!(nearest_residue(&view, &calls, 0.0), Some(0));
        // Dragging off the end of the panel must still land on the last call,
        // not on whichever peak happened to be measured first.
        assert_eq!(nearest_residue(&view, &calls, 1.0e9), Some(calls.len() - 1));
        assert_eq!(nearest_residue(&view, &calls, -1.0e9), Some(0));
    }

    #[test]
    fn gaps_in_the_row_do_not_consume_calls() {
        let t = trace();
        let mut gapped: Vec<u8> = Vec::new();
        for (i, &c) in t.calls.iter().enumerate() {
            if i == 5 {
                gapped.push(b'-');
            }
            gapped.push(c);
        }
        let ungapped: Vec<u8> = gapped.iter().copied().filter(|&c| c != b'-').collect();
        let peaks = t.peaks.clone();
        let mut view = TraceView::new(t, 0);
        view.relink(&ungapped, 1);
        // Residue 5 is still call 5, even though column 5 is a gap.
        assert_eq!(view.sample_for_residue(5), Some(peaks[5]));
        assert_eq!(nearest_residue(&view, &gapped, peaks[5] as f32), Some(5));
    }
}
