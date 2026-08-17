//! Background jobs.
//!
//! Aligning a few hundred sequences takes seconds to minutes, so it must not
//! run on the UI thread. Each job owns a clone of the input, reports progress
//! through an atomic, and can be cancelled from the UI.

use std::sync::atomic::{AtomicBool, AtomicU32, Ordering};
use std::sync::mpsc::{self, Receiver, TryRecvError};
use std::sync::{Arc, Mutex};

use tolviewer_align::{AlignParams, Progress};
use tolviewer_clean::{GblocksParams, GblocksResult};
use tolviewer_core::{Alignment, Error, Result};

/// Shared progress state. Cheap enough to poll every frame.
pub struct TaskProgress {
    /// Fraction 0..=1, stored as the bits of an f32.
    fraction: AtomicU32,
    message: Mutex<String>,
    cancel: AtomicBool,
    /// Waking the UI thread from the worker keeps the progress bar moving even
    /// when the user is not touching the mouse.
    ctx: egui::Context,
}

impl TaskProgress {
    fn new(ctx: egui::Context, message: &str) -> Self {
        TaskProgress {
            fraction: AtomicU32::new(0f32.to_bits()),
            message: Mutex::new(message.to_string()),
            cancel: AtomicBool::new(false),
            ctx,
        }
    }

    pub fn fraction(&self) -> f32 {
        f32::from_bits(self.fraction.load(Ordering::Relaxed))
    }

    pub fn message(&self) -> String {
        self.message.lock().map(|m| m.clone()).unwrap_or_default()
    }

    pub fn request_cancel(&self) {
        self.cancel.store(true, Ordering::Relaxed);
    }

    pub fn is_cancelled(&self) -> bool {
        self.cancel.load(Ordering::Relaxed)
    }
}

impl Progress for TaskProgress {
    fn tick(&self, fraction: f32, message: &str) -> bool {
        self.fraction.store(fraction.clamp(0.0, 1.0).to_bits(), Ordering::Relaxed);
        if let Ok(mut m) = self.message.lock() {
            if *m != message {
                m.clear();
                m.push_str(message);
            }
        }
        self.ctx.request_repaint();
        !self.is_cancelled()
    }
}

/// What a finished job produced.
pub enum TaskOutcome {
    /// A new alignment to install as one undoable step, with a label for the
    /// undo menu.
    Alignment { label: String, alignment: Box<Alignment> },
    /// A cleaning result: the mask is shown as a track and only applied when
    /// the user confirms.
    Clean(Box<GblocksResult>),
}

/// A job in flight.
pub struct Task {
    /// Which document it belongs to; the job is dropped if that document is
    /// closed before it finishes.
    pub doc: usize,
    pub label: String,
    pub progress: Arc<TaskProgress>,
    rx: Receiver<Result<TaskOutcome>>,
    finished: bool,
}

impl Task {
    /// Poll for completion. Returns `Some` exactly once per task.
    pub fn poll(&mut self) -> Option<Result<TaskOutcome>> {
        if self.finished {
            return None;
        }
        match self.rx.try_recv() {
            Ok(result) => {
                self.finished = true;
                Some(result)
            }
            Err(TryRecvError::Empty) => None,
            Err(TryRecvError::Disconnected) => {
                self.finished = true;
                // The worker panicked or was dropped without sending.
                Some(Err(Error::algorithm("the background job stopped unexpectedly")))
            }
        }
    }

    pub fn is_finished(&self) -> bool {
        self.finished
    }
}

fn spawn<F>(ctx: &egui::Context, doc: usize, label: &str, job: F) -> Task
where
    F: FnOnce(&TaskProgress) -> Result<TaskOutcome> + Send + 'static,
{
    let progress = Arc::new(TaskProgress::new(ctx.clone(), label));
    let (tx, rx) = mpsc::channel();
    let worker_progress = Arc::clone(&progress);
    let ctx = ctx.clone();
    std::thread::Builder::new()
        .name(format!("tolviewer-{label}"))
        .spawn(move || {
            let result = job(&worker_progress);
            // A closed receiver just means the document went away.
            let _ = tx.send(result);
            ctx.request_repaint();
        })
        .expect("the OS refused to start a worker thread");
    Task { doc, label: label.to_string(), progress, rx, finished: false }
}

/// Align every sequence in `alignment`.
pub fn align(ctx: &egui::Context, doc: usize, alignment: Alignment, params: AlignParams) -> Task {
    let label = format!("align ({})", params.engine.name());
    spawn(ctx, doc, &label, move |p| {
        let aligned = tolviewer_align::align(&alignment, &params, p)?;
        Ok(TaskOutcome::Alignment {
            label: format!("align with {}", params.engine.name()),
            alignment: Box::new(aligned),
        })
    })
}

/// Re-align only `cols`, leaving the rest of the alignment untouched.
pub fn realign_region(
    ctx: &egui::Context,
    doc: usize,
    alignment: Alignment,
    cols: std::ops::Range<usize>,
    params: AlignParams,
) -> Task {
    let label = "realign selection".to_string();
    spawn(ctx, doc, &label, move |p| {
        let aligned = tolviewer_align::realign_region(&alignment, cols, &params, p)?;
        Ok(TaskOutcome::Alignment {
            label: "realign selection".to_string(),
            alignment: Box::new(aligned),
        })
    })
}

/// Run Gblocks. Fast enough to be synchronous on small alignments, but big
/// ones deserve the same cancellable treatment as alignment.
pub fn clean(ctx: &egui::Context, doc: usize, alignment: Alignment, params: GblocksParams) -> Task {
    spawn(ctx, doc, "clean", move |p| {
        p.tick(0.1, "selecting conserved blocks");
        let result = tolviewer_clean::gblocks(&alignment, &params)?;
        p.tick(1.0, "done");
        Ok(TaskOutcome::Clean(Box::new(result)))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn progress_round_trips_and_cancels() {
        let ctx = egui::Context::default();
        let p = TaskProgress::new(ctx, "start");
        assert!(p.tick(0.25, "quarter"));
        assert!((p.fraction() - 0.25).abs() < 1e-6);
        assert_eq!(p.message(), "quarter");
        p.request_cancel();
        assert!(!p.tick(0.5, "half"), "tick must report the cancellation");
    }

    #[test]
    fn progress_fraction_is_clamped() {
        let p = TaskProgress::new(egui::Context::default(), "x");
        p.tick(5.0, "over");
        assert_eq!(p.fraction(), 1.0);
        p.tick(-1.0, "under");
        assert_eq!(p.fraction(), 0.0);
    }

    #[test]
    fn a_task_yields_its_result_exactly_once() {
        let ctx = egui::Context::default();
        let mut task = spawn(&ctx, 0, "test", |p| {
            p.tick(1.0, "done");
            Ok(TaskOutcome::Alignment {
                label: "test".into(),
                alignment: Box::new(Alignment::default()),
            })
        });
        let mut result = None;
        for _ in 0..2000 {
            if let Some(r) = task.poll() {
                result = Some(r);
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(1));
        }
        assert!(result.is_some(), "the task never finished");
        assert!(task.poll().is_none(), "a task must only report once");
        assert!(task.is_finished());
    }
}
