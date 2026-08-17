//! Core data model for TOLViewer: sequences, alignments, editing and statistics.
//!
//! Everything downstream (I/O, alignment engines, cleaning, GUI) is written
//! against the types in this crate, so it deliberately has no dependencies.

pub mod alphabet;
pub mod alignment;
pub mod edit;
pub mod error;
pub mod sequence;
pub mod stats;

pub use alphabet::{Alphabet, GAP, GAP_CHARS, MISSING};
pub use alignment::Alignment;
pub use edit::{EditOp, UndoStack};
pub use error::{Error, Result};
pub use sequence::Sequence;
pub use stats::{ColumnStats, Consensus};
