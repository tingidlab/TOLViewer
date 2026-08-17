//! Alignment cleaning for TOLViewer: a native-Rust reimplementation of the
//! Gblocks column-selection algorithm, plus the simple gap filters the GUI
//! offers alongside it.
//!
//! ```
//! use tolviewer_clean::{gblocks, GblocksParams};
//! use tolviewer_core::{Alignment, Sequence};
//!
//! let aln = Alignment::new(
//!     "demo",
//!     vec![
//!         Sequence::new("a", *b"ACGTACGTACGT"),
//!         Sequence::new("b", *b"ACGTACGTACGT"),
//!         Sequence::new("c", *b"ACGTACGTACGT"),
//!     ],
//! );
//! let result = gblocks(&aln, &GblocksParams::defaults(aln.len()))?;
//! assert_eq!(result.mask_line(), "############");
//! assert_eq!(result.apply(&aln)?.width(), result.kept);
//! # Ok::<(), tolviewer_core::Error>(())
//! ```
//!
//! # References
//!
//! * Castresana J. (2000) "Selection of conserved blocks from multiple
//!   alignments for their use in phylogenetic analysis." *Molecular Biology
//!   and Evolution* 17(4):540-552.
//! * Castresana J., *Gblocks documentation* (versions 0.91b / 1.0), which
//!   describes what the shipped program does and is the behaviour this crate
//!   matches where it differs from the paper.
//! * Talavera G. & Castresana J. (2007) "Improvement of phylogenies after
//!   removing divergent and ambiguously aligned blocks from protein sequence
//!   alignments." *Systematic Biology* 56(4):564-577, for the relaxed
//!   parameter set exposed as [`GblocksParams::relaxed`].
//!
//! The implementation is original Rust written from those descriptions.

#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod filters;
mod gblocks;
mod similarity;

pub use filters::{remove_gappy_columns, remove_gappy_sequences, trim_ends};
pub use gblocks::{gblocks, ColumnFlag, GapPolicy, GblocksParams, GblocksResult};
