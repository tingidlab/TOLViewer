#![forbid(unsafe_code)]
//! The TOLViewer project library: a tree of folders over sequence files that
//! stay where they are.
//!
//! A working phylogenetics project is dozens of files — a plate of traces per
//! locus, the alignments built from them, the supermatrix built from those —
//! and the thing a lab actually needs is somewhere to keep them straight:
//!
//! ```text
//! Lace bug project
//!   18S
//!     TL-2213_18S_F.ab1
//!     TL-2213_18S_R.ab1
//!   28S
//!     ...
//! ```
//!
//! This crate is that tree, plus the operations that go with it:
//!
//! * [`library`] — the folders, the entries, and the rule that decides whether
//!   saving an edit may overwrite the lab's file or must go to a copy.
//! * [`primer`] — mapping PCR primers onto reads and trimming back to the
//!   amplicon.
//! * [`concat`] — joining per-locus alignments into a supermatrix, matching the
//!   same specimen across loci by name.
//! * [`naming`] — the name matching that makes concatenation work.
//! * [`store`] — the `.tolvlib` file the tree is saved to.
//!
//! ## Nothing is written without being asked
//!
//! Adding a file to a library reads it once and remembers where it is. The
//! files are the sequencing facility's output and the lab's records, so an
//! edit that would land on one raises [`SaveTarget::Original`], which the GUI
//! turns into a question; answering it with [`SaveChoice::NewCopy`] diverts the
//! edit to a copy and *remembers the copy*, so the question is asked once per
//! entry rather than once per save.
//!
//! ```no_run
//! # fn main() -> tolviewer_core::Result<()> {
//! # use std::path::Path;
//! use tolviewer_io::WriteOptions;
//! use tolviewer_library::{Library, SaveChoice, SaveTarget};
//!
//! let mut library = Library::new("Lace bug project");
//! let folder = library.add_folder(None, "18S")?;
//! let read = library.add_file(Some(folder), Path::new("reads/TL-2213_18S_F.ab1"))?;
//!
//! let mut edited = library.load(read)?;
//! edited.sequences[0].residues.truncate(600);
//!
//! let choice = match library.save_target(read)? {
//!     SaveTarget::WorkingCopy(_) => SaveChoice::Overwrite, // already a copy
//!     target => SaveChoice::NewCopy(target.path().to_path_buf()), // ask first
//! };
//! library.save_entry(read, &edited, choice, &WriteOptions::default())?;
//! # Ok(()) }
//! ```

pub mod concat;
pub mod library;
pub mod naming;
pub mod primer;
pub mod store;

pub use concat::{concatenate, ConcatOptions, ConcatResult, Partition, SamplePreview};
pub use library::{
    CopyReason, Entry, EntryKind, Folder, Library, Node, NodeId, NodeKind, SaveChoice, SaveTarget,
};
pub use naming::{sample_key, MatchOptions};
pub use primer::{plan_trim, Primer, PrimerHit, PrimerSet, Strand, TrimOptions, TrimPlan};
