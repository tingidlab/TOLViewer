//! TOLViewer's application layer.
//!
//! The binary in `main.rs` is a thin wrapper around this crate so the document
//! model, selection logic and canvas geometry can be tested without opening a
//! window.

#![forbid(unsafe_code)]

pub mod app;
pub mod canvas;
pub mod document;
pub mod selection;
pub mod tasks;
pub mod theme;
pub mod ui;

pub use app::TolViewerApp;
pub use document::Document;
