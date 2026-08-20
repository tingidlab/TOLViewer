//! TOLViewer — view, edit, align and clean DNA and protein alignments.

#![forbid(unsafe_code)]
// On Windows, keep a console window from appearing behind the GUI in release
// builds while leaving it available for `--help` in debug builds.
#![cfg_attr(all(not(debug_assertions), target_os = "windows"), windows_subsystem = "windows")]

use std::path::PathBuf;

use tolviewer_app::TolViewerApp;

const HELP: &str = "\
TOLViewer — view, edit, align and clean DNA and protein alignments

USAGE:
    tolviewer [FILE]...

Any FASTA, FASTQ, PHYLIP, NEXUS, Clustal, Stockholm, MSF, GenBank or AB1 file
given on the command line is opened in a tab; an AB1 trace brings its
chromatogram with it. A .tolvlib file is opened as a project library instead.
With no arguments, TOLViewer starts empty and reopens the library you had open
last.

OPTIONS:
    -h, --help       print this message
    -V, --version    print the version
";

fn main() -> eframe::Result<()> {
    let mut paths: Vec<PathBuf> = Vec::new();
    for arg in std::env::args().skip(1) {
        match arg.as_str() {
            "-h" | "--help" => {
                println!("{HELP}");
                return Ok(());
            }
            "-V" | "--version" => {
                println!("tolviewer {}", env!("CARGO_PKG_VERSION"));
                return Ok(());
            }
            other => paths.push(PathBuf::from(other)),
        }
    }

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("TOLViewer")
            .with_inner_size([1280.0, 800.0])
            .with_min_inner_size([720.0, 460.0])
            .with_app_id("org.tingidlab.tolviewer")
            .with_icon(icon()),
        ..Default::default()
    };

    eframe::run_native(
        "TOLViewer",
        options,
        Box::new(move |cc| Ok(Box::new(TolViewerApp::new(cc, paths)))),
    )
}

/// A small procedurally drawn icon: stacked bars in the base colours, so the
/// app is recognisable in a dock or task bar without shipping a binary asset.
fn icon() -> egui::IconData {
    const SIZE: usize = 64;
    let palette = [
        [0x6F, 0xC2, 0x76, 0xFF],
        [0x6C, 0xA6, 0xE0, 0xFF],
        [0xF2, 0xC0, 0x5C, 0xFF],
        [0xE8, 0x7A, 0x6B, 0xFF],
    ];
    let mut rgba = vec![0u8; SIZE * SIZE * 4];
    for y in 0..SIZE {
        for x in 0..SIZE {
            let band = (y * 8 / SIZE).min(7);
            let color = if (4..60).contains(&x) && (4..60).contains(&y) && band.is_multiple_of(2) {
                palette[(band / 2) % palette.len()]
            } else {
                [0x1E, 0x20, 0x24, 0xFF]
            };
            let i = (y * SIZE + x) * 4;
            rgba[i..i + 4].copy_from_slice(&color);
        }
    }
    egui::IconData { rgba, width: SIZE as u32, height: SIZE as u32 }
}
