//! Residue colouring.
//!
//! The canvas paints a coloured background behind each cell and draws the
//! letter on top in a contrasting ink, so every scheme here returns a
//! background colour and lets [`ink_for`] pick the text colour.

use egui::Color32;
use tolviewer_core::{is_gap, Alphabet, ColumnStats};

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum ColorScheme {
    /// Nucleotides by base, amino acids by physicochemical class.
    Residue,
    /// Clustal X's protein scheme; falls back to `Residue` for nucleotides.
    Clustal,
    /// Colour only where a row differs from the consensus.
    Differences,
    /// Shade by how conserved the column is.
    Conservation,
    /// Shade by Phred quality where known.
    Quality,
    /// No colour; letters only.
    None,
}

impl ColorScheme {
    pub const ALL: &'static [ColorScheme] = &[
        ColorScheme::Residue,
        ColorScheme::Clustal,
        ColorScheme::Differences,
        ColorScheme::Conservation,
        ColorScheme::Quality,
        ColorScheme::None,
    ];

    pub fn name(self) -> &'static str {
        match self {
            ColorScheme::Residue => "Residue",
            ColorScheme::Clustal => "Clustal",
            ColorScheme::Differences => "Differences from consensus",
            ColorScheme::Conservation => "Conservation",
            ColorScheme::Quality => "Quality (Phred)",
            ColorScheme::None => "None",
        }
    }

    /// Needs a consensus to be computed before painting.
    pub fn needs_consensus(self) -> bool {
        matches!(self, ColorScheme::Differences | ColorScheme::Conservation)
    }
}

/// Everything the painter needs to colour one cell.
pub struct CellContext<'a> {
    pub residue: u8,
    pub alphabet: Alphabet,
    pub consensus: Option<u8>,
    pub stats: Option<&'a ColumnStats>,
    pub quality: Option<u8>,
    pub dark: bool,
}

const GAP_BG_LIGHT: Color32 = Color32::from_rgb(0xF2, 0xF2, 0xF4);
const GAP_BG_DARK: Color32 = Color32::from_rgb(0x2A, 0x2C, 0x31);

/// Background colour for a cell under `scheme`, or `None` to leave it bare.
pub fn background(scheme: ColorScheme, cx: &CellContext<'_>) -> Option<Color32> {
    let c = cx.residue.to_ascii_uppercase();
    if is_gap(cx.residue) {
        return Some(if cx.dark { GAP_BG_DARK } else { GAP_BG_LIGHT });
    }
    let base = match scheme {
        ColorScheme::None => return None,
        ColorScheme::Residue => {
            if cx.alphabet.is_nucleotide() {
                nucleotide_color(c)
            } else {
                amino_class_color(c)
            }
        }
        ColorScheme::Clustal => {
            if cx.alphabet.is_nucleotide() {
                nucleotide_color(c)
            } else {
                clustal_color(c)
            }
        }
        ColorScheme::Differences => {
            let consensus = cx.consensus?;
            if consensus.to_ascii_uppercase() == c {
                return None;
            }
            // A mismatch is what the eye should catch, so give it a warm wash
            // regardless of which residue it is.
            Color32::from_rgb(0xE8, 0x7A, 0x6B)
        }
        ColorScheme::Conservation => {
            let q = cx.stats?.identity();
            return Some(ramp(q, cx.dark));
        }
        ColorScheme::Quality => {
            let q = cx.quality? as f32;
            // Phred 0..40 is the useful range for Sanger and Illumina alike.
            return Some(ramp((q / 40.0).clamp(0.0, 1.0), cx.dark));
        }
    };
    Some(if cx.dark { dim(base) } else { base })
}

/// The colour a single base is drawn in, for the chromatogram's trace lines
/// and its letters.
///
/// This is the same palette the canvas paints cells with, so a base is the same
/// colour wherever it appears. Unlike [`background`] it returns a line colour,
/// which stays saturated in dark mode rather than being dimmed to sit behind
/// text.
pub fn base_color(residue: u8, dark: bool) -> Color32 {
    let c = nucleotide_color(residue.to_ascii_uppercase());
    if dark {
        // Lines are thin; lift them slightly so they read against a dark panel.
        Color32::from_rgb(
            c.r().saturating_add(0x18),
            c.g().saturating_add(0x18),
            c.b().saturating_add(0x18),
        )
    } else {
        c
    }
}

/// Text colour that stays readable on `bg`.
pub fn ink_for(bg: Option<Color32>, dark: bool) -> Color32 {
    match bg {
        None => {
            if dark {
                Color32::from_gray(0xE0)
            } else {
                Color32::from_gray(0x20)
            }
        }
        Some(bg) => {
            // Rec. 601 luma is good enough to choose between two inks.
            let luma = 0.299 * bg.r() as f32 + 0.587 * bg.g() as f32 + 0.114 * bg.b() as f32;
            if luma > 140.0 {
                Color32::from_gray(0x18)
            } else {
                Color32::from_gray(0xF0)
            }
        }
    }
}

fn nucleotide_color(c: u8) -> Color32 {
    match c {
        b'A' => Color32::from_rgb(0x6F, 0xC2, 0x76),
        b'C' => Color32::from_rgb(0x6C, 0xA6, 0xE0),
        b'G' => Color32::from_rgb(0xF2, 0xC0, 0x5C),
        b'T' | b'U' => Color32::from_rgb(0xE8, 0x7A, 0x6B),
        b'N' | b'?' | b'X' => Color32::from_gray(0xC8),
        // IUPAC ambiguity codes: a muted violet, distinct from the four bases.
        _ => Color32::from_rgb(0xB9, 0xA6, 0xD6),
    }
}

/// Amino acids grouped by physicochemical class.
fn amino_class_color(c: u8) -> Color32 {
    match c {
        // hydrophobic
        b'A' | b'V' | b'L' | b'I' | b'M' | b'F' | b'W' | b'C' => {
            Color32::from_rgb(0x6C, 0xA6, 0xE0)
        }
        // polar uncharged
        b'S' | b'T' | b'N' | b'Q' => Color32::from_rgb(0x6F, 0xC2, 0x76),
        // positively charged
        b'K' | b'R' | b'H' => Color32::from_rgb(0xE8, 0x7A, 0x6B),
        // negatively charged
        b'D' | b'E' => Color32::from_rgb(0xC9, 0x7B, 0xD6),
        // special
        b'G' => Color32::from_rgb(0xF2, 0xA5, 0x5C),
        b'P' => Color32::from_rgb(0xF2, 0xD9, 0x5C),
        b'Y' => Color32::from_rgb(0x5C, 0xC7, 0xC2),
        _ => Color32::from_gray(0xC8),
    }
}

/// The Clustal X protein palette.
fn clustal_color(c: u8) -> Color32 {
    match c {
        b'A' | b'I' | b'L' | b'M' | b'F' | b'W' | b'V' | b'C' => {
            Color32::from_rgb(0x80, 0xA0, 0xF0)
        }
        b'K' | b'R' => Color32::from_rgb(0xF0, 0x15, 0x05),
        b'E' | b'D' => Color32::from_rgb(0xC0, 0x48, 0xC0),
        b'N' | b'Q' | b'S' | b'T' => Color32::from_rgb(0x15, 0xC0, 0x15),
        b'G' => Color32::from_rgb(0xF0, 0x90, 0x48),
        b'P' => Color32::from_rgb(0xC0, 0xC0, 0x00),
        b'H' | b'Y' => Color32::from_rgb(0x15, 0xA4, 0xA4),
        _ => Color32::from_gray(0xC8),
    }
}

/// 0.0 (poor, red) -> 0.5 (amber) -> 1.0 (good, green).
fn ramp(t: f32, dark: bool) -> Color32 {
    let t = t.clamp(0.0, 1.0);
    let (r, g, b) = if t < 0.5 {
        let k = t * 2.0;
        (0xE0 as f32, 0x6A as f32 + k * (0xC0 as f32 - 0x6A as f32), 0x5A as f32)
    } else {
        let k = (t - 0.5) * 2.0;
        (
            0xE0 as f32 * (1.0 - k) + 0x6F as f32 * k,
            0xC0 as f32,
            0x5A as f32 * (1.0 - k) + 0x76 as f32 * k,
        )
    };
    let c = Color32::from_rgb(r as u8, g as u8, b as u8);
    if dark {
        dim(c)
    } else {
        c
    }
}

/// Darken a light-theme swatch so it works as a background in dark mode
/// without washing out the letter drawn on top.
fn dim(c: Color32) -> Color32 {
    Color32::from_rgb(
        (c.r() as f32 * 0.55) as u8,
        (c.g() as f32 * 0.55) as u8,
        (c.b() as f32 * 0.55) as u8,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cx(residue: u8) -> CellContext<'static> {
        CellContext {
            residue,
            alphabet: Alphabet::Dna,
            consensus: Some(b'A'),
            stats: None,
            quality: None,
            dark: false,
        }
    }

    #[test]
    fn bases_get_distinct_colours() {
        let a = background(ColorScheme::Residue, &cx(b'A')).unwrap();
        let c = background(ColorScheme::Residue, &cx(b'C')).unwrap();
        let g = background(ColorScheme::Residue, &cx(b'G')).unwrap();
        let t = background(ColorScheme::Residue, &cx(b'T')).unwrap();
        for (x, y) in [(a, c), (a, g), (a, t), (c, g), (c, t), (g, t)] {
            assert_ne!(x, y);
        }
    }

    #[test]
    fn lowercase_colours_like_uppercase() {
        assert_eq!(
            background(ColorScheme::Residue, &cx(b'a')),
            background(ColorScheme::Residue, &cx(b'A'))
        );
    }

    #[test]
    fn gaps_are_neutral_in_every_scheme() {
        for &s in ColorScheme::ALL {
            assert!(background(s, &cx(b'-')).is_some(), "{s:?} left a gap unpainted");
        }
    }

    #[test]
    fn differences_scheme_only_paints_mismatches() {
        assert!(background(ColorScheme::Differences, &cx(b'A')).is_none());
        assert!(background(ColorScheme::Differences, &cx(b'C')).is_some());
    }

    #[test]
    fn ink_contrasts_with_its_background() {
        let dark_bg = Color32::from_rgb(0x20, 0x20, 0x20);
        let light_bg = Color32::from_rgb(0xF0, 0xF0, 0xF0);
        assert!(ink_for(Some(dark_bg), false).r() > 0x80);
        assert!(ink_for(Some(light_bg), false).r() < 0x80);
    }
}
