#![forbid(unsafe_code)]
//! Readers and writers for the sequence and alignment formats a phylogenetics
//! lab actually has on disk: FASTA, FASTQ, PHYLIP (strict and relaxed), NEXUS,
//! Clustal, Stockholm, GCG MSF, GenBank and Applied Biosystems `.ab1` traces.
//!
//! Everything funnels through [`tolviewer_core::Alignment`]:
//!
//! ```no_run
//! # fn main() -> tolviewer_core::Result<()> {
//! use std::path::Path;
//! use tolviewer_io::{read_file, write_file, Format, WriteOptions};
//!
//! let aln = read_file(Path::new("in.fasta"))?;
//! write_file(&aln, Path::new("out.nex"), Format::Nexus, &WriteOptions::default())?;
//! # Ok(()) }
//! ```
//!
//! The readers are deliberately permissive (CRLF, blank lines, wrapped or
//! interleaved blocks, digits and whitespace inside sequence lines, lowercase
//! residues) and never panic on malformed input: they return
//! [`tolviewer_core::Error::Parse`] with a line number wherever one makes sense.
//!
//! ## What survives a write
//!
//! * FASTA and FASTQ keep descriptions; every other format keeps only the id.
//! * FASTQ is the only format that stores quality.
//! * Strict PHYLIP truncates names to 10 characters (uniqueness is preserved by
//!   replacing the tail with a counter).
//! * Stockholm `#=G*` annotation and NEXUS blocks other than `data` are dropped
//!   on read, so they cannot be written back.
//! * MSF, GenBank and AB1 are read-only ([`Format::can_write`] is false).
//!
//! ## Traces
//!
//! Sanger `.ab1` files carry a chromatogram as well as a sequence. Reading one
//! through [`read_file`] gives you the base calls as a one-row alignment;
//! [`ab1::read_file`] gives you the signal underneath them as well, which is
//! what the trace viewer needs to let an operator vet a call.

pub mod ab1;
mod clustal;
mod fasta;
mod fastq;
mod format;
mod genbank;
mod msf;
mod nexus;
mod options;
mod phylip;
mod stockholm;
mod util;

use std::fs;
use std::path::{Path, PathBuf};

use tolviewer_core::{Alignment, Error, Result};

pub use format::Format;
pub use options::{LineEnding, WriteOptions};

/// Read a file, detecting the format from content then extension.
pub fn read_file(path: &Path) -> Result<Alignment> {
    let bytes = fs::read(path)?;
    let format = detect(&bytes, path)?;
    parse(&bytes, format, &stem(path))
}

/// The format [`read_file`] would use for `path`, without reading the whole
/// file into an alignment.
///
/// Callers that need to record what a file is — a library cataloguing files it
/// does not own, say — want this rather than a second guess at the extension.
pub fn sniff_file(path: &Path) -> Result<Format> {
    let mut file = fs::File::open(path)?;
    let mut head = vec![0u8; 64 * 1024];
    let read = std::io::Read::read(&mut file, &mut head)?;
    head.truncate(read);
    detect(&head, path)
}

/// Content first, extension second, and a clear complaint if neither settles it.
fn detect(head: &[u8], path: &Path) -> Result<Format> {
    let head = &head[..head.len().min(64 * 1024)];
    Format::sniff(head).or_else(|| Format::from_path(path)).ok_or_else(|| {
        Error::parse(
            "file",
            None,
            format!(
                "cannot tell what format '{}' is: it matches no known \
                     signature and its extension is not one we recognise",
                path.display()
            ),
        )
    })
}

/// Read a file with an explicitly chosen format.
pub fn read_file_as(path: &Path, format: Format) -> Result<Alignment> {
    let bytes = fs::read(path)?;
    parse(&bytes, format, &stem(path))
}

/// Parse in-memory bytes. `name` becomes [`Alignment::name`].
pub fn parse(bytes: &[u8], format: Format, name: &str) -> Result<Alignment> {
    match format {
        Format::Fasta => fasta::parse(bytes, name),
        Format::Fastq => fastq::parse(bytes, name),
        // PHYLIP files are frequently mislabelled; if the requested variant
        // does not fit, try the other one before giving up.
        Format::Phylip => phylip::parse(bytes, name, true)
            .or_else(|e| phylip::parse(bytes, name, false).map_err(|_| e)),
        Format::PhylipRelaxed => phylip::parse(bytes, name, false)
            .or_else(|e| phylip::parse(bytes, name, true).map_err(|_| e)),
        Format::Nexus => nexus::parse(bytes, name),
        Format::Clustal => clustal::parse(bytes, name),
        Format::Stockholm => stockholm::parse(bytes, name),
        Format::Msf => msf::parse(bytes, name),
        Format::Genbank => genbank::parse(bytes, name),
        // The calls become an ordinary one-row alignment. The chromatogram is
        // dropped here; `ab1::parse` keeps it.
        Format::Ab1 => {
            let trace = ab1::parse(bytes, name)?;
            let id = if trace.sample_name.is_empty() { name } else { &trace.sample_name };
            Ok(Alignment::new(name, vec![trace.to_sequence(id)]))
        }
    }
}

/// Render an alignment in `format`.
///
/// Fails with [`Error::Format`] when the data cannot be represented, e.g. a
/// ragged (unaligned) set in a matrix format, or a format that cannot be
/// written at all.
pub fn write_string(aln: &Alignment, format: Format, opts: &WriteOptions) -> Result<String> {
    match format {
        Format::Fasta => fasta::write(aln, opts),
        Format::Fastq => fastq::write(aln, opts),
        Format::Phylip => phylip::write(aln, opts, true),
        Format::PhylipRelaxed => phylip::write(aln, opts, false),
        Format::Nexus => nexus::write(aln, opts),
        Format::Clustal => clustal::write(aln, opts),
        Format::Stockholm => stockholm::write(aln, opts),
        Format::Msf | Format::Genbank | Format::Ab1 => {
            Err(Error::format(format!("{} files can be read but not written", format.name())))
        }
    }
}

/// Write an alignment to `path`, atomically.
///
/// The bytes go to a hidden temporary file in the destination directory which
/// is then renamed over the target, so an interrupted or failing write never
/// leaves a half-written file where the user's data was.
pub fn write_file(aln: &Alignment, path: &Path, format: Format, opts: &WriteOptions) -> Result<()> {
    let text = write_string(aln, format, opts)?;
    let tmp = temp_path(path);
    match fs::write(&tmp, text.as_bytes()).and_then(|()| fs::rename(&tmp, path)) {
        Ok(()) => Ok(()),
        Err(e) => {
            let _ = fs::remove_file(&tmp);
            Err(Error::Io(e))
        }
    }
}

/// `dir/.name.tolviewer-tmp`, next to the destination so the rename is atomic.
fn temp_path(path: &Path) -> PathBuf {
    let file = path
        .file_name()
        .map(|f| f.to_string_lossy().into_owned())
        .unwrap_or_else(|| "out".to_string());
    let dir = path.parent().unwrap_or_else(|| Path::new("."));
    dir.join(format!(".{file}.tolviewer-tmp"))
}

/// The file stem, used as the alignment's display name.
fn stem(path: &Path) -> String {
    path.file_stem().map(|s| s.to_string_lossy().into_owned()).unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;
    use tolviewer_core::Sequence;

    fn aln() -> Alignment {
        Alignment::new(
            "t",
            vec![Sequence::new("alpha", *b"ACGTACGTAC"), Sequence::new("beta", *b"ACGT--GTAC")],
        )
    }

    #[test]
    fn every_writable_format_writes_and_reads_back() {
        let a = aln();
        for &f in Format::all() {
            if !f.can_write() {
                continue;
            }
            let text = write_string(&a, f, &WriteOptions::default()).unwrap();
            let back = parse(text.as_bytes(), f, "t").unwrap();
            assert_eq!(back.len(), 2, "{}", f.name());
            assert_eq!(back.sequences[0].id, "alpha", "{}", f.name());
            assert_eq!(back.sequences[1].residues, b"ACGT--GTAC", "{}", f.name());
        }
    }

    #[test]
    fn read_only_formats_refuse_to_write() {
        for f in [Format::Msf, Format::Genbank, Format::Ab1] {
            let e = write_string(&aln(), f, &WriteOptions::default()).unwrap_err();
            assert!(matches!(e, Error::Format(_)));
        }
    }

    #[test]
    fn sniff_file_agrees_with_read_file() {
        let dir = std::env::temp_dir().join("tolviewer-io-sniff");
        fs::create_dir_all(&dir).unwrap();
        for (name, bytes, expected) in [
            ("a.fasta", &b">a\nACGT\n"[..], Format::Fasta),
            // Content wins over a misleading extension.
            ("b.fasta", &b"#NEXUS\nbegin data;\n"[..], Format::Nexus),
        ] {
            let p = dir.join(name);
            fs::write(&p, bytes).unwrap();
            assert_eq!(sniff_file(&p).unwrap(), expected, "{name}");
            let _ = fs::remove_file(&p);
        }
        let p = dir.join("mystery.dat");
        fs::write(&p, b"nothing recognisable here\n").unwrap();
        assert!(sniff_file(&p).is_err());
        let _ = fs::remove_file(&p);
    }

    #[test]
    fn temp_path_sits_next_to_the_target() {
        let p = temp_path(Path::new("/data/runs/out.fasta"));
        assert_eq!(p, PathBuf::from("/data/runs/.out.fasta.tolviewer-tmp"));
    }

    #[test]
    fn unknown_content_and_extension_is_a_clear_error() {
        let dir = std::env::temp_dir().join("tolviewer-io-unknown");
        fs::create_dir_all(&dir).unwrap();
        let p = dir.join("mystery.dat");
        fs::write(&p, b"this is not a sequence file at all\n").unwrap();
        let e = read_file(&p).unwrap_err();
        match e {
            Error::Parse { message, .. } => assert!(message.contains("cannot tell"), "{message}"),
            other => panic!("wrong error: {other}"),
        }
        let _ = fs::remove_file(&p);
    }

    #[test]
    fn write_file_is_atomic_and_leaves_no_temp_file() {
        let dir = std::env::temp_dir().join("tolviewer-io-atomic");
        fs::create_dir_all(&dir).unwrap();
        let p = dir.join("out.fasta");
        fs::write(&p, b"OLD CONTENT").unwrap();
        write_file(&aln(), &p, Format::Fasta, &WriteOptions::default()).unwrap();
        let text = fs::read_to_string(&p).unwrap();
        assert!(text.starts_with(">alpha"));
        assert!(!dir.join(".out.fasta.tolviewer-tmp").exists());
        let _ = fs::remove_file(&p);
    }
}
