use std::fmt;

/// Errors produced anywhere in TOLViewer.
///
/// Kept as a single flat enum so the GUI can display one error type without
/// juggling per-crate conversions.
#[derive(Debug)]
pub enum Error {
    /// Underlying I/O failure (file missing, permissions, disk full, ...).
    Io(std::io::Error),
    /// The file could not be understood as the requested format.
    Parse { format: &'static str, line: Option<usize>, message: String },
    /// The data is well-formed but cannot be written in the requested format.
    Format(String),
    /// An operation needed an alignment (equal-length rows) but got ragged data.
    NotAligned,
    /// An index (row, column, sequence) was out of range.
    OutOfRange(String),
    /// An external tool or algorithm failed.
    Algorithm(String),
    /// The user cancelled a long-running operation.
    Cancelled,
}

pub type Result<T> = std::result::Result<T, Error>;

impl Error {
    pub fn parse(format: &'static str, line: Option<usize>, message: impl Into<String>) -> Self {
        Error::Parse { format, line, message: message.into() }
    }

    pub fn format(message: impl Into<String>) -> Self {
        Error::Format(message.into())
    }

    pub fn algorithm(message: impl Into<String>) -> Self {
        Error::Algorithm(message.into())
    }

    pub fn out_of_range(message: impl Into<String>) -> Self {
        Error::OutOfRange(message.into())
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(e) => write!(f, "I/O error: {e}"),
            Error::Parse { format, line: Some(l), message } => {
                write!(f, "{format} parse error on line {l}: {message}")
            }
            Error::Parse { format, line: None, message } => {
                write!(f, "{format} parse error: {message}")
            }
            Error::Format(m) => write!(f, "cannot write file: {m}"),
            Error::NotAligned => {
                write!(f, "this operation requires an alignment (all rows the same length)")
            }
            Error::OutOfRange(m) => write!(f, "out of range: {m}"),
            Error::Algorithm(m) => write!(f, "{m}"),
            Error::Cancelled => write!(f, "cancelled"),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Error::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Error::Io(e)
    }
}
