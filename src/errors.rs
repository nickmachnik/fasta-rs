//! Errors for this crate.

use std::error;
use std::fmt;

/// The error type for everything that can go wrong in fasta parsing.
#[derive(Debug)]
pub struct ParseError<'a> {
    kind: ErrorKind,
    message: &'a str,
}

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum ErrorKind {
    /// Index points to a line that is not a description line.
    IndexNotAtDescription,
}

impl ErrorKind {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            ErrorKind::IndexNotAtDescription => "Index points to a non-description line.",
        }
    }
}

impl<'a> ParseError<'a> {
    pub fn new(kind: ErrorKind, message: &'a str) -> ParseError<'a> {
        Self { kind, message }
    }
}

impl<'a> fmt::Display for ParseError<'a> {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(fmt, "{}", self.message)
    }
}

impl<'a> error::Error for ParseError<'a> {
    fn description(&self) -> &str {
        self.kind.as_str()
    }
}
