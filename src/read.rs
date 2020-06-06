//! Utilities for reading whole FASTA files into iterators.

use crate::helpers::open;
use flate2::bufread::MultiGzDecoder;
use std::fs::File;
use std::io::prelude::Seek;
use std::io::{BufRead, BufReader, Read, SeekFrom};
use std::path::Path;

/// An enum that wraps compressed (gz) and uncompressed files.
#[derive(Debug)]
pub enum FastaHandle {
    Compressed(MultiGzDecoder<BufReader<File>>),
    Uncompressed(BufReader<File>),
}

impl Read for FastaHandle {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            FastaHandle::Compressed(s) => s.read(buf),
            FastaHandle::Uncompressed(s) => s.read(buf),
        }
    }
}

impl Seek for FastaHandle {
    fn seek(&mut self, pos: SeekFrom) -> std::io::Result<u64> {
        match self {
            FastaHandle::Compressed(_s) => panic!("Cannot seek in gzipped file!"),
            FastaHandle::Uncompressed(s) => s.seek(pos),
        }
    }
}

impl FastaHandle {
    pub fn open_fasta(path: &Path) -> FastaHandle {
        if let Some(extension) = path.extension() {
            match extension.to_str().unwrap() {
                "gz" => {
                    let fin = File::open(path)
                        .unwrap_or_else(|_| panic!("Could not open path: {}", path.display()));
                    FastaHandle::Compressed(MultiGzDecoder::new(BufReader::new(fin)))
                }
                _ => FastaHandle::Uncompressed(BufReader::new(
                    File::open(path)
                        .unwrap_or_else(|_| panic!("Could not open path: {}", path.display())),
                )),
            }
        } else {
            FastaHandle::Uncompressed(BufReader::new(
                File::open(path)
                    .unwrap_or_else(|_| panic!("Could not open path: {}", path.display())),
            ))
        }
    }
}

/// A reader that visits entries in a FASTA file one by one.
///
/// # Examples
///
/// Iterate through a FASTA file:
/// ```
/// use fasta::read::FastaReader;
/// use std::path::Path;
///
/// let infile = Path::new("./resources/test.fasta");
/// for [description, seq] in FastaReader::new(infile) {
///     println!("{:?}", description);
///     println!("{:?}", seq);
/// }
/// ```
pub struct FastaReader {
    lines: std::io::Lines<std::io::BufReader<std::boxed::Box<dyn std::io::Read>>>,
    description: Option<String>,
    seq_buf: String,
}

impl FastaReader {
    pub fn new(path: &Path) -> Self {
        let reader = open(&path);
        let mut res = FastaReader {
            lines: BufReader::new(reader).lines(),
            description: None,
            seq_buf: String::new(),
        };

        // find first description
        while res.description == None {
            match res.lines.next() {
                Some(s) => {
                    let line = s.unwrap();
                    if line.starts_with('>') {
                        res.description = Some(line.to_string());
                    }
                }
                None => panic!("Reached EOF in FASTA parsing; No description in file?"),
            }
        }
        res
    }
}

impl Iterator for FastaReader {
    type Item = [String; 2];

    fn next(&mut self) -> Option<Self::Item> {
        self.seq_buf.clear();

        while let Some(l) = self.lines.next() {
            let line = l.unwrap();
            if line.starts_with('>') {
                let res = [self.description.clone().unwrap(), self.seq_buf.clone()];
                self.description = Some(line);
                return Some(res);
            } else {
                self.seq_buf.push_str(&line);
            }
        }

        match self.seq_buf.len() {
            0 => None,
            _ => Some([self.description.clone().unwrap(), self.seq_buf.clone()]),
        }
    }
}
