//! Structs that represent specific pieces of information
//! found in FASTA files. Useful for extracting and storing
//! these parts.

use crate::errors;
use crate::helpers::seq_id_from_description;
use crate::read::FastaReader;

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::error;
use std::fs::File;
use std::io;
use std::io::prelude::Seek;
use std::io::BufWriter;
use std::io::{BufRead, BufReader, SeekFrom, Write};
use std::path::Path;

/// A convenience struct to extract and write
/// the accessions from fasta deflines.
#[derive(Debug)]
pub struct FastaAccessions {
    pub accessions: Vec<String>,
}

impl FastaAccessions {
    pub fn from_fasta(path: &Path) -> Self {
        let reader = FastaReader::new(path);
        let mut accessions = Vec::new();
        for [header, _seq] in reader {
            accessions.push(seq_id_from_description(&header, "|", 1).to_string());
        }
        FastaAccessions { accessions }
    }

    /// Writes the accessions to json.
    pub fn to_json(&self, outpath: &Path) -> Result<(), io::Error> {
        let mut file = BufWriter::new(File::create(&outpath)?);
        serde_json::to_writer(&mut file, &self.accessions)?;
        Ok(())
    }

    /// Writes the accessions to a txt file, one per line.
    pub fn to_tsv(&self, outpath: &Path) -> Result<(), io::Error> {
        let mut file = BufWriter::new(File::create(&outpath)?);
        for id in &self.accessions {
            file.write_all(format!("{}\n", id).as_bytes())?;
        }
        Ok(())
    }
}

/// A HashMap mapping sequence ids to sequence lengths.
#[derive(Debug)]
pub struct FastaLengths {
    pub sequence_lengths: HashMap<String, usize>,
}

impl FastaLengths {
    pub fn default() -> Self {
        FastaLengths {
            sequence_lengths: HashMap::new(),
        }
    }

    pub fn from_fasta(path: &Path) -> Self {
        let reader = FastaReader::new(path);
        let mut entries: HashMap<String, usize> = HashMap::new();
        for [header, seq] in reader {
            entries.insert(
                seq_id_from_description(&header, "|", 0).to_string(),
                seq.len(),
            );
        }
        FastaLengths {
            sequence_lengths: entries,
        }
    }

    /// Writes the ID -> Sequence length mapping to .json.
    pub fn to_json(&self, outpath: &Path) -> Result<(), io::Error> {
        let mut file = BufWriter::new(File::create(&outpath)?);
        serde_json::to_writer(&mut file, &self.sequence_lengths)?;
        Ok(())
    }
}

/// A single .fasta entry with description and sequence fields.
#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct FastaEntry {
    pub description: String,
    pub sequence: String,
}

impl FastaEntry {
    pub fn from_index(data: &Path, index: u64) -> Result<Self, Box<dyn error::Error>> {
        let mut handle = BufReader::new(File::open(data)?);
        handle.seek(SeekFrom::Start(index))?;

        let mut lines = handle.lines();
        let line = lines.next().unwrap().unwrap();
        let description = if line.starts_with('>') {
            line[1..].to_string()
        } else {
            return Err(Box::new(errors::ParseError::new(
                errors::ErrorKind::IndexNotAtDescription,
                "No description line found at index when reading entry!",
            )));
        };

        let mut entry = FastaEntry {
            description,
            sequence: String::new(),
        };

        for l in lines {
            let line = l.unwrap();
            if line == "" || line.starts_with('>') {
                break;
            } else {
                entry.sequence.push_str(&line);
            }
        }

        Ok(entry)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accessions_from_fasta_short() {}

    #[test]
    fn accessions_from_fasta_long() {}
}
