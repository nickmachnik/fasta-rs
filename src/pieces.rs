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

/// A convenience struct for parsing the acession ids from FASTA description lines.
///
/// # Examples
/// Extract the accession ids from a UniProt formatted FASTA and write them to
/// a tsv file, one per line.
/// ```
/// use fasta::pieces::FastaAccessions;
/// use std::path::Path;
///
/// // parse the accessions
/// let accessions = FastaAccessions::from_fasta(Path::new("./resources/test.fasta"), "|", 1);
/// // write to tsv
/// accessions.to_tsv(Path::new("./resources/test.accessions")).expect("Dumping tsv failed");
/// ```
#[derive(Debug)]
pub struct FastaAccessions {
    pub accessions: Vec<String>,
}

impl FastaAccessions {
    pub fn from_fasta(path: &Path, separator: &str, id_index: usize) -> Self {
        let reader = FastaReader::new(path);
        let mut accessions = Vec::new();
        for [header, _seq] in reader {
            accessions.push(seq_id_from_description(&header, separator, id_index).to_string());
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

/// A convenient struct that wraps a sequence id to sequence length mapping.
///
/// # Examples
/// ```
/// use std::path::Path;
/// use fasta::pieces::FastaLengths;
///
/// // parse a fasta file
/// let lengths = FastaLengths::from_fasta(Path::new("./resources/test.fasta"), "|", 1);
/// // write to json
/// lengths.to_json(Path::new("./resources/test.accessions")).expect("JSON dump failed");
/// ```
#[derive(Debug, PartialEq)]
pub struct FastaLengths {
    pub sequence_lengths: HashMap<String, usize>,
}

impl FastaLengths {
    pub fn from_fasta(path: &Path, separator: &str, id_index: usize) -> Self {
        let reader = FastaReader::new(path);
        let mut entries: HashMap<String, usize> = HashMap::new();
        for [header, seq] in reader {
            entries.insert(
                seq_id_from_description(&header, separator, id_index).to_string(),
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
    fn accessions_from_fasta_short() {
        assert_eq!(
            FastaAccessions::from_fasta(Path::new("./resources/test_short_descr.fasta"), "|", 1)
                .accessions,
            vec!["Q2HZH0", "P93158", "H0VS30"]
        )
    }

    #[test]
    fn accessions_from_fasta_long() {
        assert_eq!(
            FastaAccessions::from_fasta(Path::new("./resources/test.fasta"), "|", 1).accessions,
            vec!["Q2HZH0", "P93158", "H0VS30"]
        )
    }

    #[test]
    fn get_single_fasta_entry() {
        let index =
            crate::index::FastaIndex::from_json(Path::new("./resources/test.index")).unwrap();
        let entry = FastaEntry::from_index(
            Path::new("./resources/test.fasta"),
            *index.id_to_offset.get("P93158").unwrap(),
        )
        .unwrap();
        let expected = FastaEntry {
            description: "tr|P93158|P93158_GOSHI Annexin (Fragment) OS=Gossypium \
            hirsutum OX=3635 GN=AnnGh2 PE=2 SV=1"
                .to_string(),
            sequence: "TLKVPVHVPSPSEDAEWQLRKAFEGWGTNEQLIIDILAHRNAAQRNSIRKVYGEAYGEDL\
            LKCLEKELTSDFERAVLLFTLDPAERDAHLANEATKKFTSSNWILMEIACSRSSHELLNV"
                .to_string(),
        };
        assert_eq!(entry, expected);
    }

    #[test]
    fn lengths_from_fasta() {
        let lengths = FastaLengths::from_fasta(Path::new("./resources/test.fasta"), "|", 1);
        let mut exp_map = HashMap::new();
        exp_map.insert("H0VS30".to_string(), 180);
        exp_map.insert("Q2HZH0".to_string(), 120);
        exp_map.insert("P93158".to_string(), 120);
        assert_eq!(
            lengths,
            FastaLengths {
                sequence_lengths: exp_map
            }
        );
        println!("{:?}", lengths);
    }
}
