//! An index that stores byte offsets of individual entries
//! in FASTA files.

use crate::helpers::seq_id_from_description;
use crate::read::FastaHandle;

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::{read_to_string, File};
use std::io::{BufRead, BufReader, BufWriter, Error};
use std::path::Path;

/// An index into FASTA files.
///
/// Wraps a sequence id -> byte offset mapping. The sequence accessions
/// provided in the description lines are used as sequence ids.
/// The description line format has to be specified when creating an index
/// by indicating the field separator and the index (0-based) of the field
/// which stores the accession. For example, the standard UniProt format
/// uses "|" as a separator and stores the accession in the first field.
///
/// # Examples
///
/// Create and index from a FASTA file and write to json and load:
/// ```
/// use fasta::index::FastaIndex;
///
/// // create the index
/// let index = FastaIndex::new(Path::new("foo.fasta"), "|", 1);
/// // write to file
/// index.to_json(Path::new("foo.index")).expect("Failed to dump json.");
/// // load from json
/// assert_eq!(index, FastaIndex::from_json(Path::new("foo.index")).unwrap());
/// ```
#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct FastaIndex {
    pub id_to_offset: HashMap<String, u64>,
}

impl FastaIndex {
    pub fn new(path: &Path, separator: &str, id_index: usize) -> Self {
        let mut res = HashMap::new();

        let fasta_handle = FastaHandle::open_fasta(path);
        if let FastaHandle::Compressed(_) = fasta_handle {
            panic!(
                "Tried to build index on non seekable compressed file: {:?}",
                path
            );
        }
        let mut reader = BufReader::new(fasta_handle);
        let mut line_buf = String::new();
        let mut global_offset: u64 = 0;

        let mut len = reader
            .read_line(&mut line_buf)
            .expect("Failed to read line!");
        while len != 0 {
            if line_buf.starts_with('>') {
                line_buf.pop();
                let key = seq_id_from_description(&line_buf, separator, id_index);
                if let Some(_old_entry) = res.insert(key.to_string(), global_offset) {
                    panic!("Multiple entries found for id: {:?}", key);
                };
            }

            global_offset += len as u64;
            line_buf.clear();
            len = reader
                .read_line(&mut line_buf)
                .expect("Failed to read line!");
        }

        FastaIndex { id_to_offset: res }
    }

    pub fn from_json(path: &Path) -> Result<Self, Error> {
        let json_file_str = read_to_string(path)?;
        let res = serde_json::from_str(&json_file_str)?;
        Ok(res)
    }

    pub fn to_json(&self, outpath: &Path) -> Result<(), Error> {
        let mut file = BufWriter::new(File::create(&outpath)?);
        serde_json::to_writer(&mut file, self)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::map::FastaMap;
    use crate::pieces::FastaEntry;

    #[test]
    fn index_building() {
        assert_eq!(
            FastaIndex::new(Path::new("./resources/test.fasta"), "|", 1),
            FastaIndex::from_json(Path::new("./resources/test.index")).unwrap()
        );
    }

    #[test]
    fn indexed_reading() {
        let index = FastaIndex::new(Path::new("./resources/test.fasta"), "|", 1);
        let fasta_map = FastaMap::from_index_with_ids(
            Path::new("./resources/test.fasta"),
            &index,
            &["P93158".to_string(), "Q2HZH0".to_string()],
        );
        assert_eq!(fasta_map.id_to_seq.len(), 2);
        assert!(fasta_map.id_to_seq.contains_key("P93158"));
        assert!(fasta_map.id_to_seq.contains_key("Q2HZH0"));
    }

    #[test]
    fn individual_entry_from_index() {
        let index = FastaIndex::new(Path::new("./resources/test.fasta"), "|", 1);
        let entry = FastaEntry::from_index(
            Path::new("./resources/test.fasta"),
            *index.id_to_offset.get("P93158").unwrap(),
        )
        .unwrap();

        let exp_entry = FastaEntry {
            description: "tr|P93158|P93158_GOSHI Annexin (Fragment) OS=Gossypium hirsutum OX=3635 GN=AnnGh2 PE=2 SV=1".to_string(),
            sequence: "TLKVPVHVPSPSEDAEWQLRKAFEGWGTNEQLIIDILAHRNAAQRNSIRKVYGEAYGEDL\
            LKCLEKELTSDFERAVLLFTLDPAERDAHLANEATKKFTSSNWILMEIACSRSSHELLNV".to_string()
        };
        assert_eq!(exp_entry, entry);
    }
}
