use crate::index::FastaIndex;
use crate::read::{FastaHandle, FastaReader};

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Seek, SeekFrom, Write};
use std::path::Path;

/// A HashMap representation of a Fasta file.
#[derive(Debug)]
pub struct FastaMap {
    pub id_to_seq: HashMap<String, String>,
}

impl FastaMap {
    pub fn default() -> Self {
        FastaMap {
            id_to_seq: HashMap::new(),
        }
    }

    pub fn from_fasta(path: &Path) -> Self {
        let reader = FastaReader::new(path);
        let mut entries: HashMap<String, String> = HashMap::new();
        for [header, seq] in reader {
            entries.insert(header, seq);
        }
        FastaMap { id_to_seq: entries }
    }

    pub fn from_index_with_ids(path: &Path, index: &FastaIndex, ids: &[String]) -> Self {
        let mut res = HashMap::new();
        let mut fasta_handle = FastaHandle::open_fasta(path);
        if let FastaHandle::Compressed(_) = fasta_handle {
            panic!(
                "Tried to use index on non seekable compressed file: {:?}",
                path
            );
        }

        for k in ids {
            if let Some(v) = index.id_to_offset.get(k) {
                let mut seq_buf = String::new();
                fasta_handle
                    .seek(SeekFrom::Start(*v))
                    .expect("File seek failed in `from_index_with_ids`.");

                let mut seen_header = false;
                for line in BufReader::new(&mut fasta_handle).lines() {
                    let lstring = line.unwrap();
                    if lstring.starts_with('>') {
                        if seen_header {
                            break;
                        } else {
                            seen_header = true;
                        }
                    } else if lstring == "" {
                        break;
                    } else {
                        seq_buf.push_str(&lstring);
                    }
                }
                res.insert((*k).to_string(), seq_buf);
            }
        }
        FastaMap { id_to_seq: res }
    }

    pub fn to_fasta(&self, path: &Path) {
        let mut f = match File::create(path) {
            Err(why) => panic!("couldn't create {:?}: {:?}", path, why),
            Ok(file) => BufWriter::new(file),
        };
        for (k, v) in self.id_to_seq.iter() {
            if let Err(why) = f.write_all(format!(">{}\n", k).as_bytes()) {
                panic!("couldn't write to {:?}: {:?}", path, why)
            };
            if let Err(why) = f.write_all(format!("{}\n\n", v).as_bytes()) {
                panic!("couldn't write to {:?}: {:?}", path, why)
            };
        }
    }
}
