extern crate bgzip;
extern crate hashbrown;
extern crate serde;

use std::io::{SeekFrom, Read, Seek};
use bgzip::read::BGzReader;
use serde::{Serialize, Deserialize};
use std::fs::{File, read_to_string};
use std::path::Path;
use std::io;
use std::io::BufRead;
use hashbrown::HashMap;
use std::error;
use std::fmt;

#[derive(Debug)]
pub enum FastaHandle {
    Compressed(BGzReader<File>),
    Uncompressed(File),
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
            FastaHandle::Compressed(s) => s.seek(pos),
            FastaHandle::Uncompressed(s) => s.seek(pos),
        }
    }
}

pub fn open_fasta(path: &Path) -> FastaHandle {
    if let Some(extension) = path.extension() {
        match extension.to_str().unwrap() {
            "gz" => {
                let fin = File::open(path)
                    .unwrap_or_else(|_| panic!("Could not open path: {}", path.display()));
                FastaHandle::Compressed(BGzReader::new(fin).unwrap())
            },
            _ => {
                FastaHandle::Uncompressed(File::open(path)
                    .unwrap_or_else(|_| panic!("Could not open path: {}", path.display())))
            }
        }
    } else {
        FastaHandle::Uncompressed(File::open(path)
            .unwrap_or_else(|_| panic!("Could not open path: {}", path.display())))
    }
}

// Open file in gz or normal mode
pub fn open(path: &Path) -> Box<dyn std::io::Read> {
    if let Some(extension) = path.extension() {
        match extension.to_str().unwrap() {
            "gz" => {
                let fin = File::open(path)
                    .unwrap_or_else(|_| panic!("Could not open path: {}", path.display()));
                Box::new(BGzReader::new(fin).unwrap())
            },
            _ => {
                Box::new(File::open(path)
                    .unwrap_or_else(|_| panic!("Could not open path: {}", path.display())))
            }
        }
    } else {
        Box::new(File::open(path)
            .unwrap_or_else(|_| panic!("Could not open path: {}", path.display())))
    }
}

pub struct FastaReader {
    lines: std::io::Lines<std::io::BufReader<std::boxed::Box<dyn std::io::Read>>>,
    header: Option<String>,
    seq_buf: String,
}

impl FastaReader {
    pub fn new(path: &Path) -> Self {
        let reader = open(&path);
        let mut res = FastaReader {
            lines: io::BufReader::new(reader).lines(),
            header: None,
            seq_buf: String::new(),
        };
        
        // find first header
        while res.header == None {
            match res.lines.next() {
                Some(s) => {
                    let line = s.unwrap();
                    if line.starts_with('>') {
                        res.header = Some(line[1..].to_string());
                    }
                },
                None => panic!("Reached EOF in FASTA parsing; No header in file?"),
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
                let res = [self.header.clone().unwrap(), self.seq_buf.clone()];
                self.header = Some(line[1..].to_string());
                return Some(res)
            } else {
                self.seq_buf.push_str(&line);
            }
        }

        match self.seq_buf.len() {
            0 => None,
            _ => Some([self.header.clone().unwrap(), self.seq_buf.clone()]),
        }
    }
}

// separate headers and seqs Vecs with indices defining the mapping
pub struct FastaSeqs {
    pub headers: Vec<String>,
    pub seqs: Vec<String>,
}

impl FastaSeqs {
    pub fn new(path: &Path) -> Self {
        let reader = FastaReader::new(path);
        let mut headers: Vec<String> = Vec::new();
        let mut seqs: Vec<String> = Vec::new();
        for [header, seq] in reader {
            headers.push(header);
            seqs.push(seq);
        }
        FastaSeqs {
            headers,
            seqs
        }
    }
}


// header -> sequence mapping
#[derive(Debug)]
pub struct FastaMap {
    pub id_to_seq: HashMap<String, String>
}

impl FastaMap {
    pub fn default() -> Self {
        FastaMap { id_to_seq: HashMap::new() }
    }

    pub fn from_fasta(path: &Path) -> Self {
        let reader = FastaReader::new(path);
        let mut entries: HashMap<String, String> = HashMap::new();
        for [header, seq] in reader {
            entries.insert(header, seq);
        }
        FastaMap { id_to_seq: entries }
    }

    pub fn from_index_with_ids (
        path: &Path,
        index: &FastaIndex,
        ids: &[&str]
    ) -> Result<Self, MissingID>
    {
        let mut res = HashMap::new();
        let mut fasta_handle = open_fasta(path);
        for k in ids {
            if let Some(v) = index.id_to_offset.get(*k) {
                let mut seq_buf = String::new();
                fasta_handle.seek(io::SeekFrom::Start(*v))
                    .expect("File seek failed in `from_index_with_ids`.");

                for line in io::BufReader::new(&mut fasta_handle).lines() {
                    let lstring = line.unwrap();
                    if lstring.starts_with('>') {
                        continue
                    } else if lstring == "\n" {
                        break
                    } else {
                        seq_buf.push_str(&lstring);
                    }
                }
                res.insert((*k).to_string(), seq_buf);
            } else {
                return Err(MissingID { id: (*k).to_string() })
            }
        }
        Ok(FastaMap{ id_to_seq: res })
    }
}

#[derive(Debug)]
pub struct MissingID { pub id: String }

impl fmt::Display for MissingID {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ID {} not found in index.", self.id)
    }
}

impl error::Error for MissingID {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        // Generic error, underlying cause isn't tracked.
        None
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FastaIndex {
    id_to_offset: HashMap<String, u64>
}

impl FastaIndex {
    pub fn new(path: &Path) -> Self {
        let mut res = HashMap::new();

        let fin = open(&path);
        let mut reader = io::BufReader::new(fin);
        let mut line_buf = String::new();
        let mut global_offset: u64 = 0;

        let mut len = reader.read_line(&mut line_buf)
                            .expect("Failed to read line!");
        while len != 0 {
            if line_buf.starts_with('>') {
                let header_split = line_buf.split('|').collect::<Vec<&str>>();
                let key = if header_split.len() > 1 {
                    header_split[1]
                } else {
                    header_split[0]
                };
                if let Some(_old_entry) = res.insert(key.to_string(), global_offset) {
                    panic!("Multiple entries found for id: {:?}", key);
                };
            }

            global_offset += len as u64;
            line_buf.clear();
            len = reader.read_line(&mut line_buf)
                        .expect("Failed to read line!");
        }

        FastaIndex{ id_to_offset: res }
    }

    pub fn from_json(path: &Path) -> Self {
        let json_file_str = read_to_string(path).expect("file not found");
        serde_json::from_str(&json_file_str).expect("error while reading json")
    }

    pub fn to_json(&self, outpath: &Path) {
        let mut file = match File::create(&outpath) {
            Err(why) => panic!("couldn't create {:?}: {:?}", outpath, why),
            Ok(file) => file,
        };
        if let Err(why) = serde_json::to_writer(&mut file, self) {
            panic!("couldn't write to {:?}: {:?}", outpath, why)
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn index_building() {
        let mut expected = HashMap::new();
        expected.insert("Q2HZH0".to_string(), 0_u64);
        expected.insert("G7PNW8".to_string(), 2764_u64);
        expected.insert("P9WNK5".to_string(), 1403_u64);
        expected.insert("H0VS30".to_string(), 774_u64);
        expected.insert("G1KTG2".to_string(), 1638_u64);
        expected.insert("Q8I5U1".to_string(), 2189_u64);
        expected.insert("P93158".to_string(), 359_u64);
        
        assert_eq!(FastaIndex::new(Path::new("./resources/test.fasta")).id_to_offset, expected);
    }

    #[test]
    fn indexed_reading() {
        let index = FastaIndex::new(Path::new("./resources/test.fasta"));
        let fasta_map = FastaMap::from_index_with_ids(
            Path::new("./resources/test.fasta"),
            &index,
            &["P9WNK5", "Q8I5U1"]
        ).unwrap();
        assert_eq!(fasta_map.id_to_seq.len(), 2);
        assert!(fasta_map.id_to_seq.contains_key("Q8I5U1"));
        assert!(fasta_map.id_to_seq.contains_key("P9WNK5"));
    }

    #[test]
    fn index_dump_and_load() {
        let index = FastaIndex::new(Path::new("./resources/test.fasta"));
        index.to_json(Path::new("./resources/test.index"));
        let loaded = FastaIndex::from_json(Path::new("./resources/test.index"));
        assert_eq!(index.id_to_offset, loaded.id_to_offset);

    }    
}
