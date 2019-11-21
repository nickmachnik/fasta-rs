extern crate flate2;

use flate2::read::GzDecoder;
use std::fs::File;
use std::path::Path;
use std::io;
use std::io::BufRead;
use std::collections::HashMap;

// Open file in gz or normal mode
pub fn open(path: &Path) -> Box<dyn std::io::Read> {
    match path.extension().expect("No path extension found!").to_str().unwrap() {
        "gz" => {
            let fin = File::open(path)
                .unwrap_or_else(|_| panic!("Could not open path: {}", path.display()));
            Box::new(GzDecoder::new(fin))
        },
        _ => {
            Box::new(File::open(path)
                .unwrap_or_else(|_| panic!("Could not open path: {}", path.display())))
        }
    }
}


pub struct FastaReader {
    lines: std::io::Lines<std::io::BufReader<std::boxed::Box<dyn std::io::Read>>>,
    header: Option<String>,
    seq_buf: String,
}

impl FastaReader {
    pub fn new(path: &Path) -> FastaReader {
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
    pub fn new(path: &Path) -> FastaSeqs {
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
pub struct FastaMap {
    pub entries: HashMap<String, String>,
}

impl FastaMap {
    pub fn new(path: &Path) -> FastaMap {
        let reader = FastaReader::new(path);
        let mut entries: HashMap<String, String> = HashMap::new();
        for [header, seq] in reader {
            entries.insert(header, seq);
        }
        FastaMap {
            entries
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
