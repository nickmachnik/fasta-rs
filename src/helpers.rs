use flate2::bufread::MultiGzDecoder;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

pub fn seq_id_from_description<'a>(line: &'a str, separator: &'a str, id_index: usize) -> &'a str {
    if line.contains(separator) {
        let fields = line.split(separator).collect::<Vec<&str>>();
        if id_index == 0 {
            &fields[id_index][1..]
        } else {
            fields[id_index]
        }
    } else {
        // remove `>`
        &line[1..]
    }
}

// Open file in gz or normal mode
pub fn open(path: &Path) -> Box<dyn std::io::Read> {
    if let Some(extension) = path.extension() {
        match extension.to_str().unwrap() {
            "gz" => {
                let fin = File::open(path)
                    .unwrap_or_else(|_| panic!("Could not open path: {}", path.display()));
                Box::new(MultiGzDecoder::new(BufReader::new(fin)))
            }
            _ => Box::new(BufReader::new(File::open(path).unwrap_or_else(|_| {
                panic!("Could not open path: {}", path.display())
            }))),
        }
    } else {
        Box::new(BufReader::new(File::open(path).unwrap_or_else(|_| {
            panic!("Could not open path: {}", path.display())
        })))
    }
}
