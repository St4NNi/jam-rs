use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

const ROW_SIZE: usize = 24;
const MERGE_FAN_IN: usize = 32;

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct RunRecord {
    pub(crate) packed_key: u64,
    pub(crate) contig_id: u32,
    pub(crate) position: u64,
    pub(crate) canonical_orientation: bool,
}

impl RunRecord {
    fn encode(self) -> [u8; ROW_SIZE] {
        let mut bytes = [0; ROW_SIZE];
        bytes[..8].copy_from_slice(&self.packed_key.to_le_bytes());
        bytes[8..12].copy_from_slice(&self.contig_id.to_le_bytes());
        bytes[12] = u8::from(self.canonical_orientation);
        bytes[16..24].copy_from_slice(&self.position.to_le_bytes());
        bytes
    }

    fn decode(bytes: [u8; ROW_SIZE]) -> io::Result<Self> {
        if bytes[12] > 1 || bytes[13..16].iter().any(|byte| *byte != 0) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid JIDX run row",
            ));
        }
        Ok(Self {
            packed_key: u64::from_le_bytes(bytes[..8].try_into().expect("run key")),
            contig_id: u32::from_le_bytes(bytes[8..12].try_into().expect("run contig")),
            canonical_orientation: bytes[12] == 1,
            position: u64::from_le_bytes(bytes[16..24].try_into().expect("run position")),
        })
    }
}

pub(crate) struct Runs {
    directory: PathBuf,
    max_records: usize,
    records: Vec<RunRecord>,
    paths: Vec<PathBuf>,
    next_file: u64,
}

impl Runs {
    pub(crate) fn new(directory: &Path, max_records: usize) -> io::Result<Self> {
        if max_records == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "JIDX run capacity must be nonzero",
            ));
        }
        if !fs::metadata(directory)?.is_dir() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "JIDX run path is not a directory",
            ));
        }
        Ok(Self {
            directory: directory.to_path_buf(),
            max_records,
            records: Vec::with_capacity(max_records),
            paths: Vec::new(),
            next_file: 0,
        })
    }

    pub(crate) fn push(&mut self, record: RunRecord) -> io::Result<()> {
        self.records.push(record);
        if self.records.len() == self.max_records {
            self.spill()?;
        }
        Ok(())
    }

    pub(crate) fn finish(mut self) -> io::Result<Merge> {
        self.spill()?;
        self.records = Vec::new();
        let original_paths = self.paths.clone();
        let mut intermediate_paths = Vec::new();
        let mut current_paths = std::mem::take(&mut self.paths);

        while current_paths.len() > MERGE_FAN_IN {
            let mut next_paths = Vec::with_capacity(current_paths.len().div_ceil(MERGE_FAN_IN));
            for chunk in current_paths.chunks(MERGE_FAN_IN) {
                let (path, file) = self.create_file("merge")?;
                merge_to_file(chunk, file)?;
                intermediate_paths.push(path.clone());
                next_paths.push(path);
            }
            current_paths = next_paths;
        }

        let merge = Merge::open(current_paths.clone(), true)?;
        for path in original_paths.iter().chain(&intermediate_paths) {
            if !current_paths.contains(path) {
                let _ = fs::remove_file(path);
            }
        }
        Ok(merge)
    }

    fn spill(&mut self) -> io::Result<()> {
        if self.records.is_empty() {
            return Ok(());
        }
        self.records.sort_unstable();
        let (path, file) = self.create_file("run")?;
        let mut writer = BufWriter::new(file);
        for record in &self.records {
            writer.write_all(&record.encode())?;
        }
        writer.flush()?;
        self.paths.push(path);
        self.records.clear();
        Ok(())
    }

    fn create_file(&mut self, kind: &str) -> io::Result<(PathBuf, File)> {
        let ordinal = self.next_file;
        self.next_file = self
            .next_file
            .checked_add(1)
            .ok_or_else(|| io::Error::other("too many JIDX run files"))?;
        let path = self.directory.join(format!(".jidx-{kind}-{ordinal}.bin"));
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)?;
        Ok((path, file))
    }
}

struct RunReader {
    reader: BufReader<File>,
    previous: Option<RunRecord>,
}

impl RunReader {
    fn open(path: &Path) -> io::Result<Self> {
        Ok(Self {
            reader: BufReader::new(File::open(path)?),
            previous: None,
        })
    }

    fn next_record(&mut self) -> io::Result<Option<RunRecord>> {
        let mut bytes = [0; ROW_SIZE];
        let read = self.reader.read(&mut bytes)?;
        if read == 0 {
            return Ok(None);
        }
        self.reader.read_exact(&mut bytes[read..])?;
        let record = RunRecord::decode(bytes)?;
        if self.previous.is_some_and(|previous| previous > record) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "unsorted JIDX run",
            ));
        }
        self.previous = Some(record);
        Ok(Some(record))
    }
}

pub(crate) struct Merge {
    readers: Vec<RunReader>,
    heap: BinaryHeap<Reverse<(RunRecord, usize)>>,
    paths: Vec<PathBuf>,
    cleanup_on_eof: bool,
    complete: bool,
}

impl Merge {
    fn open(paths: Vec<PathBuf>, cleanup_on_eof: bool) -> io::Result<Self> {
        if paths.len() > MERGE_FAN_IN {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "too many JIDX merge inputs",
            ));
        }
        let mut readers = Vec::with_capacity(paths.len());
        let mut heap = BinaryHeap::new();
        for (index, path) in paths.iter().enumerate() {
            let mut reader = RunReader::open(path)?;
            if let Some(record) = reader.next_record()? {
                heap.push(Reverse((record, index)));
            }
            readers.push(reader);
        }
        Ok(Self {
            readers,
            heap,
            paths,
            cleanup_on_eof,
            complete: false,
        })
    }

    pub(crate) fn next_record(&mut self) -> io::Result<Option<RunRecord>> {
        let Some(Reverse((record, reader_index))) = self.heap.pop() else {
            self.finish_cleanup();
            return Ok(None);
        };
        if let Some(next) = self.readers[reader_index].next_record()? {
            self.heap.push(Reverse((next, reader_index)));
        }
        Ok(Some(record))
    }

    fn finish_cleanup(&mut self) {
        if self.complete {
            return;
        }
        self.complete = true;
        if self.cleanup_on_eof {
            for path in &self.paths {
                let _ = fs::remove_file(path);
            }
        }
    }
}

fn merge_to_file(paths: &[PathBuf], file: File) -> io::Result<()> {
    let mut merge = Merge::open(paths.to_vec(), false)?;
    let mut writer = BufWriter::new(file);
    while let Some(record) = merge.next_record()? {
        writer.write_all(&record.encode())?;
    }
    writer.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn collect(mut merge: Merge) -> io::Result<Vec<RunRecord>> {
        let mut records = Vec::new();
        while let Some(record) = merge.next_record()? {
            records.push(record);
        }
        Ok(records)
    }

    #[test]
    fn forced_spills_preserve_global_tuple_order_and_duplicates() {
        let directory = tempfile::tempdir().unwrap();
        let duplicate = RunRecord {
            packed_key: 7,
            contig_id: 2,
            position: 11,
            canonical_orientation: false,
        };
        let mut expected = vec![
            duplicate,
            RunRecord {
                packed_key: 3,
                contig_id: 5,
                position: 2,
                canonical_orientation: true,
            },
            duplicate,
            RunRecord {
                packed_key: 7,
                contig_id: 2,
                position: 11,
                canonical_orientation: true,
            },
            RunRecord {
                packed_key: 3,
                contig_id: 1,
                position: 9,
                canonical_orientation: false,
            },
        ];
        let mut runs = Runs::new(directory.path(), 2).unwrap();
        for record in &expected {
            runs.push(*record).unwrap();
        }
        expected.sort_unstable();
        assert_eq!(collect(runs.finish().unwrap()).unwrap(), expected);
        assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 0);
    }

    #[test]
    fn truncated_row_is_rejected_and_retained() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("truncated.bin");
        fs::write(&path, [0; ROW_SIZE - 1]).unwrap();
        let error = match Merge::open(vec![path.clone()], true) {
            Ok(_) => panic!("truncated run unexpectedly opened"),
            Err(error) => error,
        };
        assert_eq!(error.kind(), io::ErrorKind::UnexpectedEof);
        assert!(path.exists());
    }

    #[test]
    fn more_than_thirty_two_spills_are_reduced_before_final_merge() {
        let directory = tempfile::tempdir().unwrap();
        let mut runs = Runs::new(directory.path(), 1).unwrap();
        let mut expected = Vec::new();
        for value in (0..65u64).rev() {
            let record = RunRecord {
                packed_key: value % 9,
                contig_id: value as u32,
                position: value * 3,
                canonical_orientation: value % 2 == 0,
            };
            expected.push(record);
            runs.push(record).unwrap();
        }
        expected.sort_unstable();
        let merge = runs.finish().unwrap();
        assert!(merge.readers.len() <= MERGE_FAN_IN);
        assert_eq!(collect(merge).unwrap(), expected);
        assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 0);
    }
}
