use anyhow::Result;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use crate::format::{ENTRY_SIZE, Entry};

pub fn expand_input_paths(input_paths: &[PathBuf]) -> Result<Vec<PathBuf>> {
    let mut expanded_paths = Vec::new();
    let mut rejected_paths = Vec::new();

    for path in input_paths {
        if path.is_dir() {
            let before = expanded_paths.len();
            for entry in std::fs::read_dir(path)? {
                let entry = entry?;
                let file_path = entry.path();

                if is_sequence_file(&file_path) {
                    expanded_paths.push(file_path);
                }
            }
            if expanded_paths.len() == before {
                rejected_paths.push(format!(
                    "{} (directory contains no sequence files)",
                    path.display()
                ));
            }
        } else if path.is_file() {
            if is_sequence_file(path) {
                expanded_paths.push(path.clone());
            } else {
                let content = std::fs::read_to_string(path)?;
                let base_dir = path.parent().unwrap_or_else(|| Path::new("."));
                for line in content.lines() {
                    let line = line.trim();
                    if line.is_empty() || line.starts_with('#') {
                        continue;
                    }
                    let listed_path = PathBuf::from(line);
                    let file_path = if listed_path.is_relative() {
                        base_dir.join(listed_path)
                    } else {
                        listed_path
                    };
                    if file_path.exists() && is_sequence_file(&file_path) {
                        expanded_paths.push(file_path.canonicalize().unwrap_or(file_path));
                    } else {
                        rejected_paths.push(format!(
                            "{} (listed in {}, missing or unsupported extension)",
                            file_path.display(),
                            path.display()
                        ));
                    }
                }
            }
        } else {
            rejected_paths.push(format!("{} (does not exist)", path.display()));
        }
    }

    if expanded_paths.is_empty() {
        let mut message = "No valid sequence files found in input paths".to_string();
        if !rejected_paths.is_empty() {
            message.push_str(": ");
            message.push_str(&rejected_paths.join(", "));
        }
        return Err(anyhow::anyhow!(message));
    }

    expanded_paths.sort();
    Ok(expanded_paths)
}

pub fn is_sequence_file(path: &Path) -> bool {
    if let Some(ext) = path.extension().map(|e| e.to_string_lossy().to_lowercase()) {
        if matches!(ext.as_str(), "gz" | "bz2" | "xz" | "zst" | "zstd")
            && let Some(stem_ext) = path.file_stem().and_then(|s| Path::new(s).extension())
        {
            let stem_ext = stem_ext.to_string_lossy().to_lowercase();
            return matches!(
                stem_ext.as_str(),
                "fasta" | "fa" | "fas" | "fna" | "fastq" | "fq"
            );
        }
        return matches!(
            ext.as_str(),
            "fasta" | "fa" | "fas" | "fna" | "fastq" | "fq"
        );
    }
    false
}

pub fn read_entries<P: AsRef<Path>>(path: P) -> io::Result<Vec<Entry>> {
    let file = File::open(path)?;
    let file_size = file.metadata()?.len() as usize;
    if !file_size.is_multiple_of(ENTRY_SIZE) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("entry file size {file_size} is not a multiple of {ENTRY_SIZE}"),
        ));
    }
    let entry_count = file_size / ENTRY_SIZE;

    let mut reader = BufReader::with_capacity(8 * 1024 * 1024, file);
    let mut entries = Vec::with_capacity(entry_count);

    let mut buf = [0u8; ENTRY_SIZE];
    while reader.read_exact(&mut buf).is_ok() {
        let hash = u64::from_le_bytes(buf[0..8].try_into().unwrap());
        let sample_id = u32::from_le_bytes(buf[8..12].try_into().unwrap());
        entries.push(Entry::new(hash, sample_id));
    }

    Ok(entries)
}

pub fn write_entries<P: AsRef<Path>>(path: P, entries: &[Entry]) -> io::Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::with_capacity(8 * 1024 * 1024, file);
    writer.write_all(bytemuck::cast_slice(entries))?;
    writer.flush()
}

pub struct EntryWriter {
    writer: BufWriter<File>,
    count: u64,
}

impl EntryWriter {
    pub fn new<P: AsRef<Path>>(path: P, buffer_size: usize) -> io::Result<Self> {
        let file = File::create(path)?;
        Ok(Self {
            writer: BufWriter::with_capacity(buffer_size, file),
            count: 0,
        })
    }

    pub fn write(&mut self, entry: &Entry) -> io::Result<()> {
        self.writer.write_all(bytemuck::bytes_of(entry))?;
        self.count += 1;
        Ok(())
    }

    pub fn write_batch(&mut self, entries: &[Entry]) -> io::Result<()> {
        self.writer.write_all(bytemuck::cast_slice(entries))?;
        self.count += entries.len() as u64;
        Ok(())
    }

    pub fn count(&self) -> u64 {
        self.count
    }

    pub fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}

pub fn extract_unique_hashes(entries: &[Entry]) -> Vec<u64> {
    let mut unique = Vec::with_capacity(entries.len() / 10);
    let mut last_hash: Option<u64> = None;

    for entry in entries {
        let h = entry.hash;
        if last_hash != Some(h) {
            unique.push(h);
            last_hash = Some(h);
        }
    }

    unique
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_entry_roundtrip() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("test.bin");

        let entries = vec![Entry::new(100, 1), Entry::new(200, 2), Entry::new(300, 3)];

        write_entries(&path, &entries).unwrap();
        let loaded = read_entries(&path).unwrap();
        assert_eq!(entries, loaded);
    }

    #[test]
    fn test_extract_unique_hashes() {
        let entries = vec![
            Entry::new(100, 1),
            Entry::new(100, 2),
            Entry::new(100, 3),
            Entry::new(200, 1),
            Entry::new(300, 1),
            Entry::new(300, 2),
        ];

        let unique = extract_unique_hashes(&entries);
        assert_eq!(unique, vec![100, 200, 300]);
    }

    #[test]
    fn test_entry_writer() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("test.bin");

        let mut writer = EntryWriter::new(&path, 4096).unwrap();
        writer.write(&Entry::new(10, 1)).unwrap();
        writer.write(&Entry::new(20, 2)).unwrap();
        writer.flush().unwrap();

        assert_eq!(writer.count(), 2);

        let loaded = read_entries(&path).unwrap();
        assert_eq!(loaded, vec![Entry::new(10, 1), Entry::new(20, 2)]);
    }
}
