use anyhow::Result;
use byteorder::{BigEndian, ReadBytesExt};
use crossbeam_channel::Receiver;
use heed::types::U64;
use heed::{Database, DatabaseFlags, Env, IntegerComparator, PutFlags};
use indicatif::{ProgressBar, ProgressStyle};
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use tempfile::NamedTempFile;

pub struct MergeIterator {
    heap: BinaryHeap<Reverse<(u64, u64, usize)>>,
    taken: Vec<usize>,
    readers: Vec<Option<BufReader<File>>>,
    batch_size: usize,
}

impl MergeIterator {
    pub fn new<P: AsRef<Path>>(file_paths: Vec<P>, memory_limit_gb: f64) -> Result<Self> {
        let num_files = file_paths.len();

        let memory_limit_bytes = (memory_limit_gb * 1024.0 * 1024.0 * 1024.0) as usize;
        let batch_size = std::cmp::max(1000, memory_limit_bytes / (num_files * 2 * 24));

        let mut readers: Vec<Option<BufReader<File>>> = Vec::with_capacity(num_files);
        let mut heap = BinaryHeap::new();
        let taken = vec![0; num_files];

        for (file_idx, path) in file_paths.into_iter().enumerate() {
            let file = File::open(path.as_ref())
                .map_err(|e| anyhow::anyhow!("Failed to open temp file {:?}: {}", path.as_ref(), e))?;
            let mut reader = BufReader::new(file);

            let loaded = Self::load_batch(&mut reader, &mut heap, file_idx, 2 * batch_size)?;

            if loaded > 0 {
                readers.push(Some(reader));
            } else {
                readers.push(None);
            }
        }

        Ok(MergeIterator {
            heap,
            taken,
            readers,
            batch_size,
        })
    }

    fn load_batch(
        reader: &mut BufReader<File>,
        heap: &mut BinaryHeap<Reverse<(u64, u64, usize)>>,
        file_idx: usize,
        max_count: usize,
    ) -> Result<usize> {
        let mut loaded = 0;

        for _ in 0..max_count {
            let first = match reader.read_u64::<BigEndian>() {
                Ok(val) => val,
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(anyhow::anyhow!("Failed to read hash from temp file {}: {}", file_idx, e)),
            };
            let second = reader.read_u64::<BigEndian>()
                .map_err(|e| anyhow::anyhow!("Failed to read metadata from temp file {}: {}", file_idx, e))?;
            heap.push(Reverse((first, second, file_idx)));
            loaded += 1;
        }

        Ok(loaded)
    }

    fn reload_if_needed(&mut self, file_idx: usize) -> Result<()> {
        if self.taken[file_idx] >= self.batch_size
            && let Some(ref mut reader) = self.readers[file_idx]
        {
            let loaded = Self::load_batch(reader, &mut self.heap, file_idx, self.batch_size)?;

            if loaded == 0 {
                self.readers[file_idx] = None;
            }

            self.taken[file_idx] = 0;
        }
        Ok(())
    }
}

impl Iterator for MergeIterator {
    type Item = Result<(u64, u64)>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(Reverse((first, second, file_idx))) = self.heap.pop() {
            self.taken[file_idx] += 1;

            if let Err(e) = self.reload_if_needed(file_idx) {
                return Some(Err(e));
            }

            Some(Ok((first, second)))
        } else {
            None
        }
    }
}

// Main LMDBWriter implementation
pub struct LMDBWriter {
    memory_budget: usize,
    env: Env,
    receiver: Receiver<(u64, u64)>,
    current_chunk: BinaryHeap<Reverse<(u64, u64)>>,
    current_memory_usage: usize,
    chunk_files: Vec<NamedTempFile>,
    total_entries: u64,
    temp_dir: Option<PathBuf>,
    silent: bool,
    bar: Option<indicatif::MultiProgress>,
}

impl LMDBWriter {
    pub fn new(
        memory_budget: usize,
        env: Env,
        receiver: Receiver<(u64, u64)>,
        temp_dir: Option<PathBuf>,
        silent: bool,
        bar: Option<indicatif::MultiProgress>,
    ) -> Self {
        Self {
            memory_budget,
            env,
            receiver,
            current_chunk: BinaryHeap::new(),
            current_memory_usage: 0,
            chunk_files: Vec::new(),
            total_entries: 0,
            temp_dir,
            silent,
            bar,
        }
    }

    pub fn run(mut self) -> Result<()> {
        // Process all entries from the channel
        while let Ok((hash, metadata)) = self.receiver.recv() {
            self.current_chunk.push(Reverse((hash, metadata)));
            self.current_memory_usage += 16; // 8 bytes hash + 8 bytes metadata
            self.total_entries += 1;

            // Check if we need to flush to temp file
            if self.current_memory_usage >= self.memory_budget {
                self.flush_chunk_to_tempfile()?;
            }
        }

        // Channel closed, finalize to LMDB
        self.finalize_to_lmdb()
    }

    fn flush_chunk_to_tempfile(&mut self) -> Result<()> {
        let temp_file = if let Some(ref temp_dir) = self.temp_dir {
            NamedTempFile::new_in(temp_dir)?
        } else {
            NamedTempFile::new()?
        };

        if !self.silent {
            eprintln!(
                "Creating temporary sort file: {} ({} entries, {:.1} MB)",
                temp_file.path().display(),
                self.current_chunk.len(),
                self.current_memory_usage as f64 / (1024.0 * 1024.0)
            );
        }

        let mut writer = BufWriter::with_capacity(50 * 1024 * 1024, temp_file); // 50 MiB buffer

        // Write entries in sorted order (BinaryHeap with Reverse gives us min-heap)
        while let Some(Reverse((hash, metadata))) = self.current_chunk.pop() {
            writer.write_all(&hash.to_be_bytes())?;
            writer.write_all(&metadata.to_be_bytes())?;
        }

        writer.flush()?;
        let temp_file = writer.into_inner()?;
        self.chunk_files.push(temp_file);
        self.current_memory_usage = 0;
        Ok(())
    }

    fn finalize_to_lmdb(mut self) -> Result<()> {
        // Create LMDB environment
        let mut write_txn = self.env.write_txn()?;
        let db = self
            .env
            .database_options()
            .types::<U64<BigEndian>, U64<BigEndian>>()
            .flags(DatabaseFlags::DUP_SORT | DatabaseFlags::DUP_FIXED)
            .name("HASHES")
            .key_comparator::<IntegerComparator>()
            .create(&mut write_txn)?;

        write_txn.commit()?;

        let option_bar = if let Some(ref bar) = self.bar {
            let sub_bar = ProgressBar::new(self.total_entries);
            sub_bar.set_style(
                ProgressStyle::default_bar()
                    .template("[{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} kmers processed ({percent}%)")
                    .unwrap()
                    .progress_chars("█▉▊▋▌▍▎▏  ")
            );

            Some(bar.add(sub_bar))
        } else {
            None
        };

        if self.chunk_files.is_empty() {
            // Everything fits in memory - write current chunk directly
            self.write_current_chunk_to_lmdb(&db, option_bar)?;
        } else {
            // Need to merge temp files with current chunk
            self.merge_all_to_lmdb(&db, option_bar)?;
        }

        Ok(())
    }

    fn write_current_chunk_to_lmdb(
        &mut self,
        db: &Database<U64<BigEndian>, U64<BigEndian>, IntegerComparator>,
        option_bar: Option<ProgressBar>,
    ) -> Result<()> {
        let mut wtxn = self.env.write_txn()?;
        let mut batch_count = 0;

        // Write current chunk in sorted order
        let mut prev = (0u64, 0u64);
        while let Some(Reverse((hash, metadata))) = self.current_chunk.pop() {
            if let Some(ref bar) = option_bar {
                bar.inc(1);
            }

            if prev.0 == hash && prev.1 == metadata {
                // Skip duplicates
                continue;
            }
            prev = (hash, metadata);

            db.put_with_flags(&mut wtxn, PutFlags::APPEND_DUP, &hash, &metadata)?;
            batch_count += 1;

            // Commit in batches for performance
            if batch_count >= 100000 {
                wtxn.commit()?;
                wtxn = self.env.write_txn()?;
                batch_count = 0;
            }
        }

        wtxn.commit()?;

        if let Some(ref bar) = option_bar {
            bar.finish_with_message("Data merged successfully");
        }

        Ok(())
    }

    fn merge_all_to_lmdb(
        &mut self,
        db: &Database<U64<BigEndian>, U64<BigEndian>, IntegerComparator>,
        option_bar: Option<ProgressBar>,
    ) -> Result<()> {
        // If current chunk has data, write it to a temp file first
        if !self.current_chunk.is_empty() {
            self.flush_chunk_to_tempfile()?;
        }

        let num_temp_files = self.chunk_files.len();
        if num_temp_files > 0 {
            if !self.silent {
                eprintln!(
                    "Merging {} temporary files into LMDB database...",
                    num_temp_files
                );
            }
        } else if !self.silent {
            eprintln!("Writing hash data directly to LMDB (no temporary files needed)...");
        }

        // Collect all temp file paths
        let file_paths: Vec<PathBuf> = self
            .chunk_files
            .iter()
            .map(|tf| tf.path().to_path_buf())
            .collect();

        // Calculate memory limit for MergeIterator (convert from bytes to GB)
        let memory_limit_gb = (self.memory_budget as f64) / (1024.0 * 1024.0 * 1024.0);

        let iterator = MergeIterator::new(file_paths, memory_limit_gb)?;

        let mut wtxn = self.env.write_txn()?;
        let mut batch_count = 0;
        let mut total_written = 0u64;
        let mut duplicates_skipped = 0u64;

        let mut prev = (0u64, 0u64);
        for result in iterator {
            let (hash, metadata) = result?;
            if let Some(ref bar) = option_bar {
                bar.inc(1);
            }
            if prev.0 == hash && prev.1 == metadata {
                duplicates_skipped += 1;
                continue;
            }
            prev = (hash, metadata);
            db.put_with_flags(&mut wtxn, PutFlags::APPEND_DUP, &hash, &metadata)?;
            batch_count += 1;
            total_written += 1;

            if batch_count >= 10000 {
                wtxn.commit()?;
                wtxn = self.env.write_txn()?;
                batch_count = 0;
            }
        }

        if let Some(ref bar) = option_bar {
            bar.finish_with_message("Data merged successfully");
        }

        wtxn.commit()?;

        if !self.silent {
            if num_temp_files > 0 {
                eprintln!(
                    "Merge complete! {} unique hashes written, {} duplicates skipped",
                    total_written, duplicates_skipped
                );
            } else {
                eprintln!(
                    "Direct write complete! {} unique hashes written",
                    total_written
                );
            }
        }

        Ok(())
    }
}

pub type SenderWithHandle = (
    crossbeam_channel::Sender<(u64, u64)>,
    std::thread::JoinHandle<Result<(), anyhow::Error>>,
);

// Usage example
pub fn create_lmdb_writer(
    memory_budget_gb: f64,
    env: Env,
    channel_capacity: usize,
    temp_dir: Option<PathBuf>,
    silent: bool,
    bar: Option<indicatif::MultiProgress>,
) -> SenderWithHandle {
    let memory_budget = (memory_budget_gb * 1024.0 * 1024.0 * 1024.0) as usize;
    let (sender, receiver) = crossbeam_channel::bounded(channel_capacity);

    let writer = LMDBWriter::new(memory_budget, env, receiver, temp_dir, silent, bar);

    let handle = std::thread::spawn(move || writer.run());

    (sender, handle)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    fn create_env(path: PathBuf) -> Env {
        // SAFETY: heed requires unsafe for mmap; we control file access
        unsafe {
            heed::EnvOpenOptions::new()
                .flags(
                    heed::EnvFlags::NO_SUB_DIR
                        | heed::EnvFlags::MAP_ASYNC
                        | heed::EnvFlags::NO_SYNC,
                )
                .max_dbs(3)
                .map_size(10 * 1024 * 1024 * 1024) // 10GB map size
                .open(&path)
                .unwrap()
        }
    }

    #[test]
    fn test_lmdb_writer_memory_only() {
        let temp_dir = tempdir().unwrap();
        let lmdb_path = temp_dir.path().join("test.lmdb");

        let env = create_env(lmdb_path.clone());

        let (sender, handle) = create_lmdb_writer(
            1.0, // 1GB memory budget
            env, 1000, // channel capacity
            None, // Use default temp directory
            true, // Silent mode for tests
            None, // No progress bar in tests
        );

        // Send test data
        let test_data = vec![(300, 3), (100, 1), (200, 2), (500, 5), (400, 4)];

        for (hash, metadata) in test_data {
            sender.send((hash, metadata)).unwrap();
        }

        drop(sender); // Close channel

        // Wait for completion
        handle.join().unwrap().unwrap();

        // Verify data was written (basic check)
        assert!(lmdb_path.exists());
    }

    #[test]
    fn test_lmdb_writer_with_temp_files() {
        let temp_dir = tempdir().unwrap();
        let lmdb_path = temp_dir.path().join("test.lmdb");

        let env = create_env(lmdb_path.clone());

        let (sender, handle) = create_lmdb_writer(
            0.001, // Very small memory budget to force temp files
            env, 1000, None, // Use default temp directory
            true, // Silent mode for tests
            None, // No progress bar in tests
        );

        // Send enough data to exceed memory budget
        for i in 0..1000 {
            sender.send((i * 2, i)).unwrap(); // Ensure some ordering
        }

        drop(sender);
        handle.join().unwrap().unwrap();

        assert!(lmdb_path.exists());
    }
}
