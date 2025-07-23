use anyhow::Result;
use byteorder::BigEndian;
use crossbeam_channel::Receiver;
use heed::types::U64;
use heed::{Database, DatabaseFlags, Env, IntegerComparator, PutFlags};
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use tempfile::NamedTempFile;

pub struct MergeIterator {
    heap: BinaryHeap<Reverse<(u64, u64, usize)>>,
    taken: Vec<usize>,
    readers: Vec<Option<BufReader<File>>>,
    batch_size: usize,
}

impl MergeIterator {
    pub fn new<P: AsRef<Path>>(file_paths: Vec<P>, memory_limit_gb: f64) -> Self {
        let num_files = file_paths.len();

        // Calculate batch size: memory_limit / (num_files * 2 * 24 bytes per heap entry)
        // Each heap entry is approximately 24 bytes (2*u64 + usize)
        let memory_limit_bytes = (memory_limit_gb * 1024.0 * 1024.0 * 1024.0) as usize;
        let batch_size = std::cmp::max(1000, memory_limit_bytes / (num_files * 2 * 24));

        let mut readers: Vec<Option<BufReader<File>>> = Vec::with_capacity(num_files);
        let mut heap = BinaryHeap::new();
        let taken = vec![0; num_files];

        // Initialize readers and load initial batches
        for (file_idx, path) in file_paths.into_iter().enumerate() {
            let file = File::open(path).expect("Failed to open file");
            let mut reader = BufReader::new(file);

            // Load initial 2 * batch_size entries (or all if file is smaller)
            let loaded = Self::load_batch(&mut reader, &mut heap, file_idx, 2 * batch_size);

            if loaded > 0 {
                readers.push(Some(reader));
            } else {
                readers.push(None); // Empty file
            }
        }

        MergeIterator {
            heap,
            taken,
            readers,
            batch_size,
        }
    }

    fn load_batch(
        reader: &mut BufReader<File>,
        heap: &mut BinaryHeap<Reverse<(u64, u64, usize)>>,
        file_idx: usize,
        max_count: usize,
    ) -> usize {
        let mut buffer = [0u8; 16]; // 2 * u64 = 16 bytes
        let mut loaded = 0;

        for _ in 0..max_count {
            match reader.read_exact(&mut buffer) {
                Ok(()) => {
                    let first = u64::from_be_bytes([
                        buffer[0], buffer[1], buffer[2], buffer[3], buffer[4], buffer[5],
                        buffer[6], buffer[7],
                    ]);
                    let second = u64::from_be_bytes([
                        buffer[8], buffer[9], buffer[10], buffer[11], buffer[12], buffer[13],
                        buffer[14], buffer[15],
                    ]);

                    heap.push(Reverse((first, second, file_idx)));
                    loaded += 1;
                }
                Err(_) => break, // End of file or read error
            }
        }

        loaded
    }

    fn reload_if_needed(&mut self, file_idx: usize) {
        if self.taken[file_idx] >= self.batch_size
            && let Some(ref mut reader) = self.readers[file_idx]
        {
            let loaded = Self::load_batch(reader, &mut self.heap, file_idx, self.batch_size);

            if loaded == 0 {
                // File is exhausted, remove reader
                self.readers[file_idx] = None;
            }

            // Reset taken counter
            self.taken[file_idx] = 0;
        }
    }
}

impl Iterator for MergeIterator {
    type Item = (u64, u64);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(Reverse((first, second, file_idx))) = self.heap.pop() {
            // Increment taken counter for this file
            self.taken[file_idx] += 1;

            // Check if we need to reload from this file
            self.reload_if_needed(file_idx);

            Some((first, second))
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
}

impl LMDBWriter {
    pub fn new(memory_budget: usize, env: Env, receiver: Receiver<(u64, u64)>) -> Self {
        Self {
            memory_budget,
            env,
            receiver,
            current_chunk: BinaryHeap::new(),
            current_memory_usage: 0,
            chunk_files: Vec::new(),
            total_entries: 0,
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
        let temp_file = NamedTempFile::new()?;
        let mut writer = BufWriter::new(temp_file);

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

        if self.chunk_files.is_empty() {
            // Everything fits in memory - write current chunk directly
            self.write_current_chunk_to_lmdb(&db)?;
        } else {
            // Need to merge temp files with current chunk
            self.merge_all_to_lmdb(&db)?;
        }

        Ok(())
    }

    fn write_current_chunk_to_lmdb(
        &mut self,
        db: &Database<U64<BigEndian>, U64<BigEndian>, IntegerComparator>,
    ) -> Result<()> {
        let mut wtxn = self.env.write_txn()?;
        let mut batch_count = 0;

        // Write current chunk in sorted order
        let mut prev = (0u64, 0u64);
        while let Some(Reverse((hash, metadata))) = self.current_chunk.pop() {
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
        Ok(())
    }

    fn merge_all_to_lmdb(
        &mut self,
        db: &Database<U64<BigEndian>, U64<BigEndian>, IntegerComparator>,
    ) -> Result<()> {
        // If current chunk has data, write it to a temp file first
        if !self.current_chunk.is_empty() {
            self.flush_chunk_to_tempfile()?;
        }

        // Collect all temp file paths
        let file_paths: Vec<PathBuf> = self
            .chunk_files
            .iter()
            .map(|tf| tf.path().to_path_buf())
            .collect();

        // Calculate memory limit for MergeIterator (convert from bytes to GB)
        let memory_limit_gb = (self.memory_budget as f64) / (1024.0 * 1024.0 * 1024.0);

        // Create iterator to merge all files
        let iterator = MergeIterator::new(file_paths, memory_limit_gb);

        // Write merged data to LMDB
        let mut wtxn = self.env.write_txn()?;
        let mut batch_count = 0;

        let mut prev = (0u64, 0u64);
        for (hash, metadata) in iterator {
            if prev.0 == hash && prev.1 == metadata {
                // Skip duplicates
                continue;
            }
            prev = (hash, metadata);
            db.put_with_flags(&mut wtxn, PutFlags::APPEND_DUP, &hash, &metadata)?;
            batch_count += 1;

            // Commit in batches for performance
            if batch_count >= 10000 {
                wtxn.commit()?;
                wtxn = self.env.write_txn()?;
                batch_count = 0;
            }
        }

        wtxn.commit()?;
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
) -> SenderWithHandle {
    let memory_budget = (memory_budget_gb * 1024.0 * 1024.0 * 1024.0) as usize;
    let (sender, receiver) = crossbeam_channel::bounded(channel_capacity);

    let writer = LMDBWriter::new(memory_budget, env, receiver);

    let handle = std::thread::spawn(move || writer.run());

    (sender, handle)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    fn create_env(path: PathBuf) -> Env {
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
        );

        // Send test data
        let test_data = vec![(300, 3), (100, 1), (200, 2), (500, 5), (400, 4)];

        for (hash, metadata) in test_data {
            sender.send((hash, metadata)).unwrap();
        }

        println!("Sent items");

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
            env, 1000,
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
