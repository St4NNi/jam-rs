use crate::core_utils::*;
use crate::hash_functions::ahash;
use crate::writer::create_lmdb_writer;
use anyhow::{Context, Result};
use byteorder::BigEndian;
use crossbeam_channel::Sender;
use heed::types::{Bytes, U32};
use heed::{Database, Env};
use needletail::{Sequence, parse_fastx_file, parser::SequenceRecord};
use rayon::prelude::*;
use std::borrow::Cow;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicU32, Ordering};

/// Configuration for sketching
#[derive(Debug, Clone)]
pub struct SketchConfig {
    pub kmer_size: u8,
    pub fscale: Option<u64>,
    pub nmax: Option<u64>,
    pub singleton: bool,
    pub min_entropy: f64,
    pub threads: usize,
    pub memory_budget_gb: f64,
}

impl Default for SketchConfig {
    fn default() -> Self {
        Self {
            kmer_size: 21,
            fscale: None,
            nmax: None,
            singleton: false,
            min_entropy: 1.5, // Default entropy threshold
            threads: 1,
            memory_budget_gb: 1.0,
        }
    }
}

/// Main sketching coordinator
pub struct Sketcher {
    config: SketchConfig,
    output_path: PathBuf,
    file_index_counter: Arc<AtomicU32>,
    metadata_db_path: PathBuf,
}

impl Sketcher {
    pub fn new(config: SketchConfig, output_path: PathBuf) -> Self {
        let metadata_db_path = output_path.with_extension("metadata.lmdb");

        Self {
            config,
            output_path,
            file_index_counter: Arc::new(AtomicU32::new(0)),
            metadata_db_path,
        }
    }

    /// Main entry point for sketching files
    pub fn sketch_files(&self, input_paths: &[PathBuf]) -> Result<()> {
        // Initialize metadata database
        let metadata_env = self.setup_metadata_database()?;

        // Create LMDB writer for hashes
        let (hash_sender, writer_handle) = create_lmdb_writer(
            self.config.memory_budget_gb,
            self.output_path.clone(),
            10000, // Channel capacity
        );

        // Set up rayon thread pool
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.config.threads)
            .build()
            .context("Failed to create thread pool")?;

        // Process files in parallel
        let results: Result<Vec<_>> = pool.install(|| {
            input_paths
                .par_iter()
                .map(|path| self.process_file(path, &hash_sender, &metadata_env))
                .collect()
        });

        // Close sender to signal completion
        drop(hash_sender);

        // Wait for writer to complete
        writer_handle
            .join()
            .map_err(|e| anyhow::anyhow!("Writer thread panicked: {:?}", e))??;

        // Check for any processing errors
        results?;

        println!("Sketching completed successfully!");
        Ok(())
    }

    /// Process a single file
    fn process_file(
        &self,
        file_path: &Path,
        hash_sender: &Sender<(u64, u64)>,
        metadata_env: &Env,
    ) -> Result<()> {
        let mut reader = parse_fastx_file(file_path)
            .with_context(|| format!("Failed to open file: {file_path:?}"))?;

        let mut sequence_count = 0;
        let mut total_length = 0;

        while let Some(record) = reader.next() {
            let record = record.context("Failed to parse sequence record")?;

            if self.config.singleton {
                // Each sequence gets its own file index
                let file_index = self.file_index_counter.fetch_add(1, Ordering::SeqCst);
                let sequence_length = record.seq().len();
                total_length += sequence_length;

                // Store metadata for this sequence
                self.store_sequence_metadata(
                    metadata_env,
                    file_index,
                    file_path,
                    &record,
                    sequence_length,
                    1, // Only one sequence per "file" in singleton mode
                )?;

                // Process this sequence
                self.process_sequence(&record, file_index, hash_sender)?;
                sequence_count += 1;
            } else {
                total_length += record.seq().len();
                sequence_count += 1;
            }
        }

        if !self.config.singleton && sequence_count > 0 {
            // All sequences in this file share one file index
            let file_index = self.file_index_counter.fetch_add(1, Ordering::SeqCst);

            // Store combined metadata for entire file
            self.store_file_metadata(
                metadata_env,
                file_index,
                file_path,
                sequence_count,
                total_length,
            )?;

            // Re-process file for hashing (we need the file_index)
            let mut reader = parse_fastx_file(file_path)?;
            while let Some(record) = reader.next() {
                let record = record?;
                self.process_sequence(&record, file_index, hash_sender)?;
            }
        }

        Ok(())
    }

    /// Process a single sequence record
    fn process_sequence(
        &self,
        record: &SequenceRecord,
        file_index: u32,
        hash_sender: &Sender<(u64, u64)>,
    ) -> Result<()> {
        let sequence = record.normalize(false);
        let seq_len = sequence.len();
        let gc_content = calculate_gc_content(sequence.as_ref());

        let metadata = HashMetadata::new(file_index, gc_content, seq_len);

        // Determine if we need to use a heap for nmax
        if let Some(nmax) = self.config.nmax {
            self.process_sequence_with_nmax(sequence, metadata, nmax as usize, hash_sender)
        } else {
            self.process_sequence_direct(sequence, metadata, hash_sender)
        }
    }

    /// Process sequence with nmax limit (using heap)
    fn process_sequence_with_nmax<'a>(
        &self,
        sequence: Cow<'a, [u8]>,
        metadata: HashMetadata,
        nmax: usize,
        hash_sender: &Sender<(u64, u64)>,
    ) -> Result<()> {
        let mut collector = TopNHashCollector::new(nmax);

        // Collect all valid hashes
        for (_, kmer, _) in sequence.bit_kmers(self.config.kmer_size, true) {
            // Apply entropy filter
            if !passes_entropy_filter(kmer.0, self.config.kmer_size, self.config.min_entropy) {
                continue;
            }

            let hash = ahash(kmer.0);

            // Apply FracMinHash filter if specified
            if let Some(fscale) = self.config.fscale
                && hash >= fscale
            {
                continue;
            }

            collector.add_hash(hash, metadata);
        }

        // Send the smallest nmax hashes
        for (hash, metadata) in collector.into_sorted_vec() {
            hash_sender.send((hash, metadata.pack()))?;
        }

        Ok(())
    }

    /// Process sequence without nmax limit (direct sending)
    fn process_sequence_direct<'a>(
        &self,
        sequence: Cow<'a, [u8]>,
        metadata: HashMetadata,
        hash_sender: &Sender<(u64, u64)>,
    ) -> Result<()> {
        for (_, kmer, _) in sequence.bit_kmers(self.config.kmer_size, true) {
            // Apply entropy filter
            if !passes_entropy_filter(kmer.0, self.config.kmer_size, self.config.min_entropy) {
                continue;
            }

            let hash = ahash(kmer.0);

            // Apply FracMinHash filter if specified
            if let Some(fscale) = self.config.fscale
                && hash >= fscale
            {
                continue;
            }

            hash_sender.send((hash, metadata.pack()))?;
        }

        Ok(())
    }

    /// Setup metadata database
    fn setup_metadata_database(&self) -> Result<Env> {
        std::fs::create_dir_all(self.metadata_db_path.parent().unwrap())?;

        let env = unsafe {
            heed::EnvOpenOptions::new()
                .map_size(1024 * 1024 * 1024) // 1GB for metadata
                .open(&self.metadata_db_path)?
        };

        Ok(env)
    }

    /// Store metadata for a single sequence (singleton mode)
    fn store_sequence_metadata(
        &self,
        env: &Env,
        file_index: u32,
        file_path: &Path,
        record: &SequenceRecord,
        sequence_length: usize,
        total_sequences: usize,
    ) -> Result<()> {
        let mut wtxn = env.write_txn()?;
        let db: Database<U32<BigEndian>, Bytes> =
            env.create_database(&mut wtxn, Some("METADATA"))?;

        let metadata = FileMetadata {
            filename: file_path.to_string_lossy().to_string(),
            file_size: std::fs::metadata(file_path)?.len(),
            sequence_name: String::from_utf8_lossy(record.id()).to_string(),
            sequence_length,
            total_sequences,
        };

        let metadata_json = serde_json::to_string(&metadata)?;
        db.put(&mut wtxn, &file_index, metadata_json.as_bytes())?;
        wtxn.commit()?;

        Ok(())
    }

    /// Store metadata for entire file (non-singleton mode)
    fn store_file_metadata(
        &self,
        env: &Env,
        file_index: u32,
        file_path: &Path,
        sequence_count: usize,
        total_length: usize,
    ) -> Result<()> {
        let mut wtxn = env.write_txn()?;
        let db: Database<U32<BigEndian>, Bytes> =
            env.create_database(&mut wtxn, Some("METADATA"))?;

        let metadata = FileMetadata {
            filename: file_path.to_string_lossy().to_string(),
            file_size: std::fs::metadata(file_path)?.len(),
            sequence_name: format!("{sequence_count} sequences"),
            sequence_length: total_length,
            total_sequences: sequence_count,
        };

        let metadata_json = serde_json::to_string(&metadata)?;
        db.put(&mut wtxn, &file_index, metadata_json.as_bytes())?;
        wtxn.commit()?;

        Ok(())
    }
}

/// Public interface for sketching
pub fn sketch_files(
    input_paths: &[PathBuf],
    output_path: PathBuf,
    config: SketchConfig,
) -> Result<()> {
    let sketcher = Sketcher::new(config, output_path);
    sketcher.sketch_files(input_paths)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn test_sketch_simple_sequence() {
        let temp_dir = tempdir().unwrap();
        let input_file = temp_dir.path().join("test.fasta");
        let output_file = temp_dir.path().join("test.lmdb");

        // Create a simple FASTA file
        let mut file = std::fs::File::create(&input_file).unwrap();
        writeln!(file, ">test_seq").unwrap();
        writeln!(file, "ATCGATCGATCGATCGATCGATCGATCGATCG").unwrap();

        let config = SketchConfig {
            kmer_size: 10,
            min_entropy: 0.5, // Low threshold for test
            ..Default::default()
        };

        let result = sketch_files(&[input_file], output_file.clone(), config);
        assert!(result.is_ok());
        assert!(output_file.exists());
    }
}
