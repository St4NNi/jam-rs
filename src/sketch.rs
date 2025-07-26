use crate::core_utils::*;
use crate::hash_functions::ahash;
use crate::writer::create_lmdb_writer;
use anyhow::{Context, Result};
use byteorder::BigEndian;
use crossbeam_channel::Sender;
use heed::types::{Bytes, U32};
use heed::{Database, Env};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use needletail::errors::ParseErrorKind;
use needletail::{Sequence, parse_fastx_file, parser::SequenceRecord};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicU32, Ordering};
use std::time::Instant;
use tempfile::TempDir;

/// Configuration for sketching
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SketchConfig {
    pub kmer_size: u8,
    pub fscale: u64,
    pub nmax: u64,
    pub singleton: bool,
    pub min_entropy: f64,
    pub threads: usize,
    pub memory_budget_gb: f64,
    pub temp_dir: Option<PathBuf>,
}

impl Default for SketchConfig {
    fn default() -> Self {
        Self {
            kmer_size: 21,
            fscale: u64::MAX,
            nmax: u64::MAX,
            singleton: false,
            min_entropy: 1.5, // Default entropy threshold
            threads: 1,
            memory_budget_gb: 1.0,
            temp_dir: None,
        }
    }
}

/// Main sketching coordinator
pub struct Sketcher {
    config: SketchConfig,
    output_path: PathBuf,
    file_index_counter: Arc<AtomicU32>,
    silent: bool,
}

impl Sketcher {
    pub fn new(config: SketchConfig, output_path: PathBuf, silent: bool) -> Self {
        Self {
            config,
            output_path,
            file_index_counter: Arc::new(AtomicU32::new(0)),
            silent,
        }
    }

    pub fn setup_tempdir(&self) -> Result<TempDir> {
        let tempdir = if let Some(temp_dir) = &self.config.temp_dir {
            tempfile::Builder::new()
                .prefix("temp_sketch_jam_")
                .tempdir_in(temp_dir)
                .context("Failed to create temporary directory")?
        } else {
            tempfile::Builder::new()
                .prefix("temp_sketch_jam_")
                .tempdir()
                .context("Failed to create temporary directory")?
        };
        Ok(tempdir)
    }

    /// Main entry point for sketching files
    pub fn sketch_files(&self, input_paths: &[PathBuf]) -> Result<()> {
        // Start timing the overall sketching process
        let start_time = Instant::now();

        // Initialize metadata database
        if !self.silent {
            eprintln!("Setting up LMDB environment...");
        }
        let tempdir = self.setup_tempdir()?;
        let env = self.setup_env(&tempdir)?;

        // Create progress tracking system only if not silent
        let (multi_progress, main_progress_bar, sub_progress_bar) = if !self.silent {
            let multi = MultiProgress::new();

            // Main progress bar for overall file progress
            let main_pb = multi.add(ProgressBar::new(input_paths.len() as u64));
            main_pb.set_style(
                ProgressStyle::default_bar()
                    .template("[{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} files ({percent}%) {msg}")
                    .unwrap()
                    .progress_chars("█▉▊▋▌▍▎▏  ")
            );
            main_pb.set_message("Starting...");

            // Sub-progress bar for current file progress (spinner style)
            let sub_pb = multi.add(ProgressBar::new_spinner());
            sub_pb.set_style(
                ProgressStyle::default_spinner()
                    .template("  └─ {spinner:.green} {pos} sequences processed")
                    .unwrap()
                    .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ "),
            );
            sub_pb.enable_steady_tick(std::time::Duration::from_millis(120));

            (Some(multi), Some(main_pb), Some(sub_pb))
        } else {
            (None, None, None)
        };

        let (hash_sender, writer_handle) = create_lmdb_writer(
            self.config.memory_budget_gb,
            env.clone(),
            10000, // Channel capacity
            self.config.temp_dir.clone(),
            self.silent,
            multi_progress,
        );

        // Set up rayon thread pool
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.config.threads)
            .build()
            .context("Failed to create thread pool")?;

        // Process files in parallel
        let results: Result<Vec<_>> = pool.install(|| {
            let iter = input_paths.par_iter();

            iter.map(|path| {
                let result = self.process_file(path, &hash_sender, &env, sub_progress_bar.as_ref());
                if let Some(ref pb) = main_progress_bar {
                    if let Err(ref e) = result {
                        pb.set_message(format!("Error in {}: {}", path.display(), e));
                    } else {
                        pb.inc(1);
                        pb.set_message(format!(
                            "Completed {}",
                            path.file_name().unwrap_or_default().to_string_lossy()
                        ));
                    }
                }
                result
            })
            .collect()
        });

        if let Some(pb) = main_progress_bar {
            pb.finish_with_message("File processing complete");
        }
        if let Some(pb) = sub_progress_bar {
            pb.finish_and_clear();
        }

        // Close sender to signal completion
        drop(hash_sender);

        // Wait for writer to complete
        writer_handle
            .join()
            .map_err(|e| anyhow::anyhow!("Writer thread panicked: {:?}", e))??;

        // Check for any processing errors
        results?;

        env.copy_to_path(&self.output_path, heed::CompactionOption::Enabled)
            .unwrap();

        env.prepare_for_closing().wait();

        // Display total runtime
        let elapsed = start_time.elapsed();
        if !self.silent {
            eprintln!("Sketching completed in {:.2}s", elapsed.as_secs_f64());
        }

        Ok(())
    }

    /// Process a single file
    fn process_file(
        &self,
        file_path: &Path,
        hash_sender: &Sender<(u64, u64)>,
        metadata_env: &Env,
        sub_progress_bar: Option<&ProgressBar>,
    ) -> Result<()> {
        let mut reader = match parse_fastx_file(file_path) {
            Ok(reader) => reader,
            Err(e) if e.kind == ParseErrorKind::EmptyFile => {
                eprintln!("Empty file detected: {}, skipping", file_path.display());
                // Finalize sub-progress bar for this file
                if let Some(sub_pb) = sub_progress_bar {
                    sub_pb.set_message(format!("Skipped empty file: {}", file_path.display()));
                }
                return Ok(());
            }
            Err(e) => {
                eprintln!("Error parsing file {}: {}", file_path.display(), e);
                return Err(anyhow::anyhow!(
                    "Failed to parse file: {}",
                    file_path.display()
                ));
            }
        };

        let mut sequence_count = 0;
        let mut total_length = 0;
        let mut total_hashes = 0;

        let mut file_index = self.file_index_counter.fetch_add(1, Ordering::SeqCst);

        let mut sequence_counter = 0u64;
        while let Some(record) = reader.next() {
            let record = record.context("Failed to parse sequence record")?;
            sequence_counter += 1;

            // Update sub-progress bar with current sequence info
            if let Some(sub_pb) = sub_progress_bar {
                sub_pb.inc(1);
            }

            let id = String::from_utf8_lossy(record.id()).to_string();
            let (hashes, seq_len) = self.process_sequence(record, file_index, hash_sender)?;
            total_length += seq_len;
            total_hashes += hashes;

            if self.config.singleton {
                // Each sequence gets its own file index
                let sequence_length = seq_len;
                total_length += sequence_length;
                // Store metadata for this sequence
                self.store_sequence_metadata(
                    metadata_env,
                    file_index,
                    file_path,
                    &id,
                    sequence_length,
                    1, // Only one sequence per "file" in singleton mode
                    total_hashes,
                )?;
                file_index = self.file_index_counter.fetch_add(1, Ordering::SeqCst);
                total_hashes = 0; // Reset for next sequence
            } else {
                sequence_count += 1;
            }
        }

        // Finalize sub-progress bar for this file
        if let Some(sub_pb) = sub_progress_bar {
            sub_pb.set_message(format!("Completed {} sequences", sequence_counter));
        }

        if !self.config.singleton && sequence_count > 0 {
            // Store combined metadata for entire file
            self.store_file_metadata(
                metadata_env,
                file_index,
                file_path,
                sequence_count,
                total_length,
                total_hashes,
            )?;
        }

        Ok(())
    }

    /// Process a single sequence record
    fn process_sequence(
        &self,
        record: SequenceRecord,
        file_index: u32,
        hash_sender: &Sender<(u64, u64)>,
    ) -> Result<(usize, usize)> {
        let sequence = record.normalize(false);
        let seq_len = sequence.len();
        let gc_content = calculate_gc_content(sequence.as_ref());

        let metadata = HashMetadata::new(file_index, gc_content, seq_len);

        // Determine if we need to use a heap for nmax
        if self.config.nmax != u64::MAX {
            self.process_sequence_with_nmax(
                sequence,
                metadata,
                self.config.nmax as usize,
                hash_sender,
            )
            .map(|count| (count, seq_len))
        } else {
            self.process_sequence_direct(sequence, metadata, hash_sender)
                .map(|count| (count, seq_len))
        }
    }

    /// Process sequence with nmax limit (using heap)
    fn process_sequence_with_nmax<'a>(
        &self,
        sequence: Cow<'a, [u8]>,
        metadata: HashMetadata,
        nmax: usize,
        hash_sender: &Sender<(u64, u64)>,
    ) -> Result<usize> {
        let mut collector = TopNHashCollector::new(nmax);

        let mut total_hashes = 0usize;

        // Collect all valid hashes
        for (_, kmer, _) in sequence.bit_kmers(self.config.kmer_size, true) {
            // Apply entropy filter
            if !passes_entropy_filter(kmer.0, self.config.kmer_size, self.config.min_entropy) {
                continue;
            }

            let hash = ahash(kmer.0);

            // Apply FracMinHash filter if specified
            if hash > self.config.fscale {
                continue;
            }
            collector.add_hash(hash, metadata);
        }

        // Send the smallest nmax hashes
        for (hash, metadata) in collector.into_sorted_vec() {
            total_hashes += 1;
            hash_sender.send((hash, metadata.pack()))?;
        }

        Ok(total_hashes)
    }

    /// Process sequence without nmax limit (direct sending)
    fn process_sequence_direct<'a>(
        &self,
        sequence: Cow<'a, [u8]>,
        metadata: HashMetadata,
        hash_sender: &Sender<(u64, u64)>,
    ) -> Result<usize> {
        let mut total_hashes = 0;
        for (_, kmer, _) in sequence.bit_kmers(self.config.kmer_size, true) {
            // Apply entropy filter
            if !passes_entropy_filter(kmer.0, self.config.kmer_size, self.config.min_entropy) {
                continue;
            }

            let hash = ahash(kmer.0);

            // Apply FracMinHash filter if specified
            if hash > self.config.fscale {
                continue;
            }

            total_hashes += 1;
            hash_sender.send((hash, metadata.pack()))?;
        }

        Ok(total_hashes)
    }

    /// Setup metadata database
    fn setup_env(&self, tempdir: &TempDir) -> Result<Env> {
        let mut path = tempdir.path().to_path_buf();
        path.push("sketch.lmdb");

        let env = unsafe {
            heed::EnvOpenOptions::new()
                .flags(
                    heed::EnvFlags::NO_SUB_DIR
                        | heed::EnvFlags::WRITE_MAP
                        | heed::EnvFlags::MAP_ASYNC,
                )
                .max_dbs(3)
                .map_size(10 * 1024 * 1024 * 1024 * 1024) // 10TB map size
                .open(path)?
        };

        // Create config database
        let mut wtxn = env.write_txn()?;
        let config_db: Database<Bytes, Bytes> = env.create_database(&mut wtxn, Some("CONFIG"))?;
        let config_json = serde_json::to_string(&self.config)?;
        config_db.put(&mut wtxn, b"config", config_json.as_bytes())?;
        wtxn.commit()?;

        Ok(env)
    }

    /// Store metadata for a single sequence (singleton mode)
    fn store_sequence_metadata(
        &self,
        env: &Env,
        file_index: u32,
        file_path: &Path,
        id: &str,
        sequence_length: usize,
        total_sequences: usize,
        total_hashes: usize,
    ) -> Result<()> {
        let mut wtxn = env.write_txn()?;
        let db: Database<U32<BigEndian>, Bytes> =
            env.create_database(&mut wtxn, Some("METADATA"))?;

        let metadata = FileMetadata {
            filename: file_path.to_string_lossy().to_string(),
            file_size: std::fs::metadata(file_path)?.len(),
            sequence_name: id.to_string(),
            sequence_length,
            total_sequences,
            total_hashes,
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
        total_hashes: usize,
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
            total_hashes,
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
    silent: bool,
) -> Result<()> {
    let sketcher = Sketcher::new(config, output_path, silent);
    sketcher.sketch_files(input_paths)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_sketch_simple_sequence() {
        let mut input_file = NamedTempFile::new().unwrap();
        let output_file = NamedTempFile::new().unwrap().into_temp_path().to_path_buf();

        // Create a simple FASTA file
        writeln!(input_file, ">test_seq").unwrap();
        writeln!(input_file, "ATCGATCGATCGATCGATCGATCGATCGATCG").unwrap();

        let config = SketchConfig {
            kmer_size: 10,
            min_entropy: 0.5, // Low threshold for test
            ..Default::default()
        };

        sketch_files(
            &[input_file.into_temp_path().to_path_buf()],
            output_file.clone(),
            config,
            false, // Not silent for tests
        )
        .unwrap();
        assert!(output_file.exists());
    }
}
