use crate::core_utils::*;
use anyhow::{Context, Result};
use byteorder::BigEndian;
use heed::types::{Bytes, U32, U64};
use heed::{Database, EnvFlags};
use std::collections::HashMap;
use std::path::Path;

/// Database statistics
#[derive(Debug)]
pub struct DatabaseStats {
    pub total_hashes: u64,
    pub unique_files: u32,
    pub gc_distribution: HashMap<u16, u64>,
    pub length_distribution: HashMap<u16, u64>,
    pub file_metadata: Vec<FileMetadata>,
}

impl DatabaseStats {
    /// Display short summary statistics
    pub fn print_short(&self) {
        eprintln!("Database Statistics:");
        eprintln!("Total hashes: {}", self.total_hashes);
        eprintln!("Unique files/sequences: {}", self.unique_files);
    }

    // Prints a csv summary of all sequence entries
    pub fn print_full(&self) {
        for md in &self.file_metadata {
            eprintln!(
                "{},{},{},{},{}",
                md.filename, md.file_size, md.sequence_name, md.sequence_length, md.total_sequences
            );
        }
    }

    /// Display detailed statistics
    pub fn print_detailed(&self) {
        self.print_short();

        eprintln!("\nFile/Sequence Information:");

        let mut sum_lengths: usize = 0;
        let mut sum_sequences: usize = 0;
        for md in self.file_metadata.iter() {
            sum_lengths += md.sequence_length;
            sum_sequences += md.total_sequences;
        }

        eprintln!(
            "Total files/sequences: {} ({} bp, {} sequences)",
            self.file_metadata.len(),
            sum_lengths,
            sum_sequences
        );

        for (i, metadata) in self.file_metadata.iter().take(5).enumerate() {
            eprintln!(
                "  {}. {}/{} ({} bp, {} sequences)",
                i + 1,
                metadata.filename,
                metadata.sequence_name,
                metadata.sequence_length,
                metadata.total_sequences
            );
        }

        if self.file_metadata.len() > 5 {
            eprintln!(
                "  ... and {} more files/sequences",
                self.file_metadata.len() - 5
            );
        }
    }
}

/// Statistics calculator
pub struct StatsCalculator;

impl StatsCalculator {
    /// Calculate statistics for an LMDB database
    pub fn calculate_stats(database_path: &Path) -> Result<DatabaseStats> {
        // Open hash database
        let env = unsafe {
            heed::EnvOpenOptions::new()
                .flags(EnvFlags::NO_SUB_DIR | EnvFlags::READ_ONLY)
                .max_dbs(3)
                .open(database_path)
                .with_context(|| format!("Failed to open database: {database_path:?}"))?
        };

        let rtxn = env.read_txn()?;
        let hash_db: Database<U64<BigEndian>, U64<BigEndian>> = env
            .open_database(&rtxn, Some("HASHES"))?
            .context("Hash database not found")?;

        // Open metadata database
        let metadata_db: Database<U32<BigEndian>, Bytes> = env
            .open_database(&rtxn, Some("METADATA"))?
            .context("Metadata database not found")?;

        let mut stats = DatabaseStats {
            total_hashes: hash_db.len(&rtxn)?,
            unique_files: 0,
            gc_distribution: HashMap::new(),
            length_distribution: HashMap::new(),
            file_metadata: Vec::new(),
        };

        // Load file metadata
        for item in metadata_db.iter(&rtxn)? {
            let (file_index, metadata_json) = item?;

            if let Ok(file_metadata) = serde_json::from_slice::<FileMetadata>(metadata_json) {
                stats.file_metadata.push(file_metadata);
            }

            stats.unique_files = stats.unique_files.max(file_index + 1);
        }

        Ok(stats)
    }

    /// Print statistics for a database
    pub fn print_stats(database_path: &Path, short: bool, full: bool) -> Result<()> {
        let stats = Self::calculate_stats(database_path)?;

        if full {
            stats.print_full();
            return Ok(());
        }

        if short {
            stats.print_short();
        } else {
            stats.print_detailed();
        }

        Ok(())
    }
}
