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
        println!("Database Statistics:");
        println!("  Total hashes: {}", self.total_hashes);
        println!("  Unique files/sequences: {}", self.unique_files);
    }

    /// Display detailed statistics
    pub fn print_detailed(&self) {
        self.print_short();

        println!("\nFile/Sequence Information:");
        for (i, metadata) in self.file_metadata.iter().take(5).enumerate() {
            println!(
                "  {}: {} ({} bp, {} sequences)",
                i + 1,
                metadata.sequence_name,
                metadata.sequence_length,
                metadata.total_sequences
            );
        }

        if self.file_metadata.len() > 5 {
            println!(
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
    pub fn calculate_stats(database_path: &Path, _detailed: bool) -> Result<DatabaseStats> {
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
            total_hashes: hash_db.len(&rtxn)? as u64,
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
    pub fn print_stats(database_path: &Path, short: bool) -> Result<()> {
        let stats = Self::calculate_stats(database_path, !short)?;

        if short {
            stats.print_short();
        } else {
            stats.print_detailed();
        }

        Ok(())
    }
}
