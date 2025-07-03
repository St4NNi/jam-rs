use crate::core_utils::*;
use anyhow::{Context, Result};
use byteorder::BigEndian;
use heed::Database;
use heed::types::{Bytes, U32, U64};
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

        if !self.gc_distribution.is_empty() {
            let avg_gc = self.calculate_average_gc();
            println!("  Average GC content: {:.1}%", avg_gc);
        }

        if !self.length_distribution.is_empty() {
            let avg_length = self.calculate_average_length();
            println!("  Average sequence length: {:.0} bp", avg_length);
        }
    }

    /// Display detailed statistics
    pub fn print_detailed(&self) {
        self.print_short();

        println!("\nGC Content Distribution:");
        let mut gc_entries: Vec<_> = self.gc_distribution.iter().collect();
        gc_entries.sort_by_key(|&(category, _)| category);

        for (category, count) in gc_entries.iter().take(10) {
            let gc_range = self.gc_category_to_range(**category);
            println!("  {}: {} hashes", gc_range, count);
        }

        if gc_entries.len() > 10 {
            println!("  ... and {} more categories", gc_entries.len() - 10);
        }

        println!("\nLength Distribution:");
        let mut length_entries: Vec<_> = self.length_distribution.iter().collect();
        length_entries.sort_by_key(|&(category, _)| category);

        for (category, count) in length_entries.iter().take(10) {
            let length_range = self.length_category_to_range(**category);
            println!("  {}: {} hashes", length_range, count);
        }

        if length_entries.len() > 10 {
            println!("  ... and {} more categories", length_entries.len() - 10);
        }

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

    fn calculate_average_gc(&self) -> f64 {
        let mut total_weighted_gc = 0.0;
        let mut total_count = 0;

        for (&category, &count) in &self.gc_distribution {
            let gc_midpoint = self.gc_category_to_midpoint(category);
            total_weighted_gc += gc_midpoint * count as f64;
            total_count += count;
        }

        if total_count > 0 {
            total_weighted_gc / total_count as f64
        } else {
            0.0
        }
    }

    fn calculate_average_length(&self) -> f64 {
        let mut total_weighted_length = 0.0;
        let mut total_count = 0;

        for (&category, &count) in &self.length_distribution {
            let length_midpoint = self.length_category_to_midpoint(category);
            total_weighted_length += length_midpoint * count as f64;
            total_count += count;
        }

        if total_count > 0 {
            total_weighted_length / total_count as f64
        } else {
            0.0
        }
    }

    fn gc_category_to_range(&self, category: u16) -> String {
        match category {
            0..=5 => format!("{:.0}-{:.0}%", category * 5, (category + 1) * 5),
            6..=45 => {
                let start = 30 + (category - 6);
                format!("{}-{}%", start, start + 1)
            }
            46..=51 => {
                let start = 70 + (category - 46) * 5;
                format!("{:.0}-{:.0}%", start, start + 5)
            }
            _ => format!("Unknown ({})", category),
        }
    }

    fn length_category_to_range(&self, category: u16) -> String {
        match category {
            0 => "< 1 kB".to_string(),
            1..=50 => {
                let start = 1_000 + (category as usize - 1) * 10_000;
                let end = start + 10_000;
                format!("{:.0}-{:.0} kB", start as f64 / 1000.0, end as f64 / 1000.0)
            }
            51..=60 => {
                let start = 500_000 + (category as usize - 51) * 50_000;
                let end = start + 50_000;
                format!("{:.0}-{:.0} kB", start as f64 / 1000.0, end as f64 / 1000.0)
            }
            61..=69 => {
                let start = 1_000_000 + (category as usize - 61) * 1_000_000;
                let end = start + 1_000_000;
                format!(
                    "{:.1}-{:.1} MB",
                    start as f64 / 1_000_000.0,
                    end as f64 / 1_000_000.0
                )
            }
            70 => "> 10 MB".to_string(),
            _ => format!("Unknown ({})", category),
        }
    }

    fn gc_category_to_midpoint(&self, category: u16) -> f64 {
        match category {
            0..=5 => (category * 5) as f64 + 2.5,
            6..=45 => 30.0 + (category - 6) as f64 + 0.5,
            46..=51 => 70.0 + (category - 46) as f64 * 5.0 + 2.5,
            _ => 50.0, // Default
        }
    }

    fn length_category_to_midpoint(&self, category: u16) -> f64 {
        match category {
            0 => 500.0,
            1..=50 => 1_000.0 + (category as f64 - 1.0) * 10_000.0 + 5_000.0,
            51..=60 => 500_000.0 + (category as f64 - 51.0) * 50_000.0 + 25_000.0,
            61..=69 => 1_000_000.0 + (category as f64 - 61.0) * 1_000_000.0 + 500_000.0,
            70 => 15_000_000.0, // Estimate for > 10MB
            _ => 100_000.0,     // Default
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
                .open(database_path)
                .with_context(|| format!("Failed to open database: {:?}", database_path))?
        };

        let rtxn = env.read_txn()?;
        let hash_db: Database<U64<BigEndian>, U64<BigEndian>> = env
            .open_database(&rtxn, Some("HASHES"))?
            .context("Hash database not found")?;

        // Open metadata database
        let metadata_path = database_path.with_extension("metadata.lmdb");
        let metadata_env = unsafe { heed::EnvOpenOptions::new().open(&metadata_path)? };
        let metadata_rtxn = metadata_env.read_txn()?;
        let metadata_db: Database<U32<BigEndian>, Bytes> = metadata_env
            .open_database(&metadata_rtxn, Some("METADATA"))?
            .context("Metadata database not found")?;

        let mut stats = DatabaseStats {
            total_hashes: 0,
            unique_files: 0,
            gc_distribution: HashMap::new(),
            length_distribution: HashMap::new(),
            file_metadata: Vec::new(),
        };

        // Count hashes and analyze metadata
        for item in hash_db.iter(&rtxn)? {
            let (_hash, packed_metadata) = item?;
            let metadata = HashMetadata::unpack(packed_metadata);

            stats.total_hashes += 1;

            // Count GC and length distributions
            *stats
                .gc_distribution
                .entry(metadata.gc_category)
                .or_insert(0) += 1;
            *stats
                .length_distribution
                .entry(metadata.length_category)
                .or_insert(0) += 1;
        }

        // Load file metadata
        for item in metadata_db.iter(&metadata_rtxn)? {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gc_category_ranges() {
        let stats = DatabaseStats {
            total_hashes: 0,
            unique_files: 0,
            gc_distribution: HashMap::new(),
            length_distribution: HashMap::new(),
            file_metadata: Vec::new(),
        };

        assert_eq!(stats.gc_category_to_range(0), "0-5%");
        assert_eq!(stats.gc_category_to_range(6), "30-31%");
        assert_eq!(stats.gc_category_to_range(26), "50-51%");
        assert_eq!(stats.gc_category_to_range(46), "70-75%");
    }

    #[test]
    fn test_length_category_ranges() {
        let stats = DatabaseStats {
            total_hashes: 0,
            unique_files: 0,
            gc_distribution: HashMap::new(),
            length_distribution: HashMap::new(),
            file_metadata: Vec::new(),
        };

        assert_eq!(stats.length_category_to_range(0), "< 1 kB");
        assert_eq!(stats.length_category_to_range(1), "1.0-11.0 kB");
        assert_eq!(stats.length_category_to_range(25), "240.0-250.0 kB");
        assert_eq!(stats.length_category_to_range(70), "> 10 MB");
    }
}
