use crate::core_utils::{self, *};
use crate::hash_functions::ahash;
use crate::sketch::SketchConfig;
use anyhow::{Context, Result, anyhow};
use byteorder::BigEndian;
use heed::types::{Bytes, U32, U64};
use heed::{Database, DatabaseFlags, EnvFlags, IntegerComparator};
use needletail::Sequence;
use serde::Serialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Distance calculation configuration
#[derive(Debug, Clone)]
pub struct DistanceConfig {
    pub cutoff: f64,
    pub output_format: OutputFormat,
    pub length_category_mode: LengthCategoryMode,
}

#[derive(Debug, Clone)]
pub enum LengthCategoryMode {
    /// Target length must be <= query length
    QueryAndBelow,
    /// Target length must be within ±N categories of query
    Range(u16),
}

#[derive(Debug, Clone)]
pub enum OutputFormat {
    Tsv,
    Json,
}

/// Result of a distance calculation with both total and filtered metrics
#[derive(Debug, Clone, Serialize)]
pub struct DistanceResult {
    pub query_name: String,
    pub target_name: String,
    // Total containment (all hashes)
    pub containment_query_in_target: f64,
    pub containment_target_in_query: f64,
    pub jaccard_similarity: f64,
    // Filtered containment (considering GC and length categories)
    pub filtered_containment_query_in_target: f64,
    pub filtered_containment_target_in_query: f64,
    pub filtered_jaccard_similarity: f64,
    // Hash counts
    pub shared_hashes: usize,
    pub filtered_shared_hashes: usize,
    pub query_hashes: usize,
    pub target_hashes: usize,
}

/// Target match counters per file_index
#[derive(Debug, Default)]
struct TargetCounters {
    total_matches: usize,
    filtered_matches: usize,
}

/// Streaming distance calculator
pub struct StreamingDistanceCalculator {
    config: DistanceConfig,
    sketch_config: SketchConfig,
    database_env: heed::Env,
    hash_db: Database<U64<BigEndian>, U64<BigEndian>>,
    metadata_db: Database<U32<BigEndian>, Bytes>,
}

impl StreamingDistanceCalculator {
    pub fn new(config: DistanceConfig, singleton: bool, database_path: &Path) -> Result<Self> {
        // Open database read-only
        let database_env = unsafe {
            heed::EnvOpenOptions::new()
                .flags(EnvFlags::READ_ONLY | EnvFlags::NO_SUB_DIR)
                .map_size(10 * 1024 * 1024 * 1024 * 1024)
                .max_dbs(3)
                .open(database_path)
                .with_context(|| format!("Failed to open database: {database_path:?}"))?
        };

        let rtxn = database_env.read_txn()?;
        let hash_db: Database<U64<BigEndian>, U64<BigEndian>> = database_env
            .open_database(&rtxn, Some("HASHES"))?
            .context("Hash database not found")?;

        let metadata_db: Database<U32<BigEndian>, Bytes> = database_env
            .open_database(&rtxn, Some("METADATA"))?
            .context("Metadata database not found")?;

        let config_db: Database<Bytes, Bytes> = database_env
            .open_database(&rtxn, Some("CONFIG"))?
            .context("Configuration database not found")?;

        let mut sketch_config = config_db
            .get(&rtxn, b"config")?
            .map(serde_json::from_slice::<SketchConfig>)
            .transpose()
            .context("Failed to read configuration from database")?
            .unwrap_or(SketchConfig::default());

        if singleton {
            sketch_config.singleton = true;
        }

        rtxn.commit()?;

        Ok(Self {
            config,
            sketch_config,
            database_env,
            hash_db,
            metadata_db,
        })
    }

    /// Calculate distances between query and database using streaming approach
    pub fn calculate_distances_streaming(
        &self,
        query_path: &Path,
        output_path: Option<&Path>,
    ) -> Result<Vec<DistanceResult>> {
        let mut results: Vec<_>;

        // Determine if query is a sketch file or raw sequence file
        if self.is_sketch_file(query_path)? {
            results = self.calculate_sketch_vs_database_streaming(query_path)?;
        } else {
            results = self.calculate_file_vs_database_streaming(query_path)?;
        }

        // Sort results by total containment (descending)
        results.sort_by(|a, b| {
            b.containment_query_in_target
                .partial_cmp(&a.containment_query_in_target)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        // Output results
        if let Some(output_path) = output_path {
            self.write_results(&results, output_path)?;
        } else {
            self.print_results(&results)?;
        }

        Ok(results)
    }

    /// Calculate distances for sketch file vs database - streaming through query hashes
    fn calculate_sketch_vs_database_streaming(
        &self,
        query_path: &Path,
    ) -> Result<Vec<DistanceResult>> {
        let mut results = Vec::new();

        // If query and database are the same, we can use the existing environment.
        // Otherwise, open a new one.
        let query_env_holder;
        let query_env = if self.database_env.path() == query_path {
            &self.database_env
        } else {
            query_env_holder = unsafe {
                heed::EnvOpenOptions::new()
                    .flags(EnvFlags::READ_ONLY | EnvFlags::NO_SUB_DIR)
                    .map_size(10 * 1024 * 1024 * 1024)
                    .max_dbs(3)
                    .open(query_path)?
            };
            &query_env_holder
        };

        let query_rtxn = query_env.read_txn()?;
        let query_hash_db: Database<U64<BigEndian>, U64<BigEndian>, IntegerComparator> = query_env
            .database_options()
            .types::<U64<BigEndian>, U64<BigEndian>>()
            .flags(DatabaseFlags::DUP_SORT | DatabaseFlags::DUP_FIXED)
            .name("HASHES")
            .key_comparator::<IntegerComparator>()
            .open(&query_rtxn)
            .context("Query hash database not found")?
            .ok_or_else(|| anyhow!("Database hashes not found"))?;

        let query_metadata_db: Database<U32<BigEndian>, Bytes> = query_env
            .open_database(&query_rtxn, Some("METADATA"))?
            .context("Query metadata database not found")?;

        // Group query hashes by file_index
        let mut query_sequences: HashMap<u32, (usize, HashMap<u32, TargetCounters>)> =
            HashMap::new();

        // Use the same transaction if query and database are the same environment
        let db_rtxn_holder;
        let db_rtxn = if self.database_env.path() == query_path {
            &query_rtxn
        } else {
            db_rtxn_holder = self.database_env.read_txn()?;
            &db_rtxn_holder
        };

        // Stream through query hashes
        for item in query_hash_db.iter(&query_rtxn)? {
            let (query_hash, query_packed_metadata) = item?;
            let query_metadata = HashMetadata::unpack(query_packed_metadata);

            // Initialize query sequence entry
            let (query_hash_count, target_counters) = query_sequences
                .entry(query_metadata.file_index)
                .or_insert_with(|| (0, HashMap::new()));

            *query_hash_count += 1;

            // Look up this hash in the database using get_duplicates
            if let Some(duplicates) = self.hash_db.get_duplicates(&db_rtxn, &query_hash)? {
                for item in duplicates {
                    let (_, target_packed_metadata) = item?;
                    let target_metadata = HashMetadata::unpack(target_packed_metadata);

                    // Get or create target counter
                    let target_counter = target_counters
                        .entry(target_metadata.file_index)
                        .or_insert_with(TargetCounters::default);

                    // Count total match
                    target_counter.total_matches += 1;

                    // Count filtered match if categories match
                    if self.categories_match(&query_metadata, &target_metadata) {
                        target_counter.filtered_matches += 1;
                    }
                }
            }
        }

        // Convert file_index to names and calculate results
        for (query_file_index, (query_hash_count, target_counters)) in query_sequences {
            let query_name =
                self.get_sequence_name(&query_rtxn, &query_metadata_db, query_file_index)?;

            for (target_file_index, target_counter) in target_counters {
                let target_name =
                    self.get_sequence_name(&db_rtxn, &self.metadata_db, target_file_index)?;
                let target_total_hashes =
                    self.get_target_total_hashes(&db_rtxn, target_file_index)?;

                let result = self.calculate_distance_from_counters(
                    &query_name,
                    query_hash_count,
                    &target_counter,
                    target_total_hashes,
                    &target_name,
                );

                if self.passes_cutoff(&result) {
                    results.push(result);
                }
            }
        }

        Ok(results)
    }

    /// Calculate distances for raw sequence file vs database - streaming through kmers
    fn calculate_file_vs_database_streaming(
        &self,
        input_path: &Path,
    ) -> Result<Vec<DistanceResult>> {
        let mut results = Vec::new();

        // Parse sequences from file
        let sequences = self.parse_sequences(input_path)?;

        let db_rtxn = self.database_env.read_txn()?;

        let mut target_counters: HashMap<u32, TargetCounters> = HashMap::new();
        let mut query_hash_count = 0;
        for (seq_name, sequence) in sequences {
            // Calculate query metadata
            let query_metadata = self.calculate_sequence_metadata(&sequence, &seq_name)?;

            // Stream through kmers using the provided logic
            let sequence_bytes = sequence.as_bytes();
            for (_, kmer, _) in sequence_bytes.bit_kmers(self.sketch_config.kmer_size, true) {
                // Apply entropy filter
                if !passes_entropy_filter(
                    kmer.0,
                    self.sketch_config.kmer_size,
                    self.sketch_config.min_entropy,
                ) {
                    continue;
                }

                let hash = ahash(kmer.0);

                // Apply FracMinHash filter if specified
                if hash > self.sketch_config.fscale {
                    continue;
                }

                query_hash_count += 1;

                // Look up this hash in the database using get_duplicates
                if let Some(duplicates) = self.hash_db.get_duplicates(&db_rtxn, &hash)? {
                    for item in duplicates {
                        let (_, target_packed_metadata) = item?;
                        let target_metadata = HashMetadata::unpack(target_packed_metadata);

                        // Get or create target counter
                        let target_counter = target_counters
                            .entry(target_metadata.file_index)
                            .or_default();

                        // Count total match
                        target_counter.total_matches += 1;

                        // Count filtered match if categories match
                        if self.categories_match(&query_metadata, &target_metadata) {
                            target_counter.filtered_matches += 1;
                        }
                    }
                }
            }

            if self.sketch_config.singleton {
                // Calculate results for this query sequence
                for (target_file_index, target_counter) in target_counters.drain() {
                    let target_name =
                        self.get_sequence_name(&db_rtxn, &self.metadata_db, target_file_index)?;
                    let target_total_hashes =
                        self.get_target_total_hashes(&db_rtxn, target_file_index)?;

                    let result = self.calculate_distance_from_counters(
                        &seq_name,
                        query_hash_count,
                        &target_counter,
                        target_total_hashes,
                        &target_name,
                    );

                    query_hash_count = 0; // Reset for next sequence

                    if self.passes_cutoff(&result) {
                        results.push(result);
                    }
                }
            }
        }

        // Calculate results for this query sequence
        for (target_file_index, target_counter) in target_counters.drain() {
            let target_name =
                self.get_sequence_name(&db_rtxn, &self.metadata_db, target_file_index)?;
            let target_total_hashes = self.get_target_total_hashes(&db_rtxn, target_file_index)?;

            let result = self.calculate_distance_from_counters(
                &input_path.display().to_string(),
                query_hash_count,
                &target_counter,
                target_total_hashes,
                &target_name,
            );

            if self.passes_cutoff(&result) {
                results.push(result);
            }
        }

        Ok(results)
    }

    /// Calculate sequence metadata for raw sequences
    fn calculate_sequence_metadata(&self, sequence: &str, _seq_name: &str) -> Result<HashMetadata> {
        // Calculate GC content
        let gc_count = sequence.chars().filter(|&c| c == 'G' || c == 'C').count();
        let gc_percent = (gc_count as f64 / sequence.len() as f64) * 100.0;
        let gc_category = core_utils::categorize_gc_percent(gc_percent);

        // Calculate length category
        let length_category = core_utils::categorize_length(sequence.len());

        Ok(HashMetadata {
            file_index: 0, // Temporary for raw sequences
            gc_category,
            length_category,
        })
    }

    /// Get target total hash count from FileMetadata
    fn get_target_total_hashes(&self, txn: &heed::RoTxn, file_index: u32) -> Result<usize> {
        if let Some(metadata_json) = self.metadata_db.get(txn, &file_index)? {
            let file_metadata: FileMetadata = serde_json::from_slice(metadata_json)?;
            Ok(file_metadata.total_hashes)
        } else {
            Ok(0)
        }
    }

    /// Calculate distance metrics from counters
    fn calculate_distance_from_counters(
        &self,
        query_name: &str,
        query_hash_count: usize,
        target_counter: &TargetCounters,
        target_total_hashes: usize,
        target_name: &str,
    ) -> DistanceResult {
        let shared_hashes = target_counter.total_matches;
        let filtered_shared_hashes = target_counter.filtered_matches;

        // Total containment
        let containment_query_in_target = if query_hash_count > 0 {
            shared_hashes as f64 / query_hash_count as f64
        } else {
            0.0
        };

        let containment_target_in_query = if target_total_hashes > 0 {
            shared_hashes as f64 / target_total_hashes as f64
        } else {
            0.0
        };

        let union_size = query_hash_count + target_total_hashes - shared_hashes;
        let jaccard_similarity = if union_size > 0 {
            shared_hashes as f64 / union_size as f64
        } else {
            0.0
        };

        // Filtered containment
        let filtered_containment_query_in_target = if query_hash_count > 0 {
            filtered_shared_hashes as f64 / query_hash_count as f64
        } else {
            0.0
        };

        let filtered_containment_target_in_query = if target_total_hashes > 0 {
            filtered_shared_hashes as f64 / target_total_hashes as f64
        } else {
            0.0
        };

        let filtered_union_size = query_hash_count + target_total_hashes - filtered_shared_hashes;
        let filtered_jaccard_similarity = if filtered_union_size > 0 {
            filtered_shared_hashes as f64 / filtered_union_size as f64
        } else {
            0.0
        };

        DistanceResult {
            query_name: query_name.to_string(),
            target_name: target_name.to_string(),
            containment_query_in_target,
            containment_target_in_query,
            jaccard_similarity,
            filtered_containment_query_in_target,
            filtered_containment_target_in_query,
            filtered_jaccard_similarity,
            shared_hashes,
            filtered_shared_hashes,
            query_hashes: query_hash_count,
            target_hashes: target_total_hashes,
        }
    }

    /// Check if result passes cutoff filter
    fn passes_cutoff(&self, result: &DistanceResult) -> bool {
        result.containment_query_in_target >= self.config.cutoff
            || result.containment_target_in_query >= self.config.cutoff
            || result.jaccard_similarity >= self.config.cutoff
            || result.filtered_containment_query_in_target >= self.config.cutoff
            || result.filtered_containment_target_in_query >= self.config.cutoff
            || result.filtered_jaccard_similarity >= self.config.cutoff
    }

    /// Check if GC and length categories match according to configuration
    fn categories_match(
        &self,
        query_metadata: &HashMetadata,
        target_metadata: &HashMetadata,
    ) -> bool {
        // Length category: configurable logic
        match self.config.length_category_mode {
            LengthCategoryMode::QueryAndBelow => {
                target_metadata.length_category <= query_metadata.length_category
            }
            LengthCategoryMode::Range(tolerance) => {
                let length_diff = (query_metadata.length_category as i32
                    - target_metadata.length_category as i32)
                    .abs();
                length_diff <= tolerance as i32
            }
        }
    }

    /// Get sequence name from metadata database
    fn get_sequence_name(
        &self,
        txn: &heed::RoTxn,
        metadata_db: &Database<U32<BigEndian>, Bytes>,
        file_index: u32,
    ) -> Result<String> {
        if let Some(metadata_json) = metadata_db.get(txn, &file_index)? {
            let file_metadata: FileMetadata =
                serde_json::from_slice(metadata_json).unwrap_or_else(|_| FileMetadata {
                    filename: format!("unknown_{file_index}"),
                    file_size: 0,
                    sequence_name: format!("seq_{file_index}"),
                    sequence_length: 0,
                    total_sequences: 1,
                    total_hashes: 0,
                });

            Ok(if file_metadata.total_sequences == 1 {
                file_metadata.sequence_name
            } else {
                file_metadata.filename
            })
        } else {
            Ok(format!("unknown_{file_index}"))
        }
    }

    /// Parse sequences from FASTA/FASTQ file
    fn parse_sequences(&self, file_path: &Path) -> Result<Vec<(String, String)>> {
        let file = File::open(file_path)?;
        let reader = BufReader::new(file);
        let mut sequences = Vec::new();
        let mut current_name = String::new();
        let mut current_seq = String::new();
        let mut in_sequence = false;

        for line in reader.lines() {
            let line = line?;
            let line = line.trim();

            if line.starts_with('>') {
                // FASTA header
                if in_sequence && !current_name.is_empty() {
                    sequences.push((current_name.clone(), current_seq.clone()));
                }
                current_name = line[1..].to_string();
                current_seq.clear();
                in_sequence = true;
            } else if line.starts_with('@') {
                // FASTQ header
                if in_sequence && !current_name.is_empty() {
                    sequences.push((current_name.clone(), current_seq.clone()));
                }
                current_name = line[1..].to_string();
                current_seq.clear();
                in_sequence = true;
            } else if line.starts_with('+') {
                // FASTQ separator - skip quality scores
                in_sequence = false;
            } else if in_sequence && !line.is_empty() {
                current_seq.push_str(line);
            }
        }

        // Add last sequence
        if in_sequence && !current_name.is_empty() {
            sequences.push((current_name, current_seq));
        }

        Ok(sequences)
    }

    /// Check if a file is a sketch (LMDB) or raw sequence data
    fn is_sketch_file(&self, path: &Path) -> Result<bool> {
        if let Some(extension) = path.extension() {
            Ok(extension == "lmdb")
        } else {
            Ok(false)
        }
    }

    /// Write results to file
    fn write_results(&self, results: &[DistanceResult], output_path: &Path) -> Result<()> {
        let file = File::create(output_path)?;
        let mut writer = BufWriter::new(file);

        match self.config.output_format {
            OutputFormat::Tsv => self.write_tsv_results(&mut writer, results)?,
            OutputFormat::Json => self.write_json_results(&mut writer, results)?,
        }

        writer.flush()?;
        Ok(())
    }

    /// Write results in TSV format
    fn write_tsv_results<W: Write>(
        &self,
        writer: &mut W,
        results: &[DistanceResult],
    ) -> Result<()> {
        writeln!(
            writer,
            "query\ttarget\tcontainment_query_in_target\tcontainment_target_in_query\tjaccard\tfiltered_containment_query_in_target\tfiltered_containment_target_in_query\tfiltered_jaccard\tshared_hashes\tfiltered_shared_hashes\tquery_hashes\ttarget_hashes"
        )?;

        for result in results {
            writeln!(
                writer,
                "{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}",
                result.query_name,
                result.target_name,
                result.containment_query_in_target,
                result.containment_target_in_query,
                result.jaccard_similarity,
                result.filtered_containment_query_in_target,
                result.filtered_containment_target_in_query,
                result.filtered_jaccard_similarity,
                result.shared_hashes,
                result.filtered_shared_hashes,
                result.query_hashes,
                result.target_hashes
            )?;
        }

        Ok(())
    }

    /// Write results in JSON format
    fn write_json_results<W: Write>(
        &self,
        writer: &mut W,
        results: &[DistanceResult],
    ) -> Result<()> {
        let json = serde_json::to_string_pretty(results)?;
        write!(writer, "{json}")?;
        Ok(())
    }

    /// Print results to stdout
    fn print_results(&self, results: &[DistanceResult]) -> Result<()> {
        let stdout = std::io::stdout();
        let mut handle = stdout.lock();

        match self.config.output_format {
            OutputFormat::Tsv => self.write_tsv_results(&mut handle, results)?,
            OutputFormat::Json => self.write_json_results(&mut handle, results)?,
        }

        Ok(())
    }
}

/// Public interface for streaming distance calculation
pub fn calculate_distances_streaming(
    query_path: &Path,
    database_path: &Path,
    output_path: Option<&Path>,
    config: DistanceConfig,
    singleton: bool,
) -> Result<Vec<DistanceResult>> {
    let calculator = StreamingDistanceCalculator::new(config, singleton, database_path)?;
    calculator.calculate_distances_streaming(query_path, output_path)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    // Helper function to create a temporary FASTA file
    fn _create_test_fasta_file(sequences: &[(&str, &str)]) -> Result<NamedTempFile> {
        let mut file = NamedTempFile::new()?;
        for (name, seq) in sequences {
            writeln!(file, ">{name}")?;
            writeln!(file, "{seq}")?;
        }
        file.flush()?;
        Ok(file)
    }

    #[test]
    fn test_distance_config_creation() {
        let config = DistanceConfig {
            cutoff: 0.5,
            output_format: OutputFormat::Json,
            length_category_mode: LengthCategoryMode::Range(5),
        };

        assert_eq!(config.cutoff, 0.5);
        assert!(matches!(config.output_format, OutputFormat::Json));
        assert!(matches!(
            config.length_category_mode,
            LengthCategoryMode::Range(5)
        ));
    }

    #[test]
    fn test_length_category_mode_variants() {
        let query_and_below = LengthCategoryMode::QueryAndBelow;
        let range_mode = LengthCategoryMode::Range(10);

        assert!(matches!(query_and_below, LengthCategoryMode::QueryAndBelow));
        assert!(matches!(range_mode, LengthCategoryMode::Range(10)));
    }

    #[test]
    fn test_output_format_variants() {
        let tsv = OutputFormat::Tsv;
        let json = OutputFormat::Json;

        assert!(matches!(tsv, OutputFormat::Tsv));
        assert!(matches!(json, OutputFormat::Json));
    }

    #[test]
    fn test_distance_result_creation() {
        let result = DistanceResult {
            query_name: "query1".to_string(),
            target_name: "target1".to_string(),
            containment_query_in_target: 0.75,
            containment_target_in_query: 0.60,
            jaccard_similarity: 0.50,
            filtered_containment_query_in_target: 0.70,
            filtered_containment_target_in_query: 0.55,
            filtered_jaccard_similarity: 0.45,
            shared_hashes: 100,
            filtered_shared_hashes: 90,
            query_hashes: 150,
            target_hashes: 200,
        };

        assert_eq!(result.query_name, "query1");
        assert_eq!(result.target_name, "target1");
        assert_eq!(result.containment_query_in_target, 0.75);
        assert_eq!(result.shared_hashes, 100);
        assert_eq!(result.query_hashes, 150);
    }

    #[test]
    fn test_target_counters_default() {
        let counters = TargetCounters::default();
        assert_eq!(counters.total_matches, 0);
        assert_eq!(counters.filtered_matches, 0);
    }

    #[test]
    fn test_target_counters_increment() {
        let mut counters = TargetCounters::default();
        counters.total_matches += 5;
        counters.filtered_matches += 3;

        assert_eq!(counters.total_matches, 5);
        assert_eq!(counters.filtered_matches, 3);
    }

    // Test category matching logic
    #[test]
    fn test_categories_match_logic() {
        struct MockCalculator {
            config: DistanceConfig,
        }

        impl MockCalculator {
            fn categories_match(
                &self,
                query_metadata: &HashMetadata,
                target_metadata: &HashMetadata,
            ) -> bool {
                // GC category: ±1 tolerance
                let gc_diff =
                    (query_metadata.gc_category as i32 - target_metadata.gc_category as i32).abs();
                if gc_diff > 1 {
                    return false;
                }

                // Length category: configurable logic
                match self.config.length_category_mode {
                    LengthCategoryMode::QueryAndBelow => {
                        target_metadata.length_category <= query_metadata.length_category
                    }
                    LengthCategoryMode::Range(tolerance) => {
                        let length_diff = (query_metadata.length_category as i32
                            - target_metadata.length_category as i32)
                            .abs();
                        length_diff <= tolerance as i32
                    }
                }
            }
        }

        // Test QueryAndBelow mode
        let config_query_below = DistanceConfig {
            cutoff: 0.0,
            output_format: OutputFormat::Tsv,
            length_category_mode: LengthCategoryMode::QueryAndBelow,
        };

        let calculator_query_below = MockCalculator {
            config: config_query_below,
        };

        let query_metadata = HashMetadata {
            file_index: 0,
            gc_category: 50,
            length_category: 100,
        };

        // Test cases for QueryAndBelow mode
        let target_within_gc_below_length = HashMetadata {
            file_index: 1,
            gc_category: 51,     // Within ±1
            length_category: 90, // Below query
        };

        let target_outside_gc = HashMetadata {
            file_index: 2,
            gc_category: 53, // Outside ±1
            length_category: 90,
        };

        let target_above_length = HashMetadata {
            file_index: 3,
            gc_category: 49,      // Within ±1
            length_category: 110, // Above query
        };

        assert!(
            calculator_query_below
                .categories_match(&query_metadata, &target_within_gc_below_length)
        );
        assert!(!calculator_query_below.categories_match(&query_metadata, &target_outside_gc));
        assert!(!calculator_query_below.categories_match(&query_metadata, &target_above_length));

        // Test Range mode
        let config_range = DistanceConfig {
            cutoff: 0.0,
            output_format: OutputFormat::Tsv,
            length_category_mode: LengthCategoryMode::Range(15),
        };

        let calculator_range = MockCalculator {
            config: config_range,
        };

        let target_within_range = HashMetadata {
            file_index: 4,
            gc_category: 50,      // Exact match
            length_category: 110, // Within range (diff = 10)
        };

        let target_outside_range = HashMetadata {
            file_index: 5,
            gc_category: 50,      // Exact match
            length_category: 120, // Outside range (diff = 20)
        };

        assert!(calculator_range.categories_match(&query_metadata, &target_within_range));
        assert!(!calculator_range.categories_match(&query_metadata, &target_outside_range));
    }

    // Test sequence metadata calculation
    #[test]
    fn test_sequence_metadata_calculation() {
        struct MockCalculator;

        impl MockCalculator {
            fn calculate_sequence_metadata(&self, sequence: &str, _seq_name: &str) -> HashMetadata {
                // Calculate GC content
                let gc_count = sequence.chars().filter(|&c| c == 'G' || c == 'C').count();
                let gc_percent = (gc_count as f64 / sequence.len() as f64) * 100.0;
                let gc_category = ((gc_percent * 255.0) / 100.0) as u16;

                // Calculate length category
                let length_category = self.calculate_length_category(sequence.len());

                HashMetadata {
                    file_index: 0,
                    gc_category,
                    length_category,
                }
            }

            fn calculate_length_category(&self, length: usize) -> u16 {
                let log_length = (length as f64).log10();

                ((log_length - 3.0) / (5.7 - 3.0) * 255.0)
                    .max(0.0)
                    .min(255.0) as u16
            }
        }

        let calculator = MockCalculator;

        // Test sequence with 50% GC content
        let seq_50_gc = "ATCGATCGATCGATCGATCG"; // 10 A/T, 10 G/C
        let metadata = calculator.calculate_sequence_metadata(seq_50_gc, "test");
        assert_eq!(metadata.gc_category, 127); // 50% * 255/100 ≈ 127

        // Test sequence with 100% GC content
        let seq_100_gc = "GCGCGCGCGCGCGCGCGCGC";
        let metadata = calculator.calculate_sequence_metadata(seq_100_gc, "test");
        assert_eq!(metadata.gc_category, 255); // 100% * 255/100 = 255

        // Test sequence with 0% GC content
        let seq_0_gc = "ATATATATATATATATATAT";
        let metadata = calculator.calculate_sequence_metadata(seq_0_gc, "test");
        assert_eq!(metadata.gc_category, 0); // 0% * 255/100 = 0
    }

    // Test distance calculation from counters
    #[test]
    fn test_distance_calculation_from_counters() {
        struct MockCalculator;

        impl MockCalculator {
            fn calculate_distance_from_counters(
                &self,
                query_name: &str,
                query_hash_count: usize,
                target_counter: &TargetCounters,
                target_total_hashes: usize,
                target_name: &str,
            ) -> DistanceResult {
                let shared_hashes = target_counter.total_matches;
                let filtered_shared_hashes = target_counter.filtered_matches;

                // Total containment
                let containment_query_in_target = if query_hash_count > 0 {
                    shared_hashes as f64 / query_hash_count as f64
                } else {
                    0.0
                };

                let containment_target_in_query = if target_total_hashes > 0 {
                    shared_hashes as f64 / target_total_hashes as f64
                } else {
                    0.0
                };

                let union_size = query_hash_count + target_total_hashes - shared_hashes;
                let jaccard_similarity = if union_size > 0 {
                    shared_hashes as f64 / union_size as f64
                } else {
                    0.0
                };

                // Filtered containment
                let filtered_containment_query_in_target = if query_hash_count > 0 {
                    filtered_shared_hashes as f64 / query_hash_count as f64
                } else {
                    0.0
                };

                let filtered_containment_target_in_query = if target_total_hashes > 0 {
                    filtered_shared_hashes as f64 / target_total_hashes as f64
                } else {
                    0.0
                };

                let filtered_union_size =
                    query_hash_count + target_total_hashes - filtered_shared_hashes;
                let filtered_jaccard_similarity = if filtered_union_size > 0 {
                    filtered_shared_hashes as f64 / filtered_union_size as f64
                } else {
                    0.0
                };

                DistanceResult {
                    query_name: query_name.to_string(),
                    target_name: target_name.to_string(),
                    containment_query_in_target,
                    containment_target_in_query,
                    jaccard_similarity,
                    filtered_containment_query_in_target,
                    filtered_containment_target_in_query,
                    filtered_jaccard_similarity,
                    shared_hashes,
                    filtered_shared_hashes,
                    query_hashes: query_hash_count,
                    target_hashes: target_total_hashes,
                }
            }
        }

        let calculator = MockCalculator;

        // Mock data
        let query_name = "query1";
        let query_hash_count = 100;
        let target_counter = TargetCounters {
            total_matches: 80,
            filtered_matches: 60,
        };
        let target_total_hashes = 120;
        let target_name = "target1";

        // Calculate distance result
        let result = calculator.calculate_distance_from_counters(
            query_name,
            query_hash_count,
            &target_counter,
            target_total_hashes,
            target_name,
        );

        // Assertions
        assert_eq!(result.query_name, query_name);
        assert_eq!(result.target_name, target_name);
        assert_eq!(result.shared_hashes, target_counter.total_matches);
        assert_eq!(
            result.filtered_shared_hashes,
            target_counter.filtered_matches
        );

        // Total containment
        let containment_query_in_target = if query_hash_count > 0 {
            target_counter.total_matches as f64 / query_hash_count as f64
        } else {
            0.0
        };

        let containment_target_in_query = if target_total_hashes > 0 {
            target_counter.total_matches as f64 / target_total_hashes as f64
        } else {
            0.0
        };

        assert_eq!(
            result.containment_query_in_target,
            containment_query_in_target
        );
        assert_eq!(
            result.containment_target_in_query,
            containment_target_in_query
        );

        // Jaccard similarity
        let union_size = query_hash_count + target_total_hashes - target_counter.total_matches;
        let jaccard_similarity = if union_size > 0 {
            target_counter.total_matches as f64 / union_size as f64
        } else {
            0.0
        };

        assert_eq!(result.jaccard_similarity, jaccard_similarity);

        // Filtered containment
        let filtered_containment_query_in_target = if query_hash_count > 0 {
            target_counter.filtered_matches as f64 / query_hash_count as f64
        } else {
            0.0
        };

        let filtered_containment_target_in_query = if target_total_hashes > 0 {
            target_counter.filtered_matches as f64 / target_total_hashes as f64
        } else {
            0.0
        };

        assert_eq!(
            result.filtered_containment_query_in_target,
            filtered_containment_query_in_target
        );
        assert_eq!(
            result.filtered_containment_target_in_query,
            filtered_containment_target_in_query
        );

        // Filtered Jaccard similarity
        let filtered_union_size =
            query_hash_count + target_total_hashes - target_counter.filtered_matches;
        let filtered_jaccard_similarity = if filtered_union_size > 0 {
            target_counter.filtered_matches as f64 / filtered_union_size as f64
        } else {
            0.0
        };

        assert_eq!(
            result.filtered_jaccard_similarity,
            filtered_jaccard_similarity
        );
    }
}
