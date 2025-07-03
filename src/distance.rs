use crate::core_utils::*;
use crate::sketch::{SketchConfig, sketch_files};
use anyhow::{Context, Result};
use byteorder::BigEndian;
use heed::Database;
use heed::types::{Bytes, U32, U64};
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Distance calculation configuration
#[derive(Debug, Clone)]
pub struct DistanceConfig {
    pub cutoff: f64,
    pub output_format: OutputFormat,
}

#[derive(Debug, Clone)]
pub enum OutputFormat {
    Tsv,
    Json,
}

/// Result of a distance calculation
#[derive(Debug, Clone, Serialize)]
pub struct DistanceResult {
    pub query_name: String,
    pub target_name: String,
    pub containment_query_in_target: f64,
    pub containment_target_in_query: f64,
    pub jaccard_similarity: f64,
    pub shared_hashes: usize,
    pub query_hashes: usize,
    pub target_hashes: usize,
}

/// Distance calculator
pub struct DistanceCalculator {
    config: DistanceConfig,
}

impl DistanceCalculator {
    pub fn new(config: DistanceConfig) -> Self {
        Self { config }
    }

    /// Calculate distances between query and database
    pub fn calculate_distances(
        &self,
        query_path: &Path,
        database_path: &Path,
        output_path: Option<&Path>,
    ) -> Result<Vec<DistanceResult>> {
        // Load query sketch
        let query_sketch = self.load_sketch(query_path)?;

        // Load database sketches
        let database_sketches = self.load_database(database_path)?;

        // Calculate distances
        let mut results = Vec::new();

        for (query_name, query_hashes) in &query_sketch {
            for (target_name, target_hashes) in &database_sketches {
                let distance_result = self.calculate_pairwise_distance(
                    query_name,
                    query_hashes,
                    target_name,
                    target_hashes,
                );

                // Apply cutoff filter
                if distance_result.containment_query_in_target >= self.config.cutoff
                    || distance_result.containment_target_in_query >= self.config.cutoff
                    || distance_result.jaccard_similarity >= self.config.cutoff
                {
                    results.push(distance_result);
                }
            }
        }

        // Sort results by containment (descending)
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

    /// Calculate distance between a single file/raw sequence and database
    pub fn calculate_file_vs_database(
        &self,
        input_path: &Path,
        database_path: &Path,
        output_path: Option<&Path>,
        sketch_config: SketchConfig,
    ) -> Result<Vec<DistanceResult>> {
        // Check if input is already a sketch or needs to be sketched
        let query_sketch = if self.is_sketch_file(input_path)? {
            self.load_sketch(input_path)?
        } else {
            // Sketch the input file on-the-fly
            let temp_dir = tempfile::tempdir()?;
            let temp_sketch_path = temp_dir.path().join("temp_query.lmdb");

            sketch_files(
                &[input_path.to_path_buf()],
                temp_sketch_path.clone(),
                sketch_config,
            )?;
            self.load_sketch(&temp_sketch_path)?
        };

        // Load database and calculate distances
        let database_sketches = self.load_database(database_path)?;

        let mut results = Vec::new();

        for (query_name, query_hashes) in &query_sketch {
            for (target_name, target_hashes) in &database_sketches {
                let distance_result = self.calculate_pairwise_distance(
                    query_name,
                    query_hashes,
                    target_name,
                    target_hashes,
                );

                if distance_result.containment_query_in_target >= self.config.cutoff
                    || distance_result.containment_target_in_query >= self.config.cutoff
                    || distance_result.jaccard_similarity >= self.config.cutoff
                {
                    results.push(distance_result);
                }
            }
        }

        results.sort_by(|a, b| {
            b.containment_query_in_target
                .partial_cmp(&a.containment_query_in_target)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        if let Some(output_path) = output_path {
            self.write_results(&results, output_path)?;
        } else {
            self.print_results(&results)?;
        }

        Ok(results)
    }

    /// Calculate pairwise distance between two sets of hashes
    fn calculate_pairwise_distance(
        &self,
        query_name: &str,
        query_hashes: &HashSet<u64>,
        target_name: &str,
        target_hashes: &HashSet<u64>,
    ) -> DistanceResult {
        let intersection: HashSet<_> = query_hashes.intersection(target_hashes).collect();
        let shared_hashes = intersection.len();

        let query_size = query_hashes.len();
        let target_size = target_hashes.len();

        let containment_query_in_target = if query_size > 0 {
            shared_hashes as f64 / query_size as f64
        } else {
            0.0
        };

        let containment_target_in_query = if target_size > 0 {
            shared_hashes as f64 / target_size as f64
        } else {
            0.0
        };

        let union_size = query_size + target_size - shared_hashes;
        let jaccard_similarity = if union_size > 0 {
            shared_hashes as f64 / union_size as f64
        } else {
            0.0
        };

        DistanceResult {
            query_name: query_name.to_string(),
            target_name: target_name.to_string(),
            containment_query_in_target,
            containment_target_in_query,
            jaccard_similarity,
            shared_hashes,
            query_hashes: query_size,
            target_hashes: target_size,
        }
    }

    /// Load sketch from LMDB file
    fn load_sketch(&self, sketch_path: &Path) -> Result<HashMap<String, HashSet<u64>>> {
        let env = unsafe {
            heed::EnvOpenOptions::new()
                .open(sketch_path)
                .with_context(|| format!("Failed to open sketch database: {:?}", sketch_path))?
        };

        let rtxn = env.read_txn()?;
        let hash_db: Database<U64<BigEndian>, U64<BigEndian>> = env
            .open_database(&rtxn, Some("HASHES"))?
            .context("Hash database not found")?;

        // Load metadata to get sequence/file names
        let metadata_path = sketch_path.with_extension("metadata.lmdb");
        let metadata_env = unsafe { heed::EnvOpenOptions::new().open(&metadata_path)? };
        let metadata_rtxn = metadata_env.read_txn()?;
        let metadata_db: Database<U32<BigEndian>, Bytes> = metadata_env
            .open_database(&metadata_rtxn, Some("METADATA"))?
            .context("Metadata database not found")?;

        // Group hashes by file index
        let mut sketches: HashMap<String, HashSet<u64>> = HashMap::new();

        for item in hash_db.iter(&rtxn)? {
            let (hash, packed_metadata) = item?;
            let metadata = HashMetadata::unpack(packed_metadata);

            // Get file/sequence name from metadata
            let name = if let Ok(Some(metadata_json)) =
                metadata_db.get(&metadata_rtxn, &metadata.file_index)
            {
                let file_metadata: FileMetadata = serde_json::from_slice(metadata_json)
                    .unwrap_or_else(|_| FileMetadata {
                        filename: format!("unknown_{}", metadata.file_index),
                        file_size: 0,
                        sequence_name: format!("seq_{}", metadata.file_index),
                        sequence_length: 0,
                        total_sequences: 1,
                    });

                if file_metadata.total_sequences == 1 {
                    file_metadata.sequence_name
                } else {
                    file_metadata.filename
                }
            } else {
                format!("unknown_{}", metadata.file_index)
            };

            sketches
                .entry(name)
                .or_insert_with(HashSet::new)
                .insert(hash);
        }

        Ok(sketches)
    }

    /// Load database sketches
    fn load_database(&self, database_path: &Path) -> Result<HashMap<String, HashSet<u64>>> {
        // For now, assume database is a single LMDB file
        // Could be extended to handle multiple files
        self.load_sketch(database_path)
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

    /// Write results in TSV format (BLAST-like)
    fn write_tsv_results<W: Write>(
        &self,
        writer: &mut W,
        results: &[DistanceResult],
    ) -> Result<()> {
        // Write header
        writeln!(
            writer,
            "query\ttarget\tcontainment_query_in_target\tcontainment_target_in_query\tjaccard\tshared_hashes\tquery_hashes\ttarget_hashes"
        )?;

        // Write results
        for result in results {
            writeln!(
                writer,
                "{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}",
                result.query_name,
                result.target_name,
                result.containment_query_in_target,
                result.containment_target_in_query,
                result.jaccard_similarity,
                result.shared_hashes,
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
        write!(writer, "{}", json)?;
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

/// Public interface for distance calculation
pub fn calculate_distances(
    query_path: &Path,
    database_path: &Path,
    output_path: Option<&Path>,
    config: DistanceConfig,
) -> Result<Vec<DistanceResult>> {
    let calculator = DistanceCalculator::new(config);
    calculator.calculate_distances(query_path, database_path, output_path)
}

/// Calculate distances for input file vs database
pub fn calculate_file_vs_database(
    input_path: &Path,
    database_path: &Path,
    output_path: Option<&Path>,
    distance_config: DistanceConfig,
    sketch_config: SketchConfig,
) -> Result<Vec<DistanceResult>> {
    let calculator = DistanceCalculator::new(distance_config);
    calculator.calculate_file_vs_database(input_path, database_path, output_path, sketch_config)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_distance_calculation() {
        let mut query_hashes = HashSet::new();
        query_hashes.insert(1);
        query_hashes.insert(2);
        query_hashes.insert(3);

        let mut target_hashes = HashSet::new();
        target_hashes.insert(2);
        target_hashes.insert(3);
        target_hashes.insert(4);
        target_hashes.insert(5);

        let calculator = DistanceCalculator::new(DistanceConfig {
            cutoff: 0.0,
            output_format: OutputFormat::Tsv,
        });

        let result = calculator.calculate_pairwise_distance(
            "query",
            &query_hashes,
            "target",
            &target_hashes,
        );

        assert_eq!(result.shared_hashes, 2); // Elements 2 and 3
        assert!((result.containment_query_in_target - 2.0 / 3.0).abs() < 1e-6);
        assert!((result.containment_target_in_query - 2.0 / 4.0).abs() < 1e-6);
        assert!((result.jaccard_similarity - 2.0 / 5.0).abs() < 1e-6); // 2 shared / 5 union
    }
}
