use crate::bias::{BIAS_TABLE_SERIALIZED_SIZE, BiasTable};
use crate::format::{
    BUCKET_COUNT, BUCKET_TABLE_SIZE, BucketMeta, DATA_START, ENTRY_SIZE, Entry,
    FLAG_HAS_BIAS_TABLE, HEADER_SIZE, Header, MAGIC, VERSION,
};
use crate::io::{extract_unique_hashes, read_entries, write_entries};
use crate::sketch::{SketchConfig, SketchResult};
use bytemuck;
use memmap2::MmapMut;
use rayon::prelude::*;
use std::io;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Duration;
use xorf::{BinaryFuse8, DmaSerializable};

/// BinaryFuse8 descriptor size (seed: u64 + segment_length: u32 + segment_length_mask: u32 + segment_count_length: u32)
pub const FILTER_DESCRIPTOR_SIZE: usize = 20;

/// Serialize a BinaryFuse8 filter to bytes using DMA format.
/// Format: descriptor_size(u32) + fingerprints_size(u32) + descriptor + fingerprints
pub fn serialize_filter(filter: &BinaryFuse8) -> Vec<u8> {
    let fingerprints = filter.dma_fingerprints();
    let descriptor_size = FILTER_DESCRIPTOR_SIZE as u32;
    let fingerprints_size = fingerprints.len() as u32;

    let total_size = 4 + 4 + FILTER_DESCRIPTOR_SIZE + fingerprints.len();
    let mut out = vec![0u8; total_size];

    out[0..4].copy_from_slice(&descriptor_size.to_le_bytes());
    out[4..8].copy_from_slice(&fingerprints_size.to_le_bytes());
    filter.dma_copy_descriptor_to(&mut out[8..8 + FILTER_DESCRIPTOR_SIZE]);
    out[8 + FILTER_DESCRIPTOR_SIZE..].copy_from_slice(fingerprints);

    out
}

/// Serialize sample names in order. Format: [len0][bytes0][len1][bytes1]...
/// Warns if any sample name exceeds u16::MAX bytes and will be truncated.
fn serialize_sample_names(names: &[String]) -> Vec<u8> {
    let mut buf = Vec::new();
    for (idx, name) in names.iter().enumerate() {
        let bytes = name.as_bytes();
        let len = if bytes.len() > u16::MAX as usize {
            eprintln!(
                "Warning: Sample name at index {} is {} bytes, truncating to {} bytes",
                idx,
                bytes.len(),
                u16::MAX
            );
            u16::MAX
        } else {
            bytes.len() as u16
        };
        buf.extend_from_slice(&len.to_le_bytes());
        buf.extend_from_slice(&bytes[..len as usize]);
    }
    buf
}

/// Serialize sample sizes (hash counts) as u64 array.
fn serialize_sample_sizes(sizes: &[u64]) -> Vec<u8> {
    bytemuck::cast_slice(sizes).to_vec()
}

#[derive(Clone)]
pub struct CompactConfig {
    pub num_threads: usize,
}

impl Default for CompactConfig {
    fn default() -> Self {
        Self { num_threads: 1 }
    }
}

struct ProcessedBucket {
    bucket_id: usize,
    entry_count: u64,
    unique_hash_count: u64,
    filter_bytes: Vec<u8>,
    sample_hash_counts: std::collections::HashMap<u32, u64>,
}

#[derive(Debug, Clone)]
pub struct IndexStats {
    pub total_entries: u64,
    pub unique_hashes: u64,
    pub sample_count: u32,
    pub file_size: u64,
    pub kmer_size: u8,
    pub frac_max: u64,
    pub bucket_entry_counts: [u64; BUCKET_COUNT],
}

/// Configuration for building a .jam file from input sequence files.
#[derive(Clone)]
pub struct BuildConfig {
    pub kmer_size: u8,
    pub fscale: u64,
    pub num_threads: usize,
    pub memory: usize,
    pub temp_dir_base: Option<PathBuf>,
    pub min_entropy: f64,
    pub singleton: bool,
    pub bias_table: Option<Arc<BiasTable>>,
    pub show_progress: bool,
}

impl Default for BuildConfig {
    fn default() -> Self {
        Self {
            kmer_size: 21,
            fscale: 1000,
            num_threads: 1,
            memory: 4,
            temp_dir_base: None,
            min_entropy: 0.0,
            singleton: false,
            bias_table: None,
            show_progress: false,
        }
    }
}

#[derive(Debug, thiserror::Error)]
pub enum BuildError {
    #[error("Sketch error: {0}")]
    Sketch(#[from] crate::sketch::SketchError),

    #[error("Compact error: {0}")]
    Compact(#[from] CompactError),
}

/// Build a .jam file from input sequence files.
///
/// This is the main entry point for creating a .jam database. It:
/// 1. Sketches input sequences into partitioned bucket files
/// 2. Sorts entries and builds BinaryFuse8 filters per bucket
/// 3. Writes the final .jam file with header, bucket table, entries, and filters
pub fn build(
    input_files: &[PathBuf],
    output_path: &Path,
    config: &BuildConfig,
) -> Result<IndexStats, BuildError> {
    let sketch_config = SketchConfig {
        kmer_size: config.kmer_size,
        fscale: config.fscale,
        num_threads: config.num_threads,
        memory: config.memory,
        temp_dir_base: config.temp_dir_base.clone(),
        min_entropy: config.min_entropy,
        singleton: config.singleton,
        bias_table: config.bias_table.clone(),
        send_timeout: Duration::from_millis(1),
        show_progress: config.show_progress,
    };

    let sketch_result = crate::sketch::run(input_files, &sketch_config)?;

    let compact_config = CompactConfig {
        num_threads: config.num_threads,
    };

    let stats = run(
        output_path,
        &sketch_result,
        &compact_config,
        config.kmer_size,
        config.bias_table.as_deref(),
    )?;

    Ok(stats)
}

#[derive(Debug, thiserror::Error)]
pub enum CompactError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),

    #[error("Filter construction failed for bucket {bucket}: {message}")]
    FilterConstruction { bucket: usize, message: String },
}

/// Run the two-stage compact pipeline: sort+filter, then mmap write.
pub fn run(
    output_path: &Path,
    sketch_result: &SketchResult,
    _config: &CompactConfig,
    kmer_size: u8,
    bias_table: Option<&BiasTable>,
) -> Result<IndexStats, CompactError> {
    let temp_path = sketch_result.temp_dir.path();

    // Stage 1: Parallel sort + filter construction
    let processed: Result<Vec<ProcessedBucket>, CompactError> = (0..BUCKET_COUNT)
        .into_par_iter()
        .map(|bucket_id| {
            let bucket_path = temp_path.join(format!("bucket_{bucket_id:03}.bin"));

            // Read entries from temp bucket file
            let mut entries = if bucket_path.exists() {
                read_entries(&bucket_path)?
            } else {
                Vec::new()
            };

            // Sort entries (by hash, then by sample_id - Entry derives Ord)
            entries.sort_unstable();

            // Deduplicate entries by (hash, sample_id) - Entry derives Eq
            entries.dedup();

            // Count unique hashes per sample from deduped entries
            let mut sample_hash_counts: std::collections::HashMap<u32, u64> =
                std::collections::HashMap::new();
            for entry in &entries {
                *sample_hash_counts.entry(entry.sample_id).or_insert(0) += 1;
            }

            // Extract unique hashes for filter construction
            let unique_hashes = extract_unique_hashes(&entries);
            let unique_hash_count = unique_hashes.len() as u64;

            // Build BinaryFuse8 filter from unique hashes
            let filter_bytes = if unique_hashes.is_empty() {
                // Empty filter for empty bucket
                Vec::new()
            } else {
                let filter = BinaryFuse8::try_from(&unique_hashes[..]).map_err(|e| {
                    CompactError::FilterConstruction {
                        bucket: bucket_id,
                        message: format!("{e:?}"),
                    }
                })?;
                serialize_filter(&filter)
            };

            // Write sorted (deduped) entries back to temp file
            if !entries.is_empty() {
                write_entries(&bucket_path, &entries)?;
            }

            Ok(ProcessedBucket {
                bucket_id,
                entry_count: entries.len() as u64,
                unique_hash_count,
                filter_bytes,
                sample_hash_counts,
            })
        })
        .collect();

    let processed = processed?;

    // Aggregate unique hash counts per sample from all buckets
    let mut aggregated_sample_sizes: Vec<u64> = vec![0u64; sketch_result.sample_count as usize];
    for bucket in &processed {
        for (&sample_id, &count) in &bucket.sample_hash_counts {
            if (sample_id as usize) < aggregated_sample_sizes.len() {
                aggregated_sample_sizes[sample_id as usize] += count;
            }
        }
    }

    // Stage 2: mmap write of final .jam file
    let entries_size: u64 = processed
        .iter()
        .map(|b| b.entry_count * ENTRY_SIZE as u64)
        .sum();
    let filters_size: u64 = processed.iter().map(|b| b.filter_bytes.len() as u64).sum();
    let bias_size: u64 = if bias_table.is_some() {
        BIAS_TABLE_SERIALIZED_SIZE as u64
    } else {
        0
    };
    let sample_names_bytes = serialize_sample_names(&sketch_result.sample_names);
    let sample_sizes_bytes = serialize_sample_sizes(&aggregated_sample_sizes);
    let total_size = HEADER_SIZE as u64
        + BUCKET_TABLE_SIZE as u64
        + entries_size
        + filters_size
        + bias_size
        + sample_names_bytes.len() as u64
        + sample_sizes_bytes.len() as u64;

    // Allocate output file with read+write access for mmap
    let file = std::fs::OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(true)
        .open(output_path)?;
    file.set_len(total_size)?;

    // mmap the output file (writable)
    let mut mmap = unsafe { MmapMut::map_mut(&file)? };

    // Write entries section (sequential, bucket 0..255)
    let mut entry_offset = DATA_START;
    let mut bucket_metas = vec![BucketMeta::default(); BUCKET_COUNT];
    let mut total_unique_hashes = 0u64;

    for bucket in &processed {
        let bucket_path = temp_path.join(format!("bucket_{:03}.bin", bucket.bucket_id));

        if bucket.entry_count > 0 {
            // Read sorted entries from temp file, copy directly into mmap
            let entries = read_entries(&bucket_path)?;
            let entry_bytes = bytemuck::cast_slice::<Entry, u8>(&entries);
            mmap[entry_offset..entry_offset + entry_bytes.len()].copy_from_slice(entry_bytes);
        }

        bucket_metas[bucket.bucket_id] = BucketMeta {
            entry_offset: entry_offset as u64,
            entry_count: bucket.entry_count,
            filter_offset: 0, // filled in next pass
            filter_size: 0,   // filled in next pass
        };

        entry_offset += (bucket.entry_count as usize) * ENTRY_SIZE;
        total_unique_hashes += bucket.unique_hash_count;
    }

    // Write filters section (sequential, bucket 0..255)
    let mut filter_offset = entry_offset;
    let filters_start = filter_offset;

    for bucket in &processed {
        if !bucket.filter_bytes.is_empty() {
            mmap[filter_offset..filter_offset + bucket.filter_bytes.len()]
                .copy_from_slice(&bucket.filter_bytes);
        }
        bucket_metas[bucket.bucket_id].filter_offset = filter_offset as u64;
        bucket_metas[bucket.bucket_id].filter_size = bucket.filter_bytes.len() as u64;
        filter_offset += bucket.filter_bytes.len();
    }

    // Write bias table section (after filters)
    let (bias_table_offset, bias_table_size, flags) = if let Some(bias) = bias_table {
        let bias_bytes = bias.to_bytes();
        let offset = filter_offset;
        mmap[offset..offset + bias_bytes.len()].copy_from_slice(&bias_bytes);
        (offset as u64, bias_bytes.len() as u64, FLAG_HAS_BIAS_TABLE)
    } else {
        (0, 0, 0)
    };

    // Write sample names section (after bias table if present)
    let sample_names_offset = filter_offset + bias_table_size as usize;
    mmap[sample_names_offset..sample_names_offset + sample_names_bytes.len()]
        .copy_from_slice(&sample_names_bytes);

    // Write sample sizes section
    let sample_sizes_offset = sample_names_offset + sample_names_bytes.len();
    mmap[sample_sizes_offset..sample_sizes_offset + sample_sizes_bytes.len()]
        .copy_from_slice(&sample_sizes_bytes);

    // Write bucket table at HEADER_SIZE offset
    let table_bytes = bytemuck::cast_slice::<BucketMeta, u8>(&bucket_metas);
    mmap[HEADER_SIZE..HEADER_SIZE + BUCKET_TABLE_SIZE].copy_from_slice(table_bytes);

    // Compute total entries
    let total_entries: u64 = processed.iter().map(|b| b.entry_count).sum();

    // Build header
    let header = Header {
        magic: MAGIC,
        version: VERSION,
        flags,
        entry_count: total_entries,
        unique_hash_count: total_unique_hashes,
        sample_count: sketch_result.sample_count,
        bucket_count: BUCKET_COUNT as u16,
        bucket_bits: 8,
        entry_size: ENTRY_SIZE as u8,
        hash_threshold: sketch_result.frac_max,
        kmer_size,
        _param_reserved: [0; 7],
        bucket_table_offset: HEADER_SIZE as u64,
        entries_offset: DATA_START as u64,
        filters_offset: filters_start as u64,
        bias_table_offset,
        entries_size,
        filters_size,
        bias_table_size,
        sample_names_offset: sample_names_offset as u64,
        sample_names_size: sample_names_bytes.len() as u64,
        sample_sizes_offset: sample_sizes_offset as u64,
        sample_sizes_size: sample_sizes_bytes.len() as u64,
        _padding: [0; 16],
    };

    // Write header at offset 0
    let header_bytes = bytemuck::bytes_of(&header);
    mmap[..HEADER_SIZE].copy_from_slice(header_bytes);

    // Flush
    mmap.flush()?;
    drop(mmap);

    // Build bucket entry counts array
    let mut bucket_entry_counts = [0u64; BUCKET_COUNT];
    for bucket in &processed {
        bucket_entry_counts[bucket.bucket_id] = bucket.entry_count;
    }

    Ok(IndexStats {
        total_entries,
        unique_hashes: total_unique_hashes,
        sample_count: sketch_result.sample_count,
        file_size: total_size,
        kmer_size,
        frac_max: sketch_result.frac_max,
        bucket_entry_counts,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::format::bucket_id;
    use crate::sketch::{SketchConfig, run as sketch_run};
    use std::fs::File;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_fasta(seqs: &[(&str, &str)]) -> NamedTempFile {
        let mut f = NamedTempFile::with_suffix(".fa").unwrap();
        for (name, seq) in seqs {
            writeln!(f, ">{name}").unwrap();
            writeln!(f, "{seq}").unwrap();
        }
        f
    }

    #[test]
    fn test_compact_basic() {
        let input = make_fasta(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
        let sketch_config = SketchConfig {
            kmer_size: 11,
            fscale: 1,
            num_threads: 2,
            memory: 1,
            ..Default::default()
        };

        let sketch_result = sketch_run(&[input.path().to_path_buf()], &sketch_config).unwrap();

        let output_dir = tempfile::tempdir().unwrap();
        let output_path = output_dir.path().join("test.jam");

        let compact_config = CompactConfig::default();
        let stats = run(&output_path, &sketch_result, &compact_config, 11, None).unwrap();

        assert!(stats.total_entries > 0);
        assert_eq!(stats.sample_count, 1);
        assert_eq!(stats.kmer_size, 11);
        assert!(output_path.exists());

        // Verify file size matches
        let metadata = std::fs::metadata(&output_path).unwrap();
        assert_eq!(metadata.len(), stats.file_size);
    }

    #[test]
    fn test_compact_empty_buckets() {
        // Use a very restrictive fscale to ensure most buckets are empty
        let input = make_fasta(&[("seq1", "ATCGATCGATCGATCGATCG")]);
        let sketch_config = SketchConfig {
            kmer_size: 11,
            fscale: 1_000_000, // Very restrictive - few hashes pass
            num_threads: 1,
            memory: 1,
            ..Default::default()
        };

        let sketch_result = sketch_run(&[input.path().to_path_buf()], &sketch_config).unwrap();

        let output_dir = tempfile::tempdir().unwrap();
        let output_path = output_dir.path().join("test.jam");

        let compact_config = CompactConfig::default();
        let result = run(&output_path, &sketch_result, &compact_config, 11, None);

        // Should succeed even with mostly empty buckets
        assert!(result.is_ok());
    }

    #[test]
    fn test_compact_header_validation() {
        let input = make_fasta(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
        let sketch_config = SketchConfig {
            kmer_size: 21,
            fscale: 10,
            num_threads: 1,
            memory: 1,
            ..Default::default()
        };

        let sketch_result = sketch_run(&[input.path().to_path_buf()], &sketch_config).unwrap();

        let output_dir = tempfile::tempdir().unwrap();
        let output_path = output_dir.path().join("test.jam");

        let compact_config = CompactConfig::default();
        run(&output_path, &sketch_result, &compact_config, 21, None).unwrap();

        // Read back and validate header
        let file = File::open(&output_path).unwrap();
        let mmap = unsafe { memmap2::Mmap::map(&file).unwrap() };

        let header: &Header = bytemuck::from_bytes(&mmap[..HEADER_SIZE]);
        assert!(header.validate().is_ok());
        assert_eq!(header.kmer_size, 21);
        assert_eq!(header.sample_count, 1);
    }

    #[test]
    fn test_compact_entry_sorting() {
        // Create input that will generate multiple entries per bucket
        let input = make_fasta(&[
            ("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"),
            ("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
        ]);
        let sketch_config = SketchConfig {
            kmer_size: 11,
            fscale: 1, // Keep all hashes
            num_threads: 2,
            memory: 1,
            ..Default::default()
        };

        let sketch_result = sketch_run(&[input.path().to_path_buf()], &sketch_config).unwrap();

        let output_dir = tempfile::tempdir().unwrap();
        let output_path = output_dir.path().join("test.jam");

        let compact_config = CompactConfig::default();
        run(&output_path, &sketch_result, &compact_config, 11, None).unwrap();

        // Read back and verify entries are sorted within each bucket
        let file = File::open(&output_path).unwrap();
        let mmap = unsafe { memmap2::Mmap::map(&file).unwrap() };

        let bucket_table: &[BucketMeta] =
            bytemuck::cast_slice(&mmap[HEADER_SIZE..HEADER_SIZE + BUCKET_TABLE_SIZE]);

        for (i, meta) in bucket_table.iter().enumerate() {
            if meta.entry_count == 0 {
                continue;
            }

            let start = meta.entry_offset as usize;
            let end = start + (meta.entry_count as usize) * ENTRY_SIZE;
            let entries: &[Entry] = bytemuck::cast_slice(&mmap[start..end]);

            // Verify all entries belong to this bucket
            for entry in entries {
                assert_eq!(bucket_id(entry.hash), i, "Entry in wrong bucket");
            }

            // Verify entries are sorted
            for window in entries.windows(2) {
                assert!(
                    window[0] <= window[1],
                    "Entries not sorted in bucket {i}: {:?} > {:?}",
                    window[0],
                    window[1]
                );
            }
        }
    }

    #[test]
    fn test_build_basic() {
        let input = make_fasta(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
        let output_dir = tempfile::tempdir().unwrap();
        let output_path = output_dir.path().join("test.jam");

        let config = BuildConfig {
            kmer_size: 11,
            fscale: 1,
            num_threads: 2,
            memory: 1,
            ..Default::default()
        };

        let stats = build(&[input.path().to_path_buf()], &output_path, &config).unwrap();

        assert!(stats.total_entries > 0);
        assert_eq!(stats.sample_count, 1);
        assert_eq!(stats.kmer_size, 11);
        assert!(output_path.exists());
    }

    #[test]
    fn test_build_singleton_mode() {
        let input = make_fasta(&[
            ("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
            ("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
        ]);
        let output_dir = tempfile::tempdir().unwrap();
        let output_path = output_dir.path().join("test.jam");

        let config = BuildConfig {
            kmer_size: 11,
            fscale: 1,
            singleton: true,
            num_threads: 2,
            memory: 1,
            ..Default::default()
        };

        let stats = build(&[input.path().to_path_buf()], &output_path, &config).unwrap();

        assert_eq!(
            stats.sample_count, 2,
            "Singleton mode should create one sample per sequence"
        );
    }

    #[test]
    fn test_build_multiple_files() {
        let input1 = make_fasta(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
        let input2 = make_fasta(&[("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")]);
        let output_dir = tempfile::tempdir().unwrap();
        let output_path = output_dir.path().join("test.jam");

        let config = BuildConfig {
            kmer_size: 11,
            fscale: 1,
            num_threads: 2,
            memory: 1,
            ..Default::default()
        };

        let stats = build(
            &[input1.path().to_path_buf(), input2.path().to_path_buf()],
            &output_path,
            &config,
        )
        .unwrap();

        assert_eq!(stats.sample_count, 2);
        assert!(stats.total_entries > 0);
    }

    #[test]
    fn test_build_with_entropy_filter() {
        // Low complexity sequence - should be filtered by entropy
        let input = make_fasta(&[
            ("low_complexity", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
            ("high_complexity", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
        ]);
        let output_dir = tempfile::tempdir().unwrap();
        let output_path = output_dir.path().join("test.jam");

        let config = BuildConfig {
            kmer_size: 11,
            fscale: 1,
            min_entropy: 1.5, // Filter low-complexity k-mers
            num_threads: 1,
            memory: 1,
            ..Default::default()
        };

        let stats = build(&[input.path().to_path_buf()], &output_path, &config).unwrap();

        // Should have entries from the high-complexity sequence
        assert!(stats.total_entries > 0);
    }

    #[test]
    fn test_build_with_bias_table() {
        // Create bias table from positive/negative sequences
        let pos_fasta = make_fasta(&[("pos", "ATATATATATATATATATATATATATATATAT")]);
        let neg_fasta = make_fasta(&[("neg", "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC")]);

        let bias_table =
            crate::bias::BiasTable::build(pos_fasta.path(), neg_fasta.path(), 0.5).unwrap();

        // Build JAM file with bias table - use ATATAT-rich sequence that passes bias filter
        let input = make_fasta(&[("seq1", "ATATATATATATATATATATATATATATATATATAT")]);
        let output_dir = tempfile::tempdir().unwrap();
        let output_path = output_dir.path().join("test.jam");

        let config = BuildConfig {
            kmer_size: 11,
            fscale: 1,
            num_threads: 1,
            memory: 1,
            bias_table: Some(std::sync::Arc::new(bias_table.clone())),
            ..Default::default()
        };

        build(&[input.path().to_path_buf()], &output_path, &config).unwrap();

        // Verify the bias table is embedded by reading the file back
        let reader = crate::reader::JamReader::open(&output_path).unwrap();
        assert!(reader.has_bias_table());

        let embedded_bias = reader.bias_table().unwrap();
        assert_eq!(*embedded_bias, bias_table);
    }

    #[test]
    fn test_build_without_bias_table() {
        let input = make_fasta(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
        let output_dir = tempfile::tempdir().unwrap();
        let output_path = output_dir.path().join("test.jam");

        let config = BuildConfig {
            kmer_size: 11,
            fscale: 1,
            num_threads: 1,
            memory: 1,
            bias_table: None,
            ..Default::default()
        };

        build(&[input.path().to_path_buf()], &output_path, &config).unwrap();

        // Verify no bias table is embedded
        let reader = crate::reader::JamReader::open(&output_path).unwrap();
        assert!(!reader.has_bias_table());
        assert!(reader.bias_table().is_none());
    }
}
