use crate::bias::HashBiasTable;
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

pub const FILTER_DESCRIPTOR_SIZE: usize = 20;

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

#[derive(Clone)]
pub struct BuildConfig {
    pub kmer_size: u8,
    pub fscale: u64,
    pub num_threads: usize,
    pub memory: usize,
    pub temp_dir_base: Option<PathBuf>,
    pub min_entropy: f64,
    pub singleton: bool,
    pub bias_table: Option<Arc<HashBiasTable>>,
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
        config.min_entropy,
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

pub fn run(
    output_path: &Path,
    sketch_result: &SketchResult,
    _config: &CompactConfig,
    kmer_size: u8,
    bias_table: Option<&HashBiasTable>,
    min_entropy: f64,
) -> Result<IndexStats, CompactError> {
    let temp_path = sketch_result.temp_dir.path();

    let processed: Result<Vec<ProcessedBucket>, CompactError> = (0..BUCKET_COUNT)
        .into_par_iter()
        .map(|bucket_id| {
            let bucket_path = temp_path.join(format!("bucket_{bucket_id:03}.bin"));

            let mut entries = if bucket_path.exists() {
                read_entries(&bucket_path)?
            } else {
                Vec::new()
            };

            entries.sort_unstable();
            entries.dedup();

            let mut sample_hash_counts: std::collections::HashMap<u32, u64> =
                std::collections::HashMap::new();
            for entry in &entries {
                *sample_hash_counts.entry(entry.sample_id).or_insert(0) += 1;
            }

            let unique_hashes = extract_unique_hashes(&entries);
            let unique_hash_count = unique_hashes.len() as u64;

            let filter_bytes = if unique_hashes.is_empty() {
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

    let mut aggregated_sample_sizes: Vec<u64> = vec![0u64; sketch_result.sample_count as usize];
    for bucket in &processed {
        for (&sample_id, &count) in &bucket.sample_hash_counts {
            if (sample_id as usize) < aggregated_sample_sizes.len() {
                aggregated_sample_sizes[sample_id as usize] += count;
            }
        }
    }

    use crate::format::align_to_page;

    let bias_size: u64 = bias_table.map(|b| b.to_bytes().len() as u64).unwrap_or(0);
    let sample_names_bytes = serialize_sample_names(&sketch_result.sample_names);
    let sample_sizes_bytes = serialize_sample_sizes(&aggregated_sample_sizes);

    let bucket_regions_start = align_to_page(DATA_START);
    let mut current_offset = bucket_regions_start;
    let mut bucket_offsets = Vec::with_capacity(BUCKET_COUNT);

    for bucket in &processed {
        bucket_offsets.push(current_offset);
        let bucket_size = bucket.filter_bytes.len() + (bucket.entry_count as usize) * ENTRY_SIZE;
        if bucket_size > 0 {
            current_offset = align_to_page(current_offset + bucket_size);
        }
    }

    let metadata_offset = current_offset;
    let total_size =
        metadata_offset + bias_size as usize + sample_names_bytes.len() + sample_sizes_bytes.len();

    let file = std::fs::OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(true)
        .open(output_path)?;
    file.set_len(total_size as u64)?;

    let mut mmap = unsafe { MmapMut::map_mut(&file)? };

    let mut bucket_metas = vec![BucketMeta::default(); BUCKET_COUNT];
    let mut total_unique_hashes = 0u64;
    let mut entries_size = 0u64;
    let mut filters_size = 0u64;

    for bucket in &processed {
        let bucket_offset = bucket_offsets[bucket.bucket_id];
        let filter_size = bucket.filter_bytes.len();
        let entries_bytes_len = (bucket.entry_count as usize) * ENTRY_SIZE;

        if !bucket.filter_bytes.is_empty() {
            mmap[bucket_offset..bucket_offset + filter_size].copy_from_slice(&bucket.filter_bytes);
        }

        let entry_offset = bucket_offset + filter_size;
        if bucket.entry_count > 0 {
            let bucket_path = temp_path.join(format!("bucket_{:03}.bin", bucket.bucket_id));
            let entries = read_entries(&bucket_path)?;
            let entry_bytes = bytemuck::cast_slice::<Entry, u8>(&entries);
            mmap[entry_offset..entry_offset + entry_bytes.len()].copy_from_slice(entry_bytes);
        }

        bucket_metas[bucket.bucket_id] = BucketMeta {
            filter_offset: bucket_offset as u64,
            filter_size: filter_size as u64,
            entry_offset: entry_offset as u64,
            entry_count: bucket.entry_count,
        };

        total_unique_hashes += bucket.unique_hash_count;
        entries_size += entries_bytes_len as u64;
        filters_size += filter_size as u64;
    }

    let mut meta_offset = metadata_offset;

    let (bias_table_offset, bias_table_size, flags) = if let Some(bias) = bias_table {
        let bias_bytes = bias.to_bytes();
        mmap[meta_offset..meta_offset + bias_bytes.len()].copy_from_slice(&bias_bytes);
        let offset = meta_offset;
        meta_offset += bias_bytes.len();
        (offset as u64, bias_bytes.len() as u64, FLAG_HAS_BIAS_TABLE)
    } else {
        (0, 0, 0)
    };

    let sample_names_offset = meta_offset;
    mmap[sample_names_offset..sample_names_offset + sample_names_bytes.len()]
        .copy_from_slice(&sample_names_bytes);
    meta_offset += sample_names_bytes.len();

    let sample_sizes_offset = meta_offset;
    mmap[sample_sizes_offset..sample_sizes_offset + sample_sizes_bytes.len()]
        .copy_from_slice(&sample_sizes_bytes);

    let table_bytes = bytemuck::cast_slice::<BucketMeta, u8>(&bucket_metas);
    mmap[HEADER_SIZE..HEADER_SIZE + BUCKET_TABLE_SIZE].copy_from_slice(table_bytes);

    let total_entries: u64 = processed.iter().map(|b| b.entry_count).sum();

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
        _param_reserved: [0; 3],
        min_entropy: min_entropy as f32,
        bucket_table_offset: HEADER_SIZE as u64,
        entries_offset: bucket_regions_start as u64,
        filters_offset: bucket_regions_start as u64,
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

    let header_bytes = bytemuck::bytes_of(&header);
    mmap[..HEADER_SIZE].copy_from_slice(header_bytes);

    mmap.flush()?;
    drop(mmap);

    let mut bucket_entry_counts = [0u64; BUCKET_COUNT];
    for bucket in &processed {
        bucket_entry_counts[bucket.bucket_id] = bucket.entry_count;
    }

    Ok(IndexStats {
        total_entries,
        unique_hashes: total_unique_hashes,
        sample_count: sketch_result.sample_count,
        file_size: total_size as u64,
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
        let stats = run(&output_path, &sketch_result, &compact_config, 11, None, 0.0).unwrap();

        assert!(stats.total_entries > 0);
        assert_eq!(stats.sample_count, 1);
        assert_eq!(stats.kmer_size, 11);
        assert!(output_path.exists());

        let metadata = std::fs::metadata(&output_path).unwrap();
        assert_eq!(metadata.len(), stats.file_size);
    }

    #[test]
    fn test_compact_empty_buckets() {
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
        let result = run(&output_path, &sketch_result, &compact_config, 11, None, 0.0);

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
        run(&output_path, &sketch_result, &compact_config, 21, None, 0.0).unwrap();

        let file = File::open(&output_path).unwrap();
        let mmap = unsafe { memmap2::Mmap::map(&file).unwrap() };

        let header: &Header = bytemuck::from_bytes(&mmap[..HEADER_SIZE]);
        assert!(header.validate().is_ok());
        assert_eq!(header.kmer_size, 21);
        assert_eq!(header.sample_count, 1);
    }

    #[test]
    fn test_compact_entry_sorting() {
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
        run(&output_path, &sketch_result, &compact_config, 11, None, 0.0).unwrap();

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

            for entry in entries {
                assert_eq!(bucket_id(entry.hash), i, "Entry in wrong bucket");
            }

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
        let input = make_fasta(&[
            ("low_complexity", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
            ("high_complexity", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
        ]);
        let output_dir = tempfile::tempdir().unwrap();
        let output_path = output_dir.path().join("test.jam");

        let config = BuildConfig {
            kmer_size: 11,
            fscale: 1,
            min_entropy: 1.5,
            num_threads: 1,
            memory: 1,
            ..Default::default()
        };

        let stats = build(&[input.path().to_path_buf()], &output_path, &config).unwrap();

        assert!(stats.total_entries > 0);
    }

    #[test]
    fn test_build_with_bias_table() {
        use crate::bias::{CMSConfig, HashBiasTable, RawHashCounts};

        let pos_fasta = make_fasta(&[("pos", "ATATATATATATATATATATATATATATATAT")]);
        let neg_fasta = make_fasta(&[("neg", "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC")]);

        let config = CMSConfig {
            width: 1024,
            depth: 3,
            k: 11,
            fscale: 1,
        };

        let rc = std::sync::atomic::AtomicU64::new(0);
        let hc = std::sync::atomic::AtomicU64::new(0);
        let pos_raw = RawHashCounts::build(&[pos_fasta.path()], config.clone(), &rc, &hc).unwrap();
        let neg_raw = RawHashCounts::build(&[neg_fasta.path()], config, &rc, &hc).unwrap();
        let bias_table = HashBiasTable::build(&pos_raw, &neg_raw, 1.0, Some(2.0)).unwrap();

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

        let reader = crate::reader::JamReader::open(&output_path).unwrap();
        assert!(!reader.has_bias_table());
        assert!(reader.bias_table().is_none());
    }
}
