//! Converter tool: JAM format v2 → v3
//!
//! v3 format changes:
//! - Each bucket's filter + entries are contiguous and page-aligned
//! - Allows per-bucket mmap that can be dropped after processing
//!
//! Usage: convert_v2_to_v3 <input.jam> <output.jam>

use bytemuck::{Pod, Zeroable};
use std::fs::File;
use std::io::{Read, Write};
use std::path::PathBuf;

use clap::Parser;

const MAGIC: [u8; 4] = *b"JAM\0";
const VERSION_V2: u32 = 2;
const VERSION_V3: u32 = 3;
const BUCKET_COUNT: usize = 256;
const HEADER_SIZE: usize = 160;
const BUCKET_META_SIZE: usize = 32;
const BUCKET_TABLE_SIZE: usize = BUCKET_COUNT * BUCKET_META_SIZE;
const ENTRY_SIZE: usize = 12;
const PAGE_SIZE: usize = 4096;

#[inline]
const fn align_to_page(offset: usize) -> usize {
    (offset + PAGE_SIZE - 1) & !(PAGE_SIZE - 1)
}

#[derive(Clone, Copy, Pod, Zeroable)]
#[repr(C)]
struct BucketMeta {
    entry_offset: u64,
    entry_count: u64,
    filter_offset: u64,
    filter_size: u64,
}

#[derive(Clone, Copy, Pod, Zeroable)]
#[repr(C)]
struct Header {
    magic: [u8; 4],
    version: u32,
    flags: u64,
    entry_count: u64,
    unique_hash_count: u64,
    sample_count: u32,
    bucket_count: u16,
    bucket_bits: u8,
    entry_size: u8,
    hash_threshold: u64,
    kmer_size: u8,
    _param_reserved: [u8; 7],
    bucket_table_offset: u64,
    entries_offset: u64,
    filters_offset: u64,
    bias_table_offset: u64,
    entries_size: u64,
    filters_size: u64,
    bias_table_size: u64,
    sample_names_offset: u64,
    sample_names_size: u64,
    sample_sizes_offset: u64,
    sample_sizes_size: u64,
    _padding: [u8; 16],
}

const _: () = assert!(std::mem::size_of::<Header>() == HEADER_SIZE);
const _: () = assert!(std::mem::size_of::<BucketMeta>() == BUCKET_META_SIZE);

#[derive(Debug)]
enum ConvertError {
    Io(std::io::Error),
    FileTooSmall { expected: usize, actual: usize },
    InvalidMagic([u8; 4]),
    InvalidVersion(u32),
    InvalidOffset { field: &'static str, offset: usize, size: usize, file_size: usize },
}

impl std::fmt::Display for ConvertError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(e) => write!(f, "I/O error: {}", e),
            Self::FileTooSmall { expected, actual } => {
                write!(f, "file too small: expected {} bytes, got {}", expected, actual)
            }
            Self::InvalidMagic(m) => write!(f, "invalid magic: {:?}", m),
            Self::InvalidVersion(v) => write!(f, "expected version 2, got {}", v),
            Self::InvalidOffset { field, offset, size, file_size } => {
                write!(f, "invalid {}: offset {} + size {} exceeds file size {}",
                       field, offset, size, file_size)
            }
        }
    }
}

impl From<std::io::Error> for ConvertError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}

#[derive(Parser)]
#[command(name = "convert_v2_to_v3")]
#[command(about = "Convert JAM format v2 to v3")]
struct Args {
    /// Input v2 .jam file
    input: PathBuf,
    /// Output v3 .jam file
    output: PathBuf,
}

fn convert(input_path: &PathBuf, output_path: &PathBuf) -> Result<(), ConvertError> {
    eprintln!("Reading {}...", input_path.display());

    let mut input_file = File::open(input_path)?;
    let mut input_data = Vec::new();
    input_file.read_to_end(&mut input_data)?;

    if input_data.len() < HEADER_SIZE {
        return Err(ConvertError::FileTooSmall {
            expected: HEADER_SIZE,
            actual: input_data.len(),
        });
    }

    let header: Header = *bytemuck::from_bytes(&input_data[..HEADER_SIZE]);

    if header.magic != MAGIC {
        return Err(ConvertError::InvalidMagic(header.magic));
    }

    if header.version != VERSION_V2 {
        return Err(ConvertError::InvalidVersion(header.version));
    }

    eprintln!("  Entries: {}", header.entry_count);
    eprintln!("  Samples: {}", header.sample_count);
    eprintln!("  K-mer:   {}", header.kmer_size);

    // Validate header offsets before processing (guards against corrupt files)
    let file_size = input_data.len();

    if header.bias_table_size > 0 {
        let offset = header.bias_table_offset as usize;
        let size = header.bias_table_size as usize;
        if offset.saturating_add(size) > file_size {
            return Err(ConvertError::InvalidOffset {
                field: "bias_table",
                offset,
                size,
                file_size,
            });
        }
    }

    if header.sample_names_size > 0 {
        let offset = header.sample_names_offset as usize;
        let size = header.sample_names_size as usize;
        if offset.saturating_add(size) > file_size {
            return Err(ConvertError::InvalidOffset {
                field: "sample_names",
                offset,
                size,
                file_size,
            });
        }
    }

    if header.sample_sizes_size > 0 {
        let offset = header.sample_sizes_offset as usize;
        let size = header.sample_sizes_size as usize;
        if offset.saturating_add(size) > file_size {
            return Err(ConvertError::InvalidOffset {
                field: "sample_sizes",
                offset,
                size,
                file_size,
            });
        }
    }

    let table_end = HEADER_SIZE + BUCKET_TABLE_SIZE;
    if input_data.len() < table_end {
        return Err(ConvertError::FileTooSmall {
            expected: table_end,
            actual: input_data.len(),
        });
    }

    let bucket_metas: &[BucketMeta] =
        bytemuck::cast_slice(&input_data[HEADER_SIZE..table_end]);

    // Validate bucket metadata offsets
    for (i, meta) in bucket_metas.iter().enumerate() {
        if meta.filter_size > 0 {
            let offset = meta.filter_offset as usize;
            let size = meta.filter_size as usize;
            if offset.saturating_add(size) > file_size {
                return Err(ConvertError::InvalidOffset {
                    field: "bucket filter",
                    offset,
                    size,
                    file_size,
                });
            }
        }
        if meta.entry_count > 0 {
            let offset = meta.entry_offset as usize;
            let size = (meta.entry_count as usize) * ENTRY_SIZE;
            if offset.saturating_add(size) > file_size {
                return Err(ConvertError::InvalidOffset {
                    field: "bucket entries",
                    offset,
                    size,
                    file_size,
                });
            }
        }
        let _ = i; // suppress unused warning
    }

    // Calculate v3 layout: each bucket region is page-aligned
    let bucket_regions_start = align_to_page(HEADER_SIZE + BUCKET_TABLE_SIZE);
    let mut current_offset = bucket_regions_start;
    let mut new_bucket_offsets: Vec<usize> = Vec::with_capacity(BUCKET_COUNT);

    for meta in bucket_metas {
        new_bucket_offsets.push(current_offset);
        let bucket_size = meta.filter_size as usize + (meta.entry_count as usize) * ENTRY_SIZE;
        if bucket_size > 0 {
            current_offset = align_to_page(current_offset + bucket_size);
        }
    }

    // Metadata after buckets
    let metadata_offset = current_offset;
    let total_size = metadata_offset
        + header.bias_table_size as usize
        + header.sample_names_size as usize
        + header.sample_sizes_size as usize;

    eprintln!("Converting to v3 layout...");

    let mut output_data = vec![0u8; total_size];

    // Copy header, update version and offsets
    output_data[..HEADER_SIZE].copy_from_slice(&input_data[..HEADER_SIZE]);
    output_data[4..8].copy_from_slice(&VERSION_V3.to_le_bytes());

    // Write bucket regions (filter + entries together) and update bucket table
    let mut _entries_size_total = 0u64;
    let mut _filters_size_total = 0u64;

    for (i, meta) in bucket_metas.iter().enumerate() {
        let new_offset = new_bucket_offsets[i];
        let filter_size = meta.filter_size as usize;
        let entries_size = (meta.entry_count as usize) * ENTRY_SIZE;

        // Copy filter
        if filter_size > 0 {
            let src = meta.filter_offset as usize;
            output_data[new_offset..new_offset + filter_size]
                .copy_from_slice(&input_data[src..src + filter_size]);
        }

        // Copy entries immediately after filter
        let entry_dst = new_offset + filter_size;
        if entries_size > 0 {
            let src = meta.entry_offset as usize;
            output_data[entry_dst..entry_dst + entries_size]
                .copy_from_slice(&input_data[src..src + entries_size]);
        }

        // Update bucket meta
        let new_meta = BucketMeta {
            filter_offset: new_offset as u64,
            filter_size: meta.filter_size,
            entry_offset: entry_dst as u64,
            entry_count: meta.entry_count,
        };
        let meta_offset = HEADER_SIZE + i * BUCKET_META_SIZE;
        output_data[meta_offset..meta_offset + BUCKET_META_SIZE]
            .copy_from_slice(bytemuck::bytes_of(&new_meta));

        _entries_size_total += entries_size as u64;
        _filters_size_total += filter_size as u64;
    }

    // Copy metadata sections with updated offsets
    let mut meta_dst = metadata_offset;

    if header.bias_table_size > 0 {
        let src = header.bias_table_offset as usize;
        let size = header.bias_table_size as usize;
        output_data[meta_dst..meta_dst + size]
            .copy_from_slice(&input_data[src..src + size]);
        // Update bias_table_offset in header (offset 80)
        output_data[80..88].copy_from_slice(&(meta_dst as u64).to_le_bytes());
        meta_dst += size;
    }

    if header.sample_names_size > 0 {
        let src = header.sample_names_offset as usize;
        let size = header.sample_names_size as usize;
        output_data[meta_dst..meta_dst + size]
            .copy_from_slice(&input_data[src..src + size]);
        // Update sample_names_offset in header (offset 112)
        output_data[112..120].copy_from_slice(&(meta_dst as u64).to_le_bytes());
        meta_dst += size;
    }

    if header.sample_sizes_size > 0 {
        let src = header.sample_sizes_offset as usize;
        let size = header.sample_sizes_size as usize;
        output_data[meta_dst..meta_dst + size]
            .copy_from_slice(&input_data[src..src + size]);
        // Update sample_sizes_offset in header (offset 128)
        output_data[128..136].copy_from_slice(&(meta_dst as u64).to_le_bytes());
    }

    // Update header offsets for v3 layout
    // bucket_table_offset is always immediately after header in both v2 and v3
    output_data[56..64].copy_from_slice(&(HEADER_SIZE as u64).to_le_bytes()); // bucket_table_offset
    let bucket_start = bucket_regions_start as u64;
    output_data[64..72].copy_from_slice(&bucket_start.to_le_bytes()); // entries_offset
    output_data[72..80].copy_from_slice(&bucket_start.to_le_bytes()); // filters_offset

    eprintln!("Writing {}...", output_path.display());

    let mut output_file = File::create(output_path)?;
    output_file.write_all(&output_data)?;

    let overhead = output_data.len() as i64 - input_data.len() as i64;
    let overhead_pct = (output_data.len() as f64 / input_data.len() as f64 - 1.0) * 100.0;

    eprintln!("Done!");
    eprintln!("  Input:    {} bytes", input_data.len());
    eprintln!("  Output:   {} bytes", output_data.len());
    eprintln!("  Overhead: {} bytes ({:.1}%)", overhead, overhead_pct);

    Ok(())
}

fn main() {
    let args = Args::parse();

    if let Err(e) = convert(&args.input, &args.output) {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
