//! Hash analysis tool for SLURM - finds duplicates across ALL chunks
//!
//! Each job processes one BUCKET (not chunk), reading all chunk files
//! but only keeping hashes that belong to that bucket.
//!
//! Outputs:
//! - {prefix}_duplicates_bucket_{bucket}.csv: hash,count for collisions
//! - {prefix}_distribution_bucket_{bucket}.csv: sub-bucket distribution
//!
//! Usage: analyze_hashes <bucket> <input_dir> <output_dir> --prefix <prefix> --chunks <N>

use byteorder::{LittleEndian, ReadBytesExt};
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;

const NUM_BUCKETS: usize = 2048; // Main buckets (top 11 bits) - one SLURM job per bucket
const BUCKET_SHIFT: u32 = 64 - 11;

const SUB_BUCKETS: usize = 256; // Sub-buckets within each bucket for distribution analysis
const SUB_BUCKET_SHIFT: u32 = 64 - 11 - 8; // Next 8 bits after bucket

fn main() {
    let args: Vec<String> = std::env::args().collect();

    if args.len() < 4 {
        eprintln!(
            "Usage: {} <bucket> <input_dir> <output_dir> --prefix <prefix> --chunks <N>",
            args[0]
        );
        eprintln!("  bucket:    Bucket ID (0 to 2047) - use SLURM_ARRAY_TASK_ID");
        eprintln!("  input_dir: Directory containing {{prefix}}_*.bin files");
        eprintln!("  output_dir: Output directory for CSV files");
        eprintln!("  --prefix:  Hash file prefix (default: jam_hashes)");
        eprintln!("  --chunks:  Number of chunk files (default: 2000)");
        std::process::exit(1);
    }

    let bucket: usize = args[1].parse().expect("Invalid bucket number");
    let input_dir = PathBuf::from(&args[2]);
    let output_dir = PathBuf::from(&args[3]);

    let prefix = args
        .iter()
        .position(|a| a == "--prefix")
        .and_then(|i| args.get(i + 1))
        .map(|s| s.as_str())
        .unwrap_or("jam_hashes");

    let num_chunks: usize = args
        .iter()
        .position(|a| a == "--chunks")
        .and_then(|i| args.get(i + 1))
        .and_then(|s| s.parse().ok())
        .unwrap_or(2000);

    if bucket >= NUM_BUCKETS {
        eprintln!("Bucket must be 0-{}", NUM_BUCKETS - 1);
        std::process::exit(1);
    }

    std::fs::create_dir_all(&output_dir).expect("Failed to create output dir");

    eprintln!(
        "Processing bucket {} (reading {} chunks, prefix={})",
        bucket, num_chunks, prefix
    );

    // Read ALL chunk files, keeping only hashes for our bucket
    let mut hashes: Vec<u64> = Vec::new();
    let mut sub_buckets = vec![0u64; SUB_BUCKETS];
    let mut files_read = 0;

    for chunk in 0..num_chunks {
        let path = input_dir.join(format!("{}_{}.bin", prefix, chunk));
        if !path.exists() {
            continue;
        }

        let file = match File::open(&path) {
            Ok(f) => f,
            Err(_) => continue,
        };

        let mut reader = BufReader::with_capacity(64 * 1024 * 1024, file);

        while let Ok(hash) = reader.read_u64::<LittleEndian>() {
            let hash_bucket = (hash >> BUCKET_SHIFT) as usize;
            if hash_bucket == bucket {
                // Sub-bucket for distribution (next 8 bits)
                let sub = ((hash >> SUB_BUCKET_SHIFT) & 0xFF) as usize;
                sub_buckets[sub] += 1;

                hashes.push(hash);
            }
        }

        files_read += 1;
        if chunk % 200 == 0 {
            eprintln!("  Read {}/{} chunks, {} hashes so far", chunk, num_chunks, hashes.len());
        }
    }

    eprintln!(
        "  Read {} files, collected {} hashes for bucket {}",
        files_read,
        hashes.len(),
        bucket
    );

    // Sort for duplicate detection
    eprintln!("  Sorting...");
    hashes.sort_unstable();

    // Find duplicates
    eprintln!("  Finding duplicates...");
    let mut duplicates: Vec<(u64, u32)> = Vec::new();
    let mut i = 0;
    while i < hashes.len() {
        let h = hashes[i];
        let start = i;
        while i < hashes.len() && hashes[i] == h {
            i += 1;
        }
        if i - start > 1 {
            duplicates.push((h, (i - start) as u32));
        }
    }

    let total = hashes.len();
    let dup_count: u64 = duplicates.iter().map(|(_, c)| *c as u64 - 1).sum();

    // Write distribution CSV
    let dist_path = output_dir.join(format!("{}_distribution_bucket_{}.csv", prefix, bucket));
    let mut f = BufWriter::new(File::create(&dist_path).unwrap());
    writeln!(f, "sub_bucket,count").unwrap();
    for (sb, &c) in sub_buckets.iter().enumerate() {
        writeln!(f, "{},{}", sb, c).unwrap();
    }

    // Write duplicates CSV
    let dup_path = output_dir.join(format!("{}_duplicates_bucket_{}.csv", prefix, bucket));
    let mut f = BufWriter::new(File::create(&dup_path).unwrap());
    writeln!(f, "hash,count").unwrap();
    for (h, c) in &duplicates {
        writeln!(f, "{},{}", h, c).unwrap();
    }

    println!(
        "bucket={} total_hashes={} collisions={} unique_collision_hashes={}",
        bucket,
        total,
        dup_count,
        duplicates.len()
    );
}
