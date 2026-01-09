use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use clap::Parser;
use std::env;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

#[derive(Parser)]
#[command(name = "find_duplicates_with_counts")]
#[command(about = "Find duplicate hashes and save with collision counts for analysis")]
struct Args {
    /// Input directory containing unsorted_bin_*.bin files
    input_dir: String,

    /// Output directory for duplicate files
    output_dir: String,

    /// Hash name prefix
    #[arg(long, default_value = "jam_hashes")]
    hash_name: String,

    /// Task ID (0-2047), defaults to SLURM_ARRAY_TASK_ID env var
    #[arg(long)]
    task_id: Option<usize>,

    /// Process all bins sequentially (ignores task_id)
    #[arg(long)]
    all: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    if args.all {
        // Process all bins sequentially
        println!("Processing all 2048 bins sequentially...");
        let mut total_duplicates = 0u64;
        let mut total_collisions = 0u64;

        for task_id in 0..2048 {
            let bindex = task_id / 128;
            let bin_id = task_id % 128;

            let input_file = format!(
                "{}/unsorted_bin_{}_{}_{}.bin",
                args.input_dir, args.hash_name, bindex, bin_id
            );

            let input_path = Path::new(&input_file);
            if !input_path.exists() {
                continue;
            }

            match find_duplicates_with_counts(&input_file, &args.output_dir, bindex, bin_id) {
                Ok((dups, cols)) => {
                    total_duplicates += dups;
                    total_collisions += cols;
                }
                Err(e) => eprintln!("Error processing bin {}_{}: {}", bindex, bin_id, e),
            }
        }

        println!("\n=== Summary ===");
        println!("Total unique duplicate hashes: {}", total_duplicates);
        println!("Total collision count: {}", total_collisions);
        return Ok(());
    }

    // Single bin mode
    let task_id = match args.task_id {
        Some(id) => id,
        None => env::var("SLURM_ARRAY_TASK_ID")
            .map_err(|_| "No task_id provided and SLURM_ARRAY_TASK_ID not set")?
            .parse::<usize>()?,
    };

    if task_id >= 2048 {
        return Err("Task ID must be 0-2047 (16 bindex × 128 bins = 2048 total)".into());
    }

    let bindex = task_id / 128;
    let bin_id = task_id % 128;

    let input_file = format!(
        "{}/unsorted_bin_{}_{}_{}.bin",
        args.input_dir, args.hash_name, bindex, bin_id
    );

    let input_path = Path::new(&input_file);
    if !input_path.exists() {
        eprintln!("Input file not found: {}", input_path.display());
        return Err("Input file not found".into());
    }

    println!(
        "Task {}: Processing bindex={}, bin_id={}",
        task_id, bindex, bin_id
    );

    find_duplicates_with_counts(&input_file, &args.output_dir, bindex, bin_id)?;
    Ok(())
}

/// Returns (unique_duplicates, total_collision_count)
fn find_duplicates_with_counts(
    input_file: &str,
    output_dir: &str,
    bindex: usize,
    bin_id: usize,
) -> Result<(u64, u64), Box<dyn std::error::Error>> {
    // Read all hashes into vector
    let file = File::open(input_file)?;
    let mut reader = BufReader::with_capacity(1024 * 1024 * 64, file);
    let mut hashes = Vec::new();

    while let Ok(hash) = reader.read_u64::<LittleEndian>() {
        hashes.push(hash);
    }

    let total_hashes = hashes.len();
    println!("Bin {}_{}: Read {} hashes", bindex, bin_id, total_hashes);

    // Sort hashes to group duplicates together
    hashes.sort_unstable();

    // Find duplicates with counts
    // Format: Vec<(hash, count)>
    let mut duplicates: Vec<(u64, u32)> = Vec::new();
    let mut last_hash = 0u64;
    let mut count = 0u32;

    for &current_hash in &hashes {
        if current_hash == last_hash {
            count += 1;
        } else {
            if count > 1 {
                duplicates.push((last_hash, count));
            }
            last_hash = current_hash;
            count = 1;
        }
    }

    // Handle the last group
    if count > 1 {
        duplicates.push((last_hash, count));
    }

    // Create output directory if needed
    std::fs::create_dir_all(output_dir)?;

    let unique_dups = duplicates.len() as u64;
    let total_collisions: u64 = duplicates.iter().map(|(_, c)| *c as u64).sum();

    if duplicates.is_empty() {
        println!("Bin {}_{}: No duplicates found", bindex, bin_id);
        return Ok((0, 0));
    }

    // Write duplicates with counts
    // Format: [hash: u64, count: u32] pairs (12 bytes per entry)
    let output_file = format!(
        "{}/duplicates_with_counts_{}_{}.bin",
        output_dir, bindex, bin_id
    );
    let output_path = Path::new(&output_file);

    let file = File::create(&output_path)?;
    let mut writer = BufWriter::with_capacity(1024 * 1024 * 64, file);

    for (hash, count) in &duplicates {
        writer.write_u64::<LittleEndian>(*hash)?;
        writer.write_u32::<LittleEndian>(*count)?;
    }

    writer.flush()?;

    // Also write a summary CSV for easy inspection
    let csv_file = format!("{}/duplicates_summary_{}_{}.csv", output_dir, bindex, bin_id);
    let mut csv_writer = BufWriter::new(File::create(&csv_file)?);
    writeln!(csv_writer, "hash,count")?;
    for (hash, count) in &duplicates {
        writeln!(csv_writer, "{},{}", hash, count)?;
    }
    csv_writer.flush()?;

    println!(
        "Bin {}_{}: {} total hashes, {} unique duplicates, {} total collisions",
        bindex, bin_id, total_hashes, unique_dups, total_collisions
    );

    Ok((unique_dups, total_collisions))
}
