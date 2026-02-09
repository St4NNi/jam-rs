//! Tool for analyzing hash output uniformity across all hashes (not just collisions).
//! This validates that JamHash satisfies MinHash's uniformity requirement.

use byteorder::{LittleEndian, ReadBytesExt};
use clap::Parser;
use rayon::prelude::*;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

#[derive(Parser)]
#[command(name = "hash_uniformity")]
#[command(about = "Analyze uniformity of ALL hash outputs (not just collisions)")]
struct Args {
    /// Input directory containing unsorted_bin_*.bin files
    input_dir: String,

    /// Hash name prefix
    #[arg(long, default_value = "jam_hashes")]
    hash_name: String,

    /// Number of bins for uniformity test (default: 2048)
    #[arg(long, default_value = "2048")]
    num_bins: usize,

    /// Number of threads (default: all available)
    #[arg(long)]
    threads: Option<usize>,
}

/// Statistics for a single input bin file
struct BinStats {
    total_hashes: u64,
    /// Distribution across output bins (for uniformity test)
    bin_counts: Vec<u64>,
    /// For computing hash value statistics
    sum: u128,
    sum_sq: u128,
    min_hash: u64,
    max_hash: u64,
}

impl BinStats {
    fn new(num_bins: usize) -> Self {
        Self {
            total_hashes: 0,
            bin_counts: vec![0; num_bins],
            sum: 0,
            sum_sq: 0,
            min_hash: u64::MAX,
            max_hash: 0,
        }
    }

    fn merge(&mut self, other: &BinStats) {
        self.total_hashes += other.total_hashes;
        for (i, count) in other.bin_counts.iter().enumerate() {
            self.bin_counts[i] += count;
        }
        self.sum += other.sum;
        self.sum_sq += other.sum_sq;
        self.min_hash = self.min_hash.min(other.min_hash);
        self.max_hash = self.max_hash.max(other.max_hash);
    }
}

fn process_bin_file(input_file: &str, num_bins: usize) -> Result<BinStats, Box<dyn std::error::Error + Send + Sync>> {
    let file = File::open(input_file)?;
    let mut reader = BufReader::with_capacity(1024 * 1024 * 64, file);
    let mut stats = BinStats::new(num_bins);

    while let Ok(hash) = reader.read_u64::<LittleEndian>() {
        stats.total_hashes += 1;

        // Bin by modulo for uniformity test
        let bin_idx = (hash as usize) % num_bins;
        stats.bin_counts[bin_idx] += 1;

        // Accumulate for mean/variance (use upper bits to avoid overflow issues)
        let hash_scaled = (hash >> 32) as u128;
        stats.sum += hash_scaled;
        stats.sum_sq += hash_scaled * hash_scaled;

        stats.min_hash = stats.min_hash.min(hash);
        stats.max_hash = stats.max_hash.max(hash);
    }

    Ok(stats)
}

fn chi_squared_test(observed: &[u64], expected: f64) -> (f64, f64) {
    let chi_sq: f64 = observed
        .iter()
        .map(|&o| {
            let diff = o as f64 - expected;
            (diff * diff) / expected
        })
        .sum();

    let df = observed.len() - 1;
    let p_value = chi_squared_p_value(chi_sq, df);

    (chi_sq, p_value)
}

/// Approximate chi-squared p-value using Wilson-Hilferty transformation
fn chi_squared_p_value(chi_sq: f64, df: usize) -> f64 {
    let k = df as f64;
    let z = ((chi_sq / k).powf(1.0 / 3.0) - (1.0 - 2.0 / (9.0 * k))) / (2.0 / (9.0 * k)).sqrt();
    0.5 * (1.0 + libm::erf(-z / std::f64::consts::SQRT_2))
}

fn coefficient_of_variation(counts: &[u64]) -> f64 {
    let n = counts.len() as f64;
    let mean: f64 = counts.iter().sum::<u64>() as f64 / n;
    let variance: f64 = counts.iter().map(|&x| {
        let diff = x as f64 - mean;
        diff * diff
    }).sum::<f64>() / n;
    let stddev = variance.sqrt();
    stddev / mean
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    if let Some(threads) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()?;
    }

    println!("Analyzing hash uniformity...");
    println!("Input directory: {}", args.input_dir);
    println!("Hash name: {}", args.hash_name);
    println!("Number of output bins: {}", args.num_bins);
    println!();

    // Collect all bin files
    let mut bin_files: Vec<String> = Vec::new();
    for bindex in 0..16 {
        for bin_id in 0..128 {
            let input_file = format!(
                "{}/unsorted_bin_{}_{}_{}.bin",
                args.input_dir, args.hash_name, bindex, bin_id
            );
            if Path::new(&input_file).exists() {
                bin_files.push(input_file);
            }
        }
    }

    println!("Found {} bin files to process", bin_files.len());

    // Progress counter
    let processed = AtomicU64::new(0);
    let total_files = bin_files.len() as u64;

    // Process all files in parallel
    let results: Vec<_> = bin_files
        .par_iter()
        .map(|file| {
            let result = process_bin_file(file, args.num_bins);
            let count = processed.fetch_add(1, Ordering::Relaxed) + 1;
            if count % 100 == 0 || count == total_files {
                eprintln!("Processed {}/{} files", count, total_files);
            }
            result
        })
        .collect();

    // Merge all results
    let mut combined = BinStats::new(args.num_bins);
    let mut errors = 0;
    for result in results {
        match result {
            Ok(stats) => combined.merge(&stats),
            Err(e) => {
                eprintln!("Error: {}", e);
                errors += 1;
            }
        }
    }

    if errors > 0 {
        eprintln!("\nWarning: {} files had errors", errors);
    }

    // Print results
    println!("\n{}", "=".repeat(60));
    println!("HASH UNIFORMITY ANALYSIS");
    println!("{}", "=".repeat(60));
    println!();

    println!("Total hashes analyzed:               {:>15}", format_number(combined.total_hashes));
    println!("Min hash value:                      0x{:016x}", combined.min_hash);
    println!("Max hash value:                      0x{:016x}", combined.max_hash);
    println!();

    // Bin distribution statistics
    let expected_per_bin = combined.total_hashes as f64 / args.num_bins as f64;
    let (chi_sq, p_value) = chi_squared_test(&combined.bin_counts, expected_per_bin);
    let cv = coefficient_of_variation(&combined.bin_counts);

    let min_bin_count = *combined.bin_counts.iter().min().unwrap_or(&0);
    let max_bin_count = *combined.bin_counts.iter().max().unwrap_or(&0);
    let empty_bins = combined.bin_counts.iter().filter(|&&c| c == 0).count();
    let populated_bins = args.num_bins - empty_bins;

    println!("--- Bin Distribution (mod {}) ---", args.num_bins);
    println!("Expected per bin:                    {:>15.2}", expected_per_bin);
    println!("Min bin count:                       {:>15}", format_number(min_bin_count));
    println!("Max bin count:                       {:>15}", format_number(max_bin_count));
    println!("Empty bins:                          {:>15}", empty_bins);
    println!("Populated bins:                      {:>15}", populated_bins);
    println!("Coverage:                            {:>14.2}%", (populated_bins as f64 / args.num_bins as f64) * 100.0);
    println!("Coefficient of variation:            {:>15.4}", cv);
    println!();

    println!("--- Uniformity Test (Chi-squared) ---");
    println!("Chi-squared statistic:               {:>15.2}", chi_sq);
    println!("Degrees of freedom:                  {:>15}", args.num_bins - 1);
    println!("P-value:                             {:>15.4e}", p_value);

    let interpretation = if p_value < 0.01 {
        "NOT UNIFORM (p < 0.01)"
    } else if p_value < 0.05 {
        "MARGINAL (0.01 ≤ p < 0.05)"
    } else {
        "UNIFORM (p ≥ 0.05)"
    };
    println!("Interpretation:                      {}", interpretation);
    println!();

    // Show bins with extreme counts
    println!("--- Extreme Bins ---");
    let mean_count = combined.total_hashes as f64 / args.num_bins as f64;
    let threshold_low = mean_count * 0.9;
    let threshold_high = mean_count * 1.1;

    let mut extreme_bins: Vec<(usize, u64)> = combined.bin_counts
        .iter()
        .enumerate()
        .filter(|(_, c)| (**c as f64) < threshold_low || (**c as f64) > threshold_high)
        .map(|(i, c)| (i, *c))
        .collect();

    extreme_bins.sort_by(|a, b| a.1.cmp(&b.1));

    if extreme_bins.is_empty() {
        println!("No bins deviate more than 10% from expected.");
    } else {
        println!("Bins deviating >10% from expected (showing up to 20):");
        for (i, (bin_idx, count)) in extreme_bins.iter().take(10).enumerate() {
            let deviation = (*count as f64 - mean_count) / mean_count * 100.0;
            println!("  Bin {:>4}: {:>12} ({:+.2}%)", bin_idx, format_number(*count), deviation);
            if i == 9 && extreme_bins.len() > 20 {
                println!("  ... ({} more low bins)", extreme_bins.len() / 2 - 10);
            }
        }
        if extreme_bins.len() > 10 {
            println!();
            for (bin_idx, count) in extreme_bins.iter().rev().take(10).rev() {
                let deviation = (*count as f64 - mean_count) / mean_count * 100.0;
                println!("  Bin {:>4}: {:>12} ({:+.2}%)", bin_idx, format_number(*count), deviation);
            }
        }
    }

    println!();
    println!("{}", "=".repeat(60));

    // MinHash validity conclusion
    if p_value >= 0.05 && cv < 0.2 {
        println!("✓ Hash outputs are uniformly distributed.");
        println!("✓ MinHash uniformity requirement is SATISFIED.");
    } else {
        println!("⚠ Hash outputs may not be uniformly distributed.");
        println!("⚠ MinHash uniformity requirement needs investigation.");
    }

    Ok(())
}

fn format_number(n: u64) -> String {
    let s = n.to_string();
    let mut result = String::new();
    for (i, c) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            result.push(',');
        }
        result.push(c);
    }
    result.chars().rev().collect()
}
