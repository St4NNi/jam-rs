pub mod cli;
pub mod core_utils;
pub mod distance;
pub mod hash_functions;
pub mod sketch;
pub mod stats;
pub mod writer;

use anyhow::Result;
use clap::Parser;
use cli::{Cli, Commands};
use distance::{DistanceConfig, OutputFormat};
use sketch::{SketchConfig, sketch_files};
use stats::StatsCalculator;
use std::path::{Path, PathBuf};

use crate::distance::calculate_distances_streaming;

/// Main application entry point
pub fn run() -> Result<()> {
    let cli = Cli::parse();

    // Set up rayon thread pool if threads specified
    if let Some(threads) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()?;
    }

    match cli.command {
        Commands::Sketch {
            input,
            output,
            kmer_size,
            fscale,
            nmax,
            singleton,
        } => {
            // Expand input paths to handle directories and lists
            let expanded_inputs = expand_input_paths(&input)?;

            handle_sketch_command(
                expanded_inputs,
                output,
                kmer_size,
                fscale,
                nmax,
                singleton,
                cli.threads.unwrap_or(1),
                cli.force,
                cli.silent,
            )
        }

        Commands::Dist {
            input,
            database,
            output,
            cutoff,
        } => handle_distance_command(input, database, output, cutoff, cli.silent),

        Commands::Stats { input, short } => handle_stats_command(input, short, cli.silent),
    }
}

/// Handle the sketch command
/// This is the main input for sketching where it is ok to have many arguments.
#[allow(clippy::too_many_arguments)]
fn handle_sketch_command(
    input_paths: Vec<PathBuf>,
    output_path: PathBuf,
    kmer_size: u8,
    fscale: Option<u64>,
    nmax: Option<u64>,
    singleton: bool,
    threads: usize,
    force: bool,
    silent: bool,
) -> Result<()> {
    // Ensure output has .lmdb extension
    let output_path = if output_path.extension().is_none() {
        output_path.with_extension("lmdb")
    } else {
        output_path
    };

    // Check if output file exists and force flag
    if output_path.exists() && !force {
        return Err(anyhow::anyhow!(
            "Output file {:?} already exists. Use --force to overwrite.",
            output_path
        ));
    }

    // Validate k-mer size
    if kmer_size == 0 || kmer_size >= 32 {
        return Err(anyhow::anyhow!(
            "K-mer size must be between 1 and 31, got {}",
            kmer_size
        ));
    }

    if !silent {
        println!("Starting sketching process...");
        println!("  Input files: {}", input_paths.len());
        println!("  K-mer size: {kmer_size}");
        println!("  Threads: {threads}");

        if let Some(scale) = fscale {
            println!("  FracMinHash scale: {scale}");
        }

        if let Some(max_hashes) = nmax {
            println!("  Max hashes per sequence: {max_hashes}");
        }

        if singleton {
            println!("  Mode: Singleton (separate sketch per sequence)");
        } else {
            println!("  Mode: Combined (one sketch per file)");
        }
    }

    let config = SketchConfig {
        kmer_size,
        fscale,
        nmax,
        singleton,
        min_entropy: 1.5, // Could be made configurable
        threads,
        memory_budget_gb: 2.0, // Could be made configurable
    };

    sketch_files(&input_paths, output_path.clone(), config)?;

    if !silent {
        println!("Sketching completed successfully!");
        println!("Output written to: {output_path:?}");

        // Print quick stats
        if let Ok(stats) = stats::StatsCalculator::calculate_stats(&output_path, false) {
            println!("  Total hashes: {}", stats.total_hashes);
            println!("  Files/sequences processed: {}", stats.unique_files);
        }
    }

    Ok(())
}

/// Handle the distance command  
fn handle_distance_command(
    input_path: PathBuf,
    database_path: PathBuf,
    output_path: Option<PathBuf>,
    cutoff: f64,
    silent: bool,
) -> Result<()> {
    // Validate cutoff
    if !(0.0..=1.0).contains(&cutoff) {
        return Err(anyhow::anyhow!(
            "Cutoff must be between 0.0 and 1.0, got {}",
            cutoff
        ));
    }

    // Validate input paths
    if !input_path.exists() {
        return Err(anyhow::anyhow!(
            "Input path does not exist: {:?}",
            input_path
        ));
    }

    if !database_path.exists() {
        return Err(anyhow::anyhow!(
            "Database path does not exist: {:?}",
            database_path
        ));
    }

    if !silent {
        println!("Starting distance calculation...");
        println!("  Input: {input_path:?}");
        println!("  Database: {database_path:?}");
        println!("  Cutoff: {cutoff}");

        if let Some(ref output) = output_path {
            println!("  Output: {output:?}");
        } else {
            println!("  Output: stdout");
        }
    }

    let distance_config = DistanceConfig {
        cutoff,
        length_category_mode: distance::LengthCategoryMode::QueryAndBelow,
        output_format: OutputFormat::Tsv,
    };

    // Default sketch config for on-the-fly sketching if needed
    let sketch_config = SketchConfig {
        kmer_size: 21, // Should ideally match database k-mer size
        fscale: None,
        nmax: None,
        singleton: true,
        min_entropy: 1.5,
        threads: 1,
        memory_budget_gb: 1.0,
    };

    let results = calculate_distances_streaming(
        &input_path,
        &database_path,
        output_path.as_deref(),
        distance_config,
        sketch_config,
    )?;

    if !silent {
        println!("Distance calculation completed!");
        println!("Found {} matches above cutoff threshold", results.len());

        if !results.is_empty() {
            let best_match = &results[0];
            println!(
                "Best match: {} (containment: {:.3}, Jaccard: {:.3})",
                best_match.target_name,
                best_match.containment_query_in_target,
                best_match.jaccard_similarity
            );
        }
    }

    Ok(())
}

/// Handle the stats command
fn handle_stats_command(input_path: PathBuf, short: bool, silent: bool) -> Result<()> {
    if !input_path.exists() {
        return Err(anyhow::anyhow!(
            "Input path does not exist: {:?}",
            input_path
        ));
    }

    if !silent && !short {
        println!("Calculating database statistics...");
        println!("Database: {input_path:?}");
        println!();
    }

    StatsCalculator::print_stats(&input_path, short)?;

    Ok(())
}

/// Expand input paths to handle directories and file lists
pub fn expand_input_paths(input_paths: &[PathBuf]) -> Result<Vec<PathBuf>> {
    let mut expanded_paths = Vec::new();

    for path in input_paths {
        if path.is_dir() {
            // Add all FASTA/FASTQ files in directory
            for entry in std::fs::read_dir(path)? {
                let entry = entry?;
                let file_path = entry.path();

                if is_sequence_file(&file_path) {
                    expanded_paths.push(file_path);
                }
            }
        } else if path.is_file() {
            // Check if it's a list file or sequence file
            if is_sequence_file(path) {
                expanded_paths.push(path.clone());
            } else {
                // Assume it's a list of files
                let content = std::fs::read_to_string(path)?;
                for line in content.lines() {
                    let file_path = PathBuf::from(line.trim());
                    if file_path.exists() && is_sequence_file(&file_path) {
                        expanded_paths.push(file_path);
                    }
                }
            }
        }
    }

    if expanded_paths.is_empty() {
        return Err(anyhow::anyhow!(
            "No valid sequence files found in input paths"
        ));
    }

    // Sort for consistent processing order
    expanded_paths.sort();

    Ok(expanded_paths)
}

/// Check if a file is a sequence file based on extension
fn is_sequence_file(path: &Path) -> bool {
    if let Some(extension) = path.extension() {
        let ext = extension.to_string_lossy().to_lowercase();
        matches!(
            ext.as_str(),
            "fasta" | "fa" | "fas" | "fna" | "fastq" | "fq" | "fastq.gz" | "fq.gz" | "fa.gz"
        )
    } else {
        false
    }
}
