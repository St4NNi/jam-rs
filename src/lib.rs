pub mod cli;
pub mod core_utils;
pub mod distance;
pub mod sketch;
pub mod stats;
pub mod writer;

// Re-export jamhash for backwards compatibility
pub use jamhash::jamhash_u64;

use anyhow::Result;
use clap::Parser;
use cli::{Cli, Commands};
use distance::{DistanceConfig, OutputFormat};
use sketch::{SketchConfig, sketch_files};
use stats::StatsCalculator;
use std::{
    fs::remove_file,
    path::{Path, PathBuf},
};

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
            complexity,
            singleton,
            temp_dir,
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
                cli.memory.unwrap_or(2),
                cli.force,
                cli.silent,
                complexity,
                temp_dir,
            )
        }

        Commands::Dist {
            input,
            database,
            output,
            cutoff,
            singleton,
        } => handle_distance_command(input, database, output, cutoff, singleton, cli.silent),

        Commands::Stats { input, short, full } => {
            handle_stats_command(input, short, full, cli.silent)
        }
    }
}

/// Handle the sketch command
/// This is the main input for sketching where it is ok to have many arguments.
#[allow(clippy::too_many_arguments)]
pub fn handle_sketch_command(
    input_paths: Vec<PathBuf>,
    output_path: PathBuf,
    kmer_size: u8,
    fscale: Option<u64>,
    nmax: Option<u64>,
    singleton: bool,
    threads: usize,
    memory: usize,
    force: bool,
    silent: bool,
    min_entropy: f64,
    temp_dir: Option<PathBuf>,
) -> Result<()> {
    // Validate temp directory if provided
    if let Some(ref temp_dir) = temp_dir {
        if !temp_dir.exists() {
            return Err(anyhow::anyhow!(
                "Temp directory does not exist: {:?}",
                temp_dir
            ));
        }
        if !temp_dir.is_dir() {
            return Err(anyhow::anyhow!(
                "Temp directory path is not a directory: {:?}",
                temp_dir
            ));
        }
    }

    // Check if output file exists and force flag
    if output_path.exists() {
        if !force {
            return Err(anyhow::anyhow!(
                "Output file {:?} already exists. Use --force to overwrite.",
                output_path
            ));
        }
        if !silent {
            eprintln!(
                "Warning: Overwriting existing output file: {}",
                output_path.display()
            );
        }
        if !output_path.is_file() {
            return Err(anyhow::anyhow!(
                "Output path must be a file, not a directory: {:?}",
                output_path
            ));
        }

        remove_file(&output_path)?;
    }

    // Validate k-mer size
    if kmer_size == 0 || kmer_size >= 32 {
        return Err(anyhow::anyhow!(
            "K-mer size must be between 1 and 31, got {}",
            kmer_size
        ));
    }

    if !silent {
        let mut settings = format!(
            "jam: {} files, k={}, threads={}, memory={}GB, entropy={}",
            input_paths.len(),
            kmer_size,
            threads,
            memory,
            min_entropy
        );

        if let Some(scale) = fscale {
            settings.push_str(&format!(", scale={}", scale));
        }

        if let Some(max_hashes) = nmax {
            settings.push_str(&format!(", nmax={}", max_hashes));
        }

        if singleton {
            settings.push_str(", mode=singleton");
        } else {
            settings.push_str(", mode=combined");
        }

        eprintln!("{}", settings);
    }

    let fscale = if let Some(fscale) = fscale {
        (u64::MAX as f64 / fscale as f64) as u64
    } else {
        u64::MAX
    };
    let nmax = nmax.unwrap_or(u64::MAX);

    let config = SketchConfig {
        kmer_size,
        fscale,
        nmax,
        singleton,
        min_entropy, // Could be made configurable
        threads,
        memory_budget_gb: memory as f64, // Could be made configurable
        temp_dir,
    };

    sketch_files(&input_paths, output_path.clone(), config, silent)?;

    if !silent {
        let mut completion_msg = format!("Completed: {}", output_path.display());

        // Add quick stats on the same line
        if let Ok(stats) = stats::StatsCalculator::calculate_stats(&output_path) {
            completion_msg.push_str(&format!(
                " ({} hashes, {} files/sequences)",
                stats.total_hashes, stats.unique_files
            ));
        }

        eprintln!("{}", completion_msg);
    }

    Ok(())
}

/// Handle the distance command  
pub fn handle_distance_command(
    input_path: PathBuf,
    database_path: PathBuf,
    output_path: Option<PathBuf>,
    cutoff: f64,
    singleton: bool,
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
        let output_desc = if let Some(ref output) = output_path {
            output.display().to_string()
        } else {
            "stdout".to_string()
        };

        eprintln!(
            "Distance calculation: {} vs {} (cutoff={}) -> {}",
            input_path.display(),
            database_path.display(),
            cutoff,
            output_desc
        );
    }

    let distance_config = DistanceConfig {
        cutoff,
        length_category_mode: distance::LengthCategoryMode::QueryAndBelow,
        output_format: OutputFormat::Tsv,
    };

    let results = calculate_distances_streaming(
        &input_path,
        &database_path,
        output_path.as_deref(),
        distance_config,
        singleton,
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
pub fn handle_stats_command(
    input_path: PathBuf,
    short: bool,
    full: bool,
    silent: bool,
) -> Result<()> {
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

    StatsCalculator::print_stats(&input_path, short, full)?;

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
pub fn is_sequence_file(path: &Path) -> bool {
    if let Some(ext) = path.extension().map(|e| e.to_string_lossy().to_lowercase()) {
        // Check for .gz and compound extensions
        if ext == "gz"
            && let Some(stem_ext) = path.file_stem().and_then(|s| Path::new(s).extension())
        {
            let stem_ext = stem_ext.to_string_lossy().to_lowercase();
            return matches!(
                stem_ext.as_str(),
                "fasta" | "fa" | "fas" | "fna" | "fastq" | "fq"
            );
        }
        return matches!(
            ext.as_str(),
            "fasta" | "fa" | "fas" | "fna" | "fastq" | "fq"
        );
    }
    false
}
