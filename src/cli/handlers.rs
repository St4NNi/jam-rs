use anyhow::Result;
use std::{fs::remove_file, path::PathBuf};

use crate::bias::BiasTable;
use crate::distance::{calculate_distances_streaming, DistanceConfig, LengthCategoryMode, OutputFormat};
use crate::sketch::{sketch_files, SketchConfig};
use crate::stats::StatsCalculator;

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
    bias_table_path: Option<PathBuf>,
) -> Result<()> {
    if let Some(ref temp_dir) = temp_dir {
        if !temp_dir.exists() {
            return Err(anyhow::anyhow!("Temp directory does not exist: {:?}", temp_dir));
        }
        if !temp_dir.is_dir() {
            return Err(anyhow::anyhow!("Temp directory path is not a directory: {:?}", temp_dir));
        }
    }

    if output_path.exists() {
        if !force {
            return Err(anyhow::anyhow!(
                "Output file {:?} already exists. Use --force to overwrite.",
                output_path
            ));
        }
        if !silent {
            eprintln!("Warning: Overwriting existing output file: {}", output_path.display());
        }
        if !output_path.is_file() {
            return Err(anyhow::anyhow!("Output path must be a file, not a directory: {:?}", output_path));
        }
        remove_file(&output_path)?;
    }

    if kmer_size == 0 || kmer_size >= 32 {
        return Err(anyhow::anyhow!("K-mer size must be between 1 and 31, got {}", kmer_size));
    }

    if !silent {
        let mut settings = format!(
            "jam: {} files, k={}, threads={}, memory={}GB, entropy={}",
            input_paths.len(), kmer_size, threads, memory, min_entropy
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

    let bias_table = if let Some(ref path) = bias_table_path {
        if !path.exists() {
            return Err(anyhow::anyhow!("Bias table file does not exist: {:?}", path));
        }
        let table = BiasTable::load(path)?;
        if !silent {
            table.print_stats();
        }
        Some(table)
    } else {
        None
    };

    let fscale = fscale.map(|f| (u64::MAX as f64 / f as f64) as u64).unwrap_or(u64::MAX);
    let nmax = nmax.unwrap_or(u64::MAX);

    let config = SketchConfig {
        kmer_size,
        fscale,
        nmax,
        singleton,
        min_entropy,
        threads,
        memory_budget_gb: memory as f64,
        temp_dir,
        bias_table,
    };

    sketch_files(&input_paths, output_path.clone(), config, silent)?;

    if !silent {
        let mut completion_msg = format!("Completed: {}", output_path.display());
        if let Ok(stats) = StatsCalculator::calculate_stats(&output_path) {
            completion_msg.push_str(&format!(" ({} hashes, {} files/sequences)", stats.total_hashes, stats.unique_files));
        }
        eprintln!("{}", completion_msg);
    }

    Ok(())
}

pub fn handle_distance_command(
    input_path: PathBuf,
    database_path: PathBuf,
    output_path: Option<PathBuf>,
    cutoff: f64,
    singleton: bool,
    silent: bool,
    bias_table_path: Option<PathBuf>,
) -> Result<()> {
    if !(0.0..=1.0).contains(&cutoff) {
        return Err(anyhow::anyhow!("Cutoff must be between 0.0 and 1.0, got {}", cutoff));
    }

    if !input_path.exists() {
        return Err(anyhow::anyhow!("Input path does not exist: {:?}", input_path));
    }

    if !database_path.exists() {
        return Err(anyhow::anyhow!("Database path does not exist: {:?}", database_path));
    }

    let bias_table = if let Some(ref path) = bias_table_path {
        if !path.exists() {
            return Err(anyhow::anyhow!("Bias table file does not exist: {:?}", path));
        }
        Some(BiasTable::load(path)?)
    } else {
        None
    };

    if !silent {
        let output_desc = output_path.as_ref().map(|p| p.display().to_string()).unwrap_or_else(|| "stdout".to_string());
        eprintln!("Distance calculation: {} vs {} (cutoff={}) -> {}", input_path.display(), database_path.display(), cutoff, output_desc);
    }

    let distance_config = DistanceConfig {
        cutoff,
        length_category_mode: LengthCategoryMode::QueryAndBelow,
        output_format: OutputFormat::Tsv,
    };

    let results = calculate_distances_streaming(&input_path, &database_path, output_path.as_deref(), distance_config, singleton, bias_table)?;

    if !silent {
        println!("Distance calculation completed!");
        println!("Found {} matches above cutoff threshold", results.len());
        if !results.is_empty() {
            let best_match = &results[0];
            println!("Best match: {} (containment: {:.3}, Jaccard: {:.3})", best_match.target_name, best_match.containment_query_in_target, best_match.jaccard_similarity);
        }
    }

    Ok(())
}

pub fn handle_bias_command(
    positive: PathBuf,
    negative: PathBuf,
    output: PathBuf,
    threshold: f32,
    force: bool,
    silent: bool,
) -> Result<()> {
    if !positive.exists() {
        return Err(anyhow::anyhow!("Positive file does not exist: {:?}", positive));
    }
    if !negative.exists() {
        return Err(anyhow::anyhow!("Negative file does not exist: {:?}", negative));
    }
    if output.exists() && !force {
        return Err(anyhow::anyhow!(
            "Output file {:?} already exists. Use --force to overwrite.",
            output
        ));
    }

    if !silent {
        eprintln!(
            "Building bias table: positive={}, negative={}, threshold={:.2}",
            positive.display(),
            negative.display(),
            threshold
        );
    }

    let table = BiasTable::build(&positive, &negative, threshold)?;

    if !silent {
        table.print_stats();
    }

    table.save(&output)?;

    if !silent {
        eprintln!("Bias table saved to: {}", output.display());
    }

    Ok(())
}

pub fn handle_stats_command(input_path: PathBuf, short: bool, full: bool, silent: bool) -> Result<()> {
    if !input_path.exists() {
        return Err(anyhow::anyhow!("Input path does not exist: {:?}", input_path));
    }

    if !silent && !short {
        println!("Calculating database statistics...");
        println!("Database: {input_path:?}");
        println!();
    }

    StatsCalculator::print_stats(&input_path, short, full)?;
    Ok(())
}
