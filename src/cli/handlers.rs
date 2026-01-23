use anyhow::Result;
use std::{fs::remove_file, path::PathBuf};

use crate::bias::BiasTable;
use crate::sketch::{self, SketchConfig};
use std::sync::Arc;

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

    let config = SketchConfig {
        kmer_size,
        fscale: fscale.unwrap_or(1000),
        num_threads: threads,
        memory,
        singleton,
        min_entropy,
        temp_dir_base: temp_dir,
        bias_table: bias_table.map(Arc::new),
        ..Default::default()
    };

    let result = sketch::run(&input_paths, &config)?;

    if !silent {
        let total: u64 = result.bucket_entry_counts.iter().sum();
        eprintln!(
            "Completed: {} ({} entries, {} samples)",
            output_path.display(),
            total,
            result.sample_count
        );
    }

    Ok(())
}

pub fn handle_distance_command(
    _input_path: PathBuf,
    _database_path: PathBuf,
    _output_path: Option<PathBuf>,
    _cutoff: f64,
    _singleton: bool,
    _silent: bool,
    _bias_table_path: Option<PathBuf>,
) -> Result<()> {
    unimplemented!("distance command pending storage rework")
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

pub fn handle_stats_command(
    _input_path: PathBuf,
    _short: bool,
    _full: bool,
    _silent: bool,
) -> Result<()> {
    unimplemented!("stats command pending storage rework")
}
