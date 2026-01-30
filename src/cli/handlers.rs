use anyhow::Result;
use indicatif::{ProgressBar, ProgressStyle};
use std::io::Write;
use std::{fs::remove_file, path::PathBuf};

use crate::bias::BiasTable;
use crate::query::QueryEngine;
use crate::reader::JamReader;
use crate::writer::{BuildConfig, build};
use std::sync::Arc;

#[allow(clippy::too_many_arguments)]
pub fn handle_sketch_command(
    input_paths: Vec<PathBuf>,
    output_path: PathBuf,
    kmer_size: u8,
    fscale: Option<u64>,
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
        if singleton {
            settings.push_str(", mode=singleton");
        } else {
            settings.push_str(", mode=combined");
        }
        eprintln!("{}", settings);
    }

    let bias_table = if let Some(ref path) = bias_table_path {
        if !path.exists() {
            return Err(anyhow::anyhow!(
                "Bias table file does not exist: {:?}",
                path
            ));
        }
        let table = BiasTable::load(path)?;
        if !silent {
            table.print_stats();
        }
        Some(table)
    } else {
        None
    };

    let config = BuildConfig {
        kmer_size,
        fscale: fscale.unwrap_or(1000),
        num_threads: threads,
        memory,
        singleton,
        min_entropy,
        temp_dir_base: temp_dir,
        bias_table: bias_table.map(Arc::new),
        show_progress: !silent,
    };

    let stats = build(&input_paths, &output_path, &config)?;

    if !silent {
        eprintln!(
            "Completed: {} ({} entries, {} unique hashes, {} samples)",
            output_path.display(),
            stats.total_entries,
            stats.unique_hashes,
            stats.sample_count
        );
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
) -> Result<()> {
    if !database_path.exists() {
        return Err(anyhow::anyhow!(
            "Database file does not exist: {:?}",
            database_path
        ));
    }
    if !input_path.exists() {
        return Err(anyhow::anyhow!(
            "Input file does not exist: {:?}",
            input_path
        ));
    }

    let engine = QueryEngine::open(&database_path)?;
    let db_names = engine.reader().sample_names();
    let db_sizes = engine.reader().sample_sizes();

    if !silent && engine.has_bias_table() {
        eprintln!("Using embedded bias table from database");
    }

    // Progress spinner
    let spinner = if !silent {
        let sp = ProgressBar::new_spinner();
        sp.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} [{elapsed_precise}] {msg}")
                .unwrap(),
        );
        sp.set_message("Loading query...");
        sp.enable_steady_tick(std::time::Duration::from_millis(80));
        Some(sp)
    } else {
        None
    };

    // Use QuerySketch for all inputs (handles FASTA, FASTQ, and JAM)
    // Bias filtering uses the embedded bias table from the database
    let sketch = crate::query::QuerySketch::from_inputs(
        std::slice::from_ref(&input_path),
        engine.reader(),
        singleton,
    )
    .map_err(|e| anyhow::anyhow!("{}", e))?;

    if sketch.sample_count() == 0 {
        if let Some(sp) = spinner {
            sp.finish_and_clear();
        }
        if !silent {
            eprintln!("No sequences found in input");
        }
        return Ok(());
    }

    // Update spinner with query info
    if let Some(ref sp) = spinner {
        let db_samples = engine.reader().stats().sample_count;
        sp.set_message(format!(
            "Searching {} query hashes against {} db samples...",
            sketch.total_entries(),
            db_samples
        ));
    }

    let results = engine.query_sketch(&sketch);

    // Finish spinner
    if let Some(sp) = spinner {
        sp.finish_and_clear();
    }

    // Output
    let mut writer: Box<dyn Write> = if let Some(ref out) = output_path {
        Box::new(std::fs::File::create(out)?)
    } else {
        Box::new(std::io::stdout())
    };

    writeln!(
        writer,
        "query\tdb_sample\tshared_hashes\tquery_hashes\tdb_hashes\tquery_containment\tdb_containment"
    )?;

    for (query_idx, result) in results.iter().enumerate() {
        let query_name = &sketch.sample_names[query_idx];
        let query_hashes = result.query_size;

        // Sort by containment descending
        let mut matches: Vec<_> = result
            .matches
            .iter()
            .filter(|m| m.containment >= cutoff)
            .collect();
        matches.sort_by(|a, b| b.containment.total_cmp(&a.containment));

        for m in matches {
            let db_name = db_names
                .get(m.sample_id as usize)
                .map(|s| s.as_str())
                .unwrap_or("unknown");
            let db_hashes = db_sizes.get(m.sample_id as usize).copied().unwrap_or(0);
            let shared_hashes = m.hit_count;
            let query_containment = m.containment;
            let db_containment = if db_hashes > 0 {
                shared_hashes as f64 / db_hashes as f64
            } else {
                0.0
            };

            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}",
                query_name,
                db_name,
                shared_hashes,
                query_hashes,
                db_hashes,
                query_containment,
                db_containment
            )?;
        }
    }

    if !silent {
        eprintln!("Query complete: {}", input_path.display());
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
        return Err(anyhow::anyhow!(
            "Positive file does not exist: {:?}",
            positive
        ));
    }
    if !negative.exists() {
        return Err(anyhow::anyhow!(
            "Negative file does not exist: {:?}",
            negative
        ));
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
    input_path: PathBuf,
    short: bool,
    full: bool,
    silent: bool,
) -> Result<()> {
    if !input_path.exists() {
        return Err(anyhow::anyhow!(
            "Database file does not exist: {:?}",
            input_path
        ));
    }

    let reader = JamReader::open(&input_path)?;
    let stats = reader.stats();

    if short {
        println!(
            "{}\t{}\t{}\t{}",
            stats.entry_count, stats.unique_hash_count, stats.sample_count, stats.kmer_size
        );
        return Ok(());
    }

    if !silent {
        println!("JAM Database Statistics");
        println!("=======================");
        println!("File: {}", input_path.display());
        println!("File size: {} bytes", stats.file_size);
        println!();
        println!("K-mer size: {}", stats.kmer_size);
        println!("Hash threshold: {}", stats.hash_threshold);
        println!("Sample rate: 1/{}", u64::MAX / stats.hash_threshold.max(1));
        println!(
            "Embedded bias table: {}",
            if stats.has_bias_table { "yes" } else { "no" }
        );
        println!();
        println!("Total entries: {}", stats.entry_count);
        println!("Unique hashes: {}", stats.unique_hash_count);
        println!("Sample count: {}", stats.sample_count);
    }

    if full {
        println!();
        println!("Per-Bucket Statistics");
        println!("---------------------");
        println!("bucket\tentries");
        for (i, &count) in stats.bucket_entry_counts.iter().enumerate() {
            if count > 0 {
                println!("{}\t{}", i, count);
            }
        }

        let non_empty = stats.bucket_entry_counts.iter().filter(|&&c| c > 0).count();
        let avg = if non_empty > 0 {
            stats.entry_count as f64 / non_empty as f64
        } else {
            0.0
        };
        println!();
        println!("Non-empty buckets: {}/256", non_empty);
        println!("Average entries per non-empty bucket: {:.1}", avg);
    }

    Ok(())
}
