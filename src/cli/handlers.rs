use anyhow::Result;
use needletail::Sequence;
use std::collections::HashSet;
use std::io::Write;
use std::{fs::remove_file, path::PathBuf};

use crate::bias::BiasTable;
use crate::core_utils::passes_entropy_filter;
use crate::query::QueryEngine;
use crate::reader::JamReader;
use crate::writer::{build, BuildConfig};
use jamhash::jamhash_u64;
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
    bias_table_path: Option<PathBuf>,
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
    let threshold = engine.threshold();
    let k = engine.kmer_size();

    // Determine which bias table to use:
    // 1. If --bias-table CLI arg provided, use that (with warning if different from embedded)
    // 2. Otherwise use embedded bias table if present
    // 3. Otherwise no bias filtering
    let (bias_table, bias_source) = if let Some(ref path) = bias_table_path {
        let external_bias = BiasTable::load(path)?;

        // Check if embedded bias exists and differs from external
        if let Some(ref embedded) = engine.bias_table() {
            if embedded.as_ref() != &external_bias {
                eprintln!(
                    "Warning: External bias table differs from embedded bias in database. Using external."
                );
            }
        }

        (Some(external_bias), "external")
    } else if let Some(embedded) = engine.bias_table() {
        // Use Arc's inner value - clone the BiasTable
        (Some((*embedded).clone()), "embedded")
    } else {
        (None, "none")
    };

    if !silent && bias_table.is_some() {
        eprintln!("Using bias table: {}", bias_source);
    }

    let mut writer: Box<dyn Write> = if let Some(ref out) = output_path {
        Box::new(std::fs::File::create(out)?)
    } else {
        Box::new(std::io::stdout())
    };

    writeln!(writer, "query\tsample_id\thit_count\tcontainment")?;

    let mut reader = match needletail::parse_fastx_file(&input_path) {
        Ok(reader) => reader,
        Err(e) if e.kind == needletail::errors::ParseErrorKind::EmptyFile => {
            eprintln!("Empty file detected: {}, skipping", input_path.display());
            return Ok(());
        }
        Err(e) => return Err(anyhow::anyhow!("{}", e)),
    };

    if singleton {
        while let Some(record) = reader.next() {
            let record = record.map_err(|e| anyhow::anyhow!("{}", e))?;
            let query_name = String::from_utf8_lossy(record.id()).to_string();
            let hashes = extract_hashes(&record.seq(), k, threshold, 0.0, bias_table.as_ref());

            if hashes.is_empty() {
                continue;
            }

            let result = engine.query(&hashes);
            for m in result.matches.iter().filter(|m| m.containment >= cutoff) {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{:.6}",
                    query_name, m.sample_id, m.hit_count, m.containment
                )?;
            }
        }
    } else {
        let mut all_hashes = HashSet::new();
        while let Some(record) = reader.next() {
            let record = record.map_err(|e| anyhow::anyhow!("{}", e))?;
            let hashes = extract_hashes(&record.seq(), k, threshold, 0.0, bias_table.as_ref());
            all_hashes.extend(hashes);
        }

        let hashes: Vec<u64> = all_hashes.into_iter().collect();
        if !hashes.is_empty() {
            let result = engine.query(&hashes);
            for m in result.matches.iter().filter(|m| m.containment >= cutoff) {
                writeln!(
                    writer,
                    "combined\t{}\t{}\t{:.6}",
                    m.sample_id, m.hit_count, m.containment
                )?;
            }
        }
    }

    if !silent {
        eprintln!("Query complete: {}", input_path.display());
    }

    Ok(())
}

fn extract_hashes(
    sequence: &[u8],
    k: u8,
    threshold: u64,
    min_entropy: f64,
    bias_table: Option<&BiasTable>,
) -> Vec<u64> {
    let mut hashes = HashSet::new();
    let norm = needletail::Sequence::normalize(sequence, false);

    if norm.len() < k as usize {
        return Vec::new();
    }

    for (_, kmer, _) in norm.bit_kmers(k, true) {
        let hash = jamhash_u64(kmer.0);

        if hash >= threshold {
            continue;
        }

        if min_entropy > 0.0 && !passes_entropy_filter(kmer.0, k, min_entropy) {
            continue;
        }

        if bias_table.is_some_and(|b| !b.passes_filter(kmer.0, k)) {
            continue;
        }

        hashes.insert(hash);
    }

    hashes.into_iter().collect()
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
        println!(
            "Sample rate: 1/{}",
            u64::MAX / stats.hash_threshold.max(1)
        );
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
