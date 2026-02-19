use anyhow::Result;
use indicatif::{ProgressBar, ProgressStyle};
use std::io::Write;
use std::{fs::remove_file, path::PathBuf};

use crate::bias::{BiasCreateConfig, CMSConfig, HashBiasTable};
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
        let table = HashBiasTable::load(path)?;

        if table.k() != kmer_size {
            return Err(anyhow::anyhow!(
                "Bias table k-mer size ({}) does not match sketch k-mer size ({})",
                table.k(),
                kmer_size
            ));
        }

        let sketch_fscale = fscale.unwrap_or(1000);
        if table.fscale() != sketch_fscale {
            return Err(anyhow::anyhow!(
                "Bias table fscale ({}) does not match sketch fscale ({})",
                table.fscale(),
                sketch_fscale
            ));
        }

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
    use std::time::Instant;

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

    let spinner = if !silent {
        let sp = ProgressBar::new_spinner();
        sp.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} [{elapsed_precise}] {msg}")
                .unwrap(),
        );
        sp.set_message("[1/4] Opening database...");
        sp.enable_steady_tick(std::time::Duration::from_millis(80));
        Some(sp)
    } else {
        None
    };

    let phase_start = Instant::now();
    let engine = QueryEngine::open(&database_path)?;
    let db_stats = engine.reader().stats();
    let db_names = engine.reader().sample_names();
    let db_sizes = engine.reader().sample_sizes();

    if let Some(ref sp) = spinner {
        sp.println(format!(
            "[1/4] Database opened in {:.2?}: {} samples, {} entries, {} threads",
            phase_start.elapsed(),
            db_stats.sample_count,
            db_stats.entry_count,
            rayon::current_num_threads()
        ));
        if engine.has_bias_table() {
            sp.println("      Using embedded bias table from database");
        }
        sp.set_message("[2/4] Loading query...");
    }

    let phase_start = Instant::now();
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

    if let Some(ref sp) = spinner {
        sp.println(format!(
            "[2/4] Query loaded in {:.2?}: {} samples, {} hashes",
            phase_start.elapsed(),
            sketch.sample_count(),
            sketch.total_entries()
        ));
        sp.set_message(format!(
            "[3/4] Searching {} query samples against {} db samples...",
            sketch.sample_count(),
            db_stats.sample_count
        ));
    }

    let phase_start = Instant::now();
    let results = engine.query_sketch(&sketch);

    if let Some(ref sp) = spinner {
        let total_matches: usize = results.iter().map(|r| r.matches.len()).sum();
        sp.println(format!(
            "[3/4] Search completed in {:.2?}: {} total matches",
            phase_start.elapsed(),
            total_matches
        ));
        sp.set_message("[4/4] Writing output...");
    }

    let phase_start = Instant::now();

    use rayon::prelude::*;

    let formatted_chunks: Vec<String> = results
        .par_iter()
        .enumerate()
        .map(|(query_idx, result)| {
            let query_name = &sketch.sample_names[query_idx];
            let query_hashes = result.query_size;

            let mut matches: Vec<_> = result
                .matches
                .iter()
                .filter(|m| m.containment >= cutoff)
                .collect();

            if matches.is_empty() {
                return String::new();
            }

            matches.sort_by(|a, b| b.containment.total_cmp(&a.containment));

            let mut chunk = String::with_capacity(matches.len() * 100);

            for m in &matches {
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

                use std::fmt::Write;
                let _ = writeln!(
                    chunk,
                    "{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}",
                    query_name,
                    db_name,
                    shared_hashes,
                    query_hashes,
                    db_hashes,
                    query_containment,
                    db_containment
                );
            }

            chunk
        })
        .collect();

    const WRITE_BUFFER_SIZE: usize = 64 * 1024 * 1024;
    let mut writer: Box<dyn Write> = if let Some(ref out) = output_path {
        Box::new(std::io::BufWriter::with_capacity(
            WRITE_BUFFER_SIZE,
            std::fs::File::create(out)?,
        ))
    } else {
        Box::new(std::io::BufWriter::with_capacity(
            WRITE_BUFFER_SIZE,
            std::io::stdout().lock(),
        ))
    };

    writeln!(
        writer,
        "query\tdb_sample\tshared_hashes\tquery_hashes\tdb_hashes\tquery_containment\tdb_containment"
    )?;

    for chunk in &formatted_chunks {
        if !chunk.is_empty() {
            writer.write_all(chunk.as_bytes())?;
        }
    }

    if let Some(ref sp) = spinner {
        sp.println(format!(
            "[4/4] Output written in {:.2?}",
            phase_start.elapsed()
        ));
        sp.finish_and_clear();
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn handle_bias_create_command(
    positive: Vec<PathBuf>,
    negative: Vec<PathBuf>,
    output: PathBuf,
    kmer_size: u8,
    fscale: u64,
    cms_width: usize,
    cms_depth: usize,
    alpha: f32,
    positive_fscale: Option<u64>,
    negative_fscale: Option<String>,
    unbiased_fscale: Option<u64>,
    min_positive_retention: f32,
    threads: Option<usize>,
    force: bool,
    silent: bool,
) -> Result<()> {
    use std::time::Instant;

    if positive.is_empty() {
        return Err(anyhow::anyhow!("No positive input files specified"));
    }
    if negative.is_empty() {
        return Err(anyhow::anyhow!("No negative input files specified"));
    }

    for path in positive.iter().chain(negative.iter()) {
        if !path.exists() {
            return Err(anyhow::anyhow!("Input file does not exist: {:?}", path));
        }
    }

    if output.exists() && !force {
        return Err(anyhow::anyhow!(
            "Output file {:?} already exists. Use --force to overwrite.",
            output
        ));
    }

    if !alpha.is_finite() || alpha <= 0.0 {
        return Err(anyhow::anyhow!("--alpha must be finite and > 0"));
    }

    if !min_positive_retention.is_finite() || !(0.0..=1.0).contains(&min_positive_retention) {
        return Err(anyhow::anyhow!(
            "--min-positive-retention must be in the range [0, 1]"
        ));
    }

    // Validate soft-filter fscale pair: both or neither
    match (positive_fscale.as_ref(), negative_fscale.as_ref()) {
        (Some(_), None) | (None, Some(_)) => {
            return Err(anyhow::anyhow!(
                "Both --positive-fscale and --negative-fscale must be set together"
            ));
        }
        _ => {}
    }

    let negative_fscale = match negative_fscale.as_deref() {
        Some(value) if value.eq_ignore_ascii_case("drop") => Some(u64::MAX),
        Some(value) => {
            let parsed = value
                .parse::<u64>()
                .map_err(|_| anyhow::anyhow!("--negative-fscale must be an integer or 'drop'"))?;
            Some(parsed)
        }
        None => None,
    };

    let spinner = if !silent {
        let sp = ProgressBar::new_spinner();
        sp.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} [{elapsed_precise}] {msg}")
                .unwrap(),
        );
        sp.set_message(format!(
            "Building bias table from {} positive + {} negative files...",
            positive.len(),
            negative.len()
        ));
        sp.enable_steady_tick(std::time::Duration::from_millis(80));
        Some(sp)
    } else {
        None
    };

    let start = Instant::now();

    let config = BiasCreateConfig {
        cms: CMSConfig {
            width: cms_width,
            depth: cms_depth,
            k: kmer_size,
            fscale,
        },
        alpha,
        positive_fscale,
        negative_fscale,
        unbiased_fscale,
        min_positive_retention,
    };

    let pos_paths: Vec<&std::path::Path> = positive.iter().map(|p| p.as_path()).collect();
    let neg_paths: Vec<&std::path::Path> = negative.iter().map(|p| p.as_path()).collect();

    if let Some(threads) = threads
        && threads == 0
    {
        return Err(anyhow::anyhow!("Thread count must be > 0"));
    }

    let table = if let Some(threads) = threads {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .stack_size(8 * 1024 * 1024)
            .build()?;
        pool.install(|| HashBiasTable::create(&pos_paths, &neg_paths, &config, spinner.clone()))?
    } else {
        HashBiasTable::create(&pos_paths, &neg_paths, &config, spinner.clone())?
    };

    if let Some(ref sp) = spinner {
        sp.set_message("Saving bias table...");
    }

    table.save(&output)?;

    if let Some(sp) = spinner {
        sp.finish_and_clear();
    }

    if !silent {
        eprintln!("Hash Bias Table");
        eprintln!("===============");
        eprintln!("Positive: {} files", positive.len());
        eprintln!("Negative: {} files", negative.len());
        eprintln!();
        eprintln!("Configuration:");
        eprintln!("  k-mer size:     {}", table.k());
        eprintln!("  fscale:         {}", table.fscale());
        eprintln!(
            "  CMS dimensions: {} x {}",
            table.config.width, table.config.depth
        );
        eprintln!("  Smoothing (alpha): {:.1}", alpha);
        if table.is_soft_filter() {
            eprintln!("  Filter mode: enrichment LUT");
        } else {
            eprintln!("  Filter mode: hard cutoff");
        }
        eprintln!();

        eprintln!("Calibration:");
        if table.is_soft_filter() {
            eprintln!("  positive_fscale:     {}", table.positive_fscale);
            eprintln!("  negative_fscale:     {}", table.negative_fscale_label());
            eprintln!("  min_pos_retention:   {:.2}", table.min_positive_retention);
        }
        eprintln!(
            "  positive retention:  {:.2}%",
            table.positive_retention * 100.0
        );
        eprintln!(
            "  negative retention:  {:.2}%",
            table.negative_retention * 100.0
        );
        eprintln!("  fold enrichment:     {:.2}x", table.fold_enrichment());
        if table.is_soft_filter() {
            eprintln!(
                "  eff. fscale (pos):   {:.0}",
                table.effective_fscale_on_population(table.positive_retention)
            );
            eprintln!(
                "  eff. fscale (neg):   {:.0}",
                table.effective_fscale_on_population(table.negative_retention)
            );
            eprintln!(
                "  eff. fscale (combined): {:.0}",
                table.effective_fscale_combined()
            );
            if table.unbiased_fscale > 0 {
                eprintln!("  unbiased_fscale:     {}", table.unbiased_fscale);
            }
            eprintln!("  reference points:");
            for p in table.soft_filter_reference_points() {
                eprintln!(
                    "    w={:>5.2}: eff={:>10.1}, ret={:>9.4}%, vs_base={:>7.2}x",
                    p.weight_f32, p.effective_fscale, p.retention_pct, p.vs_base
                );
            }
        }
        eprintln!(
            "  threshold:           {:.2} (quantized: {})",
            table.threshold_f32(),
            table.threshold
        );

        if table.fold_enrichment() < 1.5 {
            eprintln!();
            eprintln!(
                "Warning: Fold enrichment is very low ({:.2}x). The positive and negative \
                 sets may be too similar for effective filtering.",
                table.fold_enrichment()
            );
        }

        eprintln!();
        let (min, max, mean, std, positive_weights, above_threshold) = table.weight_stats();
        let total_cells = table.config.width * table.config.depth;
        eprintln!("Weight distribution (clamped to +/-12.70):");
        eprintln!("  min:    {:.2}", min);
        eprintln!("  max:    {:.2}", max);
        eprintln!("  mean:   {:.2}", mean);
        eprintln!("  std:    {:.2}", std);
        eprintln!(
            "  >0:             {} cells ({:.1}%)",
            positive_weights,
            positive_weights as f64 / total_cells as f64 * 100.0
        );
        eprintln!(
            "  >= threshold:   {} cells ({:.1}%)",
            above_threshold,
            above_threshold as f64 / total_cells as f64 * 100.0
        );

        if table.is_soft_filter() {
            table.print_lut_curve(min, max);
        }

        eprintln!();
        eprintln!("Saved to: {}", output.display());
        eprintln!("Built in {:.2?}", start.elapsed());
    }

    Ok(())
}

pub fn handle_bias_stats_command(
    input: PathBuf,
    output: Option<PathBuf>,
    silent: bool,
) -> Result<()> {
    if !input.exists() {
        return Err(anyhow::anyhow!("Input file does not exist: {:?}", input));
    }

    let table = HashBiasTable::load(&input)?;

    let (min, max, mean, std, positive_weights, above_threshold) = table.weight_stats();
    let total_cells = table.config.width * table.config.depth;

    let filter_mode = if table.is_soft_filter() {
        "enrichment LUT"
    } else {
        "hard cutoff"
    };

    if let Some(output_path) = output {
        let base = table.config.fscale as f64;

        let mut json = serde_json::json!({
            "file": input.display().to_string(),
            "type": "bias_v1",
            "k": table.config.k,
            "fscale": table.config.fscale,
            "cms_width": table.config.width,
            "cms_depth": table.config.depth,
            "alpha": table.alpha,
            "filter_mode": filter_mode,
            "calibration": {
                "threshold": table.threshold,
                "threshold_f32": table.threshold_f32(),
                "threshold_note": "informational only in LUT mode",
                "positive_retention": table.positive_retention,
                "negative_retention": table.negative_retention,
                "fold_enrichment": table.fold_enrichment(),
                "effective_fscale_positive": table.effective_fscale_on_population(table.positive_retention),
                "effective_fscale_negative": table.effective_fscale_on_population(table.negative_retention),
                "effective_fscale_combined": table.effective_fscale_combined(),
            },
            "soft_filter": {
                "positive_fscale": table.positive_fscale,
                "negative_fscale": table.negative_fscale,
                "negative_fscale_drop": table.negative_fscale == u64::MAX,
                "min_positive_retention": table.min_positive_retention,
                "unbiased_fscale": table.unbiased_fscale,
            },
            "weight_stats": {
                "min": min,
                "max": max,
                "mean": mean,
                "std": std,
                "positive_count": positive_weights,
                "positive_pct": positive_weights as f64 / total_cells as f64 * 100.0,
                "above_threshold_count": above_threshold,
                "above_threshold_pct": above_threshold as f64 / total_cells as f64 * 100.0,
            },
            "memory_bytes": table.memory_usage(),
        });

        if table.is_soft_filter() {
            let lut_curve: Vec<serde_json::Value> = (-127i8..=127i8)
                .map(|w| {
                    let eff = table.effective_fscale_at(w);
                    serde_json::json!({
                        "weight": w as f64 / 10.0,
                        "weight_q": w,
                        "effective_fscale": eff,
                        "retention_pct": 100.0 / eff,
                        "vs_base": base / eff,
                    })
                })
                .collect();
            json["lut_curve"] = serde_json::json!(lut_curve);

            let fscale_lut: Vec<u64> = table.fscale_lut.to_vec();
            json["soft_filter"]["fscale_lut"] = serde_json::json!(fscale_lut);
        }

        let file = std::fs::File::create(&output_path)?;
        serde_json::to_writer_pretty(file, &json)?;

        if !silent {
            eprintln!("Statistics written to: {}", output_path.display());
        }
    } else if !silent {
        table.print_stats();
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
        let base_fscale = u64::MAX / stats.hash_threshold.max(1);
        if let Some(bias) = reader.bias_table() {
            if bias.is_soft_filter() {
                let eff = bias.effective_fscale_combined().round() as u64;
                println!("Sample rate: 1/{} (effective: 1/{})", base_fscale, eff);
            } else {
                println!("Sample rate: 1/{}", base_fscale);
            }
        } else {
            println!("Sample rate: 1/{}", base_fscale);
        }
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
