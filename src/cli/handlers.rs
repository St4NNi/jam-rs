use anyhow::Result;
use indicatif::{ProgressBar, ProgressStyle};
use std::io::Write;
use std::{fs::remove_file, path::PathBuf};

use crate::bias::{
    BiasTable, CompareStats, EnrichedHexamer, RawBiasTable, idx_to_hexamer, reverse_complement_idx,
};
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

    // Immediate feedback before any loading
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
        sp.println(format!("[4/4] Output written in {:.2?}", phase_start.elapsed()));
        sp.finish_and_clear();
    }

    Ok(())
}

pub fn handle_bias_create_command(
    input: PathBuf,
    output: PathBuf,
    force: bool,
    silent: bool,
) -> Result<()> {
    use std::time::Instant;

    if !input.exists() {
        return Err(anyhow::anyhow!("Input file does not exist: {:?}", input));
    }
    if output.exists() && !force {
        return Err(anyhow::anyhow!(
            "Output file {:?} already exists. Use --force to overwrite.",
            output
        ));
    }

    let spinner = if !silent {
        let sp = ProgressBar::new_spinner();
        sp.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} [{elapsed_precise}] {msg}")
                .unwrap(),
        );
        sp.set_message("Processing input file...");
        sp.enable_steady_tick(std::time::Duration::from_millis(80));
        Some(sp)
    } else {
        None
    };

    let start = Instant::now();
    let raw = RawBiasTable::build(&input, spinner.clone())?;

    if let Some(ref sp) = spinner {
        sp.set_message("Saving raw bias table...");
    }

    raw.save(&output)?;

    if let Some(sp) = spinner {
        sp.finish_and_clear();
    }

    if !silent {
        eprintln!("Raw bias table created in {:.2?}", start.elapsed());
        eprintln!(
            "  Records: {}, Bases: {}, GC: {:.1}%",
            raw.stats.num_records,
            format_bp(raw.stats.total_bases),
            raw.stats.gc_content()
        );
        eprintln!("  Total hexamers: {}", raw.total_hexamers());
        eprintln!("Saved to: {}", output.display());
    }

    Ok(())
}

fn format_bp(bp: u64) -> String {
    if bp >= 1_000_000_000 {
        format!("{:.2} Gbp", bp as f64 / 1_000_000_000.0)
    } else if bp >= 1_000_000 {
        format!("{:.2} Mbp", bp as f64 / 1_000_000.0)
    } else if bp >= 1_000 {
        format!("{:.2} Kbp", bp as f64 / 1_000.0)
    } else {
        format!("{} bp", bp)
    }
}

pub fn handle_bias_combine_command(
    positive: PathBuf,
    negative: PathBuf,
    output: PathBuf,
    threshold: f32,
    force: bool,
    silent: bool,
) -> Result<()> {
    use std::time::Instant;

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

    let spinner = if !silent {
        let sp = ProgressBar::new_spinner();
        sp.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} [{elapsed_precise}] {msg}")
                .unwrap(),
        );
        sp.set_message("Loading raw bias tables...");
        sp.enable_steady_tick(std::time::Duration::from_millis(80));
        Some(sp)
    } else {
        None
    };

    let start = Instant::now();

    let pos_raw = RawBiasTable::load(&positive)?;
    let neg_raw = RawBiasTable::load(&negative)?;

    if let Some(ref sp) = spinner {
        sp.set_message("Computing comparison statistics...");
    }

    let compare_stats = pos_raw.compare(&neg_raw, threshold);

    if let Some(ref sp) = spinner {
        sp.set_message("Combining tables...");
    }

    let combined = pos_raw.combine(&neg_raw, threshold)?;

    if let Some(ref sp) = spinner {
        sp.set_message("Saving combined bias table...");
    }

    combined.save(&output)?;

    if let Some(sp) = spinner {
        sp.finish_and_clear();
    }

    if !silent {
        eprintln!("Bias Table Comparison");
        eprintln!("=====================");
        eprintln!(
            "Positive: {} ({} records, {}, {:.1}% GC)",
            positive.display(),
            pos_raw.stats.num_records,
            format_bp(pos_raw.stats.total_bases),
            pos_raw.stats.gc_content()
        );
        eprintln!(
            "Negative: {} ({} records, {}, {:.1}% GC)",
            negative.display(),
            neg_raw.stats.num_records,
            format_bp(neg_raw.stats.total_bases),
            neg_raw.stats.gc_content()
        );
        eprintln!();
        eprintln!("Distribution divergence:");
        eprintln!(
            "  KL divergence (pos||neg): {:.3} bits",
            compare_stats.kl_divergence
        );
        eprintln!("  GC content difference: {:+.1}%", compare_stats.gc_diff);
        eprintln!();
        print_threshold_sweep(&compare_stats);
        eprintln!();
        print_enriched_hexamers("POSITIVE", &compare_stats.top_positive_enriched);
        eprintln!();
        print_enriched_hexamers("NEGATIVE", &compare_stats.top_negative_enriched);
        eprintln!();
        eprintln!("Combined table created:");
        eprintln!("  Threshold: {:.2}", threshold);
        eprintln!(
            "  Hexamers above threshold: {} ({:.1}%)",
            compare_stats.hexamers_favoring_positive,
            compare_stats.hexamers_favoring_positive as f64 / 4096.0 * 100.0
        );
        eprintln!(
            "  Hexamers below threshold: {} ({:.1}%)",
            compare_stats.hexamers_favoring_negative,
            compare_stats.hexamers_favoring_negative as f64 / 4096.0 * 100.0
        );

        if let Some(point) = compare_stats
            .threshold_sweep
            .iter()
            .find(|p| (p.threshold - threshold).abs() < 0.01)
        {
            eprintln!();
            eprintln!("Filtering impact at this threshold:");
            eprintln!("  Expected positive k-mer retention: ~{:.1}%", point.pos_retain);
            eprintln!("  Expected negative k-mer retention: ~{:.1}%", point.neg_retain);
            eprintln!("  Selectivity ratio: {:.2}x", point.selectivity);
        }

        eprintln!();
        eprintln!("Saved to: {}", output.display());
        eprintln!("Built in {:.2?}", start.elapsed());
    }

    Ok(())
}

fn print_threshold_sweep(stats: &CompareStats) {
    eprintln!("Threshold sweep (expected retention rates):");
    eprintln!(
        "  {:>9}   {:>10}   {:>10}   {:>11}",
        "Threshold", "Pos Retain", "Neg Retain", "Selectivity"
    );
    for point in &stats.threshold_sweep {
        eprintln!(
            "  {:>9.2}   {:>9.1}%   {:>9.1}%   {:>10.2}x",
            point.threshold, point.pos_retain, point.neg_retain, point.selectivity
        );
    }
}

fn print_enriched_hexamers(label: &str, hexamers: &[EnrichedHexamer]) {
    eprintln!("Top 10 enriched in {}:", label);
    eprintln!(
        "  {:>8}   {:>8}   {:>10}   {:>8}   {:>8}   {:>6}",
        "Hexamer", "RevComp", "P(pos|hex)", "Pos Freq", "Neg Freq", "Fold"
    );
    for h in hexamers {
        let fold_str = if h.fold_change.is_infinite() {
            "inf".to_string()
        } else {
            format!("{:.1}x", h.fold_change)
        };
        eprintln!(
            "  {:>8}   {:>8}   {:>10.2}   {:>7.2}%   {:>7.2}%   {:>6}",
            h.hexamer, h.rev_comp, h.p_pos_given_hex, h.pos_freq, h.neg_freq, fold_str
        );
    }
}

pub fn handle_bias_stats_command(
    input: PathBuf,
    output: Option<PathBuf>,
    silent: bool,
) -> Result<()> {
    if !input.exists() {
        return Err(anyhow::anyhow!("Input file does not exist: {:?}", input));
    }

    let extension = input
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("");

    match extension {
        "braw" => handle_raw_stats(&input, output, silent),
        "bias" => handle_combined_stats(&input, output, silent),
        _ => Err(anyhow::anyhow!(
            "Unknown file extension: {}. Expected .braw or .bias",
            extension
        )),
    }
}

fn handle_raw_stats(input: &std::path::Path, output: Option<PathBuf>, silent: bool) -> Result<()> {
    let raw = RawBiasTable::load(input)?;

    let total_hexamers = raw.total_hexamers();
    let non_zero = raw.non_zero_count();
    let entropy = raw.shannon_entropy();
    let low_count = raw.low_count_hexamers(10);
    let top = raw.top_hexamers(10);

    if let Some(output_path) = output {
        let json = serde_json::json!({
            "file": input.display().to_string(),
            "type": "raw",
            "stats": {
                "num_records": raw.stats.num_records,
                "total_bases": raw.stats.total_bases,
                "gc_content": raw.stats.gc_content()
            },
            "hexamer_stats": {
                "total_counted": total_hexamers,
                "non_zero_count": non_zero,
                "shannon_entropy": entropy,
                "max_entropy": 12.0,
                "low_count_hexamers": low_count
            },
            "top_hexamers": top.iter().map(|(hex, count, freq)| {
                serde_json::json!({
                    "hexamer": hex,
                    "count": count,
                    "frequency": freq * 100.0
                })
            }).collect::<Vec<_>>()
        });

        let file = std::fs::File::create(&output_path)?;
        serde_json::to_writer_pretty(file, &json)?;

        if !silent {
            eprintln!("Statistics written to: {}", output_path.display());
        }
    } else if !silent {
        eprintln!("Raw Bias Table Statistics");
        eprintln!("=========================");
        eprintln!("File: {}", input.display());
        eprintln!("Records: {}", raw.stats.num_records);
        eprintln!("Total bases: {}", format_bp(raw.stats.total_bases));
        eprintln!("GC content: {:.1}%", raw.stats.gc_content());
        eprintln!();
        eprintln!("Hexamer distribution:");
        eprintln!("  Total counted: {}", total_hexamers);
        eprintln!("  Non-zero: {}/4096 ({:.1}%)", non_zero, non_zero as f64 / 40.96);
        eprintln!(
            "  Shannon entropy: {:.2} bits (max 12.0)",
            entropy
        );
        eprintln!("  Low-count (<10): {} hexamers", low_count);
        eprintln!();
        eprintln!("Top 10 most frequent:");
        eprintln!(
            "  {:>8}   {:>8}   {:>12}   {:>10}",
            "Hexamer", "RevComp", "Count", "Frequency"
        );
        for (hex, count, freq) in &top {
            let hex_idx = (0..6).fold(0usize, |acc, i| {
                let c = hex.chars().nth(i).unwrap();
                let base = match c {
                    'A' => 0,
                    'C' => 1,
                    'G' => 2,
                    'T' => 3,
                    _ => 0,
                };
                (acc << 2) | base
            });
            let rc = idx_to_hexamer(reverse_complement_idx(hex_idx));
            eprintln!(
                "  {:>8}   {:>8}   {:>12}   {:>9.2}%",
                hex, rc, count, freq * 100.0
            );
        }
    }

    Ok(())
}

fn handle_combined_stats(input: &std::path::Path, output: Option<PathBuf>, silent: bool) -> Result<()> {
    let table = BiasTable::load(input)?;

    let mut positive_count = 0usize;
    let mut negative_count = 0usize;

    for i in 0..4096 {
        let score = table.score(i);
        if score > 0.0 {
            positive_count += 1;
        } else if score < 0.0 {
            negative_count += 1;
        }
    }

    if let Some(output_path) = output {
        let json = serde_json::json!({
            "file": input.display().to_string(),
            "type": "combined",
            "threshold": table.threshold,
            "hexamers_above_threshold": positive_count,
            "hexamers_below_threshold": negative_count,
            "above_percentage": positive_count as f64 / 4096.0 * 100.0,
            "below_percentage": negative_count as f64 / 4096.0 * 100.0
        });

        let file = std::fs::File::create(&output_path)?;
        serde_json::to_writer_pretty(file, &json)?;

        if !silent {
            eprintln!("Statistics written to: {}", output_path.display());
        }
    } else if !silent {
        eprintln!("Combined Bias Table Statistics");
        eprintln!("==============================");
        eprintln!("File: {}", input.display());
        eprintln!("Threshold: {:.2}", table.threshold);
        eprintln!(
            "Hexamers above threshold: {}/4096 ({:.1}%)",
            positive_count,
            positive_count as f64 / 4096.0 * 100.0
        );
        eprintln!(
            "Hexamers below threshold: {}/4096 ({:.1}%)",
            negative_count,
            negative_count as f64 / 4096.0 * 100.0
        );
    }

    Ok(())
}

pub fn handle_bias_compare_command(
    positive: PathBuf,
    negative: PathBuf,
    output: Option<PathBuf>,
    threshold: f32,
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

    let pos_raw = RawBiasTable::load(&positive)?;
    let neg_raw = RawBiasTable::load(&negative)?;

    let stats = pos_raw.compare(&neg_raw, threshold);

    if let Some(output_path) = output {
        let json = serde_json::json!({
            "positive": {
                "file": positive.display().to_string(),
                "num_records": pos_raw.stats.num_records,
                "total_bases": pos_raw.stats.total_bases,
                "gc_content": pos_raw.stats.gc_content()
            },
            "negative": {
                "file": negative.display().to_string(),
                "num_records": neg_raw.stats.num_records,
                "total_bases": neg_raw.stats.total_bases,
                "gc_content": neg_raw.stats.gc_content()
            },
            "divergence": {
                "kl_divergence": stats.kl_divergence,
                "gc_difference": stats.gc_diff
            },
            "threshold_sweep": stats.threshold_sweep.iter().map(|p| {
                serde_json::json!({
                    "threshold": p.threshold,
                    "pos_retain": p.pos_retain,
                    "neg_retain": p.neg_retain,
                    "selectivity": p.selectivity
                })
            }).collect::<Vec<_>>(),
            "top_positive_enriched": stats.top_positive_enriched.iter().map(|h| {
                serde_json::json!({
                    "hexamer": h.hexamer,
                    "rev_comp": h.rev_comp,
                    "p_pos_given_hex": h.p_pos_given_hex,
                    "pos_freq": h.pos_freq,
                    "neg_freq": h.neg_freq,
                    "fold_change": h.fold_change
                })
            }).collect::<Vec<_>>(),
            "top_negative_enriched": stats.top_negative_enriched.iter().map(|h| {
                serde_json::json!({
                    "hexamer": h.hexamer,
                    "rev_comp": h.rev_comp,
                    "p_pos_given_hex": h.p_pos_given_hex,
                    "pos_freq": h.pos_freq,
                    "neg_freq": h.neg_freq,
                    "fold_change": h.fold_change
                })
            }).collect::<Vec<_>>(),
            "summary": {
                "threshold": threshold,
                "hexamers_favoring_positive": stats.hexamers_favoring_positive,
                "hexamers_favoring_negative": stats.hexamers_favoring_negative,
                "positive_percentage": stats.hexamers_favoring_positive as f64 / 4096.0 * 100.0,
                "negative_percentage": stats.hexamers_favoring_negative as f64 / 4096.0 * 100.0
            }
        });

        let file = std::fs::File::create(&output_path)?;
        serde_json::to_writer_pretty(file, &json)?;

        if !silent {
            eprintln!("Comparison written to: {}", output_path.display());
        }
    } else if !silent {
        eprintln!("Bias Table Comparison");
        eprintln!("=====================");
        eprintln!(
            "Positive: {} ({} records, {}, {:.1}% GC)",
            positive.display(),
            pos_raw.stats.num_records,
            format_bp(pos_raw.stats.total_bases),
            pos_raw.stats.gc_content()
        );
        eprintln!(
            "Negative: {} ({} records, {}, {:.1}% GC)",
            negative.display(),
            neg_raw.stats.num_records,
            format_bp(neg_raw.stats.total_bases),
            neg_raw.stats.gc_content()
        );
        eprintln!();
        eprintln!("Distribution divergence:");
        eprintln!(
            "  KL divergence (pos||neg): {:.3} bits",
            stats.kl_divergence
        );
        eprintln!("  GC content difference: {:+.1}%", stats.gc_diff);
        eprintln!();
        print_threshold_sweep(&stats);
        eprintln!();
        print_enriched_hexamers("POSITIVE", &stats.top_positive_enriched);
        eprintln!();
        print_enriched_hexamers("NEGATIVE", &stats.top_negative_enriched);
        eprintln!();
        eprintln!("Summary at threshold {:.2}:", threshold);
        eprintln!(
            "  Hexamers favoring positive: {} ({:.1}%)",
            stats.hexamers_favoring_positive,
            stats.hexamers_favoring_positive as f64 / 4096.0 * 100.0
        );
        eprintln!(
            "  Hexamers favoring negative: {} ({:.1}%)",
            stats.hexamers_favoring_negative,
            stats.hexamers_favoring_negative as f64 / 4096.0 * 100.0
        );
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
