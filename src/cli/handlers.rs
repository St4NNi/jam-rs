use anyhow::Result;
use indicatif::{ProgressBar, ProgressStyle};
use std::io::Write;
use std::{fs::remove_file, path::PathBuf};

use crate::bias::{BiasCreateConfig, CMSConfig, HashBiasTable};
use crate::format::VERSION;
use crate::provenance;
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
    let manifest_path = provenance::sidecar_path(&output_path);
    if manifest_path.exists() && !force {
        return Err(anyhow::anyhow!(
            "Database manifest {:?} already exists. Use --force to overwrite.",
            manifest_path
        ));
    }

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
    if !min_entropy.is_finite() || !(0.0..=2.0).contains(&min_entropy) {
        return Err(anyhow::anyhow!(
            "--complexity must be finite and between 0.0 and 2.0, got {}",
            min_entropy
        ));
    }
    if fscale == Some(0) {
        return Err(anyhow::anyhow!("--fscale must be > 0"));
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

        if fscale.is_some() {
            return Err(anyhow::anyhow!(
                "--fscale cannot be used with --bias-table. \
                 The bias table's stored fscale ({}) is used automatically.",
                table.fscale()
            ));
        }

        if !silent {
            table.print_stats();
        }
        Some(table)
    } else {
        None
    };

    let effective_fscale = match &bias_table {
        Some(table) => table.fscale(),
        None => fscale.unwrap_or(100),
    };

    let config = BuildConfig {
        kmer_size,
        fscale: effective_fscale,
        num_threads: threads,
        memory,
        singleton,
        min_entropy,
        temp_dir_base: temp_dir,
        bias_table: bias_table.map(Arc::new),
        show_progress: !silent,
    };

    let stats = build(&input_paths, &output_path, &config)?;

    let input_catalog_files = provenance::file_identities(&input_paths)?;
    let catalog_manifest_sha256 = provenance::identities_checksum(&input_catalog_files)?;
    let database_sha256 = provenance::sha256_file(&output_path)?;
    let bias = match (config.bias_table.as_deref(), bias_table_path.as_deref()) {
        (Some(table), Some(path)) => {
            let sha256 = provenance::sha256_file(path)?;
            Some(provenance::BiasDatabaseMetadata {
                table_id: format!("sha256:{sha256}"),
                source_path: path.display().to_string(),
                sha256,
                kmer_size: table.config.k,
                base_fscale: table.config.fscale,
                cms_width: table.config.width,
                cms_depth: table.config.depth,
                alpha: table.alpha,
                filter_mode: if table.is_soft_filter() {
                    "enrichment_lut".to_string()
                } else {
                    "hard_cutoff".to_string()
                },
                target_fscale: table.target_fscale,
                negative_fscale: table.negative_fscale_label(),
                unseen_fscale: table.unseen_fscale,
                positive_retention: table.positive_retention,
                negative_retention: table.negative_retention,
            })
        }
        _ => None,
    };
    let manifest = provenance::DatabaseManifest {
        schema_version: provenance::MANIFEST_SCHEMA_VERSION.to_string(),
        output_schema_version: provenance::OUTPUT_SCHEMA_VERSION.to_string(),
        database_format_version: VERSION,
        jam_rs_version: env!("CARGO_PKG_VERSION").to_string(),
        source_commit: provenance::source_commit(),
        source_dirty: provenance::source_dirty(),
        hash_id: provenance::HASH_ID.to_string(),
        hash_zero_policy: provenance::HASH_ZERO_POLICY.to_string(),
        kmer_size,
        fscale: effective_fscale,
        hash_threshold: stats.frac_max,
        entropy_threshold: min_entropy as f32 as f64,
        bias,
        input_catalog_files,
        catalog_manifest_sha256,
        sample_count: stats.sample_count,
        entry_count: stats.total_entries,
        unique_hash_count: stats.unique_hashes,
        database_file: output_path.display().to_string(),
        database_size_bytes: stats.file_size,
        database_sha256,
        creation_command: provenance::command_line(),
        creation_time_unix_seconds: provenance::unix_time_seconds(),
    };
    provenance::write_json(&manifest_path, &manifest)?;

    if !silent {
        eprintln!(
            "Completed: {} ({} entries, {} unique hashes, {} samples)",
            output_path.display(),
            stats.total_entries,
            stats.unique_hashes,
            stats.sample_count
        );
        eprintln!("Manifest: {}", manifest_path.display());
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn handle_distance_command(
    input_path: PathBuf,
    database_path: PathBuf,
    output_path: Option<PathBuf>,
    cutoff: f64,
    singleton: bool,
    force: bool,
    silent: bool,
    memory_gb: usize,
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
    if !cutoff.is_finite() {
        return Err(anyhow::anyhow!("--cutoff must be finite, got {}", cutoff));
    }
    if !(0.0..=1.0).contains(&cutoff) {
        return Err(anyhow::anyhow!(
            "--cutoff must be between 0.0 and 1.0, got {}",
            cutoff
        ));
    }
    if let Some(ref out) = output_path
        && out.exists()
    {
        if !out.is_file() {
            return Err(anyhow::anyhow!(
                "Output path must be a file, not a directory: {:?}",
                out
            ));
        }
        if !force {
            return Err(anyhow::anyhow!(
                "Output file {:?} already exists. Use --force to overwrite.",
                out
            ));
        }
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
            sp.println(
                "      Bias mode reports retained/weighted sketch evidence; uniform-hash E-values are NA",
            );
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

    let total_samples = sketch.sample_count();
    let cutoff = normalize_distance_cutoff(cutoff);
    let budget_bytes = memory_gb * 1024 * 1024 * 1024 * 7 / 10;
    let per_sample_bytes = (db_stats.sample_count as usize).max(1) * 40;
    let chunk_size = compute_distance_chunk_size(total_samples, budget_bytes, per_sample_bytes);

    if let Some(ref sp) = spinner {
        sp.set_message(format!(
            "[3/4] Searching {} query samples against {} db samples (chunk size {})...",
            total_samples, db_stats.sample_count, chunk_size,
        ));
    }

    use rayon::prelude::*;

    const WRITE_BUFFER_SIZE: usize = 1024 * 1024;
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

    let has_bias = engine.has_bias_table();
    let bias_table_id = engine
        .bias_table()
        .map(|table| format!("sha256:{}", provenance::sha256_bytes(&table.to_bytes())))
        .unwrap_or_else(|| "NA".to_string());
    if has_bias {
        writeln!(
            writer,
            "query\tdb_sample\tshared_hashes\tquery_hashes\tdb_hashes\tretained_query_containment\tretained_reference_containment\tuniform_hash_e_value\tbias_weighted_query_containment\tscore_mode\tbias_table_id"
        )?;
    } else {
        writeln!(
            writer,
            "query\tdb_sample\tshared_hashes\tquery_hashes\tdb_hashes\tquery_containment\tdb_containment\tuniform_hash_e_value"
        )?;
    }

    let flush_threshold = (memory_gb * 1024 * 1024 * 1024 / 16).clamp(1024 * 1024, 8 * 1024 * 1024);
    let mut temp_files: Vec<tempfile::NamedTempFile> = Vec::new();
    let mut pending: Vec<u8> = Vec::new();

    let phase_start = Instant::now();
    let mut total_matches: usize = 0;

    for chunk_start in (0..total_samples).step_by(chunk_size) {
        let chunk_end = (chunk_start + chunk_size).min(total_samples);

        let results =
            engine.query_sketch_chunked(&sketch, chunk_start..chunk_end, cutoff.unwrap_or(0.0));

        let formatted: Vec<String> = results
            .par_iter()
            .enumerate()
            .map(|(local_idx, result)| {
                let global_idx = chunk_start + local_idx;
                let query_name = &sketch.sample_names[global_idx];
                let query_hashes = result.query_size;

                if result.matches.is_empty() {
                    return String::new();
                }

                let mut matches: Vec<_> = result.matches.iter().collect();
                matches.sort_by(|a, b| {
                    b.containment
                        .total_cmp(&a.containment)
                        .then_with(|| b.hit_count.cmp(&a.hit_count))
                        .then_with(|| {
                            let left = db_names
                                .get(a.sample_id as usize)
                                .map(String::as_str)
                                .unwrap_or("unknown");
                            let right = db_names
                                .get(b.sample_id as usize)
                                .map(String::as_str)
                                .unwrap_or("unknown");
                            left.cmp(right)
                        })
                        .then_with(|| a.sample_id.cmp(&b.sample_id))
                });

                let mut out = String::with_capacity(matches.len() * 120);
                for m in &matches {
                    let db_name = db_names
                        .get(m.sample_id as usize)
                        .map(|s| s.as_str())
                        .unwrap_or("unknown");
                    let db_hashes = db_sizes.get(m.sample_id as usize).copied().unwrap_or(0);
                    let db_containment = if db_hashes > 0 {
                        m.hit_count as f64 / db_hashes as f64
                    } else {
                        0.0
                    };
                    use std::fmt::Write;
                    if has_bias {
                        let retained_query_containment = if query_hashes > 0 {
                            m.hit_count as f64 / query_hashes as f64
                        } else {
                            0.0
                        };
                        let _ = writeln!(
                            out,
                            "{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\tNA\t{:.6}\tbias\t{}",
                            query_name,
                            db_name,
                            m.hit_count,
                            query_hashes,
                            db_hashes,
                            retained_query_containment,
                            db_containment,
                            m.containment,
                            bias_table_id,
                        );
                    } else {
                        let _ = writeln!(
                            out,
                            "{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6e}",
                            query_name,
                            db_name,
                            m.hit_count,
                            query_hashes,
                            db_hashes,
                            m.containment,
                            db_containment,
                            m.e_value,
                        );
                    }
                }
                out
            })
            .collect();

        for s in &formatted {
            if !s.is_empty() {
                total_matches += s.lines().count();
                pending.extend_from_slice(s.as_bytes());
            }
        }

        if pending.len() > flush_threshold {
            let mut tmp = tempfile::NamedTempFile::new()?;
            tmp.write_all(&pending)?;
            temp_files.push(tmp);
            pending.clear();
        }
    }

    if let Some(ref sp) = spinner {
        sp.println(format!(
            "[3/4] Search completed in {:.2?}: {} total matches",
            phase_start.elapsed(),
            total_matches,
        ));
        sp.set_message("[4/4] Writing output...");
    }

    let phase_start = Instant::now();

    for tmp in &temp_files {
        std::io::copy(&mut std::fs::File::open(tmp.path())?, &mut writer)?;
    }
    writer.write_all(&pending)?;

    if let Some(ref sp) = spinner {
        sp.println(format!(
            "[4/4] Output written in {:.2?}",
            phase_start.elapsed()
        ));
        sp.finish_and_clear();
    }

    Ok(())
}

pub fn handle_screen_command(config: crate::screen::ScreenConfig, silent: bool) -> Result<()> {
    let stats = crate::screen::run(&config)?;
    if !silent {
        eprintln!(
            "Screened {} contigs: {} contig rows, {} assembly rows",
            stats.contig_count, stats.contig_rows, stats.summary_rows
        );
    }
    Ok(())
}

fn normalize_distance_cutoff(cutoff: f64) -> Option<f64> {
    if cutoff.is_finite() && cutoff > 0.0 {
        Some(cutoff)
    } else {
        None
    }
}

fn compute_distance_chunk_size(
    total_samples: usize,
    budget_bytes: usize,
    per_sample_bytes: usize,
) -> usize {
    debug_assert!(total_samples > 0);

    let min_chunk_size = total_samples.clamp(1, 100);
    let raw_chunk_size = budget_bytes / per_sample_bytes.max(1);

    raw_chunk_size.clamp(min_chunk_size, total_samples)
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
    target_fscale: Option<u64>,
    max_fscale: Option<String>,
    unseen_fscale: Option<u64>,
    threads: Option<usize>,
    min_positive_retention: f32,
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
    let manifest_path = provenance::sidecar_path(&output);
    if manifest_path.exists() && !force {
        return Err(anyhow::anyhow!(
            "Bias manifest {:?} already exists. Use --force to overwrite.",
            manifest_path
        ));
    }

    if !alpha.is_finite() || alpha <= 0.0 {
        return Err(anyhow::anyhow!("--alpha must be finite and > 0"));
    }
    if !(1..=31).contains(&kmer_size) {
        return Err(anyhow::anyhow!(
            "K-mer size must be between 1 and 31, got {}",
            kmer_size
        ));
    }
    if fscale == 0 {
        return Err(anyhow::anyhow!("--fscale must be > 0"));
    }
    if !min_positive_retention.is_finite() || !(0.0..=1.0).contains(&min_positive_retention) {
        return Err(anyhow::anyhow!(
            "--min-positive-retention must be finite and between 0.0 and 1.0"
        ));
    }

    match (target_fscale.as_ref(), max_fscale.as_ref()) {
        (Some(_), None) | (None, Some(_)) => {
            return Err(anyhow::anyhow!(
                "Both --target-fscale and --max-fscale must be set together"
            ));
        }
        _ => {}
    }

    let negative_fscale = match max_fscale.as_deref() {
        Some(value) if value.eq_ignore_ascii_case("drop") => Some(u64::MAX),
        Some(value) => {
            let parsed = value
                .parse::<u64>()
                .map_err(|_| anyhow::anyhow!("--max-fscale must be an integer or 'drop'"))?;
            Some(parsed)
        }
        None => None,
    };

    let unseen_fscale = unseen_fscale.or(target_fscale);

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
        target_fscale,
        negative_fscale,
        unseen_fscale,
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

    if table.positive_retention < min_positive_retention {
        return Err(anyhow::anyhow!(
            "Bias calibration retained {:.2}% of the positive set, below --min-positive-retention {:.2}%",
            table.positive_retention * 100.0,
            min_positive_retention * 100.0
        ));
    }

    if let Some(ref sp) = spinner {
        sp.set_message("Saving bias table...");
    }

    table.save(&output)?;

    let table_sha256 = provenance::sha256_file(&output)?;
    let manifest = provenance::BiasTableManifest {
        schema_version: provenance::MANIFEST_SCHEMA_VERSION.to_string(),
        jam_rs_version: env!("CARGO_PKG_VERSION").to_string(),
        source_commit: provenance::source_commit(),
        source_dirty: provenance::source_dirty(),
        hash_id: provenance::HASH_ID.to_string(),
        hash_zero_policy: provenance::HASH_ZERO_POLICY.to_string(),
        table_id: format!("sha256:{table_sha256}"),
        table_file: output.display().to_string(),
        table_size_bytes: output.metadata()?.len(),
        table_sha256,
        kmer_size: table.config.k,
        base_fscale: table.config.fscale,
        cms_width: table.config.width,
        cms_depth: table.config.depth,
        alpha: table.alpha,
        filter_mode: if table.is_soft_filter() {
            "enrichment_lut".to_string()
        } else {
            "hard_cutoff".to_string()
        },
        target_fscale: table.target_fscale,
        negative_fscale: table.negative_fscale_label(),
        unseen_fscale: table.unseen_fscale,
        positive_retention: table.positive_retention,
        negative_retention: table.negative_retention,
        minimum_positive_retention: min_positive_retention,
        positive_files: provenance::file_identities(&positive)?,
        chromosome_background_files: provenance::file_identities(&negative)?,
        creation_command: provenance::command_line(),
        creation_time_unix_seconds: provenance::unix_time_seconds(),
    };
    provenance::write_json(&manifest_path, &manifest)?;

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
            eprintln!("  requested target:    {}", table.target_fscale);
            eprintln!("  negative_fscale:     {}", table.negative_fscale_label());
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
            let target_fs = table.target_fscale as f64;
            let eff_target = table.effective_fscale_target_prior();
            let eff_pos = table.effective_fscale_on_population(table.positive_retention);
            let eff_neg = table.effective_fscale_on_population(table.negative_retention);
            let delta_pct = |achieved: f64| {
                if target_fs > 0.0 {
                    (achieved / target_fs - 1.0) * 100.0
                } else {
                    0.0
                }
            };
            eprintln!(
                "  eff. fscale (target prior): {:.0} ({:+.1}%)",
                eff_target,
                delta_pct(eff_target)
            );
            eprintln!(
                "  eff. fscale (pos):      {:.0} ({:+.1}%)",
                eff_pos,
                delta_pct(eff_pos)
            );
            eprintln!(
                "  eff. fscale (neg):      {:.0} ({:+.1}%)",
                eff_neg,
                delta_pct(eff_neg)
            );
            eprintln!(
                "  eff. fscale (combined): {:.0}",
                table.effective_fscale_combined()
            );
            if table.unseen_fscale > 0 {
                eprintln!("  unseen fscale:       {}", table.unseen_fscale);
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
                "effective_fscale_target_prior": table.effective_fscale_target_prior(),
                "effective_fscale_target_prior_delta_pct": if table.target_fscale > 0 {
                    (table.effective_fscale_target_prior() / table.target_fscale as f64 - 1.0) * 100.0
                } else {
                    0.0
                },
                "effective_fscale_combined": table.effective_fscale_combined(),
            },
            "soft_filter": {
                "target_fscale": table.target_fscale,
                "negative_fscale": table.negative_fscale,
                "negative_fscale_drop": table.negative_fscale == u64::MAX,
                "unseen_fscale": table.unseen_fscale,
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
                        "retention_pct": 100.0 * base / eff,
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
    json: bool,
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

    if json {
        let manifest = provenance::load_database_manifest(&input_path)?;
        let database_sha256 = match manifest.as_ref() {
            Some(manifest) => manifest.database_sha256.clone(),
            None => provenance::sha256_file(&input_path)?,
        };
        let bias_table_id = manifest
            .as_ref()
            .and_then(|manifest| manifest.bias.as_ref())
            .map(|bias| bias.table_id.clone())
            .or_else(|| {
                reader
                    .bias_table()
                    .map(|table| format!("sha256:{}", provenance::sha256_bytes(&table.to_bytes())))
            });
        let value = serde_json::json!({
            "schema_version": provenance::OUTPUT_SCHEMA_VERSION,
            "database_format_version": VERSION,
            "jam_rs_version": env!("CARGO_PKG_VERSION"),
            "source_commit": provenance::source_commit(),
            "hash_id": manifest.as_ref().map(|manifest| manifest.hash_id.as_str()).unwrap_or(provenance::HASH_ID),
            "hash_zero_policy": manifest.as_ref().map(|manifest| manifest.hash_zero_policy.as_str()).unwrap_or("legacy database; zero excluded at query time"),
            "database_file": input_path.display().to_string(),
            "database_sha256": database_sha256,
            "file_size_bytes": stats.file_size,
            "kmer_size": stats.kmer_size,
            "fscale": manifest.as_ref().map(|manifest| manifest.fscale).unwrap_or(u64::MAX / stats.hash_threshold.max(1)),
            "hash_threshold": stats.hash_threshold,
            "entropy_threshold": reader.min_entropy(),
            "bias_mode": stats.has_bias_table,
            "bias_table_id": bias_table_id,
            "entry_count": stats.entry_count,
            "unique_hash_count": stats.unique_hash_count,
            "sample_count": stats.sample_count,
            "bucket_entry_counts": stats.bucket_entry_counts.to_vec(),
            "manifest": manifest,
        });
        println!("{}", serde_json::to_string_pretty(&value)?);
        return Ok(());
    }

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

#[allow(clippy::too_many_arguments)]
pub fn handle_archive_command(
    input: PathBuf,
    output: PathBuf,
    block_bases: usize,
    primary_scale: u64,
    rescue_scale: Option<u64>,
    complexity: Option<f64>,
    force: bool,
    silent: bool,
) -> Result<()> {
    if !input.is_file() {
        return Err(anyhow::anyhow!(
            "Assembly input does not exist or is not a file: {}",
            input.display()
        ));
    }
    if output.exists() {
        if !output.is_file() {
            return Err(anyhow::anyhow!(
                "Archive output is not a file: {}",
                output.display()
            ));
        }
        if !force {
            return Err(anyhow::anyhow!(
                "Archive output already exists: {}. Use --force to overwrite.",
                output.display()
            ));
        }
    }
    let stats = crate::jma::builder::write_archive_from_fasta(
        &input,
        &output,
        crate::jma::builder::ArchiveBuildConfig {
            block_bases,
            k31_scale: primary_scale,
            k21_scale: rescue_scale,
            min_entropy: complexity,
            flags: crate::jma::builder::DEFAULT_FLAGS,
        },
    )?;
    if !silent {
        eprintln!(
            "Created JMA v{} archive {}: {} contigs, {} bases, {} k=31 seeds, {} k=21 seeds",
            crate::jma::JMA_FORMAT_VERSION,
            output.display(),
            stats.contig_count,
            stats.total_bases,
            stats.k31_seed_count,
            stats.k21_seed_count
        );
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn handle_trace_command(
    plasmid: PathBuf,
    database: String,
    catalog_path: PathBuf,
    output: PathBuf,
    plasmid_id: Option<String>,
    profile: crate::trace::config::SensitivityProfile,
    min_shared_hashes: u32,
    min_plasmid_containment: f64,
    min_metagenome_containment: f64,
    top_candidates: Option<usize>,
    max_alignments: usize,
    threads: usize,
    memory_gb: usize,
    cache_dir: Option<PathBuf>,
    cache_block_bytes: u64,
    request_timeout_seconds: u64,
    max_retries: u32,
    allow_full_download_fallback: bool,
    force: bool,
    silent: bool,
) -> Result<()> {
    if output.exists() {
        if !output.is_file() {
            return Err(anyhow::anyhow!(
                "Trace output is not a file: {}",
                output.display()
            ));
        }
        if !force {
            return Err(anyhow::anyhow!(
                "Trace output already exists: {}. Use --force to overwrite.",
                output.display()
            ));
        }
    }

    let resources = crate::resource::ResourceOpenOptions {
        cache_dir,
        cache_block_bytes,
        max_cache_bytes: (memory_gb as u64)
            .saturating_mul(1024 * 1024 * 1024)
            .saturating_mul(7)
            / 10,
        request_timeout_seconds,
        max_retries,
        allow_full_download_fallback,
    };
    let catalog = crate::trace::catalog::TraceCatalog::from_path(&catalog_path)?;
    let parsed_query =
        crate::trace::raw::RawAssembly::open(plasmid.to_string_lossy(), resources.clone())?;
    if parsed_query.contigs.len() != 1 {
        return Err(anyhow::anyhow!(
            "Plasmid input must contain exactly one record, got {}",
            parsed_query.contigs.len()
        ));
    }
    let record = &parsed_query.contigs[0];
    let plasmid_id = plasmid_id.unwrap_or_else(|| record.id.clone());
    let sensitivity = crate::trace::config::SensitivityConfig::for_profile(profile);
    let candidate_limit = top_candidates
        .unwrap_or_else(|| usize::try_from(sensitivity.max_candidates).unwrap_or(usize::MAX));
    let runner = crate::trace::runner::TraceRunner::new(crate::trace::runner::TraceRunnerConfig {
        sensitivity: sensitivity.clone(),
        candidates: crate::trace::screen::CandidateSearchConfig {
            min_shared_hashes,
            min_plasmid_containment,
            min_metagenome_containment,
            top_candidates: candidate_limit,
        },
        resources,
        threads,
        max_alignments_per_candidate: max_alignments,
    })?;
    let query = crate::trace::runner::TraceQuery {
        plasmid_id: plasmid_id.clone(),
        plasmid_sequence: record.sequence.clone(),
        database: database.clone(),
        catalog,
    };
    let started = provenance::unix_time_seconds();
    let run_id = format!("trace-{started}-{}", std::process::id());
    let mut result = runner.run(&query)?;
    result.set_run_id(&run_id);

    let source_commit = provenance::source_commit();
    let algorithm =
        crate::trace::config::TraceAlgorithmMetadata::for_sensitivity(sensitivity.clone());
    let header = crate::trace::model::TraceRunHeader {
        schema_version: crate::trace::TRACE_JSON_SCHEMA_VERSION.to_string(),
        run_id: run_id.clone(),
        jam_rs_version: env!("CARGO_PKG_VERSION").to_string(),
        source_commit: (source_commit != "unknown").then_some(source_commit),
        started_at_utc: format!("unix:{started}"),
        command: provenance::command_line(),
        plasmid_id,
        plasmid_length: record.sequence.len() as u64,
        sensitivity,
        algorithm,
        inputs: vec![
            trace_input("plasmid", &plasmid)?,
            result.database_input.clone(),
            trace_input("catalog", &catalog_path)?,
        ],
    };
    let resource_metrics = result
        .metagenomes
        .iter()
        .fold(result.database_metrics, |total, item| {
            add_resource_metrics(total, item.resource_metrics)
        });
    let footer = crate::trace::model::TraceRunFooter {
        schema_version: crate::trace::TRACE_JSON_SCHEMA_VERSION.to_string(),
        run_id,
        completed_at_utc: format!("unix:{}", provenance::unix_time_seconds()),
        metagenomes_total: result.metagenomes.len() as u64,
        metagenomes_with_candidates: result
            .metagenomes
            .iter()
            .filter(|item| item.candidate.is_some())
            .count() as u64,
        metagenomes_aligned: result
            .metagenomes
            .iter()
            .filter(|item| !item.alignments.is_empty())
            .count() as u64,
        metagenomes_failed: result
            .metagenomes
            .iter()
            .filter(|item| item.status == crate::trace::model::TraceStatus::Failed)
            .count() as u64,
        alignments_total: result
            .metagenomes
            .iter()
            .map(|item| item.alignments.len() as u64)
            .sum(),
        resource_metrics,
    };
    let mut writer = crate::trace::output::create(&output)?;
    writer.write_header(&header)?;
    for item in &result.metagenomes {
        writer.write_metagenome_result(item)?;
    }
    writer.write_footer(&footer)?;
    writer.finish()?;
    if !silent {
        eprintln!(
            "Trace output {}: {} candidates, {} aligned metagenomes, {} alignments",
            output.display(),
            result.search.candidates.len(),
            footer.metagenomes_aligned,
            footer.alignments_total
        );
    }
    Ok(())
}

fn trace_input(role: &str, path: &std::path::Path) -> Result<crate::trace::model::InputResource> {
    let locator = crate::resource::ResourceLocator::parse(path.to_string_lossy().as_ref())?;
    Ok(crate::trace::model::InputResource {
        role: role.to_string(),
        redacted_locator: locator.redacted(),
        sha256: Some(provenance::sha256_file(path)?),
    })
}

fn add_resource_metrics(
    left: crate::resource::ResourceMetrics,
    right: crate::resource::ResourceMetrics,
) -> crate::resource::ResourceMetrics {
    left.saturating_add(right)
}

#[cfg(test)]
mod tests {
    use super::{compute_distance_chunk_size, normalize_distance_cutoff};

    #[test]
    fn distance_chunk_size_handles_small_query_counts() {
        let chunk_size = compute_distance_chunk_size(3, 0, 1024);
        assert_eq!(chunk_size, 3);
    }

    #[test]
    fn distance_chunk_size_uses_default_minimum_for_large_queries() {
        let chunk_size = compute_distance_chunk_size(1000, 0, 1024);
        assert_eq!(chunk_size, 100);
    }

    #[test]
    fn distance_cutoff_non_positive_is_disabled() {
        assert_eq!(normalize_distance_cutoff(0.0), None);
        assert_eq!(normalize_distance_cutoff(-0.1), None);
    }

    #[test]
    fn distance_cutoff_positive_is_kept() {
        assert_eq!(normalize_distance_cutoff(0.25), Some(0.25));
    }
}
