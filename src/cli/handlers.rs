use anyhow::Result;
use indicatif::{ProgressBar, ProgressStyle};
use needletail::Sequence;
use std::io::Write;
use std::{
    fs::remove_file,
    path::{Path, PathBuf},
};

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
    sequence_policy: crate::sequence::SequenceBlockPolicy,
    sequence_codec: crate::sequence::BlockCodec,
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
            sequence_policy,
            sequence_codec,
            k31_scale: primary_scale,
            k21_scale: rescue_scale,
            min_entropy: complexity,
            flags: crate::jma::builder::DEFAULT_FLAGS,
        },
    )?;
    if !silent {
        eprintln!(
            "Created JMA format {} archive {}: {} contigs, {} bases, {} k=31 seeds, {} k=21 seeds",
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
pub fn handle_router_build_command(
    metagenomes_path: PathBuf,
    witness_metagenomes_path: Option<PathBuf>,
    output: PathBuf,
    base_scale: u32,
    tiers: &str,
    key_block_bytes: u32,
    positions_per_metagenome: usize,
    position_max_document_frequency: u32,
    force: bool,
    silent: bool,
) -> Result<()> {
    if output.exists() && !force {
        return Err(anyhow::anyhow!(
            "Router output already exists: {}. Use --force to overwrite.",
            output.display()
        ));
    }
    let catalog = crate::trace::catalog::TraceCatalog::from_path(&metagenomes_path)?;
    let metagenomes = catalog
        .entries()
        .iter()
        .map(|entry| {
            Ok(crate::router::format::MetagenomeEntry {
                id: entry.metagenome_id.clone(),
                object_uri: entry.resource_uri.clone(),
                checksum: parse_sha256(&entry.sha256)?,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let witness_metagenomes = if let Some(path) = witness_metagenomes_path {
        crate::trace::catalog::TraceCatalog::from_path(path)?
            .entries()
            .iter()
            .map(|entry| {
                Ok(crate::router::format::MetagenomeEntry {
                    id: entry.metagenome_id.clone(),
                    object_uri: entry.resource_uri.clone(),
                    checksum: parse_sha256(&entry.sha256)?,
                })
            })
            .collect::<Result<Vec<_>>>()?
    } else {
        metagenomes.clone()
    };
    let scales = tiers
        .split(',')
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(|value| {
            value
                .parse::<u32>()
                .map_err(|_| anyhow::anyhow!("invalid witness tier {value:?}"))
        })
        .collect::<Result<Vec<_>>>()?;
    let stats = crate::router::build::build_router_from_local_jma_sources(
        &witness_metagenomes,
        &metagenomes,
        &output,
        &crate::router::build::RouterCollectionBuildConfig {
            base_scale,
            available_scales: scales,
            rare_document_frequency: position_max_document_frequency,
            positions_per_metagenome,
            writer: crate::router::writer::RouterWriterOptions {
                key_blocks: crate::router::format::KeyBlockConfig {
                    target_block_bytes: key_block_bytes,
                    packed_width: crate::router::format::PackedKeyWidth::Six,
                    store_jamhash: false,
                },
                ..crate::router::writer::RouterWriterOptions::default()
            },
            ..crate::router::build::RouterCollectionBuildConfig::default()
        },
    )?;
    if !silent {
        eprintln!(
            "Created JAM Witness Router {}: {} metagenomes, {} exact keys, {} posting IDs, {} positional occurrences, {} bytes",
            output.display(),
            stats.metagenomes,
            stats.keys,
            stats.posting_ids,
            stats.positional_occurrences,
            stats.bytes
        );
    }
    Ok(())
}

pub struct IndexBuildArgs {
    pub metagenomes: PathBuf,
    pub output: PathBuf,
    pub policy: crate::jam_index::ScreenSelectionPolicy,
    pub max_part_bases: u64,
    pub max_part_signatures: u64,
    pub target_parts: usize,
    pub parallel_parts: usize,
    pub force: bool,
    pub silent: bool,
}

#[allow(clippy::too_many_arguments)]
pub fn handle_index_plan(
    metagenomes: PathBuf,
    output: PathBuf,
    parts: usize,
    fragments_per_part: usize,
    estimated_expansion: u64,
    selection_policy: crate::jam_index::ScreenSelectionPolicy,
    force: bool,
    silent: bool,
) -> Result<()> {
    if force {
        return Err(anyhow::anyhow!(
            "--force is not supported for deterministic index plans"
        ));
    }
    if output.exists() {
        return Err(anyhow::anyhow!(
            "Index plan output already exists: {}",
            output.display()
        ));
    }
    let sources = index_plan_sources(&metagenomes)?;
    let plan = crate::jam_index::plan_index(
        sources,
        provenance::sha256_file(&metagenomes)?,
        selection_policy,
        parts,
        fragments_per_part,
        estimated_expansion,
    )?;
    crate::jam_index::write_plan_atomic(&output, &plan)?;
    if !silent {
        eprintln!(
            "Jam Index plan {}: {} parts, {} fragments, {} metagenomes",
            output.display(),
            plan.parts.len(),
            plan.parts
                .iter()
                .map(|part| part.fragments.len())
                .sum::<usize>(),
            plan.parts
                .iter()
                .flat_map(|part| &part.fragments)
                .map(|fragment| fragment.sources.len())
                .sum::<usize>()
        );
    }
    Ok(())
}

pub fn handle_index_build_fragment(
    plan_path: PathBuf,
    fragment_id: u32,
    staged_metagenomes: PathBuf,
    output: PathBuf,
    force: bool,
    silent: bool,
) -> Result<()> {
    if force {
        return Err(anyhow::anyhow!(
            "--force is not supported for restartable index fragments"
        ));
    }
    let plan = crate::jam_index::load_plan(plan_path)?;
    let staged = index_staged_sources(&staged_metagenomes, &plan, fragment_id)?;
    let manifest = crate::jam_index::build_fragment(&plan, fragment_id, &staged, &output)?;
    if !silent {
        eprintln!(
            "Jam Index fragment {}: part {}, {} metagenomes, {} contigs, {} bases",
            fragment_id,
            manifest.part_id,
            manifest.metagenome_count,
            manifest.contig_count,
            manifest.total_bases
        );
    }
    Ok(())
}

pub fn handle_index_merge_part(
    plan_path: PathBuf,
    part_id: u32,
    fragments_root: PathBuf,
    output: PathBuf,
    force: bool,
    silent: bool,
) -> Result<()> {
    if force {
        return Err(anyhow::anyhow!(
            "--force is not supported for immutable merged parts"
        ));
    }
    let plan = crate::jam_index::load_plan(plan_path)?;
    let manifest = crate::jam_index::merge_part(&plan, part_id, fragments_root, &output)?;
    if !silent {
        eprintln!(
            "Jam Index merged part {}: {} fragments, {} metagenomes, {} contigs, {} bases",
            part_id,
            manifest.fragment_ids.len(),
            manifest.part.metagenome_count,
            manifest.part.contig_count,
            manifest.part.total_bases
        );
    }
    Ok(())
}

pub fn handle_index_finalize(
    plan_path: PathBuf,
    output: PathBuf,
    force: bool,
    silent: bool,
) -> Result<()> {
    if force {
        return Err(anyhow::anyhow!(
            "--force is not supported for final index publication"
        ));
    }
    let plan = crate::jam_index::load_plan(plan_path)?;
    let stats = crate::jam_index::finalize_index(&plan, &output)?;
    print_index(&output, stats, silent);
    Ok(())
}

pub fn handle_index_build(args: IndexBuildArgs) -> Result<()> {
    if args.force {
        return Err(anyhow::anyhow!(
            "--force is not supported for append-only Jam Index construction"
        ));
    }
    let sources = index_sources(&args.metagenomes)?;
    let stats = crate::jam_index::build_jam_index(
        &args.output,
        &sources,
        &crate::jam_index::JamIndexBuildConfig {
            selection_policy: args.policy,
            max_part_bases: args.max_part_bases,
            max_part_signatures: args.max_part_signatures,
            target_parts: args.target_parts,
            parallel_parts: args.parallel_parts,
            source_manifest_sha256: provenance::sha256_file(&args.metagenomes)?,
        },
    )?;
    print_index(&args.output, stats, args.silent);
    Ok(())
}

pub fn handle_index_append(args: IndexBuildArgs) -> Result<()> {
    if args.force {
        return Err(anyhow::anyhow!(
            "--force is not supported for append-only Jam Index construction"
        ));
    }
    let sources = index_sources(&args.metagenomes)?;
    let stats = crate::jam_index::append_jam_index(
        &args.output,
        &sources,
        &crate::jam_index::JamIndexBuildConfig {
            selection_policy: args.policy,
            max_part_bases: args.max_part_bases,
            max_part_signatures: args.max_part_signatures,
            target_parts: args.target_parts,
            parallel_parts: args.parallel_parts,
            source_manifest_sha256: provenance::sha256_file(&args.metagenomes)?,
        },
    )?;
    print_index(&args.output, stats, args.silent);
    Ok(())
}

fn index_sources(path: &Path) -> Result<Vec<crate::jam_index::MetagenomeSource>> {
    let catalog = crate::trace::catalog::TraceCatalog::from_path(path)?;
    catalog
        .entries()
        .iter()
        .map(|entry| {
            let locator = crate::resource::ResourceLocator::parse(entry.resource())?;
            if locator.scheme() != crate::resource::ResourceScheme::Local {
                return Err(anyhow::anyhow!(
                    "Jam Index sources must be local files: {}",
                    locator.redacted()
                ));
            }
            let sequence_path = PathBuf::from(entry.resource());
            let observed = provenance::sha256_file(&sequence_path)?;
            if observed != entry.sha256 {
                return Err(anyhow::anyhow!(
                    "source checksum mismatch for {}",
                    entry.metagenome_id
                ));
            }
            Ok(crate::jam_index::MetagenomeSource {
                metagenome_id: entry.metagenome_id.clone(),
                sequence_path,
            })
        })
        .collect()
}

fn index_source_overrides(
    path: &Path,
) -> Result<std::collections::BTreeMap<String, (PathBuf, [u8; 32])>> {
    let catalog = crate::trace::catalog::TraceCatalog::from_path(path)?;
    catalog
        .entries()
        .iter()
        .map(|entry| {
            let locator = crate::resource::ResourceLocator::parse(entry.resource())?;
            if locator.scheme() != crate::resource::ResourceScheme::Local {
                return Err(anyhow::anyhow!(
                    "Jam Index source overrides must be local files: {}",
                    locator.redacted()
                ));
            }
            Ok((
                entry.metagenome_id.clone(),
                (
                    PathBuf::from(entry.resource()),
                    crate::jam_index::distributed::parse_checksum(&entry.sha256)?,
                ),
            ))
        })
        .collect()
}

fn index_plan_sources(path: &Path) -> Result<Vec<crate::jam_index::IndexPlanSource>> {
    let catalog = crate::trace::catalog::TraceCatalog::from_path(path)?;
    catalog
        .entries()
        .iter()
        .map(|entry| {
            let locator = crate::resource::ResourceLocator::parse(entry.resource())?;
            if locator.scheme() != crate::resource::ResourceScheme::Local {
                return Err(anyhow::anyhow!(
                    "Jam Index sources must be local files: {}",
                    locator.redacted()
                ));
            }
            let source_path = PathBuf::from(entry.resource());
            let source_size = std::fs::metadata(&source_path)?.len();
            if source_size == 0
                || entry.sha256.len() != 64
                || !entry.sha256.bytes().all(|byte| byte.is_ascii_hexdigit())
            {
                return Err(anyhow::anyhow!(
                    "invalid source identity for {}",
                    entry.metagenome_id
                ));
            }
            Ok(crate::jam_index::IndexPlanSource {
                metagenome_id: entry.metagenome_id.clone(),
                source_path,
                source_size,
                source_sha256: entry.sha256.clone(),
                estimated_bases: 0,
                estimated_signatures: 0,
            })
        })
        .collect()
}

fn index_staged_sources(
    path: &Path,
    plan: &crate::jam_index::IndexBuildPlan,
    fragment_id: u32,
) -> Result<std::collections::BTreeMap<String, crate::jam_index::MetagenomeSource>> {
    let fragment = plan
        .fragment(fragment_id)
        .ok_or_else(|| anyhow::anyhow!("unknown build fragment {fragment_id}"))?;
    let planned = fragment
        .sources
        .iter()
        .map(|source| (source.metagenome_id.as_str(), source))
        .collect::<std::collections::BTreeMap<_, _>>();
    let catalog = crate::trace::catalog::TraceCatalog::from_path(path)?;
    let mut staged = std::collections::BTreeMap::new();
    for entry in catalog.entries() {
        let Some(expected) = planned.get(entry.metagenome_id.as_str()) else {
            continue;
        };
        let locator = crate::resource::ResourceLocator::parse(entry.resource())?;
        if locator.scheme() != crate::resource::ResourceScheme::Local
            || entry.sha256 != expected.source_sha256
        {
            return Err(anyhow::anyhow!(
                "staged source identity mismatch for {}",
                entry.metagenome_id
            ));
        }
        let sequence_path = PathBuf::from(entry.resource());
        if std::fs::metadata(&sequence_path)?.len() != expected.source_size {
            return Err(anyhow::anyhow!(
                "staged source size mismatch for {}",
                entry.metagenome_id
            ));
        }
        staged.insert(
            entry.metagenome_id.clone(),
            crate::jam_index::MetagenomeSource {
                metagenome_id: entry.metagenome_id.clone(),
                sequence_path,
            },
        );
    }
    if staged.len() != planned.len() {
        return Err(anyhow::anyhow!(
            "staged catalog does not contain every source for fragment {fragment_id}"
        ));
    }
    Ok(staged)
}

fn print_index(output: &Path, stats: crate::jam_index::JamIndexBuildStats, silent: bool) {
    if !silent {
        eprintln!(
            "Jam Index {}: {} new part(s), {} total metagenomes, {} total contigs, {} total bases",
            output.display(),
            stats.new_parts,
            stats.total_metagenomes,
            stats.total_contigs,
            stats.total_bases
        );
    }
}

pub struct IndexTraceArgs {
    pub query: PathBuf,
    pub query_id: Option<String>,
    pub query_kind: crate::trace::model::QueryKind,
    pub topology: crate::trace::model::TopologyRequested,
    pub index: PathBuf,
    pub output: PathBuf,
    pub profile: crate::trace::config::SensitivityProfile,
    pub min_shared: u32,
    pub min_query_windows: u32,
    pub rare_rescue_df: Option<u32>,
    pub whole_sample_min_shared: u32,
    pub top_candidates: Option<usize>,
    pub initial_contigs: usize,
    pub max_contigs: usize,
    pub max_contig_bases: u64,
    pub expansion_batch: usize,
    pub max_alignments: usize,
    pub threads: usize,
    pub topology_margin: Option<u64>,
    pub memory_gb: usize,
    pub force: bool,
    pub silent: bool,
}

pub struct IndexBatchTraceArgs {
    pub queries: Vec<PathBuf>,
    pub query_id: Option<String>,
    pub index: PathBuf,
    pub source_catalog: Option<PathBuf>,
    pub output: PathBuf,
    pub candidates: Option<PathBuf>,
    pub work: Option<PathBuf>,
    pub status: Option<PathBuf>,
    pub query_kind: crate::trace::model::QueryKind,
    pub topology: crate::trace::model::TopologyRequested,
    pub profile: crate::trace::config::SensitivityProfile,
    pub min_shared: u32,
    pub min_query_windows: u32,
    pub rare_rescue_df: Option<u32>,
    pub whole_sample_min_shared: u32,
    pub screen_only: bool,
    pub top_candidates: Option<usize>,
    pub initial_contigs: usize,
    pub max_contigs: usize,
    pub max_contig_bases: u64,
    pub expansion_batch: usize,
    pub max_alignments: usize,
    pub max_group_contig_bases: u64,
    pub fallback_contigs_per_chunk: usize,
    pub threads: usize,
    pub topology_margin: Option<u64>,
    pub memory_gb: usize,
    pub force: bool,
    pub silent: bool,
}

pub fn handle_index_trace(args: IndexTraceArgs) -> Result<()> {
    if args.output.exists() && !args.force {
        return Err(anyhow::anyhow!(
            "Trace output already exists: {}. Use --force to overwrite.",
            args.output.display()
        ));
    }
    if !args.index.is_dir() {
        return Err(anyhow::anyhow!(
            "Jam Index is not a local directory: {}",
            args.index.display()
        ));
    }
    let (record_id, sequence) = index_query(&args.query)?;
    let query_id = args.query_id.unwrap_or(record_id);
    let (sensitivity, trace_config) = index_trace_configuration(
        args.profile,
        args.min_shared,
        args.min_query_windows,
        args.rare_rescue_df,
        args.whole_sample_min_shared,
        args.top_candidates,
        args.initial_contigs,
        args.max_contigs,
        args.max_contig_bases,
        args.expansion_batch,
        args.max_alignments,
        args.threads,
        args.topology_margin,
        args.memory_gb,
        args.query_kind,
        args.topology,
        sequence.len(),
    );
    let mut result =
        crate::jam_index::trace_index(&args.index, &query_id, &sequence, &trace_config)?;
    let started = provenance::unix_time_seconds();
    let run_id = format!("trace-{started}-{}", std::process::id());
    for metagenome in &mut result.metagenomes {
        metagenome.run_id = run_id.clone();
    }
    let header = crate::trace::model::TraceRunHeader {
        schema_version: crate::trace::TRACE_JSON_SCHEMA_VERSION.to_string(),
        run_id: run_id.clone(),
        jam_rs_version: env!("CARGO_PKG_VERSION").to_string(),
        source_commit: Some(provenance::source_commit()),
        started_at_utc: format!("unix:{started}"),
        command: provenance::redacted_command_line(),
        plasmid_id: query_id,
        plasmid_length: sequence.len() as u64,
        query_kind: args.query_kind,
        topology_requested: args.topology,
        threads: args.threads,
        io_concurrency: args.threads,
        sensitivity: sensitivity.clone(),
        algorithms: crate::trace::config::algorithm_identifiers(),
        algorithm: crate::trace::config::TraceAlgorithmMetadata::for_sensitivity(sensitivity),
        inputs: vec![
            trace_input("query", &args.query)?,
            trace_input(
                "jam_index",
                &args
                    .index
                    .join(crate::jam_index::builder::JAM_INDEX_MANIFEST_FILE),
            )?,
        ],
    };
    let resources = result.metagenomes.iter().fold(
        crate::resource::ResourceMetrics::default(),
        |total, item| add_resource_metrics(total, item.resource_metrics),
    );
    let footer = crate::trace::model::TraceRunFooter {
        schema_version: crate::trace::TRACE_JSON_SCHEMA_VERSION.to_string(),
        run_id,
        completed_at_utc: format!("unix:{}", provenance::unix_time_seconds()),
        metagenomes_total: result.metagenomes.len() as u64,
        metagenomes_with_candidates: result.metagenomes.len() as u64,
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
        resource_metrics: resources,
    };
    let mut writer = crate::trace::output::create(&args.output)?;
    writer.write_header(&header)?;
    for item in &result.metagenomes {
        writer.write_metagenome_result(item)?;
    }
    writer.write_footer(&footer)?;
    writer.finish()?;
    if !args.silent {
        eprintln!(
            "Jam Index trace {}: {} candidates, {} contigs, {} selected bases",
            args.output.display(),
            result.trace_metrics.selected_candidates,
            result.trace_metrics.selected_contigs,
            result.trace_metrics.selected_contig_bases
        );
    }
    Ok(())
}

pub fn handle_index_batch_trace(args: IndexBatchTraceArgs) -> Result<()> {
    let profile_session = crate::profiling::ProfileSession::start();
    let queries = {
        let _profile = crate::profiling::scope("query_parsing");
        index_batch_queries(&args.queries, args.query_kind, args.topology)?
    };
    if queries.len() == 1
        && args.queries.len() == 1
        && args.candidates.is_none()
        && args.work.is_none()
        && args.status.is_none()
        && args.source_catalog.is_none()
        && !args.screen_only
    {
        return handle_index_trace(IndexTraceArgs {
            query: args.queries[0].clone(),
            query_id: args.query_id,
            query_kind: args.query_kind,
            topology: args.topology,
            index: args.index,
            output: args.output,
            profile: args.profile,
            min_shared: args.min_shared,
            min_query_windows: args.min_query_windows,
            rare_rescue_df: args.rare_rescue_df,
            whole_sample_min_shared: args.whole_sample_min_shared,
            top_candidates: args.top_candidates,
            initial_contigs: args.initial_contigs,
            max_contigs: args.max_contigs,
            max_contig_bases: args.max_contig_bases,
            expansion_batch: args.expansion_batch,
            max_alignments: args.max_alignments,
            threads: args.threads,
            topology_margin: args.topology_margin,
            memory_gb: args.memory_gb,
            force: args.force,
            silent: args.silent,
        });
    }
    if args.query_id.is_some() {
        return Err(anyhow::anyhow!(
            "--query-id is valid only when the combined input contains one record"
        ));
    }
    if args.force {
        return Err(anyhow::anyhow!(
            "--force is not supported for restartable batch trace outputs"
        ));
    }
    let candidates = args
        .candidates
        .clone()
        .unwrap_or_else(|| batch_sidecar(&args.output, "candidates"));
    let work = args
        .work
        .clone()
        .unwrap_or_else(|| batch_sidecar(&args.output, "work"));
    let status = args
        .status
        .clone()
        .unwrap_or_else(|| batch_sidecar(&args.output, "status"));
    let metrics_output = batch_json_sidecar(&args.output, "metrics");
    for output in [&args.output, &candidates, &work, &status, &metrics_output] {
        if output.exists() {
            return Err(anyhow::anyhow!(
                "Batch trace output already exists: {}",
                output.display()
            ));
        }
    }
    if !args.index.is_dir() {
        return Err(anyhow::anyhow!(
            "Jam Index is not a local directory: {}",
            args.index.display()
        ));
    }
    let maximum_query = queries
        .iter()
        .map(|query| query.sequence.len())
        .max()
        .unwrap_or(0);
    let (sensitivity, trace_config) = index_trace_configuration(
        args.profile,
        args.min_shared,
        args.min_query_windows,
        args.rare_rescue_df,
        args.whole_sample_min_shared,
        args.top_candidates,
        args.initial_contigs,
        args.max_contigs,
        args.max_contig_bases,
        args.expansion_batch,
        args.max_alignments,
        args.threads,
        args.topology_margin,
        args.memory_gb,
        args.query_kind,
        args.topology,
        maximum_query,
    );
    let config = crate::jam_index::JamIndexBatchConfig {
        trace: trace_config,
        max_group_contig_bases: args.max_group_contig_bases,
        fallback_contigs_per_chunk: args.fallback_contigs_per_chunk,
        group_memory_budget_bytes: (args.memory_gb as u64).saturating_mul(1024 * 1024 * 1024),
        ..crate::jam_index::JamIndexBatchConfig::default()
    };
    let mut plan = crate::jam_index::plan_batch(&args.index, &queries, &config)?;
    if let Some(source_catalog) = &args.source_catalog {
        plan.remap_sources(&index_source_overrides(source_catalog)?)?;
    }
    let mut candidate_writer = BatchTextWriter::create(&candidates)?;
    candidate_writer.write_line(&[
        "query_id",
        "metagenome_id",
        "part_id",
        "screen_policy",
        "shared_hash_count",
        "shared_spatial_signatures",
        "rare_shared_signatures",
        "shared_whole_sample_signatures",
        "supported_query_windows",
        "longest_supported_window_run",
        "weighted_hash_sum",
        "candidate_entry_reason",
        "candidate_rank",
        "candidate_truncated",
        "shared_hashes_json",
    ])?;
    for query in &plan.queries {
        for candidate in &query.candidates {
            let hashes = candidate
                .shared_hashes
                .iter()
                .map(|shared| format!("{:016x}", shared.hash))
                .collect::<Vec<_>>();
            let admission_source = serde_json::to_string(&candidate.admission_source)?;
            candidate_writer.write_line(&[
                &query.input.query_id,
                &candidate.metagenome_id,
                &candidate.part_id.to_string(),
                &candidate.screen_policy,
                &candidate.shared_hash_count.to_string(),
                &candidate.shared_spatial_signatures.to_string(),
                &candidate.rare_shared_signatures.to_string(),
                &candidate.shared_whole_sample_signatures.to_string(),
                &candidate.supported_query_windows.to_string(),
                &candidate.longest_supported_window_run.to_string(),
                &candidate.weighted_hash_sum.to_string(),
                admission_source.trim_matches('"'),
                &candidate.rank.to_string(),
                if query.candidate_limit_reached {
                    "true"
                } else {
                    "false"
                },
                &serde_json::to_string(&hashes)?,
            ])?;
        }
    }
    candidate_writer.finish()?;
    let mut work_writer = BatchTextWriter::create(&work)?;
    work_writer.write_line(&[
        "metagenome_id",
        "part_id",
        "metagenome_local_id",
        "source_assembly_path",
        "query_work_json",
        "estimated_alignment_work",
    ])?;
    for group in &plan.groups {
        let query_work = group
            .work
            .iter()
            .map(|work| {
                let query = &plan.queries[work.query_index];
                serde_json::json!({
                    "query_id": query.input.query_id,
                    "selected_contig_ids": work.contig_plan.ranked_contigs
                        .iter()
                        .map(|contig| contig.contig_id)
                        .collect::<Vec<_>>(),
                    "shared_hashes": work.candidate.shared_hashes
                        .iter()
                        .map(|shared| format!("{:016x}", shared.hash))
                        .collect::<Vec<_>>(),
                    "sequential_fallback": work.contig_plan.sequential_fallback_range.is_some(),
                    "contig_truncated": work.contig_plan.contig_truncated,
                })
            })
            .collect::<Vec<_>>();
        let estimated_work = group.work.iter().fold(0u64, |total, work| {
            let query_bases = u64::try_from(plan.queries[work.query_index].input.sequence.len())
                .unwrap_or(u64::MAX);
            let contig_bases = work
                .contig_plan
                .ranked_contigs
                .iter()
                .map(|contig| contig.base_count)
                .fold(0u64, u64::saturating_add);
            total.saturating_add(query_bases.saturating_mul(contig_bases))
        });
        work_writer.write_line(&[
            &group.metagenome_id,
            &group.part_id.to_string(),
            &group.metagenome_local_id.to_string(),
            &group.source_path.to_string_lossy(),
            &serde_json::to_string(&query_work)?,
            &estimated_work.to_string(),
        ])?;
    }
    work_writer.finish()?;
    let started = provenance::unix_time_seconds();
    let run_id = format!("trace-batch-{started}-{}", std::process::id());
    let total_query_bases = queries
        .iter()
        .map(|query| u64::try_from(query.sequence.len()).unwrap_or(u64::MAX))
        .fold(0u64, u64::saturating_add);
    let mut inputs = args
        .queries
        .iter()
        .map(|path| trace_input("queries", path))
        .chain(std::iter::once(trace_input(
            "jam_index",
            &args
                .index
                .join(crate::jam_index::builder::JAM_INDEX_MANIFEST_FILE),
        )))
        .collect::<Result<Vec<_>>>()?;
    if let Some(source_catalog) = &args.source_catalog {
        inputs.push(trace_input("staged_metagenomes", source_catalog)?);
    }
    let header = crate::trace::model::TraceRunHeader {
        schema_version: crate::trace::TRACE_JSON_SCHEMA_VERSION.to_string(),
        run_id: run_id.clone(),
        jam_rs_version: env!("CARGO_PKG_VERSION").to_string(),
        source_commit: Some(provenance::source_commit()),
        started_at_utc: format!("unix:{started}"),
        command: provenance::redacted_command_line(),
        plasmid_id: format!("batch:{}", queries.len()),
        plasmid_length: total_query_bases,
        query_kind: args.query_kind,
        topology_requested: args.topology,
        threads: args.threads,
        io_concurrency: args.threads,
        sensitivity: sensitivity.clone(),
        algorithms: crate::trace::config::algorithm_identifiers(),
        algorithm: crate::trace::config::TraceAlgorithmMetadata::for_sensitivity(sensitivity),
        inputs,
    };
    let mut trace_writer = crate::trace::output::create(&args.output)?;
    trace_writer.write_header(&header)?;
    let execution = if args.screen_only {
        index_screen_only_execution(&plan)
    } else {
        crate::jam_index::execute_batch(&plan, &config, &run_id, |result| {
            trace_writer
                .write_metagenome_result(&result)
                .map_err(|error| crate::jam_index::JamIndexBatchError::Output(error.to_string()))
        })?
    };
    let footer = crate::trace::model::TraceRunFooter {
        schema_version: crate::trace::TRACE_JSON_SCHEMA_VERSION.to_string(),
        run_id,
        completed_at_utc: format!("unix:{}", provenance::unix_time_seconds()),
        metagenomes_total: execution.metrics.results_emitted,
        metagenomes_with_candidates: execution.metrics.results_emitted,
        metagenomes_aligned: execution
            .statuses
            .iter()
            .map(|status| status.aligned_metagenomes)
            .sum(),
        metagenomes_failed: execution
            .statuses
            .iter()
            .filter(|status| status.failed)
            .count() as u64,
        alignments_total: execution.metrics.alignments_emitted,
        resource_metrics: crate::resource::ResourceMetrics::default(),
    };
    trace_writer.write_footer(&footer)?;
    trace_writer.finish()?;
    let mut status_writer = BatchTextWriter::create(&status)?;
    status_writer.write_line(&[
        "query_id",
        "original_header",
        "query_length",
        "query_sha256",
        "query_kind",
        "topology_requested",
        "status",
        "candidate_metagenomes",
        "completed_metagenomes",
        "aligned_metagenomes",
        "candidate_truncated",
        "failed",
        "error",
    ])?;
    for status in &execution.statuses {
        status_writer.write_line(&[
            &status.query_id,
            &status.original_header,
            &status.query_length.to_string(),
            &status.query_sha256,
            &format!("{:?}", status.query_kind).to_ascii_lowercase(),
            &format!("{:?}", status.topology_requested).to_ascii_lowercase(),
            &format!("{:?}", status.status).to_ascii_lowercase(),
            &status.candidate_metagenomes.to_string(),
            &status.completed_metagenomes.to_string(),
            &status.aligned_metagenomes.to_string(),
            if status.candidate_truncated {
                "true"
            } else {
                "false"
            },
            if status.failed { "true" } else { "false" },
            status.error.as_deref().unwrap_or(""),
        ])?;
    }
    status_writer.finish()?;
    let profiling = profile_session.report();
    let metrics_document = serde_json::json!({
        "schema_version": 1,
        "run_id": header.run_id,
        "binary_version": env!("CARGO_PKG_VERSION"),
        "source_identity": provenance::source_commit(),
        "query_files": args.queries,
        "index": args.index,
        "source_catalog": args.source_catalog,
        "screen_only": args.screen_only,
        "profiling": profiling,
        "execution": execution,
        "queries": plan.queries.iter().map(|query| serde_json::json!({
            "query_id": query.input.query_id,
            "query_length": query.input.sequence.len(),
            "query_sha256": query.input.sequence_sha256,
            "candidate_limit_reached": query.candidate_limit_reached,
            "screen_metrics": query.screen_metrics,
            "contig_metrics": query.contig_metrics,
            "candidate_count": query.candidates.len(),
            "contig_plan_count": query.contig_plans.len(),
        })).collect::<Vec<_>>(),
        "work_groups": plan.groups.iter().map(|group| serde_json::json!({
            "part_id": group.part_id,
            "metagenome_id": group.metagenome_id,
            "source_path": group.source_path,
            "queries": group.work.len(),
        })).collect::<Vec<_>>(),
    });
    write_json_atomic_file(&metrics_output, &metrics_document)?;
    if !args.silent {
        eprintln!(
            "Jam Index batch trace {}: {} queries, {} candidates, {} source opens, {} results, {} alignments",
            args.output.display(),
            execution.metrics.queries,
            execution.metrics.candidate_pairs,
            execution.metrics.source_open_count,
            execution.metrics.results_emitted,
            execution.metrics.alignments_emitted,
        );
    }
    Ok(())
}

fn index_screen_only_execution(
    plan: &crate::jam_index::JamIndexBatchPlan,
) -> crate::jam_index::JamIndexBatchExecution {
    let statuses = plan
        .queries
        .iter()
        .map(|query| crate::jam_index::JamIndexBatchQueryStatus {
            query_id: query.input.query_id.clone(),
            original_header: query.input.original_header.clone(),
            query_length: u64::try_from(query.input.sequence.len()).unwrap_or(u64::MAX),
            query_sha256: query.input.sequence_sha256.clone(),
            query_kind: query.input.query_kind,
            topology_requested: query.input.topology_requested,
            status: if query.candidates.is_empty() {
                crate::jam_index::JamIndexBatchStatusKind::NoCandidate
            } else {
                crate::jam_index::JamIndexBatchStatusKind::CandidateOnly
            },
            candidate_metagenomes: u64::try_from(query.candidates.len()).unwrap_or(u64::MAX),
            completed_metagenomes: 0,
            aligned_metagenomes: 0,
            candidate_truncated: query.candidate_limit_reached,
            failed: false,
            error: None,
        })
        .collect();
    crate::jam_index::JamIndexBatchExecution {
        statuses,
        metrics: crate::jam_index::JamIndexBatchMetrics {
            queries: u64::try_from(plan.queries.len()).unwrap_or(u64::MAX),
            query_bases: plan
                .queries
                .iter()
                .map(|query| u64::try_from(query.input.sequence.len()).unwrap_or(u64::MAX))
                .fold(0u64, u64::saturating_add),
            screen_parts_opened: plan.screen_parts_opened,
            contig_parts_opened: plan.contig_parts_opened,
            candidate_pairs: plan
                .queries
                .iter()
                .map(|query| u64::try_from(query.candidates.len()).unwrap_or(u64::MAX))
                .fold(0u64, u64::saturating_add),
            work_groups: u64::try_from(plan.groups.len()).unwrap_or(u64::MAX),
            ..crate::jam_index::JamIndexBatchMetrics::default()
        },
    }
}

fn index_batch_queries(
    paths: &[PathBuf],
    query_kind: crate::trace::model::QueryKind,
    topology_requested: crate::trace::model::TopologyRequested,
) -> Result<Vec<crate::jam_index::JamIndexBatchQuery>> {
    use sha2::{Digest, Sha256};

    let mut queries = Vec::new();
    for path in paths {
        let mut reader = needletail::parse_fastx_file(path)?;
        while let Some(record) = reader.next() {
            let record = record?;
            let original_header = std::str::from_utf8(record.id())?.to_string();
            let query_id = original_header
                .split_ascii_whitespace()
                .next()
                .ok_or_else(|| anyhow::anyhow!("query FASTA header is empty"))?
                .to_string();
            let sequence = record.normalize(true).into_owned();
            queries.push(crate::jam_index::JamIndexBatchQuery {
                query_id,
                original_header,
                sequence_sha256: format!("{:x}", Sha256::digest(&sequence)),
                sequence,
                query_kind,
                topology_requested,
            });
        }
    }
    if queries.is_empty() {
        return Err(anyhow::anyhow!("query batch contains no records"));
    }
    Ok(queries)
}

fn batch_sidecar(output: &Path, kind: &str) -> PathBuf {
    let mut base = output.to_string_lossy().into_owned();
    for suffix in [".zstd", ".zst", ".jsonl"] {
        if base.to_ascii_lowercase().ends_with(suffix) {
            base.truncate(base.len() - suffix.len());
        }
    }
    PathBuf::from(format!("{base}.{kind}.tsv.zst"))
}

fn batch_json_sidecar(output: &Path, kind: &str) -> PathBuf {
    let mut base = output.to_string_lossy().into_owned();
    for suffix in [".zstd", ".zst", ".jsonl"] {
        if base.to_ascii_lowercase().ends_with(suffix) {
            base.truncate(base.len() - suffix.len());
        }
    }
    PathBuf::from(format!("{base}.{kind}.json"))
}

fn write_json_atomic_file<T: serde::Serialize>(path: &Path, value: &T) -> Result<()> {
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let mut temporary = tempfile::NamedTempFile::new_in(parent)?;
    serde_json::to_writer_pretty(&mut temporary, value)?;
    temporary.write_all(b"\n")?;
    temporary.as_file_mut().sync_all()?;
    temporary.persist(path).map_err(|error| error.error)?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn index_trace_configuration(
    profile: crate::trace::config::SensitivityProfile,
    min_shared: u32,
    min_query_windows: u32,
    rare_rescue_df: Option<u32>,
    whole_sample_min_shared: u32,
    top_candidates: Option<usize>,
    initial_contigs: usize,
    max_contigs: usize,
    max_contig_bases: u64,
    expansion_batch: usize,
    max_alignments: usize,
    threads: usize,
    topology_margin: Option<u64>,
    memory_gb: usize,
    query_kind: crate::trace::model::QueryKind,
    topology: crate::trace::model::TopologyRequested,
    maximum_query_bases: usize,
) -> (
    crate::trace::config::SensitivityConfig,
    crate::jam_index::JamIndexTraceConfig,
) {
    let sensitivity = index_sensitivity(profile);
    let candidate_limit = top_candidates
        .unwrap_or_else(|| usize::try_from(sensitivity.max_candidates).unwrap_or(usize::MAX));
    let total_memory = (memory_gb as u64).saturating_mul(1024 * 1024 * 1024);
    let alignment_workers = index_workers(threads, maximum_query_bases, total_memory);
    let worker_memory = total_memory
        / u64::try_from(threads)
            .unwrap_or(u64::MAX)
            .saturating_mul(2)
            .max(1);
    let runner = crate::trace::runner::TraceRunnerConfig {
        sensitivity: sensitivity.clone(),
        candidates: crate::trace::screen::CandidateSearchConfig {
            min_shared_hashes: min_shared,
            min_plasmid_containment: 0.0,
            min_metagenome_containment: 0.0,
            top_candidates: candidate_limit,
        },
        resources: crate::resource::ResourceOpenOptions::default(),
        threads,
        io_concurrency: threads,
        max_alignments_per_candidate: max_alignments,
        query_kind,
        topology_requested: topology,
        topology_margin_bases: topology_margin.unwrap_or(sensitivity.auto_topology_margin_bases),
        memory_budget_bytes: worker_memory,
    };
    (
        sensitivity,
        crate::jam_index::JamIndexTraceConfig {
            screen: crate::jam_index::JamIndexScreenConfig {
                top_candidates: candidate_limit,
                accumulator_capacity: candidate_limit.saturating_mul(16).max(4_096),
                min_shared_hashes: min_shared,
                min_query_windows,
                rare_rescue_max_document_frequency: rare_rescue_df,
                parallel_parts: threads,
            },
            contigs: crate::jam_index::JamIndexContigSearchConfig {
                initial_contigs_per_candidate: initial_contigs,
                max_contigs_per_candidate: max_contigs,
                accumulator_capacity: max_contigs.saturating_mul(8).max(max_contigs),
                max_ranked_contig_bases: max_contig_bases,
                strong_candidate_shared_hashes: min_shared.max(4),
                min_spatial_signatures: 2,
                min_query_windows,
                rare_rescue_max_document_frequency: rare_rescue_df,
                whole_sample_min_shared_hashes: whole_sample_min_shared,
                parallel_parts: threads,
            },
            runner,
            expansion_batch,
            parallel_candidates: alignment_workers,
        },
    )
}

enum BatchTextWriter {
    Plain(std::io::BufWriter<std::fs::File>),
    Zstd(zstd::stream::write::Encoder<'static, std::io::BufWriter<std::fs::File>>),
}

impl BatchTextWriter {
    fn create(path: &Path) -> Result<Self> {
        let file = std::fs::OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(path)?;
        let writer = std::io::BufWriter::with_capacity(1024 * 1024, file);
        let compressed = path
            .extension()
            .and_then(|extension| extension.to_str())
            .is_some_and(|extension| {
                extension.eq_ignore_ascii_case("zst") || extension.eq_ignore_ascii_case("zstd")
            });
        if compressed {
            Ok(Self::Zstd(zstd::stream::write::Encoder::new(writer, 0)?))
        } else {
            Ok(Self::Plain(writer))
        }
    }

    fn write_line(&mut self, fields: &[&str]) -> Result<()> {
        let line = fields
            .iter()
            .map(|field| {
                field
                    .replace('\\', "\\\\")
                    .replace('\t', "\\t")
                    .replace('\r', "\\r")
                    .replace('\n', "\\n")
            })
            .collect::<Vec<_>>()
            .join("\t");
        match self {
            Self::Plain(writer) => writeln!(writer, "{line}")?,
            Self::Zstd(writer) => writeln!(writer, "{line}")?,
        }
        Ok(())
    }

    fn finish(self) -> Result<()> {
        match self {
            Self::Plain(mut writer) => {
                writer.flush()?;
                writer.get_ref().sync_all()?;
            }
            Self::Zstd(writer) => {
                let mut writer = writer.finish()?;
                writer.flush()?;
                writer.get_ref().sync_all()?;
            }
        }
        Ok(())
    }
}

fn index_sensitivity(
    profile: crate::trace::config::SensitivityProfile,
) -> crate::trace::config::SensitivityConfig {
    let mut sensitivity = crate::trace::config::SensitivityConfig::for_profile(profile);
    sensitivity.primary.scale = 1;
    if let Some(rescue) = sensitivity.rescue.as_mut() {
        rescue.scale = 1;
    }
    sensitivity.gap_rescue.dense_primary = None;
    sensitivity
}

fn index_query(path: &Path) -> Result<(String, Vec<u8>)> {
    let mut reader = needletail::parse_fastx_file(path)?;
    let record = reader
        .next()
        .transpose()?
        .ok_or_else(|| anyhow::anyhow!("Query input contains no FASTA/FASTQ records"))?;
    let id = std::str::from_utf8(record.id())?.to_string();
    let sequence = record.normalize(true).into_owned();
    if let Some(record) = reader.next() {
        record?;
        return Err(anyhow::anyhow!(
            "Query input must contain exactly one FASTA/FASTQ record"
        ));
    }
    Ok((id, sequence))
}

fn index_workers(threads: usize, query_length: usize, memory_bytes: u64) -> usize {
    const TARGET_BYTES: u64 = 512 * 1024 * 1024;
    const FIXED_BYTES: u64 = 192 * 1024 * 1024;
    const CELL_BYTES: u64 = 13;
    const ROW_CELLS: u64 = 129;
    const SINGLE_WORKER_BASES: usize = 95_000;

    if query_length >= SINGLE_WORKER_BASES {
        return 1;
    }
    let target = memory_bytes.min(TARGET_BYTES);
    let available = target.saturating_sub(FIXED_BYTES.min(target));
    let worker = u64::try_from(query_length)
        .unwrap_or(u64::MAX)
        .saturating_mul(CELL_BYTES)
        .saturating_mul(ROW_CELLS)
        .max(1);
    let admitted = usize::try_from(available / worker).unwrap_or(usize::MAX);
    threads.min(admitted.max(1)).max(1)
}

fn parse_sha256(value: &str) -> Result<[u8; 32]> {
    if value.len() != 64 {
        return Err(anyhow::anyhow!(
            "SHA-256 must contain 64 hexadecimal digits"
        ));
    }
    let mut bytes = [0u8; 32];
    let (chunks, remainder) = value.as_bytes().as_chunks::<2>();
    debug_assert!(remainder.is_empty());
    for (index, chunk) in chunks.iter().enumerate() {
        let text = std::str::from_utf8(chunk)?;
        bytes[index] = u8::from_str_radix(text, 16)
            .map_err(|_| anyhow::anyhow!("SHA-256 contains a non-hexadecimal digit"))?;
    }
    Ok(bytes)
}

#[allow(clippy::too_many_arguments)]
fn route_query_with_router(
    router_path: &std::path::Path,
    catalog: &crate::trace::catalog::TraceCatalog,
    sequence: &[u8],
    min_trace_length: u64,
    target_identity: f64,
    max_zero_witness_risk: f64,
    strict_witness_risk: bool,
    query_window_size: u32,
    handoff_mode: crate::router::WitnessHandoffMode,
    top_candidates: usize,
) -> Result<(
    Vec<crate::router::RoutedCandidate>,
    crate::trace::model::InputResource,
)> {
    let reader = crate::router::reader::RouterReader::open_mmap(router_path)?;
    let expected_metagenomes = catalog
        .entries()
        .iter()
        .map(|entry| {
            Ok(crate::router::format::MetagenomeEntry {
                id: entry.metagenome_id.clone(),
                object_uri: entry.resource_uri.clone(),
                checksum: parse_sha256(&entry.sha256)?,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    if reader.metagenomes().len() != expected_metagenomes.len()
        || reader
            .metagenomes()
            .iter()
            .zip(&expected_metagenomes)
            .any(|(bound, expected)| bound.id != expected.id || bound.checksum != expected.checksum)
    {
        return Err(anyhow::anyhow!(
            "router catalog binding does not match --metagenomes"
        ));
    }
    let scheme = reader
        .schemes()
        .iter()
        .find(|scheme| scheme.k == crate::router::WITNESS_K)
        .cloned()
        .ok_or_else(|| anyhow::anyhow!("router has no k=21 witness scheme"))?;
    let plan = crate::trace::query_plan::plan_witness_tier(
        crate::trace::WitnessPlanningRequest {
            min_trace_length,
            target_identity,
            max_zero_witness_risk,
            strict: strict_witness_risk,
        },
        &scheme,
    )?;
    if let Some(warning) = &plan.warning {
        eprintln!("warning: {warning}");
    }
    let nested =
        crate::router::witness::extract_nested_witnesses(sequence, &scheme, query_window_size)?;
    let mut scales = plan
        .selected_witness_tier
        .into_iter()
        .chain(plan.available_denser_tiers.iter().copied())
        .collect::<Vec<_>>();
    if scales.is_empty() {
        scales = scheme.available_scales.clone();
    }
    scales.sort_unstable_by(|left, right| right.cmp(left));
    scales.dedup();
    let mut witnesses = Vec::new();
    for witness in &nested.witnesses {
        for &scale in &scales {
            if witness.retained_at(scale) {
                witnesses.push(crate::router::search::TieredQueryWitness::new(
                    scale,
                    witness.query_witness(),
                ));
            }
        }
    }
    let max_shared_witnesses = nested.witnesses.len().clamp(1, 65_536);
    let source = crate::router::source::RouterPostingSource::new(reader, scheme)?;
    let router = crate::router::search::CandidateRouter::new(
        &source,
        crate::router::search::RouterSearchConfig {
            top_k: top_candidates,
            accumulator_capacity: catalog.len().max(top_candidates),
            max_shared_witnesses_per_candidate: max_shared_witnesses,
            max_positional_witnesses_per_candidate: max_shared_witnesses.min(4096),
            handoff_mode,
            max_positions_per_witness: 4,
            total_query_windows: Some(nested.eligible_window_count()),
            ..crate::router::search::RouterSearchConfig::default()
        },
    );
    let candidates = router.search(&witnesses)?;
    let sha256 = provenance::sha256_file(router_path)?;
    Ok((
        candidates,
        crate::trace::model::InputResource {
            role: "router".to_string(),
            redacted_locator: crate::resource::ResourceLocator::parse(
                router_path.to_string_lossy().as_ref(),
            )?
            .redacted(),
            sha256: Some(sha256),
        },
    ))
}

#[allow(clippy::too_many_arguments)]
pub fn handle_trace_command(
    query_path: PathBuf,
    query_kind: crate::trace::model::QueryKind,
    topology_requested: crate::trace::model::TopologyRequested,
    database: Option<String>,
    router: Option<PathBuf>,
    metagenomes_path: PathBuf,
    output: PathBuf,
    upload_to: Option<String>,
    query_id: Option<String>,
    profile: crate::trace::config::SensitivityProfile,
    min_shared_hashes: u32,
    min_query_containment: f64,
    min_metagenome_containment: f64,
    top_candidates: Option<usize>,
    min_trace_length: u64,
    target_identity: f64,
    max_zero_witness_risk: f64,
    strict_witness_risk: bool,
    query_window_size: u32,
    router_handoff: crate::router::WitnessHandoffMode,
    max_alignments: usize,
    threads: usize,
    io_concurrency: usize,
    topology_margin_bases: Option<u64>,
    memory_gb: usize,
    cache_dir: Option<PathBuf>,
    cache_block_bytes: u64,
    request_timeout_seconds: u64,
    max_retries: u32,
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
        expected_sha256: None,
        cache_block_bytes,
        max_cache_bytes: (memory_gb as u64)
            .saturating_mul(1024 * 1024 * 1024)
            .saturating_mul(7)
            / 10,
        request_timeout_seconds,
        max_retries,
        allow_full_download_fallback: false,
    };
    let catalog = crate::trace::catalog::TraceCatalog::from_path(&metagenomes_path)?;
    let parsed_query =
        crate::trace::raw::RawAssembly::open(query_path.to_string_lossy(), resources.clone())?;
    if parsed_query.contigs.len() != 1 {
        return Err(anyhow::anyhow!(
            "Query input must contain exactly one FASTA/FASTQ record, got {}",
            parsed_query.contigs.len()
        ));
    }
    let record = &parsed_query.contigs[0];
    let query_id = query_id.unwrap_or_else(|| record.id.clone());
    let sensitivity = crate::trace::config::SensitivityConfig::for_profile(profile);
    let topology_margin_bases =
        topology_margin_bases.unwrap_or(sensitivity.auto_topology_margin_bases);
    let candidate_limit = top_candidates
        .unwrap_or_else(|| usize::try_from(sensitivity.max_candidates).unwrap_or(usize::MAX));
    let runner = crate::trace::runner::TraceRunner::new(crate::trace::runner::TraceRunnerConfig {
        sensitivity: sensitivity.clone(),
        candidates: crate::trace::screen::CandidateSearchConfig {
            min_shared_hashes,
            min_plasmid_containment: min_query_containment,
            min_metagenome_containment,
            top_candidates: candidate_limit,
        },
        resources,
        threads,
        io_concurrency,
        max_alignments_per_candidate: max_alignments,
        query_kind,
        topology_requested,
        topology_margin_bases,
        memory_budget_bytes: (memory_gb as u64)
            .saturating_mul(1024 * 1024 * 1024)
            .saturating_mul(7)
            / 10,
    })?;
    let candidate_locator = database
        .clone()
        .or_else(|| {
            router
                .as_ref()
                .map(|path| path.to_string_lossy().into_owned())
        })
        .ok_or_else(|| anyhow::anyhow!("either --database or --router is required"))?;
    let query = crate::trace::runner::TraceQuery {
        plasmid_id: query_id.clone(),
        plasmid_sequence: record.sequence.clone(),
        database: candidate_locator,
        catalog: catalog.clone(),
    };
    let started = provenance::unix_time_seconds();
    let run_id = format!("trace-{started}-{}", std::process::id());
    let mut result = if let Some(router_path) = router.as_deref() {
        let (candidates, router_input) = route_query_with_router(
            router_path,
            &catalog,
            &record.sequence,
            min_trace_length,
            target_identity,
            max_zero_witness_risk,
            strict_witness_risk,
            query_window_size,
            router_handoff,
            candidate_limit,
        )?;
        runner.run_routed(
            &query,
            candidates,
            router_input,
            crate::resource::ResourceMetrics::default(),
        )?
    } else {
        runner.run(&query)?
    };
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
        command: provenance::redacted_command_line(),
        plasmid_id: query_id,
        plasmid_length: record.sequence.len() as u64,
        query_kind,
        topology_requested,
        threads,
        io_concurrency,
        sensitivity,
        algorithms: crate::trace::config::algorithm_identifiers(),
        algorithm,
        inputs: vec![
            trace_input("query", &query_path)?,
            result.database_input.clone(),
            trace_input("metagenomes", &metagenomes_path)?,
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
    let upload = upload_to
        .as_deref()
        .map(|locator| {
            crate::resource::upload::upload_file(
                locator,
                &output,
                crate::resource::upload::UploadOptions {
                    request_timeout_seconds,
                    max_retries,
                },
            )
        })
        .transpose()?;
    if !silent {
        eprintln!(
            "Trace output {}: {} candidates, {} aligned metagenomes, {} alignments",
            output.display(),
            result.search.candidates.len(),
            footer.metagenomes_aligned,
            footer.alignments_total
        );
        if let (Some(locator), Some(upload)) = (upload_to.as_deref(), upload.as_ref()) {
            let redacted = crate::resource::ResourceLocator::parse(locator)?.redacted();
            eprintln!(
                "Uploaded {} bytes to {} in {} attempt(s), HTTP status {}",
                upload.bytes_uploaded, redacted, upload.attempts, upload.status
            );
        }
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
    use super::{compute_distance_chunk_size, index_workers, normalize_distance_cutoff};

    #[test]
    fn index_workers_bound() {
        let memory = 2 * 1024 * 1024 * 1024;
        assert_eq!(index_workers(4, 10_000, memory), 4);
        assert_eq!(index_workers(4, 94_281, memory), 2);
        assert_eq!(index_workers(4, 97_566, memory), 1);
        assert_eq!(index_workers(4, 168_903, memory), 1);
    }

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
