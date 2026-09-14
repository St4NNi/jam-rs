use anyhow::Result;
use indicatif::{ProgressBar, ProgressStyle};
use needletail::parse_fastx_file;
use std::fs::remove_file;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::bias::{BiasCreateConfig, CMSConfig, HashBiasTable};
use crate::collection::CollectionTraceEngine;
use crate::jidx_builder::{JidxBuildConfig, build_local_jidx};
use crate::jidx_writer::sync_directory;
use crate::query::QueryEngine;
use crate::range_source::S3Config;
use crate::reader::JamReader;
use crate::trace::{TraceConfig, TraceEngine};
use crate::writer::{BuildConfig, build};
use std::sync::Arc;

pub(crate) fn handle_jidx_build_command(
    database: PathBuf,
    manifest: PathBuf,
    output: PathBuf,
    config: JidxBuildConfig,
    force: bool,
) -> Result<()> {
    if output.try_exists()? {
        if !output.is_file() {
            return Err(anyhow::anyhow!("Output path is not a file: {:?}", output));
        }
        if !force {
            return Err(anyhow::anyhow!(
                "Output file {:?} already exists. Use --force to overwrite.",
                output
            ));
        }
    }
    let parent = output
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let staged = tempfile::Builder::new()
        .prefix(".jam-jidx-")
        .tempfile_in(parent)?
        .into_temp_path();
    remove_file(&staged)?;
    build_local_jidx(database, manifest, &staged, config)?;
    if force {
        staged.persist(&output).map_err(|error| error.error)?;
    } else {
        staged
            .persist_noclobber(&output)
            .map_err(|error| error.error)?;
    }
    sync_directory(parent)?;
    Ok(())
}

pub(crate) enum TraceInput {
    Shard {
        database: PathBuf,
        index: PathBuf,
        manifest: PathBuf,
    },
    Collection(PathBuf),
    Owner(PathBuf),
    Shared {
        path: PathBuf,
        read_stats: Option<PathBuf>,
        query_topology_header: bool,
    },
}

pub(crate) struct TraceArgs {
    pub query: PathBuf,
    pub input: TraceInput,
    pub audit_index: bool,
    pub output: PathBuf,
    pub query_id: Option<String>,
    pub config: TraceConfig,
    pub s3: Option<S3Config>,
    pub force: bool,
}

pub(crate) const SHARED_BATCH_QUERY_BASES: usize = 640_000;

pub(crate) fn handle_trace_command(args: TraceArgs) -> Result<()> {
    use crate::trace_batch::{phase_cost_enabled, phase_elapsed, phase_stamp};
    let invocation_started = std::time::Instant::now();
    let invocation_phase = phase_stamp();
    let _ = crate::trace_batch::timeline_tick();
    let shared_input = matches!(&args.input, TraceInput::Shared { .. });
    let query_topology_header = matches!(
        &args.input,
        TraceInput::Shared {
            query_topology_header: true,
            ..
        }
    );
    let read_stats = match &args.input {
        TraceInput::Shared { read_stats, .. } => read_stats.clone(),
        _ => None,
    };
    if let Some(path) = &read_stats
        && (path == &args.output || path.try_exists()?)
    {
        return Err(anyhow::anyhow!(
            "Read-statistics output must be a new distinct file"
        ));
    }
    if args.output.try_exists()? {
        if !args.output.is_file() {
            return Err(anyhow::anyhow!(
                "Output path is not a file: {:?}",
                args.output
            ));
        }
        if !args.force {
            return Err(anyhow::anyhow!(
                "Output file {:?} already exists. Use --force to overwrite.",
                args.output
            ));
        }
    }
    let parent = args
        .output
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    if args.output.file_name().is_none() {
        return Err(anyhow::anyhow!("Invalid output path: {:?}", args.output));
    }

    enum Engine {
        Shard(Box<TraceEngine>),
        Collection(Box<CollectionTraceEngine>),
    }
    let engine = match args.input {
        TraceInput::Shard {
            database,
            index,
            manifest,
        } => Engine::Shard(Box::new(TraceEngine::open(
            database, index, manifest, args.s3,
        )?)),
        TraceInput::Collection(root) => {
            Engine::Collection(Box::new(CollectionTraceEngine::open(root, args.s3)?))
        }
        TraceInput::Owner(root) => Engine::Shard(Box::new(TraceEngine::open_owner(root, args.s3)?)),
        TraceInput::Shared {
            path, read_stats, ..
        } => Engine::Shard(Box::new(TraceEngine::open_shared_observed(
            path,
            args.s3,
            read_stats.is_some() && !phase_cost_enabled(),
        )?)),
    };
    if args.audit_index {
        match &engine {
            Engine::Shard(engine) => engine.verify_index()?,
            Engine::Collection(engine) => engine.verify_index()?,
        }
    }
    let batch_size = match &engine {
        Engine::Shard(_) if shared_input => 64,
        Engine::Shard(_) => rayon::current_num_threads(),
        Engine::Collection(_) => rayon::current_num_threads().clamp(1, 4),
    };
    let setup_cpu_ns = phase_elapsed(invocation_phase)[1];
    let mut parsing_cpu_ns = 0;
    let mut search_cpu_ns = 0;
    let mut output_cpu_ns = 0;
    let parsing_phase = phase_stamp();
    let startup_ended = std::time::Instant::now();
    let startup_ns = startup_ended.duration_since(invocation_started).as_nanos() as u64;
    let parsing_started = read_stats.as_ref().map(|_| startup_ended);
    let mut input = parse_fastx_file(&args.query)?;
    let mut parsing_ns = parsing_started.map_or(0, |started| started.elapsed().as_nanos() as u64);
    parsing_cpu_ns += phase_elapsed(parsing_phase)[1];
    let mut output_ns = 0u64;
    let mut search_ns = 0u64;
    let mut temporary = tempfile::Builder::new()
        .prefix(".jam-trace-")
        .tempfile_in(parent)?;
    let mut count = 0usize;
    let mut pending: Option<(String, Vec<u8>, bool)> = None;
    {
        let mut output = BufWriter::new(temporary.as_file_mut());
        loop {
            let parsing_phase = phase_stamp();
            let parsing_started = read_stats.as_ref().map(|_| std::time::Instant::now());
            let mut queries = Vec::with_capacity(batch_size);
            let mut topologies = Vec::with_capacity(batch_size);
            let mut batch_bases = 0usize;
            for _ in 0..batch_size {
                if let Some((id, sequence, circular)) = pending.take() {
                    batch_bases += sequence.len();
                    queries.push((id, sequence));
                    topologies.push(circular);
                    continue;
                }
                let Some(record) = input.next() else { break };
                let record = record?;
                if count != 0 && args.query_id.is_some() {
                    return Err(anyhow::anyhow!(
                        "--query-id requires a query file with one record"
                    ));
                }
                let id = match &args.query_id {
                    Some(id) => id.clone(),
                    None => std::str::from_utf8(
                        record
                            .id()
                            .split(|byte| byte.is_ascii_whitespace())
                            .next()
                            .unwrap_or_default(),
                    )?
                    .to_string(),
                };
                let sequence = record.seq().into_owned();
                let circular = if query_topology_header {
                    let mut tokens = record
                        .id()
                        .split(|byte| byte.is_ascii_whitespace())
                        .skip(1)
                        .filter(|token| token.starts_with(b"topology="));
                    let circular = match tokens.next() {
                        Some(b"topology=linear") => false,
                        Some(b"topology=circular") => true,
                        _ => {
                            return Err(anyhow::anyhow!(
                                "Each query header requires topology=linear or topology=circular"
                            ));
                        }
                    };
                    if tokens.next().is_some() {
                        return Err(anyhow::anyhow!("Duplicate query topology token"));
                    }
                    circular
                } else {
                    args.config.circular
                };
                count += 1;
                if shared_input
                    && !queries.is_empty()
                    && batch_bases.saturating_add(sequence.len()) > SHARED_BATCH_QUERY_BASES
                {
                    pending = Some((id, sequence, circular));
                    break;
                }
                batch_bases += sequence.len();
                queries.push((id, sequence));
                topologies.push(circular);
            }
            parsing_ns += parsing_started.map_or(0, |started| started.elapsed().as_nanos() as u64);
            parsing_cpu_ns += phase_elapsed(parsing_phase)[1];
            if queries.is_empty() {
                break;
            }
            match &engine {
                Engine::Shard(engine) => {
                    let search_phase = phase_stamp();
                    let search_started = read_stats.as_ref().map(|_| std::time::Instant::now());
                    let results = if query_topology_header {
                        engine.search_batch_topologies(&queries, args.config, &topologies)?
                    } else {
                        engine.search_batch(&queries, args.config)?
                    };
                    search_ns +=
                        search_started.map_or(0, |started| started.elapsed().as_nanos() as u64);
                    search_cpu_ns += phase_elapsed(search_phase)[1];
                    let output_phase = phase_stamp();
                    let output_started = read_stats.as_ref().map(|_| std::time::Instant::now());
                    for result in results {
                        serde_json::to_writer(&mut output, &result)?;
                        output.write_all(b"\n")?;
                    }
                    output_ns +=
                        output_started.map_or(0, |started| started.elapsed().as_nanos() as u64);
                    output_cpu_ns += phase_elapsed(output_phase)[1];
                }
                Engine::Collection(engine) => {
                    for result in engine.search_batch(&queries, args.config)? {
                        serde_json::to_writer(&mut output, &result)?;
                        output.write_all(b"\n")?;
                    }
                }
            }
        }
        output.flush()?;
    }
    if count == 0 {
        return Err(anyhow::anyhow!("Query file contains no sequence records"));
    }
    temporary.as_file().sync_all()?;
    if args.force {
        temporary
            .persist(&args.output)
            .map_err(|error| error.error)?;
    } else {
        temporary
            .persist_noclobber(&args.output)
            .map_err(|error| error.error)?;
    }
    sync_directory(parent)?;
    let publication_ns = invocation_started.elapsed().as_nanos() as u64;
    let publication_cpu_ns = phase_elapsed(invocation_phase)[1];
    if let Some(path) = read_stats {
        let Engine::Shard(engine) = &engine else {
            unreachable!()
        };
        let parent = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));
        let mut stats = tempfile::Builder::new()
            .prefix(".jam-read-stats-")
            .tempfile_in(parent)?;
        serde_json::to_writer_pretty(
            stats.as_file_mut(),
            &serde_json::json!({
                "format": "jam-shared-read-stats-v1", "index": engine.shared_read_stats(), "batch": engine.batch_stats(),
                "parsing_ns": parsing_ns, "output_ns": output_ns,
                "phase_cost_mode": phase_cost_enabled(),
                "worker_timeline": crate::trace_batch::worker_intervals(),
                "phase_cost_semantics": "phase arrays are [elapsed_ns, process_cpu_ns] at disjoint batch barriers; core lookup includes token materialization; postings is a nested contiguous subinterval of context lookup and must be subtracted for an exclusive table; downstream includes fallback postings, candidates, regions, sequence access, alignment and results; normalization is in extraction; native mode has no phase samples",
                "top_level_cpu_ns": { "setup": setup_cpu_ns, "parsing": parsing_cpu_ns, "search": search_cpu_ns, "result_serialization": output_cpu_ns, "finalization_and_other": publication_cpu_ns.saturating_sub(setup_cpu_ns + parsing_cpu_ns + search_cpu_ns + output_cpu_ns), "through_result_publication": publication_cpu_ns },
                "top_level_ns": { "setup": startup_ns, "parsing": parsing_ns, "search": search_ns, "result_serialization": output_ns, "finalization_and_other": publication_ns.saturating_sub(startup_ns + parsing_ns + search_ns + output_ns), "through_result_publication": publication_ns },
                "batch_limits": { "queries": 64, "query_bases": SHARED_BATCH_QUERY_BASES, "lookup_bytes": 134217728, "global_lookup_bytes": 268435456, "decoded_bgzf_bytes": 33554432, "concurrent_bgzf_reads": 4 },
                "histogram_semantics": "16 base-2 bins starting at one, last bin at least 32768; reuse counts per-query distinct requests for an exact context; occurrence fanout counts positions per admitted positive context; totals are per lookup batch",
                "lookup_task_semantics": "core-aligned tasks weighted by estimated directory comparisons plus requested contexts; dispatch is task construction and ordering; parallel is pool elapsed; compute and reduce sum worker spans; dispatch_to_start sums task queue delay, not worker idle time; key_lookup includes final entry ordering",
                "semantics": "top_level_ns intervals do not overlap and cover handler entry through result publication, excluding stats publication and CLI startup; batch timings are nested worker spans and may overlap; logical bytes are not physical I/O; endpoint time includes endpoint traceback"
            }),
        )?;
        stats.as_file().sync_all()?;
        stats
            .persist_noclobber(&path)
            .map_err(|error| error.error)?;
        sync_directory(parent)?;
    }
    Ok(())
}

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
                "      Bias mode reports containment on the retained/weighted k-mer subset; E-values are uniform-hash approximations",
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
    if has_bias {
        writeln!(
            writer,
            "query\tdb_sample\tshared_hashes\tquery_hashes\tdb_hashes\traw_query_containment\tdb_containment_unweighted\tuniform_hash_e_value\tbias_weighted_query_containment"
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
                matches.sort_by(|a, b| b.containment.total_cmp(&a.containment));

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
                    if has_bias && result.total_query_weight > 0.0 {
                        let raw_query_containment = if query_hashes > 0 {
                            m.hit_count as f64 / query_hashes as f64
                        } else {
                            0.0
                        };
                        let _ = writeln!(
                            out,
                            "{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6e}\t{:.6}",
                            query_name,
                            db_name,
                            m.hit_count,
                            query_hashes,
                            db_hashes,
                            raw_query_containment,
                            db_containment,
                            m.e_value,
                            m.containment,
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
    if !(1..=31).contains(&kmer_size) {
        return Err(anyhow::anyhow!(
            "K-mer size must be between 1 and 31, got {}",
            kmer_size
        ));
    }
    if fscale == 0 {
        return Err(anyhow::anyhow!("--fscale must be > 0"));
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

#[cfg(test)]
mod tests {
    use super::{
        TraceArgs, TraceInput, compute_distance_chunk_size, handle_jidx_build_command,
        handle_trace_command, normalize_distance_cutoff,
    };
    use crate::jidx_builder::JidxBuildConfig;
    use crate::trace::TraceConfig;

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

    #[test]
    fn failed_forced_jidx_build_preserves_existing_output() {
        let directory = tempfile::tempdir().unwrap();
        let output = directory.path().join("existing.jidx");
        std::fs::write(&output, b"old index").unwrap();
        assert!(
            handle_jidx_build_command(
                directory.path().join("missing.jam"),
                directory.path().join("missing.json"),
                output.clone(),
                JidxBuildConfig::default(),
                true,
            )
            .is_err()
        );
        assert_eq!(std::fs::read(output).unwrap(), b"old index");
    }

    #[test]
    fn failed_forced_trace_preserves_existing_output() {
        let directory = tempfile::tempdir().unwrap();
        let output = directory.path().join("existing.jsonl");
        std::fs::write(&output, b"old result").unwrap();
        assert!(
            handle_trace_command(TraceArgs {
                query: directory.path().join("missing.fa"),
                input: TraceInput::Shard {
                    database: directory.path().join("missing.jam"),
                    index: directory.path().join("missing.jidx"),
                    manifest: directory.path().join("missing.json"),
                },
                audit_index: false,
                output: output.clone(),
                query_id: None,
                config: TraceConfig::default(),
                s3: None,
                force: true,
            })
            .is_err()
        );
        assert_eq!(std::fs::read(output).unwrap(), b"old result");
    }
}
