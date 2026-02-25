use anyhow::{bail, Context, Result};
use bytemuck::cast_slice;
use clap::{Parser, Subcommand};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::fs::{create_dir_all, remove_dir_all, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::process::Command as ProcessCommand;
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};
use std::sync::Arc;
use std::thread;
use std::time::Instant;

#[cfg(not(target_endian = "little"))]
compile_error!("k21_pipeline requires little-endian targets");

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

const KMER_SIZE: u8 = 21;
const RAW_SPACE_SIZE: u64 = 1u64 << (2 * KMER_SIZE as u32);
const CANONICAL_SPACE_SIZE: u64 = RAW_SPACE_SIZE / 2;

const MAP_TASKS: usize = 64;
const PARTITIONS: usize = 64;
const PARTITION_BITS: u32 = 6;

const CHI_BUCKETS: usize = 2048;
const CHI_BUCKET_BITS: u32 = 11;

const IO_BUFFER_BYTES: usize = 64 * 1024 * 1024;
const DEFAULT_MAP_FILE_BUFFER_BYTES: usize = 64 * 1024 * 1024;
const DEFAULT_LOCAL_BUFFER_VALUES: usize = 32 * 1024;
const DEFAULT_MAP_IO_THREADS: usize = 16;
const DEFAULT_MAP_CHANNEL_ENTRIES: usize = 4_000_000;
const DEFAULT_LOCAL_SCRATCH_BASE: &str = "/var/scratch/sbeyvers/k21";

const FRAC_SCALES: [u64; 4] = [10, 100, 1000, 10000];

const HASH_NAMES: [&str; 4] = ["jamhash", "murmur3", "xxhash3", "wang64"];
const HASH_COUNT: usize = HASH_NAMES.len();

const U64_SPACE_F64: f64 = 18_446_744_073_709_551_616.0;

#[derive(Parser)]
#[command(name = "k21_pipeline")]
#[command(about = "Distributed k=21 canonical hash analysis pipeline")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Map stage: hash canonical k=21 space shard for this array task
    Map {
        /// Output directory root
        output_dir: PathBuf,
    },
    /// Prepare stage: combine map summaries into partition plan
    Prepare {
        /// Output directory root
        output_dir: PathBuf,
    },
    /// Sort+reduce stage: one (hash, partition) per array task
    SortReduce {
        /// Output directory root
        output_dir: PathBuf,
    },
    /// Aggregate stage per hash (array task 0..3)
    AggregateHash {
        /// Output directory root
        output_dir: PathBuf,
    },
    /// Final aggregate across all 4 hashes
    AggregateAll {
        /// Output directory root
        output_dir: PathBuf,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct MapTaskSummary {
    task_id: usize,
    raw_start: u64,
    raw_end: u64,
    canonical_count: u64,
    partition_counts: Vec<Vec<u64>>, // [hash][partition]
    elapsed_seconds: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PartitionPlan {
    map_tasks: usize,
    partitions: usize,
    chi_buckets: usize,
    total_canonical: u64,
    hash_partition_counts: Vec<Vec<u64>>,  // [hash][partition]
    hash_partition_offsets: Vec<Vec<u64>>, // [hash][partition], global rank start
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PartitionStats {
    hash_name: String,
    partition_id: usize,
    count: u64,
    min_hash: u64,
    max_hash: u64,
    observed_collision_pairs: u64,
    duplicate_hashes: u64,
    collision_elements: u64,
    max_multiplicity: u64,
    ks_statistic_max: f64,
    frac_pass_counts: Vec<u64>, // len FRAC_SCALES
    chi_bucket_counts: Vec<u64>,
    byte_hist_counts: Vec<u64>, // flattened [8][256]
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct HashMetrics {
    hash_name: String,
    total_hashes: u64,
    expected_hashes: u64,
    min_hash: u64,
    max_hash: u64,
    frac_scales: Vec<u64>,
    frac_pass_counts: Vec<u64>,
    frac_observed: Vec<f64>,
    frac_expected: Vec<f64>,
    frac_abs_error: Vec<f64>,
    frac_max_abs_error: f64,
    bit_max_abs_error: f64,
    bit_worst_index: usize,
    chi_squared_statistic: f64,
    chi_max_rel_error: f64,
    ks_statistic: f64,
    observed_collision_pairs: u64,
    expected_collision_pairs: f64,
    collision_ratio_obs_exp: f64,
    duplicate_hashes: u64,
    collision_elements: u64,
    max_multiplicity: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct CombinedReport {
    hashes: Vec<HashMetrics>,
}

#[derive(Clone)]
struct WorkerSummary {
    canonical_count: u64,
    partition_counts: [[u64; PARTITIONS]; HASH_COUNT],
}

struct WriteChunk {
    stream_idx: usize,
    values: Vec<u64>,
}

struct ScratchCleanupGuard {
    path: PathBuf,
    active: bool,
}

impl ScratchCleanupGuard {
    fn new(path: PathBuf) -> Self {
        Self { path, active: true }
    }

    fn disarm(&mut self) {
        self.active = false;
    }
}

impl Drop for ScratchCleanupGuard {
    fn drop(&mut self) {
        if !self.active {
            return;
        }

        if let Err(err) = remove_dir_all(&self.path) {
            if err.kind() != std::io::ErrorKind::NotFound {
                eprintln!(
                    "[map] warning: failed to cleanup scratch {}: {}",
                    self.path.display(),
                    err
                );
            }
        }
    }
}

impl Default for WorkerSummary {
    fn default() -> Self {
        Self {
            canonical_count: 0,
            partition_counts: [[0u64; PARTITIONS]; HASH_COUNT],
        }
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Map { output_dir } => cmd_map(&output_dir),
        Command::Prepare { output_dir } => cmd_prepare(&output_dir),
        Command::SortReduce { output_dir } => cmd_sort_reduce(&output_dir),
        Command::AggregateHash { output_dir } => cmd_aggregate_hash(&output_dir),
        Command::AggregateAll { output_dir } => cmd_aggregate_all(&output_dir),
    }
}

fn cmd_map(output_dir: &Path) -> Result<()> {
    ensure_layout(output_dir)?;
    let task_id = slurm_array_task_id_or_default(0)?;
    if task_id >= MAP_TASKS {
        bail!("map task id {} out of range 0..{}", task_id, MAP_TASKS - 1);
    }

    let total_threads = configured_threads();
    let io_threads = configured_map_io_threads(total_threads);
    let compute_threads = total_threads.saturating_sub(io_threads).max(1);
    let local_buffer_values = configured_local_buffer_values();
    let channel_entries = configured_map_channel_entries(local_buffer_values);
    let channel_capacity = div_ceil(channel_entries, local_buffer_values).max(1);
    let map_file_buffer_bytes = configured_map_file_buffer_bytes();
    let local_scratch_base = configured_local_scratch_base();
    let local_task_scratch = map_task_scratch_dir(&local_scratch_base, task_id);

    if local_task_scratch.exists() {
        remove_dir_all(&local_task_scratch).with_context(|| {
            format!(
                "failed to remove previous scratch {}",
                local_task_scratch.display()
            )
        })?;
    }
    create_dir_all(&local_task_scratch).with_context(|| {
        format!(
            "failed to create local scratch {}",
            local_task_scratch.display()
        )
    })?;
    let mut scratch_guard = ScratchCleanupGuard::new(local_task_scratch.clone());

    configure_rayon(compute_threads);

    let range = get_chunk_range(RAW_SPACE_SIZE, MAP_TASKS as u64, task_id as u64);
    let start = Instant::now();

    eprintln!(
        "[map:{}] range={}..{} compute_threads={} io_threads={} local_buf_values={} file_buf_mb={} channel_entries={} channel_chunks={} partitions={} scratch={}",
        task_id,
        range.start,
        range.end,
        compute_threads,
        io_threads,
        local_buffer_values,
        map_file_buffer_bytes / (1024 * 1024),
        channel_entries,
        channel_capacity,
        PARTITIONS,
        local_task_scratch.display()
    );

    let (io_senders, io_handles) = start_map_io_workers(
        &local_task_scratch,
        task_id,
        io_threads,
        channel_capacity,
        map_file_buffer_bytes,
    )?;
    let io_senders = Arc::new(io_senders);

    let subranges = split_range(range.clone(), compute_threads.saturating_mul(2).max(1));

    let worker_results: Vec<Result<WorkerSummary>> = subranges
        .into_par_iter()
        .map(|sub| process_subrange(sub, io_senders.clone(), io_threads, local_buffer_values))
        .collect();

    let mut combined = WorkerSummary::default();
    for result in worker_results {
        let ws = result?;
        combined.canonical_count += ws.canonical_count;
        for h in 0..HASH_COUNT {
            for p in 0..PARTITIONS {
                combined.partition_counts[h][p] += ws.partition_counts[h][p];
            }
        }
    }

    drop(io_senders);

    for handle in io_handles {
        let worker_result = handle
            .join()
            .map_err(|_| anyhow::anyhow!("map io worker panicked"))?;
        worker_result?;
    }

    copy_map_shards_to_output(&local_task_scratch, output_dir, task_id)?;

    remove_dir_all(&local_task_scratch).with_context(|| {
        format!(
            "failed to cleanup local scratch {}",
            local_task_scratch.display()
        )
    })?;
    scratch_guard.disarm();

    let partition_counts = (0..HASH_COUNT)
        .map(|h| combined.partition_counts[h].to_vec())
        .collect::<Vec<_>>();

    let summary = MapTaskSummary {
        task_id,
        raw_start: range.start,
        raw_end: range.end,
        canonical_count: combined.canonical_count,
        partition_counts,
        elapsed_seconds: start.elapsed().as_secs_f64(),
    };

    let summary_path = map_summary_path(output_dir, task_id);
    write_json_pretty(&summary_path, &summary)?;

    eprintln!(
        "[map:{}] done canonical={} elapsed={:.2}s",
        task_id, summary.canonical_count, summary.elapsed_seconds
    );

    Ok(())
}

fn cmd_prepare(output_dir: &Path) -> Result<()> {
    ensure_layout(output_dir)?;
    let mut total_canonical = 0u64;
    let mut hash_partition_counts = vec![vec![0u64; PARTITIONS]; HASH_COUNT];

    for map_task in 0..MAP_TASKS {
        let summary_path = map_summary_path(output_dir, map_task);
        let summary: MapTaskSummary = read_json(&summary_path)
            .with_context(|| format!("failed to read {}", summary_path.display()))?;

        if summary.partition_counts.len() != HASH_COUNT {
            bail!(
                "invalid partition_counts hash dimension in {}",
                summary_path.display()
            );
        }

        total_canonical += summary.canonical_count;

        for (h, row) in summary.partition_counts.iter().enumerate() {
            if row.len() != PARTITIONS {
                bail!(
                    "invalid partition_counts partition dimension in {}",
                    summary_path.display()
                );
            }
            for (p, value) in row.iter().enumerate() {
                hash_partition_counts[h][p] += *value;
            }
        }
    }

    let mut hash_partition_offsets = vec![vec![0u64; PARTITIONS]; HASH_COUNT];
    for h in 0..HASH_COUNT {
        let mut prefix = 0u64;
        for p in 0..PARTITIONS {
            hash_partition_offsets[h][p] = prefix;
            prefix += hash_partition_counts[h][p];
        }
        if prefix != total_canonical {
            bail!(
                "hash {} total mismatch: partition sum {} != canonical {}",
                HASH_NAMES[h],
                prefix,
                total_canonical
            );
        }
    }

    if total_canonical != CANONICAL_SPACE_SIZE {
        eprintln!(
            "[prepare] warning: canonical total={} expected={}",
            total_canonical, CANONICAL_SPACE_SIZE
        );
    }

    let plan = PartitionPlan {
        map_tasks: MAP_TASKS,
        partitions: PARTITIONS,
        chi_buckets: CHI_BUCKETS,
        total_canonical,
        hash_partition_counts,
        hash_partition_offsets,
    };

    let path = partition_plan_path(output_dir);
    write_json_pretty(&path, &plan)?;
    eprintln!("[prepare] wrote {}", path.display());

    Ok(())
}

fn cmd_sort_reduce(output_dir: &Path) -> Result<()> {
    ensure_layout(output_dir)?;

    let task_id = slurm_array_task_id_or_default(0)?;
    let total_tasks = HASH_COUNT * PARTITIONS;
    if task_id >= total_tasks {
        bail!(
            "sort-reduce task id {} out of range 0..{}",
            task_id,
            total_tasks - 1
        );
    }

    let hash_idx = task_id / PARTITIONS;
    let partition_id = task_id % PARTITIONS;
    let hash_name = HASH_NAMES[hash_idx];

    let plan: PartitionPlan = read_json(&partition_plan_path(output_dir))?;
    validate_plan(&plan)?;

    let expected = plan.hash_partition_counts[hash_idx][partition_id];
    let global_offset = plan.hash_partition_offsets[hash_idx][partition_id];

    let mut values = Vec::with_capacity(expected as usize);

    let read_start = Instant::now();
    for map_task in 0..plan.map_tasks {
        let path = shard_path(output_dir, hash_name, map_task, partition_id);
        read_u64_file_into_vec(&path, &mut values)
            .with_context(|| format!("failed to read {}", path.display()))?;
    }
    let read_secs = read_start.elapsed().as_secs_f64();

    if values.len() as u64 != expected {
        bail!(
            "{} part {}: loaded {} hashes, expected {}",
            hash_name,
            partition_id,
            values.len(),
            expected
        );
    }

    let threads = configured_threads();
    configure_rayon(threads);

    let sort_start = Instant::now();
    values.par_sort_unstable();
    let sort_secs = sort_start.elapsed().as_secs_f64();

    let stats_start = Instant::now();
    let stats = compute_partition_stats(
        hash_name,
        partition_id,
        &values,
        global_offset,
        plan.total_canonical,
    )?;
    let stats_secs = stats_start.elapsed().as_secs_f64();

    let sorted_path = sorted_partition_path(output_dir, hash_name, partition_id);
    write_u64_file(&sorted_path, &values)?;

    let part_stats_path = partition_stats_path(output_dir, hash_name, partition_id);
    write_json_pretty(&part_stats_path, &stats)?;

    eprintln!(
        "[sort-reduce:{}:{}] n={} read={:.2}s sort={:.2}s stats={:.2}s",
        hash_name,
        partition_id,
        values.len(),
        read_secs,
        sort_secs,
        stats_secs
    );

    Ok(())
}

fn cmd_aggregate_hash(output_dir: &Path) -> Result<()> {
    ensure_layout(output_dir)?;
    let task_id = slurm_array_task_id_or_default(0)?;
    if task_id >= HASH_COUNT {
        bail!(
            "aggregate-hash task id {} out of range 0..{}",
            task_id,
            HASH_COUNT - 1
        );
    }

    let hash_name = HASH_NAMES[task_id];
    let plan: PartitionPlan = read_json(&partition_plan_path(output_dir))?;
    validate_plan(&plan)?;

    let mut total_hashes = 0u64;
    let mut min_hash = u64::MAX;
    let mut max_hash = 0u64;

    let mut observed_pairs = 0u64;
    let mut duplicate_hashes = 0u64;
    let mut collision_elements = 0u64;
    let mut max_multiplicity = 0u64;
    let mut ks_statistic = 0.0f64;

    let mut frac_pass_counts = vec![0u64; FRAC_SCALES.len()];
    let mut chi_bucket_counts = vec![0u64; CHI_BUCKETS];
    let mut byte_hist_counts = vec![0u64; 8 * 256];

    for partition_id in 0..PARTITIONS {
        let path = partition_stats_path(output_dir, hash_name, partition_id);
        let part: PartitionStats =
            read_json(&path).with_context(|| format!("failed to read {}", path.display()))?;

        total_hashes += part.count;
        if part.count > 0 {
            min_hash = min_hash.min(part.min_hash);
            max_hash = max_hash.max(part.max_hash);
        }

        observed_pairs += part.observed_collision_pairs;
        duplicate_hashes += part.duplicate_hashes;
        collision_elements += part.collision_elements;
        max_multiplicity = max_multiplicity.max(part.max_multiplicity);
        ks_statistic = ks_statistic.max(part.ks_statistic_max);

        for (i, v) in part.frac_pass_counts.iter().enumerate() {
            frac_pass_counts[i] += *v;
        }
        for (i, v) in part.chi_bucket_counts.iter().enumerate() {
            chi_bucket_counts[i] += *v;
        }
        for (i, v) in part.byte_hist_counts.iter().enumerate() {
            byte_hist_counts[i] += *v;
        }
    }

    if total_hashes == 0 {
        min_hash = 0;
    }

    let mut frac_observed = Vec::with_capacity(FRAC_SCALES.len());
    let mut frac_expected = Vec::with_capacity(FRAC_SCALES.len());
    let mut frac_abs_error = Vec::with_capacity(FRAC_SCALES.len());
    let mut frac_max_abs_error = 0.0f64;

    for (i, scale) in FRAC_SCALES.iter().enumerate() {
        let observed = if total_hashes == 0 {
            0.0
        } else {
            frac_pass_counts[i] as f64 / total_hashes as f64
        };
        let expected = 1.0 / *scale as f64;
        let err = (observed - expected).abs();
        frac_observed.push(observed);
        frac_expected.push(expected);
        frac_abs_error.push(err);
        frac_max_abs_error = frac_max_abs_error.max(err);
    }

    let mut bit_max_abs_error = 0.0f64;
    let mut bit_worst_index = 0usize;
    if total_hashes > 0 {
        for bit in 0..64 {
            let byte_pos = bit / 8;
            let mask = 1u8 << (bit % 8);
            let mut ones = 0u64;

            for byte in 0..256usize {
                if (byte as u8 & mask) != 0 {
                    ones += byte_hist_counts[byte_pos * 256 + byte];
                }
            }

            let p1 = ones as f64 / total_hashes as f64;
            let err = (p1 - 0.5).abs();
            if err > bit_max_abs_error {
                bit_max_abs_error = err;
                bit_worst_index = bit;
            }
        }
    }

    let mut chi_squared = 0.0f64;
    let mut chi_max_rel_error = 0.0f64;
    if total_hashes > 0 {
        let expected = total_hashes as f64 / CHI_BUCKETS as f64;
        for observed in &chi_bucket_counts {
            let diff = *observed as f64 - expected;
            chi_squared += (diff * diff) / expected;
            let rel = diff.abs() / expected;
            chi_max_rel_error = chi_max_rel_error.max(rel);
        }
    }

    let expected_collision_pairs = if total_hashes > 1 {
        (total_hashes as f64) * ((total_hashes - 1) as f64) / (2.0 * U64_SPACE_F64)
    } else {
        0.0
    };

    let collision_ratio = if expected_collision_pairs > 0.0 {
        observed_pairs as f64 / expected_collision_pairs
    } else {
        0.0
    };

    let metrics = HashMetrics {
        hash_name: hash_name.to_string(),
        total_hashes,
        expected_hashes: plan.total_canonical,
        min_hash,
        max_hash,
        frac_scales: FRAC_SCALES.to_vec(),
        frac_pass_counts,
        frac_observed,
        frac_expected,
        frac_abs_error,
        frac_max_abs_error,
        bit_max_abs_error,
        bit_worst_index,
        chi_squared_statistic: chi_squared,
        chi_max_rel_error,
        ks_statistic,
        observed_collision_pairs: observed_pairs,
        expected_collision_pairs,
        collision_ratio_obs_exp: collision_ratio,
        duplicate_hashes,
        collision_elements,
        max_multiplicity,
    };

    write_json_pretty(&hash_metrics_json_path(output_dir, hash_name), &metrics)?;
    write_hash_metrics_csv(&hash_metrics_csv_path(output_dir, hash_name), &metrics)?;

    eprintln!(
        "[aggregate-hash:{}] n={} collision_ratio={:.6}",
        hash_name, metrics.total_hashes, metrics.collision_ratio_obs_exp
    );

    Ok(())
}

fn cmd_aggregate_all(output_dir: &Path) -> Result<()> {
    ensure_layout(output_dir)?;

    let mut hashes = Vec::with_capacity(HASH_COUNT);
    for hash in HASH_NAMES {
        let path = hash_metrics_json_path(output_dir, hash);
        let metrics: HashMetrics =
            read_json(&path).with_context(|| format!("failed to read {}", path.display()))?;
        hashes.push(metrics);
    }

    let report = CombinedReport { hashes };
    write_json_pretty(&combined_report_json_path(output_dir), &report)?;
    write_combined_csv(&combined_report_csv_path(output_dir), &report)?;

    eprintln!(
        "[aggregate-all] wrote {} and {}",
        combined_report_json_path(output_dir).display(),
        combined_report_csv_path(output_dir).display()
    );

    Ok(())
}

fn process_subrange(
    subrange: Range<u64>,
    io_senders: Arc<Vec<SyncSender<WriteChunk>>>,
    io_threads: usize,
    local_buffer_values: usize,
) -> Result<WorkerSummary> {
    let mut out = WorkerSummary::default();
    let mut local_buffers = (0..(HASH_COUNT * PARTITIONS))
        .map(|_| Vec::<u64>::with_capacity(local_buffer_values))
        .collect::<Vec<_>>();

    for raw in subrange {
        let rc = reverse_complement_21(raw);
        if raw > rc {
            continue;
        }

        out.canonical_count += 1;

        let hashes = [
            jamhash::jamhash_u64(raw),
            murmur3_u64(raw),
            xxhash3_u64(raw),
            wang64(raw),
        ];

        for (hash_idx, hash) in hashes.into_iter().enumerate() {
            let partition = (hash >> (64 - PARTITION_BITS)) as usize;
            let buffer_index = hash_idx * PARTITIONS + partition;
            let buf = &mut local_buffers[buffer_index];

            buf.push(hash);
            out.partition_counts[hash_idx][partition] += 1;

            if buf.len() >= local_buffer_values {
                flush_partition_buffer(
                    &io_senders,
                    io_threads,
                    buffer_index,
                    buf,
                    local_buffer_values,
                )?;
            }
        }
    }

    for (idx, buf) in local_buffers.iter_mut().enumerate() {
        flush_partition_buffer(&io_senders, io_threads, idx, buf, local_buffer_values)?;
    }

    Ok(out)
}

fn compute_partition_stats(
    hash_name: &str,
    partition_id: usize,
    values: &[u64],
    partition_offset: u64,
    total_hashes: u64,
) -> Result<PartitionStats> {
    if values.is_empty() {
        return Ok(PartitionStats {
            hash_name: hash_name.to_string(),
            partition_id,
            count: 0,
            min_hash: 0,
            max_hash: 0,
            observed_collision_pairs: 0,
            duplicate_hashes: 0,
            collision_elements: 0,
            max_multiplicity: 0,
            ks_statistic_max: 0.0,
            frac_pass_counts: vec![0u64; FRAC_SCALES.len()],
            chi_bucket_counts: vec![0u64; CHI_BUCKETS],
            byte_hist_counts: vec![0u64; 8 * 256],
        });
    }

    let mut frac_pass_counts = vec![0u64; FRAC_SCALES.len()];
    let frac_thresholds = FRAC_SCALES.map(|s| u64::MAX / s);

    let mut chi_bucket_counts = vec![0u64; CHI_BUCKETS];
    let mut byte_hist_counts = vec![0u64; 8 * 256];

    let n_f = total_hashes as f64;
    let mut ks_statistic_max = 0.0f64;

    let mut observed_collision_pairs = 0u64;
    let mut duplicate_hashes = 0u64;
    let mut collision_elements = 0u64;
    let mut max_multiplicity = 1u64;

    let mut prev = values[0];
    let mut run_len = 0u64;

    for (local_idx, &hash) in values.iter().enumerate() {
        let global_rank = partition_offset + local_idx as u64 + 1;
        let x = hash as f64 / U64_SPACE_F64;
        let rank_f = global_rank as f64;
        let d1 = (rank_f / n_f - x).abs();
        let d2 = (x - ((rank_f - 1.0) / n_f)).abs();
        ks_statistic_max = ks_statistic_max.max(d1.max(d2));

        for (i, threshold) in frac_thresholds.iter().enumerate() {
            if hash < *threshold {
                frac_pass_counts[i] += 1;
            }
        }

        let chi_bucket = (hash >> (64 - CHI_BUCKET_BITS)) as usize;
        chi_bucket_counts[chi_bucket] += 1;

        for byte_pos in 0..8 {
            let byte = ((hash >> (8 * byte_pos)) & 0xFF) as usize;
            byte_hist_counts[byte_pos * 256 + byte] += 1;
        }

        if hash == prev {
            run_len += 1;
        } else {
            finalize_collision_run(
                run_len,
                &mut observed_collision_pairs,
                &mut duplicate_hashes,
                &mut collision_elements,
                &mut max_multiplicity,
            );
            prev = hash;
            run_len = 1;
        }
    }

    finalize_collision_run(
        run_len,
        &mut observed_collision_pairs,
        &mut duplicate_hashes,
        &mut collision_elements,
        &mut max_multiplicity,
    );

    Ok(PartitionStats {
        hash_name: hash_name.to_string(),
        partition_id,
        count: values.len() as u64,
        min_hash: values[0],
        max_hash: values[values.len() - 1],
        observed_collision_pairs,
        duplicate_hashes,
        collision_elements,
        max_multiplicity,
        ks_statistic_max,
        frac_pass_counts,
        chi_bucket_counts,
        byte_hist_counts,
    })
}

fn finalize_collision_run(
    run_len: u64,
    observed_collision_pairs: &mut u64,
    duplicate_hashes: &mut u64,
    collision_elements: &mut u64,
    max_multiplicity: &mut u64,
) {
    if run_len > 1 {
        let pairs = run_len.saturating_mul(run_len.saturating_sub(1)) / 2;
        *observed_collision_pairs = observed_collision_pairs.saturating_add(pairs);
        *duplicate_hashes = duplicate_hashes.saturating_add(1);
        *collision_elements = collision_elements.saturating_add(run_len);
        *max_multiplicity = (*max_multiplicity).max(run_len);
    }
}

fn flush_partition_buffer(
    io_senders: &[SyncSender<WriteChunk>],
    io_threads: usize,
    stream_idx: usize,
    buffer: &mut Vec<u64>,
    local_buffer_values: usize,
) -> Result<()> {
    if buffer.is_empty() {
        return Ok(());
    }

    let chunk_values = std::mem::replace(buffer, Vec::with_capacity(local_buffer_values));
    let worker_idx = stream_idx % io_threads;

    io_senders[worker_idx]
        .send(WriteChunk {
            stream_idx,
            values: chunk_values,
        })
        .map_err(|_| anyhow::anyhow!("map io channel send failed for stream {}", stream_idx))?;

    Ok(())
}

fn start_map_io_workers(
    base: &Path,
    map_task: usize,
    io_threads: usize,
    channel_capacity: usize,
    map_file_buffer_bytes: usize,
) -> Result<(
    Vec<SyncSender<WriteChunk>>,
    Vec<thread::JoinHandle<Result<()>>>,
)> {
    let mut senders = Vec::with_capacity(io_threads);
    let mut handles = Vec::with_capacity(io_threads);

    for worker_idx in 0..io_threads {
        let (tx, rx) = sync_channel::<WriteChunk>(channel_capacity);
        let base = base.to_path_buf();

        let handle = thread::Builder::new()
            .name(format!("map-io-{worker_idx}"))
            .spawn(move || {
                map_io_worker_loop(
                    base,
                    map_task,
                    worker_idx,
                    io_threads,
                    map_file_buffer_bytes,
                    rx,
                )
            })
            .with_context(|| format!("failed to spawn map io worker {}", worker_idx))?;

        senders.push(tx);
        handles.push(handle);
    }

    Ok((senders, handles))
}

fn map_io_worker_loop(
    base: PathBuf,
    map_task: usize,
    worker_idx: usize,
    worker_count: usize,
    map_file_buffer_bytes: usize,
    rx: Receiver<WriteChunk>,
) -> Result<()> {
    let mut writers = (0..(HASH_COUNT * PARTITIONS))
        .map(|_| None)
        .collect::<Vec<Option<BufWriter<File>>>>();

    for stream_idx in 0..(HASH_COUNT * PARTITIONS) {
        if stream_idx % worker_count != worker_idx {
            continue;
        }

        let (hash_name, partition) = stream_to_hash_partition(stream_idx);
        let path = shard_path(&base, hash_name, map_task, partition);
        if let Some(parent) = path.parent() {
            create_dir_all(parent)
                .with_context(|| format!("failed to create {}", parent.display()))?;
        }
        let file =
            File::create(&path).with_context(|| format!("failed to create {}", path.display()))?;
        let writer = BufWriter::with_capacity(map_file_buffer_bytes, file);
        writers[stream_idx] = Some(writer);
    }

    while let Ok(chunk) = rx.recv() {
        let writer = writers[chunk.stream_idx]
            .as_mut()
            .ok_or_else(|| anyhow::anyhow!("missing writer for stream {}", chunk.stream_idx))?;
        writer.write_all(cast_slice::<u64, u8>(&chunk.values))?;
    }

    for writer in writers.iter_mut().flatten() {
        writer.flush()?;
    }

    Ok(())
}

fn stream_to_hash_partition(stream_idx: usize) -> (&'static str, usize) {
    let hash_idx = stream_idx / PARTITIONS;
    let partition = stream_idx % PARTITIONS;
    (HASH_NAMES[hash_idx], partition)
}

fn copy_map_shards_to_output(local_base: &Path, output_base: &Path, map_task: usize) -> Result<()> {
    for hash_name in HASH_NAMES {
        for partition in 0..PARTITIONS {
            let src = shard_path(local_base, hash_name, map_task, partition);
            let dst = shard_path(output_base, hash_name, map_task, partition);

            if let Some(parent) = dst.parent() {
                create_dir_all(parent)
                    .with_context(|| format!("failed to create {}", parent.display()))?;
            }

            let status = ProcessCommand::new("cpb")
                .arg(&src)
                .arg(&dst)
                .status()
                .with_context(|| {
                    format!(
                        "failed to execute cpb from {} to {}",
                        src.display(),
                        dst.display()
                    )
                })?;

            if !status.success() {
                bail!(
                    "cpb failed with status {} for {} -> {}",
                    status,
                    src.display(),
                    dst.display()
                );
            }
        }
    }

    Ok(())
}

fn read_u64_file_into_vec(path: &Path, out: &mut Vec<u64>) -> Result<()> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut reader = BufReader::with_capacity(IO_BUFFER_BYTES, file);
    let mut buf = vec![0u8; IO_BUFFER_BYTES];
    let mut carry_len = 0usize;

    loop {
        let n = reader.read(&mut buf[carry_len..])?;
        if n == 0 {
            break;
        }

        let filled = carry_len + n;
        let usable = filled / 8 * 8;

        let mut idx = 0usize;
        while idx < usable {
            let ptr = unsafe { buf.as_ptr().add(idx) as *const u64 };
            let value = unsafe { ptr.read_unaligned() };
            out.push(u64::from_le(value));
            idx += 8;
        }

        carry_len = filled - usable;
        if carry_len > 0 {
            buf.copy_within(usable..filled, 0);
        }
    }

    if carry_len != 0 {
        bail!(
            "{} has trailing bytes (not a multiple of 8)",
            path.display()
        );
    }

    Ok(())
}

fn write_u64_file(path: &Path, values: &[u64]) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = BufWriter::with_capacity(IO_BUFFER_BYTES, file);
    let bytes = cast_slice::<u64, u8>(values);

    let mut offset = 0usize;
    while offset < bytes.len() {
        let end = (offset + IO_BUFFER_BYTES).min(bytes.len());
        writer.write_all(&bytes[offset..end])?;
        offset = end;
    }
    writer.flush()?;

    Ok(())
}

fn write_hash_metrics_csv(path: &Path, metrics: &HashMetrics) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = BufWriter::with_capacity(IO_BUFFER_BYTES, file);

    writeln!(writer, "metric,value")?;
    writeln!(writer, "hash_name,{}", metrics.hash_name)?;
    writeln!(writer, "total_hashes,{}", metrics.total_hashes)?;
    writeln!(writer, "expected_hashes,{}", metrics.expected_hashes)?;
    writeln!(writer, "min_hash,{}", metrics.min_hash)?;
    writeln!(writer, "max_hash,{}", metrics.max_hash)?;
    writeln!(
        writer,
        "frac_max_abs_error,{:.12e}",
        metrics.frac_max_abs_error
    )?;
    writeln!(
        writer,
        "bit_max_abs_error,{:.12e}",
        metrics.bit_max_abs_error
    )?;
    writeln!(writer, "bit_worst_index,{}", metrics.bit_worst_index)?;
    writeln!(
        writer,
        "chi_squared_statistic,{:.12e}",
        metrics.chi_squared_statistic
    )?;
    writeln!(
        writer,
        "chi_max_rel_error,{:.12e}",
        metrics.chi_max_rel_error
    )?;
    writeln!(writer, "ks_statistic,{:.12e}", metrics.ks_statistic)?;
    writeln!(
        writer,
        "observed_collision_pairs,{}",
        metrics.observed_collision_pairs
    )?;
    writeln!(
        writer,
        "expected_collision_pairs,{:.12e}",
        metrics.expected_collision_pairs
    )?;
    writeln!(
        writer,
        "collision_ratio_obs_exp,{:.12e}",
        metrics.collision_ratio_obs_exp
    )?;
    writeln!(writer, "duplicate_hashes,{}", metrics.duplicate_hashes)?;
    writeln!(writer, "collision_elements,{}", metrics.collision_elements)?;
    writeln!(writer, "max_multiplicity,{}", metrics.max_multiplicity)?;

    for (i, scale) in metrics.frac_scales.iter().enumerate() {
        writeln!(
            writer,
            "frac_scale_{}_pass_count,{}",
            scale, metrics.frac_pass_counts[i]
        )?;
        writeln!(
            writer,
            "frac_scale_{}_observed,{:.12e}",
            scale, metrics.frac_observed[i]
        )?;
        writeln!(
            writer,
            "frac_scale_{}_expected,{:.12e}",
            scale, metrics.frac_expected[i]
        )?;
        writeln!(
            writer,
            "frac_scale_{}_abs_error,{:.12e}",
            scale, metrics.frac_abs_error[i]
        )?;
    }

    writer.flush()?;
    Ok(())
}

fn write_combined_csv(path: &Path, report: &CombinedReport) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = BufWriter::with_capacity(IO_BUFFER_BYTES, file);

    writeln!(
        writer,
        "hash,total_hashes,frac_max_abs_error,bit_max_abs_error,chi_squared_statistic,ks_statistic,observed_collision_pairs,expected_collision_pairs,collision_ratio_obs_exp"
    )?;

    for m in &report.hashes {
        writeln!(
            writer,
            "{},{},{:.12e},{:.12e},{:.12e},{:.12e},{},{:.12e},{:.12e}",
            m.hash_name,
            m.total_hashes,
            m.frac_max_abs_error,
            m.bit_max_abs_error,
            m.chi_squared_statistic,
            m.ks_statistic,
            m.observed_collision_pairs,
            m.expected_collision_pairs,
            m.collision_ratio_obs_exp
        )?;
    }

    writer.flush()?;
    Ok(())
}

fn validate_plan(plan: &PartitionPlan) -> Result<()> {
    if plan.map_tasks != MAP_TASKS {
        bail!(
            "invalid plan map_tasks={} expected={}",
            plan.map_tasks,
            MAP_TASKS
        );
    }
    if plan.partitions != PARTITIONS {
        bail!(
            "invalid plan partitions={} expected={}",
            plan.partitions,
            PARTITIONS
        );
    }
    if plan.hash_partition_counts.len() != HASH_COUNT
        || plan.hash_partition_offsets.len() != HASH_COUNT
    {
        bail!("invalid plan hash dimensions");
    }
    for h in 0..HASH_COUNT {
        if plan.hash_partition_counts[h].len() != PARTITIONS
            || plan.hash_partition_offsets[h].len() != PARTITIONS
        {
            bail!(
                "invalid plan partition dimensions for hash {}",
                HASH_NAMES[h]
            );
        }
    }
    Ok(())
}

fn configured_threads() -> usize {
    std::env::var("K21_THREADS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&n| n > 0)
        .or_else(|| {
            std::env::var("SLURM_CPUS_PER_TASK")
                .ok()
                .and_then(|v| v.parse::<usize>().ok())
                .filter(|&n| n > 0)
        })
        .unwrap_or_else(|| {
            std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1)
        })
}

fn configured_map_io_threads(total_threads: usize) -> usize {
    let requested = std::env::var("K21_IO_THREADS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&n| n > 0)
        .unwrap_or_else(|| {
            let scaled = (total_threads / 8).max(4);
            scaled.min(DEFAULT_MAP_IO_THREADS)
        });

    requested.clamp(1, total_threads.saturating_sub(1).max(1))
}

fn configured_local_buffer_values() -> usize {
    std::env::var("K21_LOCAL_BUFFER_VALUES")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&n| n > 0)
        .unwrap_or(DEFAULT_LOCAL_BUFFER_VALUES)
}

fn configured_map_channel_entries(local_buffer_values: usize) -> usize {
    if let Some(entries) = std::env::var("K21_IO_CHANNEL_ENTRIES")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&n| n > 0)
    {
        return entries;
    }

    if let Some(chunks) = std::env::var("K21_IO_CHANNEL_CAP")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&n| n > 0)
    {
        return chunks.saturating_mul(local_buffer_values);
    }

    DEFAULT_MAP_CHANNEL_ENTRIES
}

fn configured_map_file_buffer_bytes() -> usize {
    if let Some(bytes) = std::env::var("K21_MAP_FILE_BUFFER_BYTES")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&n| n > 0)
    {
        return bytes;
    }

    let mb = std::env::var("K21_MAP_FILE_BUFFER_MB")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&n| n > 0)
        .unwrap_or(DEFAULT_MAP_FILE_BUFFER_BYTES / (1024 * 1024));

    mb.saturating_mul(1024 * 1024)
}

fn configured_local_scratch_base() -> PathBuf {
    std::env::var("K21_LOCAL_SCRATCH_BASE")
        .ok()
        .filter(|v| !v.trim().is_empty())
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(DEFAULT_LOCAL_SCRATCH_BASE))
}

fn map_task_scratch_dir(base: &Path, task_id: usize) -> PathBuf {
    let job_id = std::env::var("SLURM_JOB_ID").unwrap_or_else(|_| "manual".to_string());
    base.join(format!("k21-map-{job_id}-task-{task_id:02}"))
}

fn configure_rayon(threads: usize) {
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global();
}

fn div_ceil(a: usize, b: usize) -> usize {
    if b == 0 {
        return 0;
    }
    (a + (b - 1)) / b
}

fn slurm_array_task_id_or_default(default: usize) -> Result<usize> {
    match std::env::var("SLURM_ARRAY_TASK_ID") {
        Ok(v) => v
            .parse::<usize>()
            .with_context(|| format!("invalid SLURM_ARRAY_TASK_ID: {}", v)),
        Err(_) => Ok(default),
    }
}

fn ensure_layout(base: &Path) -> Result<()> {
    create_dir_all(base).with_context(|| format!("failed to create {}", base.display()))?;
    create_dir_all(meta_map_dir(base))?;

    for hash in HASH_NAMES {
        create_dir_all(hash_shards_dir(base, hash))?;
        create_dir_all(hash_sorted_dir(base, hash))?;
        create_dir_all(hash_partition_stats_dir(base, hash))?;
    }

    Ok(())
}

fn get_chunk_range(total: u64, chunks: u64, idx: u64) -> Range<u64> {
    assert!(idx < chunks, "chunk index out of bounds");
    let base = total / chunks;
    let rem = total % chunks;

    let start = if idx < rem {
        idx * (base + 1)
    } else {
        rem * (base + 1) + (idx - rem) * base
    };
    let size = if idx < rem { base + 1 } else { base };
    start..(start + size)
}

fn split_range(range: Range<u64>, parts: usize) -> Vec<Range<u64>> {
    let total = range.end.saturating_sub(range.start);
    if total == 0 {
        return Vec::new();
    }

    let parts = parts.max(1).min(total as usize);
    let mut out = Vec::with_capacity(parts);

    let base = total / parts as u64;
    let rem = total % parts as u64;

    let mut start = range.start;
    for i in 0..parts {
        let size = if (i as u64) < rem { base + 1 } else { base };
        let end = start + size;
        out.push(start..end);
        start = end;
    }

    out
}

#[inline]
fn murmur3_u64(kmer: u64) -> u64 {
    fastmurmur3::murmur3_x64_128(&kmer.to_be_bytes(), 42) as u64
}

#[inline]
fn xxhash3_u64(kmer: u64) -> u64 {
    xxhash_rust::xxh3::xxh3_64(&kmer.to_be_bytes())
}

#[inline]
fn wang64(mut key: u64) -> u64 {
    key = (!key).wrapping_add(key << 21);
    key ^= key >> 24;
    key = key.wrapping_add(key << 3).wrapping_add(key << 8);
    key ^= key >> 14;
    key = key.wrapping_add(key << 2).wrapping_add(key << 4);
    key ^= key >> 28;
    key.wrapping_add(key << 31)
}

#[inline]
fn reverse_complement_21(sequence: u64) -> u64 {
    let mut new_kmer = sequence;
    new_kmer =
        ((new_kmer >> 2) & 0x3333_3333_3333_3333) | ((new_kmer & 0x3333_3333_3333_3333) << 2);
    new_kmer =
        ((new_kmer >> 4) & 0x0F0F_0F0F_0F0F_0F0F) | ((new_kmer & 0x0F0F_0F0F_0F0F_0F0F) << 4);
    new_kmer =
        ((new_kmer >> 8) & 0x00FF_00FF_00FF_00FF) | ((new_kmer & 0x00FF_00FF_00FF_00FF) << 8);
    new_kmer =
        ((new_kmer >> 16) & 0x0000_FFFF_0000_FFFF) | ((new_kmer & 0x0000_FFFF_0000_FFFF) << 16);
    new_kmer =
        ((new_kmer >> 32) & 0x0000_0000_FFFF_FFFF) | ((new_kmer & 0x0000_0000_FFFF_FFFF) << 32);
    new_kmer ^= 0xFFFF_FFFF_FFFF_FFFF;
    new_kmer >> 22
}

fn map_summary_path(base: &Path, task: usize) -> PathBuf {
    meta_map_dir(base).join(format!("map_task_{task:02}.json"))
}

fn partition_plan_path(base: &Path) -> PathBuf {
    base.join("_meta").join("partition_plan.json")
}

fn meta_map_dir(base: &Path) -> PathBuf {
    base.join("_meta").join("map")
}

fn hash_shards_dir(base: &Path, hash: &str) -> PathBuf {
    base.join(hash).join("shards")
}

fn hash_sorted_dir(base: &Path, hash: &str) -> PathBuf {
    base.join(hash).join("sorted")
}

fn hash_partition_stats_dir(base: &Path, hash: &str) -> PathBuf {
    base.join(hash).join("partition_stats")
}

fn shard_path(base: &Path, hash: &str, map_task: usize, partition: usize) -> PathBuf {
    hash_shards_dir(base, hash).join(format!("map_{map_task:02}/part_{partition:02}.bin"))
}

fn sorted_partition_path(base: &Path, hash: &str, partition: usize) -> PathBuf {
    hash_sorted_dir(base, hash).join(format!("part_{partition:02}.bin"))
}

fn partition_stats_path(base: &Path, hash: &str, partition: usize) -> PathBuf {
    hash_partition_stats_dir(base, hash).join(format!("part_{partition:02}.json"))
}

fn hash_metrics_json_path(base: &Path, hash: &str) -> PathBuf {
    base.join(hash).join("metrics.json")
}

fn hash_metrics_csv_path(base: &Path, hash: &str) -> PathBuf {
    base.join(hash).join("metrics.csv")
}

fn combined_report_json_path(base: &Path) -> PathBuf {
    base.join("combined_metrics.json")
}

fn combined_report_csv_path(base: &Path) -> PathBuf {
    base.join("combined_metrics.csv")
}

fn write_json_pretty<T: Serialize>(path: &Path, value: &T) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = BufWriter::with_capacity(IO_BUFFER_BYTES, file);
    serde_json::to_writer_pretty(&mut writer, value)
        .with_context(|| format!("failed to write {}", path.display()))?;
    writer.flush()?;
    Ok(())
}

fn read_json<T: for<'de> Deserialize<'de>>(path: &Path) -> Result<T> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let reader = BufReader::with_capacity(IO_BUFFER_BYTES, file);
    let value = serde_json::from_reader(reader)
        .with_context(|| format!("failed to parse {}", path.display()))?;
    Ok(value)
}
