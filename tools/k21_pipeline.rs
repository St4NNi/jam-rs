use anyhow::{bail, Context, Result};
use bytemuck::cast_slice;
use clap::{Parser, Subcommand};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::fs::{create_dir_all, remove_dir_all, File, OpenOptions};
use std::io::{BufReader, BufWriter, Read, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::process::Command as ProcessCommand;
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};
use std::sync::Arc;
use std::thread;
use std::time::{Duration, Instant};

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
const DEFAULT_MAP_MAX_OPEN_WRITERS_PER_WORKER: usize = 8;
const DEFAULT_LOCAL_SCRATCH_BASE: &str = "/var/scratch/sbeyvers/k21";
const DEFAULT_SORT_CHUNK_GIB: usize = 4;
const DEFAULT_SORT_SCRATCH_BASE: &str = "/var/scratch/sbeyvers/k21_sort";
const SORT_OUTPUT_BUFFER_VALUES: usize = 1 << 20;

const TRANSIENT_IO_RETRIES: usize = 8;
const TRANSIENT_IO_SLEEP_MS: u64 = 250;
const CPB_STATUS_RETRIES: usize = 6;

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
    /// Linear uniformity map: process one shard of the canonical k=21 space (array task)
    LinearMap {
        /// Output directory root
        output_dir: PathBuf,
    },
    /// Linear uniformity aggregate: combine all linear-map shards into final report
    LinearAggregate {
        /// Output directory root
        output_dir: PathBuf,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct MapTaskSummary {
    task_id: usize,
    raw_start: u64,
    raw_end: u64,
    #[serde(default)]
    raw_ranges: Vec<[u64; 2]>,
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

const LINEAR_TASKS: usize = 32;

#[derive(Debug, Clone, Serialize, Deserialize)]
struct LinearHashAccum {
    hash_name: String,
    sum_x: f64,
    sum_y: f64,
    sum_xy: f64,
    sum_x2: f64,
    sum_y2: f64,
    serial_n: u64,
    order_preserved: u64,
    gap_bins: Vec<u64>,
    hamming_sum: u64,
    hamming_sum_sq: u64,
    hamming_min: u32,
    hamming_max: u32,
    mutation_count: u64,
    mutation_diff_bins: Vec<u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct LinearTaskResult {
    task_id: usize,
    canonical_count: u64,
    elapsed_seconds: f64,
    hash_accums: Vec<LinearHashAccum>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct LinearHashMetrics {
    hash_name: String,
    serial_correlation: f64,
    order_preservation_frac: f64,
    gap_chi_squared: f64,
    mutation_avg_hamming: f64,
    mutation_hamming_std: f64,
    mutation_min_hamming: u32,
    mutation_max_hamming: u32,
    mutation_diff_chi_squared: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct LinearReport {
    total_canonical: u64,
    tasks: usize,
    hashes: Vec<LinearHashMetrics>,
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

struct SortStatsAccumulator {
    hash_name: String,
    partition_id: usize,
    partition_offset: u64,
    total_hashes: u64,
    count: u64,
    min_hash: u64,
    max_hash: u64,
    observed_collision_pairs: u64,
    duplicate_hashes: u64,
    collision_elements: u64,
    max_multiplicity: u64,
    ks_statistic_max: f64,
    frac_pass_counts: Vec<u64>,
    frac_thresholds: [u64; FRAC_SCALES.len()],
    chi_bucket_counts: Vec<u64>,
    byte_hist_counts: Vec<u64>,
    prev_hash: Option<u64>,
    run_len: u64,
}

impl SortStatsAccumulator {
    fn new(hash_name: &str, partition_id: usize, partition_offset: u64, total_hashes: u64) -> Self {
        Self {
            hash_name: hash_name.to_string(),
            partition_id,
            partition_offset,
            total_hashes,
            count: 0,
            min_hash: 0,
            max_hash: 0,
            observed_collision_pairs: 0,
            duplicate_hashes: 0,
            collision_elements: 0,
            max_multiplicity: 0,
            ks_statistic_max: 0.0,
            frac_pass_counts: vec![0u64; FRAC_SCALES.len()],
            frac_thresholds: FRAC_SCALES.map(|s| u64::MAX / s),
            chi_bucket_counts: vec![0u64; CHI_BUCKETS],
            byte_hist_counts: vec![0u64; 8 * 256],
            prev_hash: None,
            run_len: 0,
        }
    }

    fn observe(&mut self, hash: u64) {
        self.count += 1;
        if self.count == 1 {
            self.min_hash = hash;
        }
        self.max_hash = hash;

        if self.total_hashes > 0 {
            let n_f = self.total_hashes as f64;
            let global_rank = self.partition_offset + self.count;
            let rank_f = global_rank as f64;
            let x = hash as f64 / U64_SPACE_F64;
            let d1 = (rank_f / n_f - x).abs();
            let d2 = (x - ((rank_f - 1.0) / n_f)).abs();
            self.ks_statistic_max = self.ks_statistic_max.max(d1.max(d2));
        }

        for (i, threshold) in self.frac_thresholds.iter().enumerate() {
            if hash < *threshold {
                self.frac_pass_counts[i] += 1;
            }
        }

        let chi_bucket = (hash >> (64 - CHI_BUCKET_BITS)) as usize;
        self.chi_bucket_counts[chi_bucket] += 1;

        for byte_pos in 0..8 {
            let byte = ((hash >> (8 * byte_pos)) & 0xFF) as usize;
            self.byte_hist_counts[byte_pos * 256 + byte] += 1;
        }

        match self.prev_hash {
            Some(prev) if prev == hash => {
                self.run_len += 1;
            }
            Some(_) => {
                finalize_collision_run(
                    self.run_len,
                    &mut self.observed_collision_pairs,
                    &mut self.duplicate_hashes,
                    &mut self.collision_elements,
                    &mut self.max_multiplicity,
                );
                self.prev_hash = Some(hash);
                self.run_len = 1;
            }
            None => {
                self.prev_hash = Some(hash);
                self.run_len = 1;
            }
        }
    }

    fn finalize(mut self) -> PartitionStats {
        finalize_collision_run(
            self.run_len,
            &mut self.observed_collision_pairs,
            &mut self.duplicate_hashes,
            &mut self.collision_elements,
            &mut self.max_multiplicity,
        );

        if self.count == 0 {
            self.max_multiplicity = 0;
        }

        PartitionStats {
            hash_name: self.hash_name,
            partition_id: self.partition_id,
            count: self.count,
            min_hash: self.min_hash,
            max_hash: self.max_hash,
            observed_collision_pairs: self.observed_collision_pairs,
            duplicate_hashes: self.duplicate_hashes,
            collision_elements: self.collision_elements,
            max_multiplicity: self.max_multiplicity,
            ks_statistic_max: self.ks_statistic_max,
            frac_pass_counts: self.frac_pass_counts,
            chi_bucket_counts: self.chi_bucket_counts,
            byte_hist_counts: self.byte_hist_counts,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct MergeHeapItem {
    value: u64,
    run_idx: usize,
}

impl Ord for MergeHeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .value
            .cmp(&self.value)
            .then_with(|| other.run_idx.cmp(&self.run_idx))
    }
}

impl PartialOrd for MergeHeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

struct RunStream {
    reader: BufReader<File>,
    buf: Vec<u8>,
    pos: usize,
    usable_len: usize,
    carry_len: usize,
    eof: bool,
}

impl RunStream {
    fn open(path: &Path) -> Result<Self> {
        let file =
            File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
        Ok(Self {
            reader: BufReader::with_capacity(IO_BUFFER_BYTES, file),
            buf: vec![0u8; IO_BUFFER_BYTES],
            pos: 0,
            usable_len: 0,
            carry_len: 0,
            eof: false,
        })
    }

    fn next_u64(&mut self) -> Result<Option<u64>> {
        loop {
            if self.pos < self.usable_len {
                let ptr = unsafe { self.buf.as_ptr().add(self.pos) as *const u64 };
                let value = unsafe { ptr.read_unaligned() };
                self.pos += 8;
                return Ok(Some(u64::from_le(value)));
            }

            if self.eof {
                if self.carry_len != 0 {
                    bail!("run stream has trailing bytes (not a multiple of 8)");
                }
                return Ok(None);
            }

            let n = self.reader.read(&mut self.buf[self.carry_len..])?;
            if n == 0 {
                self.eof = true;
                continue;
            }

            let filled = self.carry_len + n;
            self.usable_len = filled / 8 * 8;
            self.pos = 0;
            self.carry_len = filled - self.usable_len;
            if self.carry_len > 0 {
                self.buf.copy_within(self.usable_len..filled, 0);
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
        Command::LinearMap { output_dir } => cmd_linear_map(&output_dir),
        Command::LinearAggregate { output_dir } => cmd_linear_aggregate(&output_dir),
    }
}

fn cmd_linear_map(output_dir: &Path) -> Result<()> {
    let linear_dir = output_dir.join("linear_test");
    create_dir_all(&linear_dir)?;

    let task_id = slurm_array_task_id_or_default(0)?;
    if task_id >= LINEAR_TASKS {
        bail!(
            "linear-map task id {} out of range 0..{}",
            task_id,
            LINEAR_TASKS - 1
        );
    }

    let start = Instant::now();
    let raw_per_task = RAW_SPACE_SIZE / LINEAR_TASKS as u64;
    let raw_start = task_id as u64 * raw_per_task;
    let raw_end = if task_id == LINEAR_TASKS - 1 {
        RAW_SPACE_SIZE
    } else {
        (task_id as u64 + 1) * raw_per_task
    };

    let hash_fns: [(&str, fn(u64) -> u64); HASH_COUNT] = [
        ("jamhash", |x| jamhash::jamhash_u64(x)),
        ("murmur3", murmur3_u64),
        ("xxhash3", xxhash3_u64),
        ("wang64", wang64),
    ];

    // Split this task's range across rayon threads
    let num_threads = rayon::current_num_threads().max(1) as u64;
    let raw_per_thread = (raw_end - raw_start) / num_threads;

    struct ThreadAccum {
        prev_hash: [Option<u64>; HASH_COUNT],
        sum_x: [f64; HASH_COUNT],
        sum_y: [f64; HASH_COUNT],
        sum_xy: [f64; HASH_COUNT],
        sum_x2: [f64; HASH_COUNT],
        sum_y2: [f64; HASH_COUNT],
        serial_n: [u64; HASH_COUNT],
        order_preserved: [u64; HASH_COUNT],
        gap_bins: [[u64; 256]; HASH_COUNT],
        hamming_sum: [u64; HASH_COUNT],
        hamming_sum_sq: [u64; HASH_COUNT],
        hamming_min: [u32; HASH_COUNT],
        hamming_max: [u32; HASH_COUNT],
        mutation_count: [u64; HASH_COUNT],
        mutation_diff_bins: [[u64; 256]; HASH_COUNT],
        canonical_count: u64,
    }

    impl ThreadAccum {
        fn new() -> Self {
            Self {
                prev_hash: [None; HASH_COUNT],
                sum_x: [0.0; HASH_COUNT],
                sum_y: [0.0; HASH_COUNT],
                sum_xy: [0.0; HASH_COUNT],
                sum_x2: [0.0; HASH_COUNT],
                sum_y2: [0.0; HASH_COUNT],
                serial_n: [0; HASH_COUNT],
                order_preserved: [0; HASH_COUNT],
                gap_bins: [[0; 256]; HASH_COUNT],
                hamming_sum: [0; HASH_COUNT],
                hamming_sum_sq: [0; HASH_COUNT],
                hamming_min: [64; HASH_COUNT],
                hamming_max: [0; HASH_COUNT],
                mutation_count: [0; HASH_COUNT],
                mutation_diff_bins: [[0; 256]; HASH_COUNT],
                canonical_count: 0,
            }
        }
    }

    let thread_results: Vec<ThreadAccum> = (0..num_threads)
        .into_par_iter()
        .map(|t| {
            let t_start = raw_start + t * raw_per_thread;
            let t_end = if t == num_threads - 1 {
                raw_end
            } else {
                raw_start + (t + 1) * raw_per_thread
            };

            let mut acc = ThreadAccum::new();

            for raw in t_start..t_end {
                let rc = reverse_complement_21(raw);
                if raw > rc {
                    continue;
                }

                for (hi, (_name, hash_fn)) in hash_fns.iter().enumerate() {
                    let h = hash_fn(raw);

                    if let Some(prev) = acc.prev_hash[hi] {
                        let x = prev as f64;
                        let y = h as f64;
                        acc.sum_x[hi] += x;
                        acc.sum_y[hi] += y;
                        acc.sum_xy[hi] += x * y;
                        acc.sum_x2[hi] += x * x;
                        acc.sum_y2[hi] += y * y;
                        acc.serial_n[hi] += 1;

                        if h > prev {
                            acc.order_preserved[hi] += 1;
                        }

                        let gap = h.wrapping_sub(prev);
                        acc.gap_bins[hi][(gap >> 56) as usize] += 1;
                    }
                    acc.prev_hash[hi] = Some(h);

                    if raw & 0x3F == 0 {
                        for pos in 0..KMER_SIZE {
                            let shift = 2 * pos as u32;
                            let current_base = (raw >> shift) & 0b11;
                            for base in 0u64..4 {
                                if base == current_base {
                                    continue;
                                }
                                let mutated = (raw & !(0b11u64 << shift)) | (base << shift);
                                let mh = hash_fn(mutated);
                                let hamming = (h ^ mh).count_ones();
                                acc.hamming_sum[hi] += hamming as u64;
                                acc.hamming_sum_sq[hi] += (hamming as u64) * (hamming as u64);
                                acc.hamming_min[hi] = acc.hamming_min[hi].min(hamming);
                                acc.hamming_max[hi] = acc.hamming_max[hi].max(hamming);
                                acc.mutation_count[hi] += 1;

                                let diff = h.wrapping_sub(mh);
                                acc.mutation_diff_bins[hi][(diff >> 56) as usize] += 1;
                            }
                        }
                    }
                }
                acc.canonical_count += 1;
            }

            acc
        })
        .collect();

    // Merge thread results into per-hash accumulators
    let mut canonical_count = 0u64;
    let mut hash_accums = Vec::with_capacity(HASH_COUNT);
    for (hi, (name, _)) in hash_fns.iter().enumerate() {
        let mut a = LinearHashAccum {
            hash_name: name.to_string(),
            sum_x: 0.0,
            sum_y: 0.0,
            sum_xy: 0.0,
            sum_x2: 0.0,
            sum_y2: 0.0,
            serial_n: 0,
            order_preserved: 0,
            gap_bins: vec![0u64; 256],
            hamming_sum: 0,
            hamming_sum_sq: 0,
            hamming_min: 64,
            hamming_max: 0,
            mutation_count: 0,
            mutation_diff_bins: vec![0u64; 256],
        };
        for t in &thread_results {
            if hi == 0 {
                canonical_count += t.canonical_count;
            }
            a.sum_x += t.sum_x[hi];
            a.sum_y += t.sum_y[hi];
            a.sum_xy += t.sum_xy[hi];
            a.sum_x2 += t.sum_x2[hi];
            a.sum_y2 += t.sum_y2[hi];
            a.serial_n += t.serial_n[hi];
            a.order_preserved += t.order_preserved[hi];
            for b in 0..256 {
                a.gap_bins[b] += t.gap_bins[hi][b];
                a.mutation_diff_bins[b] += t.mutation_diff_bins[hi][b];
            }
            a.hamming_sum += t.hamming_sum[hi];
            a.hamming_sum_sq += t.hamming_sum_sq[hi];
            a.hamming_min = a.hamming_min.min(t.hamming_min[hi]);
            a.hamming_max = a.hamming_max.max(t.hamming_max[hi]);
            a.mutation_count += t.mutation_count[hi];
        }
        hash_accums.push(a);
    }

    let result = LinearTaskResult {
        task_id,
        canonical_count,
        elapsed_seconds: start.elapsed().as_secs_f64(),
        hash_accums,
    };

    let out_path = linear_dir.join(format!("linear_task_{task_id:02}.json"));
    write_json_pretty(&out_path, &result)?;

    eprintln!(
        "[linear-map:{}] done in {:.1}s, {} canonical k-mers, wrote {}",
        task_id,
        result.elapsed_seconds,
        canonical_count,
        out_path.display()
    );

    Ok(())
}

fn cmd_linear_aggregate(output_dir: &Path) -> Result<()> {
    let linear_dir = output_dir.join("linear_test");

    // Load all task results
    let mut tasks = Vec::with_capacity(LINEAR_TASKS);
    for task_id in 0..LINEAR_TASKS {
        let path = linear_dir.join(format!("linear_task_{task_id:02}.json"));
        let result: LinearTaskResult =
            read_json(&path).with_context(|| format!("failed to read {}", path.display()))?;
        tasks.push(result);
    }

    let total_canonical: u64 = tasks.iter().map(|t| t.canonical_count).sum();

    // Aggregate per hash function
    let mut results = Vec::with_capacity(HASH_COUNT);
    for hi in 0..HASH_COUNT {
        let mut total_sx = 0.0f64;
        let mut total_sy = 0.0f64;
        let mut total_sxy = 0.0f64;
        let mut total_sx2 = 0.0f64;
        let mut total_sy2 = 0.0f64;
        let mut total_sn = 0u64;
        let mut total_order = 0u64;
        let mut total_gap_bins = [0u64; 256];
        let mut total_hamming_sum = 0u64;
        let mut total_hamming_sum_sq = 0u64;
        let mut total_hamming_min = 64u32;
        let mut total_hamming_max = 0u32;
        let mut total_mutation_count = 0u64;
        let mut total_mut_diff_bins = [0u64; 256];

        for task in &tasks {
            let a = &task.hash_accums[hi];
            total_sx += a.sum_x;
            total_sy += a.sum_y;
            total_sxy += a.sum_xy;
            total_sx2 += a.sum_x2;
            total_sy2 += a.sum_y2;
            total_sn += a.serial_n;
            total_order += a.order_preserved;
            for b in 0..256 {
                total_gap_bins[b] += a.gap_bins[b];
                total_mut_diff_bins[b] += a.mutation_diff_bins[b];
            }
            total_hamming_sum += a.hamming_sum;
            total_hamming_sum_sq += a.hamming_sum_sq;
            total_hamming_min = total_hamming_min.min(a.hamming_min);
            total_hamming_max = total_hamming_max.max(a.hamming_max);
            total_mutation_count += a.mutation_count;
        }

        let hash_name = tasks[0].hash_accums[hi].hash_name.clone();

        // Pearson correlation
        let n = total_sn as f64;
        let serial_correlation = if total_sn > 1 {
            let num = n * total_sxy - total_sx * total_sy;
            let den_x = (n * total_sx2 - total_sx * total_sx).sqrt();
            let den_y = (n * total_sy2 - total_sy * total_sy).sqrt();
            if den_x > 0.0 && den_y > 0.0 {
                num / (den_x * den_y)
            } else {
                0.0
            }
        } else {
            0.0
        };

        let order_frac = if total_sn > 0 {
            total_order as f64 / total_sn as f64
        } else {
            0.0
        };

        let gap_expected = total_sn as f64 / 256.0;
        let gap_chi: f64 = if gap_expected > 0.0 {
            total_gap_bins
                .iter()
                .map(|&c| {
                    let d = c as f64 - gap_expected;
                    d * d / gap_expected
                })
                .sum()
        } else {
            0.0
        };

        let mc = total_mutation_count as f64;
        let avg_hamming = if mc > 0.0 {
            total_hamming_sum as f64 / mc
        } else {
            0.0
        };
        let hamming_std = if mc > 1.0 {
            let var = total_hamming_sum_sq as f64 / mc - avg_hamming * avg_hamming;
            var.max(0.0).sqrt()
        } else {
            0.0
        };

        let mut_expected = total_mutation_count as f64 / 256.0;
        let mut_diff_chi: f64 = if mut_expected > 0.0 {
            total_mut_diff_bins
                .iter()
                .map(|&c| {
                    let d = c as f64 - mut_expected;
                    d * d / mut_expected
                })
                .sum()
        } else {
            0.0
        };

        results.push(LinearHashMetrics {
            hash_name,
            serial_correlation,
            order_preservation_frac: order_frac,
            gap_chi_squared: gap_chi,
            mutation_avg_hamming: avg_hamming,
            mutation_hamming_std: hamming_std,
            mutation_min_hamming: total_hamming_min,
            mutation_max_hamming: total_hamming_max,
            mutation_diff_chi_squared: mut_diff_chi,
        });
    }

    let report = LinearReport {
        total_canonical,
        tasks: LINEAR_TASKS,
        hashes: results,
    };

    let json_path = linear_dir.join("linear_test.json");
    write_json_pretty(&json_path, &report)?;

    eprintln!("\n{:=<100}", "");
    eprintln!(
        "{:<10} {:>14} {:>10} {:>12} {:>10} {:>10} {:>6} {:>6} {:>12}",
        "Hash", "SerialCorr", "OrderFrac", "GapChi²", "AvgHamm", "StdHamm", "MinH", "MaxH", "MutDiffChi²"
    );
    eprintln!("{:-<100}", "");
    for m in &report.hashes {
        eprintln!(
            "{:<10} {:>14.8} {:>10.6} {:>12.1} {:>10.4} {:>10.4} {:>6} {:>6} {:>12.1}",
            m.hash_name,
            m.serial_correlation,
            m.order_preservation_frac,
            m.gap_chi_squared,
            m.mutation_avg_hamming,
            m.mutation_hamming_std,
            m.mutation_min_hamming,
            m.mutation_max_hamming,
            m.mutation_diff_chi_squared,
        );
    }
    eprintln!("{:=<100}", "");
    eprintln!(
        "[linear-aggregate] total_canonical={}, wrote {}",
        total_canonical,
        json_path.display()
    );

    Ok(())
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
    let max_open_writers_per_worker = configured_map_max_open_writers_per_worker();
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

    let ranges = map_task_raw_ranges(task_id);
    let start = Instant::now();

    let raw_start = ranges[0].start;
    let raw_end = ranges[0].end;

    eprintln!(
        "[map:{}] raw_ranges=[{}..{},{}..{}] compute_threads={} io_threads={} local_buf_values={} file_buf_mb={} channel_entries={} channel_chunks={} max_open_writers_per_worker={} partitions={} scratch={}",
        task_id,
        ranges[0].start,
        ranges[0].end,
        ranges[1].start,
        ranges[1].end,
        compute_threads,
        io_threads,
        local_buffer_values,
        map_file_buffer_bytes / (1024 * 1024),
        channel_entries,
        channel_capacity,
        max_open_writers_per_worker,
        PARTITIONS,
        local_task_scratch.display()
    );

    let (io_senders, io_handles) = start_map_io_workers(
        &local_task_scratch,
        task_id,
        io_threads,
        channel_capacity,
        map_file_buffer_bytes,
        max_open_writers_per_worker,
    )?;
    let io_senders = Arc::new(io_senders);

    let mut subranges = Vec::new();
    for range in &ranges {
        subranges.extend(split_range(range.clone(), compute_threads.max(1)));
    }

    let worker_results: Vec<Result<WorkerSummary>> = subranges
        .into_par_iter()
        .map(|sub| process_subrange(sub, io_senders.clone(), io_threads, local_buffer_values))
        .collect();

    let mut combined = WorkerSummary::default();
    let mut compute_error: Option<anyhow::Error> = None;
    for result in worker_results {
        match result {
            Ok(ws) => {
                combined.canonical_count += ws.canonical_count;
                for h in 0..HASH_COUNT {
                    for p in 0..PARTITIONS {
                        combined.partition_counts[h][p] += ws.partition_counts[h][p];
                    }
                }
            }
            Err(err) => {
                if compute_error.is_none() {
                    compute_error = Some(err);
                }
            }
        }
    }

    drop(io_senders);

    let mut io_error: Option<anyhow::Error> = None;
    for handle in io_handles {
        match handle.join() {
            Ok(Ok(())) => {}
            Ok(Err(err)) => {
                if io_error.is_none() {
                    io_error = Some(err);
                }
            }
            Err(_) => {
                if io_error.is_none() {
                    io_error = Some(anyhow::anyhow!("map io worker panicked"));
                }
            }
        }
    }

    if let Some(err) = io_error {
        return Err(err.context("map io worker failed"));
    }
    if let Some(err) = compute_error {
        return Err(err.context("map compute worker failed"));
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
        raw_start,
        raw_end,
        raw_ranges: vec![
            [ranges[0].start, ranges[0].end],
            [ranges[1].start, ranges[1].end],
        ],
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
    let threads = configured_threads();
    configure_rayon(threads);

    let chunk_values = configured_sort_chunk_values();
    let sort_scratch_base = configured_sort_scratch_base();
    let sort_task_scratch = sort_task_scratch_dir(&sort_scratch_base, task_id);

    if sort_task_scratch.exists() {
        remove_dir_all(&sort_task_scratch).with_context(|| {
            format!(
                "failed to remove previous sort scratch {}",
                sort_task_scratch.display()
            )
        })?;
    }
    create_dir_all(&sort_task_scratch).with_context(|| {
        format!(
            "failed to create sort scratch {}",
            sort_task_scratch.display()
        )
    })?;
    let mut scratch_guard = ScratchCleanupGuard::new(sort_task_scratch.clone());

    let read_start = Instant::now();
    let mut sort_secs = 0.0f64;
    let mut loaded = 0u64;
    let mut run_idx = 0usize;
    let mut run_paths = Vec::new();
    let mut chunk = Vec::with_capacity(chunk_values);

    for map_task in 0..plan.map_tasks {
        let path = shard_path(output_dir, hash_name, map_task, partition_id);
        read_u64_file_into_runs(
            &path,
            &mut chunk,
            chunk_values,
            &mut run_paths,
            &sort_task_scratch,
            &mut run_idx,
            &mut loaded,
            &mut sort_secs,
        )
        .with_context(|| format!("failed to read {}", path.display()))?;
    }

    spill_sorted_run(
        &mut chunk,
        &mut run_paths,
        &sort_task_scratch,
        &mut run_idx,
        &mut sort_secs,
    )?;

    let read_secs = read_start.elapsed().as_secs_f64();

    if loaded != expected {
        bail!(
            "{} part {}: loaded {} hashes, expected {}",
            hash_name,
            partition_id,
            loaded,
            expected
        );
    }

    let merge_start = Instant::now();
    let sorted_path = sorted_partition_path(output_dir, hash_name, partition_id);
    let stats = merge_runs_to_sorted_output(
        &run_paths,
        &sorted_path,
        hash_name,
        partition_id,
        global_offset,
        plan.total_canonical,
    )?;
    let stats_secs = merge_start.elapsed().as_secs_f64();

    if stats.count != expected {
        bail!(
            "{} part {}: merged {} hashes, expected {}",
            hash_name,
            partition_id,
            stats.count,
            expected
        );
    }

    remove_dir_all(&sort_task_scratch).with_context(|| {
        format!(
            "failed to cleanup sort scratch {}",
            sort_task_scratch.display()
        )
    })?;
    scratch_guard.disarm();

    let part_stats_path = partition_stats_path(output_dir, hash_name, partition_id);
    write_json_pretty(&part_stats_path, &stats)?;

    eprintln!(
        "[sort-reduce:{}:{}] n={} runs={} read={:.2}s sort={:.2}s merge+stats={:.2}s chunk_values={}",
        hash_name,
        partition_id,
        loaded,
        run_paths.len(),
        read_secs,
        sort_secs,
        stats_secs,
        chunk_values,
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
    max_open_writers_per_worker: usize,
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
                    max_open_writers_per_worker,
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
    max_open_writers_per_worker: usize,
    rx: Receiver<WriteChunk>,
) -> Result<()> {
    let mut writers = (0..(HASH_COUNT * PARTITIONS))
        .map(|_| None)
        .collect::<Vec<Option<BufWriter<File>>>>();

    let mut open_order: Vec<usize> = Vec::new();
    let max_open = max_open_writers_per_worker.max(1);

    while let Ok(chunk) = rx.recv() {
        if chunk.stream_idx % worker_count != worker_idx {
            bail!(
                "received chunk for stream {} in wrong io worker {}",
                chunk.stream_idx,
                worker_idx
            );
        }

        if writers[chunk.stream_idx].is_none() {
            if open_order.len() >= max_open {
                let evict_stream = open_order.remove(0);
                if let Some(mut writer) = writers[evict_stream].take() {
                    writer.flush().with_context(|| {
                        format!("failed to flush evicted writer stream {}", evict_stream)
                    })?;
                }
            }

            let (hash_name, partition) = stream_to_hash_partition(chunk.stream_idx);
            let path = shard_path(&base, hash_name, map_task, partition);
            if let Some(parent) = path.parent() {
                create_dir_all(parent)
                    .with_context(|| format!("failed to create {}", parent.display()))?;
            }
            let file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(&path)
                .with_context(|| format!("failed to open {}", path.display()))?;
            writers[chunk.stream_idx] = Some(BufWriter::with_capacity(map_file_buffer_bytes, file));
        }

        open_order.retain(|&s| s != chunk.stream_idx);
        open_order.push(chunk.stream_idx);

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

            run_cpb_with_retries(&src, &dst)?;
        }
    }

    Ok(())
}

fn read_u64_file_into_runs(
    path: &Path,
    chunk: &mut Vec<u64>,
    chunk_values: usize,
    run_paths: &mut Vec<PathBuf>,
    scratch_dir: &Path,
    run_idx: &mut usize,
    loaded: &mut u64,
    sort_secs: &mut f64,
) -> Result<()> {
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
            chunk.push(u64::from_le(value));
            *loaded += 1;

            if chunk.len() >= chunk_values {
                spill_sorted_run(chunk, run_paths, scratch_dir, run_idx, sort_secs)?;
            }

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

fn spill_sorted_run(
    chunk: &mut Vec<u64>,
    run_paths: &mut Vec<PathBuf>,
    scratch_dir: &Path,
    run_idx: &mut usize,
    sort_secs: &mut f64,
) -> Result<()> {
    if chunk.is_empty() {
        return Ok(());
    }

    let sort_start = Instant::now();
    chunk.par_sort_unstable();
    *sort_secs += sort_start.elapsed().as_secs_f64();

    let path = scratch_dir.join(format!("run_{:04}.bin", *run_idx));
    *run_idx += 1;
    write_u64_file(&path, chunk)?;
    run_paths.push(path);
    chunk.clear();

    Ok(())
}

fn merge_runs_to_sorted_output(
    run_paths: &[PathBuf],
    sorted_path: &Path,
    hash_name: &str,
    partition_id: usize,
    partition_offset: u64,
    total_hashes: u64,
) -> Result<PartitionStats> {
    if run_paths.is_empty() {
        write_u64_file(sorted_path, &[])?;
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

    let mut streams = Vec::with_capacity(run_paths.len());
    for path in run_paths {
        streams.push(RunStream::open(path)?);
    }

    let file = File::create(sorted_path)
        .with_context(|| format!("failed to create {}", sorted_path.display()))?;
    let mut writer = BufWriter::with_capacity(IO_BUFFER_BYTES, file);
    let mut out_values = Vec::<u64>::with_capacity(SORT_OUTPUT_BUFFER_VALUES);

    let mut heap = BinaryHeap::<MergeHeapItem>::new();
    for (run_idx, stream) in streams.iter_mut().enumerate() {
        if let Some(value) = stream.next_u64()? {
            heap.push(MergeHeapItem { value, run_idx });
        }
    }

    let mut acc =
        SortStatsAccumulator::new(hash_name, partition_id, partition_offset, total_hashes);

    while let Some(item) = heap.pop() {
        out_values.push(item.value);
        acc.observe(item.value);

        if out_values.len() >= SORT_OUTPUT_BUFFER_VALUES {
            writer.write_all(cast_slice::<u64, u8>(&out_values))?;
            out_values.clear();
        }

        if let Some(next_value) = streams[item.run_idx].next_u64()? {
            heap.push(MergeHeapItem {
                value: next_value,
                run_idx: item.run_idx,
            });
        }
    }

    if !out_values.is_empty() {
        writer.write_all(cast_slice::<u64, u8>(&out_values))?;
    }
    writer.flush()?;

    Ok(acc.finalize())
}

fn run_cpb_with_retries(src: &Path, dst: &Path) -> Result<()> {
    let mut last_status = None;
    let mut last_spawn_err: Option<std::io::Error> = None;

    let max_attempts = CPB_STATUS_RETRIES.max(TRANSIENT_IO_RETRIES).max(1);

    for attempt in 1..=max_attempts {
        match ProcessCommand::new("cpb").arg(src).arg(dst).status() {
            Ok(status) if status.success() => return Ok(()),
            Ok(status) => {
                last_status = Some(status);
                if attempt < max_attempts {
                    let backoff_ms = retry_backoff_ms(attempt as u32);
                    eprintln!(
                        "[copy] cpb retry {}/{} status={} src={} dst={}",
                        attempt,
                        max_attempts,
                        status,
                        src.display(),
                        dst.display()
                    );
                    thread::sleep(Duration::from_millis(backoff_ms));
                    continue;
                }
            }
            Err(err) => {
                if is_transient_spawn_error(&err) && attempt < max_attempts {
                    let backoff_ms = retry_backoff_ms(attempt as u32);
                    eprintln!(
                        "[copy] cpb spawn retry {}/{} os_error={:?} src={} dst={}",
                        attempt,
                        max_attempts,
                        err.raw_os_error(),
                        src.display(),
                        dst.display()
                    );
                    thread::sleep(Duration::from_millis(backoff_ms));
                    last_spawn_err = Some(err);
                    continue;
                }

                return Err(err).with_context(|| {
                    format!(
                        "failed to execute cpb from {} to {}",
                        src.display(),
                        dst.display()
                    )
                });
            }
        }
    }

    if let Some(status) = last_status {
        bail!(
            "cpb failed after {} attempts with status {} for {} -> {}",
            max_attempts,
            status,
            src.display(),
            dst.display()
        );
    }

    if let Some(err) = last_spawn_err {
        return Err(err).with_context(|| {
            format!(
                "cpb spawn failed after {} attempts for {} -> {}",
                max_attempts,
                src.display(),
                dst.display()
            )
        });
    }

    bail!(
        "cpb failed after {} attempts for {} -> {}",
        max_attempts,
        src.display(),
        dst.display()
    );
}

fn retry_backoff_ms(attempt: u32) -> u64 {
    let shift = attempt.saturating_sub(1).min(6);
    TRANSIENT_IO_SLEEP_MS.saturating_mul(1u64 << shift)
}

fn is_transient_spawn_error(err: &std::io::Error) -> bool {
    matches!(err.raw_os_error(), Some(11) | Some(4) | Some(35))
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

fn configured_map_max_open_writers_per_worker() -> usize {
    std::env::var("K21_IO_MAX_OPEN_WRITERS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&n| n > 0)
        .unwrap_or(DEFAULT_MAP_MAX_OPEN_WRITERS_PER_WORKER)
}

fn configured_local_scratch_base() -> PathBuf {
    std::env::var("K21_LOCAL_SCRATCH_BASE")
        .ok()
        .filter(|v| !v.trim().is_empty())
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(DEFAULT_LOCAL_SCRATCH_BASE))
}

fn configured_sort_chunk_values() -> usize {
    if let Some(values) = std::env::var("K21_SORT_CHUNK_VALUES")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&n| n > 0)
    {
        return values;
    }

    let gib = std::env::var("K21_SORT_CHUNK_GIB")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&n| n > 0)
        .unwrap_or(DEFAULT_SORT_CHUNK_GIB);

    gib.saturating_mul(1024 * 1024 * 1024) / std::mem::size_of::<u64>()
}

fn configured_sort_scratch_base() -> PathBuf {
    std::env::var("K21_SORT_SCRATCH_BASE")
        .ok()
        .filter(|v| !v.trim().is_empty())
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(DEFAULT_SORT_SCRATCH_BASE))
}

fn map_task_scratch_dir(base: &Path, task_id: usize) -> PathBuf {
    let job_id = std::env::var("SLURM_JOB_ID").unwrap_or_else(|_| "manual".to_string());
    base.join(format!("k21-map-{job_id}-task-{task_id:02}"))
}

fn sort_task_scratch_dir(base: &Path, task_id: usize) -> PathBuf {
    let job_id = std::env::var("SLURM_JOB_ID").unwrap_or_else(|_| "manual".to_string());
    base.join(format!("k21-sort-{job_id}-task-{task_id:03}"))
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

fn map_task_raw_ranges(task_id: usize) -> [Range<u64>; 2] {
    let chunks = (MAP_TASKS * 2) as u64;
    let left_idx = task_id as u64;
    let right_idx = chunks - 1 - task_id as u64;

    [
        get_chunk_range(RAW_SPACE_SIZE, chunks, left_idx),
        get_chunk_range(RAW_SPACE_SIZE, chunks, right_idx),
    ]
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
