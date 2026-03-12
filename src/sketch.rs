use crate::bias::HashBiasTable;
use crate::core_utils::passes_entropy_filter;
use crate::format::{BUCKET_COUNT, ENTRY_SIZE, Entry, bucket_id};
use crate::io::EntryWriter;
use crossfire::mpsc;
use crossfire::{MTx, Rx};
use dashmap::DashMap;
use indicatif::{ProgressBar, ProgressStyle};
use jamhash::jamhash_u64;
use memmap2::Mmap;
use needletail::{Sequence, parse_fastx_reader};
use rayon::prelude::*;
use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicU32, AtomicU64, Ordering};
use std::time::Duration;
use tempfile::TempDir;

const WRITE_BUFFER_SIZE: usize = 8 * 1024 * 1024;
const MIN_MEMORY_GB: usize = 4;
const DEFAULT_SEND_TIMEOUT: Duration = Duration::from_millis(1);
const MIN_SPLIT_SIZE: usize = 1024 * 1024;
const MAX_CONCURRENT_MMAPS: usize = 256;
const OPTIMAL_CHANNEL_CAPACITY: usize = 512 * 1024;

type Sender = MTx<mpsc::Array<Entry>>;
type Receiver = Rx<mpsc::Array<Entry>>;

#[derive(Clone)]
pub struct SketchConfig {
    pub kmer_size: u8,
    pub fscale: u64,
    pub num_threads: usize,
    pub memory: usize,
    pub temp_dir_base: Option<PathBuf>,
    pub min_entropy: f64,
    pub singleton: bool,
    pub bias_table: Option<Arc<HashBiasTable>>,
    pub send_timeout: Duration,
    pub show_progress: bool,
}

impl Default for SketchConfig {
    fn default() -> Self {
        Self {
            kmer_size: 21,
            fscale: 1000,
            num_threads: 1,
            memory: MIN_MEMORY_GB,
            temp_dir_base: None,
            min_entropy: 0.0,
            singleton: false,
            bias_table: None,
            send_timeout: DEFAULT_SEND_TIMEOUT,
            show_progress: false,
        }
    }
}

pub struct SketchResult {
    pub sample_count: u32,
    pub bucket_entry_counts: [u64; BUCKET_COUNT],
    pub frac_max: u64,
    pub temp_dir: TempDir,
    pub sample_names: Vec<String>,
}

#[derive(Debug, thiserror::Error)]
pub enum SketchError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Parse error in {path}: {message}")]
    Parse { path: PathBuf, message: String },

    #[error("Channel send error")]
    Channel,

    #[error("Invalid configuration: {0}")]
    Config(String),
}

struct WorkUnit {
    mmap: Option<Arc<Mmap>>,
    start: usize,
    end: usize,
    sample_id: Option<u32>,
    source_path: Arc<PathBuf>,
}

struct BucketWriter {
    receiver: Receiver,
    writer: EntryWriter,
    bucket_id: usize,
}

impl BucketWriter {
    fn drain(&mut self) -> io::Result<()> {
        while let Ok(entry) = self.receiver.try_recv() {
            self.writer.write(&entry)?;
        }
        Ok(())
    }

    fn drain_until_disconnected(&mut self, timeout: Duration) -> io::Result<bool> {
        match self.receiver.recv_timeout(timeout) {
            Ok(entry) => {
                self.writer.write(&entry)?;
                Ok(true)
            }
            Err(crossfire::RecvTimeoutError::Timeout) => Ok(true),
            Err(crossfire::RecvTimeoutError::Disconnected) => Ok(false),
        }
    }
}

struct SketchContext<'a> {
    senders: &'a [Sender],
    config: &'a SketchConfig,
    sample_counter: &'a AtomicU32,
    frac_max: u64,
    sample_names: &'a DashMap<u32, String>,
}

struct MmapSliceReader {
    mmap: Arc<Mmap>,
    start: usize,
    end: usize,
    pos: usize,
}

impl MmapSliceReader {
    fn new(mmap: Arc<Mmap>, start: usize, end: usize) -> Self {
        Self {
            mmap,
            start,
            end,
            pos: 0,
        }
    }
}

impl io::Read for MmapSliceReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let current = self.start + self.pos;
        if current >= self.end {
            return Ok(0);
        }
        let remaining = &self.mmap[current..self.end];
        let n = remaining.len().min(buf.len());
        buf[..n].copy_from_slice(&remaining[..n]);
        self.pos += n;
        Ok(n)
    }
}

fn distribute_evenly(total: usize, parts: usize) -> impl Iterator<Item = usize> {
    let per_part = total / parts;
    let remainder = total % parts;
    (0..parts).map(move |i| {
        if i < remainder {
            per_part + 1
        } else {
            per_part
        }
    })
}

const GZ_MAGIC: [u8; 2] = [0x1F, 0x8B];
const BZ_MAGIC: [u8; 2] = [0x42, 0x5A];
const XZ_MAGIC: [u8; 2] = [0xFD, 0x37];
const ZST_MAGIC: [u8; 2] = [0x28, 0xB5];

#[inline]
fn is_compressed(magic: [u8; 2]) -> bool {
    matches!(magic, GZ_MAGIC | BZ_MAGIC | XZ_MAGIC | ZST_MAGIC)
}

fn validate_format(magic: [u8; 2], path: &Path) -> Result<(), SketchError> {
    if !is_compressed(magic) && !matches!(magic[0], b'>' | b'@') {
        return Err(SketchError::Parse {
            path: path.to_path_buf(),
            message: format!(
                "unrecognized format (bytes: [{:#04X}, {:#04X}])",
                magic[0], magic[1]
            ),
        });
    }
    Ok(())
}

fn validate_file_header(path: &Path) -> Result<[u8; 2], SketchError> {
    let mut file = File::open(path)?;
    let mut magic = [0u8; 2];
    use std::io::Read;
    let bytes_read = file.read(&mut magic)?;
    if bytes_read < 2 {
        return Err(SketchError::Parse {
            path: path.to_path_buf(),
            message: format!("file too small ({bytes_read} bytes) to be valid FASTA/FASTQ"),
        });
    }
    validate_format(magic, path)?;
    Ok(magic)
}

fn validate_mmap_header(mmap: &Mmap, path: &Path) -> Result<[u8; 2], SketchError> {
    if mmap.len() < 2 {
        return Err(SketchError::Parse {
            path: path.to_path_buf(),
            message: format!(
                "file too small ({} bytes) to be valid FASTA/FASTQ",
                mmap.len()
            ),
        });
    }
    let magic = [mmap[0], mmap[1]];
    validate_format(magic, path)?;
    Ok(magic)
}

fn scan_fasta_boundaries(data: &[u8]) -> Vec<usize> {
    let mut bounds = vec![0];
    bounds.extend(
        data.windows(2)
            .enumerate()
            .filter_map(|(i, w)| (w == b"\n>").then_some(i + 1)),
    );
    bounds
}

#[inline]
fn is_iupac_nucleotide(b: u8) -> bool {
    matches!(
        b | 0x20,
        b'a' | b'c'
            | b'g'
            | b't'
            | b'u'
            | b'n'
            | b'r'
            | b'y'
            | b's'
            | b'w'
            | b'k'
            | b'm'
            | b'b'
            | b'd'
            | b'h'
            | b'v'
    )
}

fn scan_fastq_boundaries(data: &[u8]) -> Vec<usize> {
    let mut bounds = vec![0];
    let mut i = 0;

    while i + 1 < data.len() {
        if data[i] == b'\n' && data[i + 1] == b'@' {
            let header_start = i + 1;
            let header_end = data[header_start..]
                .iter()
                .position(|&b| b == b'\n')
                .map(|p| header_start + p)
                .unwrap_or(data.len());

            if header_end > header_start + 1 {
                let seq_start = header_end + 1;
                if seq_start < data.len() && is_iupac_nucleotide(data[seq_start]) {
                    bounds.push(header_start);
                }
            }
        }
        i += 1;
    }

    bounds
}

fn setup_channels(
    num_threads: usize,
    memory_gb: usize,
    input_size_bytes: u64,
    temp_path: &Path,
) -> Result<(Vec<Sender>, Vec<Vec<BucketWriter>>), SketchError> {
    let capacity = compute_channel_capacity(memory_gb, input_size_bytes);

    let (senders, receivers): (Vec<_>, Vec<_>) = (0..BUCKET_COUNT)
        .map(|_| mpsc::bounded_blocking(capacity))
        .unzip();

    let bucket_threads = num_threads.min(BUCKET_COUNT);
    let chunk_sizes = distribute_evenly(BUCKET_COUNT, bucket_threads);

    let mut rx_iter = receivers.into_iter().enumerate();
    let mut bucket_writers: Vec<Vec<BucketWriter>> = chunk_sizes
        .map(|count| {
            rx_iter
                .by_ref()
                .take(count)
                .map(|(bucket_id, receiver)| {
                    let writer = EntryWriter::new(
                        temp_path.join(format!("bucket_{bucket_id:03}.bin")),
                        WRITE_BUFFER_SIZE,
                    )?;
                    Ok(BucketWriter {
                        receiver,
                        writer,
                        bucket_id,
                    })
                })
                .collect::<Result<Vec<_>, std::io::Error>>()
        })
        .collect::<Result<Vec<_>, _>>()?;

    bucket_writers.resize_with(num_threads, Vec::new);

    Ok((senders, bucket_writers))
}

fn compute_channel_capacity(memory_gb: usize, input_size_bytes: u64) -> usize {
    let memory_bytes = memory_gb as u64 * 1024 * 1024 * 1024;
    let writer_memory = BUCKET_COUNT as u64 * WRITE_BUFFER_SIZE as u64;

    let available = memory_bytes
        .saturating_sub(input_size_bytes)
        .saturating_sub(writer_memory);

    let computed = (available / (BUCKET_COUNT as u64 * ENTRY_SIZE as u64)) as usize;

    computed.clamp(1024, OPTIMAL_CHANNEL_CAPACITY)
}

fn scan_boundaries(mmap: &Mmap, magic: [u8; 2]) -> Vec<usize> {
    if mmap.len() < MIN_SPLIT_SIZE {
        return vec![0];
    }
    match magic[0] {
        b'>' => scan_fasta_boundaries(mmap),
        _ => scan_fastq_boundaries(mmap),
    }
}

fn distribute_work_units(
    positions: Vec<(Arc<Mmap>, Arc<PathBuf>, usize, usize)>,
    num_threads: usize,
    singleton: bool,
    sample_counter: &AtomicU32,
    thread_work: &mut [Vec<WorkUnit>],
) {
    if positions.is_empty() {
        return;
    }

    let mut file_sample_ids: std::collections::HashMap<PathBuf, u32> =
        std::collections::HashMap::new();

    let chunk_sizes: Vec<_> = distribute_evenly(positions.len(), num_threads).collect();
    let mut offset = 0;

    for (t, &count) in chunk_sizes.iter().enumerate() {
        for (mmap, path, start_byte, end_byte) in &positions[offset..offset + count] {
            let sample_id = (!singleton).then(|| {
                *file_sample_ids
                    .entry((**path).clone())
                    .or_insert_with(|| sample_counter.fetch_add(1, Ordering::SeqCst))
            });

            thread_work[t].push(WorkUnit {
                mmap: Some(Arc::clone(mmap)),
                start: *start_byte,
                end: *end_byte,
                sample_id,
                source_path: Arc::clone(path),
            });
        }
        offset += count;
    }
}

struct WorkUnitResult {
    thread_work: Vec<Vec<WorkUnit>>,
    total_input_bytes: u64,
}

fn build_work_units(
    input_files: &[PathBuf],
    num_threads: usize,
    singleton: bool,
    memory_gb: usize,
    sample_counter: &AtomicU32,
    show_progress: bool,
) -> Result<WorkUnitResult, SketchError> {
    let mut thread_work: Vec<Vec<WorkUnit>> = (0..num_threads).map(|_| Vec::new()).collect();
    let next_sample = || (!singleton).then(|| sample_counter.fetch_add(1, Ordering::SeqCst));

    let total_input_bytes: u64 = input_files
        .iter()
        .filter_map(|p| std::fs::metadata(p).ok())
        .map(|m| m.len())
        .sum();

    let memory_bytes = memory_gb as u64 * 1024 * 1024 * 1024;

    let skip_mmap = input_files.len() > MAX_CONCURRENT_MMAPS || total_input_bytes > memory_bytes;

    if skip_mmap {
        let validation_pb = if show_progress {
            let pb = ProgressBar::new(input_files.len() as u64);
            pb.set_style(
                ProgressStyle::default_bar()
                    .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} files validated")
                    .unwrap()
                    .progress_chars("#>-"),
            );
            Some(pb)
        } else {
            None
        };

        let validation_results: Vec<Result<(), SketchError>> = input_files
            .par_iter()
            .map(|path| {
                let result = validate_file_header(path).map(|_| ());
                if let Some(ref pb) = validation_pb {
                    pb.inc(1);
                }
                result
            })
            .collect();

        if let Some(pb) = validation_pb {
            pb.finish_with_message("validation complete");
        }

        for result in validation_results {
            result?;
        }

        for (i, path) in input_files.iter().enumerate() {
            thread_work[i % num_threads].push(WorkUnit {
                mmap: None,
                start: 0,
                end: 0,
                sample_id: next_sample(),
                source_path: Arc::new(path.clone()),
            });
        }
        return Ok(WorkUnitResult {
            thread_work,
            total_input_bytes: 0,
        });
    }

    let mut flat_positions: Vec<(Arc<Mmap>, Arc<PathBuf>, usize, usize)> = Vec::new();
    let mut compressed_files: Vec<(Arc<Mmap>, Arc<PathBuf>)> = Vec::new();

    for path in input_files {
        let file = File::open(path)?;
        let mmap = Arc::new(unsafe { Mmap::map(&file)? });
        let path = Arc::new(path.clone());
        let magic = validate_mmap_header(&mmap, &path)?;

        if is_compressed(magic) {
            compressed_files.push((mmap, path));
            continue;
        }

        let boundaries = scan_boundaries(&mmap, magic);
        for (i, &start) in boundaries.iter().enumerate() {
            let end = boundaries.get(i + 1).copied().unwrap_or(mmap.len());
            flat_positions.push((Arc::clone(&mmap), Arc::clone(&path), start, end));
        }
    }

    distribute_work_units(
        flat_positions,
        num_threads,
        singleton,
        sample_counter,
        &mut thread_work,
    );

    for (i, (mmap, path)) in compressed_files.into_iter().enumerate() {
        let end = mmap.len();
        thread_work[i % num_threads].push(WorkUnit {
            mmap: Some(mmap),
            start: 0,
            end,
            sample_id: next_sample(),
            source_path: path,
        });
    }

    Ok(WorkUnitResult {
        thread_work,
        total_input_bytes,
    })
}

pub fn run(input_files: &[PathBuf], config: &SketchConfig) -> Result<SketchResult, SketchError> {
    if config.fscale == 0 {
        return Err(SketchError::Config("fscale must be non-zero".to_string()));
    }
    if config.kmer_size == 0 || config.kmer_size > 31 {
        return Err(SketchError::Config(format!(
            "kmer_size must be between 1 and 31, got {}",
            config.kmer_size
        )));
    }

    let temp_dir = match &config.temp_dir_base {
        Some(base) => tempfile::Builder::new().prefix("jam_").tempdir_in(base)?,
        None => tempfile::Builder::new().prefix("jam_").tempdir()?,
    };

    let frac_max = match &config.bias_table {
        Some(bt) if bt.is_soft_filter() => u64::MAX / bt.min_fscale(),
        _ => u64::MAX / config.fscale,
    };
    let sample_counter = Arc::new(AtomicU32::new(0));
    let num_threads = config.num_threads.max(1);

    let WorkUnitResult {
        thread_work,
        total_input_bytes,
    } = build_work_units(
        input_files,
        num_threads,
        config.singleton,
        config.memory,
        &sample_counter,
        config.show_progress,
    )?;
    let (senders, resources) = setup_channels(
        num_threads,
        config.memory,
        total_input_bytes,
        temp_dir.path(),
    )?;

    let (result_tx, result_rx) = std::sync::mpsc::channel();

    let sample_names_map: DashMap<u32, String> = DashMap::new();

    let total_files: u64 = thread_work.iter().map(|w| w.len() as u64).sum();
    let files_processed = Arc::new(AtomicU64::new(0));

    let progress_bar = if config.show_progress {
        let pb = ProgressBar::new(total_files);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} files ({msg})")
                .unwrap()
                .progress_chars("#>-"),
        );
        pb.enable_steady_tick(Duration::from_millis(100));
        Some(pb)
    } else {
        None
    };

    rayon::scope(|s| {
        let sample_counter = &sample_counter;
        let files_processed = &files_processed;
        let progress_bar = &progress_bar;
        let sample_names_map = &sample_names_map;
        for (work, res) in thread_work.into_iter().zip(resources) {
            let thread_senders = senders.to_vec();
            let result_tx = result_tx.clone();

            s.spawn(move |_| {
                let result = process_thread_work(
                    work,
                    thread_senders,
                    res,
                    config,
                    sample_counter,
                    frac_max,
                    files_processed,
                    progress_bar,
                    sample_names_map,
                );
                let _ = result_tx.send(result);
            });
        }
        drop(senders);
        drop(result_tx);
    });

    if let Some(pb) = progress_bar {
        let final_samples = sample_counter.load(Ordering::SeqCst);
        pb.finish_with_message(format!("{} samples", final_samples));
    }

    let mut bucket_entry_counts = [0u64; BUCKET_COUNT];
    for result in result_rx {
        for (bucket_idx, count) in result? {
            bucket_entry_counts[bucket_idx] = count;
        }
    }

    let sample_count = sample_counter.load(Ordering::SeqCst);
    let mut sample_names: Vec<String> = vec![String::new(); sample_count as usize];

    for entry in sample_names_map.iter() {
        if (*entry.key() as usize) < sample_names.len() {
            sample_names[*entry.key() as usize] = entry.value().clone();
        }
    }

    Ok(SketchResult {
        sample_count,
        bucket_entry_counts,
        frac_max,
        temp_dir,
        sample_names,
    })
}

#[allow(clippy::too_many_arguments)]
fn process_thread_work(
    work_units: Vec<WorkUnit>,
    senders: Vec<Sender>,
    mut bucket_writers: Vec<BucketWriter>,
    config: &SketchConfig,
    sample_counter: &AtomicU32,
    frac_max: u64,
    files_processed: &AtomicU64,
    progress_bar: &Option<ProgressBar>,
    sample_names_map: &DashMap<u32, String>,
) -> Result<Vec<(usize, u64)>, SketchError> {
    let ctx = SketchContext {
        senders: &senders,
        config,
        sample_counter,
        frac_max,
        sample_names: sample_names_map,
    };

    for unit in &work_units {
        let reader: Box<dyn io::Read + Send> = match &unit.mmap {
            Some(mmap) => Box::new(MmapSliceReader::new(Arc::clone(mmap), unit.start, unit.end)),
            None => Box::new(io::BufReader::new(File::open(&*unit.source_path)?)),
        };

        let mut fastx = match parse_fastx_reader(reader) {
            Ok(reader) => reader,
            Err(e) if e.kind == needletail::errors::ParseErrorKind::EmptyFile => {
                eprintln!(
                    "Empty file detected: {}, skipping",
                    unit.source_path.display()
                );
                let processed = files_processed.fetch_add(1, Ordering::Relaxed) + 1;
                if let Some(pb) = progress_bar {
                    pb.set_position(processed);
                }
                continue;
            }
            Err(e) => {
                return Err(SketchError::Parse {
                    path: (*unit.source_path).clone(),
                    message: e.to_string(),
                });
            }
        };

        sketch_records(
            fastx.as_mut(),
            unit.sample_id,
            &unit.source_path,
            &ctx,
            &mut bucket_writers,
        )?;

        let processed = files_processed.fetch_add(1, Ordering::Relaxed) + 1;
        if let Some(pb) = progress_bar {
            pb.set_position(processed);
            let samples = sample_counter.load(Ordering::Relaxed);
            pb.set_message(format!("{} samples", samples));
        }
    }

    drop(senders);

    const DRAIN_TIMEOUT: Duration = Duration::from_millis(1);
    loop {
        let mut any_pending = false;
        for bw in bucket_writers.iter_mut() {
            if bw.drain_until_disconnected(DRAIN_TIMEOUT)? {
                any_pending = true;
            }
        }
        if !any_pending {
            break;
        }
    }

    bucket_writers
        .iter_mut()
        .map(|bw| {
            bw.writer.flush()?;
            Ok((bw.bucket_id, bw.writer.count()))
        })
        .collect()
}

fn sketch_records(
    reader: &mut dyn needletail::FastxReader,
    file_sample_id: Option<u32>,
    source_path: &Path,
    ctx: &SketchContext,
    bucket_writers: &mut [BucketWriter],
) -> Result<(), SketchError> {
    let k = ctx.config.kmer_size;
    let min_entropy = ctx.config.min_entropy;
    let timeout = ctx.config.send_timeout;

    while let Some(record) = reader.next() {
        let record = record.map_err(|e| SketchError::Parse {
            path: source_path.to_path_buf(),
            message: e.to_string(),
        })?;

        let sample_id =
            file_sample_id.unwrap_or_else(|| ctx.sample_counter.fetch_add(1, Ordering::SeqCst));

        if !ctx.sample_names.contains_key(&sample_id) {
            let name = if file_sample_id.is_some() {
                source_path
                    .file_name()
                    .and_then(|s| s.to_str())
                    .unwrap_or("sample")
                    .to_string()
            } else {
                String::from_utf8_lossy(record.id()).to_string()
            };
            ctx.sample_names.insert(sample_id, name);
        }

        let sequence = record.normalize(false);
        if sequence.len() < k as usize {
            continue;
        }

        for (_, kmer, _) in sequence.bit_kmers(k, true) {
            let hash = jamhash_u64(kmer.0);

            if hash >= ctx.frac_max {
                continue;
            }

            if min_entropy > 0.0 && !passes_entropy_filter(kmer.0, k, min_entropy) {
                continue;
            }

            if ctx
                .config
                .bias_table
                .as_ref()
                .is_some_and(|b| !b.passes_filter(hash))
            {
                continue;
            }

            let entry = Entry::new(hash, sample_id);
            let bucket = bucket_id(hash);

            if let Err(crossfire::SendTimeoutError::Timeout(mut entry)) =
                ctx.senders[bucket].send_timeout(entry, timeout)
            {
                const MAX_RETRIES: u32 = 10;

                for retry in 0..MAX_RETRIES {
                    for bw in bucket_writers.iter_mut() {
                        bw.drain()?;
                    }

                    let backoff_sleep = Duration::from_micros(100 << retry.min(4));
                    std::thread::sleep(backoff_sleep);

                    let backoff_timeout = timeout.saturating_mul(1 << retry.min(4));
                    match ctx.senders[bucket].send_timeout(entry, backoff_timeout) {
                        Ok(()) => break,
                        Err(crossfire::SendTimeoutError::Timeout(e)) => {
                            entry = e;
                            if retry == MAX_RETRIES - 1 && ctx.senders[bucket].send(entry).is_err()
                            {
                                return Err(SketchError::Channel);
                            }
                        }
                        Err(crossfire::SendTimeoutError::Disconnected(_)) => {
                            return Err(SketchError::Channel);
                        }
                    }
                }
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_fasta(seqs: &[(&str, &str)]) -> NamedTempFile {
        let mut f = NamedTempFile::with_suffix(".fa").unwrap();
        for (name, seq) in seqs {
            writeln!(f, ">{name}").unwrap();
            writeln!(f, "{seq}").unwrap();
        }
        f
    }

    #[test]
    fn test_sketch_basic() {
        let input = make_fasta(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
        let config = SketchConfig {
            kmer_size: 11,
            fscale: 1,
            num_threads: 2,
            memory: 1,
            ..Default::default()
        };

        let result = run(&[input.path().to_path_buf()], &config).unwrap();
        assert_eq!(result.sample_count, 1);
        assert!(result.bucket_entry_counts.iter().sum::<u64>() > 0);
    }

    #[test]
    fn test_sketch_singleton() {
        let input = make_fasta(&[
            ("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
            ("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
        ]);
        let config = SketchConfig {
            kmer_size: 11,
            fscale: 1,
            singleton: true,
            num_threads: 2,
            memory: 1,
            ..Default::default()
        };

        let result = run(&[input.path().to_path_buf()], &config).unwrap();
        assert_eq!(result.sample_count, 2);
    }

    #[test]
    fn test_sketch_fracmin_filters() {
        let input = make_fasta(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);

        let result_all = run(
            &[input.path().to_path_buf()],
            &SketchConfig {
                kmer_size: 11,
                fscale: 1,
                memory: 1,
                ..Default::default()
            },
        )
        .unwrap();

        let result_filtered = run(
            &[input.path().to_path_buf()],
            &SketchConfig {
                kmer_size: 11,
                fscale: 100,
                memory: 1,
                ..Default::default()
            },
        )
        .unwrap();

        let total_all: u64 = result_all.bucket_entry_counts.iter().sum();
        let total_filtered: u64 = result_filtered.bucket_entry_counts.iter().sum();
        assert!(total_filtered < total_all);
    }

    #[test]
    fn test_sketch_multiple_files() {
        let input1 = make_fasta(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
        let input2 = make_fasta(&[("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")]);
        let config = SketchConfig {
            kmer_size: 11,
            fscale: 1,
            num_threads: 2,
            memory: 1,
            ..Default::default()
        };

        let result = run(
            &[input1.path().to_path_buf(), input2.path().to_path_buf()],
            &config,
        )
        .unwrap();
        assert_eq!(result.sample_count, 2);
        assert!(result.bucket_entry_counts.iter().sum::<u64>() > 0);
    }

    #[test]
    fn test_sketch_backpressure() {
        let input = make_fasta(&[(
            "seq1",
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        )]);
        let config = SketchConfig {
            kmer_size: 11,
            fscale: 1,
            num_threads: 2,
            memory: MIN_MEMORY_GB,
            send_timeout: Duration::from_micros(100),
            ..Default::default()
        };

        let result = run(&[input.path().to_path_buf()], &config).unwrap();
        assert!(result.bucket_entry_counts.iter().sum::<u64>() > 0);
    }

    #[test]
    fn test_channel_capacity_calculation() {
        let cap_4gb = compute_channel_capacity(4, 0);
        assert_eq!(cap_4gb, OPTIMAL_CHANNEL_CAPACITY);

        let cap_4gb_with_input = compute_channel_capacity(4, 3 * 1024 * 1024 * 1024);
        assert!(cap_4gb_with_input < cap_4gb);

        let cap_exceeded = compute_channel_capacity(4, 10 * 1024 * 1024 * 1024);
        assert_eq!(cap_exceeded, 1024);

        assert_eq!(compute_channel_capacity(0, 0), 1024);
    }

    #[test]
    fn test_scan_fasta_boundaries() {
        let data = b">seq1\nATCG\n>seq2\nGCTA\n>seq3\nAAAA\n";
        assert_eq!(scan_fasta_boundaries(data), vec![0, 11, 22]);
    }

    #[test]
    fn test_scan_fastq_boundaries() {
        let data = b"@read1\nATCG\n+\nIIII\n@read2\nGCTA\n+\nIIII\n";
        assert_eq!(scan_fastq_boundaries(data), vec![0, 19]);
    }

    #[test]
    fn test_scan_fastq_boundaries_wrapped() {
        let data = b"@read1\nATCG\nGCTA\n+\nIIII\nJJJJ\n@read2\nAAAA\n+\nKKKK\n";
        let bounds = scan_fastq_boundaries(data);
        assert_eq!(bounds[0], 0);
        assert!(bounds.contains(&29));
    }

    #[test]
    fn test_scan_fastq_boundaries_at_in_quality() {
        let data = b"@read1\nATCG\n+\n@@@I\n@read2\nGCTA\n+\nIIII\n";
        let bounds = scan_fastq_boundaries(data);
        assert_eq!(bounds, vec![0, 19]);
    }

    #[test]
    fn test_scan_fastq_boundaries_at_followed_by_at_skipped() {
        let data = b"@read1\nATCG\n+\nIIII\n@@ambiguous\n";
        assert_eq!(scan_fastq_boundaries(data), vec![0]);
    }

    #[test]
    fn test_scan_fastq_boundaries_iupac_codes() {
        let data = b"@read1\nRYSW\n+\nIIII\n@read2\nKMBD\n+\nIIII\n";
        let bounds = scan_fastq_boundaries(data);
        assert_eq!(bounds, vec![0, 19]);
    }

    #[test]
    fn test_is_compressed() {
        assert!(is_compressed([0x1F, 0x8B]));
        assert!(is_compressed([0x42, 0x5A]));
        assert!(is_compressed([0xFD, 0x37]));
        assert!(is_compressed([0x28, 0xB5]));
        assert!(!is_compressed([b'>', b's']));
        assert!(!is_compressed([b'@', b'r']));
    }

    #[test]
    fn test_mmap_slice_reader() {
        use std::io::Read;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.bin");
        std::fs::write(&path, b"Hello, World!").unwrap();

        let file = File::open(&path).unwrap();
        let mmap = Arc::new(unsafe { Mmap::map(&file).unwrap() });

        let mut reader = MmapSliceReader::new(mmap, 7, 12);
        let mut buf = String::new();
        reader.read_to_string(&mut buf).unwrap();
        assert_eq!(buf, "World");
    }

    #[test]
    fn test_tiny_file_errors() {
        let dir = tempfile::tempdir().unwrap();

        let empty_path = dir.path().join("empty.fa");
        std::fs::write(&empty_path, b"").unwrap();

        let config = SketchConfig::default();
        let result = run(&[empty_path], &config);
        let err = match result {
            Err(e) => e,
            Ok(_) => panic!("expected error for empty file"),
        };
        assert!(err.to_string().contains("too small"));

        let tiny_path = dir.path().join("tiny.fa");
        std::fs::write(&tiny_path, b">").unwrap();

        let result = run(&[tiny_path], &config);
        let err = match result {
            Err(e) => e,
            Ok(_) => panic!("expected error for 1-byte file"),
        };
        assert!(err.to_string().contains("too small"));
    }

    #[test]
    fn test_fscale_zero_errors() {
        let input = make_fasta(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
        let config = SketchConfig {
            fscale: 0,
            ..Default::default()
        };

        let result = run(&[input.path().to_path_buf()], &config);
        let err = match result {
            Err(e) => e,
            Ok(_) => panic!("expected error for fscale=0"),
        };
        assert!(err.to_string().contains("fscale must be non-zero"));
    }

    #[test]
    fn test_kmer_size_validation() {
        let input = make_fasta(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);

        let config = SketchConfig {
            kmer_size: 0,
            memory: 1,
            ..Default::default()
        };
        let err = match run(&[input.path().to_path_buf()], &config) {
            Err(e) => e,
            Ok(_) => panic!("expected error for kmer_size=0"),
        };
        assert!(
            err.to_string()
                .contains("kmer_size must be between 1 and 31")
        );

        let config = SketchConfig {
            kmer_size: 32,
            memory: 1,
            ..Default::default()
        };
        let err = match run(&[input.path().to_path_buf()], &config) {
            Err(e) => e,
            Ok(_) => panic!("expected error for kmer_size=32"),
        };
        assert!(
            err.to_string()
                .contains("kmer_size must be between 1 and 31")
        );

        let config = SketchConfig {
            kmer_size: 31,
            fscale: 1,
            memory: 1,
            ..Default::default()
        };
        let result = run(&[input.path().to_path_buf()], &config);
        assert!(result.is_ok());
    }
}
