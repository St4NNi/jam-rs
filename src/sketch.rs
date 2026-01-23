use crate::bias::BiasTable;
use crate::core_utils::passes_entropy_filter;
use crate::format::{BUCKET_COUNT, ENTRY_SIZE, Entry, bucket_id};
use crate::io::EntryWriter;
use crossfire::mpsc;
use crossfire::{MTx, Rx};
use jamhash::jamhash_u64;
use memmap2::Mmap;
use needletail::{Sequence, parse_fastx_reader};
use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicU32, Ordering};
use std::time::Duration;
use tempfile::TempDir;

const WRITE_BUFFER_SIZE: usize = 8 * 1024 * 1024;
const MIN_MEMORY_GB: usize = 4;
const DEFAULT_SEND_TIMEOUT: Duration = Duration::from_millis(1);
const MIN_SPLIT_SIZE: usize = 1024 * 1024;
const MAX_CONCURRENT_MMAPS: usize = 256;

type Sender = MTx<mpsc::Array<Entry>>;
type Receiver = Rx<mpsc::Array<Entry>>;

#[derive(Clone)]
pub struct SketchConfig {
    pub kmer_size: u8,
    pub fscale: u64,
    pub num_threads: usize,
    pub memory: usize, // GB
    pub temp_dir_base: Option<PathBuf>,
    pub min_entropy: f64,
    pub singleton: bool,
    pub bias_table: Option<Arc<BiasTable>>,
    pub send_timeout: Duration,
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
        }
    }
}

pub struct SketchResult {
    pub sample_count: u32,
    pub bucket_entry_counts: [u64; BUCKET_COUNT],
    pub threshold: u64,
    pub temp_dir: TempDir,
}

#[derive(Debug, thiserror::Error)]
pub enum SketchError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Parse error in {path}: {message}")]
    Parse { path: PathBuf, message: String },

    #[error("Channel send error")]
    Channel,
}

struct WorkUnit {
    mmap: Option<Arc<Mmap>>,
    start: usize,
    end: usize,
    sample_id: Option<u32>,
    source_path: PathBuf,
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

// https://github.com/onecodex/needletail/blob/master/src/parser/mod.rs
const GZ_MAGIC: [u8; 2] = [0x1F, 0x8B];
const BZ_MAGIC: [u8; 2] = [0x42, 0x5A];
const XZ_MAGIC: [u8; 2] = [0xFD, 0x37];
const ZST_MAGIC: [u8; 2] = [0x28, 0xB5];

#[inline]
fn is_compressed_magic(magic: [u8; 2]) -> bool {
    matches!(magic, GZ_MAGIC | BZ_MAGIC | XZ_MAGIC | ZST_MAGIC)
}

fn scan_fasta_boundaries(data: &[u8]) -> Vec<usize> {
    let mut boundaries = vec![0];
    for i in 0..data.len() - 1 {
        if data[i] == b'\n' && data[i + 1] == b'>' {
            boundaries.push(i + 1);
        }
    }
    boundaries
}

fn scan_fastq_boundaries(data: &[u8]) -> Vec<usize> {
    let mut boundaries = vec![0];
    let mut line_in_record = 0;
    for (pos, &byte) in data.iter().enumerate() {
        if byte == b'\n' {
            line_in_record += 1;
            if line_in_record == 4 {
                if pos + 1 < data.len() {
                    boundaries.push(pos + 1);
                }
                line_in_record = 0;
            }
        }
    }
    boundaries
}

struct RecordPosition {
    mmap: Arc<Mmap>,
    offset: usize,
    path: PathBuf,
}

fn build_work_units(
    input_files: &[PathBuf],
    num_threads: usize,
    singleton: bool,
    sample_counter: &AtomicU32,
) -> Result<Vec<Vec<WorkUnit>>, SketchError> {
    // Too many files: skip mmapping, distribute round-robin for direct reading
    if input_files.len() > MAX_CONCURRENT_MMAPS {
        let mut thread_work: Vec<Vec<WorkUnit>> = (0..num_threads).map(|_| Vec::new()).collect();
        for (i, path) in input_files.iter().enumerate() {
            let t = i % num_threads;
            let sample_id = if !singleton {
                Some(sample_counter.fetch_add(1, Ordering::SeqCst))
            } else {
                None
            };
            thread_work[t].push(WorkUnit {
                mmap: None,
                start: 0,
                end: 0,
                sample_id,
                source_path: path.clone(),
            });
        }
        return Ok(thread_work);
    }

    let mut record_positions: Vec<RecordPosition> = Vec::new();
    let mut splittable: Vec<(Arc<Mmap>, PathBuf, Vec<usize>)> = Vec::new();
    let mut whole_file: Vec<(Arc<Mmap>, PathBuf)> = Vec::new();

    for path in input_files {
        let file = File::open(path)?;
        let mmap = Arc::new(unsafe { Mmap::map(&file)? });

        if mmap.len() < 2 {
            continue;
        }

        let magic = [mmap[0], mmap[1]];

        if is_compressed_magic(magic) {
            whole_file.push((mmap, path.clone()));
            continue;
        }

        if !matches!(magic[0], b'>' | b'@') {
            return Err(SketchError::Parse {
                path: path.clone(),
                message: format!(
                    "unrecognized file format (first bytes: [{:#04X}, {:#04X}])",
                    magic[0], magic[1]
                ),
            });
        }

        let is_fasta = magic[0] == b'>';
        let boundaries = if mmap.len() < MIN_SPLIT_SIZE {
            vec![0]
        } else if is_fasta {
            scan_fasta_boundaries(&mmap)
        } else {
            scan_fastq_boundaries(&mmap)
        };

        if boundaries.is_empty() {
            continue;
        }
        splittable.push((mmap, path.clone(), boundaries));
    }

    for (mmap, path, boundaries) in &splittable {
        for &offset in boundaries {
            record_positions.push(RecordPosition {
                mmap: Arc::clone(mmap),
                offset,
                path: path.clone(),
            });
        }
    }

    let mut thread_work: Vec<Vec<WorkUnit>> = (0..num_threads).map(|_| Vec::new()).collect();

    if !record_positions.is_empty() {
        let records_per_thread = record_positions.len() / num_threads;
        let extra = record_positions.len() % num_threads;

        let mut rec_idx = 0;
        for t in 0..num_threads {
            let count = records_per_thread + if t < extra { 1 } else { 0 };
            if count == 0 {
                continue;
            }

            let end_idx = rec_idx + count;
            let mut chunk_start = rec_idx;

            while chunk_start < end_idx {
                let chunk_mmap = &record_positions[chunk_start].mmap;
                let chunk_path = &record_positions[chunk_start].path;

                let mut chunk_end = chunk_start + 1;
                while chunk_end < end_idx
                    && Arc::ptr_eq(&record_positions[chunk_end].mmap, chunk_mmap)
                {
                    chunk_end += 1;
                }

                let start_byte = record_positions[chunk_start].offset;
                let end_byte = if chunk_end < record_positions.len()
                    && Arc::ptr_eq(&record_positions[chunk_end].mmap, chunk_mmap)
                {
                    record_positions[chunk_end].offset
                } else {
                    chunk_mmap.len()
                };

                thread_work[t].push(WorkUnit {
                    mmap: Some(Arc::clone(chunk_mmap)),
                    start: start_byte,
                    end: end_byte,
                    sample_id: None,
                    source_path: chunk_path.clone(),
                });

                chunk_start = chunk_end;
            }
            rec_idx = end_idx;
        }
    }

    if !singleton {
        let mut file_ids: std::collections::HashMap<PathBuf, u32> =
            std::collections::HashMap::new();
        for thread_units in &mut thread_work {
            for unit in thread_units.iter_mut() {
                let id = *file_ids
                    .entry(unit.source_path.clone())
                    .or_insert_with(|| sample_counter.fetch_add(1, Ordering::SeqCst));
                unit.sample_id = Some(id);
            }
        }
    }

    // Compressed files: round-robin, needletail handles decompression
    for (i, (mmap, path)) in whole_file.into_iter().enumerate() {
        let t = i % num_threads;
        let sample_id = if !singleton {
            Some(sample_counter.fetch_add(1, Ordering::SeqCst))
        } else {
            None
        };
        let end = mmap.len();
        thread_work[t].push(WorkUnit {
            mmap: Some(mmap),
            start: 0,
            end,
            sample_id,
            source_path: path,
        });
    }

    Ok(thread_work)
}

fn compute_channel_capacity(memory_gb: usize) -> usize {
    let memory_gb = memory_gb.max(MIN_MEMORY_GB);
    let memory_bytes = memory_gb as u64 * 1024 * 1024 * 1024;
    let writer_memory = 256u64 * WRITE_BUFFER_SIZE as u64;
    let channel_memory = memory_bytes.saturating_sub(writer_memory);
    let capacity = channel_memory / (BUCKET_COUNT as u64 * ENTRY_SIZE as u64);
    capacity.max(1024) as usize
}

pub fn run(input_files: &[PathBuf], config: &SketchConfig) -> Result<SketchResult, SketchError> {
    let temp_dir = match &config.temp_dir_base {
        Some(base) => tempfile::Builder::new().prefix("jam_").tempdir_in(base)?,
        None => tempfile::Builder::new().prefix("jam_").tempdir()?,
    };

    let threshold = u64::MAX / config.fscale;
    let sample_counter = AtomicU32::new(0);
    let num_threads = config.num_threads.max(1);

    let thread_work =
        build_work_units(input_files, num_threads, config.singleton, &sample_counter)?;

    let channel_capacity = compute_channel_capacity(config.memory);
    let mut senders: Vec<Sender> = Vec::with_capacity(BUCKET_COUNT);
    let mut receivers: Vec<Receiver> = Vec::with_capacity(BUCKET_COUNT);
    for _ in 0..BUCKET_COUNT {
        let (tx, rx) = mpsc::bounded_blocking(channel_capacity);
        senders.push(tx);
        receivers.push(rx);
    }

    let temp_path = temp_dir.path().to_path_buf();
    let receivers_per_thread = BUCKET_COUNT / num_threads;
    let extra_receivers = BUCKET_COUNT % num_threads;

    let mut receiver_groups: Vec<Vec<Receiver>> = Vec::with_capacity(num_threads);
    let mut writer_groups: Vec<Vec<EntryWriter>> = Vec::with_capacity(num_threads);
    let mut bucket_id_groups: Vec<Vec<usize>> = Vec::with_capacity(num_threads);

    let mut rx_iter = receivers.into_iter().enumerate();
    for t in 0..num_threads {
        let count = receivers_per_thread + if t < extra_receivers { 1 } else { 0 };
        let mut thread_receivers = Vec::with_capacity(count);
        let mut thread_writers = Vec::with_capacity(count);
        let mut thread_bucket_ids = Vec::with_capacity(count);

        for _ in 0..count {
            let (bucket_idx, rx) = rx_iter.next().unwrap();
            thread_receivers.push(rx);
            let path = temp_path.join(format!("bucket_{bucket_idx:03}.bin"));
            thread_writers.push(EntryWriter::new(&path, WRITE_BUFFER_SIZE)?);
            thread_bucket_ids.push(bucket_idx);
        }
        receiver_groups.push(thread_receivers);
        writer_groups.push(thread_writers);
        bucket_id_groups.push(thread_bucket_ids);
    }

    let thread_results: std::sync::Mutex<Vec<Result<Vec<(usize, u64)>, SketchError>>> =
        std::sync::Mutex::new(Vec::with_capacity(num_threads));

    {
        let mut thread_work = thread_work;
        let mut receiver_groups = receiver_groups;
        let mut writer_groups = writer_groups;
        let mut bucket_id_groups = bucket_id_groups;

        let sample_counter_ref = &sample_counter;
        let thread_results_ref = &thread_results;

        rayon::scope(|s| {
            for (((work, receivers), writers), bucket_ids) in thread_work
                .drain(..)
                .zip(receiver_groups.drain(..))
                .zip(writer_groups.drain(..))
                .zip(bucket_id_groups.drain(..))
            {
                let thread_senders: Vec<Sender> = senders.iter().map(|s| s.clone()).collect();

                s.spawn(move |_| {
                    let result = thread_work_fn(
                        work,
                        thread_senders,
                        receivers,
                        writers,
                        bucket_ids,
                        config,
                        sample_counter_ref,
                        threshold,
                    );
                    thread_results_ref.lock().unwrap().push(result);
                });
            }

            drop(senders);
        });
    }

    let mut bucket_entry_counts = [0u64; BUCKET_COUNT];
    for result in thread_results.into_inner().unwrap() {
        for (bucket_idx, count) in result? {
            bucket_entry_counts[bucket_idx] = count;
        }
    }

    let sample_count = sample_counter.load(Ordering::SeqCst);

    Ok(SketchResult {
        sample_count,
        bucket_entry_counts,
        threshold,
        temp_dir,
    })
}

fn thread_work_fn(
    work_units: Vec<WorkUnit>,
    senders: Vec<Sender>,
    receivers: Vec<Receiver>,
    mut writers: Vec<EntryWriter>,
    bucket_ids: Vec<usize>,
    config: &SketchConfig,
    sample_counter: &AtomicU32,
    threshold: u64,
) -> Result<Vec<(usize, u64)>, SketchError> {
    let send_timeout = config.send_timeout;

    for unit in &work_units {
        let fastx = match &unit.mmap {
            Some(mmap) => {
                let reader = MmapSliceReader::new(Arc::clone(mmap), unit.start, unit.end);
                parse_fastx_reader(reader)
            }
            None => {
                let file = File::open(&unit.source_path)?;
                let reader = io::BufReader::new(file);
                parse_fastx_reader(reader)
            }
        };

        let mut fastx_reader = fastx.map_err(|e| SketchError::Parse {
            path: unit.source_path.clone(),
            message: e.to_string(),
        })?;

        sketch_records(
            fastx_reader.as_mut(),
            unit.sample_id,
            &unit.source_path,
            &senders,
            &receivers,
            &mut writers,
            config,
            sample_counter,
            threshold,
            send_timeout,
        )?;
    }

    drop(senders);
    drain_until_disconnected(&receivers, &mut writers);

    for writer in &mut writers {
        writer.flush()?;
    }

    Ok(bucket_ids
        .into_iter()
        .zip(writers.iter())
        .map(|(id, w)| (id, w.count()))
        .collect())
}

#[allow(clippy::too_many_arguments)]
fn sketch_records(
    reader: &mut dyn needletail::FastxReader,
    file_sample_id: Option<u32>,
    source_path: &Path,
    senders: &[Sender],
    receivers: &[Receiver],
    writers: &mut [EntryWriter],
    config: &SketchConfig,
    sample_counter: &AtomicU32,
    threshold: u64,
    send_timeout: Duration,
) -> Result<(), SketchError> {
    let k = config.kmer_size;
    let min_entropy = config.min_entropy;

    while let Some(record) = reader.next() {
        let record = record.map_err(|e| SketchError::Parse {
            path: source_path.to_path_buf(),
            message: e.to_string(),
        })?;

        let sample_id = match file_sample_id {
            Some(id) => id,
            None => sample_counter.fetch_add(1, Ordering::SeqCst),
        };

        let sequence = record.normalize(false);
        if sequence.len() < k as usize {
            continue;
        }

        for (_, kmer, _) in sequence.bit_kmers(k, true) {
            let hash = jamhash_u64(kmer.0);

            if hash >= threshold {
                continue;
            }

            if min_entropy > 0.0 && !passes_entropy_filter(kmer.0, k, min_entropy) {
                continue;
            }

            if let Some(ref bias) = config.bias_table {
                if !bias.passes_filter(kmer.0, k) {
                    continue;
                }
            }

            let entry = Entry::new(hash, sample_id);
            let bucket = bucket_id(hash);

            match senders[bucket].send_timeout(entry, send_timeout) {
                Ok(()) => {}
                Err(crossfire::SendTimeoutError::Timeout(entry)) => {
                    drain_own_receivers(receivers, writers);
                    let _ = senders[bucket].send(entry);
                }
                Err(crossfire::SendTimeoutError::Disconnected(_)) => {
                    break;
                }
            }
        }
    }

    Ok(())
}

fn drain_own_receivers(receivers: &[Receiver], writers: &mut [EntryWriter]) {
    for (rx, writer) in receivers.iter().zip(writers.iter_mut()) {
        while let Ok(entry) = rx.try_recv() {
            let _ = writer.write(&entry);
        }
    }
}

fn drain_until_disconnected(receivers: &[Receiver], writers: &mut [EntryWriter]) {
    loop {
        let mut all_disconnected = true;
        for (rx, writer) in receivers.iter().zip(writers.iter_mut()) {
            loop {
                match rx.try_recv() {
                    Ok(entry) => {
                        let _ = writer.write(&entry);
                        all_disconnected = false;
                    }
                    Err(crossfire::TryRecvError::Empty) => {
                        all_disconnected = false;
                        break;
                    }
                    Err(crossfire::TryRecvError::Disconnected) => {
                        break;
                    }
                }
            }
        }
        if all_disconnected {
            break;
        }
        std::thread::sleep(Duration::from_micros(100));
    }
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
            ..Default::default()
        };

        let result = run(&[input.path().to_path_buf()], &config).unwrap();
        assert_eq!(result.sample_count, 1);

        let total: u64 = result.bucket_entry_counts.iter().sum();
        assert!(total > 0);
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
            ..Default::default()
        };

        let result = run(&[input.path().to_path_buf()], &config).unwrap();
        assert_eq!(result.sample_count, 2);
    }

    #[test]
    fn test_sketch_fracmin_filters() {
        let input = make_fasta(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);

        let config_all = SketchConfig {
            kmer_size: 11,
            fscale: 1,
            ..Default::default()
        };
        let result_all = run(&[input.path().to_path_buf()], &config_all).unwrap();

        let config_low = SketchConfig {
            kmer_size: 11,
            fscale: 100,
            ..Default::default()
        };
        let result_low = run(&[input.path().to_path_buf()], &config_low).unwrap();

        let total_all: u64 = result_all.bucket_entry_counts.iter().sum();
        let total_low: u64 = result_low.bucket_entry_counts.iter().sum();
        assert!(total_low < total_all);
    }

    #[test]
    fn test_sketch_multiple_files() {
        let input1 = make_fasta(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
        let input2 = make_fasta(&[("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")]);

        let config = SketchConfig {
            kmer_size: 11,
            fscale: 1,
            num_threads: 2,
            ..Default::default()
        };

        let result = run(
            &[input1.path().to_path_buf(), input2.path().to_path_buf()],
            &config,
        )
        .unwrap();
        assert_eq!(result.sample_count, 2);

        let total: u64 = result.bucket_entry_counts.iter().sum();
        assert!(total > 0);
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
        let total: u64 = result.bucket_entry_counts.iter().sum();
        assert!(total > 0);
    }

    #[test]
    fn test_channel_capacity_calculation() {
        let cap = compute_channel_capacity(4);
        assert!(cap > 100_000);

        let cap8 = compute_channel_capacity(8);
        assert!(cap8 > cap);

        let cap_low = compute_channel_capacity(1);
        assert_eq!(cap_low, cap);
    }

    #[test]
    fn test_scan_fasta_boundaries() {
        let data = b">seq1\nATCG\n>seq2\nGCTA\n>seq3\nAAAA\n";
        let bounds = scan_fasta_boundaries(data);
        assert_eq!(bounds, vec![0, 11, 22]);
    }

    #[test]
    fn test_scan_fastq_boundaries() {
        let data = b"@read1\nATCG\n+\nIIII\n@read2\nGCTA\n+\nIIII\n";
        let bounds = scan_fastq_boundaries(data);
        assert_eq!(bounds, vec![0, 19]);
    }

    #[test]
    fn test_is_compressed_magic() {
        assert!(is_compressed_magic([0x1F, 0x8B]));
        assert!(is_compressed_magic([0x42, 0x5A]));
        assert!(is_compressed_magic([0xFD, 0x37]));
        assert!(is_compressed_magic([0x28, 0xB5]));
        assert!(!is_compressed_magic([b'>', b's']));
        assert!(!is_compressed_magic([b'@', b'r']));
        assert!(!is_compressed_magic([0x00, 0x00]));
    }

    #[test]
    fn test_mmap_slice_reader() {
        use std::io::Read;
        let data = b"Hello, World!";
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.bin");
        std::fs::write(&path, data).unwrap();

        let file = File::open(&path).unwrap();
        let mmap = Arc::new(unsafe { Mmap::map(&file).unwrap() });

        let mut reader = MmapSliceReader::new(Arc::clone(&mmap), 7, 12);
        let mut buf = String::new();
        reader.read_to_string(&mut buf).unwrap();
        assert_eq!(buf, "World");
    }
}
