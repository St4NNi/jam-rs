use crate::bias::HashBiasTable;
use crate::format::{BUCKET_COUNT, bucket_id};
use crate::reader::{JamReader, ReaderError};
use jamhash::jamhash_u64;
use needletail::{Sequence, parse_fastx_file};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::Arc;

/// Query sketch with sorted (hash, sample_id) pairs across 256 buckets.
#[derive(Debug)]
pub struct QuerySketch {
    pub buckets: [Vec<(u64, u32)>; BUCKET_COUNT],
    pub sample_names: Vec<String>,
    pub query_sizes: Vec<usize>,
}

impl QuerySketch {
    pub fn new() -> Self {
        Self {
            buckets: std::array::from_fn(|_| Vec::new()),
            sample_names: Vec::new(),
            query_sizes: Vec::new(),
        }
    }

    #[inline]
    pub fn bucket(&self, idx: usize) -> &[(u64, u32)] {
        &self.buckets[idx]
    }

    #[inline]
    pub fn sample_count(&self) -> usize {
        self.sample_names.len()
    }

    #[inline]
    pub fn total_entries(&self) -> usize {
        self.buckets.iter().map(|b| b.len()).sum()
    }

    /// Build from .jam file. Fails if parameters don't match target db.
    pub fn from_jam<P: AsRef<Path>>(path: P, db: &JamReader) -> Result<Self, QueryError> {
        let source = JamReader::open(path)?;

        if source.kmer_size() != db.kmer_size() {
            return Err(QueryError::ParameterMismatch {
                parameter: "k-mer size".to_string(),
                source_value: source.kmer_size().to_string(),
                target_value: db.kmer_size().to_string(),
            });
        }

        if source.threshold() != db.threshold() {
            return Err(QueryError::ParameterMismatch {
                parameter: "hash threshold".to_string(),
                source_value: source.threshold().to_string(),
                target_value: db.threshold().to_string(),
            });
        }

        let stats = source.stats();
        let expected_sample_count = stats.sample_count as usize;

        let sample_names = source.sample_names().to_vec();
        if sample_names.len() != expected_sample_count {
            return Err(QueryError::Parse {
                path: "JAM file".to_string(),
                message: format!(
                    "sample names count ({}) doesn't match header sample_count ({})",
                    sample_names.len(),
                    expected_sample_count
                ),
            });
        }

        let stored_sizes = source.sample_sizes();
        let query_sizes: Vec<usize> = stored_sizes.iter().map(|&s| s as usize).collect();

        let mut buckets: [Vec<(u64, u32)>; BUCKET_COUNT] = std::array::from_fn(|_| Vec::new());
        for (bucket_idx, bucket) in buckets.iter_mut().enumerate() {
            for entry in source.bucket_entries(bucket_idx) {
                bucket.push((entry.hash, entry.sample_id));
            }
        }

        Ok(Self {
            buckets,
            sample_names,
            query_sizes,
        })
    }

    pub fn from_fasta<P: AsRef<Path>>(
        input: P,
        db: &JamReader,
        singleton: bool,
    ) -> Result<Self, QueryError> {
        let input_path = input.as_ref();
        let kmer_size = db.kmer_size();
        let threshold = db.threshold();
        let bias_table = db.bias_table();

        let mut reader = match parse_fastx_file(input_path) {
            Ok(reader) => reader,
            Err(e) if e.kind == needletail::errors::ParseErrorKind::EmptyFile => {
                eprintln!(
                    "Empty file detected: {}, returning empty sketch",
                    input_path.display()
                );
                return Ok(Self::new());
            }
            Err(e) => {
                return Err(QueryError::Parse {
                    path: input_path.display().to_string(),
                    message: e.to_string(),
                });
            }
        };

        let mut buckets: [Vec<(u64, u32)>; BUCKET_COUNT] = std::array::from_fn(|_| Vec::new());
        let mut sample_names: Vec<String> = Vec::new();
        let mut sample_hash_sets: Vec<HashSet<u64>> = Vec::new();
        let mut current_sample_id: u32 = 0;

        if !singleton {
            sample_names.push(
                input_path
                    .file_name()
                    .and_then(|s| s.to_str())
                    .unwrap_or("query")
                    .to_string(),
            );
            sample_hash_sets.push(HashSet::new());
        }

        while let Some(record) = reader.next() {
            let record = record.map_err(|e| QueryError::Parse {
                path: input_path.display().to_string(),
                message: e.to_string(),
            })?;

            if singleton {
                let name = std::str::from_utf8(record.id())
                    .unwrap_or("unknown")
                    .to_string();
                sample_names.push(name);
                sample_hash_sets.push(HashSet::new());
                current_sample_id = (sample_names.len() - 1) as u32;
            }

            let sequence = record.normalize(false);
            if sequence.len() < kmer_size as usize {
                continue;
            }

            for (_, kmer, _) in sequence.bit_kmers(kmer_size, true) {
                let hash = jamhash_u64(kmer.0);

                if hash >= threshold {
                    continue;
                }

                if bias_table.as_ref().is_some_and(|b| !b.passes_filter(hash)) {
                    continue;
                }

                if sample_hash_sets[current_sample_id as usize].insert(hash) {
                    buckets[bucket_id(hash)].push((hash, current_sample_id));
                }
            }
        }

        for bucket in &mut buckets {
            bucket.sort_unstable();
            bucket.dedup();
        }

        let query_sizes: Vec<usize> = sample_hash_sets.iter().map(|set| set.len()).collect();

        Ok(Self {
            buckets,
            sample_names,
            query_sizes,
        })
    }

    pub fn from_inputs(
        inputs: &[std::path::PathBuf],
        db: &JamReader,
        singleton: bool,
    ) -> Result<Self, QueryError> {
        use crate::format::MAGIC;
        use std::fs::File;
        use std::io::Read;

        if inputs.is_empty() {
            return Ok(Self::new());
        }

        let is_jam_file = |path: &std::path::PathBuf| -> bool {
            if path
                .extension()
                .is_some_and(|ext| ext.eq_ignore_ascii_case("jam"))
            {
                return true;
            }
            File::open(path)
                .ok()
                .and_then(|mut f| {
                    let mut magic = [0u8; 4];
                    f.read_exact(&mut magic).ok()?;
                    Some(magic == MAGIC)
                })
                .unwrap_or(false)
        };

        let mut combined = Self::new();

        for input in inputs {
            let sketch = if is_jam_file(input) {
                Self::from_jam(input, db)?
            } else {
                Self::from_fasta(input, db, singleton)?
            };

            let sample_offset = combined.sample_count() as u32;
            combined.sample_names.extend(sketch.sample_names);
            combined.query_sizes.extend(sketch.query_sizes);

            for (bucket_idx, bucket) in sketch.buckets.into_iter().enumerate() {
                for (hash, sample_id) in bucket {
                    combined.buckets[bucket_idx].push((hash, sample_id + sample_offset));
                }
            }
        }

        for bucket in &mut combined.buckets {
            bucket.sort_unstable();
        }

        Ok(combined)
    }
}

impl Default for QuerySketch {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Debug, thiserror::Error)]
pub enum QueryError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Database error: {0}")]
    Database(#[from] ReaderError),

    #[error("Parse error in {path}: {message}")]
    Parse { path: String, message: String },

    #[error(
        "Parameter mismatch: {parameter} - source has {source_value}, target database has {target_value}"
    )]
    ParameterMismatch {
        parameter: String,
        source_value: String,
        target_value: String,
    },
}

#[derive(Debug, Clone)]
pub struct SampleMatch {
    pub sample_id: u32,
    pub hit_count: u32,
    pub containment: f64,
}

#[derive(Debug, Clone)]
pub struct QueryResult {
    pub query_size: usize,
    pub hashes_found: usize,
    pub matches: Vec<SampleMatch>,
    /// Number of buckets that failed to mmap during query.
    /// If > 0, results may be incomplete (some hashes could not be searched).
    pub failed_bucket_count: usize,
}

impl QueryResult {
    pub fn top(&self, n: usize) -> Vec<&SampleMatch> {
        let mut sorted: Vec<_> = self.matches.iter().collect();
        sorted.sort_by(|a, b| b.containment.total_cmp(&a.containment));
        sorted.truncate(n);
        sorted
    }

    pub fn above_threshold(&self, min_containment: f64) -> Vec<&SampleMatch> {
        self.matches
            .iter()
            .filter(|m| m.containment >= min_containment)
            .collect()
    }

    pub fn has_matches(&self) -> bool {
        !self.matches.is_empty()
    }

    /// Returns true if some buckets failed during query, meaning results may be incomplete.
    pub fn is_partial(&self) -> bool {
        self.failed_bucket_count > 0
    }
}

pub struct QueryEngine {
    reader: JamReader,
    bias_table: Option<Arc<HashBiasTable>>,
}

impl QueryEngine {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self, ReaderError> {
        let reader = JamReader::open(path)?;
        let bias_table = reader.bias_table();
        Ok(Self { reader, bias_table })
    }

    pub fn threshold(&self) -> u64 {
        self.reader.threshold()
    }

    pub fn kmer_size(&self) -> u8 {
        self.reader.kmer_size()
    }

    pub fn bias_table(&self) -> Option<Arc<HashBiasTable>> {
        self.bias_table.clone()
    }

    pub fn has_bias_table(&self) -> bool {
        self.bias_table.is_some()
    }

    pub fn reader(&self) -> &JamReader {
        &self.reader
    }

    pub fn query(&self, hashes: &[u64]) -> QueryResult {
        if hashes.is_empty() {
            return QueryResult {
                query_size: 0,
                hashes_found: 0,
                matches: Vec::new(),
                failed_bucket_count: 0,
            };
        }

        let mut sorted_hashes = hashes.to_vec();
        sorted_hashes.sort_unstable_by_key(|&h| (h & 0xFF, h));

        let mut sample_hits: HashMap<u32, u32> = HashMap::new();
        let mut hashes_found = 0;

        for &hash in &sorted_hashes {
            let mut found = false;
            for sample_id in self.reader.search(hash) {
                *sample_hits.entry(sample_id).or_insert(0) += 1;
                found = true;
            }
            if found {
                hashes_found += 1;
            }
        }

        let query_size = hashes.len();
        let matches: Vec<SampleMatch> = sample_hits
            .into_iter()
            .map(|(sample_id, hit_count)| SampleMatch {
                sample_id,
                hit_count,
                containment: hit_count as f64 / query_size as f64,
            })
            .collect();

        QueryResult {
            query_size,
            hashes_found,
            matches,
            failed_bucket_count: 0,
        }
    }

    pub fn query_filtered(
        &self,
        hashes: &[u64],
        min_containment: f64,
        max_results: usize,
    ) -> QueryResult {
        let mut result = self.query(hashes);
        result.matches.retain(|m| m.containment >= min_containment);
        result
            .matches
            .sort_by(|a, b| b.containment.total_cmp(&a.containment));
        result.matches.truncate(max_results);
        result
    }

    pub fn query_batch(&self, queries: &[Vec<u64>]) -> Vec<QueryResult> {
        use rayon::prelude::*;
        queries.par_iter().map(|q| self.query(q)).collect()
    }

    pub fn query_sketch(&self, sketch: &QuerySketch) -> Vec<QueryResult> {
        use crate::format::{ENTRY_SIZE, PAGE_SIZE};
        use rayon::prelude::*;
        use std::sync::atomic::{AtomicU32, Ordering};

        let num_samples = sketch.sample_count();
        if num_samples == 0 {
            return Vec::new();
        }

        let threshold = self.reader.threshold();

        // Disable read-ahead since our access is query-driven, not sequential
        self.reader.advise_random();

        let hashes_found: Vec<AtomicU32> = (0..num_samples)
            .into_par_iter()
            .map(|_| AtomicU32::new(0))
            .collect();

        let bucket_pairs: Vec<Vec<(u32, u32)>> = (0..BUCKET_COUNT)
            .into_par_iter()
            .map(|bucket_idx| {
                let mut pairs = Vec::new();
                let query_bucket = sketch.bucket(bucket_idx);
                if query_bucket.is_empty() {
                    return pairs;
                }

                let filter = match self.reader.bucket_filter(bucket_idx) {
                    Some(f) => f,
                    None => return pairs,
                };

                let mut survivors = Vec::with_capacity(query_bucket.len() / 10);
                let mut prev_hash = u64::MAX;
                let mut prev_passed = false;

                for &(hash, sample_id) in query_bucket {
                    if hash != prev_hash {
                        prev_hash = hash;
                        prev_passed = filter.contains(&hash);
                    }
                    if prev_passed {
                        survivors.push((hash, sample_id));
                    }
                }

                // Done with filter - release those pages
                let (filter_start, filter_end) = self.reader.bucket_filter_byte_range(bucket_idx);
                self.reader.release_pages(filter_start, filter_end);

                if survivors.is_empty() {
                    return pairs;
                }

                let db_bucket = self.reader.bucket_entries(bucket_idx);
                let count = db_bucket.len();
                if count == 0 {
                    return pairs;
                }

                let (entry_start, _entry_end) = self.reader.bucket_entry_byte_range(bucket_idx);
                let mut last_released_page = entry_start & !(PAGE_SIZE - 1);

                let mut q_idx = 0;
                while q_idx < survivors.len() {
                    let q_hash = survivors[q_idx].0;

                    let est = ((q_hash as u128 * count as u128) / threshold as u128) as usize;
                    let mut d_idx = est.saturating_sub(16).min(count.saturating_sub(1));

                    while d_idx > 0 && db_bucket[d_idx].hash > q_hash {
                        d_idx -= 1;
                    }

                    while d_idx < count && db_bucket[d_idx].hash < q_hash {
                        d_idx += 1;
                    }

                    while d_idx > 0 && db_bucket[d_idx - 1].hash == q_hash {
                        d_idx -= 1;
                    }

                    // Progressive page release: release pages before current position
                    let current_byte = entry_start + d_idx * ENTRY_SIZE;
                    let current_page = current_byte & !(PAGE_SIZE - 1);
                    if current_page > last_released_page + PAGE_SIZE {
                        // Release pages we've passed (keep one page buffer for backtracking)
                        self.reader.release_pages(last_released_page, current_page - PAGE_SIZE);
                        last_released_page = current_page - PAGE_SIZE;
                    }

                    let db_start = d_idx;
                    let mut db_end = d_idx;
                    while db_end < count && db_bucket[db_end].hash == q_hash {
                        db_end += 1;
                    }
                    let has_matches = db_start < db_end;

                    let mut prev_sample = u32::MAX;
                    while q_idx < survivors.len() && survivors[q_idx].0 == q_hash {
                        let q_sample = survivors[q_idx].1;

                        if q_sample != prev_sample {
                            if has_matches {
                                for db_entry in &db_bucket[db_start..db_end] {
                                    pairs.push((q_sample, db_entry.sample_id));
                                }
                                hashes_found[q_sample as usize].fetch_add(1, Ordering::Relaxed);
                            }
                            prev_sample = q_sample;
                        }
                        q_idx += 1;
                    }
                }

                self.reader.release_bucket(bucket_idx);

                pairs
            })
            .collect();

        // Calculate offsets for parallel concatenation (prefix sum)
        let bucket_sizes: Vec<usize> = bucket_pairs.iter().map(|v| v.len()).collect();
        let total_pairs: usize = bucket_sizes.iter().sum();
        let mut bucket_offsets = Vec::with_capacity(BUCKET_COUNT + 1);
        bucket_offsets.push(0usize);
        for size in &bucket_sizes {
            bucket_offsets.push(bucket_offsets.last().unwrap() + size);
        }

        let mut all_pairs: Vec<(u32, u32)> = vec![(0, 0); total_pairs];
        bucket_pairs
            .into_par_iter()
            .enumerate()
            .for_each(|(bucket_idx, pairs)| {
                let start = bucket_offsets[bucket_idx];
                // Safety: each bucket writes to non-overlapping slice
                let dest = unsafe {
                    std::slice::from_raw_parts_mut(
                        all_pairs.as_ptr().add(start) as *mut (u32, u32),
                        pairs.len(),
                    )
                };
                dest.copy_from_slice(&pairs);
            });

        let merged_hashes_found: Vec<u32> = hashes_found
            .into_par_iter()
            .map(|a| a.load(Ordering::Relaxed))
            .collect();

        all_pairs.par_sort_unstable();

        if all_pairs.is_empty() {
            return (0..num_samples)
                .map(|i| QueryResult {
                    query_size: sketch.query_sizes[i],
                    hashes_found: merged_hashes_found[i] as usize,
                    matches: Vec::new(),
                    failed_bucket_count: 0,
                })
                .collect();
        }

        let sample_starts: Vec<usize> = (0..num_samples as u32)
            .into_par_iter()
            .map(|q_sample| {
                all_pairs
                    .partition_point(|&(qs, _)| qs < q_sample)
            })
            .collect();

        // Process each query sample's range in parallel
        let results: Vec<QueryResult> = (0..num_samples)
            .into_par_iter()
            .map(|sample_idx| {
                let q_sample = sample_idx as u32;
                let start = sample_starts[sample_idx];
                let end = if sample_idx + 1 < num_samples {
                    sample_starts[sample_idx + 1]
                } else {
                    all_pairs.len()
                };

                let mut matches = Vec::new();
                let query_size = sketch.query_sizes[sample_idx];

                // Count consecutive (q_sample, db_sample) pairs within this sample's range
                let mut i = start;
                while i < end {
                    let (_, db_sample) = all_pairs[i];
                    let mut count = 1u32;
                    while i + (count as usize) < end
                        && all_pairs[i + count as usize] == (q_sample, db_sample)
                    {
                        count += 1;
                    }
                    matches.push(SampleMatch {
                        sample_id: db_sample,
                        hit_count: count,
                        containment: if query_size > 0 {
                            count as f64 / query_size as f64
                        } else {
                            0.0
                        },
                    });
                    i += count as usize;
                }

                QueryResult {
                    query_size,
                    hashes_found: merged_hashes_found[sample_idx] as usize,
                    matches,
                    failed_bucket_count: 0,
                }
            })
            .collect();

        results
    }

    pub fn query_fasta<P: AsRef<Path>>(
        &self,
        input: P,
        singleton: bool,
    ) -> Result<Vec<QueryResult>, QueryError> {
        let sketch = QuerySketch::from_fasta(input, &self.reader, singleton)?;
        Ok(self.query_sketch(&sketch))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::writer::{BuildConfig, build};
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

    fn build_test_db(
        seqs: &[(&str, &str)],
        singleton: bool,
    ) -> (tempfile::TempDir, std::path::PathBuf) {
        let input = make_fasta(seqs);
        let output_dir = tempfile::tempdir().unwrap();
        let output_path = output_dir.path().join("test.jam");

        let config = BuildConfig {
            kmer_size: 11,
            fscale: 1,
            singleton,
            num_threads: 1,
            memory: 1,
            ..Default::default()
        };

        build(&[input.path().to_path_buf()], &output_path, &config).unwrap();
        (output_dir, output_path)
    }

    #[test]
    fn test_query_engine_open() {
        let (_dir, path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let engine = QueryEngine::open(&path).unwrap();
        assert!(engine.threshold() > 0);
        assert_eq!(engine.kmer_size(), 11);
    }

    #[test]
    fn test_query_basic() {
        let (_dir, path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let engine = QueryEngine::open(&path).unwrap();

        // Get some hashes from the database
        let reader = JamReader::open(&path).unwrap();
        let mut test_hashes = Vec::new();
        for bucket_idx in 0..256 {
            let entries = reader.bucket_entries(bucket_idx);
            for entry in entries.iter().take(5) {
                test_hashes.push(entry.hash);
            }
            if test_hashes.len() >= 10 {
                break;
            }
        }

        if !test_hashes.is_empty() {
            let result = engine.query(&test_hashes);
            assert!(result.has_matches());
            assert!(result.hashes_found > 0);
        }
    }

    #[test]
    fn test_query_empty() {
        let (_dir, path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let engine = QueryEngine::open(&path).unwrap();

        let result = engine.query(&[]);
        assert!(!result.has_matches());
        assert_eq!(result.query_size, 0);
    }

    #[test]
    fn test_query_nonexistent() {
        let (_dir, path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let engine = QueryEngine::open(&path).unwrap();

        // Query with hashes above threshold (won't exist)
        let fake_hashes: Vec<u64> = (0..10).map(|i| u64::MAX - i).collect();
        let result = engine.query(&fake_hashes);
        assert_eq!(result.hashes_found, 0);
    }

    #[test]
    fn test_query_filtered() {
        let (_dir, path) = build_test_db(
            &[
                ("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
                ("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
            ],
            true,
        );
        let engine = QueryEngine::open(&path).unwrap();

        // Get hashes from sample 0
        let reader = JamReader::open(&path).unwrap();
        let mut test_hashes = Vec::new();
        for bucket_idx in 0..256 {
            for entry in reader.bucket_entries(bucket_idx) {
                if entry.sample_id == 0 {
                    test_hashes.push(entry.hash);
                }
                if test_hashes.len() >= 20 {
                    break;
                }
            }
            if test_hashes.len() >= 20 {
                break;
            }
        }

        if !test_hashes.is_empty() {
            let result = engine.query_filtered(&test_hashes, 0.5, 10);
            assert!(result.matches.len() <= 10);
            for m in &result.matches {
                assert!(m.containment >= 0.5);
            }
        }
    }

    #[test]
    fn test_query_result_helpers() {
        let result = QueryResult {
            query_size: 100,
            hashes_found: 50,
            matches: vec![
                SampleMatch {
                    sample_id: 0,
                    hit_count: 50,
                    containment: 0.5,
                },
                SampleMatch {
                    sample_id: 1,
                    hit_count: 30,
                    containment: 0.3,
                },
                SampleMatch {
                    sample_id: 2,
                    hit_count: 80,
                    containment: 0.8,
                },
            ],
            failed_bucket_count: 0,
        };

        let top2 = result.top(2);
        assert_eq!(top2.len(), 2);
        assert_eq!(top2[0].sample_id, 2);
        assert_eq!(top2[1].sample_id, 0);

        let above_threshold = result.above_threshold(0.4);
        assert_eq!(above_threshold.len(), 2);

        assert!(result.has_matches());
        assert!(!result.is_partial());
    }

    // QuerySketch tests

    #[test]
    fn test_query_sketch_new() {
        let sketch = QuerySketch::new();

        assert_eq!(sketch.sample_count(), 0);
        assert_eq!(sketch.total_entries(), 0);
        assert_eq!(sketch.buckets.len(), 256);
        assert!(sketch.sample_names.is_empty());
        assert!(sketch.query_sizes.is_empty());
    }

    #[test]
    fn test_query_sketch_default() {
        let sketch = QuerySketch::default();

        assert_eq!(sketch.sample_count(), 0);
        assert_eq!(sketch.total_entries(), 0);
    }

    #[test]
    fn test_query_sketch_bucket_accessor() {
        let mut sketch = QuerySketch::new();

        // Add some entries to bucket 0
        sketch.buckets[0].push((100, 0));
        sketch.buckets[0].push((200, 1));

        // Add entries to bucket 255
        sketch.buckets[255].push((300, 0));

        let bucket_0 = sketch.bucket(0);
        assert_eq!(bucket_0.len(), 2);
        assert_eq!(bucket_0[0], (100, 0));
        assert_eq!(bucket_0[1], (200, 1));

        let bucket_255 = sketch.bucket(255);
        assert_eq!(bucket_255.len(), 1);
        assert_eq!(bucket_255[0], (300, 0));

        // Empty bucket
        let bucket_1 = sketch.bucket(1);
        assert!(bucket_1.is_empty());
    }

    #[test]
    fn test_query_sketch_sample_count() {
        let mut sketch = QuerySketch::new();
        assert_eq!(sketch.sample_count(), 0);

        sketch.sample_names.push("sample1".to_string());
        assert_eq!(sketch.sample_count(), 1);

        sketch.sample_names.push("sample2".to_string());
        sketch.sample_names.push("sample3".to_string());
        assert_eq!(sketch.sample_count(), 3);
    }

    #[test]
    fn test_query_sketch_total_entries() {
        let mut sketch = QuerySketch::new();
        assert_eq!(sketch.total_entries(), 0);

        // Add entries to various buckets
        sketch.buckets[0].push((100, 0));
        sketch.buckets[0].push((200, 0));
        assert_eq!(sketch.total_entries(), 2);

        sketch.buckets[50].push((300, 1));
        assert_eq!(sketch.total_entries(), 3);

        sketch.buckets[255].push((400, 0));
        sketch.buckets[255].push((500, 1));
        sketch.buckets[255].push((600, 2));
        assert_eq!(sketch.total_entries(), 6);
    }

    #[test]
    fn test_query_sketch_with_populated_fields() {
        let mut sketch = QuerySketch::new();

        // Set up sample names
        sketch.sample_names = vec!["query_sample_1".to_string(), "query_sample_2".to_string()];

        // Set up query sizes (unique hashes per sample)
        sketch.query_sizes = vec![1000, 500];

        // Populate some buckets
        for i in 0..10 {
            sketch.buckets[i].push((i as u64 * 100, 0));
            sketch.buckets[i].push((i as u64 * 100 + 1, 1));
        }

        assert_eq!(sketch.sample_count(), 2);
        assert_eq!(sketch.total_entries(), 20);
        assert_eq!(sketch.query_sizes[0], 1000);
        assert_eq!(sketch.query_sizes[1], 500);
        assert_eq!(sketch.sample_names[0], "query_sample_1");
    }

    #[test]
    #[should_panic]
    fn test_query_sketch_bucket_out_of_bounds() {
        let sketch = QuerySketch::new();
        let _ = sketch.bucket(256); // Should panic
    }

    // query_sketch tests

    #[test]
    fn test_query_sketch_empty() {
        let (_dir, path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let engine = QueryEngine::open(&path).unwrap();

        let sketch = QuerySketch::new();
        let results = engine.query_sketch(&sketch);
        assert!(results.is_empty());
    }

    #[test]
    fn test_query_sketch_single_sample() {
        // Build a database with one sample
        let (_dir, path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let engine = QueryEngine::open(&path).unwrap();
        let reader = JamReader::open(&path).unwrap();

        // Create a QuerySketch from the database entries (deduplicated)
        let mut sketch = QuerySketch::new();
        sketch.sample_names.push("query_sample".to_string());

        let mut unique_hashes = std::collections::HashSet::new();
        for bucket_idx in 0..256 {
            for entry in reader.bucket_entries(bucket_idx) {
                // Only add each unique hash once
                if unique_hashes.insert(entry.hash) {
                    sketch.buckets[bucket_idx].push((entry.hash, 0)); // sample_id 0
                }
            }
        }
        sketch.query_sizes.push(unique_hashes.len());

        let results = engine.query_sketch(&sketch);

        assert_eq!(results.len(), 1);
        assert!(results[0].has_matches());

        // Verify we got a match with db sample 0
        let db_sample_0_match = results[0].matches.iter().find(|m| m.sample_id == 0);
        assert!(db_sample_0_match.is_some(), "Should match db sample 0");

        // hit_count reflects total match count (may exceed query_size if DB has duplicates)
        // For containment correctness, verify hit_count >= query_size (all query hashes found)
        let m = db_sample_0_match.unwrap();
        assert!(
            m.hit_count >= results[0].query_size as u32,
            "Expected hit_count >= query_size, got {} vs {}",
            m.hit_count,
            results[0].query_size
        );
    }

    #[test]
    fn test_query_sketch_multiple_samples() {
        // Build a database with two different samples
        let (_dir, path) = build_test_db(
            &[
                ("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
                ("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
            ],
            true, // singleton mode - each sequence is a separate sample
        );
        let engine = QueryEngine::open(&path).unwrap();
        let reader = JamReader::open(&path).unwrap();

        // Create a QuerySketch with two query samples, each using hashes from one DB sample
        let mut sketch = QuerySketch::new();
        sketch.sample_names.push("query_0".to_string());
        sketch.sample_names.push("query_1".to_string());

        let mut hashes_per_sample: [std::collections::HashSet<u64>; 2] = Default::default();

        for bucket_idx in 0..256 {
            for entry in reader.bucket_entries(bucket_idx) {
                let query_sample_id = entry.sample_id; // Map db sample to query sample
                if (query_sample_id as usize) < 2 {
                    hashes_per_sample[query_sample_id as usize].insert(entry.hash);
                    sketch.buckets[bucket_idx].push((entry.hash, query_sample_id));
                }
            }
        }

        sketch.query_sizes.push(hashes_per_sample[0].len());
        sketch.query_sizes.push(hashes_per_sample[1].len());

        let results = engine.query_sketch(&sketch);

        assert_eq!(results.len(), 2);

        // Each query sample should have high containment with its corresponding DB sample
        for (query_idx, result) in results.iter().enumerate() {
            assert!(result.has_matches());
            // Find match with corresponding DB sample
            let self_match = result
                .matches
                .iter()
                .find(|m| m.sample_id == query_idx as u32);
            if let Some(m) = self_match {
                // Should have high containment with self
                assert!(
                    m.containment >= 0.9,
                    "Query {} should have high containment with DB sample {}, got {}",
                    query_idx,
                    query_idx,
                    m.containment
                );
            }
        }
    }

    #[test]
    fn test_query_sketch_no_matches() {
        let (_dir, path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let engine = QueryEngine::open(&path).unwrap();

        // Create a QuerySketch with hashes that won't be in the database
        let mut sketch = QuerySketch::new();
        sketch.sample_names.push("fake_sample".to_string());
        sketch.query_sizes.push(10);

        // Add fake hashes above threshold (won't exist in DB)
        for i in 0..10 {
            let fake_hash = u64::MAX - i;
            let bucket_idx = (fake_hash & 0xFF) as usize;
            sketch.buckets[bucket_idx].push((fake_hash, 0));
        }

        let results = engine.query_sketch(&sketch);

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].hashes_found, 0);
        assert!(results[0].matches.is_empty());
    }

    #[test]
    fn test_query_sketch_containment_calculation() {
        let (_dir, path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let engine = QueryEngine::open(&path).unwrap();
        let reader = JamReader::open(&path).unwrap();

        // Get half the hashes from the database
        let mut sketch = QuerySketch::new();
        sketch.sample_names.push("half_sample".to_string());

        let mut all_hashes = Vec::new();
        for bucket_idx in 0..256 {
            for entry in reader.bucket_entries(bucket_idx) {
                all_hashes.push((entry.hash, bucket_idx));
            }
        }

        // Take every other hash
        let selected_hashes: Vec<_> = all_hashes.iter().step_by(2).collect();
        sketch.query_sizes.push(selected_hashes.len());

        for &(hash, bucket_idx) in &selected_hashes {
            sketch.buckets[*bucket_idx].push((*hash, 0));
        }

        let results = engine.query_sketch(&sketch);

        assert_eq!(results.len(), 1);
        assert!(results[0].has_matches());
        // With half the hashes, we should get ~100% containment (all queried hashes exist)
        // but hit count should be about half
        let top = results[0].top(1);
        assert!(!top.is_empty());
        // All queried hashes should be found, so containment should be ~1.0
        assert!(top[0].containment >= 0.9);
    }

    // QuerySketch::from_fasta tests

    #[test]
    fn test_from_fasta_non_singleton() {
        // Build a database
        let (_dir, db_path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();

        // Create a FASTA file with the same sequence
        let query_fasta = make_fasta(&[("query_seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);

        // Build QuerySketch from FASTA (non-singleton mode)
        let sketch = QuerySketch::from_fasta(query_fasta.path(), &db, false).unwrap();

        // Should have one sample (entire file)
        assert_eq!(sketch.sample_count(), 1);
        assert!(!sketch.sample_names[0].is_empty());

        // Should have hashes
        assert!(sketch.total_entries() > 0);
        assert!(sketch.query_sizes[0] > 0);

        // query_sizes should match total_entries for single sample
        assert_eq!(sketch.query_sizes[0], sketch.total_entries());
    }

    #[test]
    fn test_from_fasta_singleton() {
        // Build a database
        let (_dir, db_path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();

        // Create a FASTA file with two sequences
        let query_fasta = make_fasta(&[
            ("query_seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
            ("query_seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
        ]);

        // Build QuerySketch from FASTA (singleton mode)
        let sketch = QuerySketch::from_fasta(query_fasta.path(), &db, true).unwrap();

        // Should have two samples (one per sequence)
        assert_eq!(sketch.sample_count(), 2);
        assert_eq!(sketch.sample_names[0], "query_seq1");
        assert_eq!(sketch.sample_names[1], "query_seq2");

        // Each sample should have hashes
        assert!(sketch.query_sizes[0] > 0);
        assert!(sketch.query_sizes[1] > 0);

        // query_sizes should sum to approximately total_entries
        // (may be less due to shared hashes between samples)
        let total_unique: usize = sketch.query_sizes.iter().sum();
        assert!(total_unique <= sketch.total_entries() + sketch.sample_count());
    }

    #[test]
    fn test_from_fasta_uses_db_parameters() {
        // Build a database with specific k-mer size and scale
        let input = make_fasta(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
        let output_dir = tempfile::tempdir().unwrap();
        let db_path = output_dir.path().join("test.jam");

        let config = BuildConfig {
            kmer_size: 15, // Different k-mer size
            fscale: 10,    // Scale of 10 means ~10% of hashes kept
            singleton: false,
            num_threads: 1,
            memory: 1,
            ..Default::default()
        };

        build(&[input.path().to_path_buf()], &db_path, &config).unwrap();
        let db = JamReader::open(&db_path).unwrap();

        // Verify DB parameters
        assert_eq!(db.kmer_size(), 15);

        // Create a query FASTA
        let query_fasta = make_fasta(&[("query", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);

        // Build QuerySketch - should use DB's parameters
        let sketch = QuerySketch::from_fasta(query_fasta.path(), &db, false).unwrap();

        // Should have hashes (using DB's threshold, so filtered)
        assert!(sketch.sample_count() == 1);

        // The sketch should have fewer hashes than with fscale=1
        // because it uses the DB's threshold
    }

    #[test]
    fn test_from_fasta_deduplication() {
        // Build a database
        let (_dir, db_path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();

        // Create a FASTA file with repeated sequence (should produce duplicate k-mers)
        let query_fasta = make_fasta(&[(
            "query",
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        )]);

        let sketch = QuerySketch::from_fasta(query_fasta.path(), &db, false).unwrap();

        // query_sizes should be unique hashes (deduplicated)
        // total_entries should also be deduplicated (one entry per unique hash)
        assert_eq!(sketch.query_sizes[0], sketch.total_entries());
    }

    #[test]
    fn test_from_fasta_bucketization() {
        // Build a database
        let (_dir, db_path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();

        let query_fasta = make_fasta(&[(
            "query",
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        )]);

        let sketch = QuerySketch::from_fasta(query_fasta.path(), &db, false).unwrap();

        // Verify bucketization: each entry's hash should be in the correct bucket
        for (bucket_idx, bucket) in sketch.buckets.iter().enumerate() {
            for &(hash, _sample_id) in bucket {
                assert_eq!(
                    bucket_id(hash),
                    bucket_idx,
                    "Hash {} should be in bucket {}, not {}",
                    hash,
                    bucket_id(hash),
                    bucket_idx
                );
            }
        }
    }

    #[test]
    fn test_from_fasta_sorted_buckets() {
        // Build a database
        let (_dir, db_path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();

        let query_fasta = make_fasta(&[
            (
                "query1",
                "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
            ),
            (
                "query2",
                "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA",
            ),
        ]);

        let sketch = QuerySketch::from_fasta(query_fasta.path(), &db, true).unwrap();

        // Verify each bucket is sorted by (hash, sample_id)
        for bucket in &sketch.buckets {
            for window in bucket.windows(2) {
                assert!(
                    window[0] <= window[1],
                    "Bucket not sorted: {:?} > {:?}",
                    window[0],
                    window[1]
                );
            }
        }
    }

    #[test]
    fn test_from_fasta_short_sequences_skipped() {
        // Build a database with k=11
        let (_dir, db_path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();
        assert_eq!(db.kmer_size(), 11);

        // Create a FASTA with sequences shorter than k
        let query_fasta = make_fasta(&[
            ("short", "ATCGATCG"),                        // 8 bp, < 11
            ("long", "ATCGATCGATCGATCGATCGATCGATCGATCG"), // 32 bp, > 11
        ]);

        let sketch = QuerySketch::from_fasta(query_fasta.path(), &db, true).unwrap();

        // Should have two sample names (both sequences processed)
        assert_eq!(sketch.sample_count(), 2);

        // First sample (short) should have 0 hashes
        assert_eq!(sketch.query_sizes[0], 0);

        // Second sample (long) should have hashes
        assert!(sketch.query_sizes[1] > 0);
    }

    #[test]
    fn test_from_fasta_file_not_found() {
        let (_dir, db_path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();

        let result = QuerySketch::from_fasta("/nonexistent/path.fasta", &db, false);
        assert!(result.is_err());

        if let Err(QueryError::Parse { path, message: _ }) = result {
            assert!(path.contains("nonexistent"));
        } else {
            panic!("Expected Parse error");
        }
    }

    #[test]
    fn test_from_fasta_integration_with_query_engine() {
        // Build a database
        let (_dir, db_path) =
            build_test_db(&[("db_seq", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();
        let engine = QueryEngine::open(&db_path).unwrap();

        // Create a query FASTA with the same sequence
        let query_fasta = make_fasta(&[("query_seq", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);

        // Build QuerySketch from FASTA
        let sketch = QuerySketch::from_fasta(query_fasta.path(), &db, false).unwrap();

        // Query using the engine
        let results = engine.query_sketch(&sketch);

        assert_eq!(results.len(), 1);
        assert!(results[0].has_matches());

        // Should have high containment (same sequence)
        let top = results[0].top(1);
        assert!(!top.is_empty());
        assert!(
            top[0].containment >= 0.9,
            "Expected high containment, got {}",
            top[0].containment
        );
    }

    // QuerySketch::from_jam tests

    fn build_test_db_with_params(
        seqs: &[(&str, &str)],
        kmer_size: u8,
        fscale: u64,
        singleton: bool,
    ) -> (tempfile::TempDir, std::path::PathBuf) {
        let input = make_fasta(seqs);
        let output_dir = tempfile::tempdir().unwrap();
        let output_path = output_dir.path().join("test.jam");

        let config = BuildConfig {
            kmer_size,
            fscale,
            singleton,
            num_threads: 1,
            memory: 1,
            ..Default::default()
        };

        build(&[input.path().to_path_buf()], &output_path, &config).unwrap();
        (output_dir, output_path)
    }

    #[test]
    fn test_from_jam_success() {
        // Create two compatible .jam files with same parameters
        let (_dir1, db_path) = build_test_db_with_params(
            &[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            11,
            1,
            false,
        );
        let (_dir2, query_path) = build_test_db_with_params(
            &[("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")],
            11,
            1,
            false,
        );

        let db = JamReader::open(&db_path).unwrap();
        let sketch = QuerySketch::from_jam(&query_path, &db).unwrap();

        // Verify sketch was loaded
        assert_eq!(sketch.sample_count(), 1);
        assert!(sketch.total_entries() > 0);
        // Sample name comes from the stored names in the JAM file (FASTA filename in combined mode)
        assert!(!sketch.sample_names[0].is_empty());
        assert!(sketch.query_sizes[0] > 0);
    }

    #[test]
    fn test_from_jam_multiple_samples() {
        // Create query .jam with multiple samples (singleton mode)
        let (_dir1, db_path) = build_test_db_with_params(
            &[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            11,
            1,
            false,
        );
        let (_dir2, query_path) = build_test_db_with_params(
            &[
                ("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
                ("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
            ],
            11,
            1,
            true, // singleton mode - each sequence is a sample
        );

        let db = JamReader::open(&db_path).unwrap();
        let sketch = QuerySketch::from_jam(&query_path, &db).unwrap();

        assert_eq!(sketch.sample_count(), 2);
        // Sample names come from the stored names in the JAM file (record IDs in singleton mode)
        assert_eq!(sketch.sample_names[0], "seq1");
        assert_eq!(sketch.sample_names[1], "seq2");
        assert_eq!(sketch.query_sizes.len(), 2);
        assert!(sketch.query_sizes[0] > 0);
        assert!(sketch.query_sizes[1] > 0);
    }

    #[test]
    fn test_from_jam_kmer_size_mismatch() {
        // Create two .jam files with different k-mer sizes
        let (_dir1, db_path) = build_test_db_with_params(
            &[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            11, // k=11
            1,
            false,
        );
        let (_dir2, query_path) = build_test_db_with_params(
            &[("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")],
            21, // k=21 (different!)
            1,
            false,
        );

        let db = JamReader::open(&db_path).unwrap();
        let result = QuerySketch::from_jam(&query_path, &db);

        assert!(result.is_err());
        let err = result.unwrap_err();
        match err {
            QueryError::ParameterMismatch {
                parameter,
                source_value,
                target_value,
            } => {
                assert!(parameter.contains("k-mer"));
                assert_eq!(source_value, "21");
                assert_eq!(target_value, "11");
            }
            _ => panic!("Expected ParameterMismatch error, got {:?}", err),
        }
    }

    #[test]
    fn test_from_jam_threshold_mismatch() {
        // Create two .jam files with different fscale (threshold)
        let (_dir1, db_path) = build_test_db_with_params(
            &[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            11,
            1, // fscale=1 -> threshold = u64::MAX
            false,
        );
        let (_dir2, query_path) = build_test_db_with_params(
            &[("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")],
            11,
            1000, // fscale=1000 -> different threshold!
            false,
        );

        let db = JamReader::open(&db_path).unwrap();
        let result = QuerySketch::from_jam(&query_path, &db);

        assert!(result.is_err());
        let err = result.unwrap_err();
        match err {
            QueryError::ParameterMismatch {
                parameter,
                source_value,
                target_value,
            } => {
                assert!(parameter.contains("threshold"));
                // Just verify they're different - exact values depend on u64::MAX / fscale
                assert_ne!(source_value, target_value);
            }
            _ => panic!("Expected ParameterMismatch error, got {:?}", err),
        }
    }

    #[test]
    fn test_from_jam_preserves_bucketization() {
        // Verify that entries are correctly distributed to buckets
        let (_dir1, db_path) = build_test_db_with_params(
            &[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            11,
            1,
            false,
        );
        let (_dir2, query_path) = build_test_db_with_params(
            &[("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")],
            11,
            1,
            false,
        );

        let db = JamReader::open(&db_path).unwrap();
        let sketch = QuerySketch::from_jam(&query_path, &db).unwrap();

        // Verify all entries are in the correct bucket based on hash & 0xFF
        for bucket_idx in 0..BUCKET_COUNT {
            for &(hash, _sample_id) in sketch.bucket(bucket_idx) {
                assert_eq!(
                    bucket_id(hash),
                    bucket_idx,
                    "Entry with hash {} is in wrong bucket",
                    hash
                );
            }
        }
    }

    #[test]
    fn test_from_jam_query_sizes_correct() {
        // Create a query with known multiple samples
        let (_dir1, db_path) = build_test_db_with_params(
            &[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            11,
            1,
            false,
        );
        let (_dir2, query_path) = build_test_db_with_params(
            &[
                ("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
                ("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
            ],
            11,
            1,
            true, // singleton mode
        );

        let db = JamReader::open(&db_path).unwrap();
        let sketch = QuerySketch::from_jam(&query_path, &db).unwrap();

        // Verify query_sizes matches the actual unique hash count per sample
        for (sample_id, &expected_size) in sketch.query_sizes.iter().enumerate() {
            let mut unique_hashes = std::collections::HashSet::new();
            for bucket_idx in 0..BUCKET_COUNT {
                for &(hash, sid) in sketch.bucket(bucket_idx) {
                    if sid as usize == sample_id {
                        unique_hashes.insert(hash);
                    }
                }
            }
            assert_eq!(
                unique_hashes.len(),
                expected_size,
                "query_sizes[{}] should match actual unique hash count",
                sample_id
            );
        }
    }

    #[test]
    fn test_from_jam_empty_source() {
        // Test with a source that has minimal/no entries due to high fscale
        // Note: With very high fscale, most hashes are filtered out
        let (_dir1, db_path) = build_test_db_with_params(
            &[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            11,
            1_000_000, // Very high fscale
            false,
        );
        let (_dir2, query_path) = build_test_db_with_params(
            &[("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")],
            11,
            1_000_000, // Same high fscale
            false,
        );

        let db = JamReader::open(&db_path).unwrap();
        let result = QuerySketch::from_jam(&query_path, &db);

        // Should succeed even with few/no entries
        assert!(result.is_ok());
    }

    // QuerySketch::from_inputs tests

    #[test]
    fn test_from_inputs_empty() {
        let (_dir, db_path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();

        let sketch = QuerySketch::from_inputs(&[], &db, false).unwrap();

        assert_eq!(sketch.sample_count(), 0);
        assert_eq!(sketch.total_entries(), 0);
    }

    #[test]
    fn test_from_inputs_single_fasta() {
        let (_dir, db_path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();

        let query_fasta = make_fasta(&[("query_seq", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);

        let sketch =
            QuerySketch::from_inputs(&[query_fasta.path().to_path_buf()], &db, false).unwrap();

        assert_eq!(sketch.sample_count(), 1);
        assert!(sketch.total_entries() > 0);
    }

    #[test]
    fn test_from_inputs_single_jam() {
        let (_dir1, db_path) = build_test_db_with_params(
            &[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            11,
            1,
            false,
        );
        let (_dir2, query_jam) = build_test_db_with_params(
            &[("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")],
            11,
            1,
            false,
        );

        let db = JamReader::open(&db_path).unwrap();

        let sketch = QuerySketch::from_inputs(&[query_jam], &db, false).unwrap();

        assert_eq!(sketch.sample_count(), 1);
        assert!(sketch.total_entries() > 0);
    }

    #[test]
    fn test_from_inputs_multiple_fasta() {
        let (_dir, db_path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();

        let fasta1 = make_fasta(&[("query1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
        let fasta2 = make_fasta(&[("query2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")]);

        let sketch = QuerySketch::from_inputs(
            &[fasta1.path().to_path_buf(), fasta2.path().to_path_buf()],
            &db,
            false,
        )
        .unwrap();

        // Two FASTA files in non-singleton mode = 2 samples
        assert_eq!(sketch.sample_count(), 2);
        assert!(sketch.total_entries() > 0);
        assert_eq!(sketch.query_sizes.len(), 2);
    }

    #[test]
    fn test_from_inputs_multiple_fasta_singleton() {
        let (_dir, db_path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();

        // Each FASTA has 2 sequences
        let fasta1 = make_fasta(&[
            ("seq1a", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
            ("seq1b", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
        ]);
        let fasta2 = make_fasta(&[
            ("seq2a", "TATATATATATATATATATATATATATATATA"),
            ("seq2b", "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"),
        ]);

        let sketch = QuerySketch::from_inputs(
            &[fasta1.path().to_path_buf(), fasta2.path().to_path_buf()],
            &db,
            true, // singleton mode
        )
        .unwrap();

        // Two FASTA files with 2 seqs each in singleton mode = 4 samples
        assert_eq!(sketch.sample_count(), 4);
        assert_eq!(sketch.sample_names.len(), 4);
        assert_eq!(sketch.sample_names[0], "seq1a");
        assert_eq!(sketch.sample_names[1], "seq1b");
        assert_eq!(sketch.sample_names[2], "seq2a");
        assert_eq!(sketch.sample_names[3], "seq2b");
    }

    #[test]
    fn test_from_inputs_mixed_fasta_and_jam() {
        // Create a database
        let (_dir1, db_path) = build_test_db_with_params(
            &[("db_seq", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            11,
            1,
            false,
        );
        let db = JamReader::open(&db_path).unwrap();

        // Create a FASTA query
        let query_fasta = make_fasta(&[("fasta_query", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);

        // Create a .jam query with compatible parameters
        let (_dir2, query_jam) = build_test_db_with_params(
            &[("jam_query", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")],
            11,
            1,
            false,
        );

        let sketch =
            QuerySketch::from_inputs(&[query_fasta.path().to_path_buf(), query_jam], &db, false)
                .unwrap();

        // One FASTA + one .jam = 2 samples
        assert_eq!(sketch.sample_count(), 2);
        assert!(sketch.total_entries() > 0);
    }

    #[test]
    fn test_from_inputs_sample_id_renumbering() {
        // Verify sample IDs are renumbered to be contiguous
        let (_dir1, db_path) = build_test_db_with_params(
            &[("db_seq", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            11,
            1,
            false,
        );
        let db = JamReader::open(&db_path).unwrap();

        // Create two .jam files with singleton mode (multiple samples each)
        let (_dir2, jam1) = build_test_db_with_params(
            &[
                ("seq1a", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
                ("seq1b", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
            ],
            11,
            1,
            true, // singleton mode
        );
        let (_dir3, jam2) = build_test_db_with_params(
            &[
                ("seq2a", "TATATATATATATATATATATATATATATATA"),
                ("seq2b", "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"),
            ],
            11,
            1,
            true, // singleton mode
        );

        let sketch = QuerySketch::from_inputs(&[jam1, jam2], &db, false).unwrap();

        // 2 + 2 = 4 samples total
        assert_eq!(sketch.sample_count(), 4);

        // Verify sample IDs in buckets are in range [0, 4)
        for bucket in &sketch.buckets {
            for &(_hash, sample_id) in bucket {
                assert!(sample_id < 4, "Sample ID {} should be < 4", sample_id);
            }
        }

        // Verify we have entries for all sample IDs
        let mut seen_samples = std::collections::HashSet::new();
        for bucket in &sketch.buckets {
            for &(_hash, sample_id) in bucket {
                seen_samples.insert(sample_id);
            }
        }
        // All 4 samples should have at least some entries
        assert_eq!(seen_samples.len(), 4, "All samples should have entries");
    }

    #[test]
    fn test_from_inputs_buckets_sorted() {
        let (_dir, db_path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();

        let fasta1 = make_fasta(&[("q1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
        let fasta2 = make_fasta(&[("q2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")]);

        let sketch = QuerySketch::from_inputs(
            &[fasta1.path().to_path_buf(), fasta2.path().to_path_buf()],
            &db,
            false,
        )
        .unwrap();

        // Verify each bucket is sorted by (hash, sample_id)
        for bucket in &sketch.buckets {
            for window in bucket.windows(2) {
                assert!(
                    window[0] <= window[1],
                    "Bucket not sorted: {:?} > {:?}",
                    window[0],
                    window[1]
                );
            }
        }
    }

    #[test]
    fn test_from_inputs_query_sizes_preserved() {
        let (_dir, db_path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();

        let fasta1 = make_fasta(&[("q1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
        let fasta2 = make_fasta(&[("q2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")]);

        // Get individual sketches to compare query_sizes
        let sketch1 = QuerySketch::from_fasta(fasta1.path(), &db, false).unwrap();
        let sketch2 = QuerySketch::from_fasta(fasta2.path(), &db, false).unwrap();

        // Re-create FASTAs (NamedTempFile might have been moved)
        let fasta1_new = make_fasta(&[("q1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
        let fasta2_new = make_fasta(&[("q2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")]);

        let combined = QuerySketch::from_inputs(
            &[
                fasta1_new.path().to_path_buf(),
                fasta2_new.path().to_path_buf(),
            ],
            &db,
            false,
        )
        .unwrap();

        // Query sizes should match individual sketches
        assert_eq!(combined.query_sizes.len(), 2);
        assert_eq!(combined.query_sizes[0], sketch1.query_sizes[0]);
        assert_eq!(combined.query_sizes[1], sketch2.query_sizes[0]);
    }

    #[test]
    fn test_from_inputs_jam_detection_by_extension() {
        let (_dir1, db_path) = build_test_db_with_params(
            &[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            11,
            1,
            false,
        );
        let db = JamReader::open(&db_path).unwrap();

        // Create a .jam file - it should be detected by extension
        let (_dir2, jam_path) = build_test_db_with_params(
            &[("jam_seq", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")],
            11,
            1,
            false,
        );

        // Verify the file has .jam extension
        assert_eq!(jam_path.extension().unwrap(), "jam");

        let sketch = QuerySketch::from_inputs(&[jam_path], &db, false).unwrap();

        // Should load as .jam (sample names now come from the stored names in the JAM file)
        assert_eq!(sketch.sample_count(), 1);
        // In combined mode (singleton=false), the sample name is the FASTA filename
        assert!(!sketch.sample_names[0].is_empty());
    }

    #[test]
    fn test_from_inputs_propagates_errors() {
        let (_dir, db_path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();

        // Try to load a non-existent file
        let result = QuerySketch::from_inputs(
            &[std::path::PathBuf::from("/nonexistent/file.fasta")],
            &db,
            false,
        );

        assert!(result.is_err());
    }

    #[test]
    fn test_from_inputs_jam_parameter_mismatch_propagates() {
        let (_dir1, db_path) = build_test_db_with_params(
            &[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            11, // k=11
            1,
            false,
        );
        let db = JamReader::open(&db_path).unwrap();

        // Create a .jam with different k-mer size
        let (_dir2, jam_path) = build_test_db_with_params(
            &[("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")],
            21, // k=21 (mismatch!)
            1,
            false,
        );

        let result = QuerySketch::from_inputs(&[jam_path], &db, false);

        assert!(result.is_err());
        match result.unwrap_err() {
            QueryError::ParameterMismatch { parameter, .. } => {
                assert!(parameter.contains("k-mer"));
            }
            e => panic!("Expected ParameterMismatch error, got {:?}", e),
        }
    }

    #[test]
    fn test_from_inputs_integration_with_query_engine() {
        // Create a database with a known sequence
        let (_dir1, db_path) = build_test_db_with_params(
            &[("db_seq", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            11,
            1,
            false,
        );
        let db = JamReader::open(&db_path).unwrap();
        let engine = QueryEngine::open(&db_path).unwrap();

        // Create query inputs: one FASTA with same sequence, one .jam
        let query_fasta = make_fasta(&[("same_seq", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
        let (_dir2, query_jam) = build_test_db_with_params(
            &[("different_seq", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")],
            11,
            1,
            false,
        );

        let sketch =
            QuerySketch::from_inputs(&[query_fasta.path().to_path_buf(), query_jam], &db, false)
                .unwrap();

        assert_eq!(sketch.sample_count(), 2);

        // Query using the engine
        let results = engine.query_sketch(&sketch);

        assert_eq!(results.len(), 2);

        // First query (same sequence as DB) should have high containment
        assert!(results[0].has_matches());
        let top0 = results[0].top(1);
        assert!(!top0.is_empty());
        assert!(
            top0[0].containment >= 0.9,
            "Same sequence should have high containment, got {}",
            top0[0].containment
        );

        // Second query (different sequence) should have low/no containment
        // (Different sequences may still share some k-mers by chance)
    }

    // QueryEngine::query_fasta tests

    #[test]
    fn test_query_fasta_non_singleton() {
        // Build a database with a known sequence
        let (_dir, db_path) =
            build_test_db(&[("db_seq", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let engine = QueryEngine::open(&db_path).unwrap();

        // Create a query FASTA with the same sequence
        let query_fasta = make_fasta(&[("query_seq", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);

        // Query in non-singleton mode (entire file as one sample)
        let results = engine.query_fasta(query_fasta.path(), false).unwrap();

        // Should have exactly one result (one sample)
        assert_eq!(results.len(), 1);
        assert!(results[0].has_matches());

        // Should have high containment (same sequence)
        let top = results[0].top(1);
        assert!(!top.is_empty());
        assert!(
            top[0].containment >= 0.9,
            "Expected high containment, got {}",
            top[0].containment
        );
    }

    #[test]
    fn test_query_fasta_singleton() {
        // Build a database with two sequences (singleton mode)
        let (_dir, db_path) = build_test_db(
            &[
                ("db_seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
                ("db_seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
            ],
            true,
        );
        let engine = QueryEngine::open(&db_path).unwrap();

        // Create a query FASTA with two sequences
        let query_fasta = make_fasta(&[
            ("query1", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
            ("query2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
        ]);

        // Query in singleton mode (each sequence as a separate sample)
        let results = engine.query_fasta(query_fasta.path(), true).unwrap();

        // Should have two results (one per sequence)
        assert_eq!(results.len(), 2);

        // Both should have matches
        assert!(results[0].has_matches());
        assert!(results[1].has_matches());

        // Each query should have high containment with its corresponding DB sample
        for (i, result) in results.iter().enumerate() {
            let self_match = result.matches.iter().find(|m| m.sample_id == i as u32);
            if let Some(m) = self_match {
                assert!(
                    m.containment >= 0.9,
                    "Query {} should have high containment with DB sample {}, got {}",
                    i,
                    i,
                    m.containment
                );
            }
        }
    }

    #[test]
    fn test_query_fasta_file_not_found() {
        // Build a database
        let (_dir, db_path) =
            build_test_db(&[("db_seq", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let engine = QueryEngine::open(&db_path).unwrap();

        // Try to query a non-existent file
        let result = engine.query_fasta("/nonexistent/path.fasta", false);

        assert!(result.is_err());
        match result.unwrap_err() {
            QueryError::Parse { path, message: _ } => {
                assert!(path.contains("nonexistent"));
            }
            e => panic!("Expected Parse error, got {:?}", e),
        }
    }
}
