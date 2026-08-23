use crate::bias::HashBiasTable;
use crate::core_utils::passes_entropy_filter;
use crate::format::{BUCKET_COUNT, bucket_id};
use crate::jamhash_u64_v1;
use crate::reader::{JamReader, ReaderError, lower_bound_hash};
use needletail::{Sequence, parse_fastx_file};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::Arc;

#[derive(Debug)]
pub struct QuerySketch {
    pub buckets: [Vec<(u64, u32)>; BUCKET_COUNT],
    pub sample_names: Vec<String>,
    pub query_sizes: Vec<usize>,
    pub query_weight_sums: Vec<f64>,
}

impl QuerySketch {
    pub fn new() -> Self {
        Self {
            buckets: std::array::from_fn(|_| Vec::new()),
            sample_names: Vec::new(),
            query_sizes: Vec::new(),
            query_weight_sums: Vec::new(),
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

        if source.min_entropy().to_bits() != db.min_entropy().to_bits() {
            return Err(QueryError::ParameterMismatch {
                parameter: "minimum entropy".to_string(),
                source_value: source.min_entropy().to_string(),
                target_value: db.min_entropy().to_string(),
            });
        }

        match (source.bias_table(), db.bias_table()) {
            (None, None) => {}
            (Some(source_bias), Some(db_bias)) if source_bias.as_ref() == db_bias.as_ref() => {}
            (Some(_), Some(_)) => {
                return Err(QueryError::ParameterMismatch {
                    parameter: "bias table".to_string(),
                    source_value: "different bias table".to_string(),
                    target_value: "database bias table".to_string(),
                });
            }
            (Some(_), None) => {
                return Err(QueryError::ParameterMismatch {
                    parameter: "bias table".to_string(),
                    source_value: "present".to_string(),
                    target_value: "absent".to_string(),
                });
            }
            (None, Some(_)) => {
                return Err(QueryError::ParameterMismatch {
                    parameter: "bias table".to_string(),
                    source_value: "absent".to_string(),
                    target_value: "present".to_string(),
                });
            }
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

        let mut buckets: [Vec<(u64, u32)>; BUCKET_COUNT] = std::array::from_fn(|_| Vec::new());
        let mut query_sizes = vec![0usize; expected_sample_count];
        for (bucket_idx, bucket) in buckets.iter_mut().enumerate() {
            for entry in source.bucket_entries(bucket_idx) {
                if entry.hash == 0 {
                    continue;
                }
                bucket.push((entry.hash, entry.sample_id));
                query_sizes[entry.sample_id as usize] += 1;
            }
        }

        let query_weight_sums = if let Some(ref bt) = db.bias_table() {
            let mut sums = vec![0.0f64; query_sizes.len()];
            for bucket in &buckets {
                for &(hash, sample_id) in bucket {
                    sums[sample_id as usize] += bt.effective_fscale_at(bt.weight(hash));
                }
            }
            sums
        } else {
            vec![0.0; query_sizes.len()]
        };
        Ok(Self {
            buckets,
            sample_names,
            query_sizes,
            query_weight_sums,
        })
    }

    pub fn from_fasta<P: AsRef<Path>>(
        input: P,
        db: &JamReader,
        singleton: bool,
    ) -> Result<Self, QueryError> {
        let input_path = input.as_ref();
        let kmer_size = db.kmer_size();
        let threshold = match db.bias_table() {
            Some(ref bt) if bt.is_soft_filter() => u64::MAX / bt.min_fscale(),
            _ => db.threshold(),
        };
        let min_entropy = db.min_entropy();
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
        let mut query_sizes: Vec<usize> = Vec::new();
        let mut combined_hashes: HashSet<u64> = HashSet::new();
        let mut weight_sums: Vec<f64> = Vec::new();
        let mut current_sample_id: u32 = 0;

        if !singleton {
            sample_names.push(
                input_path
                    .file_name()
                    .and_then(|s| s.to_str())
                    .unwrap_or("query")
                    .to_string(),
            );
            weight_sums.push(0.0);
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
                weight_sums.push(0.0);
                current_sample_id = u32::try_from(sample_names.len() - 1).map_err(|_| {
                    QueryError::Config("query sample count exceeds u32::MAX".to_string())
                })?;
            }

            let sequence = record.normalize(false);
            let mut record_hashes = HashSet::new();
            if sequence.len() < kmer_size as usize {
                if singleton {
                    query_sizes.push(0);
                }
                continue;
            }

            for (_, kmer, _) in sequence.bit_kmers(kmer_size, true) {
                let hash = jamhash_u64_v1(kmer.0);

                if hash == 0 || hash >= threshold {
                    continue;
                }

                if min_entropy > 0.0 && !passes_entropy_filter(kmer.0, kmer_size, min_entropy) {
                    continue;
                }

                if bias_table.as_ref().is_some_and(|b| !b.passes_filter(hash)) {
                    continue;
                }

                let hashes = if singleton {
                    &mut record_hashes
                } else {
                    &mut combined_hashes
                };
                if hashes.insert(hash) {
                    buckets[bucket_id(hash)].push((hash, current_sample_id));
                    if let Some(ref bt) = bias_table {
                        weight_sums[current_sample_id as usize] +=
                            bt.effective_fscale_at(bt.weight(hash));
                    }
                }
            }
            if singleton {
                query_sizes.push(record_hashes.len());
            }
        }

        for bucket in &mut buckets {
            bucket.sort_unstable();
            bucket.dedup();
        }

        if !singleton {
            query_sizes.push(combined_hashes.len());
        }

        Ok(Self {
            buckets,
            sample_names,
            query_sizes,
            query_weight_sums: weight_sums,
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

            let sample_offset = u32::try_from(combined.sample_count()).map_err(|_| {
                QueryError::Config("combined query sample count exceeds u32::MAX".to_string())
            })?;
            combined.sample_names.extend(sketch.sample_names);
            combined.query_sizes.extend(sketch.query_sizes);
            combined.query_weight_sums.extend(sketch.query_weight_sums);

            for (bucket_idx, bucket) in sketch.buckets.into_iter().enumerate() {
                for (hash, sample_id) in bucket {
                    let combined_sample_id =
                        sample_id.checked_add(sample_offset).ok_or_else(|| {
                            QueryError::Config(
                                "combined query sample ID exceeds u32::MAX".to_string(),
                            )
                        })?;
                    combined.buckets[bucket_idx].push((hash, combined_sample_id));
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

    #[error("Invalid query configuration: {0}")]
    Config(String),
}

#[derive(Debug, Clone)]
pub struct SampleMatch {
    pub sample_id: u32,
    pub hit_count: u32,
    pub containment: f64,
    pub hit_weight: f64,
    pub e_value: f64,
}

#[derive(Debug, Clone)]
pub struct QueryResult {
    pub query_size: usize,
    pub hashes_found: usize,
    pub matches: Vec<SampleMatch>,
    pub failed_bucket_count: usize,
    pub total_query_weight: f64,
}

type SampleHitMap = HashMap<u32, (u32, f64)>;
type QueryBucketHitMap = HashMap<usize, SampleHitMap>;

impl QueryResult {
    pub fn top(&self, n: usize) -> Vec<&SampleMatch> {
        let mut sorted: Vec<_> = self.matches.iter().collect();
        sorted.sort_by(|a, b| {
            b.containment
                .total_cmp(&a.containment)
                .then_with(|| b.hit_count.cmp(&a.hit_count))
                .then_with(|| a.sample_id.cmp(&b.sample_id))
        });
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

    pub fn is_partial(&self) -> bool {
        self.failed_bucket_count > 0
    }
}

fn compute_e_value(
    hit_count: u32,
    query_size: usize,
    db_size: u64,
    threshold: u64,
    num_db_samples: u32,
) -> f64 {
    use statrs::distribution::{DiscreteCDF, Poisson};

    if hit_count == 0 {
        return num_db_samples as f64;
    }
    if query_size == 0 || threshold == 0 {
        return 0.0;
    }

    let lambda = query_size as f64 * db_size as f64 / threshold as f64;
    if lambda == 0.0 {
        return 0.0;
    }

    let p = match Poisson::new(lambda) {
        Ok(pois) => 1.0 - pois.cdf(hit_count.saturating_sub(1) as u64),
        Err(_) => 1.0,
    };

    p * num_db_samples as f64
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
                total_query_weight: 0.0,
            };
        }

        let mut sorted_hashes = hashes.to_vec();
        sorted_hashes.retain(|&hash| hash != 0);
        sorted_hashes.sort_unstable_by_key(|&h| (h & 0xFF, h));
        sorted_hashes.dedup();

        let mut sample_hits: HashMap<u32, (u32, f64)> = HashMap::new();
        let mut hashes_found = 0;

        for &hash in &sorted_hashes {
            let hash_weight = self
                .bias_table
                .as_ref()
                .map(|bt| bt.effective_fscale_at(bt.weight(hash)))
                .unwrap_or(0.0);
            let mut found = false;
            for sample_id in self.reader.search(hash) {
                let entry = sample_hits.entry(sample_id).or_insert((0, 0.0));
                entry.0 += 1;
                entry.1 += hash_weight;
                found = true;
            }
            if found {
                hashes_found += 1;
            }
        }

        let query_size = sorted_hashes.len();
        let threshold = self.reader.threshold();
        let sample_sizes = self.reader.sample_sizes();
        let num_db_samples = self.reader.stats().sample_count;
        let total_query_weight = self
            .bias_table
            .as_ref()
            .map(|bt| {
                sorted_hashes
                    .iter()
                    .map(|&h| bt.effective_fscale_at(bt.weight(h)))
                    .sum()
            })
            .unwrap_or(0.0);

        let matches: Vec<SampleMatch> = sample_hits
            .into_iter()
            .map(|(sample_id, (hit_count, hit_weight))| SampleMatch {
                sample_id,
                hit_count,
                containment: if total_query_weight > 0.0 && hit_weight > 0.0 {
                    hit_weight / total_query_weight
                } else if query_size > 0 {
                    hit_count as f64 / query_size as f64
                } else {
                    0.0
                },
                hit_weight,
                e_value: compute_e_value(
                    hit_count,
                    query_size,
                    sample_sizes.get(sample_id as usize).copied().unwrap_or(0),
                    threshold,
                    num_db_samples,
                ),
            })
            .collect();

        QueryResult {
            query_size,
            hashes_found,
            matches,
            failed_bucket_count: 0,
            total_query_weight,
        }
    }

    pub fn query_filtered(&self, hashes: &[u64], min_containment: f64) -> QueryResult {
        let mut result = self.query(hashes);
        result.matches.retain(|m| m.containment >= min_containment);
        result.matches.sort_by(|a, b| {
            b.containment
                .total_cmp(&a.containment)
                .then_with(|| b.hit_count.cmp(&a.hit_count))
                .then_with(|| a.sample_id.cmp(&b.sample_id))
        });
        result
    }

    pub fn query_batch(&self, queries: &[Vec<u64>]) -> Vec<QueryResult> {
        use rayon::prelude::*;
        queries.par_iter().map(|q| self.query(q)).collect()
    }

    pub fn query_sketch(&self, sketch: &QuerySketch) -> Vec<QueryResult> {
        self.query_sketch_chunked(sketch, 0..sketch.sample_count(), 0.0)
    }

    pub fn query_sketch_chunked(
        &self,
        sketch: &QuerySketch,
        sample_range: std::ops::Range<usize>,
        min_containment: f64,
    ) -> Vec<QueryResult> {
        use crate::format::{ENTRY_SIZE, PAGE_SIZE};
        use rayon::prelude::*;
        use std::collections::HashMap;
        use std::sync::atomic::{AtomicU32, Ordering};

        let chunk_len = sample_range.len();
        if chunk_len == 0 {
            return Vec::new();
        }

        let range_start = sample_range.start;
        let threshold = self.reader.threshold();

        self.reader.advise_random();

        let hashes_found: Vec<AtomicU32> = (0..chunk_len).map(|_| AtomicU32::new(0)).collect();

        let bias_table = &self.bias_table;

        let bucket_maps: Vec<QueryBucketHitMap> = (0..BUCKET_COUNT)
            .into_par_iter()
            .map(|bucket_idx| {
                let mut local_maps: QueryBucketHitMap = HashMap::new();

                let query_bucket = sketch.bucket(bucket_idx);
                if query_bucket.is_empty() {
                    return local_maps;
                }

                let filter = match self.reader.bucket_filter(bucket_idx) {
                    Some(f) => f,
                    None => return local_maps,
                };

                let mut survivors: Vec<(u64, u32)> = Vec::with_capacity(query_bucket.len() / 10);
                let mut prev_hash = u64::MAX;
                let mut prev_passed = false;

                for &(hash, sample_id) in query_bucket {
                    if hash != prev_hash {
                        prev_hash = hash;
                        prev_passed = filter.contains(&hash);
                    }
                    if prev_passed
                        && (sample_id as usize) >= range_start
                        && (sample_id as usize) < range_start + chunk_len
                    {
                        survivors.push((hash, sample_id));
                    }
                }

                if survivors.is_empty() {
                    return local_maps;
                }

                let db_bucket = self.reader.bucket_entries(bucket_idx);
                let count = db_bucket.len();
                if count == 0 {
                    return local_maps;
                }

                let (entry_start, entry_end) = self.reader.bucket_entry_byte_range(bucket_idx);
                let mut last_released_page = entry_start & !(PAGE_SIZE - 1);

                let mut q_idx = 0;
                while q_idx < survivors.len() {
                    let q_hash = survivors[q_idx].0;

                    let d_idx = lower_bound_hash(db_bucket, q_hash);

                    let current_byte = entry_start + d_idx * ENTRY_SIZE;
                    let current_page = current_byte & !(PAGE_SIZE - 1);
                    if current_page > last_released_page + PAGE_SIZE {
                        self.reader
                            .release_pages(last_released_page, current_page - PAGE_SIZE);
                        last_released_page = current_page - PAGE_SIZE;
                    }

                    let db_start = d_idx;
                    let mut db_end = d_idx;
                    while db_end < count && db_bucket[db_end].hash == q_hash {
                        db_end += 1;
                    }
                    let has_matches = db_start < db_end;

                    let hash_weight = bias_table
                        .as_ref()
                        .map(|bt| bt.effective_fscale_at(bt.weight(q_hash)))
                        .unwrap_or(0.0);

                    let mut prev_sample = u32::MAX;
                    while q_idx < survivors.len() && survivors[q_idx].0 == q_hash {
                        let q_sample = survivors[q_idx].1;

                        if q_sample != prev_sample {
                            if has_matches {
                                let local_idx = (q_sample as usize) - range_start;
                                let sample_map = local_maps.entry(local_idx).or_default();
                                for db_entry in &db_bucket[db_start..db_end] {
                                    let entry =
                                        sample_map.entry(db_entry.sample_id).or_insert((0, 0.0));
                                    entry.0 += 1;
                                    entry.1 += hash_weight;
                                }
                                hashes_found[local_idx].fetch_add(1, Ordering::Relaxed);
                            }
                            prev_sample = q_sample;
                        }
                        q_idx += 1;
                    }
                }

                self.reader.release_pages(last_released_page, entry_end);

                local_maps
            })
            .collect();

        let merged: Vec<SampleHitMap> = (0..chunk_len)
            .into_par_iter()
            .map(|local_idx| {
                let mut combined: SampleHitMap = HashMap::new();
                for per_bucket in &bucket_maps {
                    if let Some(local_map) = per_bucket.get(&local_idx) {
                        for (&db_id, &(count, weight)) in local_map {
                            let entry = combined.entry(db_id).or_insert((0, 0.0));
                            entry.0 += count;
                            entry.1 += weight;
                        }
                    }
                }
                combined
            })
            .collect();

        let sample_sizes = self.reader.sample_sizes();
        let num_db_samples = self.reader.stats().sample_count;

        (0..chunk_len)
            .into_par_iter()
            .map(|local_idx| {
                let global_idx = range_start + local_idx;
                let query_size = sketch.query_sizes[global_idx];

                let matches: Vec<SampleMatch> = merged[local_idx]
                    .iter()
                    .filter_map(|(&db_id, &(count, weight))| {
                        let containment = {
                            let total_w = sketch
                                .query_weight_sums
                                .get(global_idx)
                                .copied()
                                .unwrap_or(0.0);
                            if total_w > 0.0 && weight > 0.0 {
                                weight / total_w
                            } else if query_size > 0 {
                                count as f64 / query_size as f64
                            } else {
                                0.0
                            }
                        };
                        (containment >= min_containment).then(|| SampleMatch {
                            sample_id: db_id,
                            hit_count: count,
                            containment,
                            hit_weight: weight,
                            e_value: compute_e_value(
                                count,
                                query_size,
                                sample_sizes.get(db_id as usize).copied().unwrap_or(0),
                                threshold,
                                num_db_samples,
                            ),
                        })
                    })
                    .collect();

                QueryResult {
                    query_size,
                    hashes_found: hashes_found[local_idx].load(Ordering::Relaxed) as usize,
                    matches,
                    failed_bucket_count: 0,
                    total_query_weight: sketch
                        .query_weight_sums
                        .get(global_idx)
                        .copied()
                        .unwrap_or(0.0),
                }
            })
            .collect()
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
    use crate::bias::{BiasCreateConfig, CMSConfig, HashBiasTable};
    use crate::writer::{BuildConfig, build};
    use std::io::Write;
    use std::sync::Arc;
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
            let result = engine.query_filtered(&test_hashes, 0.5);
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
                    hit_weight: 0.0,
                    e_value: 0.0,
                },
                SampleMatch {
                    sample_id: 1,
                    hit_count: 30,
                    containment: 0.3,
                    hit_weight: 0.0,
                    e_value: 0.0,
                },
                SampleMatch {
                    sample_id: 2,
                    hit_count: 80,
                    containment: 0.8,
                    hit_weight: 0.0,
                    e_value: 0.0,
                },
            ],
            failed_bucket_count: 0,
            total_query_weight: 0.0,
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

        sketch.buckets[0].push((100, 0));
        sketch.buckets[0].push((200, 1));

        sketch.buckets[255].push((300, 0));

        let bucket_0 = sketch.bucket(0);
        assert_eq!(bucket_0.len(), 2);
        assert_eq!(bucket_0[0], (100, 0));
        assert_eq!(bucket_0[1], (200, 1));

        let bucket_255 = sketch.bucket(255);
        assert_eq!(bucket_255.len(), 1);
        assert_eq!(bucket_255[0], (300, 0));

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

        sketch.sample_names = vec!["query_sample_1".to_string(), "query_sample_2".to_string()];

        sketch.query_sizes = vec![1000, 500];

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
        let _ = sketch.bucket(256);
    }

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
        let (_dir, path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let engine = QueryEngine::open(&path).unwrap();
        let reader = JamReader::open(&path).unwrap();

        let mut sketch = QuerySketch::new();
        sketch.sample_names.push("query_sample".to_string());

        let mut unique_hashes = std::collections::HashSet::new();
        for bucket_idx in 0..256 {
            for entry in reader.bucket_entries(bucket_idx) {
                if unique_hashes.insert(entry.hash) {
                    sketch.buckets[bucket_idx].push((entry.hash, 0));
                }
            }
        }
        sketch.query_sizes.push(unique_hashes.len());

        let results = engine.query_sketch(&sketch);

        assert_eq!(results.len(), 1);
        assert!(results[0].has_matches());

        let db_sample_0_match = results[0].matches.iter().find(|m| m.sample_id == 0);
        assert!(db_sample_0_match.is_some(), "Should match db sample 0");

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
        let (_dir, path) = build_test_db(
            &[
                ("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
                ("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
            ],
            true,
        );
        let engine = QueryEngine::open(&path).unwrap();
        let reader = JamReader::open(&path).unwrap();

        let mut sketch = QuerySketch::new();
        sketch.sample_names.push("query_0".to_string());
        sketch.sample_names.push("query_1".to_string());

        let mut hashes_per_sample: [std::collections::HashSet<u64>; 2] = Default::default();

        for bucket_idx in 0..256 {
            for entry in reader.bucket_entries(bucket_idx) {
                let query_sample_id = entry.sample_id;
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

        for (query_idx, result) in results.iter().enumerate() {
            assert!(result.has_matches());
            let self_match = result
                .matches
                .iter()
                .find(|m| m.sample_id == query_idx as u32);
            if let Some(m) = self_match {
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

        let mut sketch = QuerySketch::new();
        sketch.sample_names.push("fake_sample".to_string());
        sketch.query_sizes.push(10);

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

        let mut sketch = QuerySketch::new();
        sketch.sample_names.push("half_sample".to_string());

        let mut all_hashes = Vec::new();
        for bucket_idx in 0..256 {
            for entry in reader.bucket_entries(bucket_idx) {
                all_hashes.push((entry.hash, bucket_idx));
            }
        }

        let selected_hashes: Vec<_> = all_hashes.iter().step_by(2).collect();
        sketch.query_sizes.push(selected_hashes.len());

        for &(hash, bucket_idx) in &selected_hashes {
            sketch.buckets[*bucket_idx].push((*hash, 0));
        }

        let results = engine.query_sketch(&sketch);

        assert_eq!(results.len(), 1);
        assert!(results[0].has_matches());
        let top = results[0].top(1);
        assert!(!top.is_empty());
        assert!(top[0].containment >= 0.9);
    }

    #[test]
    fn test_from_fasta_non_singleton() {
        let (_dir, db_path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();

        let query_fasta = make_fasta(&[("query_seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);

        let sketch = QuerySketch::from_fasta(query_fasta.path(), &db, false).unwrap();

        assert_eq!(sketch.sample_count(), 1);
        assert!(!sketch.sample_names[0].is_empty());

        assert!(sketch.total_entries() > 0);
        assert!(sketch.query_sizes[0] > 0);

        assert_eq!(sketch.query_sizes[0], sketch.total_entries());
    }

    #[test]
    fn test_from_fasta_singleton() {
        let (_dir, db_path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();

        let query_fasta = make_fasta(&[
            ("query_seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
            ("query_seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
        ]);

        let sketch = QuerySketch::from_fasta(query_fasta.path(), &db, true).unwrap();

        assert_eq!(sketch.sample_count(), 2);
        assert_eq!(sketch.sample_names[0], "query_seq1");
        assert_eq!(sketch.sample_names[1], "query_seq2");

        assert!(sketch.query_sizes[0] > 0);
        assert!(sketch.query_sizes[1] > 0);

        let total_unique: usize = sketch.query_sizes.iter().sum();
        assert!(total_unique <= sketch.total_entries() + sketch.sample_count());
    }

    #[test]
    fn test_from_fasta_uses_db_parameters() {
        let input = make_fasta(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
        let output_dir = tempfile::tempdir().unwrap();
        let db_path = output_dir.path().join("test.jam");

        let config = BuildConfig {
            kmer_size: 15,
            fscale: 10,
            singleton: false,
            num_threads: 1,
            memory: 1,
            ..Default::default()
        };

        build(&[input.path().to_path_buf()], &db_path, &config).unwrap();
        let db = JamReader::open(&db_path).unwrap();

        assert_eq!(db.kmer_size(), 15);

        let query_fasta = make_fasta(&[("query", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);

        let sketch = QuerySketch::from_fasta(query_fasta.path(), &db, false).unwrap();

        assert!(sketch.sample_count() == 1);
    }

    #[test]
    fn test_from_fasta_deduplication() {
        let (_dir, db_path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();

        let query_fasta = make_fasta(&[(
            "query",
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        )]);

        let sketch = QuerySketch::from_fasta(query_fasta.path(), &db, false).unwrap();

        assert_eq!(sketch.query_sizes[0], sketch.total_entries());
    }

    #[test]
    fn test_from_fasta_bucketization() {
        let (_dir, db_path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();

        let query_fasta = make_fasta(&[(
            "query",
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        )]);

        let sketch = QuerySketch::from_fasta(query_fasta.path(), &db, false).unwrap();

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
        let (_dir, db_path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();
        assert_eq!(db.kmer_size(), 11);

        let query_fasta = make_fasta(&[
            ("short", "ATCGATCG"),
            ("long", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
        ]);

        let sketch = QuerySketch::from_fasta(query_fasta.path(), &db, true).unwrap();

        assert_eq!(sketch.sample_count(), 2);

        assert_eq!(sketch.query_sizes[0], 0);

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
        let (_dir, db_path) =
            build_test_db(&[("db_seq", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();
        let engine = QueryEngine::open(&db_path).unwrap();

        let query_fasta = make_fasta(&[("query_seq", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);

        let sketch = QuerySketch::from_fasta(query_fasta.path(), &db, false).unwrap();

        let results = engine.query_sketch(&sketch);

        assert_eq!(results.len(), 1);
        assert!(results[0].has_matches());

        let top = results[0].top(1);
        assert!(!top.is_empty());
        assert!(
            top[0].containment >= 0.9,
            "Expected high containment, got {}",
            top[0].containment
        );
    }

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

    fn build_test_db_with_config(
        seqs: &[(&str, &str)],
        config: BuildConfig,
    ) -> (tempfile::TempDir, std::path::PathBuf) {
        let input = make_fasta(seqs);
        let output_dir = tempfile::tempdir().unwrap();
        let output_path = output_dir.path().join("test.jam");
        build(&[input.path().to_path_buf()], &output_path, &config).unwrap();
        (output_dir, output_path)
    }

    fn make_bias_table(k: u8, fscale: u64) -> HashBiasTable {
        let pos = make_fasta(&[("pos", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
        let neg = make_fasta(&[("neg", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")]);
        let config = BiasCreateConfig {
            cms: CMSConfig {
                width: 1024,
                depth: 3,
                k,
                fscale,
            },
            alpha: 1.0,
            target_fscale: None,
            negative_fscale: None,
            unseen_fscale: None,
        };
        HashBiasTable::create(&[pos.path()], &[neg.path()], &config, None).unwrap()
    }

    #[test]
    fn test_query_deduplicates_direct_hashes() {
        let (_dir, path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let reader = JamReader::open(&path).unwrap();
        let hash = (0..BUCKET_COUNT)
            .find_map(|bucket_idx| {
                reader
                    .bucket_entries(bucket_idx)
                    .first()
                    .map(|entry| entry.hash)
            })
            .expect("test database should contain at least one hash");

        let engine = QueryEngine::open(&path).unwrap();
        let result = engine.query(&[hash, hash, hash]);

        assert_eq!(result.query_size, 1);
        assert_eq!(result.hashes_found, 1);
        let m = result.matches.iter().find(|m| m.sample_id == 0).unwrap();
        assert_eq!(m.hit_count, 1);
    }

    #[test]
    fn test_from_jam_min_entropy_mismatch() {
        let base_config = BuildConfig {
            kmer_size: 11,
            fscale: 1,
            num_threads: 1,
            memory: 1,
            ..Default::default()
        };
        let (_source_dir, source_path) = build_test_db_with_config(
            &[("source", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            base_config.clone(),
        );

        let mut target_config = base_config;
        target_config.min_entropy = 1.0;
        let (_target_dir, target_path) = build_test_db_with_config(
            &[("target", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            target_config,
        );

        let db = JamReader::open(&target_path).unwrap();
        let err = QuerySketch::from_jam(&source_path, &db).unwrap_err();
        match err {
            QueryError::ParameterMismatch { parameter, .. } => {
                assert!(parameter.contains("entropy"));
            }
            e => panic!("Expected ParameterMismatch error, got {:?}", e),
        }
    }

    #[test]
    fn test_from_jam_bias_presence_mismatch() {
        let bias_table = make_bias_table(11, 1);
        let (_source_dir, source_path) = build_test_db_with_params(
            &[("source", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            11,
            1,
            false,
        );

        let target_config = BuildConfig {
            kmer_size: 11,
            fscale: 1,
            num_threads: 1,
            memory: 1,
            bias_table: Some(Arc::new(bias_table)),
            ..Default::default()
        };
        let (_target_dir, target_path) = build_test_db_with_config(
            &[("target", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            target_config,
        );

        let db = JamReader::open(&target_path).unwrap();
        let err = QuerySketch::from_jam(&source_path, &db).unwrap_err();
        match err {
            QueryError::ParameterMismatch { parameter, .. } => {
                assert!(parameter.contains("bias"));
            }
            e => panic!("Expected ParameterMismatch error, got {:?}", e),
        }
    }

    #[test]
    fn test_from_jam_success() {
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

        assert_eq!(sketch.sample_count(), 1);
        assert!(sketch.total_entries() > 0);
        assert!(!sketch.sample_names[0].is_empty());
        assert!(sketch.query_sizes[0] > 0);
    }

    #[test]
    fn test_from_jam_multiple_samples() {
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
            true,
        );

        let db = JamReader::open(&db_path).unwrap();
        let sketch = QuerySketch::from_jam(&query_path, &db).unwrap();

        assert_eq!(sketch.sample_count(), 2);
        assert_eq!(sketch.sample_names[0], "seq1");
        assert_eq!(sketch.sample_names[1], "seq2");
        assert_eq!(sketch.query_sizes.len(), 2);
        assert!(sketch.query_sizes[0] > 0);
        assert!(sketch.query_sizes[1] > 0);
    }

    #[test]
    fn test_from_jam_kmer_size_mismatch() {
        let (_dir1, db_path) = build_test_db_with_params(
            &[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            11,
            1,
            false,
        );
        let (_dir2, query_path) = build_test_db_with_params(
            &[("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")],
            21,
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
        let (_dir1, db_path) = build_test_db_with_params(
            &[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            11,
            1,
            false,
        );
        let (_dir2, query_path) = build_test_db_with_params(
            &[("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")],
            11,
            1000,
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
                assert_ne!(source_value, target_value);
            }
            _ => panic!("Expected ParameterMismatch error, got {:?}", err),
        }
    }

    #[test]
    fn test_from_jam_preserves_bucketization() {
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
            true,
        );

        let db = JamReader::open(&db_path).unwrap();
        let sketch = QuerySketch::from_jam(&query_path, &db).unwrap();

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
        let (_dir1, db_path) = build_test_db_with_params(
            &[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            11,
            1_000_000,
            false,
        );
        let (_dir2, query_path) = build_test_db_with_params(
            &[("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")],
            11,
            1_000_000,
            false,
        );

        let db = JamReader::open(&db_path).unwrap();
        let result = QuerySketch::from_jam(&query_path, &db);

        assert!(result.is_ok());
    }

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

        assert_eq!(sketch.sample_count(), 2);
        assert!(sketch.total_entries() > 0);
        assert_eq!(sketch.query_sizes.len(), 2);
    }

    #[test]
    fn test_from_inputs_multiple_fasta_singleton() {
        let (_dir, db_path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();

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
            true,
        )
        .unwrap();

        assert_eq!(sketch.sample_count(), 4);
        assert_eq!(sketch.sample_names.len(), 4);
        assert_eq!(sketch.sample_names[0], "seq1a");
        assert_eq!(sketch.sample_names[1], "seq1b");
        assert_eq!(sketch.sample_names[2], "seq2a");
        assert_eq!(sketch.sample_names[3], "seq2b");
    }

    #[test]
    fn test_from_inputs_mixed_fasta_and_jam() {
        let (_dir1, db_path) = build_test_db_with_params(
            &[("db_seq", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            11,
            1,
            false,
        );
        let db = JamReader::open(&db_path).unwrap();

        let query_fasta = make_fasta(&[("fasta_query", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);

        let (_dir2, query_jam) = build_test_db_with_params(
            &[("jam_query", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")],
            11,
            1,
            false,
        );

        let sketch =
            QuerySketch::from_inputs(&[query_fasta.path().to_path_buf(), query_jam], &db, false)
                .unwrap();

        assert_eq!(sketch.sample_count(), 2);
        assert!(sketch.total_entries() > 0);
    }

    #[test]
    fn test_from_inputs_sample_id_renumbering() {
        let (_dir1, db_path) = build_test_db_with_params(
            &[("db_seq", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            11,
            1,
            false,
        );
        let db = JamReader::open(&db_path).unwrap();

        let (_dir2, jam1) = build_test_db_with_params(
            &[
                ("seq1a", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
                ("seq1b", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
            ],
            11,
            1,
            true,
        );
        let (_dir3, jam2) = build_test_db_with_params(
            &[
                ("seq2a", "TATATATATATATATATATATATATATATATA"),
                ("seq2b", "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"),
            ],
            11,
            1,
            true,
        );

        let sketch = QuerySketch::from_inputs(&[jam1, jam2], &db, false).unwrap();

        assert_eq!(sketch.sample_count(), 4);

        for bucket in &sketch.buckets {
            for &(_hash, sample_id) in bucket {
                assert!(sample_id < 4, "Sample ID {} should be < 4", sample_id);
            }
        }

        let mut seen_samples = std::collections::HashSet::new();
        for bucket in &sketch.buckets {
            for &(_hash, sample_id) in bucket {
                seen_samples.insert(sample_id);
            }
        }
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

        let sketch1 = QuerySketch::from_fasta(fasta1.path(), &db, false).unwrap();
        let sketch2 = QuerySketch::from_fasta(fasta2.path(), &db, false).unwrap();

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

        let (_dir2, jam_path) = build_test_db_with_params(
            &[("jam_seq", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")],
            11,
            1,
            false,
        );

        assert_eq!(jam_path.extension().unwrap(), "jam");

        let sketch = QuerySketch::from_inputs(&[jam_path], &db, false).unwrap();

        assert_eq!(sketch.sample_count(), 1);
        assert!(!sketch.sample_names[0].is_empty());
    }

    #[test]
    fn test_from_inputs_propagates_errors() {
        let (_dir, db_path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let db = JamReader::open(&db_path).unwrap();

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
            11,
            1,
            false,
        );
        let db = JamReader::open(&db_path).unwrap();

        let (_dir2, jam_path) = build_test_db_with_params(
            &[("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")],
            21,
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

    fn normalize_results(results: &mut [QueryResult]) {
        for r in results.iter_mut() {
            r.matches.sort_by_key(|m| m.sample_id);
        }
    }

    #[test]
    fn test_chunked_full_range_matches_query_sketch() {
        let (_dir, db_path) = build_test_db_with_params(
            &[
                ("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
                ("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
            ],
            11,
            1,
            true,
        );
        let engine = QueryEngine::open(&db_path).unwrap();
        let reader = JamReader::open(&db_path).unwrap();
        let sketch = QuerySketch::from_jam(&db_path, &reader).unwrap();
        let n = sketch.sample_count();

        let mut expected = engine.query_sketch(&sketch);
        let mut actual = engine.query_sketch_chunked(&sketch, 0..n, 0.0);

        normalize_results(&mut expected);
        normalize_results(&mut actual);

        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.iter().zip(expected.iter()) {
            assert_eq!(a.query_size, e.query_size, "query_size mismatch");
            assert_eq!(a.hashes_found, e.hashes_found, "hashes_found mismatch");
            assert_eq!(a.matches.len(), e.matches.len(), "match count mismatch");
            for (am, em) in a.matches.iter().zip(e.matches.iter()) {
                assert_eq!(am.sample_id, em.sample_id);
                assert_eq!(am.hit_count, em.hit_count);
            }
        }
    }

    #[test]
    fn test_chunked_by_one_sample_matches_query_sketch() {
        let (_dir, db_path) = build_test_db_with_params(
            &[
                ("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
                ("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
            ],
            11,
            1,
            true,
        );
        let engine = QueryEngine::open(&db_path).unwrap();
        let reader = JamReader::open(&db_path).unwrap();
        let sketch = QuerySketch::from_jam(&db_path, &reader).unwrap();
        let n = sketch.sample_count();

        let mut expected = engine.query_sketch(&sketch);
        normalize_results(&mut expected);

        let mut chunked_results: Vec<QueryResult> = Vec::new();
        for i in 0..n {
            let mut chunk = engine.query_sketch_chunked(&sketch, i..i + 1, 0.0);
            normalize_results(&mut chunk);
            chunked_results.extend(chunk);
        }

        assert_eq!(chunked_results.len(), expected.len());
        for (a, e) in chunked_results.iter().zip(expected.iter()) {
            assert_eq!(a.query_size, e.query_size);
            assert_eq!(a.hashes_found, e.hashes_found);
            assert_eq!(a.matches.len(), e.matches.len());
            for (am, em) in a.matches.iter().zip(e.matches.iter()) {
                assert_eq!(am.sample_id, em.sample_id);
                assert_eq!(am.hit_count, em.hit_count);
            }
        }
    }

    #[test]
    fn test_chunked_cutoff_removes_low_containment() {
        let (_dir, db_path) = build_test_db_with_params(
            &[
                ("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
                ("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
            ],
            11,
            1,
            true,
        );
        let engine = QueryEngine::open(&db_path).unwrap();
        let reader = JamReader::open(&db_path).unwrap();
        let sketch = QuerySketch::from_jam(&db_path, &reader).unwrap();
        let n = sketch.sample_count();

        let all_results = engine.query_sketch_chunked(&sketch, 0..n, 0.0);
        let filtered_results = engine.query_sketch_chunked(&sketch, 0..n, 0.9);

        for result in &filtered_results {
            for m in &result.matches {
                assert!(
                    m.containment >= 0.9,
                    "containment {} below cutoff",
                    m.containment
                );
            }
        }

        for (filtered, full) in filtered_results.iter().zip(all_results.iter()) {
            assert!(filtered.matches.len() <= full.matches.len());
        }
    }

    #[test]
    fn test_chunked_two_halves_match_full_range() {
        let (_dir, db_path) = build_test_db_with_params(
            &[
                ("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
                ("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
                ("seq3", "TATATATATATATATATATATATATATATATA"),
                ("seq4", "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"),
            ],
            11,
            1,
            true,
        );
        let engine = QueryEngine::open(&db_path).unwrap();
        let reader = JamReader::open(&db_path).unwrap();
        let sketch = QuerySketch::from_jam(&db_path, &reader).unwrap();
        let n = sketch.sample_count();

        let mut expected = engine.query_sketch(&sketch);
        normalize_results(&mut expected);

        let mut first_half = engine.query_sketch_chunked(&sketch, 0..n / 2, 0.0);
        let mut second_half = engine.query_sketch_chunked(&sketch, n / 2..n, 0.0);
        normalize_results(&mut first_half);
        normalize_results(&mut second_half);

        let mut combined = first_half;
        combined.extend(second_half);

        assert_eq!(combined.len(), expected.len());
        for (a, e) in combined.iter().zip(expected.iter()) {
            assert_eq!(a.query_size, e.query_size);
            assert_eq!(a.hashes_found, e.hashes_found);
            assert_eq!(a.matches.len(), e.matches.len());
            for (am, em) in a.matches.iter().zip(e.matches.iter()) {
                assert_eq!(am.sample_id, em.sample_id);
                assert_eq!(am.hit_count, em.hit_count);
            }
        }
    }

    #[test]
    fn test_from_inputs_integration_with_query_engine() {
        let (_dir1, db_path) = build_test_db_with_params(
            &[("db_seq", "ATCGATCGATCGATCGATCGATCGATCGATCG")],
            11,
            1,
            false,
        );
        let db = JamReader::open(&db_path).unwrap();
        let engine = QueryEngine::open(&db_path).unwrap();

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

        let results = engine.query_sketch(&sketch);

        assert_eq!(results.len(), 2);

        assert!(results[0].has_matches());
        let top0 = results[0].top(1);
        assert!(!top0.is_empty());
        assert!(
            top0[0].containment >= 0.9,
            "Same sequence should have high containment, got {}",
            top0[0].containment
        );
    }

    #[test]
    fn test_query_fasta_non_singleton() {
        let (_dir, db_path) =
            build_test_db(&[("db_seq", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let engine = QueryEngine::open(&db_path).unwrap();

        let query_fasta = make_fasta(&[("query_seq", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);

        let results = engine.query_fasta(query_fasta.path(), false).unwrap();

        assert_eq!(results.len(), 1);
        assert!(results[0].has_matches());

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
        let (_dir, db_path) = build_test_db(
            &[
                ("db_seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
                ("db_seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
            ],
            true,
        );
        let engine = QueryEngine::open(&db_path).unwrap();

        let query_fasta = make_fasta(&[
            ("query1", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
            ("query2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
        ]);

        let results = engine.query_fasta(query_fasta.path(), true).unwrap();

        assert_eq!(results.len(), 2);

        assert!(results[0].has_matches());
        assert!(results[1].has_matches());

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
        let (_dir, db_path) =
            build_test_db(&[("db_seq", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let engine = QueryEngine::open(&db_path).unwrap();

        let result = engine.query_fasta("/nonexistent/path.fasta", false);

        assert!(result.is_err());
        match result.unwrap_err() {
            QueryError::Parse { path, message: _ } => {
                assert!(path.contains("nonexistent"));
            }
            e => panic!("Expected Parse error, got {:?}", e),
        }
    }

    #[test]
    fn test_compute_e_value_no_hits() {
        let e = compute_e_value(0, 100, 200, u64::MAX, 10);
        assert_eq!(e, 10.0);
    }

    #[test]
    fn test_compute_e_value_zero_query_with_hits() {
        let e = compute_e_value(5, 0, 0, u64::MAX, 10);
        assert_eq!(e, 0.0);
    }

    #[test]
    fn test_compute_e_value_significant_hit() {
        let e = compute_e_value(100, 1000, 1000, u64::MAX, 100);
        assert!(e < 1.0, "Expected significant e-value, got {e}");
    }

    #[test]
    fn test_compute_e_value_insignificant_hit() {
        let threshold = 1_000_000u64;
        let e = compute_e_value(1, 1000, 1000, threshold, 100);
        assert!(e > 1.0, "Expected non-significant e-value, got {e}");
    }

    #[test]
    fn test_query_result_has_weight_fields() {
        let (_dir, path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let engine = QueryEngine::open(&path).unwrap();
        let result = engine.query(&[]);
        assert_eq!(result.total_query_weight, 0.0);
    }

    #[test]
    fn test_e_value_present_in_query_results() {
        let (_dir, path) = build_test_db(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")], false);
        let engine = QueryEngine::open(&path).unwrap();

        let reader = JamReader::open(&path).unwrap();
        let mut test_hashes = Vec::new();
        for bucket_idx in 0..256 {
            for entry in reader.bucket_entries(bucket_idx) {
                test_hashes.push(entry.hash);
            }
        }

        if !test_hashes.is_empty() {
            let result = engine.query(&test_hashes);
            assert!(result.has_matches());
            for m in &result.matches {
                assert!(m.e_value >= 0.0, "e_value should be non-negative");
            }
        }
    }
}
