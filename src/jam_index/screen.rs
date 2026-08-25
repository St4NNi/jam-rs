//! Bounded parallel Stage-1 screening across pure `.jam` part shards.

use super::builder::{JamIndexBuildError, load_manifest};
use super::manifest::{JamIndexManifest, JamIndexPart, ScreenSelectionPolicy};
use crate::format::{Entry, bucket_id};
use crate::jamhash_u64_v1;
use crate::reader::{JamReader, ReaderError, lower_bound_hash_counted};
use crate::trace::model::CandidateAdmissionSource;
use needletail::Sequence;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeMap, BTreeSet, BinaryHeap, HashMap, HashSet};
use std::path::Path;
use thiserror::Error;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct QueryHashOccurrence {
    pub packed_kmer: u64,
    pub query_position: u64,
    pub query_reverse: bool,
    pub query_window_id: u32,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct QueryHashEvidence {
    occurrences: Vec<QueryHashOccurrence>,
    windows: Vec<u32>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PreparedJamIndexQuery {
    pub query_length: u64,
    pub query_window_bases: u32,
    pub total_query_windows: u32,
    hashes: BTreeMap<u64, QueryHashEvidence>,
}

impl PreparedJamIndexQuery {
    #[must_use]
    pub fn unique_hash_count(&self) -> usize {
        self.hashes.len()
    }

    pub fn hashes(&self) -> impl Iterator<Item = u64> + '_ {
        self.hashes.keys().copied()
    }

    #[must_use]
    pub fn occurrences(&self, hash: u64) -> Option<&[QueryHashOccurrence]> {
        self.hashes
            .get(&hash)
            .map(|evidence| evidence.occurrences.as_slice())
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SharedScreenHash {
    pub hash: u64,
    pub entry_ordinal: u64,
    pub document_frequency: u32,
    pub occurrences: Vec<QueryHashOccurrence>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct JamIndexCandidate {
    pub part_id: u32,
    pub metagenome_local_id: u32,
    pub metagenome_id: String,
    pub rank: u32,
    pub shared_hash_count: u32,
    pub supported_query_windows: u32,
    pub longest_supported_window_run: u32,
    pub weighted_hash_sum: f64,
    pub shared_spatial_signatures: u32,
    pub rare_shared_signatures: u32,
    pub shared_whole_sample_signatures: u32,
    pub screen_policy: String,
    pub admission_source: CandidateAdmissionSource,
    pub shared_hashes: Vec<SharedScreenHash>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct JamIndexScreenConfig {
    pub top_candidates: usize,
    pub accumulator_capacity: usize,
    pub min_shared_hashes: u32,
    pub min_query_windows: u32,
    pub rare_rescue_max_document_frequency: Option<u32>,
    pub parallel_parts: usize,
}

impl Default for JamIndexScreenConfig {
    fn default() -> Self {
        Self {
            top_candidates: 250,
            accumulator_capacity: 4_096,
            min_shared_hashes: 1,
            min_query_windows: 1,
            rare_rescue_max_document_frequency: None,
            parallel_parts: 1,
        }
    }
}

impl JamIndexScreenConfig {
    fn validate(self) -> Result<Self, JamIndexScreenError> {
        if self.top_candidates == 0
            || self.accumulator_capacity < self.top_candidates
            || self.min_shared_hashes == 0
            || self.min_query_windows == 0
            || self.rare_rescue_max_document_frequency == Some(0)
            || self.parallel_parts == 0
        {
            return Err(JamIndexScreenError::InvalidConfig);
        }
        Ok(self)
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub struct JamIndexScreenMetrics {
    pub parts_searched: u64,
    pub query_hashes: u64,
    pub hash_postings_seen: u64,
    pub maximum_accumulator_entries: u64,
    pub accumulator_evictions: u64,
    pub exact_candidates_recounted: u64,
    pub candidate_limit_hits: u64,
    pub selected_candidates: u64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct JamIndexScreenResult {
    pub candidates: Vec<JamIndexCandidate>,
    pub metrics: JamIndexScreenMetrics,
    pub candidate_limit_reached: bool,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct JamIndexBatchScreenResult {
    pub queries: Vec<JamIndexScreenResult>,
    pub parts_opened: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct CandidateCore {
    part_id: u32,
    sample_id: u32,
    shared_hash_count: u32,
    supported_query_windows: u32,
    longest_supported_window_run: u32,
    weighted_hash_micros: i64,
    admission_source: CandidateAdmissionSource,
}

#[derive(Clone, Copy, Debug)]
struct BatchQueryHashReference {
    query_index: usize,
    window_bits: u64,
}

#[derive(Clone, Debug)]
struct BatchBucketHash {
    hash: u64,
    queries: Vec<BatchQueryHashReference>,
}

#[derive(Clone, Copy, Debug)]
struct ApproxBatchEvidence {
    query_index: usize,
    candidate_key: u64,
    count: u32,
    windows: u64,
}

#[derive(Debug)]
struct FirstBucketResult {
    evidence: Vec<ApproxBatchEvidence>,
    document_frequencies: Vec<(u64, u32)>,
    query_postings: Vec<(usize, u64)>,
    query_evictions: Vec<(usize, u64)>,
    maximum_entries: u64,
    comparisons: u64,
    used_merge: bool,
}

#[derive(Clone, Copy, Debug)]
struct ExactBatchHit {
    query_index: usize,
    candidate_key: u64,
    hash: u64,
    entry_ordinal: u64,
}

#[derive(Debug)]
struct ExactBucketResult {
    hits: Vec<ExactBatchHit>,
    comparisons: u64,
    used_merge: bool,
}

#[derive(Clone, Debug, Default)]
struct ExactBatchCounter {
    hashes: BTreeMap<u64, u64>,
    windows: BTreeSet<u32>,
}

impl CandidateCore {
    fn weighted_hash_sum(self) -> f64 {
        self.weighted_hash_micros as f64 / 1_000_000.0
    }
}

#[derive(Clone, Copy, Debug)]
struct ExactCounter {
    shared_hash_count: u32,
    weighted_hash_micros: i64,
    minimum_document_frequency: u32,
}

impl Default for ExactCounter {
    fn default() -> Self {
        Self {
            shared_hash_count: 0,
            weighted_hash_micros: 0,
            minimum_document_frequency: u32::MAX,
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct ApproxEntry {
    count: u32,
    windows: u64,
    generation: u64,
}

#[derive(Debug)]
struct BoundedAccumulator {
    capacity: usize,
    generation: u64,
    evictions: u64,
    entries: HashMap<u64, ApproxEntry>,
    minimum: BinaryHeap<Reverse<(u32, u32, u64, u64)>>,
}

impl BoundedAccumulator {
    fn new(capacity: usize) -> Self {
        Self {
            capacity,
            generation: 0,
            evictions: 0,
            entries: HashMap::new(),
            minimum: BinaryHeap::new(),
        }
    }

    fn add(&mut self, sample_id: u64, window_bits: u64) {
        self.generation = self.generation.wrapping_add(1);
        if let Some(entry) = self.entries.get_mut(&sample_id) {
            entry.count = entry.count.saturating_add(1);
            entry.windows |= window_bits;
            entry.generation = self.generation;
            self.minimum.push(Reverse((
                entry.count,
                entry.windows.count_ones(),
                entry.generation,
                sample_id,
            )));
            return;
        }
        let base = if self.entries.len() < self.capacity {
            0
        } else {
            self.evictions = self.evictions.saturating_add(1);
            self.pop_minimum().map_or(0, |(_, entry)| entry.count)
        };
        let entry = ApproxEntry {
            count: base.saturating_add(1),
            windows: window_bits,
            generation: self.generation,
        };
        self.entries.insert(sample_id, entry);
        self.minimum.push(Reverse((
            entry.count,
            entry.windows.count_ones(),
            entry.generation,
            sample_id,
        )));
    }

    fn add_count(&mut self, sample_id: u64, window_bits: u64, count: u32) {
        for _ in 0..count {
            self.add(sample_id, window_bits);
        }
    }

    fn pop_minimum(&mut self) -> Option<(u64, ApproxEntry)> {
        while let Some(Reverse((count, windows, generation, sample_id))) = self.minimum.pop() {
            let Some(entry) = self.entries.get(&sample_id).copied() else {
                continue;
            };
            if (entry.count, entry.windows.count_ones(), entry.generation)
                != (count, windows, generation)
            {
                continue;
            }
            self.entries.remove(&sample_id);
            return Some((sample_id, entry));
        }
        None
    }

    fn sample_ids(&self) -> HashSet<u64> {
        self.entries.keys().copied().collect()
    }
}

pub fn prepare_screen_query(
    sequence: &[u8],
    policy: &ScreenSelectionPolicy,
) -> Result<PreparedJamIndexQuery, JamIndexScreenError> {
    policy
        .validate()
        .map_err(|_| JamIndexScreenError::InvalidConfig)?;
    if sequence.is_empty() {
        return Err(JamIndexScreenError::EmptyQuery);
    }
    let query_length = u64::try_from(sequence.len()).map_err(|_| JamIndexScreenError::Overflow)?;
    let total_query_windows = u32::try_from(
        query_length.saturating_add(u64::from(policy.query_window_bases) - 1)
            / u64::from(policy.query_window_bases),
    )
    .map_err(|_| JamIndexScreenError::Overflow)?;
    let normalized = sequence.normalize(false);
    let mut generated = Vec::<(u64, QueryHashOccurrence)>::new();
    {
        let _profile = crate::profiling::scope("query_k21_generation");
        for (position, kmer, reverse) in normalized.bit_kmers(policy.k, true) {
            let hash = jamhash_u64_v1(kmer.0);
            if hash == 0 {
                continue;
            }
            let position = u64::try_from(position).map_err(|_| JamIndexScreenError::Overflow)?;
            let first_window = u32::try_from(position / u64::from(policy.query_window_bases))
                .map_err(|_| JamIndexScreenError::Overflow)?;
            let last_position = position
                .checked_add(u64::from(policy.k - 1))
                .ok_or(JamIndexScreenError::Overflow)?;
            let last_window = u32::try_from(last_position / u64::from(policy.query_window_bases))
                .map_err(|_| JamIndexScreenError::Overflow)?;
            for window in first_window..=last_window {
                generated.push((
                    hash,
                    QueryHashOccurrence {
                        packed_kmer: kmer.0,
                        query_position: position,
                        query_reverse: reverse,
                        query_window_id: window,
                    },
                ));
            }
        }
    }
    let mut hashes = BTreeMap::<u64, QueryHashEvidence>::new();
    {
        let _profile = crate::profiling::scope("query_hash_sort_deduplication");
        generated.sort_unstable_by_key(|(hash, occurrence)| {
            (
                *hash,
                occurrence.query_position,
                occurrence.query_window_id,
                occurrence.packed_kmer,
                occurrence.query_reverse,
            )
        });
        generated.dedup();
        for (hash, occurrence) in generated {
            let evidence = hashes.entry(hash).or_insert_with(|| QueryHashEvidence {
                occurrences: Vec::new(),
                windows: Vec::new(),
            });
            evidence.windows.push(occurrence.query_window_id);
            evidence.occurrences.push(occurrence);
        }
        for evidence in hashes.values_mut() {
            evidence.windows.sort_unstable();
            evidence.windows.dedup();
        }
    }
    Ok(PreparedJamIndexQuery {
        query_length,
        query_window_bases: policy.query_window_bases,
        total_query_windows,
        hashes,
    })
}

pub fn search_jam_index(
    root: impl AsRef<Path>,
    query: &PreparedJamIndexQuery,
    config: JamIndexScreenConfig,
) -> Result<JamIndexScreenResult, JamIndexScreenError> {
    let config = config.validate()?;
    let root = root.as_ref();
    let manifest = load_manifest(root)?;
    if query.query_window_bases != manifest.selection_policy.query_window_bases {
        return Err(JamIndexScreenError::PolicyMismatch);
    }
    let screens = open_screens(root, &manifest)?;
    search_open_index(&manifest, &screens, query, config)
}

pub fn search_jam_index_batch(
    root: impl AsRef<Path>,
    queries: &[PreparedJamIndexQuery],
    config: JamIndexScreenConfig,
) -> Result<JamIndexBatchScreenResult, JamIndexScreenError> {
    let config = config.validate()?;
    let root = root.as_ref();
    let manifest = load_manifest(root)?;
    if queries.is_empty() {
        return Err(JamIndexScreenError::EmptyBatch);
    }
    if queries
        .iter()
        .any(|query| query.query_window_bases != manifest.selection_policy.query_window_bases)
    {
        return Err(JamIndexScreenError::PolicyMismatch);
    }
    let screens = open_screens(root, &manifest)?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.parallel_parts)
        .build()
        .map_err(|_| JamIndexScreenError::InvalidConfig)?;
    let results = {
        let _profile = crate::profiling::process_scope("screen_part_candidate_search");
        search_open_batch(&manifest, &screens, queries, config, &pool)?
    };
    Ok(JamIndexBatchScreenResult {
        queries: results,
        parts_opened: u64::try_from(screens.len()).unwrap_or(u64::MAX),
    })
}

fn open_screens<'a>(
    root: &Path,
    manifest: &'a JamIndexManifest,
) -> Result<Vec<(&'a JamIndexPart, JamReader)>, JamIndexScreenError> {
    let screens = {
        let _profile = crate::profiling::scope("screen_part_open");
        manifest
            .parts
            .iter()
            .map(|part| Ok((part, open_screen(root, part)?)))
            .collect::<Result<Vec<_>, JamIndexScreenError>>()?
    };
    crate::profiling::add_counter(
        "screen_jam_mapped_bytes",
        screens
            .iter()
            .map(|(_, reader)| reader.stats().file_size)
            .fold(0u64, u64::saturating_add),
    );
    Ok(screens)
}

fn search_open_batch(
    manifest: &JamIndexManifest,
    screens: &[(&JamIndexPart, JamReader)],
    queries: &[PreparedJamIndexQuery],
    config: JamIndexScreenConfig,
    pool: &rayon::ThreadPool,
) -> Result<Vec<JamIndexScreenResult>, JamIndexScreenError> {
    let buckets = batch_query_buckets(queries);
    let work = screens
        .iter()
        .enumerate()
        .flat_map(|(part_index, _)| {
            buckets
                .iter()
                .enumerate()
                .filter(|(_, hashes)| !hashes.is_empty())
                .map(move |(bucket_index, _)| (part_index, bucket_index))
        })
        .collect::<Vec<_>>();
    crate::profiling::add_counter(
        "database_buckets_processed",
        u64::try_from(work.len()).unwrap_or(u64::MAX),
    );
    crate::profiling::add_counter(
        "database_bucket_passes",
        u64::try_from(work.len())
            .unwrap_or(u64::MAX)
            .saturating_mul(2),
    );
    crate::profiling::add_counter(
        "query_hashes_processed",
        queries
            .iter()
            .map(|query| u64::try_from(query.unique_hash_count()).unwrap_or(u64::MAX))
            .fold(0u64, u64::saturating_add),
    );
    let mut accumulators = (0..queries.len())
        .map(|_| BoundedAccumulator::new(config.accumulator_capacity))
        .collect::<Vec<_>>();
    let mut metrics = vec![JamIndexScreenMetrics::default(); queries.len()];
    let mut document_frequencies = HashMap::<u64, u32>::new();
    let wave_size = config.parallel_parts.saturating_mul(2).max(1);
    for wave in work.chunks(wave_size) {
        let results = pool.install(|| {
            wave.par_iter()
                .map(|(part_index, bucket_index)| {
                    let _worker = crate::profiling::worker("screen_part_candidate_search");
                    let (part, reader) = screens
                        .get(*part_index)
                        .ok_or(JamIndexScreenError::Overflow)?;
                    let _profile = crate::profiling::scope("bucket_search");
                    search_first_bucket(
                        part.part_id,
                        reader,
                        *bucket_index,
                        &buckets[*bucket_index],
                        config.accumulator_capacity,
                    )
                })
                .collect::<Result<Vec<_>, JamIndexScreenError>>()
        })?;
        let _accumulation_profile = crate::profiling::scope("candidate_accumulation");
        for result in results {
            record_bucket_method(result.used_merge);
            record_signature_comparisons(result.comparisons);
            for (hash, frequency) in result.document_frequencies {
                let total = document_frequencies.entry(hash).or_default();
                *total = total.saturating_add(frequency);
            }
            for (query_index, postings) in result.query_postings {
                metrics[query_index].hash_postings_seen = metrics[query_index]
                    .hash_postings_seen
                    .saturating_add(postings);
            }
            for (query_index, evictions) in result.query_evictions {
                metrics[query_index].accumulator_evictions = metrics[query_index]
                    .accumulator_evictions
                    .saturating_add(evictions);
            }
            for evidence in result.evidence {
                accumulators[evidence.query_index].add_count(
                    evidence.candidate_key,
                    evidence.windows,
                    evidence.count,
                );
            }
            for query_metrics in &mut metrics {
                query_metrics.maximum_accumulator_entries = query_metrics
                    .maximum_accumulator_entries
                    .max(result.maximum_entries);
            }
        }
        drop(_accumulation_profile);
    }
    crate::profiling::add_counter(
        "candidate_maps_created",
        u64::try_from(queries.len()).unwrap_or(u64::MAX),
    );
    let selected = accumulators
        .iter()
        .map(BoundedAccumulator::sample_ids)
        .collect::<Vec<_>>();
    for (query_index, accumulator) in accumulators.iter().enumerate() {
        metrics[query_index].maximum_accumulator_entries = metrics[query_index]
            .maximum_accumulator_entries
            .max(u64::try_from(accumulator.entries.len()).unwrap_or(u64::MAX));
        metrics[query_index].accumulator_evictions = metrics[query_index]
            .accumulator_evictions
            .saturating_add(accumulator.evictions);
        metrics[query_index].exact_candidates_recounted =
            u64::try_from(selected[query_index].len()).unwrap_or(u64::MAX);
        crate::profiling::add_counter(
            "candidate_map_entries",
            u64::try_from(accumulator.entries.len()).unwrap_or(u64::MAX),
        );
    }

    let mut exact = vec![BTreeMap::<u64, ExactBatchCounter>::new(); queries.len()];
    {
        let _profile = crate::profiling::process_scope("exact_candidate_recount");
        for wave in work.chunks(wave_size) {
            let results = pool.install(|| {
                wave.par_iter()
                    .map(|(part_index, bucket_index)| {
                        let _worker = crate::profiling::worker("exact_candidate_recount");
                        let (part, reader) = screens
                            .get(*part_index)
                            .ok_or(JamIndexScreenError::Overflow)?;
                        search_exact_bucket(
                            part.part_id,
                            reader,
                            *bucket_index,
                            &buckets[*bucket_index],
                            &selected,
                        )
                    })
                    .collect::<Result<Vec<_>, JamIndexScreenError>>()
            })?;
            for result in results {
                record_bucket_method(result.used_merge);
                record_signature_comparisons(result.comparisons);
                for hit in result.hits {
                    let query = queries
                        .get(hit.query_index)
                        .ok_or(JamIndexScreenError::Overflow)?;
                    let evidence = query
                        .hashes
                        .get(&hit.hash)
                        .ok_or(JamIndexScreenError::Overflow)?;
                    let counter = exact[hit.query_index].entry(hit.candidate_key).or_default();
                    counter.hashes.insert(hit.hash, hit.entry_ordinal);
                    counter.windows.extend(evidence.windows.iter().copied());
                    metrics[hit.query_index].hash_postings_seen = metrics[hit.query_index]
                        .hash_postings_seen
                        .saturating_add(1);
                }
            }
        }
    }

    let total_metagenomes =
        u32::try_from(manifest.total_metagenomes).map_err(|_| JamIndexScreenError::Overflow)?;
    let mut results = Vec::with_capacity(queries.len());
    for (query_index, query) in queries.iter().enumerate() {
        let _profile = crate::profiling::scope("candidate_merge");
        let mut cores = exact[query_index]
            .iter()
            .filter_map(|(candidate_key, counter)| {
                let shared_hash_count = u32::try_from(counter.hashes.len()).unwrap_or(u32::MAX);
                let supported_query_windows =
                    u32::try_from(counter.windows.len()).unwrap_or(u32::MAX);
                let minimum_document_frequency = counter
                    .hashes
                    .keys()
                    .filter_map(|hash| document_frequencies.get(hash).copied())
                    .min()
                    .unwrap_or(u32::MAX);
                let weighted_hash_micros = counter.hashes.keys().fold(0i64, |total, hash| {
                    total.saturating_add(rarity_weight_micros(
                        total_metagenomes,
                        document_frequencies.get(hash).copied().unwrap_or(0),
                    ))
                });
                let exact_counter = ExactCounter {
                    shared_hash_count,
                    weighted_hash_micros,
                    minimum_document_frequency,
                };
                candidate_admission(exact_counter, supported_query_windows, config).map(
                    |admission_source| CandidateCore {
                        part_id: candidate_part(*candidate_key),
                        sample_id: candidate_sample(*candidate_key),
                        shared_hash_count,
                        supported_query_windows,
                        longest_supported_window_run: longest_run(&counter.windows),
                        weighted_hash_micros,
                        admission_source,
                    },
                )
            })
            .collect::<Vec<_>>();
        cores.sort_by(core_cmp);
        let candidate_limit_hit = cores.len() > config.top_candidates;
        cores.truncate(config.top_candidates);
        let mut candidates = cores
            .into_iter()
            .map(|core| {
                let (_, reader) = screens
                    .iter()
                    .find(|(part, _)| part.part_id == core.part_id)
                    .ok_or(JamIndexScreenError::UnknownPart(core.part_id))?;
                let metagenome_id = reader.sample_name(core.sample_id).ok_or(
                    JamIndexScreenError::UnknownSample {
                        part_id: core.part_id,
                        sample_id: core.sample_id,
                    },
                )?;
                let counter = exact[query_index]
                    .get(&candidate_key(core.part_id, core.sample_id))
                    .ok_or(JamIndexScreenError::Overflow)?;
                let shared_hashes = counter
                    .hashes
                    .iter()
                    .map(|(hash, entry_ordinal)| {
                        let evidence = query
                            .hashes
                            .get(hash)
                            .ok_or(JamIndexScreenError::Overflow)?;
                        Ok(SharedScreenHash {
                            hash: *hash,
                            entry_ordinal: *entry_ordinal,
                            document_frequency: document_frequencies
                                .get(hash)
                                .copied()
                                .unwrap_or(0),
                            occurrences: evidence.occurrences.clone(),
                        })
                    })
                    .collect::<Result<Vec<_>, JamIndexScreenError>>()?;
                let rare_shared_signatures =
                    rare_shared_count(&shared_hashes, config.rare_rescue_max_document_frequency);
                Ok(JamIndexCandidate {
                    part_id: core.part_id,
                    metagenome_local_id: core.sample_id,
                    metagenome_id: metagenome_id.to_string(),
                    rank: 0,
                    shared_hash_count: core.shared_hash_count,
                    supported_query_windows: core.supported_query_windows,
                    longest_supported_window_run: core.longest_supported_window_run,
                    weighted_hash_sum: core.weighted_hash_sum(),
                    shared_spatial_signatures: 0,
                    rare_shared_signatures,
                    shared_whole_sample_signatures: 0,
                    screen_policy: manifest.selection_policy.policy_id.clone(),
                    admission_source: core.admission_source,
                    shared_hashes,
                })
            })
            .collect::<Result<Vec<_>, JamIndexScreenError>>()?;
        candidates.sort_by(candidate_cmp);
        for (rank, candidate) in candidates.iter_mut().enumerate() {
            candidate.rank = u32::try_from(rank + 1).map_err(|_| JamIndexScreenError::Overflow)?;
        }
        metrics[query_index].parts_searched = u64::try_from(screens.len()).unwrap_or(u64::MAX);
        metrics[query_index].query_hashes = u64::try_from(query.hashes.len()).unwrap_or(u64::MAX);
        metrics[query_index].selected_candidates =
            u64::try_from(candidates.len()).unwrap_or(u64::MAX);
        metrics[query_index].candidate_limit_hits = u64::from(candidate_limit_hit);
        results.push(JamIndexScreenResult {
            candidates,
            metrics: metrics[query_index],
            candidate_limit_reached: candidate_limit_hit
                || metrics[query_index].accumulator_evictions != 0,
        });
    }
    Ok(results)
}

fn batch_query_buckets(queries: &[PreparedJamIndexQuery]) -> Vec<Vec<BatchBucketHash>> {
    let mut buckets = (0..crate::format::BUCKET_COUNT)
        .map(|_| BTreeMap::<u64, Vec<BatchQueryHashReference>>::new())
        .collect::<Vec<_>>();
    for (query_index, query) in queries.iter().enumerate() {
        for (hash, evidence) in &query.hashes {
            buckets[bucket_id(*hash)]
                .entry(*hash)
                .or_default()
                .push(BatchQueryHashReference {
                    query_index,
                    window_bits: window_bits(&evidence.windows),
                });
        }
    }
    buckets
        .into_iter()
        .map(|bucket| {
            bucket
                .into_iter()
                .map(|(hash, queries)| BatchBucketHash { hash, queries })
                .collect()
        })
        .collect()
}

fn search_first_bucket(
    part_id: u32,
    reader: &JamReader,
    bucket_index: usize,
    hashes: &[BatchBucketHash],
    accumulator_capacity: usize,
) -> Result<FirstBucketResult, JamIndexScreenError> {
    let entries = reader.bucket_entries(bucket_index);
    let (intersections, comparisons, used_merge) = bucket_intersections(entries, hashes);
    let mut accumulators = BTreeMap::<usize, BoundedAccumulator>::new();
    let mut query_postings = BTreeMap::<usize, u64>::new();
    let mut document_frequencies = Vec::with_capacity(intersections.len());
    for (hash_index, start, end) in intersections {
        let hash = hashes
            .get(hash_index)
            .ok_or(JamIndexScreenError::Overflow)?;
        let frequency =
            u32::try_from(end.saturating_sub(start)).map_err(|_| JamIndexScreenError::Overflow)?;
        document_frequencies.push((hash.hash, frequency));
        for query in &hash.queries {
            let accumulator = accumulators
                .entry(query.query_index)
                .or_insert_with(|| BoundedAccumulator::new(accumulator_capacity));
            for entry in &entries[start..end] {
                accumulator.add(candidate_key(part_id, entry.sample_id), query.window_bits);
            }
            let postings = query_postings.entry(query.query_index).or_default();
            *postings = postings.saturating_add(u64::from(frequency));
        }
    }
    let maximum_entries = accumulators
        .values()
        .map(|accumulator| u64::try_from(accumulator.entries.len()).unwrap_or(u64::MAX))
        .max()
        .unwrap_or(0);
    let mut evidence = accumulators
        .iter()
        .flat_map(|(query_index, accumulator)| {
            accumulator
                .entries
                .iter()
                .map(|(candidate_key, entry)| ApproxBatchEvidence {
                    query_index: *query_index,
                    candidate_key: *candidate_key,
                    count: entry.count,
                    windows: entry.windows,
                })
        })
        .collect::<Vec<_>>();
    evidence.sort_by_key(|item| (item.query_index, item.candidate_key));
    let query_evictions = accumulators
        .iter()
        .map(|(query_index, accumulator)| (*query_index, accumulator.evictions))
        .collect();
    Ok(FirstBucketResult {
        evidence,
        document_frequencies,
        query_postings: query_postings.into_iter().collect(),
        query_evictions,
        maximum_entries,
        comparisons,
        used_merge,
    })
}

fn search_exact_bucket(
    part_id: u32,
    reader: &JamReader,
    bucket_index: usize,
    hashes: &[BatchBucketHash],
    selected: &[HashSet<u64>],
) -> Result<ExactBucketResult, JamIndexScreenError> {
    let entries = reader.bucket_entries(bucket_index);
    let (intersections, comparisons, used_merge) = bucket_intersections(entries, hashes);
    let ordinal_start = reader.bucket_entry_ordinal_start(bucket_index);
    let mut hits = Vec::new();
    for (hash_index, start, end) in intersections {
        let hash = hashes
            .get(hash_index)
            .ok_or(JamIndexScreenError::Overflow)?;
        for query in &hash.queries {
            let selected = selected
                .get(query.query_index)
                .ok_or(JamIndexScreenError::Overflow)?;
            for (offset, entry) in entries[start..end].iter().enumerate() {
                let candidate_key = candidate_key(part_id, entry.sample_id);
                if selected.contains(&candidate_key) {
                    hits.push(ExactBatchHit {
                        query_index: query.query_index,
                        candidate_key,
                        hash: hash.hash,
                        entry_ordinal: ordinal_start
                            .saturating_add(u64::try_from(start).unwrap_or(u64::MAX))
                            .saturating_add(u64::try_from(offset).unwrap_or(u64::MAX)),
                    });
                }
            }
        }
    }
    Ok(ExactBucketResult {
        hits,
        comparisons,
        used_merge,
    })
}

fn bucket_intersections(
    entries: &[Entry],
    hashes: &[BatchBucketHash],
) -> (Vec<(usize, usize, usize)>, u64, bool) {
    // One lower bound costs approximately log2(database entries); an ordered
    // merge costs one pass over both sides. The rule is deterministic and is
    // evaluated once for the complete multi-query bucket, never per query.
    let binary_steps = usize::BITS
        .saturating_sub(entries.len().max(1).leading_zeros())
        .saturating_add(1) as usize;
    let binary_cost = hashes.len().saturating_mul(binary_steps);
    let merge_cost = hashes.len().saturating_add(entries.len());
    let used_merge = merge_cost < binary_cost;
    let mut intersections = Vec::new();
    let mut comparisons = 0u64;
    if used_merge {
        let mut entry_index = 0usize;
        let mut hash_index = 0usize;
        while entry_index < entries.len() && hash_index < hashes.len() {
            comparisons = comparisons.saturating_add(1);
            let entry_hash = entries[entry_index].hash;
            match entry_hash.cmp(&hashes[hash_index].hash) {
                Ordering::Less => entry_index += 1,
                Ordering::Greater => hash_index += 1,
                Ordering::Equal => {
                    let start = entry_index;
                    while entry_index < entries.len()
                        && entries[entry_index].hash == hashes[hash_index].hash
                    {
                        comparisons = comparisons.saturating_add(1);
                        entry_index += 1;
                    }
                    intersections.push((hash_index, start, entry_index));
                    hash_index += 1;
                }
            }
        }
    } else {
        for (hash_index, hash) in hashes.iter().enumerate() {
            let (start, probes) = lower_bound_hash_counted(entries, hash.hash);
            comparisons = comparisons.saturating_add(probes);
            let mut end = start;
            while end < entries.len() {
                comparisons = comparisons.saturating_add(1);
                if entries[end].hash != hash.hash {
                    break;
                }
                end += 1;
            }
            if end != start {
                intersections.push((hash_index, start, end));
            }
        }
    }
    (intersections, comparisons, used_merge)
}

fn candidate_key(part_id: u32, sample_id: u32) -> u64 {
    (u64::from(part_id) << 32) | u64::from(sample_id)
}

fn candidate_part(candidate_key: u64) -> u32 {
    u32::try_from(candidate_key >> 32).unwrap_or(u32::MAX)
}

fn candidate_sample(candidate_key: u64) -> u32 {
    u32::try_from(candidate_key & u64::from(u32::MAX)).unwrap_or(u32::MAX)
}

fn record_bucket_method(used_merge: bool) {
    crate::profiling::add_counter(
        if used_merge {
            "screen_merge_buckets"
        } else {
            "screen_binary_buckets"
        },
        1,
    );
}

fn search_open_index(
    manifest: &JamIndexManifest,
    screens: &[(&JamIndexPart, JamReader)],
    query: &PreparedJamIndexQuery,
    config: JamIndexScreenConfig,
) -> Result<JamIndexScreenResult, JamIndexScreenError> {
    let mut document_frequencies = BTreeMap::new();
    let mut df_postings = 0u64;
    for (_, reader) in screens {
        let (part_frequencies, part_postings) = part_document_frequencies(reader, query)?;
        for (hash, count) in part_frequencies {
            *document_frequencies.entry(hash).or_insert(0u32) = document_frequencies
                .get(&hash)
                .copied()
                .unwrap_or(0)
                .saturating_add(count);
        }
        df_postings = df_postings.saturating_add(part_postings);
    }
    let total_metagenomes =
        u32::try_from(manifest.total_metagenomes).map_err(|_| JamIndexScreenError::Overflow)?;
    let mut selected = Vec::new();
    let mut rank_metrics = JamIndexScreenMetrics::default();
    for (part, reader) in screens {
        let (part_candidates, part_metrics) = rank_part(
            part,
            reader,
            query,
            &document_frequencies,
            total_metagenomes,
            config,
        )?;
        selected.extend(part_candidates);
        rank_metrics = add_metrics(rank_metrics, part_metrics);
    }
    let candidates = {
        let _profile = crate::profiling::scope("candidate_merge");
        selected.sort_by(core_cmp);
        if selected.len() > config.top_candidates {
            rank_metrics.candidate_limit_hits = rank_metrics.candidate_limit_hits.saturating_add(1);
            selected.truncate(config.top_candidates);
        }
        let mut candidates = materialize_selected(
            screens,
            query,
            &document_frequencies,
            selected,
            &manifest.selection_policy.policy_id,
            config.rare_rescue_max_document_frequency,
        )?;
        candidates.sort_by(candidate_cmp);
        for (index, candidate) in candidates.iter_mut().enumerate() {
            candidate.rank = u32::try_from(index + 1).map_err(|_| JamIndexScreenError::Overflow)?;
        }
        candidates
    };
    Ok(JamIndexScreenResult {
        candidate_limit_reached: rank_metrics.candidate_limit_hits != 0
            || rank_metrics.accumulator_evictions != 0,
        metrics: JamIndexScreenMetrics {
            parts_searched: u64::try_from(manifest.parts.len()).unwrap_or(u64::MAX),
            query_hashes: u64::try_from(query.hashes.len()).unwrap_or(u64::MAX),
            hash_postings_seen: df_postings.saturating_add(rank_metrics.hash_postings_seen),
            selected_candidates: u64::try_from(candidates.len()).unwrap_or(u64::MAX),
            ..rank_metrics
        },
        candidates,
    })
}

fn part_document_frequencies(
    reader: &JamReader,
    query: &PreparedJamIndexQuery,
) -> Result<(BTreeMap<u64, u32>, u64), JamIndexScreenError> {
    let mut frequencies = BTreeMap::new();
    let mut postings = 0u64;
    let mut comparisons = 0u64;
    for hash in query.hashes.keys().copied() {
        let mut hits = reader.search_entries_profiled(hash);
        let count = hits.by_ref().count();
        comparisons = comparisons.saturating_add(hits.comparisons());
        if count != 0 {
            frequencies.insert(
                hash,
                u32::try_from(count).map_err(|_| JamIndexScreenError::Overflow)?,
            );
            postings = postings.saturating_add(u64::try_from(count).unwrap_or(u64::MAX));
        }
    }
    record_signature_comparisons(comparisons);
    Ok((frequencies, postings))
}

fn rank_part(
    part: &JamIndexPart,
    reader: &JamReader,
    query: &PreparedJamIndexQuery,
    document_frequencies: &BTreeMap<u64, u32>,
    total_metagenomes: u32,
    config: JamIndexScreenConfig,
) -> Result<(Vec<CandidateCore>, JamIndexScreenMetrics), JamIndexScreenError> {
    let mut accumulator = BoundedAccumulator::new(config.accumulator_capacity);
    crate::profiling::add_counter("candidate_maps_created", 1);
    let mut postings = 0u64;
    let mut comparisons = 0u64;
    for (hash, evidence) in &query.hashes {
        let bits = window_bits(&evidence.windows);
        let mut hits = reader.search_entries_profiled(*hash);
        for hit in hits.by_ref() {
            accumulator.add(u64::from(hit.sample_id), bits);
            postings = postings.saturating_add(1);
        }
        comparisons = comparisons.saturating_add(hits.comparisons());
    }
    record_signature_comparisons(comparisons);
    crate::profiling::add_counter(
        "candidate_map_entries",
        u64::try_from(accumulator.entries.len()).unwrap_or(u64::MAX),
    );
    let shortlist = accumulator.sample_ids();
    let exact = {
        let _profile = crate::profiling::scope("exact_candidate_recount");
        let mut exact = BTreeMap::<u32, (ExactCounter, BTreeSet<u32>)>::new();
        let mut comparisons = 0u64;
        for (hash, evidence) in &query.hashes {
            let document_frequency = document_frequencies.get(hash).copied().unwrap_or(0);
            let weight = rarity_weight_micros(total_metagenomes, document_frequency);
            let mut hits = reader.search_entries_profiled(*hash);
            for hit in hits.by_ref() {
                let sample_id = hit.sample_id;
                if !shortlist.contains(&u64::from(sample_id)) {
                    continue;
                }
                let entry = exact.entry(sample_id).or_default();
                entry.0.shared_hash_count = entry.0.shared_hash_count.saturating_add(1);
                entry.0.weighted_hash_micros = entry.0.weighted_hash_micros.saturating_add(weight);
                entry.0.minimum_document_frequency =
                    entry.0.minimum_document_frequency.min(document_frequency);
                entry.1.extend(evidence.windows.iter().copied());
            }
            comparisons = comparisons.saturating_add(hits.comparisons());
        }
        record_signature_comparisons(comparisons);
        exact
    };
    let mut candidates = exact
        .into_iter()
        .filter_map(|(sample_id, (counter, windows))| {
            let supported = u32::try_from(windows.len()).unwrap_or(u32::MAX);
            let admission_source = candidate_admission(counter, supported, config);
            admission_source.map(|admission_source| CandidateCore {
                part_id: part.part_id,
                sample_id,
                shared_hash_count: counter.shared_hash_count,
                supported_query_windows: supported,
                longest_supported_window_run: longest_run(&windows),
                weighted_hash_micros: counter.weighted_hash_micros,
                admission_source,
            })
        })
        .collect::<Vec<_>>();
    candidates.sort_by(core_cmp);
    let candidate_limit_hit = candidates.len() > config.top_candidates;
    candidates.truncate(config.top_candidates);
    Ok((
        candidates,
        JamIndexScreenMetrics {
            hash_postings_seen: postings,
            maximum_accumulator_entries: u64::try_from(accumulator.entries.len())
                .unwrap_or(u64::MAX),
            accumulator_evictions: accumulator.evictions,
            exact_candidates_recounted: u64::try_from(shortlist.len()).unwrap_or(u64::MAX),
            candidate_limit_hits: u64::from(candidate_limit_hit),
            ..JamIndexScreenMetrics::default()
        },
    ))
}

fn candidate_admission(
    counter: ExactCounter,
    supported_query_windows: u32,
    config: JamIndexScreenConfig,
) -> Option<CandidateAdmissionSource> {
    if counter.shared_hash_count >= config.min_shared_hashes {
        Some(CandidateAdmissionSource::Standard)
    } else if supported_query_windows >= config.min_query_windows {
        Some(CandidateAdmissionSource::WindowSpread)
    } else if counter.shared_hash_count == 1
        && config
            .rare_rescue_max_document_frequency
            .is_some_and(|threshold| counter.minimum_document_frequency <= threshold)
    {
        Some(CandidateAdmissionSource::RareRescue)
    } else {
        None
    }
}

fn materialize_selected(
    screens: &[(&JamIndexPart, JamReader)],
    query: &PreparedJamIndexQuery,
    document_frequencies: &BTreeMap<u64, u32>,
    selected: Vec<CandidateCore>,
    screen_policy: &str,
    rare_rescue_max_document_frequency: Option<u32>,
) -> Result<Vec<JamIndexCandidate>, JamIndexScreenError> {
    let mut grouped = BTreeMap::<u32, Vec<CandidateCore>>::new();
    for candidate in selected {
        grouped
            .entry(candidate.part_id)
            .or_default()
            .push(candidate);
    }
    let mut nested = Vec::new();
    for (part_id, candidates) in grouped {
        let (part, screen) = screens
            .iter()
            .find(|(part, _)| part.part_id == part_id)
            .ok_or(JamIndexScreenError::UnknownPart(part_id))?;
        nested.push(materialize_part(
            part,
            screen,
            query,
            document_frequencies,
            &candidates,
            screen_policy,
            rare_rescue_max_document_frequency,
        )?);
    }
    Ok(nested.into_iter().flatten().collect())
}

fn materialize_part(
    part: &JamIndexPart,
    screen: &JamReader,
    query: &PreparedJamIndexQuery,
    document_frequencies: &BTreeMap<u64, u32>,
    selected: &[CandidateCore],
    screen_policy: &str,
    rare_rescue_max_document_frequency: Option<u32>,
) -> Result<Vec<JamIndexCandidate>, JamIndexScreenError> {
    let mut by_sample = selected
        .iter()
        .enumerate()
        .map(|(index, candidate)| (candidate.sample_id, index))
        .collect::<HashMap<_, _>>();
    let mut details = vec![Vec::<SharedScreenHash>::new(); selected.len()];
    let mut comparisons = 0u64;
    for (hash, evidence) in &query.hashes {
        let mut hits = screen.search_entries_profiled(*hash);
        for hit in hits.by_ref() {
            let sample_id = hit.sample_id;
            let Some(index) = by_sample.get(&sample_id).copied() else {
                continue;
            };
            details[index].push(SharedScreenHash {
                hash: *hash,
                entry_ordinal: hit.ordinal,
                document_frequency: document_frequencies.get(hash).copied().unwrap_or(0),
                occurrences: evidence.occurrences.clone(),
            });
        }
        comparisons = comparisons.saturating_add(hits.comparisons());
    }
    record_signature_comparisons(comparisons);
    by_sample.clear();
    selected
        .iter()
        .copied()
        .zip(details)
        .map(|(candidate, mut shared_hashes)| {
            shared_hashes.sort_by_key(|shared| shared.hash);
            let rare_shared_signatures =
                rare_shared_count(&shared_hashes, rare_rescue_max_document_frequency);
            let metagenome_id = screen.sample_name(candidate.sample_id).ok_or(
                JamIndexScreenError::UnknownSample {
                    part_id: part.part_id,
                    sample_id: candidate.sample_id,
                },
            )?;
            Ok(JamIndexCandidate {
                part_id: candidate.part_id,
                metagenome_local_id: candidate.sample_id,
                metagenome_id: metagenome_id.to_string(),
                rank: 0,
                shared_hash_count: candidate.shared_hash_count,
                supported_query_windows: candidate.supported_query_windows,
                longest_supported_window_run: candidate.longest_supported_window_run,
                weighted_hash_sum: candidate.weighted_hash_sum(),
                shared_spatial_signatures: 0,
                rare_shared_signatures,
                shared_whole_sample_signatures: 0,
                screen_policy: screen_policy.to_string(),
                admission_source: candidate.admission_source,
                shared_hashes,
            })
        })
        .collect()
}

fn rare_shared_count(shared: &[SharedScreenHash], threshold: Option<u32>) -> u32 {
    u32::try_from(
        shared
            .iter()
            .filter(|item| threshold.is_some_and(|limit| item.document_frequency <= limit))
            .count(),
    )
    .unwrap_or(u32::MAX)
}

fn record_signature_comparisons(comparisons: u64) {
    crate::profiling::add_counter("signature_comparisons", comparisons);
    crate::profiling::add_counter(
        "screen_jam_bytes_read",
        comparisons.saturating_mul(crate::format::ENTRY_SIZE as u64),
    );
}

fn open_screen(root: &Path, part: &JamIndexPart) -> Result<JamReader, JamIndexScreenError> {
    let reader = JamReader::open(root.join(&part.directory).join(&part.screen_file))?;
    if reader.kmer_size() != 21
        || reader.threshold() != u64::MAX
        || u32::try_from(reader.sample_names().len()).ok() != Some(part.metagenome_count)
    {
        return Err(JamIndexScreenError::PartBinding(part.part_id));
    }
    Ok(reader)
}

fn core_cmp(left: &CandidateCore, right: &CandidateCore) -> Ordering {
    right
        .weighted_hash_micros
        .cmp(&left.weighted_hash_micros)
        .then_with(|| {
            right
                .supported_query_windows
                .cmp(&left.supported_query_windows)
        })
        .then_with(|| {
            right
                .longest_supported_window_run
                .cmp(&left.longest_supported_window_run)
        })
        .then_with(|| right.shared_hash_count.cmp(&left.shared_hash_count))
        .then_with(|| left.part_id.cmp(&right.part_id))
        .then_with(|| left.sample_id.cmp(&right.sample_id))
}

fn candidate_cmp(left: &JamIndexCandidate, right: &JamIndexCandidate) -> Ordering {
    right
        .weighted_hash_sum
        .total_cmp(&left.weighted_hash_sum)
        .then_with(|| {
            right
                .supported_query_windows
                .cmp(&left.supported_query_windows)
        })
        .then_with(|| {
            right
                .longest_supported_window_run
                .cmp(&left.longest_supported_window_run)
        })
        .then_with(|| right.shared_hash_count.cmp(&left.shared_hash_count))
        .then_with(|| left.part_id.cmp(&right.part_id))
        .then_with(|| left.metagenome_local_id.cmp(&right.metagenome_local_id))
}

fn rarity_weight_micros(total_metagenomes: u32, document_frequency: u32) -> i64 {
    (((f64::from(total_metagenomes) + 1.0) / (f64::from(document_frequency) + 1.0)).ln()
        * 1_000_000.0)
        .round() as i64
}

fn window_bits(windows: &[u32]) -> u64 {
    windows
        .iter()
        .fold(0u64, |bits, window| bits | (1u64 << (window % 64)))
}

fn longest_run(windows: &BTreeSet<u32>) -> u32 {
    let mut longest = 0u32;
    let mut current = 0u32;
    let mut previous = None;
    for window in windows {
        current = if previous.is_some_and(|old| old + 1 == *window) {
            current.saturating_add(1)
        } else {
            1
        };
        longest = longest.max(current);
        previous = Some(*window);
    }
    longest
}

fn add_metrics(left: JamIndexScreenMetrics, right: JamIndexScreenMetrics) -> JamIndexScreenMetrics {
    JamIndexScreenMetrics {
        parts_searched: left.parts_searched.saturating_add(right.parts_searched),
        query_hashes: left.query_hashes.saturating_add(right.query_hashes),
        hash_postings_seen: left
            .hash_postings_seen
            .saturating_add(right.hash_postings_seen),
        maximum_accumulator_entries: left
            .maximum_accumulator_entries
            .max(right.maximum_accumulator_entries),
        accumulator_evictions: left
            .accumulator_evictions
            .saturating_add(right.accumulator_evictions),
        exact_candidates_recounted: left
            .exact_candidates_recounted
            .saturating_add(right.exact_candidates_recounted),
        candidate_limit_hits: left
            .candidate_limit_hits
            .saturating_add(right.candidate_limit_hits),
        selected_candidates: left
            .selected_candidates
            .saturating_add(right.selected_candidates),
    }
}

#[derive(Debug, Error)]
pub enum JamIndexScreenError {
    #[error("invalid Jam Index screen config")]
    InvalidConfig,
    #[error("Jam Index query is empty")]
    EmptyQuery,
    #[error("Jam Index query batch is empty")]
    EmptyBatch,
    #[error("Jam Index query policy does not match the manifest")]
    PolicyMismatch,
    #[error("Jam Index screen coordinate overflow")]
    Overflow,
    #[error("unknown Jam Index part {0}")]
    UnknownPart(u32),
    #[error("unknown sample {sample_id} in Jam Index part {part_id}")]
    UnknownSample { part_id: u32, sample_id: u32 },
    #[error("Jam Index part {0} screen/data binding mismatch")]
    PartBinding(u32),
    #[error(transparent)]
    Build(#[from] JamIndexBuildError),
    #[error(transparent)]
    Screen(#[from] ReaderError),
    #[error(transparent)]
    Part(#[from] super::part::JamIndexPartError),
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counter(shared: u32, document_frequency: u32) -> ExactCounter {
        ExactCounter {
            shared_hash_count: shared,
            minimum_document_frequency: document_frequency,
            ..ExactCounter::default()
        }
    }

    fn config(rare_rescue: Option<u32>) -> JamIndexScreenConfig {
        JamIndexScreenConfig {
            min_shared_hashes: 2,
            min_query_windows: 2,
            rare_rescue_max_document_frequency: rare_rescue,
            ..JamIndexScreenConfig::default()
        }
    }

    #[test]
    fn candidate_admission_uses_standard_window_and_rare_rules() {
        assert_eq!(
            candidate_admission(counter(2, 20), 1, config(Some(4))),
            Some(CandidateAdmissionSource::Standard)
        );
        assert_eq!(
            candidate_admission(counter(1, 20), 2, config(Some(4))),
            Some(CandidateAdmissionSource::WindowSpread)
        );
        assert_eq!(
            candidate_admission(counter(1, 4), 1, config(Some(4))),
            Some(CandidateAdmissionSource::RareRescue)
        );
        assert_eq!(candidate_admission(counter(1, 5), 1, config(Some(4))), None);
        assert_eq!(candidate_admission(counter(1, 1), 1, config(None)), None);
    }

    #[test]
    fn query_preparation_retains_every_unique_nonzero_canonical_k21() {
        let sequence = (0..1_024)
            .map(|index| b"ACGTTGCA"[(index * 5 + index / 7) % 8])
            .collect::<Vec<_>>();
        let normalized = sequence.normalize(false);
        let expected = normalized
            .bit_kmers(21, true)
            .map(|(_, kmer, _)| jamhash_u64_v1(kmer.0))
            .filter(|hash| *hash != 0)
            .collect::<BTreeSet<_>>();
        let prepared =
            prepare_screen_query(&sequence, &ScreenSelectionPolicy::spatial_256(512)).unwrap();
        assert_eq!(prepared.hashes().collect::<BTreeSet<_>>(), expected);
    }
}
