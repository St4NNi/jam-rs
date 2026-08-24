//! Bounded parallel Stage-1 screening across pure `.jam` part shards.

use super::builder::{JamIndexBuildError, load_manifest};
use super::manifest::{JamIndexManifest, JamIndexPart, ScreenSelectionPolicy};
use super::part::ExternalPartReader;
use crate::jamhash_u64_v1;
use crate::reader::{JamReader, ReaderError};
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
    pub shared_hashes: Vec<SharedScreenHash>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct JamIndexScreenConfig {
    pub top_candidates: usize,
    pub accumulator_capacity: usize,
    pub min_shared_hashes: u32,
    pub min_query_windows: u32,
    pub parallel_parts: usize,
}

impl Default for JamIndexScreenConfig {
    fn default() -> Self {
        Self {
            top_candidates: 250,
            accumulator_capacity: 4_096,
            min_shared_hashes: 1,
            min_query_windows: 1,
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
    pub exact_candidates_recounted: u64,
    pub selected_candidates: u64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct JamIndexScreenResult {
    pub candidates: Vec<JamIndexCandidate>,
    pub metrics: JamIndexScreenMetrics,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct CandidateCore {
    part_id: u32,
    sample_id: u32,
    shared_hash_count: u32,
    supported_query_windows: u32,
    longest_supported_window_run: u32,
    weighted_hash_micros: i64,
}

impl CandidateCore {
    fn weighted_hash_sum(self) -> f64 {
        self.weighted_hash_micros as f64 / 1_000_000.0
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct ExactCounter {
    shared_hash_count: u32,
    weighted_hash_micros: i64,
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
    entries: HashMap<u32, ApproxEntry>,
    minimum: BinaryHeap<Reverse<(u32, u32, u64, u32)>>,
}

impl BoundedAccumulator {
    fn new(capacity: usize) -> Self {
        Self {
            capacity,
            generation: 0,
            entries: HashMap::with_capacity(capacity),
            minimum: BinaryHeap::with_capacity(capacity.saturating_mul(2)),
        }
    }

    fn add(&mut self, sample_id: u32, window_bits: u64) {
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

    fn pop_minimum(&mut self) -> Option<(u32, ApproxEntry)> {
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

    fn sample_ids(&self) -> HashSet<u32> {
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
    let mut hashes = BTreeMap::<u64, QueryHashEvidence>::new();
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
        let evidence = hashes.entry(hash).or_insert_with(|| QueryHashEvidence {
            occurrences: Vec::new(),
            windows: Vec::new(),
        });
        for window in first_window..=last_window {
            evidence.occurrences.push(QueryHashOccurrence {
                packed_kmer: kmer.0,
                query_position: position,
                query_reverse: reverse,
                query_window_id: window,
            });
            if !evidence.windows.contains(&window) {
                evidence.windows.push(window);
            }
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
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.parallel_parts)
        .build()
        .map_err(|_| JamIndexScreenError::InvalidConfig)?;
    let (document_frequencies, df_postings) = pool.install(|| {
        manifest
            .parts
            .par_iter()
            .map(|part| part_document_frequencies(root, part, query))
            .try_reduce(
                || (BTreeMap::new(), 0u64),
                |(mut left, left_postings), (right, right_postings)| {
                    for (hash, count) in right {
                        *left.entry(hash).or_insert(0u32) =
                            left.get(&hash).copied().unwrap_or(0).saturating_add(count);
                    }
                    Ok((left, left_postings.saturating_add(right_postings)))
                },
            )
    })?;
    let total_metagenomes =
        u32::try_from(manifest.total_metagenomes).map_err(|_| JamIndexScreenError::Overflow)?;
    let (selected, rank_metrics) = pool.install(|| {
        manifest
            .parts
            .par_iter()
            .map(|part| {
                rank_part(
                    root,
                    part,
                    query,
                    &document_frequencies,
                    total_metagenomes,
                    config,
                )
            })
            .try_reduce(
                || (Vec::new(), JamIndexScreenMetrics::default()),
                |(left, left_metrics), (right, right_metrics)| {
                    Ok((
                        merge_top_candidates(left, right, config.top_candidates),
                        add_metrics(left_metrics, right_metrics),
                    ))
                },
            )
    })?;
    let mut candidates = materialize_selected(
        root,
        &manifest,
        query,
        &document_frequencies,
        selected,
        config.parallel_parts,
    )?;
    candidates.sort_by(candidate_cmp);
    for (index, candidate) in candidates.iter_mut().enumerate() {
        candidate.rank = u32::try_from(index + 1).map_err(|_| JamIndexScreenError::Overflow)?;
    }
    Ok(JamIndexScreenResult {
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
    root: &Path,
    part: &JamIndexPart,
    query: &PreparedJamIndexQuery,
) -> Result<(BTreeMap<u64, u32>, u64), JamIndexScreenError> {
    let reader = open_screen(root, part)?;
    let mut frequencies = BTreeMap::new();
    let mut postings = 0u64;
    for hash in query.hashes.keys().copied() {
        let count = reader.search(hash).count();
        if count != 0 {
            frequencies.insert(
                hash,
                u32::try_from(count).map_err(|_| JamIndexScreenError::Overflow)?,
            );
            postings = postings.saturating_add(u64::try_from(count).unwrap_or(u64::MAX));
        }
    }
    Ok((frequencies, postings))
}

fn rank_part(
    root: &Path,
    part: &JamIndexPart,
    query: &PreparedJamIndexQuery,
    document_frequencies: &BTreeMap<u64, u32>,
    total_metagenomes: u32,
    config: JamIndexScreenConfig,
) -> Result<(Vec<CandidateCore>, JamIndexScreenMetrics), JamIndexScreenError> {
    let reader = open_screen(root, part)?;
    let mut accumulator = BoundedAccumulator::new(config.accumulator_capacity);
    let mut postings = 0u64;
    for (hash, evidence) in &query.hashes {
        let bits = window_bits(&evidence.windows);
        for sample_id in reader.search(*hash) {
            accumulator.add(sample_id, bits);
            postings = postings.saturating_add(1);
        }
    }
    let shortlist = accumulator.sample_ids();
    let mut exact = BTreeMap::<u32, (ExactCounter, BTreeSet<u32>)>::new();
    for (hash, evidence) in &query.hashes {
        let document_frequency = document_frequencies.get(hash).copied().unwrap_or(0);
        let weight = rarity_weight_micros(total_metagenomes, document_frequency);
        for sample_id in reader.search(*hash) {
            if !shortlist.contains(&sample_id) {
                continue;
            }
            let entry = exact.entry(sample_id).or_default();
            entry.0.shared_hash_count = entry.0.shared_hash_count.saturating_add(1);
            entry.0.weighted_hash_micros = entry.0.weighted_hash_micros.saturating_add(weight);
            entry.1.extend(evidence.windows.iter().copied());
        }
    }
    let mut candidates = exact
        .into_iter()
        .filter_map(|(sample_id, (counter, windows))| {
            let supported = u32::try_from(windows.len()).unwrap_or(u32::MAX);
            (counter.shared_hash_count >= config.min_shared_hashes
                && supported >= config.min_query_windows)
                .then_some(CandidateCore {
                    part_id: part.part_id,
                    sample_id,
                    shared_hash_count: counter.shared_hash_count,
                    supported_query_windows: supported,
                    longest_supported_window_run: longest_run(&windows),
                    weighted_hash_micros: counter.weighted_hash_micros,
                })
        })
        .collect::<Vec<_>>();
    candidates.sort_by(core_cmp);
    candidates.truncate(config.top_candidates);
    Ok((
        candidates,
        JamIndexScreenMetrics {
            hash_postings_seen: postings,
            maximum_accumulator_entries: u64::try_from(accumulator.entries.len())
                .unwrap_or(u64::MAX),
            exact_candidates_recounted: u64::try_from(shortlist.len()).unwrap_or(u64::MAX),
            ..JamIndexScreenMetrics::default()
        },
    ))
}

fn materialize_selected(
    root: &Path,
    manifest: &JamIndexManifest,
    query: &PreparedJamIndexQuery,
    document_frequencies: &BTreeMap<u64, u32>,
    selected: Vec<CandidateCore>,
    parallel_parts: usize,
) -> Result<Vec<JamIndexCandidate>, JamIndexScreenError> {
    let mut grouped = BTreeMap::<u32, Vec<CandidateCore>>::new();
    for candidate in selected {
        grouped
            .entry(candidate.part_id)
            .or_default()
            .push(candidate);
    }
    let groups = grouped.into_iter().collect::<Vec<_>>();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(parallel_parts)
        .build()
        .map_err(|_| JamIndexScreenError::InvalidConfig)?;
    let nested = pool.install(|| {
        groups
            .par_iter()
            .map(|(part_id, candidates)| {
                let part = manifest
                    .parts
                    .get(usize::try_from(*part_id).map_err(|_| JamIndexScreenError::Overflow)?)
                    .ok_or(JamIndexScreenError::UnknownPart(*part_id))?;
                materialize_part(root, part, query, document_frequencies, candidates)
            })
            .collect::<Result<Vec<_>, JamIndexScreenError>>()
    })?;
    Ok(nested.into_iter().flatten().collect())
}

fn materialize_part(
    root: &Path,
    part: &JamIndexPart,
    query: &PreparedJamIndexQuery,
    document_frequencies: &BTreeMap<u64, u32>,
    selected: &[CandidateCore],
) -> Result<Vec<JamIndexCandidate>, JamIndexScreenError> {
    let screen = open_screen(root, part)?;
    let data = ExternalPartReader::open(root.join(&part.directory).join(&part.data_file))?;
    let mut by_sample = selected
        .iter()
        .enumerate()
        .map(|(index, candidate)| (candidate.sample_id, index))
        .collect::<HashMap<_, _>>();
    let mut details = vec![Vec::<SharedScreenHash>::new(); selected.len()];
    for (hash, evidence) in &query.hashes {
        for hit in screen.search_entries(*hash) {
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
    }
    by_sample.clear();
    selected
        .iter()
        .copied()
        .zip(details)
        .map(|(candidate, mut shared_hashes)| {
            shared_hashes.sort_by_key(|shared| shared.hash);
            let metagenome = data
                .metagenomes()
                .get(
                    usize::try_from(candidate.sample_id)
                        .map_err(|_| JamIndexScreenError::Overflow)?,
                )
                .ok_or(JamIndexScreenError::UnknownSample {
                    part_id: part.part_id,
                    sample_id: candidate.sample_id,
                })?;
            if screen.sample_name(candidate.sample_id) != Some(metagenome.metagenome_id.as_str()) {
                return Err(JamIndexScreenError::PartBinding(part.part_id));
            }
            Ok(JamIndexCandidate {
                part_id: candidate.part_id,
                metagenome_local_id: candidate.sample_id,
                metagenome_id: metagenome.metagenome_id.clone(),
                rank: 0,
                shared_hash_count: candidate.shared_hash_count,
                supported_query_windows: candidate.supported_query_windows,
                longest_supported_window_run: candidate.longest_supported_window_run,
                weighted_hash_sum: candidate.weighted_hash_sum(),
                shared_hashes,
            })
        })
        .collect()
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

fn merge_top_candidates(
    mut left: Vec<CandidateCore>,
    right: Vec<CandidateCore>,
    top_candidates: usize,
) -> Vec<CandidateCore> {
    left.extend(right);
    left.sort_by(core_cmp);
    left.truncate(top_candidates);
    left
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
        exact_candidates_recounted: left
            .exact_candidates_recounted
            .saturating_add(right.exact_candidates_recounted),
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
