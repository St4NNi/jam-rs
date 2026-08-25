//! Position-free Stage-2 mapping of selected screen hashes to bounded contigs.

use super::builder::{JamIndexBuildError, load_manifest};
use super::part::JamIndexPartReader;
use super::screen::{JamIndexCandidate, PreparedJamIndexQuery};
use crate::trace::model::CandidateAdmissionSource;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeMap, BTreeSet, BinaryHeap, HashMap, HashSet};
use std::ops::Range;
use std::path::Path;
use thiserror::Error;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct JamIndexContigSearchConfig {
    pub initial_contigs_per_candidate: usize,
    pub max_contigs_per_candidate: usize,
    pub accumulator_capacity: usize,
    pub max_ranked_contig_bases: u64,
    pub strong_candidate_shared_hashes: u32,
    pub min_spatial_signatures: u32,
    pub min_query_windows: u32,
    pub rare_rescue_max_document_frequency: Option<u32>,
    pub whole_sample_min_shared_hashes: u32,
    pub parallel_parts: usize,
}

impl Default for JamIndexContigSearchConfig {
    fn default() -> Self {
        Self {
            initial_contigs_per_candidate: 8,
            max_contigs_per_candidate: 64,
            accumulator_capacity: 512,
            max_ranked_contig_bases: 64 * 1024 * 1024,
            strong_candidate_shared_hashes: 4,
            min_spatial_signatures: 2,
            min_query_windows: 2,
            rare_rescue_max_document_frequency: None,
            whole_sample_min_shared_hashes: 2,
            parallel_parts: 1,
        }
    }
}

impl JamIndexContigSearchConfig {
    fn validate(self) -> Result<Self, JamIndexContigSearchError> {
        if self.initial_contigs_per_candidate == 0
            || self.max_contigs_per_candidate < self.initial_contigs_per_candidate
            || self.accumulator_capacity < self.max_contigs_per_candidate
            || self.max_ranked_contig_bases == 0
            || self.strong_candidate_shared_hashes == 0
            || self.min_spatial_signatures == 0
            || self.min_query_windows == 0
            || self.rare_rescue_max_document_frequency == Some(0)
            || self.whole_sample_min_shared_hashes == 0
            || self.parallel_parts == 0
        {
            return Err(JamIndexContigSearchError::InvalidConfig);
        }
        Ok(self)
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RankedJamIndexContig {
    pub contig_id: u32,
    pub contig_name: String,
    pub base_count: u64,
    pub shared_hash_count: u32,
    pub supported_query_windows: u32,
    pub longest_supported_window_run: u32,
    pub weighted_hash_sum: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct JamIndexContigPlan {
    pub candidate_rank: u32,
    pub part_id: u32,
    pub metagenome_local_id: u32,
    pub metagenome_id: String,
    pub ranked_contigs: Vec<RankedJamIndexContig>,
    pub initial_contig_count: u32,
    pub sequential_fallback_range: Option<Range<u32>>,
    pub contig_truncated: bool,
    pub shared_spatial_signatures: u32,
    pub rare_shared_signatures: u32,
    pub shared_whole_sample_signatures: u32,
    pub candidate_entry_reason: CandidateAdmissionSource,
}

impl JamIndexContigPlan {
    #[must_use]
    pub fn initial_contigs(&self) -> &[RankedJamIndexContig] {
        &self.ranked_contigs[..usize::try_from(self.initial_contig_count)
            .unwrap_or(self.ranked_contigs.len())
            .min(self.ranked_contigs.len())]
    }

    #[must_use]
    pub fn weaker_contigs(&self) -> &[RankedJamIndexContig] {
        &self.ranked_contigs[usize::try_from(self.initial_contig_count)
            .unwrap_or(self.ranked_contigs.len())
            .min(self.ranked_contigs.len())..]
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub struct JamIndexContigSearchMetrics {
    pub candidates_processed: u64,
    pub posting_contigs_seen: u64,
    pub maximum_accumulator_entries: u64,
    pub accumulator_evictions: u64,
    pub ranked_contigs: u64,
    pub ranked_contig_bases: u64,
    pub contig_limit_hits: u64,
    pub base_limit_hits: u64,
    pub sequential_fallback_candidates: u64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct JamIndexContigSearchResult {
    pub plans: Vec<JamIndexContigPlan>,
    pub metrics: JamIndexContigSearchMetrics,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct JamIndexBatchContigSearchResult {
    pub queries: Vec<JamIndexContigSearchResult>,
    pub parts_opened: u64,
}

#[derive(Clone, Copy, Debug)]
struct ApproxEntry {
    count: u32,
    generation: u64,
}

#[derive(Debug)]
struct ContigAccumulator {
    capacity: usize,
    generation: u64,
    evictions: u64,
    entries: HashMap<u32, ApproxEntry>,
    minimum: BinaryHeap<Reverse<(u32, u64, u32)>>,
}

impl ContigAccumulator {
    fn new(capacity: usize) -> Self {
        Self {
            capacity,
            generation: 0,
            evictions: 0,
            entries: HashMap::with_capacity(capacity),
            minimum: BinaryHeap::with_capacity(capacity.saturating_mul(2)),
        }
    }

    fn add(&mut self, contig_id: u32) {
        self.generation = self.generation.wrapping_add(1);
        if let Some(entry) = self.entries.get_mut(&contig_id) {
            entry.count = entry.count.saturating_add(1);
            entry.generation = self.generation;
            self.minimum
                .push(Reverse((entry.count, entry.generation, contig_id)));
            return;
        }
        let base = if self.entries.len() < self.capacity {
            0
        } else {
            self.evictions = self.evictions.saturating_add(1);
            self.pop_minimum().map_or(0, |entry| entry.count)
        };
        let entry = ApproxEntry {
            count: base.saturating_add(1),
            generation: self.generation,
        };
        self.entries.insert(contig_id, entry);
        self.minimum
            .push(Reverse((entry.count, entry.generation, contig_id)));
    }

    fn pop_minimum(&mut self) -> Option<ApproxEntry> {
        while let Some(Reverse((count, generation, contig_id))) = self.minimum.pop() {
            let Some(entry) = self.entries.get(&contig_id).copied() else {
                continue;
            };
            if (entry.count, entry.generation) != (count, generation) {
                continue;
            }
            self.entries.remove(&contig_id);
            return Some(entry);
        }
        None
    }

    fn shortlist(&self) -> HashSet<u32> {
        self.entries.keys().copied().collect()
    }
}

#[derive(Clone, Debug, Default)]
struct ExactContigCounter {
    hashes: BTreeSet<u64>,
    windows: BTreeSet<u32>,
    weighted_hash_micros: i64,
}

pub fn select_candidate_contigs(
    root: impl AsRef<Path>,
    query: &PreparedJamIndexQuery,
    candidates: &[JamIndexCandidate],
    config: JamIndexContigSearchConfig,
) -> Result<JamIndexContigSearchResult, JamIndexContigSearchError> {
    let mut batch = select_candidate_contigs_batch(
        root,
        std::slice::from_ref(query),
        std::slice::from_ref(&candidates.to_vec()),
        config,
    )?;
    batch
        .queries
        .pop()
        .ok_or(JamIndexContigSearchError::EmptyBatch)
}

pub fn select_candidate_contigs_batch(
    root: impl AsRef<Path>,
    queries: &[PreparedJamIndexQuery],
    candidates: &[Vec<JamIndexCandidate>],
    config: JamIndexContigSearchConfig,
) -> Result<JamIndexBatchContigSearchResult, JamIndexContigSearchError> {
    select_candidate_contigs_batch_with_readers(root, queries, candidates, config)
        .map(|(result, _)| result)
}

pub(crate) fn select_candidate_contigs_batch_with_readers(
    root: impl AsRef<Path>,
    queries: &[PreparedJamIndexQuery],
    candidates: &[Vec<JamIndexCandidate>],
    config: JamIndexContigSearchConfig,
) -> Result<
    (
        JamIndexBatchContigSearchResult,
        BTreeMap<u32, JamIndexPartReader>,
    ),
    JamIndexContigSearchError,
> {
    let _profile = crate::profiling::process_scope("contig_posting_lookup");
    let config = config.validate()?;
    let root = root.as_ref();
    let manifest = load_manifest(root)?;
    if queries.is_empty() || queries.len() != candidates.len() {
        return Err(JamIndexContigSearchError::EmptyBatch);
    }
    let requested_parts = candidates
        .iter()
        .flatten()
        .map(|candidate| candidate.part_id)
        .collect::<BTreeSet<_>>();
    let readers = requested_parts
        .iter()
        .map(|part_id| {
            let part = manifest
                .parts
                .iter()
                .find(|part| part.part_id == *part_id)
                .ok_or(JamIndexContigSearchError::UnknownPart(*part_id))?;
            Ok((
                *part_id,
                JamIndexPartReader::open(root.join(&part.directory).join(&part.data_file))?,
            ))
        })
        .collect::<Result<BTreeMap<_, _>, JamIndexContigSearchError>>()?;
    crate::profiling::add_counter(
        "part_bin_bytes_read",
        readers
            .values()
            .map(JamIndexPartReader::object_size)
            .fold(0u64, u64::saturating_add),
    );
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.parallel_parts)
        .build()
        .map_err(|_| JamIndexContigSearchError::InvalidConfig)?;
    let results = pool.install(|| {
        queries
            .par_iter()
            .zip(candidates.par_iter())
            .map(|(query, candidates)| {
                let _worker = crate::profiling::worker("contig_posting_lookup");
                select_open_index(
                    &readers,
                    manifest.total_metagenomes,
                    query,
                    candidates,
                    config,
                )
            })
            .collect::<Result<Vec<_>, JamIndexContigSearchError>>()
    })?;
    Ok((
        JamIndexBatchContigSearchResult {
            queries: results,
            parts_opened: u64::try_from(readers.len()).unwrap_or(u64::MAX),
        },
        readers,
    ))
}

fn select_open_index(
    readers: &BTreeMap<u32, JamIndexPartReader>,
    total_metagenomes: u64,
    query: &PreparedJamIndexQuery,
    candidates: &[JamIndexCandidate],
    config: JamIndexContigSearchConfig,
) -> Result<JamIndexContigSearchResult, JamIndexContigSearchError> {
    let mut grouped = BTreeMap::<u32, Vec<&JamIndexCandidate>>::new();
    for candidate in candidates {
        grouped
            .entry(candidate.part_id)
            .or_default()
            .push(candidate);
    }
    let mut nested = Vec::new();
    for (part_id, candidates) in grouped {
        let reader = readers
            .get(&part_id)
            .ok_or(JamIndexContigSearchError::UnknownPart(part_id))?;
        nested.push(
            candidates
                .iter()
                .map(|candidate| {
                    let metagenome = reader
                        .metagenomes()
                        .get(
                            usize::try_from(candidate.metagenome_local_id)
                                .map_err(|_| JamIndexContigSearchError::Overflow)?,
                        )
                        .ok_or(JamIndexContigSearchError::CandidateBinding)?;
                    if metagenome.metagenome_id != candidate.metagenome_id {
                        return Err(JamIndexContigSearchError::CandidateBinding);
                    }
                    rank_candidate(reader, query, candidate, config, total_metagenomes)
                })
                .collect::<Result<Vec<_>, JamIndexContigSearchError>>()?,
        );
    }
    let mut plans_and_metrics = nested.into_iter().flatten().collect::<Vec<_>>();
    plans_and_metrics
        .sort_by_key(|(plan, _)| plan.as_ref().map_or(u32::MAX, |plan| plan.candidate_rank));
    let mut metrics = JamIndexContigSearchMetrics::default();
    let mut plans = Vec::with_capacity(plans_and_metrics.len());
    for (plan, part_metrics) in plans_and_metrics {
        metrics = add_metrics(metrics, part_metrics);
        if let Some(plan) = plan {
            plans.push(plan);
        }
    }
    Ok(JamIndexContigSearchResult { plans, metrics })
}

fn rank_candidate(
    reader: &JamIndexPartReader,
    _query: &PreparedJamIndexQuery,
    candidate: &JamIndexCandidate,
    config: JamIndexContigSearchConfig,
    total_metagenomes: u64,
) -> Result<(Option<JamIndexContigPlan>, JamIndexContigSearchMetrics), JamIndexContigSearchError> {
    let mut spatial_hashes = BTreeSet::new();
    let mut whole_hashes = BTreeSet::new();
    let mut spatial_windows = BTreeSet::new();
    let mut rare_spatial_hashes = BTreeSet::new();
    for shared in &candidate.shared_hashes {
        let kind = reader.posting_kind(shared.entry_ordinal)?;
        if kind.spatial {
            spatial_hashes.insert(shared.hash);
            spatial_windows.extend(
                shared
                    .occurrences
                    .iter()
                    .map(|occurrence| occurrence.query_window_id),
            );
            if config
                .rare_rescue_max_document_frequency
                .is_some_and(|limit| shared.document_frequency <= limit)
            {
                rare_spatial_hashes.insert(shared.hash);
            }
        }
        if kind.whole_sample {
            whole_hashes.insert(shared.hash);
        }
    }
    let shared_spatial_signatures = u32::try_from(spatial_hashes.len()).unwrap_or(u32::MAX);
    let shared_whole_sample_signatures = u32::try_from(whole_hashes.len()).unwrap_or(u32::MAX);
    let rare_shared_signatures = u32::try_from(rare_spatial_hashes.len()).unwrap_or(u32::MAX);
    let candidate_entry_reason = if shared_spatial_signatures >= config.min_spatial_signatures {
        Some(CandidateAdmissionSource::Standard)
    } else if u32::try_from(spatial_windows.len()).unwrap_or(u32::MAX) >= config.min_query_windows {
        Some(CandidateAdmissionSource::WindowSpread)
    } else if rare_shared_signatures != 0 {
        Some(CandidateAdmissionSource::RareRescue)
    } else if shared_whole_sample_signatures >= config.whole_sample_min_shared_hashes {
        Some(CandidateAdmissionSource::WholeSampleFallback)
    } else {
        None
    };
    let Some(candidate_entry_reason) = candidate_entry_reason else {
        return Ok((
            None,
            JamIndexContigSearchMetrics {
                candidates_processed: 1,
                ..JamIndexContigSearchMetrics::default()
            },
        ));
    };
    let contig_range = reader.metagenome_contigs(candidate.metagenome_local_id)?;
    let mut accumulator = ContigAccumulator::new(config.accumulator_capacity);
    let mut posting_contigs_seen = 0u64;
    for shared in &candidate.shared_hashes {
        for contig_id in reader.posting(shared.entry_ordinal)? {
            if in_range(contig_id, &contig_range) {
                accumulator.add(contig_id);
                posting_contigs_seen = posting_contigs_seen.saturating_add(1);
            }
        }
    }
    let shortlist = accumulator.shortlist();
    let mut exact = BTreeMap::<u32, ExactContigCounter>::new();
    for shared in &candidate.shared_hashes {
        let weight = rarity_weight_micros(total_metagenomes, shared.document_frequency);
        let windows = shared
            .occurrences
            .iter()
            .map(|occurrence| occurrence.query_window_id)
            .collect::<BTreeSet<_>>();
        for contig_id in reader.posting(shared.entry_ordinal)? {
            if !shortlist.contains(&contig_id) || !in_range(contig_id, &contig_range) {
                continue;
            }
            let counter = exact.entry(contig_id).or_default();
            if counter.hashes.insert(shared.hash) {
                counter.weighted_hash_micros = counter.weighted_hash_micros.saturating_add(weight);
            }
            counter.windows.extend(windows.iter().copied());
        }
    }
    let mut ranked = exact
        .into_iter()
        .map(|(contig_id, counter)| {
            let contig = reader
                .contigs()
                .get(usize::try_from(contig_id).map_err(|_| JamIndexContigSearchError::Overflow)?)
                .ok_or(JamIndexContigSearchError::UnknownContig(contig_id))?;
            Ok(RankedJamIndexContig {
                contig_id,
                contig_name: contig.name.clone(),
                base_count: contig.base_count,
                shared_hash_count: u32::try_from(counter.hashes.len()).unwrap_or(u32::MAX),
                supported_query_windows: u32::try_from(counter.windows.len()).unwrap_or(u32::MAX),
                longest_supported_window_run: longest_run(&counter.windows),
                weighted_hash_sum: counter.weighted_hash_micros as f64 / 1_000_000.0,
            })
        })
        .collect::<Result<Vec<_>, JamIndexContigSearchError>>()?;
    ranked.sort_by(contig_cmp);
    let mut retained = Vec::new();
    let mut retained_bases = 0u64;
    let mut contig_limit_hit = false;
    let mut base_limit_hit = false;
    for contig in ranked {
        if retained.len() == config.max_contigs_per_candidate {
            contig_limit_hit = true;
            break;
        }
        if !retained.is_empty()
            && retained_bases.saturating_add(contig.base_count) > config.max_ranked_contig_bases
        {
            base_limit_hit = true;
            break;
        }
        retained_bases = retained_bases.saturating_add(contig.base_count);
        retained.push(contig);
    }
    let sequential_fallback_range = (retained.is_empty()
        && candidate.shared_hash_count >= config.strong_candidate_shared_hashes)
        .then_some(contig_range);
    let initial_contig_count =
        u32::try_from(retained.len().min(config.initial_contigs_per_candidate))
            .map_err(|_| JamIndexContigSearchError::Overflow)?;
    let metrics = JamIndexContigSearchMetrics {
        candidates_processed: 1,
        posting_contigs_seen,
        maximum_accumulator_entries: u64::try_from(accumulator.entries.len()).unwrap_or(u64::MAX),
        accumulator_evictions: accumulator.evictions,
        ranked_contigs: u64::try_from(retained.len()).unwrap_or(u64::MAX),
        ranked_contig_bases: retained_bases,
        contig_limit_hits: u64::from(contig_limit_hit),
        base_limit_hits: u64::from(base_limit_hit),
        sequential_fallback_candidates: u64::from(sequential_fallback_range.is_some()),
    };
    Ok((
        Some(JamIndexContigPlan {
            candidate_rank: candidate.rank,
            part_id: candidate.part_id,
            metagenome_local_id: candidate.metagenome_local_id,
            metagenome_id: candidate.metagenome_id.clone(),
            ranked_contigs: retained,
            initial_contig_count,
            sequential_fallback_range,
            contig_truncated: accumulator.evictions != 0 || contig_limit_hit || base_limit_hit,
            shared_spatial_signatures,
            rare_shared_signatures,
            shared_whole_sample_signatures,
            candidate_entry_reason,
        }),
        metrics,
    ))
}

fn in_range(contig_id: u32, range: &Range<u32>) -> bool {
    contig_id >= range.start && contig_id < range.end
}

fn contig_cmp(left: &RankedJamIndexContig, right: &RankedJamIndexContig) -> Ordering {
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
        .then_with(|| left.base_count.cmp(&right.base_count))
        .then_with(|| left.contig_id.cmp(&right.contig_id))
}

fn rarity_weight_micros(total_metagenomes: u64, document_frequency: u32) -> i64 {
    ((((total_metagenomes as f64) + 1.0) / (f64::from(document_frequency) + 1.0)).ln()
        * 1_000_000.0)
        .round() as i64
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

fn add_metrics(
    left: JamIndexContigSearchMetrics,
    right: JamIndexContigSearchMetrics,
) -> JamIndexContigSearchMetrics {
    JamIndexContigSearchMetrics {
        candidates_processed: left
            .candidates_processed
            .saturating_add(right.candidates_processed),
        posting_contigs_seen: left
            .posting_contigs_seen
            .saturating_add(right.posting_contigs_seen),
        maximum_accumulator_entries: left
            .maximum_accumulator_entries
            .max(right.maximum_accumulator_entries),
        accumulator_evictions: left
            .accumulator_evictions
            .saturating_add(right.accumulator_evictions),
        ranked_contigs: left.ranked_contigs.saturating_add(right.ranked_contigs),
        ranked_contig_bases: left
            .ranked_contig_bases
            .saturating_add(right.ranked_contig_bases),
        contig_limit_hits: left
            .contig_limit_hits
            .saturating_add(right.contig_limit_hits),
        base_limit_hits: left.base_limit_hits.saturating_add(right.base_limit_hits),
        sequential_fallback_candidates: left
            .sequential_fallback_candidates
            .saturating_add(right.sequential_fallback_candidates),
    }
}

#[derive(Debug, Error)]
pub enum JamIndexContigSearchError {
    #[error("invalid Jam Index contig search config")]
    InvalidConfig,
    #[error("Jam Index contig search batch is empty or mismatched")]
    EmptyBatch,
    #[error("Jam Index screen candidate does not match its data part")]
    CandidateBinding,
    #[error("unknown Jam Index part {0}")]
    UnknownPart(u32),
    #[error("unknown Jam Index contig {0}")]
    UnknownContig(u32),
    #[error("Jam Index contig search coordinate overflow")]
    Overflow,
    #[error(transparent)]
    Build(#[from] JamIndexBuildError),
    #[error(transparent)]
    Part(#[from] super::part::JamIndexPartError),
}
