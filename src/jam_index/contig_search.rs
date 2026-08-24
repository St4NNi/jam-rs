//! Position-free Stage-2 mapping of selected screen hashes to bounded contigs.

use super::builder::{JamIndexBuildError, load_manifest};
use super::part::{JamIndexPartReader, SignatureHit};
use super::screen::{JamIndexCandidate, PreparedJamIndexQuery};
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
    pub signature_hits_seen: u64,
    pub maximum_accumulator_entries: u64,
    pub ranked_contigs: u64,
    pub ranked_contig_bases: u64,
    pub sequential_fallback_candidates: u64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct JamIndexContigSearchResult {
    pub plans: Vec<JamIndexContigPlan>,
    pub metrics: JamIndexContigSearchMetrics,
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
    entries: HashMap<u32, ApproxEntry>,
    minimum: BinaryHeap<Reverse<(u32, u64, u32)>>,
}

impl ContigAccumulator {
    fn new(capacity: usize) -> Self {
        Self {
            capacity,
            generation: 0,
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
    let config = config.validate()?;
    let root = root.as_ref();
    let manifest = load_manifest(root)?;
    let mut grouped = BTreeMap::<u32, Vec<&JamIndexCandidate>>::new();
    for candidate in candidates {
        grouped
            .entry(candidate.part_id)
            .or_default()
            .push(candidate);
    }
    let groups = grouped.into_iter().collect::<Vec<_>>();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.parallel_parts)
        .build()
        .map_err(|_| JamIndexContigSearchError::InvalidConfig)?;
    let nested = pool.install(|| {
        groups
            .par_iter()
            .map(|(part_id, candidates)| {
                let part = manifest
                    .parts
                    .get(
                        usize::try_from(*part_id)
                            .map_err(|_| JamIndexContigSearchError::Overflow)?,
                    )
                    .ok_or(JamIndexContigSearchError::UnknownPart(*part_id))?;
                let reader =
                    JamIndexPartReader::open(root.join(&part.directory).join(&part.data_file))?;
                candidates
                    .iter()
                    .map(|candidate| {
                        rank_candidate(
                            &reader,
                            query,
                            candidate,
                            config,
                            manifest.total_metagenomes,
                        )
                    })
                    .collect::<Result<Vec<_>, JamIndexContigSearchError>>()
            })
            .collect::<Result<Vec<_>, JamIndexContigSearchError>>()
    })?;
    let mut plans_and_metrics = nested.into_iter().flatten().collect::<Vec<_>>();
    plans_and_metrics.sort_by_key(|(plan, _)| plan.candidate_rank);
    let mut metrics = JamIndexContigSearchMetrics::default();
    let mut plans = Vec::with_capacity(plans_and_metrics.len());
    for (plan, part_metrics) in plans_and_metrics {
        metrics = add_metrics(metrics, part_metrics);
        plans.push(plan);
    }
    Ok(JamIndexContigSearchResult { plans, metrics })
}

fn rank_candidate(
    reader: &JamIndexPartReader,
    _query: &PreparedJamIndexQuery,
    candidate: &JamIndexCandidate,
    config: JamIndexContigSearchConfig,
    total_metagenomes: u64,
) -> Result<(JamIndexContigPlan, JamIndexContigSearchMetrics), JamIndexContigSearchError> {
    let contig_range = reader.contig_ids_for_metagenome(candidate.metagenome_local_id)?;
    let mut accumulator = ContigAccumulator::new(config.accumulator_capacity);
    let mut signature_hits_seen = 0u64;
    for shared in &candidate.shared_hashes {
        reader.visit_signature_hits(shared.hash, &mut |hit| {
            if in_range(hit, &contig_range) {
                accumulator.add(hit.contig_id);
                signature_hits_seen = signature_hits_seen.saturating_add(1);
            }
        })?;
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
        reader.visit_signature_hits(shared.hash, &mut |hit| {
            if !shortlist.contains(&hit.contig_id) || !in_range(hit, &contig_range) {
                return;
            }
            let counter = exact.entry(hit.contig_id).or_default();
            if counter.hashes.insert(shared.hash) {
                counter.weighted_hash_micros = counter.weighted_hash_micros.saturating_add(weight);
            }
            counter.windows.extend(windows.iter().copied());
        })?;
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
    for contig in ranked {
        if retained.len() == config.max_contigs_per_candidate {
            break;
        }
        if !retained.is_empty()
            && retained_bases.saturating_add(contig.base_count) > config.max_ranked_contig_bases
        {
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
        signature_hits_seen,
        maximum_accumulator_entries: u64::try_from(accumulator.entries.len()).unwrap_or(u64::MAX),
        ranked_contigs: u64::try_from(retained.len()).unwrap_or(u64::MAX),
        ranked_contig_bases: retained_bases,
        sequential_fallback_candidates: u64::from(sequential_fallback_range.is_some()),
    };
    Ok((
        JamIndexContigPlan {
            candidate_rank: candidate.rank,
            part_id: candidate.part_id,
            metagenome_local_id: candidate.metagenome_local_id,
            metagenome_id: candidate.metagenome_id.clone(),
            ranked_contigs: retained,
            initial_contig_count,
            sequential_fallback_range,
        },
        metrics,
    ))
}

fn in_range(hit: SignatureHit, range: &Range<u32>) -> bool {
    hit.contig_id >= range.start && hit.contig_id < range.end
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
        signature_hits_seen: left
            .signature_hits_seen
            .saturating_add(right.signature_hits_seen),
        maximum_accumulator_entries: left
            .maximum_accumulator_entries
            .max(right.maximum_accumulator_entries),
        ranked_contigs: left.ranked_contigs.saturating_add(right.ranked_contigs),
        ranked_contig_bases: left
            .ranked_contig_bases
            .saturating_add(right.ranked_contig_bases),
        sequential_fallback_candidates: left
            .sequential_fallback_candidates
            .saturating_add(right.sequential_fallback_candidates),
    }
}

#[derive(Debug, Error)]
pub enum JamIndexContigSearchError {
    #[error("invalid Jam Index contig search config")]
    InvalidConfig,
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
