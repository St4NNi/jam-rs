//! Native multi-query Jam Index planning and grouped execution.

use super::builder::{JamIndexBuildError, load_manifest};
use super::contig_search::{
    JamIndexBatchContigSearchResult, JamIndexContigPlan, JamIndexContigSearchError,
    JamIndexContigSearchMetrics, select_candidate_contigs_batch_with_readers,
};
use super::part::{JamIndexPartError, JamIndexPartReader};
use super::screen::{
    JamIndexCandidate, JamIndexScreenError, JamIndexScreenMetrics, PreparedJamIndexQuery,
    prepare_screen_query, search_jam_index_batch,
};
use super::trace::{
    JamIndexTraceConfig, JamIndexTraceError, JamIndexTraceMetrics, candidate_result, complete,
    run_archive, run_plan_loaded,
};
use crate::archive::ArchiveError;
use crate::trace::memory::TraceMemorySemaphore;
use crate::trace::model::{QueryKind, TopologyRequested, TraceMetagenomeResult};
use crate::trace::runner::{RunnerError, candidate_memory_estimate, merge_index_results};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::mpsc::sync_channel;
use thiserror::Error;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct JamIndexBatchQuery {
    pub query_id: String,
    pub original_header: String,
    pub sequence: Vec<u8>,
    pub sequence_sha256: String,
    pub query_kind: QueryKind,
    pub topology_requested: TopologyRequested,
}

#[derive(Clone, Debug)]
pub struct JamIndexBatchConfig {
    pub trace: JamIndexTraceConfig,
    pub max_group_contig_bases: u64,
    pub fallback_contigs_per_chunk: usize,
    pub group_memory_budget_bytes: u64,
    pub max_total_query_bases: u64,
    pub max_prepared_query_hashes: u64,
}

impl Default for JamIndexBatchConfig {
    fn default() -> Self {
        Self {
            trace: JamIndexTraceConfig::default(),
            max_group_contig_bases: 4 * 1024 * 1024 * 1024,
            fallback_contigs_per_chunk: 8,
            group_memory_budget_bytes: 8 * 1024 * 1024 * 1024,
            max_total_query_bases: 100_000_000,
            max_prepared_query_hashes: 100_000_000,
        }
    }
}

impl JamIndexBatchConfig {
    fn validate(&self) -> Result<(), JamIndexBatchError> {
        self.trace.validate()?;
        if self.max_group_contig_bases == 0
            || self.fallback_contigs_per_chunk == 0
            || self.group_memory_budget_bytes == 0
            || self.max_total_query_bases == 0
            || self.max_prepared_query_hashes == 0
        {
            return Err(JamIndexBatchError::InvalidInput(
                "batch group limits must be positive".to_string(),
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Debug)]
pub struct JamIndexBatchWork {
    pub query_index: usize,
    pub candidate: JamIndexCandidate,
    pub contig_plan: JamIndexContigPlan,
}

#[derive(Clone, Debug)]
pub struct JamIndexBatchWorkGroup {
    pub part_id: u32,
    pub metagenome_local_id: u32,
    pub metagenome_id: String,
    pub source_path: PathBuf,
    pub work: Vec<JamIndexBatchWork>,
}

#[derive(Clone, Debug)]
pub struct JamIndexPlannedQuery {
    pub input: JamIndexBatchQuery,
    pub prepared: PreparedJamIndexQuery,
    pub candidates: Vec<JamIndexCandidate>,
    pub candidate_limit_reached: bool,
    pub screen_metrics: JamIndexScreenMetrics,
    pub contig_plans: Vec<JamIndexContigPlan>,
    pub contig_metrics: JamIndexContigSearchMetrics,
}

pub struct JamIndexBatchPlan {
    pub queries: Vec<JamIndexPlannedQuery>,
    pub groups: Vec<JamIndexBatchWorkGroup>,
    pub screen_parts_opened: u64,
    pub contig_parts_opened: u64,
    readers: BTreeMap<u32, JamIndexPartReader>,
}

impl JamIndexBatchPlan {
    pub fn remap_sources(
        &mut self,
        sources: &BTreeMap<String, (PathBuf, [u8; 32])>,
    ) -> Result<(), JamIndexBatchError> {
        for group in &mut self.groups {
            let (source_path, source_sha256) =
                sources.get(&group.metagenome_id).ok_or_else(|| {
                    JamIndexBatchError::InvalidInput(format!(
                        "source override is missing candidate metagenome {}",
                        group.metagenome_id
                    ))
                })?;
            let reader = self
                .readers
                .get_mut(&group.part_id)
                .ok_or(JamIndexBatchError::UnknownPart(group.part_id))?;
            reader.remap_source(&group.metagenome_id, source_path, *source_sha256)?;
            group.source_path = source_path.clone();
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum JamIndexBatchStatusKind {
    NoCandidate,
    CandidateOnly,
    AlignedTrace,
    PartialFailure,
    Failed,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct JamIndexBatchQueryStatus {
    pub query_id: String,
    pub original_header: String,
    pub query_length: u64,
    pub query_sha256: String,
    pub query_kind: QueryKind,
    pub topology_requested: TopologyRequested,
    pub status: JamIndexBatchStatusKind,
    pub candidate_metagenomes: u64,
    pub completed_metagenomes: u64,
    pub aligned_metagenomes: u64,
    pub candidate_truncated: bool,
    pub failed: bool,
    pub error: Option<String>,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub struct JamIndexBatchMetrics {
    pub queries: u64,
    pub query_bases: u64,
    pub screen_parts_opened: u64,
    pub contig_parts_opened: u64,
    pub candidate_pairs: u64,
    pub work_groups: u64,
    pub source_open_count: u64,
    pub sequential_source_open_count: u64,
    pub source_bytes_read: u64,
    pub maximum_loaded_contig_bases: u64,
    pub results_emitted: u64,
    pub alignments_emitted: u64,
    pub failed_groups: u64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct JamIndexBatchExecution {
    pub statuses: Vec<JamIndexBatchQueryStatus>,
    pub metrics: JamIndexBatchMetrics,
}

struct JamIndexBatchGroupExecution {
    outcomes: Vec<(usize, Result<TraceMetagenomeResult, String>)>,
    metrics: JamIndexBatchMetrics,
}

pub fn plan_batch(
    root: impl AsRef<std::path::Path>,
    queries: &[JamIndexBatchQuery],
    config: &JamIndexBatchConfig,
) -> Result<JamIndexBatchPlan, JamIndexBatchError> {
    config.validate()?;
    if queries.is_empty() {
        return Err(JamIndexBatchError::InvalidInput(
            "batch requires at least one query".to_string(),
        ));
    }
    let mut query_ids = BTreeSet::new();
    if queries.iter().any(|query| {
        query.query_id.trim().is_empty()
            || query.sequence.is_empty()
            || query.sequence_sha256 != format!("{:x}", Sha256::digest(&query.sequence))
            || !query_ids.insert(query.query_id.clone())
    }) {
        return Err(JamIndexBatchError::InvalidInput(
            "query IDs, sequences, and checksums must be valid and unique".to_string(),
        ));
    }
    let root = root.as_ref();
    let manifest = load_manifest(root)?;
    let total_query_bases = queries
        .iter()
        .map(|query| u64::try_from(query.sequence.len()).unwrap_or(u64::MAX))
        .fold(0u64, u64::saturating_add);
    if total_query_bases > config.max_total_query_bases {
        return Err(JamIndexBatchError::QueryBatchLimit {
            observed: total_query_bases,
            limit: config.max_total_query_bases,
            unit: "bases",
        });
    }
    let prepared = {
        let _profile = crate::profiling::scope("query_k21_preparation");
        queries
            .iter()
            .map(|query| prepare_screen_query(&query.sequence, &manifest.selection_policy))
            .collect::<Result<Vec<_>, JamIndexScreenError>>()?
    };
    let prepared_hashes = prepared
        .iter()
        .map(|query| u64::try_from(query.unique_hash_count()).unwrap_or(u64::MAX))
        .fold(0u64, u64::saturating_add);
    if prepared_hashes > config.max_prepared_query_hashes {
        return Err(JamIndexBatchError::QueryBatchLimit {
            observed: prepared_hashes,
            limit: config.max_prepared_query_hashes,
            unit: "hashes",
        });
    }
    let screen = search_jam_index_batch(root, &prepared, config.trace.screen)?;
    let candidate_vectors = screen
        .queries
        .iter()
        .map(|result| result.candidates.clone())
        .collect::<Vec<_>>();
    let (contigs, readers): (
        JamIndexBatchContigSearchResult,
        BTreeMap<u32, JamIndexPartReader>,
    ) = select_candidate_contigs_batch_with_readers(
        root,
        &prepared,
        &candidate_vectors,
        config.trace.contigs,
    )?;
    let mut planned = Vec::with_capacity(queries.len());
    for (query_index, input) in queries.iter().cloned().enumerate() {
        let screen_result = screen
            .queries
            .get(query_index)
            .ok_or(JamIndexBatchError::Overflow)?;
        let contig_result = contigs
            .queries
            .get(query_index)
            .ok_or(JamIndexBatchError::Overflow)?;
        let mut candidates = screen_result.candidates.clone();
        candidates.retain_mut(|candidate| {
            let Some(plan) = contig_result
                .plans
                .iter()
                .find(|plan| plan.candidate_rank == candidate.rank)
            else {
                return false;
            };
            candidate.shared_spatial_signatures = plan.shared_spatial_signatures;
            candidate.rare_shared_signatures = plan.rare_shared_signatures;
            candidate.shared_whole_sample_signatures = plan.shared_whole_sample_signatures;
            candidate.admission_source = plan.candidate_entry_reason;
            true
        });
        planned.push(JamIndexPlannedQuery {
            input,
            prepared: prepared
                .get(query_index)
                .cloned()
                .ok_or(JamIndexBatchError::Overflow)?,
            candidates,
            candidate_limit_reached: screen_result.candidate_limit_reached,
            screen_metrics: screen_result.metrics,
            contig_plans: contig_result.plans.clone(),
            contig_metrics: contig_result.metrics,
        });
    }
    let groups = {
        let _profile = crate::profiling::scope("work_group_construction");
        let mut grouped = BTreeMap::<(u32, u32), Vec<JamIndexBatchWork>>::new();
        for (query_index, query) in planned.iter().enumerate() {
            for contig_plan in &query.contig_plans {
                let candidate = query
                    .candidates
                    .iter()
                    .find(|candidate| candidate.rank == contig_plan.candidate_rank)
                    .cloned()
                    .ok_or(JamIndexBatchError::CandidateBinding)?;
                grouped
                    .entry((candidate.part_id, candidate.metagenome_local_id))
                    .or_default()
                    .push(JamIndexBatchWork {
                        query_index,
                        candidate,
                        contig_plan: contig_plan.clone(),
                    });
            }
        }
        let mut groups = Vec::with_capacity(grouped.len());
        for ((part_id, metagenome_local_id), mut work) in grouped {
            work.sort_by_key(|item| item.query_index);
            let reader = readers
                .get(&part_id)
                .ok_or(JamIndexBatchError::UnknownPart(part_id))?;
            let metagenome = reader
                .metagenomes()
                .get(
                    usize::try_from(metagenome_local_id)
                        .map_err(|_| JamIndexBatchError::Overflow)?,
                )
                .ok_or(JamIndexBatchError::CandidateBinding)?;
            if work
                .iter()
                .any(|item| item.candidate.metagenome_id != metagenome.metagenome_id)
            {
                return Err(JamIndexBatchError::CandidateBinding);
            }
            groups.push(JamIndexBatchWorkGroup {
                part_id,
                metagenome_local_id,
                metagenome_id: metagenome.metagenome_id.clone(),
                source_path: metagenome.source_path.clone(),
                work,
            });
        }
        groups
    };
    Ok(JamIndexBatchPlan {
        queries: planned,
        groups,
        screen_parts_opened: screen.parts_opened,
        contig_parts_opened: contigs.parts_opened,
        readers,
    })
}

pub fn execute_batch<F>(
    plan: &JamIndexBatchPlan,
    config: &JamIndexBatchConfig,
    run_id: &str,
    mut emit: F,
) -> Result<JamIndexBatchExecution, JamIndexBatchError>
where
    F: FnMut(TraceMetagenomeResult) -> Result<(), JamIndexBatchError>,
{
    config.validate()?;
    let mut statuses = plan
        .queries
        .iter()
        .map(|query| JamIndexBatchQueryStatus {
            query_id: query.input.query_id.clone(),
            original_header: query.input.original_header.clone(),
            query_length: u64::try_from(query.input.sequence.len()).unwrap_or(u64::MAX),
            query_sha256: query.input.sequence_sha256.clone(),
            query_kind: query.input.query_kind,
            topology_requested: query.input.topology_requested,
            status: if query.candidates.is_empty() {
                JamIndexBatchStatusKind::NoCandidate
            } else {
                JamIndexBatchStatusKind::CandidateOnly
            },
            candidate_metagenomes: u64::try_from(query.candidates.len()).unwrap_or(u64::MAX),
            completed_metagenomes: 0,
            aligned_metagenomes: 0,
            candidate_truncated: query.candidate_limit_reached,
            failed: false,
            error: None,
        })
        .collect::<Vec<_>>();
    let mut metrics = JamIndexBatchMetrics {
        queries: u64::try_from(plan.queries.len()).unwrap_or(u64::MAX),
        query_bases: plan
            .queries
            .iter()
            .map(|query| u64::try_from(query.input.sequence.len()).unwrap_or(u64::MAX))
            .fold(0u64, u64::saturating_add),
        screen_parts_opened: plan.screen_parts_opened,
        contig_parts_opened: plan.contig_parts_opened,
        candidate_pairs: plan
            .queries
            .iter()
            .map(|query| u64::try_from(query.candidates.len()).unwrap_or(u64::MAX))
            .fold(0u64, u64::saturating_add),
        work_groups: u64::try_from(plan.groups.len()).unwrap_or(u64::MAX),
        ..JamIndexBatchMetrics::default()
    };
    let worker_count = config.trace.runner.threads.min(plan.groups.len()).max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(worker_count)
        .build()
        .map_err(|error| JamIndexBatchError::InvalidInput(error.to_string()))?;
    let memory = TraceMemorySemaphore::new(config.group_memory_budget_bytes);
    let (sender, receiver) = sync_channel(worker_count.saturating_mul(2).max(1));
    let mut output_error = None;
    let producer_panicked = std::thread::scope(|scope| {
        let producer = scope.spawn(|| {
            pool.install(|| {
                plan.groups.par_iter().enumerate().for_each_with(
                    sender,
                    |sender, (group_index, group)| {
                        let execution = (|| {
                            let required = group_memory_bytes(plan, group, config)?;
                            let _reservation = memory.acquire_bytes(required).map_err(|error| {
                                JamIndexBatchError::InvalidInput(format!(
                                    "group memory admission failed for {}: {error}",
                                    group.metagenome_id
                                ))
                            })?;
                            let _profile = crate::profiling::scope("work_group_execution");
                            execute_group(plan, group, config, run_id)
                        })();
                        let _ = sender.send((group_index, execution));
                    },
                );
            });
        });
        for (group_index, execution) in receiver {
            let group = plan
                .groups
                .get(group_index)
                .ok_or(JamIndexBatchError::Overflow)?;
            match execution {
                Ok(execution) => {
                    add_group_metrics(&mut metrics, execution.metrics);
                    for (query_index, outcome) in execution.outcomes {
                        let status = statuses
                            .get_mut(query_index)
                            .ok_or(JamIndexBatchError::Overflow)?;
                        match outcome {
                            Ok(result) => {
                                let aligned = !result.alignments.is_empty();
                                if output_error.is_none()
                                    && let Err(error) = emit(result)
                                {
                                    output_error = Some(error);
                                }
                                status.completed_metagenomes =
                                    status.completed_metagenomes.saturating_add(1);
                                if aligned {
                                    status.aligned_metagenomes =
                                        status.aligned_metagenomes.saturating_add(1);
                                }
                            }
                            Err(error) => {
                                status.failed = true;
                                status.error = Some(error);
                            }
                        }
                    }
                }
                Err(error) => {
                    metrics.failed_groups = metrics.failed_groups.saturating_add(1);
                    for work in &group.work {
                        let status = statuses
                            .get_mut(work.query_index)
                            .ok_or(JamIndexBatchError::Overflow)?;
                        status.failed = true;
                        status.error = Some(error.to_string());
                    }
                }
            }
        }
        Ok::<_, JamIndexBatchError>(producer.join().is_err())
    })?;
    if producer_panicked {
        return Err(JamIndexBatchError::InvalidInput(
            "work-group worker pool panicked".to_string(),
        ));
    }
    if let Some(error) = output_error {
        return Err(error);
    }
    for status in &mut statuses {
        status.status = if status.failed {
            if status.completed_metagenomes == 0 {
                JamIndexBatchStatusKind::Failed
            } else {
                JamIndexBatchStatusKind::PartialFailure
            }
        } else if status.aligned_metagenomes != 0 {
            JamIndexBatchStatusKind::AlignedTrace
        } else if status.candidate_metagenomes == 0 {
            JamIndexBatchStatusKind::NoCandidate
        } else {
            JamIndexBatchStatusKind::CandidateOnly
        };
    }
    Ok(JamIndexBatchExecution { statuses, metrics })
}

fn group_memory_bytes(
    plan: &JamIndexBatchPlan,
    group: &JamIndexBatchWorkGroup,
    config: &JamIndexBatchConfig,
) -> Result<u64, JamIndexBatchError> {
    let reader = plan
        .readers
        .get(&group.part_id)
        .ok_or(JamIndexBatchError::UnknownPart(group.part_id))?;
    let metagenome = reader
        .metagenomes()
        .get(usize::try_from(group.metagenome_local_id).map_err(|_| JamIndexBatchError::Overflow)?)
        .ok_or(JamIndexBatchError::CandidateBinding)?;
    let fallback = group
        .work
        .iter()
        .any(|work| work.contig_plan.sequential_fallback_range.is_some());
    let selected_bases = if fallback {
        metagenome.total_bases
    } else {
        let selected = group
            .work
            .iter()
            .flat_map(|work| {
                work.contig_plan
                    .ranked_contigs
                    .iter()
                    .map(|contig| contig.contig_id)
            })
            .collect::<BTreeSet<_>>();
        selected.iter().try_fold(0u64, |total, contig_id| {
            let contig = reader
                .contigs()
                .get(usize::try_from(*contig_id).map_err(|_| JamIndexBatchError::Overflow)?)
                .ok_or(JamIndexBatchError::CandidateBinding)?;
            Ok::<_, JamIndexBatchError>(total.saturating_add(contig.base_count))
        })?
    };
    let query_work = group
        .work
        .iter()
        .map(|work| {
            plan.queries
                .get(work.query_index)
                .map(|query| {
                    candidate_memory_estimate(&config.trace.runner, query.input.sequence.len())
                        .total_bytes()
                })
                .ok_or(JamIndexBatchError::Overflow)
        })
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .max()
        .unwrap_or(0);
    Ok(metagenome
        .source_size
        .saturating_add(selected_bases)
        .saturating_add(query_work)
        .saturating_add(
            u64::try_from(group.work.len())
                .unwrap_or(u64::MAX)
                .saturating_mul(64 * 1024),
        ))
}

fn add_group_metrics(total: &mut JamIndexBatchMetrics, group: JamIndexBatchMetrics) {
    total.source_open_count = total
        .source_open_count
        .saturating_add(group.source_open_count);
    total.sequential_source_open_count = total
        .sequential_source_open_count
        .saturating_add(group.sequential_source_open_count);
    total.source_bytes_read = total
        .source_bytes_read
        .saturating_add(group.source_bytes_read);
    total.maximum_loaded_contig_bases = total
        .maximum_loaded_contig_bases
        .max(group.maximum_loaded_contig_bases);
    total.results_emitted = total.results_emitted.saturating_add(group.results_emitted);
    total.alignments_emitted = total
        .alignments_emitted
        .saturating_add(group.alignments_emitted);
    total.failed_groups = total.failed_groups.saturating_add(group.failed_groups);
}

fn execute_group(
    plan: &JamIndexBatchPlan,
    group: &JamIndexBatchWorkGroup,
    config: &JamIndexBatchConfig,
    run_id: &str,
) -> Result<JamIndexBatchGroupExecution, JamIndexBatchError> {
    let mut metrics = JamIndexBatchMetrics::default();
    let mut results = BTreeMap::<usize, TraceMetagenomeResult>::new();
    let reader = plan
        .readers
        .get(&group.part_id)
        .ok_or(JamIndexBatchError::UnknownPart(group.part_id))?;
    let fallback = group
        .work
        .iter()
        .filter(|work| work.contig_plan.sequential_fallback_range.is_some())
        .map(|work| work.query_index)
        .collect::<BTreeSet<_>>();
    let selected_ids = group
        .work
        .iter()
        .filter(|work| work.contig_plan.sequential_fallback_range.is_none())
        .flat_map(|work| {
            work.contig_plan
                .ranked_contigs
                .iter()
                .map(|contig| contig.contig_id)
        })
        .collect::<BTreeSet<_>>();
    let selected_bases = selected_ids.iter().try_fold(0u64, |total, contig_id| {
        let contig = reader
            .contigs()
            .get(usize::try_from(*contig_id).map_err(|_| JamIndexBatchError::Overflow)?)
            .ok_or(JamIndexBatchError::CandidateBinding)?;
        Ok::<_, JamIndexBatchError>(total.saturating_add(contig.base_count))
    })?;
    if selected_bases > config.max_group_contig_bases {
        return Err(JamIndexBatchError::GroupTooLarge {
            metagenome_id: group.metagenome_id.clone(),
            selected_bases,
            limit: config.max_group_contig_bases,
        });
    }
    crate::profiling::add_counter("selected_contig_bases", selected_bases);
    let mut loaded = BTreeMap::<u32, Arc<[u8]>>::new();
    let mut fallback_results = BTreeMap::<usize, TraceMetagenomeResult>::new();
    let mut query_failures = BTreeMap::<usize, String>::new();
    let _source_profile = crate::profiling::scope("source_staging");
    let source_bytes = if fallback.is_empty() {
        let read = reader.read_contigs(
            group.metagenome_local_id,
            &selected_ids.iter().copied().collect::<Vec<_>>(),
        )?;
        let loaded_bases = read
            .contigs
            .values()
            .map(|sequence| u64::try_from(sequence.len()).unwrap_or(u64::MAX))
            .fold(0u64, u64::saturating_add);
        metrics.maximum_loaded_contig_bases = metrics.maximum_loaded_contig_bases.max(loaded_bases);
        loaded.extend(
            read.contigs
                .into_iter()
                .map(|(contig_id, sequence)| (contig_id, Arc::<[u8]>::from(sequence))),
        );
        read.source_bytes
    } else {
        metrics.sequential_source_open_count =
            metrics.sequential_source_open_count.saturating_add(1);
        reader.stream_contigs::<JamIndexBatchError, _>(
            group.metagenome_local_id,
            config.fallback_contigs_per_chunk,
            |chunk| {
                let chunk = chunk
                    .into_iter()
                    .map(|(contig_id, sequence)| (contig_id, Arc::<[u8]>::from(sequence)))
                    .collect::<BTreeMap<_, _>>();
                let active_ids = loaded
                    .keys()
                    .chain(chunk.keys())
                    .copied()
                    .collect::<BTreeSet<_>>();
                let active_bases = active_ids.iter().try_fold(0u64, |total, contig_id| {
                    let contig = reader
                        .contigs()
                        .get(
                            usize::try_from(*contig_id)
                                .map_err(|_| JamIndexBatchError::Overflow)?,
                        )
                        .ok_or(JamIndexBatchError::CandidateBinding)?;
                    Ok::<_, JamIndexBatchError>(total.saturating_add(contig.base_count))
                })?;
                if active_bases > config.max_group_contig_bases {
                    return Err(JamIndexBatchError::GroupTooLarge {
                        metagenome_id: group.metagenome_id.clone(),
                        selected_bases: active_bases,
                        limit: config.max_group_contig_bases,
                    });
                }
                metrics.maximum_loaded_contig_bases =
                    metrics.maximum_loaded_contig_bases.max(active_bases);
                for (contig_id, sequence) in &chunk {
                    if selected_ids.contains(contig_id) {
                        loaded.insert(*contig_id, Arc::clone(sequence));
                    }
                }
                for work in group
                    .work
                    .iter()
                    .filter(|work| fallback.contains(&work.query_index))
                {
                    if query_failures.contains_key(&work.query_index) {
                        continue;
                    }
                    let query = plan
                        .queries
                        .get(work.query_index)
                        .ok_or(JamIndexBatchError::Overflow)?;
                    if fallback_results
                        .get(&work.query_index)
                        .is_some_and(|result| complete(result, query.input.sequence.len()))
                    {
                        continue;
                    }
                    let archive = super::archive::JamIndexArchive::from_loaded(
                        reader,
                        group.metagenome_local_id,
                        &chunk,
                        chunk.keys().copied(),
                        0,
                    )?;
                    let candidate = candidate_result(
                        reader,
                        &work.candidate,
                        u64::try_from(query.prepared.unique_hash_count())
                            .map_err(|_| JamIndexBatchError::Overflow)?,
                    )?;
                    let result = match run_archive(
                        &query.input.query_id,
                        &query.input.sequence,
                        candidate,
                        &archive,
                        &config.trace.runner,
                    ) {
                        Ok(result) => result,
                        Err(error) => {
                            query_failures.insert(work.query_index, error.to_string());
                            continue;
                        }
                    };
                    let merged = match fallback_results.remove(&work.query_index) {
                        Some(previous) => match merge_index_results(
                            previous,
                            result,
                            query.input.sequence.len(),
                            &config.trace.runner,
                        ) {
                            Ok(merged) => merged,
                            Err(error) => {
                                query_failures.insert(work.query_index, error.to_string());
                                continue;
                            }
                        },
                        None => result,
                    };
                    fallback_results.insert(work.query_index, merged);
                }
                Ok(())
            },
        )?
    };
    metrics.source_open_count = metrics.source_open_count.saturating_add(1);
    metrics.source_bytes_read = metrics.source_bytes_read.saturating_add(source_bytes);
    crate::profiling::add_counter("gzip_source_bytes_read", source_bytes);
    drop(_source_profile);
    for work in &group.work {
        if query_failures.contains_key(&work.query_index) {
            continue;
        }
        let query = plan
            .queries
            .get(work.query_index)
            .ok_or(JamIndexBatchError::Overflow)?;
        let outcome = if fallback.contains(&work.query_index) {
            let result = fallback_results
                .remove(&work.query_index)
                .ok_or_else(|| JamIndexBatchError::NoResult(query.input.query_id.clone()));
            result.map(|result| {
                (
                    result,
                    JamIndexTraceMetrics {
                        selected_candidates: 1,
                        sequential_scans: 1,
                        ..JamIndexTraceMetrics::default()
                    },
                )
            })
        } else {
            run_plan_loaded(
                reader,
                &query.input.query_id,
                &query.input.sequence,
                u64::try_from(query.prepared.unique_hash_count())
                    .map_err(|_| JamIndexBatchError::Overflow)?,
                &work.candidate,
                &work.contig_plan,
                &config.trace.runner,
                config.trace.expansion_batch,
                &loaded,
            )
            .map_err(JamIndexBatchError::from)
        };
        let (mut result, trace_metrics) = match outcome {
            Ok(result) => result,
            Err(error) => {
                query_failures.insert(work.query_index, error.to_string());
                continue;
            }
        };
        result.run_id = run_id.to_string();
        let alignment_count = u64::try_from(result.alignments.len()).unwrap_or(u64::MAX);
        results.insert(work.query_index, result);
        metrics.results_emitted = metrics.results_emitted.saturating_add(1);
        metrics.alignments_emitted = metrics.alignments_emitted.saturating_add(alignment_count);
        let _ = trace_metrics;
    }
    if !query_failures.is_empty() {
        metrics.failed_groups = metrics.failed_groups.saturating_add(1);
    }
    let mut outcomes = Vec::with_capacity(group.work.len());
    for work in &group.work {
        let outcome = if let Some(result) = results.remove(&work.query_index) {
            Ok(result)
        } else {
            Err(query_failures.remove(&work.query_index).unwrap_or_else(|| {
                JamIndexBatchError::NoResult(plan.queries.get(work.query_index).map_or_else(
                    || "unknown".to_string(),
                    |query| query.input.query_id.clone(),
                ))
                .to_string()
            }))
        };
        outcomes.push((work.query_index, outcome));
    }
    Ok(JamIndexBatchGroupExecution { outcomes, metrics })
}

#[derive(Debug, Error)]
pub enum JamIndexBatchError {
    #[error("invalid Jam Index batch input: {0}")]
    InvalidInput(String),
    #[error("unknown Jam Index part {0}")]
    UnknownPart(u32),
    #[error("Jam Index batch candidate and contig plan do not match")]
    CandidateBinding,
    #[error("Jam Index batch coordinate overflow")]
    Overflow,
    #[error("Jam Index query batch has {observed} {unit}, above limit {limit}")]
    QueryBatchLimit {
        observed: u64,
        limit: u64,
        unit: &'static str,
    },
    #[error(
        "Jam Index batch group {metagenome_id} selected {selected_bases} bases above limit {limit}"
    )]
    GroupTooLarge {
        metagenome_id: String,
        selected_bases: u64,
        limit: u64,
    },
    #[error("Jam Index batch produced no result for query {0}")]
    NoResult(String),
    #[error("Jam Index batch output failed: {0}")]
    Output(String),
    #[error(transparent)]
    Build(#[from] JamIndexBuildError),
    #[error(transparent)]
    Screen(#[from] JamIndexScreenError),
    #[error(transparent)]
    Contig(#[from] JamIndexContigSearchError),
    #[error(transparent)]
    Part(#[from] JamIndexPartError),
    #[error(transparent)]
    Archive(#[from] ArchiveError),
    #[error(transparent)]
    Trace(#[from] JamIndexTraceError),
    #[error(transparent)]
    Runner(#[from] RunnerError),
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::jam_index::{
        JamIndexBuildConfig, JamIndexContigPlan, MetagenomeSource, build_jam_index,
    };
    use crate::trace::config::SensitivityProfile;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use sha2::{Digest, Sha256};
    use std::io::Write;

    fn sequence(length: usize) -> Vec<u8> {
        let mut state = 0x1234_5678_9abc_def0_u64;
        (0..length)
            .map(|_| {
                state ^= state << 7;
                state ^= state >> 9;
                state ^= state << 8;
                b"ACGT"[(state as usize) & 3]
            })
            .collect()
    }

    #[test]
    fn shared_fallback_stream_opens_ordinary_gzip_once() {
        let sources = tempfile::Builder::new()
            .prefix("batch-fallback-source-")
            .tempdir_in("target")
            .unwrap();
        let outputs = tempfile::Builder::new()
            .prefix("batch-fallback-output-")
            .tempdir_in("target")
            .unwrap();
        let query_sequence = sequence(160);
        let source = sources.path().join("sample.fa.gz");
        let file = std::fs::File::create(&source).unwrap();
        let mut writer = GzEncoder::new(file, Compression::fast());
        writeln!(writer, ">noise\n{}", "A".repeat(200)).unwrap();
        writeln!(
            writer,
            ">target\n{}",
            String::from_utf8_lossy(&query_sequence)
        )
        .unwrap();
        writeln!(writer, ">unused\n{}", "C".repeat(200)).unwrap();
        writer.finish().unwrap();
        let root = outputs.path().join("index");
        build_jam_index(
            &root,
            &[MetagenomeSource {
                metagenome_id: "sample".to_string(),
                sequence_path: source,
            }],
            &JamIndexBuildConfig {
                source_manifest_sha256: "a".repeat(64),
                ..JamIndexBuildConfig::default()
            },
        )
        .unwrap();
        let manifest = load_manifest(&root).unwrap();
        let part = &manifest.parts[0];
        let reader =
            JamIndexPartReader::open(root.join(&part.directory).join(&part.data_file)).unwrap();
        let prepared = prepare_screen_query(&query_sequence, &manifest.selection_policy).unwrap();
        let candidate = JamIndexCandidate {
            part_id: 0,
            metagenome_local_id: 0,
            metagenome_id: "sample".to_string(),
            rank: 1,
            shared_hash_count: 4,
            supported_query_windows: 1,
            longest_supported_window_run: 1,
            weighted_hash_sum: 1.0,
            shared_spatial_signatures: 4,
            rare_shared_signatures: 0,
            shared_whole_sample_signatures: 0,
            screen_policy: "test".to_string(),
            admission_source: crate::trace::model::CandidateAdmissionSource::Standard,
            shared_hashes: Vec::new(),
        };
        let contig_plan = JamIndexContigPlan {
            candidate_rank: 1,
            part_id: 0,
            metagenome_local_id: 0,
            metagenome_id: "sample".to_string(),
            ranked_contigs: Vec::new(),
            initial_contig_count: 0,
            sequential_fallback_range: Some(0..3),
            contig_truncated: false,
            shared_spatial_signatures: 4,
            rare_shared_signatures: 0,
            shared_whole_sample_signatures: 0,
            candidate_entry_reason: crate::trace::model::CandidateAdmissionSource::Standard,
        };
        let input = JamIndexBatchQuery {
            query_id: "query".to_string(),
            original_header: "query retained".to_string(),
            sequence_sha256: format!("{:x}", Sha256::digest(&query_sequence)),
            sequence: query_sequence.clone(),
            query_kind: QueryKind::Plasmid,
            topology_requested: TopologyRequested::Linear,
        };
        let planned = JamIndexPlannedQuery {
            input,
            prepared,
            candidates: vec![candidate.clone()],
            candidate_limit_reached: false,
            screen_metrics: JamIndexScreenMetrics::default(),
            contig_plans: vec![contig_plan.clone()],
            contig_metrics: JamIndexContigSearchMetrics::default(),
        };
        let mut readers = BTreeMap::new();
        readers.insert(0, reader);
        let plan = JamIndexBatchPlan {
            queries: vec![planned],
            groups: vec![JamIndexBatchWorkGroup {
                part_id: 0,
                metagenome_local_id: 0,
                metagenome_id: "sample".to_string(),
                source_path: readers[&0].metagenomes()[0].source_path.clone(),
                work: vec![JamIndexBatchWork {
                    query_index: 0,
                    candidate,
                    contig_plan,
                }],
            }],
            screen_parts_opened: 1,
            contig_parts_opened: 1,
            readers,
        };
        let config = JamIndexBatchConfig {
            trace: JamIndexTraceConfig {
                runner: crate::trace::runner::TraceRunnerConfig {
                    sensitivity: super::super::trace::dense_profile(SensitivityProfile::Sensitive),
                    topology_requested: TopologyRequested::Linear,
                    query_kind: QueryKind::Plasmid,
                    ..crate::trace::runner::TraceRunnerConfig::default()
                },
                expansion_batch: 1,
                ..JamIndexTraceConfig::default()
            },
            max_group_contig_bases: 1024 * 1024,
            fallback_contigs_per_chunk: 1,
            ..JamIndexBatchConfig::default()
        };
        let mut results = Vec::new();
        let execution = execute_batch(&plan, &config, "fallback-run", |result| {
            results.push(result);
            Ok(())
        })
        .unwrap();
        assert_eq!(execution.metrics.source_open_count, 1);
        assert_eq!(execution.metrics.sequential_source_open_count, 1);
        assert_eq!(results.len(), 1);
        assert_eq!(
            results[0].coverage.as_ref().unwrap().supported_bases,
            query_sequence.len() as u64
        );
    }
}
