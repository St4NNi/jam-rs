//! Three-stage local Jam Index trace orchestration.

use super::archive::{ExactContigMatch, JamIndexArchive};
use super::builder::{JamIndexBuildError, load_manifest};
use super::contig_search::{
    JamIndexContigPlan, JamIndexContigSearchConfig, JamIndexContigSearchError,
    JamIndexContigSearchMetrics, select_candidate_contigs,
};
use super::part::JamIndexPartReader;
use super::screen::{
    JamIndexCandidate, JamIndexScreenConfig, JamIndexScreenError, JamIndexScreenMetrics,
    prepare_screen_query, search_jam_index,
};
use crate::trace::config::SensitivityProfile;
use crate::trace::model::{
    AlignmentRole, BaseAlignment, BaseInterval, CandidateResult, EditOperation, EditRun,
    SeedEvidence, TopologyRequested, TraceMetagenomeResult,
};
use crate::trace::runner::{
    RunnerError, TraceRunnerConfig, merge_index_results, run_exact_index, run_index_archive,
};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::Path;
use thiserror::Error;

#[derive(Clone, Debug)]
pub struct JamIndexTraceConfig {
    pub screen: JamIndexScreenConfig,
    pub contigs: JamIndexContigSearchConfig,
    pub runner: TraceRunnerConfig,
    pub expansion_batch: usize,
    pub parallel_candidates: usize,
}

impl Default for JamIndexTraceConfig {
    fn default() -> Self {
        Self {
            screen: JamIndexScreenConfig::default(),
            contigs: JamIndexContigSearchConfig::default(),
            runner: TraceRunnerConfig::default(),
            expansion_batch: 8,
            parallel_candidates: 1,
        }
    }
}

impl JamIndexTraceConfig {
    fn validate(&self) -> Result<(), JamIndexTraceError> {
        if self.expansion_batch == 0 || self.parallel_candidates == 0 {
            return Err(JamIndexTraceError::InvalidConfig);
        }
        self.runner.validate()?;
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub struct JamIndexTraceMetrics {
    pub selected_candidates: u64,
    pub selected_contigs: u64,
    pub selected_contig_bases: u64,
    pub alignment_passes: u64,
    pub expanded_candidates: u64,
    pub sequential_scans: u64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct JamIndexTraceOutput {
    pub screen_metrics: JamIndexScreenMetrics,
    pub contig_metrics: JamIndexContigSearchMetrics,
    pub trace_metrics: JamIndexTraceMetrics,
    pub metagenomes: Vec<TraceMetagenomeResult>,
}

pub fn trace_index(
    root: impl AsRef<Path>,
    query_id: &str,
    query: &[u8],
    config: &JamIndexTraceConfig,
) -> Result<JamIndexTraceOutput, JamIndexTraceError> {
    config.validate()?;
    let root = root.as_ref();
    let manifest = load_manifest(root)?;
    let prepared = prepare_screen_query(query, &manifest.selection_policy)?;
    let screen = search_jam_index(root, &prepared, config.screen)?;
    let contigs = select_candidate_contigs(root, &prepared, &screen.candidates, config.contigs)?;
    let candidates = screen
        .candidates
        .iter()
        .map(|candidate| (candidate.rank, candidate))
        .collect::<BTreeMap<_, _>>();
    let mut runner = config.runner.clone();
    runner.sensitivity = dense_profile(runner.sensitivity.profile);
    runner.topology_margin_bases = runner.sensitivity.auto_topology_margin_bases;
    runner.validate()?;
    let query_hashes =
        u64::try_from(prepared.unique_hash_count()).map_err(|_| JamIndexTraceError::Overflow)?;
    let mut metrics = JamIndexTraceMetrics {
        selected_candidates: u64::try_from(screen.candidates.len()).unwrap_or(u64::MAX),
        ..JamIndexTraceMetrics::default()
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.parallel_candidates)
        .build()
        .map_err(|_| JamIndexTraceError::InvalidConfig)?;
    let completed = pool.install(|| {
        contigs
            .plans
            .par_iter()
            .map(|plan| {
                let candidate = candidates
                    .get(&plan.candidate_rank)
                    .copied()
                    .ok_or(JamIndexTraceError::CandidateBinding)?;
                let part = manifest
                    .parts
                    .get(usize::try_from(plan.part_id).map_err(|_| JamIndexTraceError::Overflow)?)
                    .ok_or(JamIndexTraceError::CandidateBinding)?;
                let reader =
                    JamIndexPartReader::open(root.join(&part.directory).join(&part.data_file))?;
                run_plan(
                    &reader,
                    query_id,
                    query,
                    query_hashes,
                    candidate,
                    plan,
                    &runner,
                    config.expansion_batch,
                )
            })
            .collect::<Result<Vec<_>, JamIndexTraceError>>()
    })?;
    let mut metagenomes = Vec::with_capacity(completed.len());
    for (result, delta) in completed {
        metrics.selected_contigs = metrics
            .selected_contigs
            .saturating_add(delta.selected_contigs);
        metrics.selected_contig_bases = metrics
            .selected_contig_bases
            .saturating_add(delta.selected_contig_bases);
        metrics.alignment_passes = metrics
            .alignment_passes
            .saturating_add(delta.alignment_passes);
        metrics.expanded_candidates = metrics
            .expanded_candidates
            .saturating_add(delta.expanded_candidates);
        metrics.sequential_scans = metrics
            .sequential_scans
            .saturating_add(delta.sequential_scans);
        metagenomes.push(result);
    }
    metagenomes.sort_by(|left, right| left.metagenome_id.cmp(&right.metagenome_id));
    Ok(JamIndexTraceOutput {
        screen_metrics: screen.metrics,
        contig_metrics: contigs.metrics,
        trace_metrics: metrics,
        metagenomes,
    })
}

#[allow(clippy::too_many_arguments)]
fn run_plan(
    reader: &JamIndexPartReader,
    query_id: &str,
    query: &[u8],
    query_hashes: u64,
    candidate: &JamIndexCandidate,
    plan: &JamIndexContigPlan,
    runner: &TraceRunnerConfig,
    expansion_batch: usize,
) -> Result<(TraceMetagenomeResult, JamIndexTraceMetrics), JamIndexTraceError> {
    let mut selected = plan
        .initial_contigs()
        .iter()
        .map(|contig| contig.contig_id)
        .collect::<Vec<_>>();
    if selected.is_empty() {
        if let Some(range) = &plan.sequential_fallback_range {
            return run_fallback(
                reader,
                query_id,
                query,
                query_hashes,
                candidate,
                range.clone(),
                runner,
                expansion_batch,
            );
        }
        return Err(JamIndexTraceError::NoContigs(plan.metagenome_id.clone()));
    }
    let candidate_result = candidate_result(reader, candidate, query_hashes)?;
    let mut next =
        usize::try_from(plan.initial_contig_count).map_err(|_| JamIndexTraceError::Overflow)?;
    let mut passes = 0u64;
    let mut expanded = false;
    let result = loop {
        let archive = JamIndexArchive::load(reader, plan.metagenome_local_id, selected.clone())?;
        let result = run_archive(query_id, query, candidate_result.clone(), &archive, runner)?;
        passes = passes.saturating_add(1);
        if complete(&result, query.len()) || next >= plan.ranked_contigs.len() {
            break result;
        }
        let end = next
            .saturating_add(expansion_batch)
            .min(plan.ranked_contigs.len());
        selected.extend(
            plan.ranked_contigs[next..end]
                .iter()
                .map(|contig| contig.contig_id),
        );
        next = end;
        expanded = true;
    };
    let selected_contig_bases = selected.iter().try_fold(0u64, |total, contig_id| {
        let contig = reader
            .contigs()
            .get(usize::try_from(*contig_id).map_err(|_| JamIndexTraceError::Overflow)?)
            .ok_or(JamIndexTraceError::CandidateBinding)?;
        Ok::<_, JamIndexTraceError>(total.saturating_add(contig.base_count))
    })?;
    Ok((
        result,
        JamIndexTraceMetrics {
            selected_candidates: 1,
            selected_contigs: u64::try_from(selected.len()).unwrap_or(u64::MAX),
            selected_contig_bases,
            alignment_passes: passes,
            expanded_candidates: u64::from(expanded),
            sequential_scans: 0,
        },
    ))
}

#[allow(clippy::too_many_arguments)]
fn run_fallback(
    reader: &JamIndexPartReader,
    query_id: &str,
    query: &[u8],
    query_hashes: u64,
    candidate: &JamIndexCandidate,
    range: std::ops::Range<u32>,
    runner: &TraceRunnerConfig,
    batch_size: usize,
) -> Result<(TraceMetagenomeResult, JamIndexTraceMetrics), JamIndexTraceError> {
    let candidate_result = candidate_result(reader, candidate, query_hashes)?;
    let batch_size = u32::try_from(batch_size).unwrap_or(u32::MAX);
    let mut start = range.start;
    let mut merged = None;
    let mut passes = 0u64;
    let mut selected_contigs = 0u64;
    let mut selected_contig_bases = 0u64;
    while start < range.end {
        let end = start.saturating_add(batch_size).min(range.end);
        let contig_ids = start..end;
        for contig_id in contig_ids.clone() {
            let contig = reader
                .contigs()
                .get(usize::try_from(contig_id).map_err(|_| JamIndexTraceError::Overflow)?)
                .ok_or(JamIndexTraceError::CandidateBinding)?;
            selected_contig_bases = selected_contig_bases.saturating_add(contig.base_count);
        }
        selected_contigs = selected_contigs.saturating_add(u64::from(end - start));
        let archive = JamIndexArchive::load(reader, candidate.metagenome_local_id, contig_ids)?;
        let result = run_archive(query_id, query, candidate_result.clone(), &archive, runner)?;
        passes = passes.saturating_add(1);
        let result = match merged.take() {
            Some(previous) => merge_index_results(previous, result, query.len(), runner)?,
            None => result,
        };
        let finished = complete(&result, query.len());
        merged = Some(result);
        if finished {
            break;
        }
        start = end;
    }
    let result =
        merged.ok_or_else(|| JamIndexTraceError::NoContigs(candidate.metagenome_id.clone()))?;
    Ok((
        result,
        JamIndexTraceMetrics {
            selected_candidates: 1,
            selected_contigs,
            selected_contig_bases,
            alignment_passes: passes,
            expanded_candidates: 0,
            sequential_scans: 1,
        },
    ))
}

fn run_archive(
    query_id: &str,
    query: &[u8],
    candidate: CandidateResult,
    archive: &JamIndexArchive,
    runner: &TraceRunnerConfig,
) -> Result<TraceMetagenomeResult, JamIndexTraceError> {
    let allow_wrap = !matches!(runner.topology_requested, TopologyRequested::Linear);
    let exact = archive.exact_matches(query, allow_wrap);
    if exact.is_empty() {
        return Ok(run_index_archive(
            query_id, query, candidate, archive, runner,
        )?);
    }
    let alignments = exact
        .into_iter()
        .map(|item| {
            exact_alignment(
                query_id,
                &candidate.metagenome_id,
                item,
                runner.sensitivity.alignment.match_score,
            )
        })
        .collect::<Vec<_>>();
    let metrics = crate::archive::TraceArchive::metrics(archive);
    let contigs_considered = u64::try_from(
        crate::archive::TraceArchive::metadata(archive)?
            .contigs
            .len(),
    )
    .unwrap_or(u64::MAX);
    Ok(run_exact_index(
        query_id,
        query.len(),
        candidate,
        alignments,
        metrics,
        contigs_considered,
        runner,
    )?)
}

fn exact_alignment(
    query_id: &str,
    metagenome_id: &str,
    exact: ExactContigMatch,
    match_score: i32,
) -> BaseAlignment {
    let query_segments = if exact.query_start == 0 {
        vec![BaseInterval {
            start: 0,
            end: exact.length,
        }]
    } else {
        vec![
            BaseInterval {
                start: exact.query_start,
                end: exact.length,
            },
            BaseInterval {
                start: 0,
                end: exact.query_start,
            },
        ]
    };
    let score =
        i64::from(match_score).saturating_mul(i64::try_from(exact.length).unwrap_or(i64::MAX));
    BaseAlignment {
        alignment_id: String::new(),
        plasmid_id: query_id.to_string(),
        metagenome_id: metagenome_id.to_string(),
        contig_id: exact.contig_id,
        strand: exact.strand,
        query_segments,
        target_interval: BaseInterval {
            start: exact.target_start,
            end: exact.target_start.saturating_add(exact.length),
        },
        query_length: exact.length,
        target_length: exact.length,
        origin_crossing: exact.query_start != 0,
        score,
        matches: exact.length,
        substitutions: 0,
        insertions: 0,
        deletions: 0,
        cigar: format!("{}=", exact.length),
        edit_script: equal_runs(exact.length),
        chain_score: score,
        identity: 1.0,
        seed_evidence: SeedEvidence::default(),
        primary_supported_bases: 0,
        secondary_supported_bases: 0,
        newly_supported_bases: 0,
        role: AlignmentRole::AlternativeMapping,
        primary: false,
    }
}

fn equal_runs(mut length: u64) -> Vec<EditRun> {
    let mut runs = Vec::new();
    while length != 0 {
        let chunk = length.min(u64::from(u32::MAX));
        runs.push(EditRun {
            operation: EditOperation::Equal,
            length: u32::try_from(chunk).unwrap_or(u32::MAX),
        });
        length -= chunk;
    }
    runs
}

fn candidate_result(
    reader: &JamIndexPartReader,
    candidate: &JamIndexCandidate,
    query_hashes: u64,
) -> Result<CandidateResult, JamIndexTraceError> {
    let metagenome = reader
        .metagenomes()
        .get(
            usize::try_from(candidate.metagenome_local_id)
                .map_err(|_| JamIndexTraceError::Overflow)?,
        )
        .ok_or(JamIndexTraceError::CandidateBinding)?;
    Ok(CandidateResult {
        metagenome_id: candidate.metagenome_id.clone(),
        shared_hashes: u64::from(candidate.shared_hash_count),
        plasmid_hashes: query_hashes,
        metagenome_hashes: u64::from(metagenome.screen_hash_count),
        plasmid_containment: if query_hashes == 0 {
            0.0
        } else {
            f64::from(candidate.shared_hash_count) / query_hashes as f64
        },
        metagenome_containment: if metagenome.screen_hash_count == 0 {
            0.0
        } else {
            f64::from(candidate.shared_hash_count) / f64::from(metagenome.screen_hash_count)
        },
        rank: candidate.rank,
        score_mode: "jam_index".to_string(),
        bias_weighted_plasmid_containment: None,
        uniform_hash_e_value: None,
    })
}

fn dense_profile(profile: SensitivityProfile) -> crate::trace::config::SensitivityConfig {
    let mut sensitivity = crate::trace::config::SensitivityConfig::for_profile(profile);
    sensitivity.primary.scale = 1;
    if let Some(rescue) = sensitivity.rescue.as_mut() {
        rescue.scale = 1;
    }
    sensitivity.gap_rescue.dense_primary = None;
    sensitivity
}

fn complete(result: &TraceMetagenomeResult, query_length: usize) -> bool {
    result
        .coverage
        .as_ref()
        .is_some_and(|coverage| coverage.supported_bases >= query_length as u64)
}

#[derive(Debug, Error)]
pub enum JamIndexTraceError {
    #[error("invalid Jam Index trace config")]
    InvalidConfig,
    #[error("Jam Index candidate and part binding mismatch")]
    CandidateBinding,
    #[error("Jam Index candidate {0} has no selected contigs")]
    NoContigs(String),
    #[error("Jam Index trace coordinate overflow")]
    Overflow,
    #[error(transparent)]
    Build(#[from] JamIndexBuildError),
    #[error(transparent)]
    Screen(#[from] JamIndexScreenError),
    #[error(transparent)]
    Contig(#[from] JamIndexContigSearchError),
    #[error(transparent)]
    Part(#[from] super::part::JamIndexPartError),
    #[error(transparent)]
    Archive(#[from] crate::archive::ArchiveError),
    #[error(transparent)]
    Runner(#[from] RunnerError),
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::jam_index::{JamIndexBuildConfig, MetagenomeSource, build_jam_index};
    use crate::trace::config::SensitivityProfile;
    use std::fs;
    use tempfile::Builder;

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
    fn fallback_scans_batches() {
        let sources = Builder::new()
            .prefix("index-fallback-source-")
            .tempdir_in("target")
            .unwrap();
        let outputs = Builder::new()
            .prefix("index-fallback-output-")
            .tempdir_in("target")
            .unwrap();
        let query = sequence(160);
        let source = sources.path().join("sample.fasta");
        fs::write(
            &source,
            format!(
                ">noise\n{}\n>target\n{}\n>unused\n{}\n",
                "A".repeat(200),
                String::from_utf8_lossy(&query),
                "C".repeat(200)
            ),
        )
        .unwrap();
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
        let candidate = JamIndexCandidate {
            part_id: 0,
            metagenome_local_id: 0,
            metagenome_id: "sample".to_string(),
            rank: 1,
            shared_hash_count: 4,
            supported_query_windows: 1,
            longest_supported_window_run: 1,
            weighted_hash_sum: 1.0,
            shared_hashes: Vec::new(),
        };
        let plan = JamIndexContigPlan {
            candidate_rank: 1,
            part_id: 0,
            metagenome_local_id: 0,
            metagenome_id: "sample".to_string(),
            ranked_contigs: Vec::new(),
            initial_contig_count: 0,
            sequential_fallback_range: Some(0..3),
        };
        let runner = TraceRunnerConfig {
            sensitivity: dense_profile(SensitivityProfile::Sensitive),
            ..TraceRunnerConfig::default()
        };
        let (result, metrics) =
            run_plan(&reader, "query", &query, 4, &candidate, &plan, &runner, 1).unwrap();
        assert_eq!(metrics.sequential_scans, 1);
        assert_eq!(metrics.alignment_passes, 2);
        assert_eq!(metrics.selected_contigs, 2);
        assert_eq!(
            result.coverage.as_ref().unwrap().supported_bases,
            query.len() as u64
        );
    }
}
