//! Three-stage local Jam Index trace orchestration.

use super::archive::JamIndexArchive;
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
use crate::trace::model::{CandidateResult, TraceMetagenomeResult};
use crate::trace::runner::{RunnerError, TraceRunnerConfig, run_index_archive};
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
}

impl Default for JamIndexTraceConfig {
    fn default() -> Self {
        Self {
            screen: JamIndexScreenConfig::default(),
            contigs: JamIndexContigSearchConfig::default(),
            runner: TraceRunnerConfig::default(),
            expansion_batch: 8,
        }
    }
}

impl JamIndexTraceConfig {
    fn validate(&self) -> Result<(), JamIndexTraceError> {
        if self.expansion_batch == 0 {
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
    let mut metagenomes = Vec::with_capacity(contigs.plans.len());
    for plan in &contigs.plans {
        let candidate = candidates
            .get(&plan.candidate_rank)
            .copied()
            .ok_or(JamIndexTraceError::CandidateBinding)?;
        let part = manifest
            .parts
            .get(usize::try_from(plan.part_id).map_err(|_| JamIndexTraceError::Overflow)?)
            .ok_or(JamIndexTraceError::CandidateBinding)?;
        let reader = JamIndexPartReader::open(root.join(&part.directory).join(&part.data_file))?;
        let (result, delta) = run_plan(
            &reader,
            query_id,
            query,
            query_hashes,
            candidate,
            plan,
            &runner,
            config.expansion_batch,
        )?;
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
    let mut sequential = false;
    if selected.is_empty()
        && let Some(range) = &plan.sequential_fallback_range
    {
        selected.extend(range.clone());
        sequential = true;
    }
    if selected.is_empty() {
        return Err(JamIndexTraceError::NoContigs(plan.metagenome_id.clone()));
    }
    let candidate_result = candidate_result(reader, candidate, query_hashes)?;
    let mut next =
        usize::try_from(plan.initial_contig_count).map_err(|_| JamIndexTraceError::Overflow)?;
    let mut passes = 0u64;
    let mut expanded = false;
    let result = loop {
        let archive = JamIndexArchive::load(reader, plan.metagenome_local_id, selected.clone())?;
        let result =
            run_index_archive(query_id, query, candidate_result.clone(), &archive, runner)?;
        passes = passes.saturating_add(1);
        if complete(&result, query.len()) || next >= plan.ranked_contigs.len() || sequential {
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
            sequential_scans: u64::from(sequential),
        },
    ))
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
