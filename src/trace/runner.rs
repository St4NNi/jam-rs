//! Candidate-to-sequence trace runner.
//!
//! The runner owns orchestration and bounded resource use.  It queries the
//! existing `.jam` candidate index once, resolves candidate resources through
//! the catalog, retrieves indexed windows when JMA is available, and invokes
//! the alignment kernel.  Its output remains evidence for downstream
//! confirmation; it never upgrades a candidate to a biological presence call.

#[path = "cascade/mod.rs"]
pub mod cascade;

use crate::archive::{
    ArchiveError, ArchiveMetrics, NativeJmaArchive, SeedKey, SeedSchemeDescriptor, SeedSchemeId,
    SequenceRequest, TraceArchive,
};
use crate::jma::{ArchiveReader, JmaError, SequenceRange};
use crate::resource::{RangeReader, ResourceError, ResourceMetrics, ResourceOpenOptions};
use crate::router::RoutedCandidate;
use crate::trace::TraceMemoryEstimate;
use crate::trace::alignment::exact_blocks::extend_and_merge_anchors;
use crate::trace::alignment::{
    AlignmentCorridor, AlignmentError, AlignmentOptions, AlignmentRetryMetadata,
    AlignmentWorkspace, CorridorAnchor,
};
use crate::trace::anchors::{
    Anchor, AnchorOccurrenceEvidence, AnchorSet, SeedOccurrenceGroup, anchor_strand,
    summarize_anchor_evidence,
};
use crate::trace::catalog::{CatalogEntry, CatalogError, TraceCatalog};
use crate::trace::chain::{
    AnchorChain, AnchorClass, ChainConfig, ChainError, WeightedAnchor, WeightedAnchorChain,
    chain_weighted_anchors,
};
use crate::trace::config::{SeedSensitivity, SensitivityConfig};
use crate::trace::memory::TraceMemorySemaphore;
use crate::trace::model::{
    AlignmentRole, BaseAlignment, BaseInterval, CandidatePerformanceCounters, CandidateResult,
    CoordinateModel, CoverageSummary, FragmentMosaicSummary, InputResource, QueryKind,
    RescueRoundMetrics, SeedEvidence, Strand, TopologyAssessment, TopologyEvidence,
    TopologyRequested, TraceFailure, TraceMetagenomeResult, TraceStageMetrics, TraceStatus,
};
use crate::trace::mosaic::{MosaicError, assess_topology, assign_alignment_ids};
use crate::trace::raw::{AssemblyResource, RawAssembly, RawError, open_resource};
use crate::trace::screen::{
    CandidateError, CandidateScoreMode, CandidateSearchConfig, CandidateSearchResult,
    CandidateSearcher, RankedCandidate,
};
use crate::trace::seeds::gear::{
    ExactFragment, FragmentOrientation, FragmentationMode, GearConfig, GearTableKind,
    VerifiedFragment, fragment_bytes, fragment_sequence, merge_exact_runs, verify_exact_fragment,
};
use crate::trace::seeds::{
    QuerySeed, SeedError, extract_seed_level, extract_seed_level_in_intervals,
};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;
use thiserror::Error;

/// Runner configuration.  All resource and candidate limits are explicit so
/// a run can be reproduced from its JSON metadata.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TraceRunnerConfig {
    pub sensitivity: SensitivityConfig,
    pub candidates: CandidateSearchConfig,
    pub resources: ResourceOpenOptions,
    pub threads: usize,
    pub io_concurrency: usize,
    pub max_alignments_per_candidate: usize,
    pub query_kind: QueryKind,
    pub topology_requested: TopologyRequested,
    pub topology_margin_bases: u64,
    #[serde(default = "default_memory_budget_bytes")]
    pub memory_budget_bytes: u64,
}

const fn default_memory_budget_bytes() -> u64 {
    1_400 * 1024 * 1024
}

impl Default for TraceRunnerConfig {
    fn default() -> Self {
        let resources = ResourceOpenOptions {
            allow_full_download_fallback: false,
            ..ResourceOpenOptions::default()
        };
        Self {
            sensitivity: SensitivityConfig::default(),
            candidates: CandidateSearchConfig::default(),
            resources,
            threads: 1,
            io_concurrency: 1,
            max_alignments_per_candidate: 256,
            query_kind: QueryKind::Unknown,
            topology_requested: TopologyRequested::Auto,
            topology_margin_bases: SensitivityConfig::default().auto_topology_margin_bases,
            memory_budget_bytes: default_memory_budget_bytes(),
        }
    }
}

impl TraceRunnerConfig {
    pub fn validate(&self) -> Result<(), RunnerError> {
        self.sensitivity
            .validate()
            .map_err(|error| RunnerError::InvalidConfig(error.to_string()))?;
        self.candidates.validate()?;
        if self.threads == 0 {
            return Err(RunnerError::InvalidConfig(
                "threads must be greater than zero".to_string(),
            ));
        }
        if self.io_concurrency == 0 {
            return Err(RunnerError::InvalidConfig(
                "io_concurrency must be greater than zero".to_string(),
            ));
        }
        if self.max_alignments_per_candidate == 0 {
            return Err(RunnerError::InvalidConfig(
                "max_alignments_per_candidate must be greater than zero".to_string(),
            ));
        }
        if self.topology_margin_bases == 0 {
            return Err(RunnerError::InvalidConfig(
                "topology_margin_bases must be greater than zero".to_string(),
            ));
        }
        if self.memory_budget_bytes == 0 {
            return Err(RunnerError::InvalidConfig(
                "memory_budget_bytes must be greater than zero".to_string(),
            ));
        }
        if self.resources.cache_block_bytes == 0 || self.resources.max_cache_bytes == 0 {
            return Err(RunnerError::InvalidConfig(
                "resource cache limits must be greater than zero".to_string(),
            ));
        }
        Ok(())
    }
}

/// Input to one trace run.  The sequence is exactly one plasmid record.
#[derive(Clone, Debug)]
pub struct TraceQuery {
    pub plasmid_id: String,
    pub plasmid_sequence: Vec<u8>,
    pub database: String,
    pub catalog: TraceCatalog,
}

impl TraceQuery {
    pub fn from_fasta(
        plasmid_id: impl Into<String>,
        path: impl AsRef<Path>,
        database: impl AsRef<Path>,
        catalog: TraceCatalog,
        resources: ResourceOpenOptions,
    ) -> Result<Self, RunnerError> {
        let assembly = RawAssembly::open(path.as_ref().to_string_lossy(), resources)?;
        if assembly.contigs.len() != 1 {
            return Err(RunnerError::InvalidQuery(format!(
                "plasmid FASTA must contain exactly one record, got {}",
                assembly.contigs.len()
            )));
        }
        Ok(Self {
            plasmid_id: plasmid_id.into(),
            plasmid_sequence: assembly.contigs[0].sequence.clone(),
            database: database.as_ref().to_string_lossy().into_owned(),
            catalog,
        })
    }

    pub fn validate(&self) -> Result<(), RunnerError> {
        if self.plasmid_id.trim().is_empty() {
            return Err(RunnerError::InvalidQuery(
                "plasmid_id must not be empty".to_string(),
            ));
        }
        if self.plasmid_sequence.is_empty() {
            return Err(RunnerError::InvalidQuery(
                "plasmid sequence must not be empty".to_string(),
            ));
        }
        if self.database.trim().is_empty() {
            return Err(RunnerError::InvalidQuery(
                "candidate database path must not be empty".to_string(),
            ));
        }
        Ok(())
    }
}

/// Bounded counters from candidate search and positional processing.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub struct TracePerformanceCounters {
    pub candidate_queries: u64,
    pub candidates_considered: u64,
    pub candidates_processed: u64,
    pub contigs_considered: u64,
    pub windows_retrieved: u64,
    pub alignments_attempted: u64,
    pub alignments_succeeded: u64,
    pub failures: u64,
    pub elapsed_millis: u64,
}

impl TracePerformanceCounters {
    fn add_assign(&mut self, other: Self) {
        self.candidates_processed = self
            .candidates_processed
            .saturating_add(other.candidates_processed);
        self.contigs_considered = self
            .contigs_considered
            .saturating_add(other.contigs_considered);
        self.windows_retrieved = self
            .windows_retrieved
            .saturating_add(other.windows_retrieved);
        self.alignments_attempted = self
            .alignments_attempted
            .saturating_add(other.alignments_attempted);
        self.alignments_succeeded = self
            .alignments_succeeded
            .saturating_add(other.alignments_succeeded);
        self.failures = self.failures.saturating_add(other.failures);
    }

    fn candidate_snapshot(self) -> CandidatePerformanceCounters {
        CandidatePerformanceCounters {
            candidates_processed: self.candidates_processed,
            contigs_considered: self.contigs_considered,
            windows_retrieved: self.windows_retrieved,
            alignments_attempted: self.alignments_attempted,
            alignments_succeeded: self.alignments_succeeded,
            failures: self.failures,
            elapsed_millis: self.elapsed_millis,
        }
    }
}

fn candidate_memory_estimate(
    config: &TraceRunnerConfig,
    query_length: usize,
) -> TraceMemoryEstimate {
    let anchors = u64::from(config.sensitivity.max_anchors_per_candidate);
    let chains = u64::from(config.sensitivity.max_chains_per_candidate);
    let query_length_u64 = u64::try_from(query_length).unwrap_or(u64::MAX);
    let band = u64::from(alignment_half_width_cap(config, query_length_u64));
    let row_width = band.saturating_mul(2).saturating_add(1);
    let query_rows = u64::try_from(query_length)
        .unwrap_or(u64::MAX)
        .saturating_add(1);
    TraceMemoryEstimate {
        candidate_accumulator_bytes: u64::from(config.sensitivity.max_candidates)
            .saturating_mul(256),
        router_key_page_bytes: 0,
        router_posting_bytes: 0,
        jma_key_page_bytes: config.resources.cache_block_bytes.saturating_mul(4),
        jma_occurrence_bytes: anchors.saturating_mul(24),
        anchor_bytes: anchors
            .saturating_mul(u64::try_from(std::mem::size_of::<Anchor>()).unwrap_or(u64::MAX)),
        chain_bytes: anchors
            .saturating_mul(
                u64::try_from(std::mem::size_of::<WeightedAnchor>()).unwrap_or(u64::MAX),
            )
            .saturating_add(chains.saturating_mul(256)),
        sequence_block_bytes: config.sensitivity.max_alignment_window_bases,
        alignment_score_bytes: row_width.saturating_mul(2).saturating_mul(32),
        traceback_bytes: query_rows.saturating_mul(row_width).saturating_mul(32),
        output_bytes: u64::try_from(config.max_alignments_per_candidate)
            .unwrap_or(u64::MAX)
            .saturating_mul(8 * 1024),
    }
}

fn alignment_half_width_cap(config: &TraceRunnerConfig, query_length: u64) -> u32 {
    const BUDGET_BYTES_PER_BANDED_POSITION: u64 = 40;

    let row_budget = config.memory_budget_bytes
        / query_length
            .saturating_add(1)
            .saturating_mul(BUDGET_BYTES_PER_BANDED_POSITION)
            .max(1);
    u32::try_from(row_budget.saturating_sub(1) / 2)
        .unwrap_or(u32::MAX)
        .clamp(64, 2048)
}

/// Complete in-memory orchestration result.  The JSON writer consumes
/// `metagenomes` incrementally; callers should not treat this container as a
/// replacement for streaming output on very large runs.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TraceRunOutput {
    pub search: CandidateSearchResult,
    pub metagenomes: Vec<TraceMetagenomeResult>,
    pub counters: TracePerformanceCounters,
    pub database_input: InputResource,
    pub database_metrics: ResourceMetrics,
}

impl TraceRunOutput {
    /// Attach the caller's run identifier before records are handed to the
    /// JSONL writer.  The runner does not invent identifiers because the
    /// header owns the run-level provenance.
    pub fn set_run_id(&mut self, run_id: impl Into<String>) {
        let run_id = run_id.into();
        for result in &mut self.metagenomes {
            result.run_id = run_id.clone();
        }
    }
}

#[derive(Debug)]
struct ResolvedCandidateDatabase {
    path: PathBuf,
    input: InputResource,
    metrics: ResourceMetrics,
}

#[derive(Debug, Deserialize, Serialize)]
struct CandidateDatabaseCacheManifest {
    schema_version: String,
    redacted_locator: String,
    version: Option<String>,
    size: u64,
    sha256: String,
}

fn resolve_candidate_database(
    locator: &str,
    options: &ResourceOpenOptions,
) -> Result<ResolvedCandidateDatabase, RunnerError> {
    let resource = open_resource(locator, options.clone())?;
    let redacted_locator = resource.locator().redacted();
    match resource {
        AssemblyResource::Local(resource) => {
            let path = resource.path().to_path_buf();
            let sha256 = crate::provenance::sha256_file(&path).map_err(|error| {
                RunnerError::DatabaseCache {
                    locator: redacted_locator.clone(),
                    message: format!("unable to checksum local database: {error}"),
                }
            })?;
            Ok(ResolvedCandidateDatabase {
                path,
                input: InputResource {
                    role: "database".to_string(),
                    redacted_locator,
                    sha256: Some(sha256),
                },
                metrics: resource.metrics(),
            })
        }
        AssemblyResource::Object(resource) => {
            let metadata = resource.metadata()?;
            let metrics = resource.metrics();
            resolve_remote_candidate_database(locator, redacted_locator, metadata, metrics, options)
        }
    }
}

fn resolve_remote_candidate_database(
    locator: &str,
    redacted_locator: String,
    metadata: crate::resource::ResourceMetadata,
    mut metrics: ResourceMetrics,
    options: &ResourceOpenOptions,
) -> Result<ResolvedCandidateDatabase, RunnerError> {
    let cache_dir = options
        .cache_dir
        .as_deref()
        .ok_or_else(|| RunnerError::DatabaseCache {
            locator: redacted_locator.clone(),
            message: "remote .jam input requires --cache-dir so it can be memory mapped"
                .to_string(),
        })?;
    fs::create_dir_all(cache_dir).map_err(|error| RunnerError::DatabaseCache {
        locator: redacted_locator.clone(),
        message: format!(
            "unable to create cache directory {}: {error}",
            cache_dir.display()
        ),
    })?;
    let version = database_version(&metadata);
    if let Some(version) = version.as_deref() {
        let stem = database_cache_key(&redacted_locator, version, metadata.size);
        if let Some(resolved) = cached_database(
            cache_dir,
            &stem,
            &redacted_locator,
            Some(version),
            metadata.size,
            metrics,
        )? {
            return Ok(resolved);
        }
    }

    let mut last_error = None;
    let mut downloaded = None;
    for attempt in 0..options.max_retries.saturating_add(1) {
        let resource = match open_resource(locator, options.clone()) {
            Ok(resource) => resource,
            Err(error) => {
                last_error = Some(error.to_string());
                if attempt < options.max_retries {
                    metrics.retries = metrics.retries.saturating_add(1);
                    continue;
                }
                break;
            }
        };
        let mut temporary = tempfile::NamedTempFile::new_in(cache_dir).map_err(|error| {
            RunnerError::DatabaseCache {
                locator: redacted_locator.clone(),
                message: format!("unable to create temporary cache file: {error}"),
            }
        })?;
        let result = resource
            .stream()
            .and_then(|mut stream| {
                io::copy(&mut stream, temporary.as_file_mut()).map_err(|error| ResourceError::Io {
                    locator: redacted_locator.clone(),
                    message: error.to_string(),
                })
            })
            .and_then(|copied| {
                temporary
                    .as_file_mut()
                    .flush()
                    .map_err(|error| ResourceError::Io {
                        locator: redacted_locator.clone(),
                        message: error.to_string(),
                    })?;
                if copied != metadata.size {
                    return Err(ResourceError::Transport {
                        locator: redacted_locator.clone(),
                        message: format!(
                            "downloaded {copied} bytes but metadata declared {}",
                            metadata.size
                        ),
                    });
                }
                Ok(())
            });
        metrics = add_resource_metrics(metrics, resource.metrics());
        match result {
            Ok(()) => {
                downloaded = Some(temporary);
                break;
            }
            Err(error) => {
                last_error = Some(error.to_string());
                if attempt < options.max_retries {
                    metrics.retries = metrics.retries.saturating_add(1);
                }
            }
        }
    }
    let temporary = downloaded.ok_or_else(|| RunnerError::DatabaseCache {
        locator: redacted_locator.clone(),
        message: last_error.unwrap_or_else(|| "remote database download failed".to_string()),
    })?;

    let verification = open_resource(locator, options.clone())?;
    let verified_metadata = verification.metadata()?;
    metrics = add_resource_metrics(metrics, verification.metrics());
    if verified_metadata.size != metadata.size || database_version(&verified_metadata) != version {
        return Err(RunnerError::DatabaseCache {
            locator: redacted_locator,
            message: "remote database identity changed during download".to_string(),
        });
    }
    let sha256 = crate::provenance::sha256_file(temporary.path()).map_err(|error| {
        RunnerError::DatabaseCache {
            locator: redacted_locator.clone(),
            message: format!("unable to checksum downloaded database: {error}"),
        }
    })?;
    let stem = version.as_deref().map_or_else(
        || format!("content-{sha256}"),
        |version| database_cache_key(&redacted_locator, version, metadata.size),
    );
    if let Some(resolved) = cached_database(
        cache_dir,
        &stem,
        &redacted_locator,
        version.as_deref(),
        metadata.size,
        metrics,
    )? {
        return Ok(resolved);
    }
    let path = cache_dir.join(format!("candidate-{stem}.jam"));
    temporary
        .persist_noclobber(&path)
        .map_err(|error| RunnerError::DatabaseCache {
            locator: redacted_locator.clone(),
            message: format!("unable to persist cached database: {}", error.error),
        })?;
    let manifest = CandidateDatabaseCacheManifest {
        schema_version: "1.0.0".to_string(),
        redacted_locator: redacted_locator.clone(),
        version,
        size: metadata.size,
        sha256: sha256.clone(),
    };
    write_cache_manifest(cache_dir, &stem, &manifest)?;
    Ok(ResolvedCandidateDatabase {
        path,
        input: InputResource {
            role: "database".to_string(),
            redacted_locator,
            sha256: Some(sha256),
        },
        metrics,
    })
}

fn cached_database(
    cache_dir: &Path,
    stem: &str,
    redacted_locator: &str,
    version: Option<&str>,
    size: u64,
    metrics: ResourceMetrics,
) -> Result<Option<ResolvedCandidateDatabase>, RunnerError> {
    let path = cache_dir.join(format!("candidate-{stem}.jam"));
    let manifest_path = cache_dir.join(format!("candidate-{stem}.jam.cache.json"));
    match (path.exists(), manifest_path.exists()) {
        (false, false) => return Ok(None),
        (true, true) => {}
        _ => {
            return Err(RunnerError::DatabaseCache {
                locator: redacted_locator.to_string(),
                message: format!("incomplete cache entry for {}", path.display()),
            });
        }
    }
    let bytes = fs::read(&manifest_path).map_err(|error| RunnerError::DatabaseCache {
        locator: redacted_locator.to_string(),
        message: format!("unable to read {}: {error}", manifest_path.display()),
    })?;
    let manifest: CandidateDatabaseCacheManifest =
        serde_json::from_slice(&bytes).map_err(|error| RunnerError::DatabaseCache {
            locator: redacted_locator.to_string(),
            message: format!(
                "invalid cache manifest {}: {error}",
                manifest_path.display()
            ),
        })?;
    if manifest.schema_version != "1.0.0"
        || manifest.redacted_locator != redacted_locator
        || manifest.version.as_deref() != version
        || manifest.size != size
    {
        return Err(RunnerError::DatabaseCache {
            locator: redacted_locator.to_string(),
            message: format!("cache identity mismatch in {}", manifest_path.display()),
        });
    }
    let file_size = fs::metadata(&path)
        .map_err(|error| RunnerError::DatabaseCache {
            locator: redacted_locator.to_string(),
            message: format!("unable to inspect {}: {error}", path.display()),
        })?
        .len();
    let sha256 =
        crate::provenance::sha256_file(&path).map_err(|error| RunnerError::DatabaseCache {
            locator: redacted_locator.to_string(),
            message: format!("unable to checksum cached database: {error}"),
        })?;
    if file_size != size || sha256 != manifest.sha256 {
        return Err(RunnerError::DatabaseCache {
            locator: redacted_locator.to_string(),
            message: format!(
                "cached database failed size or checksum validation: {}",
                path.display()
            ),
        });
    }
    Ok(Some(ResolvedCandidateDatabase {
        path,
        input: InputResource {
            role: "database".to_string(),
            redacted_locator: redacted_locator.to_string(),
            sha256: Some(sha256),
        },
        metrics,
    }))
}

fn write_cache_manifest(
    cache_dir: &Path,
    stem: &str,
    manifest: &CandidateDatabaseCacheManifest,
) -> Result<(), RunnerError> {
    let path = cache_dir.join(format!("candidate-{stem}.jam.cache.json"));
    let mut temporary =
        tempfile::NamedTempFile::new_in(cache_dir).map_err(|error| RunnerError::DatabaseCache {
            locator: manifest.redacted_locator.clone(),
            message: format!("unable to create temporary cache manifest: {error}"),
        })?;
    serde_json::to_writer_pretty(temporary.as_file_mut(), manifest).map_err(|error| {
        RunnerError::DatabaseCache {
            locator: manifest.redacted_locator.clone(),
            message: format!("unable to encode cache manifest: {error}"),
        }
    })?;
    temporary
        .as_file_mut()
        .write_all(b"\n")
        .map_err(|error| RunnerError::DatabaseCache {
            locator: manifest.redacted_locator.clone(),
            message: format!("unable to finish cache manifest: {error}"),
        })?;
    temporary
        .persist_noclobber(path)
        .map_err(|error| RunnerError::DatabaseCache {
            locator: manifest.redacted_locator.clone(),
            message: format!("unable to persist cache manifest: {}", error.error),
        })?;
    Ok(())
}

fn database_version(metadata: &crate::resource::ResourceMetadata) -> Option<String> {
    metadata
        .etag
        .as_ref()
        .map(|value| format!("etag:{value}"))
        .or_else(|| {
            metadata
                .last_modified
                .as_ref()
                .map(|value| format!("last-modified:{value}"))
        })
}

fn database_cache_key(locator: &str, version: &str, size: u64) -> String {
    let mut digest = Sha256::new();
    digest.update(locator.as_bytes());
    digest.update([0]);
    digest.update(version.as_bytes());
    digest.update([0]);
    digest.update(size.to_le_bytes());
    format!("{:x}", digest.finalize())
}

fn add_resource_metrics(left: ResourceMetrics, right: ResourceMetrics) -> ResourceMetrics {
    left.saturating_add(right)
}

#[derive(Clone, Copy, Debug, Default)]
struct ArchiveMetricDelta {
    seed_buckets_read: u64,
    sequence_blocks_read: u64,
    bytes_read: u64,
}

fn subtract_archive_metrics(after: ArchiveMetrics, before: ArchiveMetrics) -> ArchiveMetricDelta {
    ArchiveMetricDelta {
        seed_buckets_read: after
            .resource
            .seed_buckets_read
            .saturating_sub(before.resource.seed_buckets_read),
        sequence_blocks_read: after
            .resource
            .sequence_blocks_read
            .saturating_sub(before.resource.sequence_blocks_read),
        bytes_read: after
            .resource
            .returned_bytes
            .saturating_sub(before.resource.returned_bytes),
    }
}

fn archive_error(error: ArchiveError) -> RunnerError {
    match error {
        ArchiveError::Resource(error) => RunnerError::Resource(error),
        ArchiveError::InvalidSequenceRange { start, end } => {
            RunnerError::Jma(JmaError::InvalidSequenceRange { start, end })
        }
        ArchiveError::UnknownContig(contig) => RunnerError::Jma(JmaError::UnknownContig(contig)),
        ArchiveError::ChecksumMismatch(section) => {
            RunnerError::Jma(JmaError::ChecksumMismatch(section))
        }
        ArchiveError::CorruptMetadata(message) | ArchiveError::Backend(message) => {
            RunnerError::Jma(JmaError::CorruptSection(message))
        }
        ArchiveError::UnknownSeedScheme(scheme) => RunnerError::Jma(JmaError::CorruptSection(
            format!("archive seed scheme {scheme} is unavailable"),
        )),
    }
}

/// Query an advertised exact Gear-fragment scheme without constructing a
/// target-wide index. Every returned occurrence is range-read and byte-
/// verified before an exact run can become a direct `=` alignment.
fn exact_gear_fragment_fast_path<A: TraceArchive>(
    archive: &A,
    query: &[u8],
    descriptor: SeedSchemeDescriptor,
    max_occurrences: u32,
) -> Result<Vec<BaseAlignment>, RunnerError> {
    let fragments = fragment_sequence(
        query,
        0,
        GearConfig::fine(GearTableKind::SingleBase, 0),
        FragmentationMode::StrandSymmetric,
    )
    .map_err(|error| RunnerError::InvalidConfig(format!("Gear fragmentation failed: {error}")))?;
    let metadata = archive.metadata().map_err(archive_error)?;
    let mut pending = Vec::new();
    let mut fragment_pairs = Vec::new();
    for fragment in fragments {
        let Some(bytes) = fragment_bytes(query, fragment) else {
            continue;
        };
        let lookup = archive
            .lookup_seeds(
                SeedSchemeId(descriptor.scheme_id),
                &[SeedKey {
                    digest: fragment.digest,
                    verification: bytes,
                }],
            )
            .map_err(archive_error)?;
        for seed_match in lookup.matches {
            for occurrence in seed_match
                .occurrences
                .into_iter()
                .take(usize::try_from(max_occurrences).unwrap_or(usize::MAX))
            {
                let Some(contig) = metadata
                    .contigs
                    .iter()
                    .find(|contig| contig.id == occurrence.contig_id)
                else {
                    continue;
                };
                let end = occurrence
                    .position
                    .checked_add(u64::from(occurrence.span))
                    .ok_or_else(|| {
                        RunnerError::InvalidConfig("Gear occurrence overflow".to_string())
                    })?;
                if end > contig.length || occurrence.span == 0 {
                    continue;
                }
                let request =
                    SequenceRequest::new(occurrence.contig_id, occurrence.position, end, false)
                        .map_err(archive_error)?;
                pending.push(request);
                fragment_pairs.push((fragment, occurrence, contig.length, end));
            }
        }
    }
    if pending.is_empty() {
        return Ok(Vec::new());
    }
    let slices = archive.read_sequences(&pending).map_err(archive_error)?;
    let mut verified = Vec::new();
    for ((query_fragment, occurrence, contig_length, occurrence_end), slice) in
        fragment_pairs.into_iter().zip(slices)
    {
        let target_fragment = ExactFragment {
            contig_id: occurrence.contig_id,
            start: occurrence.position,
            length: occurrence.span.into(),
            orientation: if occurrence.reverse {
                FragmentOrientation::Reverse
            } else {
                FragmentOrientation::Forward
            },
            digest: query_fragment.digest,
        };
        if !verify_exact_fragment(
            query,
            query_fragment,
            &slice.bases,
            ExactFragment {
                start: 0,
                ..target_fragment
            },
        ) {
            continue;
        }
        let relative_orientation =
            relative_fragment_orientation(query_fragment.orientation, target_fragment.orientation);
        let target_axis_start = if relative_orientation == FragmentOrientation::Reverse {
            contig_length.checked_sub(occurrence_end).ok_or_else(|| {
                RunnerError::InvalidConfig("reverse Gear occurrence exceeds contig".to_string())
            })?
        } else {
            occurrence.position
        };
        verified.push(VerifiedFragment {
            query_start: query_fragment.start,
            target_axis_start,
            length: u64::from(query_fragment.length),
            contig_id: occurrence.contig_id,
            orientation: relative_orientation,
        });
    }
    let runs = merge_exact_runs(&verified)
        .map_err(|error| RunnerError::InvalidConfig(format!("Gear run merge failed: {error}")))?;
    let mut alignments = Vec::new();
    for run in runs.into_iter().filter(|run| run.direct_cigar().is_some()) {
        let Some(contig) = metadata
            .contigs
            .iter()
            .find(|contig| contig.id == run.contig_id)
        else {
            continue;
        };
        let (target_start, target_end) = if run.orientation == FragmentOrientation::Reverse {
            (
                contig.length.checked_sub(run.target_end).ok_or_else(|| {
                    RunnerError::InvalidConfig("reverse Gear run exceeds contig".to_string())
                })?,
                contig.length.checked_sub(run.target_start).ok_or_else(|| {
                    RunnerError::InvalidConfig("reverse Gear run exceeds contig".to_string())
                })?,
            )
        } else {
            (run.target_start, run.target_end)
        };
        let query_interval = BaseInterval {
            start: run.query_start,
            end: run.query_end,
        };
        let exact_score = i64::try_from(run.verified_bases).map_err(|_| {
            RunnerError::InvalidConfig("exact Gear run score exceeds i64".to_string())
        })?;
        alignments.push(BaseAlignment {
            alignment_id: String::new(),
            plasmid_id: String::new(),
            metagenome_id: String::new(),
            contig_id: contig.name.clone(),
            strand: if run.orientation == FragmentOrientation::Reverse {
                Strand::Reverse
            } else {
                Strand::Forward
            },
            query_segments: vec![query_interval],
            target_interval: BaseInterval {
                start: target_start,
                end: target_end,
            },
            query_length: run.verified_bases,
            target_length: run.verified_bases,
            origin_crossing: false,
            score: exact_score,
            matches: run.verified_bases,
            substitutions: 0,
            insertions: 0,
            deletions: 0,
            cigar: run.direct_cigar().unwrap_or_default(),
            edit_script: equal_edit_script(run.verified_bases),
            chain_score: exact_score,
            identity: 1.0,
            seed_evidence: SeedEvidence::default(),
            primary_supported_bases: 0,
            secondary_supported_bases: 0,
            newly_supported_bases: 0,
            role: AlignmentRole::AlternativeMapping,
            primary: false,
        });
    }
    Ok(alignments)
}

fn equal_edit_script(mut length: u64) -> Vec<crate::trace::model::EditRun> {
    let mut runs = Vec::new();
    while length != 0 {
        let chunk = length.min(u64::from(u32::MAX));
        runs.push(crate::trace::model::EditRun {
            operation: crate::trace::model::EditOperation::Equal,
            length: u32::try_from(chunk).expect("chunk is bounded by u32::MAX"),
        });
        length -= chunk;
    }
    runs
}

const fn relative_fragment_orientation(
    query: FragmentOrientation,
    target: FragmentOrientation,
) -> FragmentOrientation {
    if matches!(query, FragmentOrientation::Reverse)
        ^ matches!(target, FragmentOrientation::Reverse)
    {
        FragmentOrientation::Reverse
    } else {
        FragmentOrientation::Forward
    }
}

/// Executes candidate retrieval and bounded positional processing.
pub struct TraceRunner {
    config: TraceRunnerConfig,
}

impl TraceRunner {
    pub fn new(config: TraceRunnerConfig) -> Result<Self, RunnerError> {
        config.validate()?;
        Ok(Self { config })
    }

    #[must_use]
    pub fn config(&self) -> &TraceRunnerConfig {
        &self.config
    }

    pub fn run(&self, query: &TraceQuery) -> Result<TraceRunOutput, RunnerError> {
        query.validate()?;
        let started = Instant::now();
        let database = resolve_candidate_database(&query.database, &self.config.resources)?;
        let searcher = CandidateSearcher::open(&database.path)?;
        let positional_candidate_limit =
            usize::try_from(self.config.sensitivity.max_candidates).unwrap_or(usize::MAX);
        let mut candidate_config = self.config.candidates.clone();
        candidate_config.top_candidates = candidate_config
            .top_candidates
            .min(positional_candidate_limit);
        let mut search = searcher.search_sequence(&query.plasmid_sequence, &candidate_config)?;
        // CandidateSearchConfig bounds the sketch result, while the
        // sensitivity profile bounds positional work.  Apply both limits
        // before constructing the worker input so a pathological search
        // result cannot create an unbounded queue or result vector.
        search.candidates.truncate(positional_candidate_limit);
        self.run_selected(
            query,
            search,
            database.input,
            database.metrics,
            None,
            started,
        )
    }

    pub fn run_routed(
        &self,
        query: &TraceQuery,
        routed: Vec<RoutedCandidate>,
        router_input: InputResource,
        router_metrics: ResourceMetrics,
    ) -> Result<TraceRunOutput, RunnerError> {
        query.validate()?;
        let started = Instant::now();
        let positional_candidate_limit =
            usize::try_from(self.config.sensitivity.max_candidates).unwrap_or(usize::MAX);
        let mut routed = routed;
        routed.truncate(positional_candidate_limit);
        let candidates = routed
            .iter()
            .enumerate()
            .map(|(index, candidate)| RankedCandidate {
                candidate: CandidateResult {
                    metagenome_id: candidate.metagenome_id.clone(),
                    shared_hashes: u64::from(
                        candidate
                            .rare_shared_witnesses
                            .saturating_add(candidate.common_shared_witnesses),
                    ),
                    plasmid_hashes: candidate
                        .window_evidence
                        .total_shared_witness_count
                        .max(1)
                        .into(),
                    metagenome_hashes: 0,
                    plasmid_containment: candidate.estimated_query_containment,
                    metagenome_containment: 0.0,
                    rank: u32::try_from(index + 1).unwrap_or(u32::MAX),
                    score_mode: CandidateScoreMode::Witness.as_str().to_string(),
                    bias_weighted_plasmid_containment: None,
                    uniform_hash_e_value: None,
                },
                score_mode: CandidateScoreMode::Witness,
                bias_weighted_plasmid_containment: None,
                uniform_hash_e_value: None,
            })
            .collect::<Vec<_>>();
        let search = CandidateSearchResult {
            query_hashes: candidates
                .iter()
                .map(|candidate| candidate.candidate.plasmid_hashes)
                .max()
                .unwrap_or(0),
            hashes_found: candidates
                .iter()
                .map(|candidate| candidate.candidate.shared_hashes)
                .max()
                .unwrap_or(0),
            candidates,
            score_mode: CandidateScoreMode::Witness,
        };
        let routed = routed
            .into_iter()
            .map(|candidate| (candidate.metagenome_id.clone(), candidate))
            .collect::<BTreeMap<_, _>>();
        self.run_selected(
            query,
            search,
            router_input,
            router_metrics,
            Some(&routed),
            started,
        )
    }

    fn run_selected(
        &self,
        query: &TraceQuery,
        search: CandidateSearchResult,
        database_input: InputResource,
        database_metrics: ResourceMetrics,
        routed: Option<&BTreeMap<String, RoutedCandidate>>,
        started: Instant,
    ) -> Result<TraceRunOutput, RunnerError> {
        let mut counters = TracePerformanceCounters {
            candidate_queries: 1,
            candidates_considered: search.candidates.len() as u64,
            ..TracePerformanceCounters::default()
        };

        let run_config = &self.config;
        let plasmid_id = &query.plasmid_id;
        let plasmid_sequence = &query.plasmid_sequence;
        let worker_count = self
            .config
            .threads
            .min(self.config.io_concurrency)
            .min(self.config.sensitivity.max_concurrent_candidates as usize)
            .max(1);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(worker_count)
            .build()
            .map_err(|error| RunnerError::InvalidConfig(error.to_string()))?;
        let memory = TraceMemorySemaphore::new(self.config.memory_budget_bytes);
        let memory_estimate = candidate_memory_estimate(&self.config, plasmid_sequence.len());
        let mut work = pool.install(|| {
            search
                .candidates
                .par_iter()
                .map(|candidate| {
                    let reservation = memory.acquire(memory_estimate).map_err(|error| {
                        RunnerError::InvalidConfig(format!(
                            "candidate memory admission failed: {error}"
                        ))
                    })?;
                    let candidate_started = Instant::now();
                    let entry = query.catalog.get(candidate.metagenome_id());
                    let router_candidate =
                        routed.and_then(|candidates| candidates.get(candidate.metagenome_id()));
                    let mut work = process_candidate(
                        plasmid_id,
                        plasmid_sequence,
                        candidate,
                        entry,
                        router_candidate,
                        run_config,
                    );
                    work.counters.failures = work.result.failures.len() as u64;
                    work.counters.elapsed_millis = candidate_started
                        .elapsed()
                        .as_millis()
                        .try_into()
                        .unwrap_or(u64::MAX);
                    work.result.performance_counters = work.counters.candidate_snapshot();
                    // Account a complete candidate serialization before the
                    // permit is released. The public in-memory API writes the
                    // same record later with its final run identifier.
                    let serialized = serde_json::to_vec(&work.result).map_err(|error| {
                        RunnerError::InvalidConfig(format!(
                            "candidate serialization failed: {error}"
                        ))
                    })?;
                    drop(serialized);
                    drop(reservation);
                    Ok(work)
                })
                .collect::<Result<Vec<_>, RunnerError>>()
        })?;
        work.sort_by(|left, right| {
            left.result
                .metagenome_id
                .cmp(&right.result.metagenome_id)
                .then_with(|| {
                    left.result
                        .candidate
                        .as_ref()
                        .map(|candidate| candidate.rank)
                        .cmp(
                            &right
                                .result
                                .candidate
                                .as_ref()
                                .map(|candidate| candidate.rank),
                        )
                })
        });
        let mut metagenomes = Vec::with_capacity(work.len());
        for item in work {
            counters.add_assign(item.counters);
            metagenomes.push(item.result);
        }
        counters.failures = metagenomes
            .iter()
            .map(|result| result.failures.len() as u64)
            .sum();
        counters.elapsed_millis = started.elapsed().as_millis().try_into().unwrap_or(u64::MAX);
        Ok(TraceRunOutput {
            search,
            metagenomes,
            counters,
            database_input,
            database_metrics,
        })
    }
}

struct CandidateWork {
    result: TraceMetagenomeResult,
    counters: TracePerformanceCounters,
}

#[derive(Debug)]
struct TraceProcessingFailure {
    error: Box<RunnerError>,
    metrics: Box<ResourceMetrics>,
}

impl TraceProcessingFailure {
    fn new(error: RunnerError, metrics: ResourceMetrics) -> Self {
        Self {
            error: Box::new(error),
            metrics: Box::new(metrics),
        }
    }
}

fn process_candidate(
    plasmid_id: &str,
    plasmid: &[u8],
    candidate: &RankedCandidate,
    entry: Option<&CatalogEntry>,
    router_candidate: Option<&RoutedCandidate>,
    config: &TraceRunnerConfig,
) -> CandidateWork {
    let mut counters = TracePerformanceCounters {
        candidates_processed: 1,
        ..TracePerformanceCounters::default()
    };
    let Some(entry) = entry else {
        let mut result = failed_result(
            plasmid_id,
            candidate,
            config,
            "catalog",
            "missing_catalog_entry",
            "candidate is absent from the resource catalog",
            false,
            None,
        );
        result.router_candidate = router_candidate.cloned();
        return CandidateWork { result, counters };
    };

    let resource = entry.resource();
    let mut result = match process_jma(
        plasmid_id,
        plasmid,
        candidate,
        resource,
        &entry.sha256,
        router_candidate,
        config,
        &mut counters,
    ) {
        Ok(processed) => processed.result,
        Err(error) => failed_result_with_metrics(
            plasmid_id,
            candidate,
            config,
            "jma",
            error.error.code(),
            &error.error.to_string(),
            error.error.retryable(),
            Some(resource),
            *error.metrics,
        ),
    };
    result.router_candidate = router_candidate.cloned();
    CandidateWork { result, counters }
}

struct ProcessedJma {
    result: TraceMetagenomeResult,
}

fn positional_round_configs(sensitivity: &SensitivityConfig) -> Vec<Vec<SeedSensitivity>> {
    let mut rounds = vec![vec![sensitivity.primary]];
    let mut primary_for_mixed = sensitivity.primary;
    if let Some(dense_primary) = sensitivity.gap_rescue.dense_primary
        && (dense_primary.k, dense_primary.scale)
            != (sensitivity.primary.k, sensitivity.primary.scale)
    {
        primary_for_mixed = dense_primary;
        rounds.push(vec![dense_primary]);
    }
    if let Some(rescue) = sensitivity.rescue {
        rounds.push(vec![primary_for_mixed, rescue]);
    }
    rounds.truncate(usize::from(sensitivity.gap_rescue.max_rounds));
    rounds
}

#[allow(clippy::too_many_arguments)]
fn process_jma(
    plasmid_id: &str,
    plasmid: &[u8],
    candidate: &RankedCandidate,
    locator: &str,
    expected_sha256: &str,
    router_candidate: Option<&RoutedCandidate>,
    config: &TraceRunnerConfig,
    counters: &mut TracePerformanceCounters,
) -> Result<ProcessedJma, TraceProcessingFailure> {
    let archive_options = ResourceOpenOptions {
        allow_full_download_fallback: false,
        expected_sha256: Some(expected_sha256.to_string()),
        ..config.resources.clone()
    };
    let resource =
        Arc::new(open_resource(locator, archive_options).map_err(|error| {
            TraceProcessingFailure::new(error.into(), ResourceMetrics::default())
        })?);
    let archive = NativeJmaArchive::from_resource(Arc::clone(&resource))
        .map_err(|error| TraceProcessingFailure::new(archive_error(error), resource.metrics()))?;
    if archive.reader().header().algorithm_id.as_deref() != Some(crate::trace::TRACE_ALGORITHM_ID)
        || archive.reader().header().algorithm_version
            != Some(crate::trace::TRACE_ALGORITHM_VERSION)
    {
        return Err(TraceProcessingFailure::new(
            RunnerError::from(JmaError::CorruptSection(format!(
                "JMA algorithm metadata is missing or incompatible; expected {} version {}",
                crate::trace::TRACE_ALGORITHM_ID,
                crate::trace::TRACE_ALGORITHM_VERSION
            ))),
            archive.metrics().resource,
        ));
    }
    process_indexed_archive(
        plasmid_id,
        plasmid,
        candidate,
        &archive,
        router_candidate,
        config,
        counters,
    )
}

#[allow(clippy::too_many_arguments)]
fn process_indexed_archive<A: TraceArchive>(
    plasmid_id: &str,
    plasmid: &[u8],
    candidate: &RankedCandidate,
    archive: &A,
    router_candidate: Option<&RoutedCandidate>,
    config: &TraceRunnerConfig,
    counters: &mut TracePerformanceCounters,
) -> Result<ProcessedJma, TraceProcessingFailure> {
    let mut alignments = Vec::new();
    let mut workspace = AlignmentWorkspace::new();
    let mut seen_contigs = HashSet::new();
    let mut seen_windows = HashSet::new();
    let mut seen_seed_keys = HashSet::new();
    let mut retry_metadata = Vec::new();
    let mut stage_metrics = Vec::new();
    let descriptors = archive.available_seed_schemes().map_err(|error| {
        TraceProcessingFailure::new(archive_error(error), archive.metrics().resource)
    })?;
    if let Some(descriptor) = cascade::advertised_exact_gear_scheme(&descriptors) {
        let started = Instant::now();
        let exact = exact_gear_fragment_fast_path(
            archive,
            plasmid,
            descriptor,
            config.sensitivity.primary.max_occurrences,
        )
        .map_err(|error| TraceProcessingFailure::new(error, archive.metrics().resource))?;
        let exact_count = exact.len() as u64;
        for mut alignment in exact {
            alignment.plasmid_id = plasmid_id.to_string();
            alignment.metagenome_id = candidate.metagenome_id().to_string();
            alignment.seed_evidence.primary_anchor_count = 1;
            alignments.push(alignment);
        }
        stage_metrics.push(TraceStageMetrics {
            stage: 0,
            name: cascade::CascadeStageKind::ExactGearFragments
                .name()
                .to_string(),
            anchors_retained: exact_count,
            chains_retained: exact_count,
            // Direct '=' runs intentionally perform no DP alignment.
            alignment_attempts: 0,
            new_query_bases_supported: alignments.iter().map(|alignment| alignment.matches).sum(),
            wall_micros: started.elapsed().as_micros().try_into().unwrap_or(u64::MAX),
            ..TraceStageMetrics::default()
        });
    }
    let chain_model = chaining_coordinate_model(config.topology_requested);
    let round_configs = positional_round_configs(&config.sensitivity);

    let query_length = plasmid.len() as u64;
    let mut rescue_rounds = Vec::new();
    let mut target_gaps = vec![BaseInterval {
        start: 0,
        end: query_length,
    }];
    let mut supported_before = 0_u64;
    for (round_index, seed_configs) in round_configs.into_iter().enumerate() {
        if round_index > 0 {
            target_gaps.retain(|gap| gap.len() >= config.sensitivity.gap_rescue.min_gap_bases);
            if target_gaps.is_empty() {
                break;
            }
        }
        let round_target_gaps = target_gaps.clone();

        let started = Instant::now();
        let before_metrics = archive.metrics();
        let round = chains_for_configs(
            archive,
            plasmid,
            &seed_configs,
            &config.sensitivity,
            chain_model,
            (round_index > 0).then_some(target_gaps.as_slice()),
            router_candidate,
            &mut seen_seed_keys,
        )
        .map_err(|error| TraceProcessingFailure::new(error, archive.metrics().resource))?;
        let LevelChainResult {
            chains,
            seed_keys_tested,
            anchors_created,
            chains_accepted,
            occurrences_decoded,
        } = round;
        let protected_alignments = alignments.len();
        let alignment_windows_attempted = align_indexed_chains(
            archive,
            plasmid_id,
            plasmid,
            candidate,
            config,
            chains,
            if round_index == 0 {
                u64::MAX
            } else {
                u64::from(
                    config
                        .sensitivity
                        .gap_rescue
                        .max_alignment_windows_per_round,
                )
            },
            if round_index == 0 {
                u64::MAX
            } else {
                u64::from(config.sensitivity.gap_rescue.max_sequence_blocks_per_round)
            },
            &mut workspace,
            &mut seen_contigs,
            &mut seen_windows,
            &mut alignments,
            &mut retry_metadata,
            counters,
        )
        .map_err(|error| TraceProcessingFailure::new(error, archive.metrics().resource))?;
        retain_new_alignments(
            &mut alignments,
            protected_alignments,
            config.max_alignments_per_candidate,
        );

        let topology = assess_topology(
            query_length,
            config.topology_requested,
            config.topology_margin_bases,
            &alignments,
        )
        .map_err(|error| TraceProcessingFailure::new(error.into(), archive.metrics().resource))?;
        let selected = selected_mosaic(&topology);
        let supported_after = selected.base_covered_bases;
        target_gaps = selected
            .unsupported_gaps
            .iter()
            .map(|gap| gap.interval)
            .collect();
        let after_metrics = archive.metrics();
        let metric_delta = subtract_archive_metrics(after_metrics, before_metrics);
        let new_query_bases_supported = supported_after.saturating_sub(supported_before);
        let reported_seed = seed_configs
            .last()
            .copied()
            .unwrap_or(config.sensitivity.primary);
        rescue_rounds.push(RescueRoundMetrics {
            round: u8::try_from(round_index + 1).unwrap_or(u8::MAX),
            seed_k: reported_seed.k,
            seed_scale: reported_seed.scale,
            target_gaps: round_target_gaps,
            seed_buckets_requested: metric_delta.seed_buckets_read,
            seed_keys_tested,
            anchors_created,
            chains_accepted,
            sequence_blocks_fetched: metric_delta.sequence_blocks_read,
            alignment_windows_attempted,
            new_query_bases_supported,
            elapsed_millis: started.elapsed().as_millis().try_into().unwrap_or(u64::MAX),
        });
        stage_metrics.push(TraceStageMetrics {
            stage: u8::try_from(round_index + 1).unwrap_or(u8::MAX),
            name: if round_index == 0 {
                "sparse_k31".to_string()
            } else if seed_configs.len() > 1 {
                "mixed_gap_k31_k21".to_string()
            } else if reported_seed.k == crate::trace::seeds::PRIMARY_K {
                "dense_gap_k31".to_string()
            } else {
                "rescue_k21".to_string()
            },
            seed_pages_read: metric_delta.seed_buckets_read,
            keys_tested: seed_keys_tested,
            occurrences_decoded,
            anchors_retained: anchors_created,
            chains_retained: chains_accepted,
            sequence_blocks_read: metric_delta.sequence_blocks_read,
            alignment_attempts: alignment_windows_attempted,
            new_query_bases_supported,
            wall_micros: started.elapsed().as_micros().try_into().unwrap_or(u64::MAX),
            cpu_micros: 0,
            bytes_read: metric_delta.bytes_read,
        });
        if (round_index > 0 && new_query_bases_supported == 0) || supported_after == query_length {
            break;
        }
        supported_before = supported_after;
    }

    let finalized = finalize_evidence(query_length, config, &mut alignments)
        .map_err(|error| TraceProcessingFailure::new(error, archive.metrics().resource))?;
    let warnings = evidence_warnings(&finalized);
    let result = TraceMetagenomeResult {
        schema_version: crate::trace::TRACE_JSON_SCHEMA_VERSION.to_string(),
        run_id: String::new(),
        plasmid_id: plasmid_id.to_string(),
        metagenome_id: candidate.metagenome_id().to_string(),
        query_kind: config.query_kind,
        topology_requested: config.topology_requested,
        coordinate_model: finalized.topology.coordinate_model,
        topology_evidence: finalized.topology.topology_evidence,
        algorithms: crate::trace::config::algorithm_identifiers(),
        algorithm: crate::trace::config::TraceAlgorithmMetadata::for_sensitivity(
            config.sensitivity.clone(),
        ),
        status: TraceStatus::Complete,
        candidate: Some(candidate.candidate.clone()),
        router_candidate: None,
        alignments,
        primary_fragment_mosaic: Some(finalized.primary_mosaic),
        topology: Some(finalized.topology),
        rescue_rounds,
        stages: stage_metrics,
        alignment_retries: retry_metadata,
        performance_counters: CandidatePerformanceCounters::default(),
        coverage: Some(finalized.coverage),
        warnings,
        failures: Vec::new(),
        archive_metrics: Some(archive.metrics()),
        resource_metrics: archive.metrics().resource,
    };
    Ok(ProcessedJma { result })
}

#[allow(dead_code)]
pub(crate) fn run_index_archive<A: TraceArchive>(
    query_id: &str,
    query: &[u8],
    candidate: CandidateResult,
    archive: &A,
    config: &TraceRunnerConfig,
) -> Result<TraceMetagenomeResult, RunnerError> {
    let ranked = RankedCandidate {
        candidate,
        score_mode: CandidateScoreMode::Witness,
        bias_weighted_plasmid_containment: None,
        uniform_hash_e_value: None,
    };
    let mut counters = TracePerformanceCounters {
        candidates_processed: 1,
        ..TracePerformanceCounters::default()
    };
    let started = Instant::now();
    let mut result = process_indexed_archive(
        query_id,
        query,
        &ranked,
        archive,
        None,
        config,
        &mut counters,
    )
    .map_err(|failure| *failure.error)?
    .result;
    counters.failures = result.failures.len() as u64;
    counters.elapsed_millis = started.elapsed().as_millis().try_into().unwrap_or(u64::MAX);
    result.performance_counters = counters.candidate_snapshot();
    Ok(result)
}

pub(crate) fn run_exact_index(
    query_id: &str,
    query_length: usize,
    candidate: CandidateResult,
    mut alignments: Vec<BaseAlignment>,
    archive_metrics: ArchiveMetrics,
    contigs_considered: u64,
    config: &TraceRunnerConfig,
) -> Result<TraceMetagenomeResult, RunnerError> {
    let alignments_succeeded = u64::try_from(alignments.len()).unwrap_or(u64::MAX);
    let finalized = finalize_evidence(query_length as u64, config, &mut alignments)?;
    let warnings = evidence_warnings(&finalized);
    let supported_bases = finalized.coverage.supported_bases;
    Ok(TraceMetagenomeResult {
        schema_version: crate::trace::TRACE_JSON_SCHEMA_VERSION.to_string(),
        run_id: String::new(),
        plasmid_id: query_id.to_string(),
        metagenome_id: candidate.metagenome_id.clone(),
        query_kind: config.query_kind,
        topology_requested: config.topology_requested,
        coordinate_model: finalized.topology.coordinate_model,
        topology_evidence: finalized.topology.topology_evidence,
        algorithms: crate::trace::config::algorithm_identifiers(),
        algorithm: crate::trace::config::TraceAlgorithmMetadata::for_sensitivity(
            config.sensitivity.clone(),
        ),
        status: TraceStatus::Complete,
        candidate: Some(candidate),
        router_candidate: None,
        alignments,
        primary_fragment_mosaic: Some(finalized.primary_mosaic),
        topology: Some(finalized.topology),
        rescue_rounds: Vec::new(),
        stages: vec![TraceStageMetrics {
            stage: 0,
            name: "exact_contig".to_string(),
            new_query_bases_supported: supported_bases,
            ..TraceStageMetrics::default()
        }],
        alignment_retries: Vec::new(),
        performance_counters: CandidatePerformanceCounters {
            candidates_processed: 1,
            contigs_considered,
            windows_retrieved: contigs_considered,
            alignments_succeeded,
            ..CandidatePerformanceCounters::default()
        },
        coverage: Some(finalized.coverage),
        warnings,
        failures: Vec::new(),
        archive_metrics: Some(archive_metrics),
        resource_metrics: archive_metrics.resource,
    })
}

pub(crate) fn merge_index_results(
    mut merged: TraceMetagenomeResult,
    mut next: TraceMetagenomeResult,
    query_length: usize,
    config: &TraceRunnerConfig,
) -> Result<TraceMetagenomeResult, RunnerError> {
    if merged.plasmid_id != next.plasmid_id || merged.metagenome_id != next.metagenome_id {
        return Err(RunnerError::InvalidConfig(
            "cannot merge trace results from different queries or metagenomes".to_string(),
        ));
    }
    merged.alignments.append(&mut next.alignments);
    retain_best_alignments(&mut merged.alignments, config.max_alignments_per_candidate);
    let finalized = finalize_evidence(query_length as u64, config, &mut merged.alignments)?;
    merged.warnings = evidence_warnings(&finalized);
    merged.coordinate_model = finalized.topology.coordinate_model;
    merged.topology_evidence = finalized.topology.topology_evidence;
    merged.primary_fragment_mosaic = Some(finalized.primary_mosaic);
    merged.topology = Some(finalized.topology);
    merged.coverage = Some(finalized.coverage);
    merged.rescue_rounds.append(&mut next.rescue_rounds);
    merged.stages.append(&mut next.stages);
    merged.alignment_retries.append(&mut next.alignment_retries);
    merged.failures.append(&mut next.failures);
    merged.performance_counters =
        merge_candidate_counters(merged.performance_counters, next.performance_counters);
    merged.archive_metrics = merge_archive_metrics(merged.archive_metrics, next.archive_metrics);
    merged.resource_metrics = merged
        .resource_metrics
        .saturating_add(next.resource_metrics);
    Ok(merged)
}

fn merge_candidate_counters(
    left: CandidatePerformanceCounters,
    right: CandidatePerformanceCounters,
) -> CandidatePerformanceCounters {
    CandidatePerformanceCounters {
        candidates_processed: left
            .candidates_processed
            .saturating_add(right.candidates_processed),
        contigs_considered: left
            .contigs_considered
            .saturating_add(right.contigs_considered),
        windows_retrieved: left
            .windows_retrieved
            .saturating_add(right.windows_retrieved),
        alignments_attempted: left
            .alignments_attempted
            .saturating_add(right.alignments_attempted),
        alignments_succeeded: left
            .alignments_succeeded
            .saturating_add(right.alignments_succeeded),
        failures: left.failures.saturating_add(right.failures),
        elapsed_millis: left.elapsed_millis.saturating_add(right.elapsed_millis),
    }
}

fn merge_archive_metrics(
    left: Option<ArchiveMetrics>,
    right: Option<ArchiveMetrics>,
) -> Option<ArchiveMetrics> {
    match (left, right) {
        (Some(left), Some(right)) => Some(ArchiveMetrics {
            resource: left.resource.saturating_add(right.resource),
            mapped_bytes: left.mapped_bytes.saturating_add(right.mapped_bytes),
            resident_bytes: left.resident_bytes.saturating_add(right.resident_bytes),
            metadata_bytes_read: left
                .metadata_bytes_read
                .saturating_add(right.metadata_bytes_read),
            seed_bytes_read: left.seed_bytes_read.saturating_add(right.seed_bytes_read),
            sequence_bytes_read: left
                .sequence_bytes_read
                .saturating_add(right.sequence_bytes_read),
            decoded_sequence_bases: left
                .decoded_sequence_bases
                .saturating_add(right.decoded_sequence_bases),
            coalesced_range_requests: left
                .coalesced_range_requests
                .saturating_add(right.coalesced_range_requests),
        }),
        (Some(metrics), None) | (None, Some(metrics)) => Some(metrics),
        (None, None) => None,
    }
}

#[allow(clippy::too_many_arguments)]
fn align_window(
    plasmid_id: &str,
    candidate: &RankedCandidate,
    workspace: &mut AlignmentWorkspace,
    plasmid: &[u8],
    target: &[u8],
    target_start: u64,
    contig_id: &str,
    config: &TraceRunnerConfig,
    chain_score: i64,
    seed_evidence: SeedEvidence,
    strand_hint: Option<Strand>,
    query_start: u64,
    query_span: u64,
    diagonal_offset: isize,
    corridor: Option<&AlignmentCorridor>,
    alignments: &mut Vec<BaseAlignment>,
    retry_metadata: &mut Vec<AlignmentRetryMetadata>,
    counters: &mut TracePerformanceCounters,
) -> Result<(), RunnerError> {
    let band = usize::try_from(config.sensitivity.alignment.band_width)
        .map_err(|_| RunnerError::InvalidConfig("alignment band exceeds platform".to_string()))?;
    let query_span_usize = usize::try_from(query_span)
        .map_err(|_| RunnerError::InvalidConfig("query span exceeds platform".to_string()))?;
    if query_span_usize == 0 {
        return Err(RunnerError::InvalidConfig(
            "query span must be greater than zero".to_string(),
        ));
    }
    let retry_band = corridor
        .map(|corridor| corridor.required_half_width())
        .unwrap_or(config.sensitivity.alignment.band_width)
        .max(config.sensitivity.alignment.band_width)
        .min(alignment_half_width_cap(config, query_span));
    let max_cells = query_span_usize.saturating_add(1).saturating_mul(
        usize::try_from(retry_band)
            .unwrap_or(band)
            .saturating_mul(2)
            .saturating_add(1),
    );
    let effective_diagonal_offset = corridor.map_or(diagonal_offset, |_| 0);
    let options = AlignmentOptions::new(config.sensitivity.alignment)
        .with_diagonal_offset(effective_diagonal_offset)
        .with_max_cells(max_cells.max(1));
    let both_strands = [Strand::Forward, Strand::Reverse];
    let hinted_strand = [strand_hint.unwrap_or(Strand::Forward)];
    let strands = if strand_hint.is_some() {
        &hinted_strand[..]
    } else {
        &both_strands[..]
    };
    for &strand in strands {
        counters.alignments_attempted = counters.alignments_attempted.saturating_add(1);
        let result = if let Some(corridor) = corridor {
            let result = workspace.align_circular_with_retries(
                plasmid,
                query_start,
                query_span,
                target,
                target_start,
                strand,
                corridor,
                options,
            );
            match result {
                Ok(result)
                    if result.retry_metadata.band_edge_touched
                        || result.retry_metadata.predicted_drift
                        || result.retry_metadata.chain_span_under_explained =>
                {
                    workspace
                        .align_circular_anchored_semiglobal(
                            plasmid,
                            query_start,
                            query_span,
                            target,
                            target_start,
                            strand,
                            corridor,
                            options,
                        )
                        .or(Ok(result))
                }
                other => other,
            }
        } else {
            workspace.align_circular(
                plasmid,
                query_start,
                query_span,
                target,
                target_start,
                strand,
                options,
            )
        };
        match result {
            Ok(result) => {
                counters.alignments_succeeded = counters.alignments_succeeded.saturating_add(1);
                retry_metadata.push(result.retry_metadata.clone());
                let mut alignment = result.into_base_alignment(
                    plasmid_id,
                    candidate.metagenome_id(),
                    contig_id,
                    chain_score,
                    true,
                );
                alignment.seed_evidence = seed_evidence.clone();
                if alignment.matches > 0 {
                    alignments.push(alignment);
                }
            }
            Err(AlignmentError::NoAlignment) => {}
            Err(error) => return Err(RunnerError::Alignment(error)),
        }
    }
    Ok(())
}

struct EvidencedChain {
    chain: WeightedAnchorChain,
    seed_evidence: SeedEvidence,
}

struct LevelChainResult {
    chains: Vec<EvidencedChain>,
    seed_keys_tested: u64,
    occurrences_decoded: u64,
    anchors_created: u64,
    chains_accepted: u64,
}

#[allow(clippy::too_many_arguments)]
fn chains_for_configs<A: TraceArchive>(
    archive: &A,
    plasmid: &[u8],
    seed_configs: &[SeedSensitivity],
    sensitivity: &SensitivityConfig,
    coordinate_model: CoordinateModel,
    intervals: Option<&[BaseInterval]>,
    router_candidate: Option<&RoutedCandidate>,
    seen_seed_keys: &mut HashSet<(u8, u64, u64, u64)>,
) -> Result<LevelChainResult, RunnerError> {
    let descriptors = archive.available_seed_schemes().map_err(archive_error)?;
    let mut groups = Vec::new();
    let mut seed_keys_tested = 0_u64;
    let mut occurrences_decoded = 0_u64;
    let mut remaining_bucket_budget = intervals.map_or(usize::MAX, |_| {
        usize::try_from(sensitivity.gap_rescue.max_seed_buckets_per_round).unwrap_or(usize::MAX)
    });

    for seed_config in seed_configs {
        let effective_seed_config = if seed_config.k == crate::trace::seeds::RESCUE_K {
            router_candidate.map_or(*seed_config, |candidate| SeedSensitivity {
                scale: u64::from(candidate.candidate_tier),
                ..*seed_config
            })
        } else {
            *seed_config
        };
        let scheme = select_seed_scheme(&descriptors, effective_seed_config)?;
        let mut level = match intervals {
            Some(intervals) => extract_seed_level_in_intervals(
                plasmid,
                effective_seed_config,
                intervals,
                sensitivity.gap_rescue.flank_bases,
            )?,
            None => extract_seed_level(plasmid, effective_seed_config)?,
        };
        if seed_config.k == crate::trace::seeds::RESCUE_K
            && let Some(router_candidate) = router_candidate
        {
            let handed = router_candidate
                .shared_witnesses
                .iter()
                .map(|witness| {
                    (
                        witness.key.jamhash,
                        witness.key.packed,
                        witness.query_position,
                        witness.query_reverse,
                    )
                })
                .collect::<BTreeSet<_>>();
            level.seeds.retain(|seed| {
                handed.contains(&(seed.hash, seed.canonical_kmer, seed.position, seed.reverse))
            });
        }
        let reuse_primary_for_mixed =
            seed_configs.len() > 1 && seed_config.k == crate::trace::seeds::PRIMARY_K;
        level.seeds.retain(|seed| {
            reuse_primary_for_mixed
                || seen_seed_keys.insert((level.k, seed.hash, seed.canonical_kmer, seed.position))
        });
        if intervals.is_some() {
            level.seeds.truncate(remaining_bucket_budget);
            remaining_bucket_budget = remaining_bucket_budget.saturating_sub(level.seeds.len());
        }
        seed_keys_tested = seed_keys_tested.saturating_add(level.seeds.len() as u64);
        let mut handed_positions =
            BTreeMap::<(u64, u64, u64, bool), Vec<crate::jma::SeedOccurrence>>::new();
        if seed_config.k == crate::trace::seeds::RESCUE_K
            && let Some(router_candidate) = router_candidate
        {
            for position in &router_candidate.positional_witnesses {
                if position.scheme_id != scheme.scheme_id {
                    return Err(RunnerError::InvalidConfig(
                        "router/JMA witness scheme mismatch".to_string(),
                    ));
                }
                handed_positions
                    .entry((
                        position.witness.key.jamhash,
                        position.witness.key.packed,
                        position.witness.query_position,
                        position.witness.query_reverse,
                    ))
                    .or_default()
                    .push(crate::jma::SeedOccurrence {
                        contig_id: position.contig_id,
                        position: position.position,
                        reverse: position.reverse,
                    });
            }
        }
        let mut by_key = BTreeMap::<SeedKey, Vec<QuerySeed>>::new();
        for seed in level.seeds {
            if let Some(occurrences) =
                handed_positions.get(&(seed.hash, seed.canonical_kmer, seed.position, seed.reverse))
            {
                groups.push(SeedOccurrenceGroup {
                    seed,
                    k: level.k,
                    occurrences: occurrences.clone(),
                });
                continue;
            }
            by_key
                .entry(seed_key(seed, level.k))
                .or_default()
                .push(seed);
        }
        if by_key.is_empty() {
            continue;
        }
        let keys = by_key.keys().cloned().collect::<Vec<_>>();
        let lookup = archive
            .lookup_seeds_bounded(
                SeedSchemeId(scheme.scheme_id),
                &keys,
                Some(seed_config.max_occurrences),
            )
            .map_err(archive_error)?;
        let mut matches = BTreeMap::new();
        for seed_match in lookup.matches {
            if let Some(key) = keys.get(seed_match.key_index) {
                occurrences_decoded =
                    occurrences_decoded.saturating_add(seed_match.occurrences.len() as u64);
                let limited = seed_match
                    .occurrences
                    .into_iter()
                    .take(
                        usize::try_from(seed_config.max_occurrences)
                            .unwrap_or(usize::MAX)
                            .saturating_add(1),
                    )
                    .map(|occurrence| crate::jma::SeedOccurrence {
                        contig_id: occurrence.contig_id,
                        position: occurrence.position,
                        reverse: occurrence.reverse,
                    })
                    .collect::<Vec<_>>();
                matches.insert(key.clone(), limited);
            }
        }
        for (key, seeds) in by_key {
            let Some(occurrences) = matches.get(&key) else {
                continue;
            };
            if occurrences.is_empty() {
                continue;
            }
            for seed in seeds {
                groups.push(SeedOccurrenceGroup {
                    seed,
                    k: level.k,
                    occurrences: occurrences.clone(),
                });
            }
        }
    }

    // `generate_anchors` sorts and truncates only after constructing every
    // occurrence. Keep the cap in this runner-owned path so a repeat-heavy
    // candidate cannot allocate an unbounded temporary anchor vector.
    let occurrence_limits = seed_configs
        .iter()
        .map(|config| (config.k, config.max_occurrences))
        .collect::<BTreeMap<_, _>>();
    let anchors = bounded_generate_anchors(
        &groups,
        &occurrence_limits,
        sensitivity.max_anchors_per_candidate,
    );
    let anchor_evidence = summarize_anchor_evidence(
        &anchors.anchors,
        &groups,
        seed_configs
            .iter()
            .map(|config| config.max_occurrences)
            .max()
            .unwrap_or(0),
        sensitivity.common_seed_candidate_occurrence_threshold,
    );
    let weighted_anchors = anchors
        .anchors
        .iter()
        .copied()
        .map(|anchor| {
            let class = anchor_evidence
                .iter()
                .find(|item| item.anchor == anchor)
                .map_or_else(
                    || class_for_anchor(anchor, false, false),
                    |item| class_for_anchor(anchor, item.seed.is_common, item.seed.is_repetitive),
                );
            WeightedAnchor::new(anchor, class)
        })
        .collect::<Vec<_>>();
    let chains = chain_weighted_anchors(
        &weighted_anchors,
        plasmid.len() as u64,
        ChainConfig::from_sensitivity(sensitivity).with_coordinate_model(coordinate_model),
    )
    .map_err(RunnerError::from)?;
    let chains_accepted = chains.len() as u64;
    let chains = chains
        .into_iter()
        .map(|chain| EvidencedChain {
            seed_evidence: evidence_for_weighted_chain(&chain, &anchor_evidence),
            chain,
        })
        .collect();
    Ok(LevelChainResult {
        chains,
        seed_keys_tested,
        occurrences_decoded,
        anchors_created: anchors.anchors.len() as u64,
        chains_accepted,
    })
}

fn select_seed_scheme(
    descriptors: &[SeedSchemeDescriptor],
    config: SeedSensitivity,
) -> Result<SeedSchemeDescriptor, RunnerError> {
    let mut matching = descriptors
        .iter()
        .copied()
        .filter(|descriptor| {
            descriptor.algorithm_id == 1
                && descriptor.span == u16::from(config.k)
                && descriptor.informative_bases == descriptor.span
                && descriptor.key_encoding == 1
                && descriptor.occurrence_encoding == 1
        })
        .collect::<Vec<_>>();
    matching.sort_by_key(|descriptor| (descriptor.density_parameter, descriptor.scheme_id));
    let Some(scheme) = matching.first().copied() else {
        return Err(RunnerError::SeedLevelMismatch {
            k: config.k,
            required_scale: config.scale,
            available_scale: None,
        });
    };
    let available_scale = u64::from(scheme.density_parameter);
    if available_scale > config.scale {
        return Err(RunnerError::SeedLevelMismatch {
            k: config.k,
            required_scale: config.scale,
            available_scale: Some(available_scale),
        });
    }
    Ok(scheme)
}

fn seed_key(seed: QuerySeed, k: u8) -> SeedKey {
    let width = usize::from(k).div_ceil(4);
    let bytes = seed.canonical_kmer.to_be_bytes();
    SeedKey {
        digest: seed.hash,
        verification: bytes[bytes.len().saturating_sub(width)..].to_vec(),
    }
}

fn class_for_anchor(anchor: Anchor, common: bool, repetitive: bool) -> AnchorClass {
    if repetitive || common {
        AnchorClass::Repetitive
    } else if anchor.k == crate::trace::seeds::PRIMARY_K {
        AnchorClass::SpecificK31
    } else if anchor.k == crate::trace::seeds::RESCUE_K {
        AnchorClass::SpecificK21
    } else {
        AnchorClass::SpecificPaired
    }
}

fn bounded_generate_anchors(
    groups: &[SeedOccurrenceGroup],
    occurrence_limits: &BTreeMap<u8, u32>,
    max_anchors: u32,
) -> AnchorSet {
    let max_anchors = usize::try_from(max_anchors).unwrap_or(usize::MAX);
    let mut anchors = Vec::with_capacity(max_anchors.min(groups.len()));
    let mut seen = HashSet::new();
    let mut repetitive_seeds = 0_u64;
    let mut skipped_occurrences = 0_u64;
    let mut truncated_anchors = 0_u64;
    for group in groups {
        if group.seed.hash == 0 {
            skipped_occurrences = skipped_occurrences
                .saturating_add(u64::try_from(group.occurrences.len()).unwrap_or(u64::MAX));
            continue;
        }
        let occurrence_limit = occurrence_limits.get(&group.k).copied().unwrap_or(0);
        if group.occurrences.len() > usize::try_from(occurrence_limit).unwrap_or(usize::MAX) {
            repetitive_seeds = repetitive_seeds.saturating_add(1);
            skipped_occurrences = skipped_occurrences
                .saturating_add(u64::try_from(group.occurrences.len()).unwrap_or(u64::MAX));
            continue;
        }
        for occurrence in &group.occurrences {
            let anchor = Anchor {
                query_position: group.seed.position,
                target_position: occurrence.position,
                contig_id: occurrence.contig_id,
                strand: anchor_strand(group.seed.reverse, occurrence.reverse),
                k: group.k,
                hash: group.seed.hash,
                canonical_kmer: group.seed.canonical_kmer,
                query_reverse: group.seed.reverse,
                target_reverse: occurrence.reverse,
            };
            let key = (
                anchor.query_position,
                anchor.target_position,
                anchor.contig_id,
                anchor.k,
                anchor.hash,
                anchor.canonical_kmer,
                matches!(anchor.strand, Strand::Reverse),
                anchor.query_reverse,
                anchor.target_reverse,
            );
            if seen.insert(key) {
                if anchors.len() < max_anchors {
                    anchors.push(anchor);
                } else {
                    truncated_anchors = truncated_anchors.saturating_add(1);
                }
            }
        }
    }
    anchors.sort_by_key(|anchor| {
        (
            anchor.contig_id,
            matches!(anchor.strand, Strand::Reverse),
            anchor.query_position,
            anchor.target_position,
            anchor.hash,
            anchor.canonical_kmer,
            anchor.k,
            anchor.query_reverse,
            anchor.target_reverse,
        )
    });
    AnchorSet {
        anchors,
        repetitive_seeds,
        skipped_occurrences,
        truncated_anchors,
    }
}

fn evidence_for_weighted_chain(
    chain: &WeightedAnchorChain,
    evidence: &[AnchorOccurrenceEvidence],
) -> SeedEvidence {
    let mut result = SeedEvidence::default();
    for weighted in &chain.anchors {
        let anchor = weighted.anchor;
        match weighted.class {
            AnchorClass::SpecificK31 => {
                result.primary_anchor_count = result.primary_anchor_count.saturating_add(1)
            }
            AnchorClass::SpecificK21 => {
                result.rescue_anchor_count = result.rescue_anchor_count.saturating_add(1)
            }
            _ => {}
        }
        if let Some(item) = evidence.iter().find(|item| item.anchor == anchor) {
            result.query_occurrence_max = result
                .query_occurrence_max
                .max(item.seed.query_occurrence_count);
            result.candidate_occurrence_max = result
                .candidate_occurrence_max
                .max(item.seed.candidate_occurrence_count);
            result.collection_document_frequency_max = match (
                result.collection_document_frequency_max,
                item.seed.collection_document_frequency,
            ) {
                (Some(left), Some(right)) => Some(left.max(right)),
                (None, value) | (value, None) => value,
            };
            if item.seed.is_repetitive {
                result.repetitive_seed_count = result.repetitive_seed_count.saturating_add(1);
            } else if item.seed.is_common {
                result.common_anchor_count = result.common_anchor_count.saturating_add(1);
            } else {
                result.nonrepetitive_anchor_count =
                    result.nonrepetitive_anchor_count.saturating_add(1);
            }
        }
    }
    result
}

#[allow(clippy::too_many_arguments)]
fn align_indexed_chains<A: TraceArchive>(
    archive: &A,
    query_id: &str,
    query: &[u8],
    candidate: &RankedCandidate,
    config: &TraceRunnerConfig,
    chains: Vec<EvidencedChain>,
    max_windows: u64,
    max_sequence_blocks: u64,
    workspace: &mut AlignmentWorkspace,
    seen_contigs: &mut HashSet<u32>,
    seen_windows: &mut HashSet<(u32, u64, u64, u8, u64, u64, isize)>,
    alignments: &mut Vec<BaseAlignment>,
    retry_metadata: &mut Vec<AlignmentRetryMetadata>,
    counters: &mut TracePerformanceCounters,
) -> Result<u64, RunnerError> {
    let mut attempted = 0_u64;
    let blocks_before = archive.metrics().resource.sequence_blocks_read;
    let metadata = archive.metadata().map_err(archive_error)?;
    for evidenced in chains {
        if attempted >= max_windows
            || archive
                .metrics()
                .resource
                .sequence_blocks_read
                .saturating_sub(blocks_before)
                >= max_sequence_blocks
        {
            break;
        }
        let weighted_chain = evidenced.chain;
        let chain = weighted_to_anchor_chain(&weighted_chain)?;
        let contig = metadata
            .contigs
            .iter()
            .find(|contig| contig.id == chain.contig_id)
            .ok_or(JmaError::UnknownContig(chain.contig_id))?;
        if seen_contigs.insert(contig.id) {
            counters.contigs_considered = counters.contigs_considered.saturating_add(1);
        }
        let range = chain_window(
            &chain,
            contig.length,
            query.len(),
            config.sensitivity.max_alignment_window_bases,
        )?;
        let coordinates = chain_alignment_coordinates(
            &chain,
            range,
            query.len() as u64,
            chaining_coordinate_model(config.topology_requested),
        )?;
        let strand_rank = match chain.strand {
            Strand::Forward => 0,
            Strand::Reverse => 1,
        };
        if !seen_windows.insert((
            contig.id,
            range.start,
            range.end,
            strand_rank,
            coordinates.query_start,
            coordinates.query_span,
            coordinates.diagonal_offset,
        )) {
            continue;
        }
        let request = SequenceRequest::new(contig.id, range.start, range.end, false)
            .map_err(archive_error)?;
        let sequence = archive
            .read_sequences(&[request])
            .map_err(archive_error)?
            .into_iter()
            .next()
            .ok_or_else(|| {
                RunnerError::InvalidConfig("archive returned no sequence slice".to_string())
            })?
            .bases;
        counters.windows_retrieved = counters.windows_retrieved.saturating_add(1);
        if let Some(alignment) = exact_full_query_alignment(
            query_id,
            candidate,
            query,
            &sequence,
            range,
            &contig.name,
            &weighted_chain,
            &coordinates,
            evidenced.seed_evidence.clone(),
            config.sensitivity.alignment.match_score,
        )? {
            counters.alignments_succeeded = counters.alignments_succeeded.saturating_add(1);
            alignments.push(alignment);
            continue;
        }
        let corridor = build_chain_corridor(
            &weighted_chain,
            range,
            &coordinates,
            query.len() as u64,
            config
                .sensitivity
                .alignment
                .band_width
                .min(crate::trace::alignment::DEFAULT_RETRY_WIDTHS[0]),
            alignment_half_width_cap(config, coordinates.query_span),
        )?;
        attempted = attempted.saturating_add(1);
        align_window(
            query_id,
            candidate,
            workspace,
            query,
            &sequence,
            range.start,
            &contig.name,
            config,
            chain.score,
            evidenced.seed_evidence,
            Some(chain.strand),
            coordinates.query_start,
            coordinates.query_span,
            coordinates.diagonal_offset,
            Some(&corridor),
            alignments,
            retry_metadata,
            counters,
        )?;
    }
    Ok(attempted)
}

#[allow(clippy::too_many_arguments)]
fn exact_full_query_alignment(
    query_id: &str,
    candidate: &RankedCandidate,
    query: &[u8],
    target: &[u8],
    target_range: SequenceRange,
    contig_id: &str,
    chain: &WeightedAnchorChain,
    coordinates: &ChainAlignmentCoordinates,
    seed_evidence: SeedEvidence,
    match_score: i32,
) -> Result<Option<BaseAlignment>, RunnerError> {
    if chain.origin_crossing
        || coordinates.query_start != 0
        || coordinates.query_span != query.len() as u64
    {
        return Ok(None);
    }
    let mut anchors = Vec::with_capacity(chain.anchors.len());
    for weighted in &chain.anchors {
        let Some(target_position) = weighted
            .anchor
            .target_position
            .checked_sub(target_range.start)
        else {
            return Ok(None);
        };
        if target_position >= target_range.len() {
            return Ok(None);
        }
        anchors.push(Anchor {
            target_position,
            ..weighted.anchor
        });
    }
    let blocks = extend_and_merge_anchors(query, target, &anchors)
        .map_err(|error| RunnerError::InvalidConfig(error.to_string()))?;
    let Some(block) = blocks.into_iter().find(|block| {
        block.query_interval.start == 0 && block.query_interval.end == query.len() as u64
    }) else {
        return Ok(None);
    };
    let target_interval = BaseInterval::new(
        block
            .target_interval
            .start
            .checked_add(target_range.start)
            .ok_or_else(|| RunnerError::InvalidConfig("exact block target overflow".to_string()))?,
        block
            .target_interval
            .end
            .checked_add(target_range.start)
            .ok_or_else(|| RunnerError::InvalidConfig("exact block target overflow".to_string()))?,
    )
    .map_err(|error| RunnerError::InvalidConfig(error.to_string()))?;
    let matches = block.len();
    let score = i64::from(match_score).saturating_mul(i64::try_from(matches).unwrap_or(i64::MAX));
    Ok(Some(BaseAlignment {
        alignment_id: String::new(),
        plasmid_id: query_id.to_string(),
        metagenome_id: candidate.metagenome_id().to_string(),
        contig_id: contig_id.to_string(),
        strand: block.strand,
        query_segments: vec![block.query_interval],
        target_interval,
        query_length: matches,
        target_length: matches,
        origin_crossing: false,
        score,
        matches,
        substitutions: 0,
        insertions: 0,
        deletions: 0,
        cigar: block.cigar(),
        edit_script: block
            .edit_runs_checked()
            .map_err(|error| RunnerError::InvalidConfig(error.to_string()))?,
        chain_score: chain.score,
        identity: 1.0,
        seed_evidence,
        primary_supported_bases: 0,
        secondary_supported_bases: 0,
        newly_supported_bases: 0,
        role: AlignmentRole::AlternativeMapping,
        primary: false,
    }))
}

fn weighted_to_anchor_chain(chain: &WeightedAnchorChain) -> Result<AnchorChain, RunnerError> {
    let Some(first) = chain.anchors.first() else {
        return Err(RunnerError::InvalidConfig(
            "accepted weighted chain has no anchors".to_string(),
        ));
    };
    Ok(AnchorChain {
        contig_id: chain.contig_id,
        strand: chain.strand,
        k: first.anchor.k,
        anchors: chain.anchors.iter().map(|anchor| anchor.anchor).collect(),
        score: chain.score,
        query_interval: chain.query_interval,
        query_segments: chain.query_segments.clone(),
        target_interval: chain.target_interval,
        origin_crossing: chain.origin_crossing,
        linear_query_start: chain.linear_query_start,
        linear_query_end: chain.linear_query_end,
    })
}

fn build_chain_corridor(
    chain: &WeightedAnchorChain,
    range: SequenceRange,
    coordinates: &ChainAlignmentCoordinates,
    query_length: u64,
    safety_margin: u32,
    max_half_width: u32,
) -> Result<AlignmentCorridor, RunnerError> {
    let mut anchors = chain
        .anchors
        .iter()
        .filter_map(|weighted| {
            let span = u32::try_from(weighted.span).ok()?;
            let mut query_position = weighted.anchor.query_position;
            while query_position < coordinates.query_start {
                query_position = query_position.checked_add(query_length)?;
            }
            let query_position = query_position.checked_sub(coordinates.query_start)?;
            if query_position >= coordinates.query_span {
                return None;
            }
            let target_position = match chain.strand {
                Strand::Forward => weighted.anchor.target_position.checked_sub(range.start)?,
                Strand::Reverse => range
                    .end
                    .checked_sub(weighted.anchor.target_position.checked_add(weighted.span)?)?,
            };
            Some(CorridorAnchor::new(query_position, target_position, span))
        })
        .collect::<Vec<_>>();
    anchors.sort_by_key(|anchor| (anchor.query_position, anchor.target_position));
    anchors.dedup_by_key(|anchor| (anchor.query_position, anchor.target_position));
    AlignmentCorridor::new(anchors, safety_margin)
        .map(|corridor| corridor.with_max_half_width(max_half_width))
        .map_err(RunnerError::Alignment)
}

const fn chaining_coordinate_model(topology: TopologyRequested) -> CoordinateModel {
    match topology {
        TopologyRequested::Linear => CoordinateModel::Linear,
        TopologyRequested::Circular | TopologyRequested::Auto | TopologyRequested::Unknown => {
            CoordinateModel::Wrap
        }
    }
}

fn chain_window(
    chain: &AnchorChain,
    contig_length: u64,
    query_length: usize,
    max_window_bases: u64,
) -> Result<SequenceRange, RunnerError> {
    let query_length = u64::try_from(query_length)
        .map_err(|_| RunnerError::InvalidConfig("query length exceeds u64".to_string()))?;
    if query_length > max_window_bases {
        return Err(RunnerError::WindowTooLarge {
            query_length: usize::try_from(query_length).unwrap_or(usize::MAX),
            max_window: usize::try_from(max_window_bases).unwrap_or(usize::MAX),
        });
    }
    if chain.target_interval.end > contig_length {
        return Err(RunnerError::InvalidConfig(
            "seed chain exceeds the declared contig length".to_string(),
        ));
    }
    let span = max_window_bases.min(contig_length);
    let desired_start = chain.target_interval.start.saturating_sub(query_length);
    let desired_end = chain
        .target_interval
        .end
        .saturating_add(query_length)
        .min(contig_length);
    if desired_end.saturating_sub(desired_start) <= span {
        return SequenceRange::new(desired_start, desired_end).map_err(RunnerError::from);
    }
    let midpoint = chain
        .target_interval
        .start
        .saturating_add(chain.target_interval.len() / 2);
    let mut start = midpoint.saturating_sub(span / 2);
    start = start.min(contig_length.saturating_sub(span));
    SequenceRange::new(start, start.saturating_add(span)).map_err(RunnerError::from)
}

struct ChainAlignmentCoordinates {
    query_start: u64,
    query_span: u64,
    diagonal_offset: isize,
}

/// Translate chain coordinates into the linearized query/target window used
/// by the banded aligner.  JMA target positions are contig-local, while the
/// alignment target is the retrieved window and the query is circular.  The
/// query origin is backed up from the first chain anchor by the matching
/// target-window prefix, with reverse chains measured in the reversed target
/// orientation.
fn chain_alignment_coordinates(
    chain: &AnchorChain,
    range: SequenceRange,
    query_length: u64,
    coordinate_model: CoordinateModel,
) -> Result<ChainAlignmentCoordinates, RunnerError> {
    if query_length == 0 || range.is_empty() {
        return Err(RunnerError::InvalidConfig(
            "chain alignment requires non-empty query and target windows".to_string(),
        ));
    }
    let chain_span = chain
        .linear_query_end
        .checked_sub(chain.linear_query_start)
        .ok_or_else(|| {
            RunnerError::InvalidConfig("chain query coordinates are reversed".to_string())
        })?;
    if chain_span == 0 || chain_span > query_length {
        return Err(RunnerError::InvalidConfig(
            "chain query span is outside the plasmid".to_string(),
        ));
    }
    let target_offset = match chain.strand {
        Strand::Forward => chain
            .target_interval
            .start
            .checked_sub(range.start)
            .ok_or_else(|| {
                RunnerError::InvalidConfig(
                    "chain target starts before its retrieved window".to_string(),
                )
            })?,
        Strand::Reverse => range
            .end
            .checked_sub(chain.target_interval.end)
            .ok_or_else(|| {
                RunnerError::InvalidConfig(
                    "reverse chain ends after its retrieved window".to_string(),
                )
            })?,
    };
    if matches!(
        coordinate_model,
        CoordinateModel::Linear | CoordinateModel::Undetermined
    ) {
        let query_anchor = chain.linear_query_start % query_length;
        let diagonal_offset = isize::try_from(i128::from(target_offset) - i128::from(query_anchor))
            .map_err(|_| {
                RunnerError::InvalidConfig(
                    "linear chain diagonal exceeds alignment coordinate range".to_string(),
                )
            })?;
        return Ok(ChainAlignmentCoordinates {
            query_start: 0,
            query_span: query_length,
            diagonal_offset,
        });
    }
    let left_context = target_offset.min(query_length - chain_span);
    let query_origin = chain.linear_query_start % query_length;
    let query_start = if left_context <= query_origin {
        query_origin - left_context
    } else {
        query_length - (left_context - query_origin)
    };
    let required_span = left_context
        .checked_add(chain_span)
        .ok_or_else(|| RunnerError::InvalidConfig("chain query span overflowed".to_string()))?;
    let query_span = range.len().max(required_span).min(query_length);
    let diagonal_offset = isize::try_from(target_offset - left_context).map_err(|_| {
        RunnerError::InvalidConfig(
            "chain target offset exceeds alignment coordinate range".to_string(),
        )
    })?;
    Ok(ChainAlignmentCoordinates {
        query_start,
        query_span,
        diagonal_offset,
    })
}

struct FinalEvidence {
    topology: TopologyAssessment,
    primary_mosaic: FragmentMosaicSummary,
    coverage: CoverageSummary,
}

fn finalize_evidence(
    query_length: u64,
    config: &TraceRunnerConfig,
    alignments: &mut [BaseAlignment],
) -> Result<FinalEvidence, RunnerError> {
    assign_alignment_ids(alignments);
    let topology = assess_topology(
        query_length,
        config.topology_requested,
        config.topology_margin_bases,
        alignments,
    )?;
    let primary_mosaic = selected_mosaic(&topology).clone();
    let mut selected_by_id = primary_mosaic
        .alignment_evidence
        .iter()
        .map(|evidence| (evidence.alignment_id.as_str(), evidence))
        .collect::<BTreeMap<_, _>>();
    if let Some(wrap) = topology.wrap_model.as_ref() {
        for evidence in &wrap.mosaic.alignment_evidence {
            selected_by_id
                .entry(evidence.alignment_id.as_str())
                .or_insert(evidence);
        }
    }
    for alignment in alignments.iter_mut() {
        alignment.primary = false;
        alignment.primary_supported_bases = 0;
        alignment.secondary_supported_bases = 0;
        alignment.newly_supported_bases = 0;
        alignment.role = AlignmentRole::AlternativeMapping;
        if let Some(evidence) = selected_by_id.get(alignment.alignment_id.as_str()) {
            alignment.primary_supported_bases = evidence.primary_supported_bases;
            alignment.secondary_supported_bases = evidence.secondary_supported_bases;
            alignment.newly_supported_bases = evidence.newly_supported_bases;
            alignment.role = evidence.role;
            alignment.primary = evidence.primary_supported_bases > 0;
        }
    }
    let coverage = coverage_from_mosaic(&primary_mosaic, query_length);
    Ok(FinalEvidence {
        topology,
        primary_mosaic,
        coverage,
    })
}

fn selected_mosaic(topology: &TopologyAssessment) -> &FragmentMosaicSummary {
    if topology.coordinate_model == CoordinateModel::Wrap {
        topology
            .wrap_model
            .as_ref()
            .map_or(&topology.linear_model.mosaic, |model| &model.mosaic)
    } else {
        &topology.linear_model.mosaic
    }
}

fn coverage_from_mosaic(mosaic: &FragmentMosaicSummary, query_length: u64) -> CoverageSummary {
    let mut largest_gap = mosaic
        .unsupported_gaps
        .iter()
        .map(|gap| gap.length)
        .max()
        .unwrap_or(0);
    if mosaic.coordinate_model == CoordinateModel::Wrap
        && mosaic.unsupported_gaps.len() >= 2
        && let (Some(first), Some(last)) = (
            mosaic.unsupported_gaps.first(),
            mosaic.unsupported_gaps.last(),
        )
        && first.interval.start == 0
        && last.interval.end == query_length
    {
        largest_gap = largest_gap.max(first.length.saturating_add(last.length));
    }
    CoverageSummary {
        plasmid_length: query_length,
        supported_bases: mosaic.base_covered_bases,
        supported_fraction: mosaic.base_coverage_fraction,
        primary_intervals: mosaic.covered_intervals.clone(),
        secondary_intervals: Vec::new(),
        gaps: mosaic.unsupported_gaps.clone(),
        largest_gap,
    }
}

fn evidence_warnings(finalized: &FinalEvidence) -> Vec<String> {
    let mut warnings = Vec::new();
    if finalized.topology.coordinate_model == CoordinateModel::Undetermined {
        warnings.push(
            "coordinate model is undetermined; wrapped coordinates do not establish biological topology"
                .to_string(),
        );
    }
    if finalized.primary_mosaic.supporting_contigs.len() > 1 {
        warnings.push(
            "supporting contigs are independent sequence observations and are not physically joined"
                .to_string(),
        );
    }
    if finalized.primary_mosaic.repeat_only_supported_bases > 0
        && finalized.primary_mosaic.nonrepetitive_supported_bases == 0
    {
        warnings.push("all retained support is common or repetitive sequence evidence".to_string());
    }
    warnings
}

fn retain_best_alignments(alignments: &mut Vec<BaseAlignment>, limit: usize) {
    alignments.sort_by(compare_alignments);
    alignments.dedup_by(|right, left| same_alignment(right, left));
    alignments.truncate(limit);
}

fn retain_new_alignments(
    alignments: &mut Vec<BaseAlignment>,
    protected_count: usize,
    limit: usize,
) {
    let protected_count = protected_count.min(alignments.len());
    let mut newly_added = alignments.split_off(protected_count);
    newly_added.retain(|alignment| {
        !alignments
            .iter()
            .any(|existing| same_alignment(alignment, existing))
    });
    retain_best_alignments(&mut newly_added, limit);
    alignments.extend(newly_added);
}

fn same_alignment(left: &BaseAlignment, right: &BaseAlignment) -> bool {
    left.contig_id == right.contig_id
        && left.strand == right.strand
        && left.target_interval == right.target_interval
        && left.query_segments == right.query_segments
}

fn compare_alignments(left: &BaseAlignment, right: &BaseAlignment) -> Ordering {
    right
        .score
        .cmp(&left.score)
        .then_with(|| right.matches.cmp(&left.matches))
        .then_with(|| right.chain_score.cmp(&left.chain_score))
        .then_with(|| left.contig_id.cmp(&right.contig_id))
        .then_with(|| left.target_interval.start.cmp(&right.target_interval.start))
        .then_with(|| left.target_interval.end.cmp(&right.target_interval.end))
        .then_with(|| strand_key(left).cmp(&strand_key(right)))
}

fn strand_key(alignment: &BaseAlignment) -> u8 {
    match alignment.strand {
        crate::trace::model::Strand::Forward => 0,
        crate::trace::model::Strand::Reverse => 1,
    }
}

#[allow(clippy::too_many_arguments)]
fn failed_result(
    plasmid_id: &str,
    candidate: &RankedCandidate,
    config: &TraceRunnerConfig,
    stage: &str,
    code: &str,
    message: &str,
    retryable: bool,
    resource: Option<&str>,
) -> TraceMetagenomeResult {
    failed_result_with_metrics(
        plasmid_id,
        candidate,
        config,
        stage,
        code,
        message,
        retryable,
        resource,
        ResourceMetrics::default(),
    )
}

#[allow(clippy::too_many_arguments)]
fn failed_result_with_metrics(
    plasmid_id: &str,
    candidate: &RankedCandidate,
    config: &TraceRunnerConfig,
    stage: &str,
    code: &str,
    message: &str,
    retryable: bool,
    resource: Option<&str>,
    resource_metrics: ResourceMetrics,
) -> TraceMetagenomeResult {
    TraceMetagenomeResult {
        schema_version: crate::trace::TRACE_JSON_SCHEMA_VERSION.to_string(),
        run_id: String::new(),
        plasmid_id: plasmid_id.to_string(),
        metagenome_id: candidate.metagenome_id().to_string(),
        query_kind: config.query_kind,
        topology_requested: config.topology_requested,
        coordinate_model: match config.topology_requested {
            TopologyRequested::Linear => CoordinateModel::Linear,
            TopologyRequested::Circular => CoordinateModel::Wrap,
            TopologyRequested::Auto | TopologyRequested::Unknown => CoordinateModel::Undetermined,
        },
        topology_evidence: if config.topology_requested == TopologyRequested::Unknown {
            TopologyEvidence::Undetermined
        } else {
            TopologyEvidence::Insufficient
        },
        algorithms: crate::trace::config::algorithm_identifiers(),
        algorithm: crate::trace::config::TraceAlgorithmMetadata::for_sensitivity(
            config.sensitivity.clone(),
        ),
        status: TraceStatus::Failed,
        candidate: Some(candidate.candidate.clone()),
        router_candidate: None,
        alignments: Vec::new(),
        primary_fragment_mosaic: None,
        topology: None,
        rescue_rounds: Vec::new(),
        stages: Vec::new(),
        alignment_retries: Vec::new(),
        performance_counters: CandidatePerformanceCounters::default(),
        coverage: None,
        warnings: Vec::new(),
        failures: vec![TraceFailure {
            stage: stage.to_string(),
            code: code.to_string(),
            message: message.to_string(),
            resource: resource.map(redact_resource),
            retryable,
        }],
        archive_metrics: None,
        resource_metrics,
    }
}

fn redact_resource(resource: &str) -> String {
    crate::resource::ResourceLocator::parse(resource)
        .map(|locator| locator.redacted())
        .unwrap_or_else(|_| "<invalid-resource>".to_string())
}

fn raw_error_code(error: &RawError) -> &'static str {
    match error {
        RawError::Resource(error) => resource_error_code(error),
        RawError::Parse { .. } => "raw_parse_error",
        RawError::Io { .. } => "raw_io_error",
    }
}

fn resource_error_code(error: &ResourceError) -> &'static str {
    match error {
        ResourceError::HttpStatus {
            status: 401 | 403, ..
        } => "access_denied",
        ResourceError::HttpStatus { status: 404, .. } => "missing_object",
        ResourceError::HttpStatus { status, .. }
            if matches!(*status, 408 | 425 | 429 | 500..=599) =>
        {
            "retryable_server_error"
        }
        ResourceError::HttpStatus { .. } => "http_status_error",
        ResourceError::Timeout { .. } => "timeout",
        ResourceError::Transport { .. } => "retryable_resource_error",
        ResourceError::CacheIdentityChanged(_) => "stale_cache",
        ResourceError::RangeUnsupported(_) => "range_unsupported",
        ResourceError::InvalidLocator(_) => "invalid_resource",
        ResourceError::UnsupportedScheme(_) => "unsupported_resource",
        ResourceError::RangeOverflow { .. } | ResourceError::RangeOutOfBounds { .. } => {
            "resource_range_error"
        }
        ResourceError::Io { .. } => "resource_io_error",
    }
}

fn resource_error_retryable(error: &ResourceError) -> bool {
    match error {
        ResourceError::HttpStatus { status, .. } => {
            matches!(*status, 408 | 425 | 429 | 500..=599)
        }
        ResourceError::Timeout { .. } | ResourceError::Transport { .. } => true,
        ResourceError::InvalidLocator(_)
        | ResourceError::UnsupportedScheme(_)
        | ResourceError::RangeOverflow { .. }
        | ResourceError::RangeOutOfBounds { .. }
        | ResourceError::CacheIdentityChanged(_)
        | ResourceError::RangeUnsupported(_)
        | ResourceError::Io { .. } => false,
    }
}

impl RunnerError {
    fn code(&self) -> &'static str {
        match self {
            Self::Candidate(_) => "candidate_error",
            Self::Catalog(_) => "catalog_error",
            Self::Resource(error) => resource_error_code(error),
            Self::Raw(error) => raw_error_code(error),
            Self::Jma(error) => match error {
                JmaError::Resource(error) => resource_error_code(error),
                JmaError::ChecksumMismatch(_) => "jma_checksum_mismatch",
                JmaError::CorruptSection(_) => "jma_corrupt_section",
                JmaError::UnsupportedVersion(_) => "jma_unsupported_version",
                JmaError::InvalidMagic => "jma_invalid_magic",
                _ => "jma_error",
            },
            Self::Alignment(_) => "alignment_error",
            Self::Seed(_) => "seed_error",
            Self::Chain(_) => "chain_error",
            Self::Mosaic(_) => "coverage_error",
            Self::SeedLevelMismatch { .. } => "seed_level_mismatch",
            Self::WindowTooLarge { .. } => "alignment_window_too_large",
            Self::DatabaseCache { .. } => "candidate_database_cache_error",
            Self::InvalidQuery(_) | Self::InvalidConfig(_) => "invalid_configuration",
        }
    }

    fn retryable(&self) -> bool {
        match self {
            Self::Resource(error) => resource_error_retryable(error),
            Self::Raw(RawError::Resource(error)) => resource_error_retryable(error),
            Self::Jma(JmaError::Resource(error)) => resource_error_retryable(error),
            _ => false,
        }
    }
}

#[derive(Debug, Error)]
pub enum RunnerError {
    #[error("candidate search failed: {0}")]
    Candidate(#[from] CandidateError),
    #[error("catalog failed: {0}")]
    Catalog(#[from] CatalogError),
    #[error("resource access failed: {0}")]
    Resource(#[from] ResourceError),
    #[error("raw assembly failed: {0}")]
    Raw(#[from] RawError),
    #[error("JMA processing failed: {0}")]
    Jma(#[from] JmaError),
    #[error("alignment failed: {0}")]
    Alignment(#[from] AlignmentError),
    #[error("seed extraction failed: {0}")]
    Seed(#[from] SeedError),
    #[error("seed chaining failed: {0}")]
    Chain(#[from] ChainError),
    #[error("coverage selection failed: {0}")]
    Mosaic(#[from] MosaicError),
    #[error(
        "JMA seed level k={k} cannot serve required scale {required_scale}; available densest scale is {available_scale:?}"
    )]
    SeedLevelMismatch {
        k: u8,
        required_scale: u64,
        available_scale: Option<u64>,
    },
    #[error("alignment window of {query_length} bases exceeds configured maximum {max_window}")]
    WindowTooLarge {
        query_length: usize,
        max_window: usize,
    },
    #[error("invalid trace query: {0}")]
    InvalidQuery(String),
    #[error("invalid trace runner configuration: {0}")]
    InvalidConfig(String),
    #[error("candidate database cache failed for {locator}: {message}")]
    DatabaseCache { locator: String, message: String },
}

#[cfg(test)]
mod coordinate_tests {
    use super::*;
    use crate::trace::model::BaseInterval;

    fn dna(length: usize) -> Vec<u8> {
        let mut state = 0xD1B54A32D192ED03u64;
        (0..length)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                b"ACGT"[((state >> 62) & 3) as usize]
            })
            .collect()
    }

    fn chain(strand: Strand, target_start: u64) -> AnchorChain {
        AnchorChain {
            contig_id: 7,
            strand,
            k: 31,
            anchors: Vec::new(),
            score: 100,
            query_interval: BaseInterval::new(18_856, 18_887).unwrap(),
            query_segments: vec![BaseInterval::new(18_856, 18_887).unwrap()],
            target_interval: BaseInterval::new(target_start, target_start + 31).unwrap(),
            origin_crossing: false,
            linear_query_start: 18_856,
            linear_query_end: 18_887,
        }
    }

    #[test]
    fn split_chain_translation_uses_query_start_and_local_target_offset() {
        let range = SequenceRange::new(0, 20_856).unwrap();
        let coordinates = chain_alignment_coordinates(
            &chain(Strand::Forward, 5_000),
            range,
            94_281,
            CoordinateModel::Wrap,
        )
        .unwrap();
        assert_eq!(coordinates.query_start, 13_856);
        assert_eq!(coordinates.query_span, 20_856);
        assert_eq!(coordinates.diagonal_offset, 0);
    }

    #[test]
    fn reverse_chain_translation_uses_reversed_window_origin() {
        let range = SequenceRange::new(10_000, 30_856).unwrap();
        let coordinates = chain_alignment_coordinates(
            &chain(Strand::Reverse, 15_000),
            range,
            94_281,
            CoordinateModel::Wrap,
        )
        .unwrap();
        assert_eq!(coordinates.query_start, 3_031);
        assert_eq!(coordinates.query_span, 20_856);
        assert_eq!(coordinates.diagonal_offset, 0);
    }

    #[test]
    fn translated_split_chain_recovers_global_query_interval() {
        let query = dna(400);
        let target = &query[120..320];
        let chain = AnchorChain {
            contig_id: 1,
            strand: Strand::Forward,
            k: 31,
            anchors: Vec::new(),
            score: 100,
            query_interval: BaseInterval::new(150, 181).unwrap(),
            query_segments: vec![BaseInterval::new(150, 181).unwrap()],
            target_interval: BaseInterval::new(30, 61).unwrap(),
            origin_crossing: false,
            linear_query_start: 150,
            linear_query_end: 181,
        };
        let coordinates = chain_alignment_coordinates(
            &chain,
            SequenceRange::new(0, 200).unwrap(),
            400,
            CoordinateModel::Wrap,
        )
        .unwrap();
        let mut workspace = AlignmentWorkspace::new();
        let result = workspace
            .align_circular(
                &query,
                coordinates.query_start,
                coordinates.query_span,
                target,
                0,
                Strand::Forward,
                AlignmentOptions::new(SensitivityConfig::default().alignment)
                    .with_diagonal_offset(coordinates.diagonal_offset),
            )
            .unwrap();
        assert_eq!(
            result.query_segments,
            vec![BaseInterval::new(120, 320).unwrap()]
        );
        assert_eq!(result.matches, 200);
    }
}

#[cfg(test)]
mod failure_tests {
    use super::*;

    fn alignment(id: &str, score: i64, target_start: u64) -> BaseAlignment {
        BaseAlignment {
            alignment_id: id.to_string(),
            plasmid_id: "query".to_string(),
            metagenome_id: "sample".to_string(),
            contig_id: "contig".to_string(),
            strand: Strand::Forward,
            query_segments: vec![BaseInterval {
                start: target_start,
                end: target_start + 4,
            }],
            target_interval: BaseInterval {
                start: target_start,
                end: target_start + 4,
            },
            query_length: 100,
            target_length: 4,
            origin_crossing: false,
            score,
            matches: 4,
            substitutions: 0,
            insertions: 0,
            deletions: 0,
            cigar: "4=".to_string(),
            edit_script: Vec::new(),
            chain_score: score,
            identity: 1.0,
            seed_evidence: SeedEvidence::default(),
            primary_supported_bases: 0,
            secondary_supported_bases: 0,
            newly_supported_bases: 0,
            role: crate::trace::model::AlignmentRole::AlternativeMapping,
            primary: false,
        }
    }

    #[test]
    fn resource_failure_codes_and_retryability_are_stable() {
        let cases = [
            (
                ResourceError::HttpStatus {
                    locator: "https://example.test/object".to_string(),
                    status: 403,
                },
                "access_denied",
                false,
            ),
            (
                ResourceError::HttpStatus {
                    locator: "https://example.test/object".to_string(),
                    status: 404,
                },
                "missing_object",
                false,
            ),
            (
                ResourceError::Timeout {
                    locator: "https://example.test/object".to_string(),
                },
                "timeout",
                true,
            ),
            (
                ResourceError::HttpStatus {
                    locator: "https://example.test/object".to_string(),
                    status: 503,
                },
                "retryable_server_error",
                true,
            ),
            (
                ResourceError::CacheIdentityChanged("https://example.test/object".to_string()),
                "stale_cache",
                false,
            ),
        ];
        for (error, expected_code, expected_retryable) in cases {
            let runner_error = RunnerError::Resource(error);
            assert_eq!(runner_error.code(), expected_code);
            assert_eq!(runner_error.retryable(), expected_retryable);
        }
    }

    #[test]
    fn jma_corruption_failure_codes_remain_specific() {
        assert_eq!(
            RunnerError::Jma(JmaError::ChecksumMismatch("header".to_string())).code(),
            "jma_checksum_mismatch"
        );
        assert_eq!(
            RunnerError::Jma(JmaError::UnsupportedVersion(99)).code(),
            "jma_unsupported_version"
        );
        assert_eq!(
            RunnerError::Jma(JmaError::InvalidMagic).code(),
            "jma_invalid_magic"
        );
        assert_eq!(
            RunnerError::Jma(JmaError::CorruptSection("truncated".to_string())).code(),
            "jma_corrupt_section"
        );
    }

    #[test]
    fn rescue_retention_preserves_earlier_alignments_and_caps_new_evidence() {
        let earlier = alignment("earlier-exact", 1, 0);
        let later_low = alignment("later-low", 10, 10);
        let later_high = alignment("later-high", 20, 20);
        let mut alignments = vec![earlier, later_low, later_high];

        retain_new_alignments(&mut alignments, 1, 1);

        assert_eq!(alignments.len(), 2);
        assert_eq!(alignments[0].alignment_id, "earlier-exact");
        assert_eq!(alignments[1].alignment_id, "later-high");
    }
}

#[cfg(test)]
mod memory_cap_tests {
    use super::*;

    #[test]
    fn long_query_retry_width_is_capped_before_matrix_allocation() {
        let config = TraceRunnerConfig {
            threads: 4,
            io_concurrency: 2,
            memory_budget_bytes: 4 * 1024 * 1024 * 1024,
            ..TraceRunnerConfig::default()
        };
        assert!(alignment_half_width_cap(&config, 94_281) >= 512);
        assert!(alignment_half_width_cap(&config, 94_281) < 1024);
        assert!(alignment_half_width_cap(&config, 10_000) >= 512);
        assert!(candidate_memory_estimate(&config, 94_281).total_bytes() <= 4 * 1024 * 1024 * 1024);
    }
}

#[cfg(test)]
mod cascade_plan_tests {
    use super::*;
    use crate::trace::config::{SensitivityConfig, SensitivityProfile};

    #[test]
    fn sensitive_profile_does_not_repeat_same_k31_density_before_mixed_rescue() {
        let rounds = positional_round_configs(&SensitivityConfig::for_profile(
            SensitivityProfile::Sensitive,
        ));
        assert_eq!(rounds.len(), 2);
        assert_eq!(rounds[0].len(), 1);
        assert_eq!(rounds[1].len(), 2);
        assert_eq!(rounds[1][0].k, crate::trace::seeds::PRIMARY_K);
        assert_eq!(rounds[1][1].k, crate::trace::seeds::RESCUE_K);
    }

    #[test]
    fn exact_gear_strand_is_relative_to_query_and_target_canonicalization() {
        assert_eq!(
            relative_fragment_orientation(
                FragmentOrientation::Forward,
                FragmentOrientation::Reverse,
            ),
            FragmentOrientation::Reverse
        );
        assert_eq!(
            relative_fragment_orientation(
                FragmentOrientation::Reverse,
                FragmentOrientation::Reverse,
            ),
            FragmentOrientation::Forward
        );
    }
}
