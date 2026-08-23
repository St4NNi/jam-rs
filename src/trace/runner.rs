//! Candidate-to-sequence trace runner.
//!
//! The runner owns orchestration and bounded resource use.  It queries the
//! existing `.jam` candidate index once, resolves candidate resources through
//! the catalog, retrieves indexed windows when JMA is available, and invokes
//! the alignment kernel.  Its output remains evidence for downstream
//! confirmation; it never upgrades a candidate to a biological presence call.

use crate::jma::reader::JmaReader;
use crate::jma::{ArchiveReader, JmaError, SequenceRange};
use crate::resource::{RangeReader, ResourceError, ResourceMetrics, ResourceOpenOptions};
use crate::trace::alignment::{AlignmentError, AlignmentOptions, AlignmentWorkspace};
use crate::trace::anchors::{
    AnchorOccurrenceEvidence, SeedOccurrenceGroup, generate_anchors, summarize_anchor_evidence,
};
use crate::trace::catalog::{CatalogEntry, CatalogError, TraceCatalog};
use crate::trace::chain::{AnchorChain, ChainConfig, ChainError, chain_anchors};
use crate::trace::config::{SeedSensitivity, SensitivityConfig};
use crate::trace::model::{
    AlignmentRole, BaseAlignment, BaseInterval, CandidatePerformanceCounters, CoordinateModel,
    CoverageSummary, FragmentMosaicSummary, InputResource, QueryKind, RescueRoundMetrics,
    SeedEvidence, Strand, TopologyAssessment, TopologyEvidence, TopologyRequested, TraceFailure,
    TraceMetagenomeResult, TraceStatus,
};
use crate::trace::mosaic::{MosaicError, assess_topology, assign_alignment_ids};
use crate::trace::raw::{AssemblyResource, RawAssembly, RawError, open_resource};
use crate::trace::screen::{
    CandidateError, CandidateSearchConfig, CandidateSearchResult, CandidateSearcher,
    RankedCandidate,
};
use crate::trace::seeds::{SeedError, extract_seed_level, extract_seed_level_in_intervals};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashSet};
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
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
}

impl Default for TraceRunnerConfig {
    fn default() -> Self {
        Self {
            sensitivity: SensitivityConfig::default(),
            candidates: CandidateSearchConfig::default(),
            resources: ResourceOpenOptions::default(),
            threads: 1,
            io_concurrency: 1,
            max_alignments_per_candidate: 256,
            query_kind: QueryKind::Unknown,
            topology_requested: TopologyRequested::Auto,
            topology_margin_bases: SensitivityConfig::default().auto_topology_margin_bases,
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

fn subtract_resource_metrics(after: ResourceMetrics, before: ResourceMetrics) -> ResourceMetrics {
    ResourceMetrics {
        metadata_requests: after
            .metadata_requests
            .saturating_sub(before.metadata_requests),
        head_requests: after.head_requests.saturating_sub(before.head_requests),
        get_requests: after.get_requests.saturating_sub(before.get_requests),
        range_requests: after.range_requests.saturating_sub(before.range_requests),
        stream_requests: after.stream_requests.saturating_sub(before.stream_requests),
        requested_bytes: after.requested_bytes.saturating_sub(before.requested_bytes),
        returned_bytes: after.returned_bytes.saturating_sub(before.returned_bytes),
        decoded_bytes: after.decoded_bytes.saturating_sub(before.decoded_bytes),
        remote_bytes: after.remote_bytes.saturating_sub(before.remote_bytes),
        cache_bytes: after.cache_bytes.saturating_sub(before.cache_bytes),
        cache_hits: after.cache_hits.saturating_sub(before.cache_hits),
        cache_misses: after.cache_misses.saturating_sub(before.cache_misses),
        cache_evictions: after.cache_evictions.saturating_sub(before.cache_evictions),
        stale_cache_rejections: after
            .stale_cache_rejections
            .saturating_sub(before.stale_cache_rejections),
        retries: after.retries.saturating_sub(before.retries),
        full_object_fallbacks: after
            .full_object_fallbacks
            .saturating_sub(before.full_object_fallbacks),
        seed_buckets_read: after
            .seed_buckets_read
            .saturating_sub(before.seed_buckets_read),
        sequence_blocks_read: after
            .sequence_blocks_read
            .saturating_sub(before.sequence_blocks_read),
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
        let mut work = pool.install(|| {
            search
                .candidates
                .par_iter()
                .map(|candidate| {
                    let candidate_started = Instant::now();
                    let entry = query.catalog.get(candidate.metagenome_id());
                    let mut work = process_candidate(
                        plasmid_id,
                        plasmid_sequence,
                        candidate,
                        entry,
                        run_config,
                    );
                    work.counters.failures = work.result.failures.len() as u64;
                    work.counters.elapsed_millis = candidate_started
                        .elapsed()
                        .as_millis()
                        .try_into()
                        .unwrap_or(u64::MAX);
                    work.result.performance_counters = work.counters.candidate_snapshot();
                    work
                })
                .collect::<Vec<_>>()
        });
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
            database_input: database.input,
            database_metrics: database.metrics,
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
    config: &TraceRunnerConfig,
) -> CandidateWork {
    let mut counters = TracePerformanceCounters {
        candidates_processed: 1,
        ..TracePerformanceCounters::default()
    };
    let Some(entry) = entry else {
        return CandidateWork {
            result: failed_result(
                plasmid_id,
                candidate,
                config,
                "catalog",
                "missing_catalog_entry",
                "candidate is absent from the resource catalog",
                false,
                None,
            ),
            counters,
        };
    };

    let mut failures = Vec::new();
    let indexed = entry.jma.as_deref();
    let mut result = if let Some(jma) = indexed {
        if entry.jma_index.is_none() && !config.resources.allow_full_download_fallback {
            failed_result(
                plasmid_id,
                candidate,
                config,
                "jma_index",
                "jma_index_required",
                "indexed JMA access requires a catalog jma_index; enable full-download fallback explicitly to use an archive without one",
                false,
                Some(jma),
            )
        } else {
            match process_jma(
                plasmid_id,
                plasmid,
                candidate,
                jma,
                entry.jma_index.as_deref(),
                config,
                &mut counters,
            ) {
                Ok(result) => result,
                Err(error)
                    if entry.fallback_raw().is_some()
                        && !matches!(
                            error.error.as_ref(),
                            RunnerError::SeedLevelMismatch { .. }
                        ) =>
                {
                    failures.push(failure_for_error("jma", jma, error.error.as_ref()));
                    let mut fallback = process_raw(
                        plasmid_id,
                        plasmid,
                        candidate,
                        entry.raw.as_deref().unwrap_or_default(),
                        config,
                        &mut counters,
                    );
                    fallback.resource_metrics =
                        add_resource_metrics(*error.metrics, fallback.resource_metrics);
                    fallback
                }
                Err(error) => failed_result_with_metrics(
                    plasmid_id,
                    candidate,
                    config,
                    "jma",
                    error.error.code(),
                    &error.error.to_string(),
                    error.error.retryable(),
                    Some(jma),
                    *error.metrics,
                ),
            }
        }
    } else if let Some(raw) = entry.raw.as_deref() {
        process_raw(plasmid_id, plasmid, candidate, raw, config, &mut counters)
    } else {
        failed_result(
            plasmid_id,
            candidate,
            config,
            "catalog",
            "missing_resource",
            "catalog row has neither a JMA nor raw assembly resource",
            false,
            None,
        )
    };
    if indexed.is_some()
        && entry.jma_index.is_none()
        && config.resources.allow_full_download_fallback
    {
        result.failures.push(TraceFailure {
            stage: "jma_index".to_string(),
            code: "jma_full_download_fallback".to_string(),
            message: "catalog row has no jma_index; archive was opened eagerly because full-download fallback is enabled".to_string(),
            resource: indexed.map(redact_resource),
            retryable: false,
        });
        result.resource_metrics.full_object_fallbacks = result
            .resource_metrics
            .full_object_fallbacks
            .saturating_add(1);
        if result.status == TraceStatus::Complete {
            result.status = TraceStatus::Partial;
        }
    }
    result.failures.splice(0..0, failures);
    if !result.failures.is_empty() && result.status == TraceStatus::Complete {
        result.status = TraceStatus::Partial;
    }
    CandidateWork { result, counters }
}

fn process_jma(
    plasmid_id: &str,
    plasmid: &[u8],
    candidate: &RankedCandidate,
    locator: &str,
    index_locator: Option<&str>,
    config: &TraceRunnerConfig,
    counters: &mut TracePerformanceCounters,
) -> Result<TraceMetagenomeResult, TraceProcessingFailure> {
    let resource = open_resource(locator, config.resources.clone())
        .map_err(|error| TraceProcessingFailure::new(error.into(), ResourceMetrics::default()))?;
    let reader = match index_locator {
        Some(index_locator) => {
            let index =
                open_resource(index_locator, config.resources.clone()).map_err(|error| {
                    TraceProcessingFailure::new(error.into(), ResourceMetrics::default())
                })?;
            JmaReader::open_indexed(resource, index).map_err(|error| {
                TraceProcessingFailure::new(error.into(), ResourceMetrics::default())
            })?
        }
        None => JmaReader::open(resource).map_err(|error| {
            TraceProcessingFailure::new(error.into(), ResourceMetrics::default())
        })?,
    };
    if reader.header().algorithm_id.as_deref() != Some(crate::trace::TRACE_ALGORITHM_ID)
        || reader.header().algorithm_version != Some(crate::trace::TRACE_ALGORITHM_VERSION)
    {
        return Err(TraceProcessingFailure::new(
            RunnerError::from(JmaError::CorruptSection(format!(
                "JMA algorithm metadata is missing or incompatible; expected {} version {}",
                crate::trace::TRACE_ALGORITHM_ID,
                crate::trace::TRACE_ALGORITHM_VERSION
            ))),
            reader.metrics(),
        ));
    }
    let mut alignments = Vec::new();
    let mut workspace = AlignmentWorkspace::new();
    let mut seen_contigs = HashSet::new();
    let mut seen_windows = HashSet::new();
    let mut seen_seed_keys = HashSet::new();
    let chain_model = chaining_coordinate_model(config.topology_requested);
    let mut round_configs = vec![config.sensitivity.primary];
    if let Some(dense_primary) = config.sensitivity.gap_rescue.dense_primary
        && dense_primary != config.sensitivity.primary
    {
        round_configs.push(dense_primary);
    }
    if let Some(rescue) = config.sensitivity.rescue {
        round_configs.push(rescue);
    }
    round_configs.truncate(usize::from(config.sensitivity.gap_rescue.max_rounds));

    let query_length = plasmid.len() as u64;
    let mut rescue_rounds = Vec::new();
    let mut target_gaps = vec![BaseInterval {
        start: 0,
        end: query_length,
    }];
    let mut supported_before = 0_u64;
    let mut terminal_seed_evidence = false;

    for (round_index, seed_config) in round_configs.into_iter().enumerate() {
        if round_index > 0 {
            target_gaps.retain(|gap| gap.len() >= config.sensitivity.gap_rescue.min_gap_bases);
            if matches!(config.topology_requested, TopologyRequested::Auto)
                && !terminal_seed_evidence
            {
                target_gaps.retain(|gap| gap.start != 0 && gap.end != query_length);
            }
            if target_gaps.is_empty() {
                break;
            }
        }
        let round_target_gaps = target_gaps.clone();

        let started = Instant::now();
        let before_metrics = reader.metrics();
        let round = chains_for_level(
            &reader,
            plasmid,
            seed_config,
            &config.sensitivity,
            chain_model,
            (round_index > 0).then_some(target_gaps.as_slice()),
            &mut seen_seed_keys,
        )
        .map_err(|error| TraceProcessingFailure::new(error, reader.metrics()))?;
        let LevelChainResult {
            chains,
            seed_keys_tested,
            anchors_created,
            chains_accepted,
            terminal_seed_evidence: round_terminal_seed_evidence,
        } = round;
        terminal_seed_evidence |= round_terminal_seed_evidence;
        let alignment_windows_attempted = align_indexed_chains(
            &reader,
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
            counters,
        )
        .map_err(|error| TraceProcessingFailure::new(error, reader.metrics()))?;

        let topology = assess_topology(
            query_length,
            config.topology_requested,
            config.topology_margin_bases,
            &alignments,
        )
        .map_err(|error| TraceProcessingFailure::new(error.into(), reader.metrics()))?;
        let selected = selected_mosaic(&topology);
        let supported_after = selected.base_covered_bases;
        target_gaps = selected
            .unsupported_gaps
            .iter()
            .map(|gap| gap.interval)
            .collect();
        let after_metrics = reader.metrics();
        let metric_delta = subtract_resource_metrics(after_metrics, before_metrics);
        rescue_rounds.push(RescueRoundMetrics {
            round: u8::try_from(round_index + 1).unwrap_or(u8::MAX),
            seed_k: seed_config.k,
            seed_scale: seed_config.scale,
            target_gaps: round_target_gaps,
            seed_buckets_requested: metric_delta.seed_buckets_read,
            seed_keys_tested,
            anchors_created,
            chains_accepted,
            sequence_blocks_fetched: metric_delta.sequence_blocks_read,
            alignment_windows_attempted,
            new_query_bases_supported: supported_after.saturating_sub(supported_before),
            elapsed_millis: started.elapsed().as_millis().try_into().unwrap_or(u64::MAX),
        });
        supported_before = supported_after;
        if supported_after == query_length {
            break;
        }
    }

    let finalized = finalize_evidence(query_length, config, &mut alignments)
        .map_err(|error| TraceProcessingFailure::new(error, reader.metrics()))?;
    let warnings = evidence_warnings(&finalized);
    Ok(TraceMetagenomeResult {
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
        alignments,
        primary_fragment_mosaic: Some(finalized.primary_mosaic),
        topology: Some(finalized.topology),
        rescue_rounds,
        performance_counters: CandidatePerformanceCounters::default(),
        coverage: Some(finalized.coverage),
        warnings,
        failures: Vec::new(),
        resource_metrics: reader.metrics(),
    })
}

fn process_raw(
    plasmid_id: &str,
    plasmid: &[u8],
    candidate: &RankedCandidate,
    locator: &str,
    config: &TraceRunnerConfig,
    counters: &mut TracePerformanceCounters,
) -> TraceMetagenomeResult {
    let resource = match RawAssembly::open_with_metrics(locator, config.resources.clone()) {
        Ok(resource) => resource,
        Err(failure) => {
            return failed_result_with_metrics(
                plasmid_id,
                candidate,
                config,
                "raw",
                raw_error_code(failure.error.as_ref()),
                &failure.error.to_string(),
                raw_error_retryable(failure.error.as_ref()),
                Some(locator),
                *failure.metrics,
            );
        }
    };
    let mut alignments = Vec::new();
    let mut workspace = AlignmentWorkspace::new();
    for contig in &resource.contigs {
        counters.contigs_considered = counters.contigs_considered.saturating_add(1);
        let windows = match window_ranges(
            contig.sequence.len(),
            plasmid.len(),
            config.sensitivity.max_alignment_window_bases,
        ) {
            Ok(windows) => windows,
            Err(error) => {
                return failed_result_with_metrics(
                    plasmid_id,
                    candidate,
                    config,
                    "raw",
                    error.code(),
                    &error.to_string(),
                    error.retryable(),
                    Some(locator),
                    resource.metrics,
                );
            }
        };
        for (start, end) in windows {
            counters.windows_retrieved = counters.windows_retrieved.saturating_add(1);
            let sequence = &contig.sequence[start..end];
            if let Err(error) = align_window(
                plasmid_id,
                candidate,
                &mut workspace,
                plasmid,
                sequence,
                start as u64,
                &contig.id,
                config,
                0,
                SeedEvidence::default(),
                None,
                0,
                plasmid.len() as u64,
                0,
                &mut alignments,
                counters,
            ) {
                return failed_result_with_metrics(
                    plasmid_id,
                    candidate,
                    config,
                    "alignment",
                    error.code(),
                    &error.to_string(),
                    error.retryable(),
                    Some(locator),
                    resource.metrics,
                );
            }
            retain_best_alignments(&mut alignments, config.max_alignments_per_candidate);
        }
    }
    let finalized = match finalize_evidence(plasmid.len() as u64, config, &mut alignments) {
        Ok(finalized) => finalized,
        Err(error) => {
            return failed_result_with_metrics(
                plasmid_id,
                candidate,
                config,
                "coverage",
                error.code(),
                &error.to_string(),
                false,
                Some(locator),
                resource.metrics,
            );
        }
    };
    let mut warnings = evidence_warnings(&finalized);
    warnings
        .push("raw assembly mode does not provide indexed gap-directed seed rescue".to_string());
    TraceMetagenomeResult {
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
        alignments,
        primary_fragment_mosaic: Some(finalized.primary_mosaic),
        topology: Some(finalized.topology),
        rescue_rounds: Vec::new(),
        performance_counters: CandidatePerformanceCounters::default(),
        coverage: Some(finalized.coverage),
        warnings,
        failures: Vec::new(),
        resource_metrics: resource.metrics,
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
    alignments: &mut Vec<BaseAlignment>,
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
    let max_cells = query_span_usize
        .saturating_add(1)
        .saturating_mul(band.saturating_mul(2).saturating_add(1));
    let options = AlignmentOptions::new(config.sensitivity.alignment)
        .with_diagonal_offset(diagonal_offset)
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
        match workspace.align_circular(
            plasmid,
            query_start,
            query_span,
            target,
            target_start,
            strand,
            options,
        ) {
            Ok(result) => {
                counters.alignments_succeeded = counters.alignments_succeeded.saturating_add(1);
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
    chain: AnchorChain,
    seed_evidence: SeedEvidence,
}

struct LevelChainResult {
    chains: Vec<EvidencedChain>,
    seed_keys_tested: u64,
    anchors_created: u64,
    chains_accepted: u64,
    terminal_seed_evidence: bool,
}

fn chains_for_level<R: RangeReader>(
    reader: &JmaReader<R>,
    plasmid: &[u8],
    seed_config: SeedSensitivity,
    sensitivity: &SensitivityConfig,
    coordinate_model: CoordinateModel,
    intervals: Option<&[BaseInterval]>,
    seen_seed_keys: &mut HashSet<(u8, u64, u64, u64)>,
) -> Result<LevelChainResult, RunnerError> {
    let available_scale = reader
        .header()
        .seed_levels
        .iter()
        .filter(|level| level.k == seed_config.k)
        .map(|level| level.scale)
        .min()
        .ok_or(RunnerError::SeedLevelMismatch {
            k: seed_config.k,
            required_scale: seed_config.scale,
            available_scale: None,
        })?;
    if available_scale > seed_config.scale {
        return Err(RunnerError::SeedLevelMismatch {
            k: seed_config.k,
            required_scale: seed_config.scale,
            available_scale: Some(available_scale),
        });
    }
    let mut level = match intervals {
        Some(intervals) => extract_seed_level_in_intervals(
            plasmid,
            seed_config,
            intervals,
            sensitivity.gap_rescue.flank_bases,
        )?,
        None => extract_seed_level(plasmid, seed_config)?,
    };
    level.seeds.retain(|seed| {
        seen_seed_keys.insert((level.k, seed.hash, seed.canonical_kmer, seed.position))
    });
    if intervals.is_some() {
        level.seeds.truncate(
            usize::try_from(sensitivity.gap_rescue.max_seed_buckets_per_round)
                .unwrap_or(usize::MAX),
        );
    }
    let seed_keys_tested = level.seeds.len() as u64;
    let mut groups = Vec::with_capacity(level.seeds.len());
    for seed in level.seeds {
        let occurrences = reader.seed_occurrences(seed.query(level.k))?;
        if !occurrences.is_empty() {
            groups.push(SeedOccurrenceGroup {
                seed,
                k: level.k,
                occurrences,
            });
        }
    }
    let anchors = generate_anchors(
        &groups,
        seed_config.max_occurrences,
        sensitivity.max_anchors_per_candidate,
    );
    let anchor_evidence = summarize_anchor_evidence(
        &anchors.anchors,
        &groups,
        seed_config.max_occurrences,
        sensitivity.common_seed_candidate_occurrence_threshold,
    );
    let terminal_seed_evidence = has_terminal_seed_pair(
        &anchor_evidence,
        plasmid.len() as u64,
        sensitivity
            .gap_rescue
            .flank_bases
            .max(u64::from(seed_config.k)),
    );
    let chains = chain_anchors(
        &anchors.anchors,
        plasmid.len() as u64,
        ChainConfig::from_sensitivity(sensitivity).with_coordinate_model(coordinate_model),
    )
    .map_err(RunnerError::from)?;
    let chains_accepted = chains.len() as u64;
    let chains = chains
        .into_iter()
        .map(|chain| EvidencedChain {
            seed_evidence: evidence_for_chain(&chain, &anchor_evidence),
            chain,
        })
        .collect();
    Ok(LevelChainResult {
        chains,
        seed_keys_tested,
        anchors_created: anchors.anchors.len() as u64,
        chains_accepted,
        terminal_seed_evidence,
    })
}

fn evidence_for_chain(chain: &AnchorChain, evidence: &[AnchorOccurrenceEvidence]) -> SeedEvidence {
    let mut result = SeedEvidence::default();
    for anchor in &chain.anchors {
        if anchor.k == crate::trace::seeds::PRIMARY_K {
            result.primary_anchor_count = result.primary_anchor_count.saturating_add(1);
        } else if anchor.k == crate::trace::seeds::RESCUE_K {
            result.rescue_anchor_count = result.rescue_anchor_count.saturating_add(1);
        }
        if let Some(item) = evidence.iter().find(|item| item.anchor == *anchor) {
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

fn has_terminal_seed_pair(
    evidence: &[AnchorOccurrenceEvidence],
    query_length: u64,
    terminal_width: u64,
) -> bool {
    let mut by_target = BTreeMap::<(u32, u8), (bool, bool)>::new();
    for item in evidence {
        if item.seed.is_repetitive {
            continue;
        }
        let strand = match item.anchor.strand {
            Strand::Forward => 0,
            Strand::Reverse => 1,
        };
        let flags = by_target
            .entry((item.anchor.contig_id, strand))
            .or_default();
        flags.0 |= item.anchor.query_position <= terminal_width;
        flags.1 |= item
            .anchor
            .query_position
            .saturating_add(u64::from(item.anchor.k))
            >= query_length.saturating_sub(terminal_width);
    }
    by_target.values().any(|(start, end)| *start && *end)
}

#[allow(clippy::too_many_arguments)]
fn align_indexed_chains<R: RangeReader>(
    reader: &JmaReader<R>,
    query_id: &str,
    query: &[u8],
    candidate: &RankedCandidate,
    config: &TraceRunnerConfig,
    chains: Vec<EvidencedChain>,
    max_windows: u64,
    max_sequence_blocks: u64,
    workspace: &mut AlignmentWorkspace,
    seen_contigs: &mut HashSet<u32>,
    seen_windows: &mut HashSet<(u32, u64, u64, u8)>,
    alignments: &mut Vec<BaseAlignment>,
    counters: &mut TracePerformanceCounters,
) -> Result<u64, RunnerError> {
    let mut attempted = 0_u64;
    let blocks_before = reader.metrics().sequence_blocks_read;
    for evidenced in chains {
        if attempted >= max_windows
            || reader
                .metrics()
                .sequence_blocks_read
                .saturating_sub(blocks_before)
                >= max_sequence_blocks
        {
            break;
        }
        let chain = evidenced.chain;
        let contig = reader
            .contigs()
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
        if !seen_windows.insert((contig.id, range.start, range.end, strand_rank)) {
            continue;
        }
        let sequence = reader.read_sequence(contig.id, range)?;
        counters.windows_retrieved = counters.windows_retrieved.saturating_add(1);
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
            alignments,
            counters,
        )?;
        retain_best_alignments(alignments, config.max_alignments_per_candidate);
    }
    Ok(attempted)
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

fn window_ranges(
    contig_length: usize,
    query_length: usize,
    max_window_bases: u64,
) -> Result<Vec<(usize, usize)>, RunnerError> {
    if contig_length == 0 {
        return Ok(Vec::new());
    }
    let max_window = usize::try_from(max_window_bases)
        .map_err(|_| RunnerError::InvalidConfig("alignment window exceeds platform".to_string()))?;
    if max_window == 0 {
        return Err(RunnerError::InvalidConfig(
            "max_alignment_window_bases must be greater than zero".to_string(),
        ));
    }
    if query_length > max_window {
        return Err(RunnerError::WindowTooLarge {
            query_length,
            max_window,
        });
    }
    let span = contig_length.min(max_window);
    let overlap = query_length.min(span.saturating_sub(1)).max(1);
    let step = span.saturating_sub(overlap).max(1);
    let mut ranges = Vec::new();
    let mut start = 0usize;
    loop {
        let end = start.saturating_add(span).min(contig_length);
        ranges.push((start, end));
        if end == contig_length {
            break;
        }
        let next = start.saturating_add(step);
        if next <= start {
            break;
        }
        start = next.min(contig_length - 1);
    }
    Ok(ranges)
}

fn retain_best_alignments(alignments: &mut Vec<BaseAlignment>, limit: usize) {
    alignments.sort_by(compare_alignments);
    alignments.dedup_by(|right, left| {
        right.contig_id == left.contig_id
            && right.strand == left.strand
            && right.target_interval == left.target_interval
            && right.query_segments == left.query_segments
    });
    alignments.truncate(limit);
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
        alignments: Vec::new(),
        primary_fragment_mosaic: None,
        topology: None,
        rescue_rounds: Vec::new(),
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
        resource_metrics,
    }
}

fn failure_for_error(stage: &str, resource: &str, error: &RunnerError) -> TraceFailure {
    TraceFailure {
        stage: stage.to_string(),
        code: error.code().to_string(),
        message: error.to_string(),
        resource: Some(redact_resource(resource)),
        retryable: error.retryable(),
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

fn raw_error_retryable(error: &RawError) -> bool {
    match error {
        RawError::Resource(error) => resource_error_retryable(error),
        RawError::Parse { .. } | RawError::Io { .. } => false,
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
}
