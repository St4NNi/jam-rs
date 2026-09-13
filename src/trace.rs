use crate::alignment::{
    AlignmentConfig, AlignmentError, Interval, Strand, TraceAlignmentWorkspace,
};
use crate::bgzf::{BgzfError, BgzfReader};
use crate::bgzf_cache::{BgzfBlockCache, DEFAULT_BATCH_BGZF_CACHE_BYTES};
use crate::jidx::{JidxError, RESCUE_K15_TAG, sha256, sha256_reader};
use crate::jidx_reader::{
    ContigId, JidxReader, JidxReaderError, MetagenomeId, SEED_LOOKUP_BATCH_KEYS,
};
use crate::mosaic::{Fragment, Mosaic, MosaicError, build_mosaic};
use crate::query::{QueryEngine, QueryError, QuerySketch};
use crate::range_source::S3Config;
use crate::reader::ReaderError;
use crate::trace_batch::{
    SharedCoreLookups, SharedSeedLookups, TraceBatch, lookup_budget, lookup_bytes, phase_elapsed,
    phase_stamp, prepare_cores, prepare_lookup_with_cores, prepare_screened_cores,
};
use crate::trace_index::{
    TraceCacheIdentity, TraceDocument as SeedDocument, TraceIndex, TraceSeed,
};
use needletail::Sequence;
use rayon::prelude::*;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fmt::Write as _;
use std::fs::File;
use std::io::{self, BufReader};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;
use thiserror::Error;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TraceConfig {
    pub min_containment: f64,
    pub max_metagenomes: usize,
    pub use_sketch: bool,
    pub min_seed_hits: u32,
    pub diagonal_bin_bases: u64,
    pub flank_bases: u64,
    pub min_identity: f64,
    pub min_aligned_bases: u64,
    pub endpoint_bases: usize,
    pub circular: bool,
    pub verify_resources: bool,
    pub alignment: AlignmentConfig,
}

impl Default for TraceConfig {
    fn default() -> Self {
        Self {
            min_containment: 0.01,
            max_metagenomes: usize::MAX,
            use_sketch: true,
            min_seed_hits: 2,
            diagonal_bin_bases: 64,
            flank_bases: 256,
            min_identity: 0.8,
            min_aligned_bases: 40,
            endpoint_bases: 256,
            circular: true,
            verify_resources: false,
            alignment: AlignmentConfig {
                max_cells: 1 << 24,
                ..AlignmentConfig::default()
            },
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct TraceResult {
    pub query_id: String,
    pub query_length: u64,
    pub index: TraceIndexIdentity,
    pub completion: SearchCompletion,
    pub candidates_screened: u32,
    pub metagenomes: Vec<MetagenomeTrace>,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct TraceIndexIdentity {
    pub manifest_sha256: String,
    pub body_sha256: String,
    pub seed_k: u8,
    pub rescue_k15: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(tag = "status", rename_all = "snake_case")]
pub enum SearchCompletion {
    Complete,
    CandidateBudgetExceeded { candidates_omitted: u32 },
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct MetagenomeTrace {
    pub metagenome_id: MetagenomeId,
    pub name: String,
    pub shared_hashes: u32,
    pub containment: f64,
    pub exact_seed_hits: u64,
    pub compressed_bytes_read: u64,
    pub range_requests: u64,
    pub bgzf_blocks_decoded: u64,
    pub contigs: Vec<TraceContig>,
    pub mosaic: Mosaic,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct TraceContig {
    pub id: ContigId,
    pub name: String,
}

pub struct TraceEngine {
    jam_path: Option<PathBuf>,
    screen: Option<QueryEngine>,
    index: TraceIndex,
    sample_to_metagenome: Vec<MetagenomeId>,
    s3: Option<S3Config>,
    batch_stats: Mutex<TraceBatchStats>,
    pub(crate) observed: bool,
    pub(crate) region_policy: crate::trace_islands::RegionPolicy,
    /// Directory for per-query region support ledgers (diagnostic runs only).
    pub(crate) region_ledger: Option<PathBuf>,
}

#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct TraceBatchStats {
    pub timings_observed: bool,
    // Candidate ordering, geometry, windows, inclusive circular retries, result reduction.
    // Retry CPU overlaps the alignment kernel samples; it is never added to them.
    pub downstream_thread_cpu_ns: [u64; 5],
    pub alignment_work: crate::alignment::AlignmentWork,
    pub alignment_tasks: u64,
    pub identical_alignment_tasks: u64,
    pub task_signature_unavailable_queries: u64,
    pub alignment_workspaces: u64,
    pub alignment_query_bytes: u64,
    pub alignment_target_bytes: u64,
    pub alignment_query_max: u64,
    pub alignment_target_max: u64,
    pub circular_retries: u64,
    pub copied_query_bytes: u64,
    pub returned_alignments: u64,
    pub rejected_alignments: u64,
    pub region_predecessor_tests: u64,
    pub geometric_hits: u64,
    pub fragments_before_dedup: u64,
    pub fragments_after_dedup: u64,
    pub seed_generation_ns: u64,
    pub key_lookup_ns: u64,
    pub membership_access_ns: u64,
    pub position_access_ns: u64,
    pub posting_member_tasks: u64,
    pub posting_plan_ns: u64,
    pub posting_position_tasks: u64,
    pub posting_plan_hash: u64,
    pub posting_admitted_member_rows: u64,
    pub posting_admitted_position_rows: u64,
    pub posting_peak_parallel_tasks: u64,
    pub posting_scratch_bytes: u64,
    pub posting_member_copies: u64,
    pub posting_initialized_position_bytes: u64,
    pub posting_member_worker_elapsed_ns: u64,
    pub posting_position_worker_elapsed_ns: u64,
    pub posting_member_worker_cpu_ns: u64,
    pub posting_position_worker_cpu_ns: u64,
    pub posting_cpu_unavailable_tasks: u64,
    pub posting_unused_member_lists: u64,
    pub posting_unused_position_rows: u64,
    pub posting_usage_unavailable_batches: u64,
    pub candidate_routing_ns: u64,
    pub region_formation_ns: u64,
    pub sequence_read_ns: u64,
    pub bgzf_decode_and_handling_ns: u64,
    pub alignment_ns: u64,
    pub traceback_and_cigar_ns: u64,
    pub endpoint_completion_ns: u64,
    pub search_critical_ns: u64,
    pub batches: u64,
    pub unique_keys: u64,
    pub distinct_cores: u64,
    pub split_core_resolutions: u64,
    pub phase_extraction_ns: [u64; 2],
    pub phase_core_lookup_ns: [u64; 2],
    pub phase_context_generation_ns: [u64; 2],
    pub phase_context_lookup_ns: [u64; 2],
    pub phase_postings_ns: [u64; 2],
    pub phase_downstream_ns: [u64; 2],
    pub query_core_occurrences: u64,
    pub extraction_occurrence_probes: u64,
    pub extraction_covered_probes: u64,
    pub extraction_uncovered_probes: u64,
    pub tokens_rejected_before_sort: u64,
    pub surviving_tokens_sorted: u64,
    pub surviving_distinct_query_cores: u64,
    pub query_context_associations: u64,
    pub query_distinct_context_requests: u64,
    pub query_distinct_cores: u64,
    pub query_executed_context_associations: u64,
    pub nested_context_calls: u64,
    pub core_lookup_tasks: u64,
    pub core_requests_attempted: u64,
    pub core_requests_covered: u64,
    pub core_requests_rejected: u64,
    pub core_requests_uncovered: u64,
    pub core_requests_prescreened: u64,
    pub core_requests_retained: u64,
    pub core_requests_planned: u64,
    pub core_lookup_fallbacks: u64,
    pub core_lookup_ns: u64,
    pub core_lookup_peak_bytes: usize,
    pub core_lookup_retained_bytes: usize,
    pub numeric_metadata_capacity_bound: usize,
    pub physical_geometry_capacity_bound: usize,
    pub geometry_phase_reserved_bytes: usize,
    pub executed_anchor_associations: u64,
    pub query_sequence_capacity_bytes: u64,
    pub query_core_capacity_bytes: u64,
    pub query_token_capacity_bytes: u64,
    pub query_token_reserved_peak_bytes: u64,
    pub query_core_conversion_capacity_bound: u64,
    pub query_compact_core_queries: u64,
    pub query_wide_core_queries: u64,
    pub query_nested_capacity_bytes: u64,
    pub query_directory_capacity_bytes: u64,
    pub emitted_anchor_associations: u64,
    pub restored_query_context_requests: u64,
    pub context_reuse_histogram_log2: [u64; 16],
    pub context_occurrence_histogram_log2: [u64; 16],
    pub lookup_tasks: u64,
    pub lookup_plan_hash: u64,
    pub lookup_dispatch_ns: u64,
    pub lookup_parallel_ns: u64,
    pub lookup_compute_ns: u64,
    pub lookup_reduce_ns: u64,
    pub lookup_dispatch_to_start_ns: u64,
    pub cached_groups: u64,
    pub cached_positions: u64,
    pub lookup_peak_bytes: usize,
    pub lookup_retained_bytes: usize,
    pub lookup_reserved_bytes: usize,
    pub bgzf_cache_hits: u64,
    pub bgzf_blocks_decoded: u64,
    pub bgzf_evictions: u64,
    pub bgzf_peak_bytes: usize,
    pub bgzf_io_limit: usize,
    pub bgzf_io_peak: usize,
}

pub(crate) struct PreparedQuery {
    pub(crate) query_id: String,
    pub(crate) query_length: u64,
    query: Vec<u8>,
    positions_by_key: QueryPositions,
    lookup_identity: [u8; 32],
    batch_ordinal: usize,
    shared_cores: Option<Arc<SharedCoreLookups>>,
}

impl PreparedQuery {
    fn find_seeds(
        &self,
        index: &TraceIndex,
        keys: &[u64],
    ) -> Result<Vec<Option<TraceSeed>>, TraceError> {
        if let Some(cores) = &self.shared_cores {
            index.find_seeds_in_cores(keys, cores)
        } else {
            index.find_seeds_batch(keys)
        }
    }
}

pub(crate) struct TraceCensus {
    pub(crate) candidates: Vec<Candidate>,
    pub(crate) frequencies: Vec<(u64, u32)>,
    pub(crate) lookups: Option<CachedSeedLookups>,
}

pub(crate) const LOOKUP_CACHE_BYTES: usize = 256 * 1024 * 1024;
pub(crate) static LOOKUP_CACHE_AVAILABLE: AtomicUsize = AtomicUsize::new(LOOKUP_CACHE_BYTES);

pub(crate) struct CacheReservation<'a> {
    available: &'a AtomicUsize,
    pub(crate) bytes: usize,
}

impl<'a> CacheReservation<'a> {
    pub(crate) fn acquire(available: &'a AtomicUsize, bytes: usize) -> Option<Self> {
        available
            .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |remaining| {
                remaining.checked_sub(bytes)
            })
            .ok()
            .map(|_| Self { available, bytes })
    }

    pub(crate) fn retain(&mut self, bytes: usize) {
        assert!(bytes <= self.bytes);
        self.available
            .fetch_add(self.bytes - bytes, Ordering::Relaxed);
        self.bytes = bytes;
    }
}

impl Drop for CacheReservation<'_> {
    fn drop(&mut self) {
        self.available.fetch_add(self.bytes, Ordering::Relaxed);
    }
}

pub(crate) struct CachedSeedLookups {
    header_sha256: [u8; 32],
    query_identity: [u8; 32],
    file_identity: TraceCacheIdentity,
    through: Option<u64>,
    frozen: bool,
    groups: Vec<CachedDocumentGroup>,
    documents: Vec<SeedDocument>,
    _reservation: CacheReservation<'static>,
}

#[derive(Clone, Copy)]
struct CachedDocumentGroup {
    seed: TraceSeed,
    document_start: usize,
    document_count: usize,
}

#[derive(Clone, Copy)]
struct CachedSeedLookup<'a> {
    seed: TraceSeed,
    documents: &'a [SeedDocument],
}

impl CachedSeedLookups {
    fn new(
        header_sha256: [u8; 32],
        query_identity: [u8; 32],
        file_identity: TraceCacheIdentity,
        bytes: usize,
        query_keys: usize,
        document_count: u32,
    ) -> Option<Self> {
        let fixed = std::mem::size_of::<Option<Self>>().checked_add(4096)?;
        let documents_per_group = usize::try_from(document_count).ok()?;
        let bytes_per_group = std::mem::size_of::<CachedDocumentGroup>()
            .checked_add(documents_per_group.checked_mul(std::mem::size_of::<SeedDocument>())?)?;
        let group_limit = bytes
            .checked_sub(fixed)?
            .checked_div(bytes_per_group)?
            .min(query_keys);
        if group_limit == 0 {
            return None;
        }
        let document_limit = group_limit.checked_mul(documents_per_group)?;
        let reservation = CacheReservation::acquire(&LOOKUP_CACHE_AVAILABLE, bytes)?;
        let mut groups = Vec::new();
        let mut documents = Vec::new();
        groups.try_reserve_exact(group_limit).ok()?;
        documents.try_reserve_exact(document_limit).ok()?;
        let allocated = fixed
            .checked_add(
                groups
                    .capacity()
                    .checked_mul(std::mem::size_of::<CachedDocumentGroup>())?,
            )?
            .checked_add(
                documents
                    .capacity()
                    .checked_mul(std::mem::size_of::<SeedDocument>())?,
            )?;
        if allocated > reservation.bytes {
            return None;
        }
        Some(Self {
            header_sha256,
            query_identity,
            file_identity,
            through: None,
            frozen: false,
            groups,
            documents,
            _reservation: reservation,
        })
    }

    fn get(&self, key: u64) -> Option<Option<CachedSeedLookup<'_>>> {
        self.through.is_some_and(|through| key <= through).then(|| {
            self.groups
                .binary_search_by_key(&key, |group| group.seed.packed_key())
                .ok()
                .map(|index| {
                    let group = self.groups[index];
                    CachedSeedLookup {
                        seed: group.seed,
                        documents: &self.documents
                            [group.document_start..group.document_start + group.document_count],
                    }
                })
        })
    }

    fn cache_negative(&mut self, key: u64) {
        if !self.frozen {
            self.through = Some(key);
        }
    }

    fn cache_group(&mut self, key: u64, seed: TraceSeed, documents: &[SeedDocument]) {
        if self.frozen {
            return;
        }
        let Some(document_end) = self.documents.len().checked_add(documents.len()) else {
            self.frozen = true;
            return;
        };
        if self.groups.len() == self.groups.capacity() || document_end > self.documents.capacity() {
            self.frozen = true;
            return;
        }
        let document_start = self.documents.len();
        self.documents.extend_from_slice(documents);
        self.groups.push(CachedDocumentGroup {
            seed,
            document_start,
            document_count: documents.len(),
        });
        self.through = Some(key);
    }
}

/// Internal region-policy study switch; release builds always use parent envelopes.
fn study_region_policy() -> Result<crate::trace_islands::RegionPolicy, TraceError> {
    #[cfg(feature = "bench-internals")]
    match std::env::var("JAM_REGION_POLICY").as_deref() {
        Err(std::env::VarError::NotPresent) | Ok("parent") => {}
        Ok("islands") => return Ok(crate::trace_islands::RegionPolicy::Islands),
        _ => return Err(TraceError::Invalid("JAM_REGION_POLICY")),
    }
    Ok(crate::trace_islands::RegionPolicy::Parent)
}

fn study_region_ledger() -> Option<PathBuf> {
    #[cfg(feature = "bench-internals")]
    if let Some(path) = std::env::var_os("JAM_REGION_LEDGER") {
        return Some(path.into());
    }
    None
}

impl TraceEngine {
    pub fn open(
        jam: impl AsRef<Path>,
        jidx: impl AsRef<Path>,
        manifest: impl AsRef<Path>,
        s3: Option<S3Config>,
    ) -> Result<Self, TraceError> {
        let jam = jam.as_ref();
        let index = JidxReader::open(jidx)?;
        let manifest_sha256 = sha256_reader(BufReader::new(File::open(manifest)?))?;
        if manifest_sha256 != index.header().manifest_sha256 {
            return Err(TraceError::Invalid(
                "JIDX belongs to a different root manifest",
            ));
        }
        let screen = QueryEngine::open(jam)?;
        let mut by_name = HashMap::new();
        for id in 0..index.header().document_count {
            let name = index
                .metagenome_name(id)?
                .ok_or(TraceError::Invalid("missing JIDX metagenome"))?;
            if by_name.insert(name.to_string(), id).is_some() {
                return Err(TraceError::Invalid("duplicate JIDX metagenome"));
            }
        }
        let mut sample_to_metagenome = Vec::with_capacity(screen.reader().sample_names().len());
        for name in screen.reader().sample_names() {
            sample_to_metagenome.push(
                *by_name
                    .get(name)
                    .ok_or(TraceError::Invalid("JAM and JIDX names differ"))?,
            );
        }
        if sample_to_metagenome.len() != by_name.len() {
            return Err(TraceError::Invalid("JAM and JIDX names differ"));
        }
        Ok(Self {
            jam_path: Some(jam.to_path_buf()),
            screen: Some(screen),
            index: TraceIndex::Shard(Box::new(index)),
            sample_to_metagenome,
            s3,
            batch_stats: Mutex::default(),
            observed: false,
            region_policy: Default::default(),
            region_ledger: None,
        })
    }

    pub fn open_owner(root: impl AsRef<Path>, s3: Option<S3Config>) -> Result<Self, TraceError> {
        let index = crate::owner_reader::OwnerReader::open(root)?;
        if !index.is_complete() {
            return Err(TraceError::Invalid("owner generation is incomplete"));
        }
        Ok(Self {
            jam_path: None,
            screen: None,
            index: TraceIndex::Owner(index),
            sample_to_metagenome: Vec::new(),
            s3,
            batch_stats: Mutex::default(),
            observed: false,
            region_policy: Default::default(),
            region_ledger: None,
        })
    }

    pub fn open_shared(path: impl AsRef<Path>, s3: Option<S3Config>) -> Result<Self, TraceError> {
        Self::open_shared_observed(path, s3, false)
    }

    pub(crate) fn open_shared_observed(
        path: impl AsRef<Path>,
        s3: Option<S3Config>,
        observed: bool,
    ) -> Result<Self, TraceError> {
        let index = if observed {
            crate::shared_reader::SharedReader::open_observed(path)?
        } else {
            crate::shared_reader::SharedReader::open(path)?
        };
        Ok(Self {
            jam_path: None,
            screen: None,
            index: TraceIndex::Shared(Box::new(index)),
            sample_to_metagenome: Vec::new(),
            s3,
            batch_stats: Mutex::default(),
            observed,
            region_policy: study_region_policy()?,
            region_ledger: study_region_ledger(),
        })
    }

    pub fn verify_index(&self) -> Result<(), TraceError> {
        if let Some(jam_path) = &self.jam_path
            && sha256_reader(BufReader::new(File::open(jam_path)?))?
                != self.index.shard()?.header().jam_sha256
        {
            return Err(TraceError::Invalid(
                "JIDX belongs to a different JAM database",
            ));
        }
        self.index.verify_checksum()?;
        Ok(())
    }

    pub(crate) fn index(&self) -> Result<&JidxReader, TraceError> {
        self.index.shard()
    }

    pub fn search(
        &self,
        query_id: impl Into<String>,
        sequence: &[u8],
        config: TraceConfig,
    ) -> Result<TraceResult, TraceError> {
        let prepared = self.prepare(query_id, sequence, config)?;
        if self.index.is_shared() {
            let batch = self.prepare_batch(std::slice::from_ref(&prepared))?;
            let result = self.search_prepared(prepared, config, Some(&batch));
            self.record_batch(&batch);
            result
        } else {
            self.search_prepared(prepared, config, None)
        }
    }

    #[cfg(test)]
    pub(crate) fn search_without_batch(
        &self,
        query_id: &str,
        sequence: &[u8],
        config: TraceConfig,
    ) -> Result<TraceResult, TraceError> {
        let prepared = self.prepare(query_id, sequence, config)?;
        self.search_prepared(prepared, config, None)
    }

    #[cfg(feature = "bench-internals")]
    pub fn benchmark_core_planning(
        &self,
        keys: &[u32],
        early: bool,
    ) -> Result<(impl Sized, [usize; 7]), TraceError> {
        let result = if early {
            prepare_cores(&self.index, keys.iter().copied(), keys.len(), false)?
        } else {
            crate::trace_batch::prepare_cores_late(
                &self.index,
                keys.iter().copied(),
                keys.len(),
                false,
            )?
        };
        let counts = result.as_ref().map_or([0; 7], |r| {
            [
                r.requests.attempted,
                r.requests.covered,
                r.requests.rejected,
                r.requests.uncovered,
                r.requests.retained,
                r.requests.planned,
                r.tasks,
            ]
        });
        Ok((result, counts))
    }

    #[cfg(feature = "bench-internals")]
    pub fn benchmark_prepare(
        &self,
        sequence: &[u8],
        config: TraceConfig,
    ) -> Result<impl Sized, TraceError> {
        self.prepare("benchmark", sequence, config)
    }

    #[cfg(feature = "bench-internals")]
    pub fn benchmark_prepare_directory(
        &self,
        sequence: &[u8],
        config: TraceConfig,
    ) -> Result<impl Sized, TraceError> {
        self.prepare_directory(sequence, config, true)
    }

    #[cfg(any(test, feature = "bench-internals"))]
    fn prepare_directory(
        &self,
        sequence: &[u8],
        config: TraceConfig,
        compact: bool,
    ) -> Result<PreparedQuery, TraceError> {
        let scratch = compact
            .then(|| reserve_query_tokens([sequence.len()], lookup_budget(&self.index)))
            .flatten();
        let mut prepared =
            prepare_query_with_tokens("benchmark", sequence, config, 15, false, scratch.is_some())?;
        self.prepare_shared_queries(
            std::slice::from_mut(&mut prepared),
            &[config.circular],
            scratch,
            None,
        )?;
        Ok(prepared)
    }

    #[cfg(test)]
    pub(crate) fn assert_directory_associations(&self, sequence: &[u8], config: TraceConfig) {
        let actual = self.prepare("benchmark", sequence, config).unwrap();
        if let Some(operation) = self.extraction_operation(true).unwrap() {
            let mut fallback = prepare_query_screened(
                "benchmark",
                sequence,
                config,
                15,
                false,
                false,
                Some(&operation),
                true,
            )
            .unwrap();
            self.prepare_shared_queries(
                std::slice::from_mut(&mut fallback),
                &[config.circular],
                None,
                Some(&operation),
            )
            .unwrap();
            operation.finish().unwrap();
            assert_eq!(
                actual.positions_by_key.iter().collect::<Vec<_>>(),
                fallback.positions_by_key.iter().collect::<Vec<_>>()
            );
        }
        for compact in [true, false] {
            let expected = self.prepare_directory(sequence, config, compact).unwrap();
            assert_eq!(
                actual.positions_by_key.core_occurrences,
                expected.positions_by_key.core_occurrences
            );
            assert_eq!(
                actual.positions_by_key.core_distinct,
                expected.positions_by_key.core_distinct
            );
            assert_eq!(
                actual.positions_by_key.iter().collect::<Vec<_>>(),
                expected.positions_by_key.iter().collect::<Vec<_>>()
            );
        }
    }

    #[cfg(feature = "bench-internals")]
    pub fn benchmark_posting_preparation(
        &self,
        keys: &[u64],
        parallel: bool,
        observed: bool,
    ) -> Result<(impl Sized, [u64; 12]), TraceError> {
        let TraceIndex::Shared(reader) = &self.index else {
            return Err(TraceError::Invalid("shared posting benchmark"));
        };
        let seeds = self.index.find_seeds_batch(keys)?;
        let mut entries = keys
            .iter()
            .copied()
            .zip(seeds)
            .filter(|(_, seed)| seed.is_some())
            .collect::<Vec<_>>();
        entries.sort_unstable_by_key(|entry| entry.0);
        entries.dedup_by_key(|entry| entry.0);
        let mut postings = (0..entries.len()).map(|_| None).collect::<Vec<_>>();
        let stats = crate::trace_postings::benchmark_postings(
            reader,
            &entries,
            &mut postings,
            crate::trace_batch::lookup_budget(&self.index),
            observed,
            parallel,
        )?
        .ok_or(TraceError::Invalid("posting benchmark admission"))?;
        Ok((
            postings,
            [
                stats.member_tasks as u64,
                stats.position_tasks as u64,
                stats.membership_ns,
                stats.position_ns,
                stats.member_worker_cpu_ns,
                stats.position_worker_cpu_ns,
                stats.admitted_member_rows,
                stats.admitted_position_rows,
                stats.scratch_bytes as u64,
                stats.retained_bytes as u64,
                stats.peak_parallel_tasks as u64,
                stats.task_hash,
            ],
        ))
    }

    #[cfg(feature = "bench-internals")]
    pub fn benchmark_lookup(
        &self,
        keys: &[u64],
        observed: bool,
    ) -> Result<(impl Sized, [u64; 3]), TraceError> {
        let requests = keys.iter().map(|&key| (key, 0)).collect();
        let lookups = prepare_lookup_with_cores(&self.index, requests, 1, observed, None)?;
        let stats = lookups.as_ref().map_or([0; 3], |lookup| {
            [
                lookup.lookup_tasks as u64,
                lookup.lookup_plan_hash,
                lookup.entries.len() as u64,
            ]
        });
        Ok((lookups, stats))
    }

    fn prepare(
        &self,
        query_id: impl Into<String>,
        sequence: &[u8],
        config: TraceConfig,
    ) -> Result<PreparedQuery, TraceError> {
        let started = self.observed.then(Instant::now);
        let scratch = self
            .index
            .is_shared()
            .then(|| reserve_query_tokens([sequence.len()], lookup_budget(&self.index)))
            .flatten();
        let operation = self.extraction_operation(scratch.is_some())?;
        let mut prepared = prepare_query_screened(
            query_id,
            sequence,
            config,
            self.index.k(),
            self.index.rescue_k15(),
            scratch.is_some(),
            operation.as_ref(),
            self.observed,
        )?;
        if self.index.is_shared() {
            self.prepare_shared_queries(
                std::slice::from_mut(&mut prepared),
                &[config.circular],
                scratch,
                operation.as_ref(),
            )?;
            if let Some(operation) = operation {
                operation.finish()?;
            }
        } else {
            self.record_prepared(&prepared, config.circular);
        }
        if let Some(started) = started {
            self.batch_stats.lock().unwrap().seed_generation_ns +=
                started.elapsed().as_nanos() as u64;
        }
        Ok(prepared)
    }

    fn extraction_operation(
        &self,
        admitted: bool,
    ) -> Result<Option<crate::shared_reader::SharedCoreOperation<'_>>, TraceError> {
        match &self.index {
            TraceIndex::Shared(reader) if admitted => Ok(reader.core_operation()?),
            _ => Ok(None),
        }
    }

    fn prepare_shared_queries(
        &self,
        prepared: &mut [PreparedQuery],
        circular: &[bool],
        mut scratch: Option<CacheReservation<'static>>,
        operation: Option<&crate::shared_reader::SharedCoreOperation<'_>>,
    ) -> Result<(), TraceError> {
        let core_phase = phase_stamp();
        if let Some(operation) = operation {
            for query in prepared.iter() {
                operation.record(query.positions_by_key.occurrence_probes);
            }
        }
        if self.observed {
            let mut stats = self.batch_stats.lock().unwrap();
            stats.query_token_reserved_peak_bytes = stats.query_token_reserved_peak_bytes.max(
                scratch
                    .as_ref()
                    .map_or(0, |reservation| reservation.bytes as u64),
            );
        }
        if let Some(scratch) = &mut scratch {
            let bytes = prepared
                .iter()
                .try_fold(0usize, |sum, query| {
                    let bytes = query.positions_by_key.tokens.as_ref().map_or(0, |tokens| {
                        tokens.capacity() * std::mem::size_of::<u64>()
                            + tokens.len() * std::mem::size_of::<QuerySeed>()
                    });
                    sum.checked_add(bytes)
                })
                .ok_or(TraceError::Invalid("query token reservation"))?;
            scratch.retain(bytes);
        }
        if self.observed {
            let mut stats = self.batch_stats.lock().unwrap();
            for query in prepared.iter() {
                let positions = &query.positions_by_key;
                let probes = positions.occurrence_probes;
                stats.extraction_occurrence_probes += probes.attempted as u64;
                stats.extraction_covered_probes += probes.covered as u64;
                stats.extraction_uncovered_probes += probes.uncovered as u64;
                stats.tokens_rejected_before_sort += probes.rejected as u64;
                stats.surviving_tokens_sorted +=
                    positions.tokens.as_ref().map_or(0, Vec::len) as u64;
                stats.surviving_distinct_query_cores += positions.directory.len() as u64;
            }
        }
        let count = prepared
            .iter()
            .try_fold(0usize, |sum, query| {
                sum.checked_add(query.positions_by_key.len())
            })
            .ok_or(TraceError::Invalid("core request count"))?;
        let keys = prepared
            .iter()
            .flat_map(|query| query.positions_by_key.keys().map(|&key| key as u32));
        let cores = match operation {
            Some(operation) => {
                prepare_screened_cores(&self.index, keys, count, self.observed, operation)?
            }
            None => prepare_cores(&self.index, keys, count, self.observed)?,
        };
        if self.observed && cores.is_none() {
            self.batch_stats.lock().unwrap().core_lookup_fallbacks += 1;
        }
        if self.observed
            && let Some(cores) = &cores
        {
            let mut stats = self.batch_stats.lock().unwrap();
            stats.core_lookup_tasks += cores.tasks as u64;
            stats.core_requests_attempted += cores.requests.attempted as u64;
            stats.core_requests_covered += cores.requests.covered as u64;
            stats.core_requests_rejected += cores.requests.rejected as u64;
            stats.core_requests_uncovered += cores.requests.uncovered as u64;
            stats.core_requests_prescreened += cores.requests.prescreened as u64;
            stats.core_requests_retained += cores.requests.retained as u64;
            stats.core_requests_planned += cores.requests.planned as u64;
            stats.core_lookup_ns += cores.lookup_ns;
            stats.core_lookup_peak_bytes =
                stats.core_lookup_peak_bytes.max(cores.peak_capacity_bound);
            stats.core_lookup_retained_bytes =
                stats.core_lookup_retained_bytes.max(cores.capacity_bytes());
        }
        prepared.par_iter_mut().for_each(|query| {
            query.shared_cores = cores.clone();
            if let Some(cores) = &cores {
                query.positions_by_key.directory.retain(|entry| {
                    cores
                        .groups
                        .binary_search_by_key(&(entry.key as u32), |group| group.key().core)
                        .is_ok()
                });
            }
        });
        let mut conversion_bytes = 0usize;
        let mut capacity_bound = 0usize;
        for query in prepared.iter() {
            let positions = &query.positions_by_key;
            let token_bytes = positions
                .tokens
                .as_ref()
                .map_or(0, |tokens| tokens.capacity() * std::mem::size_of::<u64>());
            let core_bytes = if positions.tokens.is_some() {
                positions
                    .directory
                    .iter()
                    .map(|entry| entry.count)
                    .sum::<usize>()
                    * std::mem::size_of::<QuerySeed>()
            } else {
                positions.core.capacity() * std::mem::size_of::<QuerySeed>()
            };
            if positions.tokens.is_some() {
                conversion_bytes = conversion_bytes
                    .checked_add(token_bytes + core_bytes)
                    .ok_or(TraceError::Invalid("query core conversion size"))?;
            }
            capacity_bound = capacity_bound
                .checked_add(
                    token_bytes
                        + core_bytes
                        + positions.directory.capacity() * std::mem::size_of::<QueryKeyRange>()
                        + query.query.capacity(),
                )
                .ok_or(TraceError::Invalid("query core capacity bound"))?;
            if self.observed {
                let mut stats = self.batch_stats.lock().unwrap();
                stats.query_token_capacity_bytes += token_bytes as u64;
                stats.query_compact_core_queries += u64::from(positions.tokens.is_some());
                stats.query_wide_core_queries += u64::from(positions.tokens.is_none());
            }
        }
        if let Some(scratch) = &mut scratch {
            scratch.retain(conversion_bytes);
        }
        if self.observed {
            let mut stats = self.batch_stats.lock().unwrap();
            stats.query_core_conversion_capacity_bound = stats
                .query_core_conversion_capacity_bound
                .max(capacity_bound as u64);
        }
        prepared
            .par_iter_mut()
            .try_for_each(|query| query.positions_by_key.materialize_tokens())?;
        drop(scratch);
        if core_phase.is_some() {
            let elapsed = phase_elapsed(core_phase);
            let mut stats = self.batch_stats.lock().unwrap();
            for (total, value) in stats.phase_core_lookup_ns.iter_mut().zip(elapsed) {
                *total += value;
            }
        }
        let context_phase = phase_stamp();
        let result = prepared
            .par_iter_mut()
            .zip(circular)
            .try_for_each(|(query, &circular)| {
                let count: usize = query
                    .positions_by_key
                    .directory
                    .iter()
                    .map(|entry| entry.count)
                    .sum();
                let scratch = count
                    .checked_mul(std::mem::size_of::<usize>())
                    .and_then(|bytes| CacheReservation::acquire(&LOOKUP_CACHE_AVAILABLE, bytes));
                let mut positions = Vec::new();
                if scratch.is_some() {
                    positions
                        .try_reserve_exact(count)
                        .map_err(|_| TraceError::Invalid("surviving core positions"))?;
                    for entry in &query.positions_by_key.directory {
                        positions.extend(entry.first..entry.first + entry.count);
                    }
                    positions.sort_unstable_by_key(|&ordinal| {
                        query.positions_by_key.core[ordinal].position
                    });
                }
                let mut nested = Vec::new();
                let mut ordered = positions.iter().copied();
                let mut fallback = query
                    .positions_by_key
                    .directory
                    .iter()
                    .flat_map(|entry| entry.first..entry.first + entry.count);
                for ordinal in std::iter::from_fn(|| {
                    if scratch.is_some() {
                        ordered.next()
                    } else {
                        fallback.next()
                    }
                }) {
                    let seed = query.positions_by_key.core[ordinal];
                    let context = crate::shared_seed::context_seed(
                        &query.query,
                        seed.position as usize,
                        seed.packed_key as u32,
                        seed.canonical_orientation,
                        circular,
                    )
                    .ok_or(TraceError::Invalid("query core context"))?;
                    for length in [21, 31] {
                        if let Some(key) = context.key(length) {
                            let flank = u64::from((length - 15) / 2);
                            let position = if circular {
                                (seed.position + query.query_length - flank) % query.query_length
                            } else {
                                seed.position
                                    .checked_sub(flank)
                                    .ok_or(TraceError::Invalid("query context start"))?
                            };
                            nested.push(QuerySeed {
                                packed_key: key
                                    .packed()
                                    .ok_or(TraceError::Invalid("query context key"))?,
                                position,
                                canonical_orientation: seed.canonical_orientation,
                            });
                        }
                    }
                }
                drop(positions);
                drop(scratch);
                query
                    .positions_by_key
                    .add_nested(nested, query.query_length);
                query.positions_by_key.directory.shrink_to_fit();
                query.lookup_identity =
                    sha256(&[query.lookup_identity.as_slice(), b"shared-contexts-v1"].concat());
                self.record_prepared(query, circular);
                Ok(())
            });
        if context_phase.is_some() {
            let elapsed = phase_elapsed(context_phase);
            let mut stats = self.batch_stats.lock().unwrap();
            for (total, value) in stats.phase_context_generation_ns.iter_mut().zip(elapsed) {
                *total += value;
            }
        }
        result
    }

    fn record_prepared(&self, prepared: &PreparedQuery, circular: bool) {
        if !self.observed {
            return;
        }
        let executed = prepared
            .positions_by_key
            .values()
            .map(|positions| positions.len() as u64)
            .sum::<u64>();
        let associations = if self.index.is_shared() {
            prepared.positions_by_key.core_occurrences as u64
                + possible_context_associations(&prepared.query, circular)
        } else {
            executed
        };
        let mut stats = self.batch_stats.lock().unwrap();
        stats.query_context_associations += associations;
        stats.query_executed_context_associations += executed;
        stats.query_distinct_context_requests += prepared.positions_by_key.len() as u64;
        stats.query_sequence_capacity_bytes += prepared.query.capacity() as u64;
        stats.query_core_capacity_bytes +=
            (prepared.positions_by_key.core.capacity() * std::mem::size_of::<QuerySeed>()) as u64;
        stats.query_nested_capacity_bytes +=
            (prepared.positions_by_key.nested.capacity() * std::mem::size_of::<QuerySeed>()) as u64;
        stats.query_directory_capacity_bytes += (prepared.positions_by_key.directory.capacity()
            * std::mem::size_of::<QueryKeyRange>())
            as u64;
        if self.index.is_shared() {
            stats.query_core_occurrences += prepared.positions_by_key.core_occurrences as u64;
            stats.query_distinct_cores += prepared.positions_by_key.core_distinct as u64;
            stats.nested_context_calls += prepared
                .positions_by_key
                .iter()
                .filter(|(key, _)| **key >> 62 == 0)
                .map(|(_, positions)| positions.len() as u64)
                .sum::<u64>();
        }
    }

    pub(crate) fn search_batch(
        &self,
        queries: &[(String, Vec<u8>)],
        config: TraceConfig,
    ) -> Result<Vec<TraceResult>, TraceError> {
        self.search_batch_topologies(queries, config, &vec![config.circular; queries.len()])
    }

    pub(crate) fn search_batch_topologies(
        &self,
        queries: &[(String, Vec<u8>)],
        config: TraceConfig,
        circular: &[bool],
    ) -> Result<Vec<TraceResult>, TraceError> {
        if queries.len() != circular.len() {
            return Err(TraceError::Invalid("query topology count"));
        }
        let started = self.observed.then(Instant::now);
        if queries.len() < 2 && !self.index.is_shared() {
            return queries
                .par_iter()
                .zip(circular)
                .map(|((id, sequence), &circular)| {
                    self.search(id.as_str(), sequence, TraceConfig { circular, ..config })
                })
                .collect();
        }
        let extraction_phase = phase_stamp();
        let scratch = self
            .index
            .is_shared()
            .then(|| {
                reserve_query_tokens(
                    queries.iter().map(|(_, sequence)| sequence.len()),
                    lookup_budget(&self.index),
                )
            })
            .flatten();
        let operation = self.extraction_operation(scratch.is_some())?;
        let mut prepared = queries
            .par_iter()
            .zip(circular)
            .enumerate()
            .map(|(ordinal, ((id, sequence), &circular))| {
                let mut prepared = if self.index.is_shared() {
                    prepare_query_screened(
                        id.as_str(),
                        sequence,
                        TraceConfig { circular, ..config },
                        15,
                        false,
                        scratch.is_some(),
                        operation.as_ref(),
                        self.observed,
                    )?
                } else {
                    self.prepare(id.as_str(), sequence, TraceConfig { circular, ..config })?
                };
                prepared.batch_ordinal = ordinal;
                Ok::<_, TraceError>(prepared)
            })
            .collect::<Result<Vec<_>, _>>()?;
        if extraction_phase.is_some() {
            let elapsed = phase_elapsed(extraction_phase);
            let mut stats = self.batch_stats.lock().unwrap();
            for (total, value) in stats.phase_extraction_ns.iter_mut().zip(elapsed) {
                *total += value;
            }
        }
        if self.index.is_shared() {
            self.prepare_shared_queries(&mut prepared, circular, scratch, operation.as_ref())?;
            if let Some(operation) = operation {
                operation.finish()?;
            }
            if let Some(started) = started {
                self.batch_stats.lock().unwrap().seed_generation_ns +=
                    started.elapsed().as_nanos() as u64;
            }
        }
        let lookup_phase = phase_stamp();
        let batch = self.prepare_batch(&prepared)?;
        if lookup_phase.is_some() {
            let elapsed = phase_elapsed(lookup_phase);
            let mut stats = self.batch_stats.lock().unwrap();
            for (total, value) in stats.phase_context_lookup_ns.iter_mut().zip(elapsed) {
                *total += value;
            }
        }
        let downstream_phase = phase_stamp();
        let results = prepared
            .into_par_iter()
            .zip(circular)
            .map(|(prepared, &circular)| {
                self.search_prepared(prepared, TraceConfig { circular, ..config }, Some(&batch))
            })
            .collect();
        self.record_batch(&batch);
        if downstream_phase.is_some() {
            let elapsed = phase_elapsed(downstream_phase);
            let mut stats = self.batch_stats.lock().unwrap();
            for (total, value) in stats.phase_downstream_ns.iter_mut().zip(elapsed) {
                *total += value;
            }
        }
        if let Some(started) = started {
            self.batch_stats.lock().unwrap().search_critical_ns +=
                started.elapsed().as_nanos() as u64;
        }
        results
    }

    fn prepare_batch(&self, prepared: &[PreparedQuery]) -> Result<TraceBatch, TraceError> {
        let started = self.observed.then(Instant::now);
        let count = prepared
            .iter()
            .try_fold(0usize, |sum, query| {
                sum.checked_add(query.positions_by_key.len())
            })
            .ok_or(TraceError::Invalid("batch query key count"))?;
        let mut keys = Vec::new();
        let lookups = if lookup_bytes(&self.index, count, prepared.len())
            .is_some_and(|bytes| bytes <= lookup_budget(&self.index))
        {
            keys.try_reserve_exact(count)
                .map_err(|_| TraceError::Invalid("batch query key allocation"))?;
            for (ordinal, query) in prepared.iter().enumerate() {
                keys.extend(query.positions_by_key.keys().map(|&key| (key, ordinal)));
            }
            let setup_ns = started.map_or(0, |started| started.elapsed().as_nanos() as u64);
            let cores = prepared
                .first()
                .and_then(|query| query.shared_cores.as_deref());
            let mut lookups =
                prepare_lookup_with_cores(&self.index, keys, prepared.len(), self.observed, cores)?;
            if let Some(lookups) = &mut lookups {
                lookups.lookup_ns += setup_ns;
            }
            lookups
        } else {
            None
        };
        Ok(TraceBatch {
            lookups,
            sequence: Arc::new(
                BgzfBlockCache::new(DEFAULT_BATCH_BGZF_CACHE_BYTES)
                    .ok_or(TraceError::Invalid("batch sequence cache budget"))?,
            ),
        })
    }

    fn record_batch(&self, batch: &TraceBatch) {
        let cache = batch.sequence.stats();
        let mut stats = self
            .batch_stats
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        stats.batches += 1;
        stats.timings_observed = self.observed;
        stats.bgzf_cache_hits += cache.hits;
        stats.bgzf_blocks_decoded += cache.blocks_decoded;
        stats.bgzf_evictions += cache.evictions;
        stats.bgzf_peak_bytes = stats.bgzf_peak_bytes.max(cache.peak_accounted_bytes);
        stats.bgzf_io_limit = cache.max_loading_blocks;
        stats.bgzf_io_peak = stats.bgzf_io_peak.max(cache.peak_loading_blocks);
        if let Some(lookups) = &batch.lookups {
            stats.unique_keys += lookups.attempted_keys as u64;
            stats.distinct_cores += lookups.distinct_cores;
            stats.split_core_resolutions += lookups.split_core_resolutions;
            stats.restored_query_context_requests += lookups.query_entries.len() as u64;
            stats.lookup_tasks += lookups.lookup_tasks as u64;
            stats.lookup_plan_hash =
                stats.lookup_plan_hash.wrapping_mul(0x100_0000_01b3) ^ lookups.lookup_plan_hash;
            stats.lookup_dispatch_ns += lookups.lookup_dispatch_ns;
            stats.lookup_parallel_ns += lookups.lookup_parallel_ns;
            stats.lookup_compute_ns += lookups.lookup_compute_ns;
            stats.lookup_reduce_ns += lookups.lookup_reduce_ns;
            stats.lookup_dispatch_to_start_ns += lookups.lookup_dispatch_to_start_ns;
            for bucket in 0..16 {
                stats.context_reuse_histogram_log2[bucket] +=
                    lookups.context_reuse_histogram_log2[bucket];
                stats.context_occurrence_histogram_log2[bucket] +=
                    lookups.context_occurrence_histogram_log2[bucket];
            }
            stats.cached_groups += lookups.postings.iter().flatten().count() as u64;
            stats.cached_positions += lookups
                .postings
                .iter()
                .flatten()
                .filter_map(|p| p.occurrences.as_ref())
                .flatten()
                .map(|positions| positions.len() as u64)
                .sum::<u64>();
            stats.lookup_peak_bytes = stats.lookup_peak_bytes.max(lookups.peak_capacity_bound);
            stats.lookup_retained_bytes = stats.lookup_retained_bytes.max(lookups.capacity_bytes);
            stats.lookup_reserved_bytes =
                stats.lookup_reserved_bytes.max(lookups._reservation.bytes);
            stats.key_lookup_ns += lookups.lookup_ns;
            for (total, value) in stats
                .phase_postings_ns
                .iter_mut()
                .zip(lookups.phase_postings_ns)
            {
                *total += value;
            }
            stats.membership_access_ns += lookups.membership_ns;
            stats.position_access_ns += lookups.position_ns;
            let posting = &lookups.posting_execution;
            stats.posting_plan_ns += posting.plan_ns;
            stats.posting_member_tasks += posting.member_tasks as u64;
            stats.posting_position_tasks += posting.position_tasks as u64;
            stats.posting_plan_hash =
                stats.posting_plan_hash.wrapping_mul(0x100_0000_01b3) ^ posting.task_hash;
            stats.posting_admitted_member_rows += posting.admitted_member_rows;
            stats.posting_admitted_position_rows += posting.admitted_position_rows;
            stats.posting_peak_parallel_tasks = stats
                .posting_peak_parallel_tasks
                .max(posting.peak_parallel_tasks as u64);
            stats.posting_scratch_bytes = stats
                .posting_scratch_bytes
                .max(posting.scratch_bytes as u64);
            stats.posting_member_copies += posting.member_copies;
            stats.posting_initialized_position_bytes += posting.initialized_position_bytes as u64;
            stats.posting_member_worker_elapsed_ns += posting.member_worker_elapsed_ns;
            stats.posting_position_worker_elapsed_ns += posting.position_worker_elapsed_ns;
            stats.posting_member_worker_cpu_ns += posting.member_worker_cpu_ns;
            stats.posting_position_worker_cpu_ns += posting.position_worker_cpu_ns;
            stats.posting_cpu_unavailable_tasks += posting.cpu_unavailable as u64;
            if let Some((lists, rows)) = lookups.unused_positions() {
                stats.posting_unused_member_lists += lists;
                stats.posting_unused_position_rows += rows;
            } else {
                stats.posting_usage_unavailable_batches += 1;
            }
        }
    }

    pub fn batch_stats(&self) -> TraceBatchStats {
        *self
            .batch_stats
            .lock()
            .unwrap_or_else(|error| error.into_inner())
    }

    pub fn shared_read_stats(&self) -> Option<crate::shared_reader::SharedReadStats> {
        match &self.index {
            TraceIndex::Shared(index) => Some(index.stats()),
            _ => None,
        }
    }

    fn downstream_cpu(&self, stage: usize, start: Option<u64>) {
        if start.is_some() {
            let elapsed = crate::alignment::elapsed_cpu(start);
            self.batch_stats.lock().unwrap().downstream_thread_cpu_ns[stage] += elapsed;
        }
    }

    fn record_alignment_time(
        &self,
        started: Option<Instant>,
        workspace: &crate::alignment::AlignmentWorkspace,
        before: (u64, u64),
    ) {
        if let Some(started) = started {
            let elapsed = started.elapsed().as_nanos() as u64;
            let traceback = workspace.traceback_nanoseconds().saturating_sub(before.0);
            let endpoint = workspace.endpoint_nanoseconds().saturating_sub(before.1);
            let mut stats = self.batch_stats.lock().unwrap();
            let mut work = workspace.work;
            work.capacity_bytes = workspace.retained_bytes() as u64;
            stats.alignment_work.add(work);
            stats.alignment_ns += elapsed.saturating_sub(traceback).saturating_sub(endpoint);
            stats.traceback_and_cigar_ns += traceback;
            stats.endpoint_completion_ns += endpoint;
        }
    }

    fn search_prepared(
        &self,
        prepared: PreparedQuery,
        config: TraceConfig,
        batch: Option<&TraceBatch>,
    ) -> Result<TraceResult, TraceError> {
        let cpu = crate::alignment::observed_cpu(self.observed);
        let shared = batch.and_then(|batch| batch.lookups.as_ref());
        let cache_bytes = LOOKUP_CACHE_BYTES / rayon::current_num_threads().max(1);
        let census = self.candidate_census_with_shared(&prepared, config, cache_bytes, shared)?;
        let completion = candidate_completion(census.candidates.len(), config.max_metagenomes)?;
        let mut candidates = census.candidates;
        candidates.truncate(config.max_metagenomes);
        let key_order = if let Some(shared) = shared {
            let mut ordinals =
                shared.query_entries[shared.query_ranges[prepared.batch_ordinal].clone()].to_vec();
            ordinals.sort_unstable_by_key(|&ordinal| {
                let (key, seed) = shared.entries[ordinal as usize];
                (
                    std::cmp::Reverse(if self.index.is_shared() { key >> 62 } else { 0 }),
                    seed.unwrap().document_frequency(),
                    key,
                )
            });
            ordinals
        } else {
            let mut frequencies = census.frequencies;
            frequencies.sort_unstable_by_key(|&(key, frequency)| {
                (
                    std::cmp::Reverse(if self.index.is_shared() { key >> 62 } else { 0 }),
                    frequency,
                    key,
                )
            });
            frequencies
                .into_iter()
                .map(|(key, _)| key)
                .collect::<Vec<_>>()
        };
        let candidates_screened =
            u32::try_from(candidates.len()).map_err(|_| TraceError::Invalid("candidate count"))?;
        self.downstream_cpu(0, cpu);
        let metagenomes = self.trace_selected_with_batch(
            &prepared,
            candidates,
            &key_order,
            census.lookups.as_ref(),
            config,
            batch,
        )?;
        Ok(TraceResult {
            query_id: prepared.query_id,
            query_length: prepared.query_length,
            index: TraceIndexIdentity {
                manifest_sha256: digest_hex(self.index.manifest_sha256()),
                body_sha256: digest_hex(self.index.body_sha256()),
                seed_k: self.index.k(),
                rescue_k15: self.index.rescue_k15(),
            },
            completion,
            candidates_screened,
            metagenomes,
        })
    }

    pub(crate) fn candidate_census(
        &self,
        prepared: &PreparedQuery,
        config: TraceConfig,
        cache_bytes: usize,
    ) -> Result<TraceCensus, TraceError> {
        self.candidate_census_with_shared(prepared, config, cache_bytes, None)
    }

    fn candidate_census_with_shared(
        &self,
        prepared: &PreparedQuery,
        config: TraceConfig,
        cache_bytes: usize,
        shared: Option<&SharedSeedLookups>,
    ) -> Result<TraceCensus, TraceError> {
        let started = self.observed.then(Instant::now);
        validate_config(config)?;
        if let Some(shared) = shared
            && self.index.cache_file_identity()? != Some(shared.identity)
        {
            return Err(TraceError::Invalid("shared seed lookup identity"));
        }
        let sketch_candidates =
            self.screen_candidates(&prepared.query_id, &prepared.query, config)?;
        let mut candidates = BTreeMap::new();
        for candidate in sketch_candidates {
            candidates.insert(candidate.id, candidate);
        }
        self.index
            .verify_query_filter_pages(prepared.positions_by_key.keys())?;
        let mut frequencies = Vec::new();
        let header_sha256 = self.index.header_sha256()?;
        let mut lookups = if shared.is_some_and(|shared| shared.postings_complete) {
            None
        } else {
            self.index
                .cache_file_identity()
                .ok()
                .flatten()
                .and_then(|file_identity| {
                    CachedSeedLookups::new(
                        header_sha256,
                        prepared.lookup_identity,
                        file_identity,
                        cache_bytes,
                        prepared.positions_by_key.len(),
                        self.index.document_count(),
                    )
                })
        };
        let mut entries = prepared
            .positions_by_key
            .iter()
            .map(|(&key, seeds)| (key, seeds));
        let mut matching_keys = shared.map(|shared| {
            shared.query_entries[shared.query_ranges[prepared.batch_ordinal].clone()].iter()
        });
        let mut resolved_seeds = Vec::new();
        let mut chunk = Vec::new();
        chunk
            .try_reserve_exact(SEED_LOOKUP_BATCH_KEYS)
            .map_err(|_| TraceError::Invalid("query seed batch"))?;
        let mut packed_keys = Vec::new();
        packed_keys
            .try_reserve_exact(SEED_LOOKUP_BATCH_KEYS)
            .map_err(|_| TraceError::Invalid("query seed batch"))?;
        loop {
            chunk.clear();
            packed_keys.clear();
            resolved_seeds.clear();
            let posting_ordinals = matching_keys.as_ref().map(|keys| keys.as_slice());
            for _ in 0..SEED_LOOKUP_BATCH_KEYS {
                let next = if let Some(matching_keys) = &mut matching_keys {
                    matching_keys.next().map(|&ordinal| {
                        let (key, seed) = shared.unwrap().entries[ordinal as usize];
                        resolved_seeds.push(seed);
                        (key, &prepared.positions_by_key[&key])
                    })
                } else {
                    entries.next()
                };
                let Some((packed_key, query_seeds)) = next else {
                    break;
                };
                chunk.push((packed_key, query_seeds));
                packed_keys.push(packed_key);
            }
            if chunk.is_empty() {
                break;
            }

            if shared.is_none() {
                resolved_seeds = prepared.find_seeds(&self.index, &packed_keys)?;
            }
            let index_seeds = &resolved_seeds;
            for (seed_ordinal, ((packed_key, query_seeds), &index_seed)) in
                chunk.iter().copied().zip(index_seeds).enumerate()
            {
                let Some(index_seed) = index_seed else {
                    if let Some(lookups) = &mut lookups {
                        lookups.cache_negative(packed_key);
                    }
                    continue;
                };
                frequencies.push((packed_key, index_seed.document_frequency()));
                let query_positions = u64::try_from(query_seeds.len())
                    .map_err(|_| TraceError::Invalid("query seed count"))?;
                let decoded_documents;
                let posting_ordinal =
                    posting_ordinals.map(|ordinals| ordinals[seed_ordinal] as usize);
                let posting = shared
                    .map(|s| s.posting(packed_key, posting_ordinal))
                    .transpose()?
                    .flatten();
                let documents = if let Some(posting) = posting {
                    &posting.documents
                } else {
                    decoded_documents = self.index.seed_documents(index_seed)?;
                    &decoded_documents
                };
                for &document in documents {
                    let hits = document
                        .occurrence_count()
                        .checked_mul(query_positions)
                        .ok_or(TraceError::Invalid("exact seed hit count"))?;
                    if let Some(candidate) = candidates.get_mut(&document.metagenome_id()) {
                        candidate.exact_seed_hits = candidate
                            .exact_seed_hits
                            .checked_add(hits)
                            .ok_or(TraceError::Invalid("exact seed hit count"))?;
                        continue;
                    }
                    let name = self
                        .index
                        .metagenome_name(document.metagenome_id())?
                        .ok_or(TraceError::Invalid("JIDX metagenome ID"))?
                        .to_string();
                    candidates.insert(
                        document.metagenome_id(),
                        Candidate {
                            id: document.metagenome_id(),
                            name,
                            shared_hashes: 0,
                            containment: 0.0,
                            exact_seed_hits: hits,
                        },
                    );
                }
                if let Some(lookups) = &mut lookups {
                    lookups.cache_group(packed_key, index_seed, documents);
                }
            }
        }
        let mut candidates = candidates.into_values().collect::<Vec<_>>();
        candidates.sort_by(compare_candidates);
        if let Some(started) = started {
            self.batch_stats.lock().unwrap().candidate_routing_ns +=
                started.elapsed().as_nanos() as u64;
        }
        Ok(TraceCensus {
            candidates,
            frequencies,
            lookups,
        })
    }

    pub(crate) fn trace_selected(
        &self,
        prepared: &PreparedQuery,
        candidates: Vec<Candidate>,
        key_order: &[u64],
        lookups: Option<&CachedSeedLookups>,
        config: TraceConfig,
    ) -> Result<Vec<MetagenomeTrace>, TraceError> {
        self.trace_selected_with_batch(prepared, candidates, key_order, lookups, config, None)
    }

    fn trace_selected_with_batch(
        &self,
        prepared: &PreparedQuery,
        candidates: Vec<Candidate>,
        key_order: &[u64],
        lookups: Option<&CachedSeedLookups>,
        config: TraceConfig,
        batch: Option<&TraceBatch>,
    ) -> Result<Vec<MetagenomeTrace>, TraceError> {
        let started = self.observed.then(Instant::now);
        let geometry_cpu = crate::alignment::observed_cpu(self.observed);
        validate_config(config)?;
        self.index.enable_selected_front_metadata();
        // Contig metadata in the hit loop is read without per-record identity checks; the
        // file identity is checked here and again after the loop.
        self.index.cache_file_identity()?;
        let lookups = if let Some(lookups) = lookups {
            if lookups.header_sha256 != self.index.header_sha256()?
                || lookups.query_identity != prepared.lookup_identity
            {
                return Err(TraceError::Invalid("cached seed lookup identity"));
            }
            match self.index.cache_file_identity() {
                Ok(Some(identity)) if identity == lookups.file_identity => Some(lookups),
                Ok(Some(_)) => return Err(TraceError::Invalid("cached seed lookup identity")),
                _ => None,
            }
        } else {
            None
        };
        let candidate_ids = candidates
            .iter()
            .map(|candidate| candidate.id)
            .collect::<HashSet<_>>();
        let mut region_hits = BTreeMap::<RegionKey, RegionHits>::new();
        let metadata_reservation = CacheReservation::acquire(&LOOKUP_CACHE_AVAILABLE, 1024 * 1024);
        let metadata_limit = metadata_reservation
            .as_ref()
            .map_or(0, |reservation| (reservation.bytes - 4096) / 256);
        let mut numeric_contigs = BTreeMap::new();
        type SharedGeometry = (
            crate::shared_reader::SharedOccurrenceStorage,
            u64,
            u64,
            u64,
            bool,
        );
        let geometry_reservation = self
            .index
            .is_shared()
            .then(|| CacheReservation::acquire(&LOOKUP_CACHE_AVAILABLE, 1024 * 1024))
            .flatten();
        let geometry_row_bytes = 4 * std::mem::size_of::<SharedGeometry>();
        let geometry_limit = geometry_reservation.as_ref().map_or(0, |reservation| {
            (reservation.bytes - 4096) / geometry_row_bytes
        });
        let mut geometries = BTreeSet::<SharedGeometry>::new();
        // Ledger runs record 21/31 contexts per core pair, so repeated geometries are replayed.
        let mut region_ledger = self.region_ledger.as_ref().map(|_| RegionLedger::default());
        let mut routing_reserved_bytes = 4096usize;
        let batch_lookups = batch.and_then(|batch| batch.lookups.as_ref());
        if let Some(shared) = batch_lookups
            && self.index.cache_file_identity()? != Some(shared.identity)
        {
            return Err(TraceError::Invalid("batch lookup identity"));
        }
        let mut packed_keys = Vec::new();
        packed_keys
            .try_reserve_exact(SEED_LOOKUP_BATCH_KEYS)
            .map_err(|_| TraceError::Invalid("query seed batch"))?;
        let mut filter_keys = Vec::new();
        filter_keys
            .try_reserve_exact(SEED_LOOKUP_BATCH_KEYS)
            .map_err(|_| TraceError::Invalid("query seed batch"))?;
        let mut emitted_anchor_associations = 0u64;
        let mut executed_anchor_associations = 0u64;
        for key_chunk in key_order.chunks(SEED_LOOKUP_BATCH_KEYS) {
            packed_keys.clear();
            for &packed_key in key_chunk {
                if let Some(shared) = batch_lookups {
                    packed_keys.push(shared.entries[packed_key as usize].0);
                } else if prepared.positions_by_key.contains_key(&packed_key) {
                    packed_keys.push(packed_key);
                }
            }
            if packed_keys.is_empty() {
                continue;
            }
            filter_keys.clear();
            filter_keys.extend(packed_keys.iter().copied().filter(|&key| {
                batch_lookups.is_none()
                    && lookups
                        .is_none_or(|lookups| lookups.through.is_none_or(|through| key > through))
            }));
            filter_keys.sort_unstable();
            filter_keys.dedup();
            self.index.verify_query_filter_pages(filter_keys.iter())?;
            let uncached_seeds = if filter_keys.is_empty() {
                Vec::new()
            } else {
                prepared.find_seeds(&self.index, &filter_keys)?
            };
            let index_seeds = packed_keys
                .iter()
                .enumerate()
                .map(|(ordinal, &key)| {
                    if let Some(shared) = batch_lookups {
                        return Ok((shared.entries[key_chunk[ordinal] as usize].1, None));
                    }
                    Ok(match lookups.and_then(|lookups| lookups.get(key)) {
                        Some(Some(cached)) => (Some(cached.seed), Some(cached.documents)),
                        Some(None) => (None, None),
                        None => (
                            uncached_seeds
                                [filter_keys.binary_search(&key).expect("uncached query key")],
                            None,
                        ),
                    })
                })
                .collect::<Result<Vec<_>, TraceError>>()?;

            for (seed_ordinal, (&packed_key, (index_seed, cached_documents))) in
                packed_keys.iter().zip(index_seeds).enumerate()
            {
                let query_seeds = prepared
                    .positions_by_key
                    .get(&packed_key)
                    .expect("present query key");
                let Some(index_seed) = index_seed else {
                    continue;
                };
                let seed_k = self.index.seed_length(packed_key)?;
                let decoded_documents;
                let posting = batch_lookups
                    .map(|s| s.posting(packed_key, Some(key_chunk[seed_ordinal] as usize)))
                    .transpose()?
                    .flatten();
                let documents = if let Some(posting) = posting {
                    &posting.documents
                } else if let Some(documents) = cached_documents {
                    documents
                } else {
                    decoded_documents = self.index.seed_documents(index_seed)?;
                    &decoded_documents
                };
                for (document_ordinal, document) in documents
                    .iter()
                    .copied()
                    .enumerate()
                    .filter(|(_, document)| candidate_ids.contains(&document.metagenome_id()))
                {
                    let storage = match document {
                        SeedDocument::Shared(member) => Some(member.occurrence_storage_identity()),
                        _ => None,
                    };
                    let mut occurrence_start = 0u64;
                    let mut visit = |occurrences: &[crate::jidx_reader::SeedOccurrence]| {
                        let start = occurrence_start;
                        occurrence_start += occurrences.len() as u64;
                        if self.observed {
                            emitted_anchor_associations +=
                                (query_seeds.len() as u64) * (occurrences.len() as u64);
                        }
                        for occurrence in occurrences {
                            if numeric_contigs.len() < metadata_limit
                                && let std::collections::btree_map::Entry::Vacant(slot) =
                                    numeric_contigs.entry(occurrence.contig_id)
                            {
                                let contig = self
                                    .index
                                    .numeric_contig_unchecked(occurrence.contig_id)?
                                    .ok_or(TraceError::Invalid("missing occurrence contig"))?;
                                slot.insert(contig);
                            }
                        }
                        for seed in query_seeds {
                            let (query_position, region_k) = if self.index.is_shared() {
                                (
                                    (seed.position + u64::from((seed_k - 15) / 2))
                                        % prepared.query_length,
                                    15,
                                )
                            } else {
                                (seed.position, seed_k)
                            };
                            if let Some(storage) = storage {
                                let geometry = (
                                    storage,
                                    start,
                                    occurrences.len() as u64,
                                    query_position,
                                    seed.canonical_orientation,
                                );
                                if geometries.len() < geometry_limit {
                                    if !geometries.insert(geometry) && region_ledger.is_none() {
                                        continue;
                                    }
                                } else if geometries.contains(&geometry) && region_ledger.is_none()
                                {
                                    continue;
                                }
                            }
                            if self.observed {
                                executed_anchor_associations += occurrences.len() as u64;
                            }
                            for occurrence in occurrences {
                                let contig = if let Some(contig) =
                                    numeric_contigs.get(&occurrence.contig_id)
                                {
                                    *contig
                                } else {
                                    self.index
                                        .numeric_contig_unchecked(occurrence.contig_id)?
                                        .ok_or(TraceError::Invalid("missing occurrence contig"))?
                                };
                                let strand = if seed.canonical_orientation
                                    == occurrence.canonical_orientation
                                {
                                    Strand::Forward
                                } else {
                                    Strand::Reverse
                                };
                                let oriented_position = match strand {
                                    Strand::Forward => occurrence.position,
                                    Strand::Reverse => contig
                                        .length
                                        .checked_sub(
                                            occurrence
                                                .position
                                                .checked_add(u64::from(region_k))
                                                .ok_or(TraceError::Invalid(
                                                    "occurrence position",
                                                ))?,
                                        )
                                        .ok_or(TraceError::Invalid("occurrence position"))?,
                                };
                                let diagonal =
                                    i128::from(oriented_position) - i128::from(query_position);
                                let region = region_hits
                                    .entry(RegionKey {
                                        metagenome_id: contig.metagenome_id,
                                        contig_id: contig.id,
                                        strand,
                                        k: region_k,
                                    })
                                    .or_default();
                                let hit = SeedHit {
                                    query: query_position,
                                    target: oriented_position,
                                    diagonal,
                                };
                                if let Some(ledger) = region_ledger.as_mut() {
                                    let context = match seed_k {
                                        21 => 2,
                                        31 => 4,
                                        _ => 1,
                                    };
                                    let key = RegionKey {
                                        metagenome_id: contig.metagenome_id,
                                        contig_id: contig.id,
                                        strand,
                                        k: region_k,
                                    };
                                    *ledger
                                        .contexts
                                        .entry((key, query_position, oriented_position))
                                        .or_default() |= context;
                                }
                                if self.index.is_shared() {
                                    region.push_unique(hit, &mut routing_reserved_bytes)?;
                                } else {
                                    region.push(hit);
                                }
                            }
                        }
                        Ok(())
                    };
                    if let Some(positions) = posting.and_then(|p| p.occurrences.as_ref()) {
                        if let Some(shared) = batch_lookups {
                            shared.mark_positions_used(
                                key_chunk[seed_ordinal] as usize,
                                document_ordinal,
                            );
                        }
                        visit(&positions[document_ordinal])?;
                    } else {
                        self.index.visit_occurrences(index_seed, document, visit)?;
                    }
                }
            }
        }
        #[cfg(test)]
        tests::run_hit_loop_end_hook(&prepared.query_id);
        self.index.cache_file_identity()?;
        let mut predecessor_tests = 0;
        let geometric_hits = if self.observed {
            region_hits
                .values()
                .map(|hits| match hits {
                    RegionHits::Empty => 0,
                    RegionHits::One(_) => 1,
                    RegionHits::Many(hits) => hits.len(),
                    RegionHits::Unique(hits) => hits.len(),
                })
                .sum::<usize>()
        } else {
            0
        };
        // Island planning and ledgers keep each region's ordered anchors; charged by the
        // per-pair routing reservation, which covers later SeedHit storage.
        let retain_support = self.region_policy != crate::trace_islands::RegionPolicy::Parent
            || region_ledger.is_some();
        let mut support = Vec::new();
        let regions = if self.observed || retain_support {
            form_regions_observed(
                region_hits,
                config.diagonal_bin_bases,
                Some(&mut predecessor_tests),
                retain_support.then_some(&mut support),
            )
        } else {
            form_regions(region_hits, config.diagonal_bin_bases)
        };
        self.downstream_cpu(1, geometry_cpu);
        let task_cpu = crate::alignment::observed_cpu(self.observed);

        let (tasks, split_parents) = self.tasks(
            regions,
            &support,
            prepared.query_length,
            config,
            region_ledger.as_mut(),
        )?;
        if let Some(started) = started {
            let mut stats = self.batch_stats.lock().unwrap();
            stats.region_predecessor_tests += predecessor_tests;
            stats.geometric_hits += geometric_hits as u64;
            stats.region_formation_ns += started.elapsed().as_nanos() as u64;
            stats.emitted_anchor_associations += emitted_anchor_associations;
            stats.executed_anchor_associations += executed_anchor_associations;
            stats.physical_geometry_capacity_bound = stats
                .physical_geometry_capacity_bound
                .max(4096 + geometries.len() * geometry_row_bytes);
            stats.geometry_phase_reserved_bytes = stats.geometry_phase_reserved_bytes.max(
                metadata_reservation
                    .as_ref()
                    .map_or(0, |reservation| reservation.bytes)
                    + geometry_reservation
                        .as_ref()
                        .map_or(0, |reservation| reservation.bytes),
            );
            stats.numeric_metadata_capacity_bound = stats
                .numeric_metadata_capacity_bound
                .max(4096 + numeric_contigs.len() * 256);
        }
        drop(numeric_contigs);
        drop(metadata_reservation);
        drop(geometries);
        drop(geometry_reservation);
        let (loaded, mut reads) = self.load_ranges(
            &tasks,
            config.verify_resources,
            batch.map(|batch| &batch.sequence),
        )?;
        self.downstream_cpu(2, task_cpu);
        let outcomes = self.align_tasks(&prepared.query, &tasks, &loaded, config)?;
        drop(loaded);
        let fragments = if split_parents.is_empty() && region_ledger.is_none() {
            outcomes
                .into_iter()
                .filter_map(|outcome| outcome.fragment)
                .collect::<Vec<_>>()
        } else {
            self.finish_island_tasks(
                prepared,
                &tasks,
                outcomes,
                &split_parents,
                &support,
                &mut reads,
                region_ledger,
                config,
                batch,
            )?
        };
        let result_cpu = crate::alignment::observed_cpu(self.observed);
        let mut before_dedup = 0u64;
        let mut after_dedup = 0u64;
        let mut by_metagenome = BTreeMap::<MetagenomeId, Vec<Fragment>>::new();
        for (metagenome_id, fragment) in fragments {
            by_metagenome
                .entry(metagenome_id)
                .or_default()
                .push(fragment);
        }

        let mut metagenomes = Vec::with_capacity(candidates.len());
        for candidate in candidates {
            let fragments = by_metagenome.remove(&candidate.id).unwrap_or_default();
            let mosaic = build_mosaic(prepared.query_length, &fragments)?;
            before_dedup += fragments.len() as u64;
            after_dedup += (mosaic.primary.len() + mosaic.alternatives.len()) as u64;
            let mut contig_ids = BTreeSet::new();
            for fragment in mosaic
                .primary
                .iter()
                .map(|selected| &selected.fragment)
                .chain(mosaic.alternatives.iter())
            {
                contig_ids.insert(fragment.contig_id);
            }
            let mut contigs = Vec::with_capacity(contig_ids.len());
            for id in contig_ids {
                let contig = self
                    .index
                    .contig(id)?
                    .ok_or(TraceError::Invalid("missing result contig"))?;
                contigs.push(TraceContig {
                    id,
                    name: contig.name.to_string(),
                });
            }
            let (stats, bgzf_blocks_decoded) =
                reads.get(&candidate.id).copied().unwrap_or_default();
            metagenomes.push(MetagenomeTrace {
                metagenome_id: candidate.id,
                name: candidate.name,
                shared_hashes: candidate.shared_hashes,
                containment: candidate.containment,
                exact_seed_hits: candidate.exact_seed_hits,
                compressed_bytes_read: stats.bytes_read,
                range_requests: stats.read_requests,
                bgzf_blocks_decoded,
                contigs,
                mosaic,
            });
        }
        self.downstream_cpu(4, result_cpu);
        if self.observed {
            let mut stats = self.batch_stats.lock().unwrap();
            stats.fragments_before_dedup += before_dedup;
            stats.fragments_after_dedup += after_dedup;
        }
        Ok(metagenomes)
    }

    fn screen_candidates(
        &self,
        query_id: &str,
        query: &[u8],
        config: TraceConfig,
    ) -> Result<Vec<Candidate>, TraceError> {
        if !config.use_sketch {
            return Ok(Vec::new());
        }
        let Some(screen) = &self.screen else {
            return Ok(Vec::new());
        };
        let sketch = QuerySketch::from_sequence(query_id, query, screen.reader())?;
        let result = screen
            .query_sketch(&sketch)
            .into_iter()
            .next()
            .ok_or(TraceError::Invalid("missing JAM query result"))?;
        result
            .matches
            .into_iter()
            .filter(|candidate| candidate.containment >= config.min_containment)
            .map(|candidate| {
                let id = *self
                    .sample_to_metagenome
                    .get(candidate.sample_id as usize)
                    .ok_or(TraceError::Invalid("JAM sample ID"))?;
                let name = self
                    .index
                    .metagenome_name(id)?
                    .ok_or(TraceError::Invalid("JIDX metagenome ID"))?
                    .to_string();
                Ok(Candidate {
                    id,
                    name,
                    shared_hashes: candidate.hit_count,
                    containment: candidate.containment,
                    exact_seed_hits: 0,
                })
            })
            .collect::<Result<Vec<_>, TraceError>>()
    }

    /// Plans alignment tasks. Under the island policy a split parent contributes its islands to
    /// the returned tasks and its full envelope to the second list, which holds fallbacks only.
    fn tasks(
        &self,
        regions: Vec<(RegionKey, RegionAccumulator)>,
        support: &[SeedHit],
        query_length: u64,
        config: TraceConfig,
        mut ledger: Option<&mut RegionLedger>,
    ) -> Result<(Vec<AlignmentTask>, Vec<AlignmentTask>), TraceError> {
        let mut tasks = Vec::new();
        let mut split_parents = Vec::new();
        let mut parent = 0u32;
        for (key, region) in regions {
            let minimum = if self.index.is_shared() {
                config.min_seed_hits
            } else {
                minimum_region_hits(key.k, config.min_seed_hits)
            };
            if region.hits < minimum {
                continue;
            }
            let contig = self
                .index
                .contig(key.contig_id)?
                .ok_or(TraceError::Invalid("missing region contig"))?;
            let envelope = fragment_envelope(&region, key, query_length, contig.length, config)?;
            let task = |envelope: &FragmentEnvelope, island| AlignmentTask {
                metagenome_id: key.metagenome_id,
                contig_id: key.contig_id,
                strand: key.strand,
                query_start: envelope.query_start,
                query_span: envelope.query_span,
                target_start: envelope.target_start,
                target_end: envelope.target_end,
                diagonal_offset: envelope.diagonal_offset,
                parent,
                island,
            };
            let anchors = &support[region.support.clone()];
            let islands = match self.region_policy {
                crate::trace_islands::RegionPolicy::Parent => None,
                crate::trace_islands::RegionPolicy::Islands => crate::trace_islands::plan_islands(
                    &region,
                    &envelope,
                    anchors,
                    key,
                    query_length,
                    contig.length,
                    config,
                )?,
            };
            if let Some(ledger) = ledger.as_deref_mut() {
                let contexts = anchors
                    .iter()
                    .map(|hit| {
                        ledger
                            .contexts
                            .get(&(key, hit.query, hit.target))
                            .copied()
                            .unwrap_or(0)
                    })
                    .collect::<Vec<_>>();
                ledger.support.push(region.support.clone());
                ledger.parents.push(crate::trace_islands::ParentLedger {
                    query_id: String::new(),
                    query_length,
                    circular: config.circular,
                    policy: self.region_policy.name(),
                    parent,
                    metagenome_id: key.metagenome_id,
                    contig_id: key.contig_id,
                    contig_length: contig.length,
                    strand: key.strand,
                    support: crate::trace_islands::SupportLedger::new(anchors, &contexts),
                    parent_window: crate::trace_islands::WindowLedger::new(&envelope, region.hits),
                    islands: islands
                        .iter()
                        .flatten()
                        .map(|island| {
                            crate::trace_islands::WindowLedger::new(
                                &island.envelope,
                                island.region.hits,
                            )
                        })
                        .collect(),
                    fallback: false,
                    tasks: Vec::new(),
                });
            }
            if let Some(islands) = islands {
                tasks.extend(
                    islands
                        .iter()
                        .map(|island| task(&island.envelope, Some(island.edges))),
                );
                split_parents.push(task(&envelope, None));
            } else {
                tasks.push(task(&envelope, None));
            }
            parent = parent
                .checked_add(1)
                .ok_or(TraceError::Invalid("region count"))?;
        }
        tasks.sort_unstable_by_key(|task| {
            (
                task.metagenome_id,
                task.contig_id,
                task.strand,
                task.target_start,
                task.query_start,
            )
        });
        Ok((tasks, split_parents))
    }

    /// Applies the island fallback rule and returns fragments in parent task order. Any inner-edge
    /// contact of an island reruns its full parent, whose result replaces all of its islands.
    #[allow(clippy::too_many_arguments)]
    fn finish_island_tasks(
        &self,
        prepared: &PreparedQuery,
        tasks: &[AlignmentTask],
        outcomes: Vec<TaskOutcome>,
        split_parents: &[AlignmentTask],
        support: &[SeedHit],
        reads: &mut HashMap<MetagenomeId, (crate::range_source::RangeStats, u64)>,
        mut ledger: Option<RegionLedger>,
        config: TraceConfig,
        batch: Option<&TraceBatch>,
    ) -> Result<Vec<(MetagenomeId, Fragment)>, TraceError> {
        let fallback_parents = tasks
            .iter()
            .zip(&outcomes)
            .filter(|(_, outcome)| outcome.contact)
            .map(|(task, _)| task.parent)
            .collect::<BTreeSet<_>>();
        let fallback = split_parents
            .iter()
            .filter(|task| fallback_parents.contains(&task.parent))
            .cloned()
            .collect::<Vec<_>>();
        let fallback_outcomes = if fallback.is_empty() {
            Vec::new()
        } else {
            let (loaded, fallback_reads) = self.load_ranges(
                &fallback,
                config.verify_resources,
                batch.map(|batch| &batch.sequence),
            )?;
            for (metagenome_id, (stats, blocks)) in fallback_reads {
                let entry = reads.entry(metagenome_id).or_default();
                entry.0.metadata_requests += stats.metadata_requests;
                entry.0.read_requests += stats.read_requests;
                entry.0.bytes_read += stats.bytes_read;
                entry.0.read_nanoseconds = match (entry.0.read_nanoseconds, stats.read_nanoseconds)
                {
                    (Some(first), Some(second)) => Some(first + second),
                    (first, second) => first.or(second),
                };
                entry.1 += blocks;
            }
            self.align_tasks(&prepared.query, &fallback, &loaded, config)?
        };
        let query_length = prepared.query_length;
        let mut keyed = Vec::new();
        let roles = tasks
            .iter()
            .map(|task| {
                if task.island.is_some() {
                    "island"
                } else {
                    "parent"
                }
            })
            .chain(fallback.iter().map(|_| "fallback"));
        for ((role, task), outcome) in roles
            .zip(tasks.iter().chain(&fallback))
            .zip(outcomes.into_iter().chain(fallback_outcomes))
        {
            if let Some(RegionLedger {
                support: ranges,
                parents,
                ..
            }) = ledger.as_mut()
            {
                let parent = &mut parents[task.parent as usize];
                parent.fallback |= role == "fallback";
                let anchors = &support[ranges[task.parent as usize].clone()];
                let observation = outcome.ledger.as_deref();
                parent.tasks.push(crate::trace_islands::TaskLedger {
                    role,
                    window: crate::trace_islands::WindowLedger::new(
                        &FragmentEnvelope {
                            query_start: task.query_start,
                            query_span: task.query_span,
                            target_start: task.target_start,
                            target_end: task.target_end,
                            diagonal_offset: task.diagonal_offset,
                        },
                        0,
                    ),
                    band_cells: observation.map_or(0, |observed| observed.band_cells),
                    work: observation
                        .map(|observed| observed.work)
                        .unwrap_or_default(),
                    core: observation.and_then(|observed| {
                        observed.core.as_ref().map(|alignment| {
                            crate::trace_islands::AlignmentLedger::new(
                                alignment,
                                task.query_start,
                                anchors,
                                query_length,
                            )
                        })
                    }),
                    selected: observation.and_then(|observed| {
                        observed.selected.as_ref().map(|alignment| {
                            crate::trace_islands::AlignmentLedger::new(
                                alignment,
                                observed.query_start,
                                anchors,
                                query_length,
                            )
                        })
                    }),
                    accepted: outcome.fragment.is_some(),
                    contact: outcome.contact,
                });
            }
            if role == "island" && fallback_parents.contains(&task.parent) {
                continue;
            }
            if let Some(fragment) = outcome.fragment {
                let key = (
                    task.metagenome_id,
                    task.contig_id,
                    task.strand,
                    task.target_start,
                    task.query_start,
                );
                keyed.push((key, fragment));
            }
        }
        keyed.sort_by_key(|(key, _)| *key);
        if let (Some(ledger), Some(directory)) = (ledger, &self.region_ledger) {
            use std::io::Write as _;
            let name = format!(
                "{}-{}.jsonl",
                &digest_hex(sha256(prepared.query_id.as_bytes()))[..16],
                self.region_policy.name()
            );
            let file = std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(directory.join(name))?;
            let mut writer = io::BufWriter::new(file);
            for mut parent in ledger.parents {
                parent.query_id.clone_from(&prepared.query_id);
                serde_json::to_writer(&mut writer, &parent)
                    .map_err(|_| TraceError::Invalid("region ledger serialization"))?;
                writer.write_all(b"\n")?;
            }
            writer
                .into_inner()
                .map_err(|error| TraceError::Io(error.into_error()))?
                .sync_all()?;
        }
        Ok(keyed.into_iter().map(|(_, fragment)| fragment).collect())
    }

    fn load_ranges(
        &self,
        tasks: &[AlignmentTask],
        verify: bool,
        cache: Option<&Arc<BgzfBlockCache>>,
    ) -> Result<LoadedRanges, TraceError> {
        let mut spans = BTreeMap::<(MetagenomeId, ContigId), Vec<(u64, u64)>>::new();
        for task in tasks {
            spans
                .entry((task.metagenome_id, task.contig_id))
                .or_default()
                .push((task.target_start, task.target_end));
        }
        for contig_spans in spans.values_mut() {
            coalesce_spans(contig_spans);
        }
        // Metagenomes are read independently; the shared block cache bounds concurrent decodes.
        let metagenomes = tasks
            .iter()
            .map(|task| task.metagenome_id)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let per_metagenome = metagenomes
            .par_iter()
            .map(|&metagenome_id| {
                let source = self
                    .index
                    .metagenome(metagenome_id)?
                    .ok_or(TraceError::Invalid("missing source metagenome"))?;
                let mut reader = if let Some(cache) = cache {
                    BgzfReader::open_with_cache(
                        source,
                        self.s3.as_ref(),
                        verify,
                        Arc::clone(cache),
                    )?
                } else {
                    BgzfReader::open(source, self.s3.as_ref(), verify)?
                };
                if self.observed {
                    reader.enable_timing();
                }
                let mut loaded = Vec::new();
                for (&(_, contig_id), contig_spans) in
                    spans.range((metagenome_id, 0)..=(metagenome_id, u32::MAX))
                {
                    let contig = self
                        .index
                        .contig(contig_id)?
                        .ok_or(TraceError::Invalid("missing source contig"))?;
                    let mut loaded_ranges = Vec::with_capacity(contig_spans.len());
                    for &(start, end) in contig_spans {
                        loaded_ranges.push(LoadedRange {
                            offset: start,
                            end,
                            sequence: reader.read_contig_range(contig, start, end)?,
                        });
                    }
                    loaded.push(((metagenome_id, contig_id), loaded_ranges));
                }
                if self.observed {
                    let mut stats = self.batch_stats.lock().unwrap();
                    stats.sequence_read_ns += reader.range_stats().read_nanoseconds.unwrap_or(0);
                    stats.bgzf_decode_and_handling_ns +=
                        reader.decompression_nanoseconds().unwrap_or(0);
                }
                Ok((
                    metagenome_id,
                    loaded,
                    (reader.range_stats(), reader.blocks_decoded()),
                ))
            })
            .collect::<Result<Vec<_>, TraceError>>()?;
        let mut loaded = BTreeMap::new();
        let mut reads = HashMap::new();
        for (metagenome_id, ranges, read) in per_metagenome {
            loaded.extend(ranges);
            reads.insert(metagenome_id, read);
        }
        Ok((loaded, reads))
    }

    fn align_tasks(
        &self,
        query: &[u8],
        tasks: &[AlignmentTask],
        loaded: &BTreeMap<(MetagenomeId, ContigId), Vec<LoadedRange>>,
        config: TraceConfig,
    ) -> Result<Vec<TaskOutcome>, TraceError> {
        let ledger = self.region_ledger.is_some();
        let mut geometry = [0usize; 4];
        for task in tasks {
            let query_bases = usize::try_from(task.query_span)
                .map_err(|_| TraceError::Invalid("query window"))?;
            let target_bases = task
                .target_end
                .checked_sub(task.target_start)
                .and_then(|span| usize::try_from(span).ok())
                .ok_or(TraceError::Invalid("loaded range"))?;
            let path_bases = query_bases
                .checked_add(target_bases)
                .ok_or(TraceError::Invalid("task window"))?;
            let cells = task_local_cells(task, target_bases, query.len(), config)?;
            for (maximum, value) in
                geometry
                    .iter_mut()
                    .zip([query_bases, target_bases, path_bases, cells])
            {
                *maximum = (*maximum).max(value);
            }
        }
        let [
            max_query_bases,
            max_target_bases,
            max_path_bases,
            max_local_cells,
        ] = geometry;
        let workspace_bytes = crate::alignment::trace_alignment_bytes(
            max_query_bases,
            max_target_bases,
            max_path_bases,
            max_local_cells,
            config.endpoint_bases,
            config.alignment,
        );
        if self.observed {
            let diagnostic = tasks
                .len()
                .checked_mul(128)
                .and_then(|bytes| bytes.checked_add(4096))
                .and_then(|bytes| CacheReservation::acquire(&LOOKUP_CACHE_AVAILABLE, bytes));
            // Query identity, topology, index generation and configuration are fixed for this call.
            let duplicates = diagnostic.as_ref().map(|_| {
                let unique: BTreeSet<_> = tasks
                    .iter()
                    .map(|t| {
                        (
                            t.metagenome_id,
                            t.contig_id,
                            t.strand,
                            t.query_start,
                            t.query_span,
                            t.target_start,
                            t.target_end,
                            t.diagonal_offset,
                        )
                    })
                    .collect();
                (tasks.len() - unique.len()) as u64
            });
            let mut stats = self.batch_stats.lock().unwrap();
            stats.alignment_tasks += tasks.len() as u64;
            if let Some(duplicates) = duplicates {
                stats.identical_alignment_tasks += duplicates;
            } else {
                stats.task_signature_unavailable_queries += 1;
            }
        }
        let workspace = || {
            if self.observed {
                self.batch_stats.lock().unwrap().alignment_workspaces += 1;
            }
            workspace_bytes.and_then(TraceAlignmentWorkspace::acquire)
        };
        tasks
            .par_iter()
            .map_init(workspace, |workspace, task| {
                let workspace = workspace
                    .as_mut()
                    .map_err(|error| TraceError::AlignmentAdmission(error.to_string()))?
                    .workspace_mut();
                if self.observed {
                    workspace.enable_timing();
                    workspace.work = Default::default();
                }
                let started = self.observed.then(Instant::now);
                let before = (
                    workspace.traceback_nanoseconds(),
                    workspace.endpoint_nanoseconds(),
                );
                let window_cpu = crate::alignment::observed_cpu(self.observed);
                let query_window =
                    linearize_query(query, task.query_start, task.query_span, config.circular)?;
                let loaded = loaded
                    .get(&(task.metagenome_id, task.contig_id))
                    .and_then(|ranges| {
                        ranges.iter().find(|range| {
                            range.offset <= task.target_start && range.end >= task.target_end
                        })
                    })
                    .ok_or(TraceError::Invalid("missing loaded range"))?;
                let start = usize::try_from(task.target_start - loaded.offset)
                    .map_err(|_| TraceError::Invalid("loaded range"))?;
                let end = usize::try_from(task.target_end - loaded.offset)
                    .map_err(|_| TraceError::Invalid("loaded range"))?;
                let target = loaded
                    .sequence
                    .get(start..end)
                    .ok_or(TraceError::Invalid("loaded range"))?;
                self.downstream_cpu(2, window_cpu);
                if self.observed {
                    let mut stats = self.batch_stats.lock().unwrap();
                    stats.alignment_query_bytes += query_window.len() as u64;
                    stats.alignment_target_bytes += target.len() as u64;
                    stats.alignment_query_max =
                        stats.alignment_query_max.max(query_window.len() as u64);
                    stats.alignment_target_max =
                        stats.alignment_target_max.max(target.len() as u64);
                    stats.copied_query_bytes += query_window.len() as u64;
                }
                let mut alignment_config = config.alignment;
                alignment_config.diagonal_offset = task.diagonal_offset;
                let margin = config.endpoint_bases as u64;
                let Some(initial) = align_task_window(
                    workspace,
                    &query_window,
                    target,
                    task.target_start,
                    task.strand,
                    alignment_config,
                    config,
                )?
                else {
                    self.record_alignment_time(started, workspace, before);
                    if self.observed {
                        self.batch_stats.lock().unwrap().rejected_alignments += 1;
                    }
                    return Ok(TaskOutcome {
                        fragment: None,
                        contact: false,
                        ledger: ledger.then(|| {
                            Box::new(TaskObservation {
                                work: workspace.work,
                                band_cells: task_local_cells(
                                    task,
                                    target.len(),
                                    query.len(),
                                    config,
                                )
                                .unwrap_or(usize::MAX)
                                    as u64,
                                core: None,
                                selected: None,
                                query_start: task.query_start,
                            })
                        }),
                    });
                };
                let mut contact =
                    crate::trace_islands::touches_inner_edge(task, &initial.core, margin);
                let mut alignment = initial.selected;
                let initial_accepted = alignment_accepted(&alignment, config);
                let mut query_start = task.query_start;
                let retry = if config.circular && task.query_span == query.len() as u64 {
                    let retry = circular_retry(
                        &initial.core,
                        task,
                        u64::try_from(query.len())
                            .map_err(|_| TraceError::Invalid("query length"))?,
                    )?;
                    if retry.is_some() {
                        retry
                    } else {
                        circular_retry(
                            &alignment,
                            task,
                            u64::try_from(query.len())
                                .map_err(|_| TraceError::Invalid("query length"))?,
                        )?
                    }
                } else {
                    None
                };
                if let Some((retry_start, retry_diagonal)) = retry {
                    let retry_cpu = crate::alignment::observed_cpu(self.observed);
                    if self.observed {
                        let mut stats = self.batch_stats.lock().unwrap();
                        stats.circular_retries += 1;
                        stats.copied_query_bytes += task.query_span;
                    }
                    let retry_query = linearize_query(query, retry_start, task.query_span, true)?;
                    alignment_config.diagonal_offset = retry_diagonal;
                    if let Some(retry) = align_task_window(
                        workspace,
                        &retry_query,
                        target,
                        task.target_start,
                        task.strand,
                        alignment_config,
                        config,
                    )? && alignment_accepted(&retry.selected, config)
                        && retry_improves(&alignment, &retry.selected, initial_accepted)
                    {
                        alignment = retry.selected;
                        query_start = retry_start;
                    }
                    self.downstream_cpu(3, retry_cpu);
                }
                contact |= crate::trace_islands::touches_inner_edge(task, &alignment, margin);
                let observation =
                    |workspace: &crate::alignment::AlignmentWorkspace,
                     alignment: &crate::alignment::Alignment| {
                        ledger.then(|| {
                            Box::new(TaskObservation {
                                work: workspace.work,
                                band_cells: task_local_cells(
                                    task,
                                    target.len(),
                                    query.len(),
                                    config,
                                )
                                .unwrap_or(usize::MAX)
                                    as u64,
                                core: Some(initial.core.clone()),
                                selected: Some(alignment.clone()),
                                query_start,
                            })
                        })
                    };
                if !alignment_accepted(&alignment, config) {
                    self.record_alignment_time(started, workspace, before);
                    if self.observed {
                        self.batch_stats.lock().unwrap().rejected_alignments += 1;
                    }
                    return Ok(TaskOutcome {
                        fragment: None,
                        contact,
                        ledger: observation(workspace, &alignment),
                    });
                }
                let projection_cpu = crate::alignment::observed_cpu(self.observed);
                let query_segments = query_segments(
                    query_start,
                    alignment.query_interval,
                    u64::try_from(query.len()).map_err(|_| TraceError::Invalid("query length"))?,
                    config.circular,
                )?;
                self.downstream_cpu(4, projection_cpu);
                self.record_alignment_time(started, workspace, before);
                if self.observed {
                    self.batch_stats.lock().unwrap().returned_alignments += 1;
                }
                let ledger = observation(workspace, &alignment);
                Ok(TaskOutcome {
                    fragment: Some((
                        task.metagenome_id,
                        Fragment {
                            contig_id: task.contig_id,
                            query_segments,
                            alignment,
                        },
                    )),
                    contact,
                    ledger,
                })
            })
            .collect()
    }
}

/// Diagnostic region support for one query: 21/31 context bits per core pair, each admitted
/// parent's anchor range and its ledger record.
#[derive(Default)]
struct RegionLedger {
    contexts: BTreeMap<(RegionKey, u64, u64), u8>,
    support: Vec<std::ops::Range<usize>>,
    parents: Vec<crate::trace_islands::ParentLedger>,
}

/// Result of one alignment task. `contact` is set only for island tasks.
struct TaskOutcome {
    fragment: Option<(MetagenomeId, Fragment)>,
    contact: bool,
    ledger: Option<Box<TaskObservation>>,
}

/// Per-task work and alignments retained for the region support ledger.
struct TaskObservation {
    work: crate::alignment::AlignmentWork,
    band_cells: u64,
    core: Option<crate::alignment::Alignment>,
    selected: Option<crate::alignment::Alignment>,
    query_start: u64,
}

struct WindowAlignment {
    core: crate::alignment::Alignment,
    selected: crate::alignment::Alignment,
}

fn alignment_accepted(alignment: &crate::alignment::Alignment, config: TraceConfig) -> bool {
    alignment.identity() >= config.min_identity
        && alignment.query_interval.len() >= config.min_aligned_bases
}

fn retry_improves(
    initial: &crate::alignment::Alignment,
    retry: &crate::alignment::Alignment,
    initial_accepted: bool,
) -> bool {
    !initial_accepted
        || (retry.query_interval.len() >= initial.query_interval.len()
            && retry.score >= initial.score
            && (retry.query_interval.len() > initial.query_interval.len()
                || retry.score > initial.score))
}

fn align_task_window(
    workspace: &mut crate::alignment::AlignmentWorkspace,
    query: &[u8],
    target: &[u8],
    target_start: u64,
    strand: Strand,
    alignment_config: AlignmentConfig,
    config: TraceConfig,
) -> Result<Option<WindowAlignment>, TraceError> {
    let core = match workspace.align_oriented(query, target, target_start, strand, alignment_config)
    {
        Ok(alignment) => alignment,
        Err(AlignmentError::NoAlignment) => return Ok(None),
        Err(error) => return Err(error.into()),
    };
    let completion = workspace.complete_endpoints(
        core.clone(),
        query,
        target,
        target_start,
        config.endpoint_bases,
        alignment_config,
    )?;
    #[cfg(feature = "bench-internals")]
    if workspace.work.local_passes > 0 {
        retain_alignment_fixture(
            query,
            target,
            target_start,
            alignment_config,
            config,
            &core,
            &completion,
            workspace.work,
        )?;
    }
    let completed = completion.alignment;
    let selected = if completed.identity() >= config.min_identity {
        completed
    } else {
        core.clone()
    };
    Ok(Some(WindowAlignment { core, selected }))
}

#[cfg(feature = "bench-internals")]
#[allow(clippy::too_many_arguments)]
fn retain_alignment_fixture(
    query: &[u8],
    target: &[u8],
    target_start: u64,
    alignment: AlignmentConfig,
    config: TraceConfig,
    core: &crate::alignment::Alignment,
    completed: &crate::alignment::EndpointCompletion,
    work: crate::alignment::AlignmentWork,
) -> Result<(), TraceError> {
    static ROOT: std::sync::OnceLock<Option<std::path::PathBuf>> = std::sync::OnceLock::new();
    static COUNT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let Some(root) =
        ROOT.get_or_init(|| std::env::var_os("JAM_ALIGNMENT_FIXTURES").map(Into::into))
    else {
        return Ok(());
    };
    let ordinal = COUNT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    if ordinal < 256 {
        let selected = if completed.alignment.identity() >= config.min_identity {
            &completed.alignment
        } else {
            core
        };
        let file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(root.join(format!("local-task-{ordinal:03}.json")))?;
        serde_json::to_writer(&file, &serde_json::json!({
            "query_length": query.len(), "target_length": target.len(),
            "strand": core.strand, "band": alignment.band_width,
            "diagonal": alignment.diagonal_offset,
            "scoring": [alignment.match_score, alignment.mismatch_score, alignment.gap_open_score, alignment.gap_extend_score],
            "accepted": alignment_accepted(selected, config), "work": work,
        })).map_err(|_| TraceError::Invalid("local task summary"))?;
        file.sync_all()?;
    }
    if ordinal >= 16 {
        return Ok(());
    }
    let bytes = query
        .len()
        .checked_add(target.len())
        .and_then(|bases| bases.checked_mul(64))
        .and_then(|bytes| bytes.checked_add(65536))
        .ok_or(TraceError::Invalid("alignment fixture budget"))?;
    let _reservation = CacheReservation::acquire(&LOOKUP_CACHE_AVAILABLE, bytes)
        .ok_or(TraceError::Invalid("alignment fixture budget"))?;
    let path = root.join(format!("task-{ordinal:02}.json"));
    let file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)?;
    serde_json::to_writer(&file, &serde_json::json!({
        "query": query, "target": target, "target_start": target_start,
        "scoring": [alignment.match_score, alignment.mismatch_score, alignment.gap_open_score, alignment.gap_extend_score],
        "band_width": alignment.band_width, "diagonal_offset": alignment.diagonal_offset,
        "max_cells": alignment.max_cells, "endpoint_bases": config.endpoint_bases,
        "circular": config.circular, "min_identity": config.min_identity,
        "min_aligned_bases": config.min_aligned_bases, "core": core,
        "completed": completed.alignment,
        "metrics": [completed.metrics.left_query_bases, completed.metrics.left_target_bases,
                    completed.metrics.right_query_bases, completed.metrics.right_target_bases,
                    completed.metrics.matrix_cells],
    })).map_err(|_| TraceError::Invalid("alignment fixture serialization"))?;
    file.sync_all()?;
    Ok(())
}

fn circular_retry(
    alignment: &crate::alignment::Alignment,
    task: &AlignmentTask,
    query_length: u64,
) -> Result<Option<(u64, i64)>, TraceError> {
    let suffix = alignment.query_interval.end == query_length && alignment.query_interval.start > 0;
    let prefix = alignment.query_interval.start == 0 && alignment.query_interval.end < query_length;
    let continuation = match (suffix, prefix, task.strand) {
        (true, false, Strand::Forward) => alignment.target_interval.end < task.target_end,
        (true, false, Strand::Reverse) => alignment.target_interval.start > task.target_start,
        (false, true, Strand::Forward) => alignment.target_interval.start > task.target_start,
        (false, true, Strand::Reverse) => alignment.target_interval.end < task.target_end,
        _ => false,
    };
    if !continuation {
        return Ok(None);
    }
    let boundary = if suffix {
        alignment.query_interval.start
    } else {
        alignment.query_interval.end
    };
    let retry_start = (u128::from(task.query_start) + u128::from(boundary))
        .checked_rem(u128::from(query_length))
        .and_then(|value| u64::try_from(value).ok())
        .ok_or(TraceError::Invalid("circular retry position"))?;
    let retry_query_offset = if suffix { 0 } else { query_length - boundary };
    let oriented_target_start = match task.strand {
        Strand::Forward => alignment
            .target_interval
            .start
            .checked_sub(task.target_start),
        Strand::Reverse => task.target_end.checked_sub(alignment.target_interval.end),
    }
    .ok_or(TraceError::Invalid("circular retry target"))?;
    let diagonal =
        i64::try_from(i128::from(oriented_target_start) - i128::from(retry_query_offset))
            .map_err(|_| TraceError::Invalid("circular retry diagonal"))?;
    Ok(Some((retry_start, diagonal)))
}

type LoadedRanges = (
    BTreeMap<(MetagenomeId, ContigId), Vec<LoadedRange>>,
    HashMap<MetagenomeId, (crate::range_source::RangeStats, u64)>,
);

struct LoadedRange {
    offset: u64,
    end: u64,
    sequence: Vec<u8>,
}

pub(crate) struct Candidate {
    pub(crate) id: MetagenomeId,
    pub(crate) name: String,
    pub(crate) shared_hashes: u32,
    pub(crate) containment: f64,
    pub(crate) exact_seed_hits: u64,
}

pub(crate) fn compare_candidates(left: &Candidate, right: &Candidate) -> std::cmp::Ordering {
    right
        .exact_seed_hits
        .cmp(&left.exact_seed_hits)
        .then_with(|| right.containment.total_cmp(&left.containment))
        .then_with(|| right.shared_hashes.cmp(&left.shared_hashes))
        .then_with(|| left.name.cmp(&right.name))
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct QuerySeed {
    packed_key: u64,
    position: u64,
    canonical_orientation: bool,
}

impl QuerySeed {
    fn token(self) -> Option<u64> {
        let position = u32::try_from(self.position).ok()?;
        (self.packed_key < 1 << 30).then_some(
            (self.packed_key << 33)
                | (u64::from(position) << 1)
                | u64::from(self.canonical_orientation),
        )
    }

    fn from_token(token: u64) -> Self {
        Self {
            packed_key: token >> 33,
            position: u64::from((token >> 1) as u32),
            canonical_orientation: token & 1 != 0,
        }
    }
}

struct QueryKeyRange {
    key: u64,
    first: usize,
    count: usize,
}

struct QueryPositions {
    core: Vec<QuerySeed>,
    tokens: Option<Vec<u64>>,
    core_occurrences: usize,
    core_distinct: usize,
    occurrence_probes: crate::shared_reader::CoreRequestCounts,
    nested: Vec<QuerySeed>,
    directory: Vec<QueryKeyRange>,
}

impl QueryPositions {
    fn new(mut core: Vec<QuerySeed>) -> Self {
        core.sort_unstable_by_key(|seed| {
            (seed.packed_key, seed.position, seed.canonical_orientation)
        });
        let mut first = 0;
        let directory: Vec<_> = core
            .chunk_by(|left, right| left.packed_key == right.packed_key)
            .map(|group| {
                let entry = QueryKeyRange {
                    key: group[0].packed_key,
                    first,
                    count: group.len(),
                };
                first += group.len();
                entry
            })
            .collect();
        Self {
            core_occurrences: core.len(),
            core_distinct: directory.len(),
            occurrence_probes: Default::default(),
            core,
            tokens: None,
            nested: Vec::new(),
            directory,
        }
    }

    #[cfg(test)]
    fn compact(query: &[u8], circular: bool) -> Option<Self> {
        Self::compact_screened(query, circular, None, false)
    }

    fn compact_screened(
        query: &[u8],
        circular: bool,
        operation: Option<&crate::shared_reader::SharedCoreOperation<'_>>,
        observed: bool,
    ) -> Option<Self> {
        u32::try_from(query.len().saturating_sub(1)).ok()?;
        let mut tokens = Vec::new();
        tokens.try_reserve_exact(query.len()).ok()?;
        if tokens.capacity() > query.len() {
            return None;
        }
        let mut original_keys = Vec::<u32>::new();
        let _diagnostic_reservation = if observed && operation.is_some() {
            let bytes = query.len().checked_mul(std::mem::size_of::<u32>())?;
            let reservation = CacheReservation::acquire(&LOOKUP_CACHE_AVAILABLE, bytes)?;
            original_keys.try_reserve_exact(query.len()).ok()?;
            if original_keys.capacity() > query.len() {
                return None;
            }
            Some(reservation)
        } else {
            None
        };
        let mut occurrences = 0;
        let mut probes = crate::shared_reader::CoreRequestCounts::default();
        let mut retain = |core: u32| {
            occurrences += 1;
            let Some(operation) = operation else {
                return Some(true);
            };
            if observed {
                original_keys.push(core);
            }
            let (covered, keep) = operation.screen(core).ok()?;
            probes.attempted += 1;
            probes.covered += usize::from(covered);
            probes.uncovered += usize::from(!covered);
            probes.rejected += usize::from(!keep);
            probes.retained += usize::from(keep);
            Some(keep)
        };
        if query.len() >= 15 {
            for (position, key, orientation) in query.bit_kmers(15, true) {
                if !retain(key.0 as u32)? {
                    continue;
                }
                tokens.push(
                    QuerySeed {
                        packed_key: key.0,
                        position: position as u64,
                        canonical_orientation: orientation,
                    }
                    .token()?,
                );
            }
            if circular {
                let start = query.len() - 14;
                let mut boundary = [0; 28];
                boundary[..14].copy_from_slice(&query[start..]);
                boundary[14..].copy_from_slice(&query[..14]);
                for (position, key, orientation) in boundary.bit_kmers(15, true) {
                    if !retain(key.0 as u32)? {
                        continue;
                    }
                    tokens.push(
                        QuerySeed {
                            packed_key: key.0,
                            position: (start + position) as u64,
                            canonical_orientation: orientation,
                        }
                        .token()?,
                    );
                }
            }
        }
        tokens.sort_unstable();
        let mut first = 0;
        let directory: Vec<_> = tokens
            .chunk_by(|left, right| left >> 33 == right >> 33)
            .map(|group| {
                let entry = QueryKeyRange {
                    key: group[0] >> 33,
                    first,
                    count: group.len(),
                };
                first += group.len();
                entry
            })
            .collect();
        let core_distinct = if observed && operation.is_some() {
            original_keys.sort_unstable();
            original_keys.dedup();
            original_keys.len()
        } else if operation.is_none() {
            directory.len()
        } else {
            0
        };
        Some(Self {
            core: Vec::new(),
            core_occurrences: occurrences,
            core_distinct,
            occurrence_probes: probes,
            tokens: Some(tokens),
            nested: Vec::new(),
            directory,
        })
    }

    fn materialize_tokens(&mut self) -> Result<(), TraceError> {
        let Some(tokens) = self.tokens.take() else {
            return Ok(());
        };
        let count = self.directory.iter().map(|entry| entry.count).sum();
        self.core
            .try_reserve_exact(count)
            .map_err(|_| TraceError::Invalid("query core associations"))?;
        if self.core.capacity() > count {
            return Err(TraceError::Invalid("query core association capacity"));
        }
        for entry in &mut self.directory {
            let first = self.core.len();
            self.core.extend(
                tokens[entry.first..entry.first + entry.count]
                    .iter()
                    .copied()
                    .map(QuerySeed::from_token),
            );
            entry.first = first;
        }
        Ok(())
    }

    fn add_nested(&mut self, mut nested: Vec<QuerySeed>, query_length: u64) {
        nested.sort_unstable_by_key(|seed| {
            let flank = if seed.packed_key >> 62 == 1 { 3 } else { 8 };
            let position = seed.position + flank;
            let core_position = if position >= query_length {
                position - query_length
            } else {
                position
            };
            (seed.packed_key, core_position, seed.canonical_orientation)
        });
        let mut first = self.core.len();
        self.directory.extend(
            nested
                .chunk_by(|left, right| left.packed_key == right.packed_key)
                .map(|group| {
                    let entry = QueryKeyRange {
                        key: group[0].packed_key,
                        first,
                        count: group.len(),
                    };
                    first += group.len();
                    entry
                }),
        );
        self.nested = nested;
    }

    fn positions(&self, range: &QueryKeyRange) -> &[QuerySeed] {
        if range.first < self.core.len() {
            &self.core[range.first..range.first + range.count]
        } else {
            let first = range.first - self.core.len();
            &self.nested[first..first + range.count]
        }
    }

    fn len(&self) -> usize {
        self.directory.len()
    }

    fn keys(&self) -> impl Iterator<Item = &u64> {
        self.directory.iter().map(|entry| &entry.key)
    }

    fn values(&self) -> impl Iterator<Item = &[QuerySeed]> {
        self.directory.iter().map(|entry| self.positions(entry))
    }

    fn iter(&self) -> impl Iterator<Item = (&u64, &[QuerySeed])> {
        self.directory
            .iter()
            .map(|entry| (&entry.key, self.positions(entry)))
    }

    fn get(&self, key: &u64) -> Option<&[QuerySeed]> {
        self.directory
            .binary_search_by_key(key, |entry| entry.key)
            .ok()
            .map(|ordinal| self.positions(&self.directory[ordinal]))
    }

    fn contains_key(&self, key: &u64) -> bool {
        self.get(key).is_some()
    }
}

impl std::ops::Index<&u64> for QueryPositions {
    type Output = [QuerySeed];

    fn index(&self, key: &u64) -> &Self::Output {
        self.get(key).expect("query key")
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct RegionKey {
    pub(crate) metagenome_id: MetagenomeId,
    pub(crate) contig_id: ContigId,
    pub(crate) strand: Strand,
    pub(crate) k: u8,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct SeedHit {
    pub(crate) query: u64,
    pub(crate) target: u64,
    pub(crate) diagonal: i128,
}

#[derive(Default)]
enum RegionHits {
    #[default]
    Empty,
    One(SeedHit),
    Many(Vec<SeedHit>),
    Unique(BTreeSet<(u64, u64)>),
}

impl RegionHits {
    fn push(&mut self, hit: SeedHit) {
        match self {
            Self::Empty => *self = Self::One(hit),
            Self::One(first) => {
                let first = *first;
                let mut hits = Vec::new();
                hits.extend([first, hit]);
                *self = Self::Many(hits);
            }
            Self::Many(hits) => hits.push(hit),
            Self::Unique(_) => unreachable!("shared anchors require bounded admission"),
        }
    }

    fn push_unique(&mut self, hit: SeedHit, reserved: &mut usize) -> Result<(), TraceError> {
        let pair = (hit.query, hit.target);
        let next = reserved
            .checked_add(512)
            .filter(|&bytes| bytes <= 64 * 1024 * 1024);
        if let Self::Unique(pairs) = self
            && let Some(next) = next
        {
            if pairs.insert(pair) {
                *reserved = next;
            }
            return Ok(());
        }
        if match self {
            Self::One(first) => (first.query, first.target) == pair,
            Self::Unique(pairs) => pairs.contains(&pair),
            _ => false,
        } {
            return Ok(());
        }
        // Covers sparse BTree nodes, the outer region node and later SeedHit Vec growth.
        let next = next.ok_or(TraceError::Invalid(
            "shared anchor workspace exceeds byte budget",
        ))?;
        match self {
            Self::Empty => *self = Self::One(hit),
            Self::One(first) => {
                *self = Self::Unique(BTreeSet::from([(first.query, first.target), pair]))
            }
            Self::Unique(pairs) => {
                pairs.insert(pair);
            }
            Self::Many(_) => return Err(TraceError::Invalid("mixed shared anchor storage")),
        }
        *reserved = next;
        Ok(())
    }
}

pub(crate) struct RegionAccumulator {
    pub(crate) query_start: u64,
    pub(crate) query_end: u64,
    pub(crate) target_start: u64,
    pub(crate) target_end: u64,
    pub(crate) diagonal_min: i128,
    pub(crate) diagonal_max: i128,
    pub(crate) hits: u32,
    /// Ordered anchors of this region in the query's support array, when retained.
    pub(crate) support: std::ops::Range<usize>,
}

impl RegionAccumulator {
    pub(crate) fn new(hit: SeedHit) -> Self {
        Self {
            query_start: hit.query,
            query_end: hit.query,
            target_start: hit.target,
            target_end: hit.target,
            diagonal_min: hit.diagonal,
            diagonal_max: hit.diagonal,
            hits: 1,
            support: 0..0,
        }
    }

    fn accepts(&self, hit: SeedHit, max_diagonal_drift: u64) -> bool {
        hit.query > self.query_end
            && hit.target > self.target_end
            && self.diagonal_min.min(hit.diagonal) + i128::from(max_diagonal_drift)
                >= self.diagonal_max.max(hit.diagonal)
    }

    pub(crate) fn add(&mut self, hit: SeedHit) {
        self.query_end = hit.query;
        self.target_end = hit.target;
        self.diagonal_min = self.diagonal_min.min(hit.diagonal);
        self.diagonal_max = self.diagonal_max.max(hit.diagonal);
        self.hits = self.hits.saturating_add(1);
    }
}

fn form_regions(
    hits_by_contig: BTreeMap<RegionKey, RegionHits>,
    max_diagonal_drift: u64,
) -> Vec<(RegionKey, RegionAccumulator)> {
    form_regions_observed(hits_by_contig, max_diagonal_drift, None, None)
}

/// Forms regions. With `support`, each region's anchors are appended to one ordered array and
/// the region keeps its range; anchors stay in query order within the range.
fn form_regions_observed(
    hits_by_contig: BTreeMap<RegionKey, RegionHits>,
    max_diagonal_drift: u64,
    mut predecessor_tests: Option<&mut u64>,
    mut support: Option<&mut Vec<SeedHit>>,
) -> Vec<(RegionKey, RegionAccumulator)> {
    let mut output = Vec::new();
    for (key, hits) in hits_by_contig {
        let ordered = matches!(hits, RegionHits::Unique(_));
        let mut hits = match hits {
            RegionHits::Empty => continue,
            RegionHits::One(hit) => {
                output.push((key, RegionAccumulator::new(hit)));
                continue;
            }
            RegionHits::Many(hits) => hits,
            RegionHits::Unique(pairs) => pairs
                .into_iter()
                .map(|(query, target)| SeedHit {
                    query,
                    target,
                    diagonal: i128::from(target) - i128::from(query),
                })
                .collect(),
        };
        if !ordered {
            hits.sort_unstable_by_key(|hit| (hit.query, hit.target));
        }
        let mut regions = Vec::<RegionAccumulator>::new();
        let mut owners = Vec::new();
        for &hit in &hits {
            if let Some((ordinal, region)) =
                regions.iter_mut().enumerate().rev().find(|(_, region)| {
                    if let Some(tests) = predecessor_tests.as_deref_mut() {
                        *tests += 1;
                    }
                    region.accepts(hit, max_diagonal_drift)
                })
            {
                region.add(hit);
                if support.is_some() {
                    owners.push(ordinal);
                }
            } else {
                if support.is_some() {
                    owners.push(regions.len());
                }
                regions.push(RegionAccumulator::new(hit));
            }
        }
        if let Some(support) = support.as_deref_mut() {
            let mut next = vec![0usize; regions.len()];
            for &owner in &owners {
                next[owner] += 1;
            }
            let mut offset = support.len();
            for (region, next) in regions.iter_mut().zip(&mut next) {
                region.support = offset..offset + *next;
                *next = offset;
                offset = region.support.end;
            }
            support.resize(offset, SeedHit::default());
            for (&hit, &owner) in hits.iter().zip(&owners) {
                support[next[owner]] = hit;
                next[owner] += 1;
            }
        }
        output.extend(regions.into_iter().map(|region| (key, region)));
    }
    output
}

#[derive(Clone)]
pub(crate) struct AlignmentTask {
    metagenome_id: MetagenomeId,
    contig_id: ContigId,
    strand: Strand,
    pub(crate) query_start: u64,
    pub(crate) query_span: u64,
    pub(crate) target_start: u64,
    pub(crate) target_end: u64,
    diagonal_offset: i64,
    /// Ordinal of the admitted parent region within this query.
    parent: u32,
    /// Inner island sides; `None` for a full parent envelope.
    pub(crate) island: Option<crate::trace_islands::InnerEdges>,
}

pub(crate) const SHORT_CONTIG_ENVELOPE_BYTES: u64 = 64 * 1024;
const LONG_CONTIG_ENVELOPE_FLANK_LIMIT: u64 = 16 * 1024;
const ALIGNMENT_CELL_RESERVE_DIVISOR: usize = 16;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct FragmentEnvelope {
    pub(crate) query_start: u64,
    pub(crate) query_span: u64,
    pub(crate) target_start: u64,
    pub(crate) target_end: u64,
    pub(crate) diagonal_offset: i64,
}

pub(crate) fn fragment_envelope(
    region: &RegionAccumulator,
    key: RegionKey,
    query_length: u64,
    contig_length: u64,
    config: TraceConfig,
) -> Result<FragmentEnvelope, TraceError> {
    let k = u64::from(key.k);
    let query_seed_end = region
        .query_end
        .checked_add(k)
        .ok_or(TraceError::Invalid("query seed interval"))?;
    let target_seed_end = region
        .target_end
        .checked_add(k)
        .ok_or(TraceError::Invalid("target seed interval"))?;
    if region.query_start >= query_length
        || (!config.circular && query_seed_end > query_length)
        || target_seed_end > contig_length
    {
        return Err(TraceError::Invalid("seed interval outside sequence"));
    }

    let chain_span = query_seed_end
        .saturating_sub(region.query_start)
        .max(target_seed_end.saturating_sub(region.target_start));
    let diagonal_spread = u64::try_from(region.diagonal_max - region.diagonal_min)
        .map_err(|_| TraceError::Invalid("region diagonal range"))?;
    let identity_allowance = if config.min_identity == 0.0 {
        chain_span
    } else {
        ((chain_span as f64 * (1.0 - config.min_identity) / config.min_identity).ceil() as u64)
            .min(LONG_CONTIG_ENVELOPE_FLANK_LIMIT)
    };
    let minimum_extension = config.flank_bases.saturating_mul(4);
    let extension = minimum_extension
        .max(
            chain_span
                .saturating_add(diagonal_spread)
                .saturating_add(identity_allowance),
        )
        .min(LONG_CONTIG_ENVELOPE_FLANK_LIMIT);

    let bounded_target = (
        region.target_start.saturating_sub(extension),
        target_seed_end.saturating_add(extension).min(contig_length),
    );
    let projection = |(start, end)| {
        projected_query_window(
            i128::from(start) - region.diagonal_max - i128::from(identity_allowance),
            i128::from(end) - region.diagonal_min + i128::from(identity_allowance),
            region.query_start,
            query_seed_end.saturating_sub(region.query_start),
            query_length,
            config.circular,
        )
    };
    let full_query = (contig_length <= SHORT_CONTIG_ENVELOPE_BYTES)
        .then(|| projection(bounded_target))
        .transpose()?;
    // The admitted geometry uses the same final diagonal that alignment receives.
    let full = match full_query {
        Some((query_start, query_span)) => {
            let diagonal = envelope_diagonal(region, query_start, 0, query_length, config)?;
            envelope_fits_workspace(query_span, contig_length, diagonal, config)?.then_some((
                0,
                contig_length,
                query_start,
                query_span,
                diagonal,
            ))
        }
        None => None,
    };
    let (oriented_start, oriented_end, query_start, query_span, diagonal_offset) =
        if let Some(full) = full {
            full
        } else {
            let (query_start, query_span) = projection(bounded_target)?;
            let target_span = bounded_target.1 - bounded_target.0;
            let diagonal =
                envelope_diagonal(region, query_start, bounded_target.0, query_length, config)?;
            // Longer local tasks run in resident chunks; only endpoint completion must fit.
            if !endpoint_fits_workspace(query_span, target_span, config)? {
                return Err(TraceError::Invalid(
                    "fragment envelope exceeds alignment workspace",
                ));
            }
            (
                bounded_target.0,
                bounded_target.1,
                query_start,
                query_span,
                diagonal,
            )
        };
    let (target_start, target_end) = match key.strand {
        Strand::Forward => (oriented_start, oriented_end),
        Strand::Reverse => (contig_length - oriented_end, contig_length - oriented_start),
    };
    Ok(FragmentEnvelope {
        query_start,
        query_span,
        target_start,
        target_end,
        diagonal_offset,
    })
}

/// Local diagonal of the region's first anchor inside an oriented task window.
fn envelope_diagonal(
    region: &RegionAccumulator,
    query_start: u64,
    oriented_start: u64,
    query_length: u64,
    config: TraceConfig,
) -> Result<i64, TraceError> {
    let target_relative = region
        .target_start
        .checked_sub(oriented_start)
        .ok_or(TraceError::Invalid("task target position"))?;
    let query_relative = if region.query_start >= query_start {
        region.query_start - query_start
    } else if config.circular {
        region
            .query_start
            .checked_add(query_length)
            .and_then(|position| position.checked_sub(query_start))
            .ok_or(TraceError::Invalid("task query position"))?
    } else {
        return Err(TraceError::Invalid("task query position"));
    };
    i64::try_from(i128::from(target_relative) - i128::from(query_relative))
        .map_err(|_| TraceError::Invalid("task diagonal"))
}

fn envelope_fits_workspace(
    query_span: u64,
    target_span: u64,
    diagonal_offset: i64,
    config: TraceConfig,
) -> Result<bool, TraceError> {
    let query = usize::try_from(query_span).map_err(|_| TraceError::Invalid("query window"))?;
    let target = usize::try_from(target_span).map_err(|_| TraceError::Invalid("target window"))?;
    let local_cells =
        crate::alignment::band_cells(query, target, diagonal_offset, config.alignment.band_width)?;
    let reserve = (config.alignment.max_cells / ALIGNMENT_CELL_RESERVE_DIVISOR).max(1);
    let usable = config
        .alignment
        .max_cells
        .checked_sub(reserve)
        .ok_or(TraceError::Invalid("alignment workspace reserve"))?;
    Ok(local_cells <= usable && endpoint_fits_workspace(query_span, target_span, config)?)
}

fn endpoint_fits_workspace(
    query_span: u64,
    target_span: u64,
    config: TraceConfig,
) -> Result<bool, TraceError> {
    let query = usize::try_from(query_span).map_err(|_| TraceError::Invalid("query window"))?;
    let target = usize::try_from(target_span).map_err(|_| TraceError::Invalid("target window"))?;
    let endpoint_cells = query
        .min(config.endpoint_bases)
        .checked_add(1)
        .and_then(|rows| rows.checked_mul(target.min(config.endpoint_bases).saturating_add(1)))
        .ok_or(TraceError::Invalid("endpoint workspace"))?;
    Ok(endpoint_cells <= config.alignment.max_cells)
}

/// Local cells a task can store on any diagonal it runs. A circular whole-query task may retry
/// on another diagonal, so it uses row and column bounds that hold for every diagonal.
fn task_local_cells(
    task: &AlignmentTask,
    target: usize,
    query_length: usize,
    config: TraceConfig,
) -> Result<usize, TraceError> {
    let query =
        usize::try_from(task.query_span).map_err(|_| TraceError::Invalid("query window"))?;
    if !(config.circular && query == query_length) {
        return Ok(crate::alignment::band_cells(
            query,
            target,
            task.diagonal_offset,
            config.alignment.band_width,
        )?);
    }
    let band = usize::try_from(config.alignment.band_width)
        .ok()
        .and_then(|band| band.checked_mul(2)?.checked_add(1))
        .ok_or(TraceError::Invalid("alignment band"))?;
    let rows = query
        .checked_add(1)
        .and_then(|rows| rows.checked_mul(target.saturating_add(1).min(band)));
    let columns = target
        .checked_add(1)
        .and_then(|columns| columns.checked_mul(query.saturating_add(1).min(band)));
    rows.zip(columns)
        .map(|(rows, columns)| rows.min(columns))
        .ok_or(TraceError::Invalid("task local cells"))
}

fn projected_query_window(
    low: i128,
    high: i128,
    anchor: u64,
    chain_span: u64,
    query_length: u64,
    circular: bool,
) -> Result<(u64, u64), TraceError> {
    if query_length == 0 || low >= high {
        return Err(TraceError::Invalid("projected query interval"));
    }
    if !circular {
        let start = low.clamp(0, i128::from(query_length));
        let end = high.clamp(0, i128::from(query_length));
        if start >= end || i128::from(anchor) < start || i128::from(anchor) >= end {
            return Err(TraceError::Invalid("projected query interval"));
        }
        return Ok((start as u64, (end - start) as u64));
    }
    let width = u64::try_from(high - low).unwrap_or(u64::MAX);
    if width >= query_length {
        let flank = (query_length - chain_span.min(query_length)) / 2;
        let start =
            (i128::from(anchor) - i128::from(flank)).rem_euclid(i128::from(query_length)) as u64;
        return Ok((start, query_length));
    }
    let length = i128::from(query_length);
    let start = u64::try_from(low.rem_euclid(length))
        .map_err(|_| TraceError::Invalid("projected query interval"))?;
    let anchor_offset = (i128::from(anchor) - i128::from(start)).rem_euclid(length);
    if anchor_offset >= i128::from(width) {
        return Err(TraceError::Invalid("projected query interval"));
    }
    Ok((start, width))
}

pub(crate) fn prepare_query(
    query_id: impl Into<String>,
    sequence: &[u8],
    config: TraceConfig,
    k: u8,
    rescue_k15: bool,
) -> Result<PreparedQuery, TraceError> {
    prepare_query_with_tokens(query_id, sequence, config, k, rescue_k15, false)
}

fn reserve_query_tokens(
    lengths: impl IntoIterator<Item = usize>,
    budget: usize,
) -> Option<CacheReservation<'static>> {
    let bytes = lengths
        .into_iter()
        .filter(|&length| u32::try_from(length.saturating_sub(1)).is_ok())
        .try_fold(0usize, |sum, length| {
            sum.checked_add(
                length
                    .checked_mul(std::mem::size_of::<u64>() + std::mem::size_of::<QuerySeed>())?,
            )
        })?;
    (bytes <= budget)
        .then(|| CacheReservation::acquire(&LOOKUP_CACHE_AVAILABLE, bytes))
        .flatten()
}

fn prepare_query_with_tokens(
    query_id: impl Into<String>,
    sequence: &[u8],
    config: TraceConfig,
    k: u8,
    rescue_k15: bool,
    compact: bool,
) -> Result<PreparedQuery, TraceError> {
    prepare_query_screened(
        query_id, sequence, config, k, rescue_k15, compact, None, false,
    )
}

#[allow(clippy::too_many_arguments)]
fn prepare_query_screened(
    query_id: impl Into<String>,
    sequence: &[u8],
    config: TraceConfig,
    k: u8,
    rescue_k15: bool,
    compact: bool,
    operation: Option<&crate::shared_reader::SharedCoreOperation<'_>>,
    observed: bool,
) -> Result<PreparedQuery, TraceError> {
    validate_config(config)?;
    let query_id = query_id.into();
    if query_id.is_empty()
        || query_id
            .bytes()
            .any(|byte| matches!(byte, 0 | b'\n' | b'\r'))
    {
        return Err(TraceError::Invalid("query ID"));
    }
    let query = sequence.normalize(false).into_owned();
    if query.is_empty() {
        return Err(TraceError::Invalid("empty query"));
    }
    let query_length =
        u64::try_from(query.len()).map_err(|_| TraceError::Invalid("query length"))?;
    let positions_by_key = if compact && k == 15 && !rescue_k15 {
        QueryPositions::compact_screened(&query, config.circular, operation, observed)
    } else {
        None
    };
    let positions_by_key = match positions_by_key {
        Some(positions) => positions,
        None => {
            let mut positions =
                QueryPositions::new(query_seeds(&query, k, rescue_k15, config.circular)?);
            if let Some(operation) = operation {
                let mut counts = crate::shared_reader::CoreRequestCounts::default();
                let mut error = None;
                positions
                    .directory
                    .retain(|entry| match operation.screen(entry.key as u32) {
                        Ok((covered, keep)) => {
                            counts.attempted += 1;
                            counts.covered += usize::from(covered);
                            counts.uncovered += usize::from(!covered);
                            counts.rejected += usize::from(!keep);
                            counts.retained += usize::from(keep);
                            keep
                        }
                        Err(failure) => {
                            error = Some(failure);
                            true
                        }
                    });
                if let Some(error) = error {
                    return Err(error.into());
                }
                operation.record(counts);
            }
            positions
        }
    };
    let mut lookup_identity = [0; 35];
    lookup_identity[..32].copy_from_slice(&sha256(&query));
    lookup_identity[32..].copy_from_slice(&[k, u8::from(rescue_k15), u8::from(config.circular)]);
    Ok(PreparedQuery {
        batch_ordinal: 0,
        shared_cores: None,
        query_id,
        query_length,
        query,
        positions_by_key,
        lookup_identity: sha256(&lookup_identity),
    })
}

fn query_seeds(
    query: &[u8],
    k: u8,
    rescue_k15: bool,
    circular: bool,
) -> Result<Vec<QuerySeed>, TraceError> {
    let mut seeds = extract_query_seeds(query, k, circular)?;
    if rescue_k15 {
        let mut rescue = extract_query_seeds(query, 15, circular)?;
        for seed in &mut rescue {
            seed.packed_key |= RESCUE_K15_TAG;
        }
        seeds.extend(rescue);
    }
    Ok(seeds)
}

fn possible_context_associations(query: &[u8], circular: bool) -> u64 {
    let mut run = 0usize;
    let mut prefix = 0usize;
    let mut count = 0u64;
    for (position, base) in query.iter().enumerate() {
        if matches!(base, b'A' | b'C' | b'G' | b'T') {
            run += 1;
            if run == position + 1 {
                prefix = run;
            }
            count += u64::from(run >= 21) + u64::from(run >= 31);
        } else {
            run = 0;
        }
    }
    if circular {
        for length in [21, 31] {
            if query.len() >= length {
                count += (prefix.min(length - 1) + run.min(length - 1)).saturating_sub(length - 1)
                    as u64;
            }
        }
    }
    count
}

fn extract_query_seeds(query: &[u8], k: u8, circular: bool) -> Result<Vec<QuerySeed>, TraceError> {
    if query.len() < usize::from(k) {
        return Ok(Vec::new());
    }
    let mut seeds = Vec::new();
    for (position, kmer, orientation) in query.bit_kmers(k, true) {
        let seed = QuerySeed {
            packed_key: kmer.0,
            position: u64::try_from(position)
                .map_err(|_| TraceError::Invalid("query seed position"))?,
            canonical_orientation: orientation,
        };
        seeds.push(seed);
    }
    if circular {
        let flank = usize::from(k) - 1;
        let start = query.len() - flank;
        let mut boundary = Vec::with_capacity(2 * flank);
        boundary.extend_from_slice(&query[start..]);
        boundary.extend_from_slice(&query[..flank]);
        for (position, kmer, orientation) in boundary.bit_kmers(k, true) {
            seeds.push(QuerySeed {
                packed_key: kmer.0,
                position: (start + position) as u64,
                canonical_orientation: orientation,
            });
        }
    }
    Ok(seeds)
}

fn minimum_region_hits(k: u8, configured: u32) -> u32 {
    if k == 15 {
        configured.max(3)
    } else {
        configured
    }
}

pub(crate) fn digest_hex(digest: [u8; 32]) -> String {
    let mut output = String::with_capacity(64);
    for byte in digest {
        write!(&mut output, "{byte:02x}").expect("writing to a string cannot fail");
    }
    output
}

pub(crate) fn candidate_completion(
    candidate_count: usize,
    max_metagenomes: usize,
) -> Result<SearchCompletion, TraceError> {
    let omitted = candidate_count.saturating_sub(max_metagenomes);
    if omitted == 0 {
        Ok(SearchCompletion::Complete)
    } else {
        Ok(SearchCompletion::CandidateBudgetExceeded {
            candidates_omitted: u32::try_from(omitted)
                .map_err(|_| TraceError::Invalid("candidate count"))?,
        })
    }
}

fn coalesce_spans(spans: &mut Vec<(u64, u64)>) {
    spans.sort_unstable();
    let mut merged = Vec::with_capacity(spans.len());
    for &(start, end) in spans.iter() {
        if let Some((_, previous_end)) = merged.last_mut()
            && start <= *previous_end
        {
            *previous_end = (*previous_end).max(end);
        } else {
            merged.push((start, end));
        }
    }
    *spans = merged;
}

fn linearize_query(
    query: &[u8],
    start: u64,
    span: u64,
    circular: bool,
) -> Result<Vec<u8>, TraceError> {
    let start = usize::try_from(start).map_err(|_| TraceError::Invalid("query window"))?;
    let span = usize::try_from(span).map_err(|_| TraceError::Invalid("query window"))?;
    if !circular {
        return query
            .get(start..start + span)
            .map(<[u8]>::to_vec)
            .ok_or(TraceError::Invalid("query window"));
    }
    Ok((0..span)
        .map(|offset| query[(start + offset) % query.len()])
        .collect())
}

fn query_segments(
    window_start: u64,
    local: Interval,
    query_length: u64,
    circular: bool,
) -> Result<Vec<Interval>, TraceError> {
    let start = window_start
        .checked_add(local.start)
        .ok_or(TraceError::Invalid("query coordinates"))?;
    let end = window_start
        .checked_add(local.end)
        .ok_or(TraceError::Invalid("query coordinates"))?;
    if !circular || end <= query_length {
        return Ok(vec![Interval::new(start, end)?]);
    }
    if start >= query_length {
        return Ok(vec![Interval::new(
            start - query_length,
            end - query_length,
        )?]);
    }
    let wrapped = end % query_length;
    let mut segments = vec![Interval::new(start, query_length)?];
    if wrapped != 0 {
        segments.push(Interval::new(0, wrapped)?);
    }
    Ok(segments)
}

fn validate_config(config: TraceConfig) -> Result<(), TraceError> {
    if !config.min_containment.is_finite()
        || !(0.0..=1.0).contains(&config.min_containment)
        || config.max_metagenomes == 0
        || config.min_seed_hits == 0
        || config.diagonal_bin_bases == 0
        || !config.min_identity.is_finite()
        || !(0.0..=1.0).contains(&config.min_identity)
        || config.min_aligned_bases == 0
    {
        return Err(TraceError::Invalid("trace configuration"));
    }
    Ok(())
}

#[derive(Debug, Error)]
pub enum TraceError {
    #[error("trace I/O failed: {0}")]
    Io(#[from] io::Error),
    #[error(transparent)]
    Database(#[from] ReaderError),
    #[error(transparent)]
    Query(#[from] QueryError),
    #[error(transparent)]
    Jidx(#[from] JidxReaderError),
    #[error(transparent)]
    JidxFormat(#[from] JidxError),
    #[error("owner index failed: {0}")]
    Owner(#[from] crate::owner_format::OwnerReaderError),
    #[error(transparent)]
    Shared(#[from] crate::shared_format::SharedError),
    #[error(transparent)]
    Bgzf(#[from] BgzfError),
    #[error(transparent)]
    Alignment(#[from] AlignmentError),
    #[error("alignment workspace admission failed: {0}")]
    AlignmentAdmission(String),
    #[error(transparent)]
    Mosaic(#[from] MosaicError),
    #[error("invalid trace input: {0}")]
    Invalid(&'static str),
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::{Alignment, AlignmentWorkspace};
    use crate::cli::handlers::{TraceArgs, TraceInput, handle_trace_command};
    use crate::jidx_builder::{JidxBuildConfig, build_local_jidx};
    use crate::writer::{BuildConfig, build};
    use noodles_bgzf::{self as bgzf, gzi};
    use std::io::Write;

    type HitLoopEndHook = (String, Box<dyn FnOnce() + Send>);

    /// Runs once after the hit loop of the named query and before the identity check that ends
    /// the phase, so a test can change the index file identity inside the phase.
    static HIT_LOOP_END_HOOK: Mutex<Option<HitLoopEndHook>> = Mutex::new(None);

    pub(super) fn run_hit_loop_end_hook(query_id: &str) {
        let mut hook = HIT_LOOP_END_HOOK
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        if hook.as_ref().is_some_and(|(id, _)| id == query_id) {
            let (_, run) = hook.take().unwrap();
            drop(hook);
            run();
        }
    }

    #[test]
    fn lookup_cache_reservations_share_and_release_the_byte_limit() {
        let available = AtomicUsize::new(10);
        let mut first = CacheReservation::acquire(&available, 6).unwrap();
        assert!(CacheReservation::acquire(&available, 5).is_none());
        assert_eq!(available.load(Ordering::Relaxed), 4);
        let second = CacheReservation::acquire(&available, 4).unwrap();
        assert_eq!(available.load(Ordering::Relaxed), 0);
        first.retain(2);
        assert_eq!(available.load(Ordering::Relaxed), 4);
        assert_eq!(first.bytes, 2);
        assert!(CacheReservation::acquire(&available, 5).is_none());
        drop(first);
        assert_eq!(available.load(Ordering::Relaxed), 6);
        drop(second);
        assert_eq!(available.load(Ordering::Relaxed), 10);
    }

    #[test]
    fn unique_anchor_admission_keeps_duplicates_at_the_byte_cap() {
        let mut hits = RegionHits::Unique(BTreeSet::from([(1, 2), (2, 3)]));
        let mut reserved = 64 * 1024 * 1024 - 512;
        let hit = |query| SeedHit {
            query,
            target: query + 1,
            diagonal: 1,
        };
        hits.push_unique(hit(1), &mut reserved).unwrap();
        assert_eq!(reserved, 64 * 1024 * 1024 - 512);
        hits.push_unique(hit(3), &mut reserved).unwrap();
        assert_eq!(reserved, 64 * 1024 * 1024);
        hits.push_unique(hit(1), &mut reserved).unwrap();
        assert!(hits.push_unique(hit(4), &mut reserved).is_err());
        let RegionHits::Unique(pairs) = hits else {
            panic!("unique anchors");
        };
        assert_eq!(pairs, BTreeSet::from([(1, 2), (2, 3), (3, 4)]));
    }

    #[test]
    fn compact_query_tokens_preserve_checked_boundaries_and_filtered_associations() {
        let mut seeds = Vec::new();
        for packed_key in [0, (1 << 30) - 1] {
            for position in [0, u64::from(u32::MAX)] {
                for canonical_orientation in [false, true] {
                    let seed = QuerySeed {
                        packed_key,
                        position,
                        canonical_orientation,
                    };
                    assert_eq!(QuerySeed::from_token(seed.token().unwrap()), seed);
                    seeds.push(seed);
                }
            }
        }
        let mut tokens = seeds
            .iter()
            .rev()
            .map(|seed| seed.token().unwrap())
            .collect::<Vec<_>>();
        tokens.sort_unstable();
        assert_eq!(
            tokens
                .into_iter()
                .map(QuerySeed::from_token)
                .collect::<Vec<_>>(),
            seeds
        );
        for seed in [
            QuerySeed {
                packed_key: 1 << 30,
                position: 0,
                canonical_orientation: false,
            },
            QuerySeed {
                packed_key: 1,
                position: u64::from(u32::MAX) + 1,
                canonical_orientation: true,
            },
        ] {
            assert!(seed.token().is_none());
            let mut wide = QueryPositions::new(vec![seed]);
            wide.materialize_tokens().unwrap();
            assert_eq!(wide.core, [seed]);
        }
        for query in [b"".as_slice(), b"ACGT", b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"] {
            let mut compact = QueryPositions::compact(query, true).unwrap();
            compact.materialize_tokens().unwrap();
            assert_eq!(compact.core_occurrences, 0);
            assert!(compact.core.is_empty() && compact.directory.is_empty());
        }
        assert!(
            prepare_query_with_tokens("empty", b"", TraceConfig::default(), 15, false, true)
                .is_err()
        );
        let query = sequence();
        let wide = QueryPositions::new(extract_query_seeds(query.as_bytes(), 15, true).unwrap());
        for parity in [0, 1, 2] {
            let mut compact = QueryPositions::compact(query.as_bytes(), true).unwrap();
            compact.directory.retain(|entry| entry.key % 2 == parity);
            compact.materialize_tokens().unwrap();
            assert_eq!(compact.core_occurrences, wide.core.len());
            assert_eq!(compact.core_distinct, wide.directory.len());
            let expected = wide
                .core
                .iter()
                .copied()
                .filter(|seed| seed.packed_key % 2 == parity)
                .collect::<Vec<_>>();
            assert_eq!(compact.core, expected);
            for (key, positions) in compact.iter() {
                assert_eq!(positions, &wide[key]);
            }
        }
    }

    #[test]
    fn flat_query_associations_match_vector_groups_and_circular_extension() {
        for length in [2_000, 64_000, 250_000] {
            for variant in 0..3 {
                let mut state = 7u64;
                let input = (0..length)
                    .map(|position| {
                        state ^= state << 13;
                        state ^= state >> 7;
                        state ^= state << 17;
                        match variant {
                            0 => b"ACGTTGCA"[position % 8],
                            1 => b"ACGT"[(state & 3) as usize],
                            _ if position % 97 < 20 => b'N',
                            _ => b"acgt"[(state & 3) as usize],
                        }
                    })
                    .collect::<Vec<_>>();
                for circular in [false, true] {
                    let config = TraceConfig {
                        circular,
                        ..TraceConfig::default()
                    };
                    let wide = prepare_query("flat", &input, config, 15, false).unwrap();
                    let mut prepared =
                        prepare_query_with_tokens("flat", &input, config, 15, false, true).unwrap();
                    assert!(prepared.positions_by_key.tokens.is_some());
                    prepared.positions_by_key.materialize_tokens().unwrap();
                    assert_eq!(prepared.positions_by_key.core, wide.positions_by_key.core);
                    drop(wide);
                    let mut extended = prepared.query.clone();
                    if circular {
                        extended.extend_from_slice(&prepared.query[..14]);
                    }
                    let mut seeds = extended
                        .bit_kmers(15, true)
                        .take_while(|(position, _, _)| *position < input.len())
                        .map(|(position, key, reverse)| QuerySeed {
                            packed_key: key.0,
                            position: position as u64,
                            canonical_orientation: reverse,
                        })
                        .collect::<Vec<_>>();
                    seeds.sort_by_key(|seed| seed.packed_key);
                    let mut reference = seeds
                        .chunk_by(|left, right| left.packed_key == right.packed_key)
                        .map(|group| (group[0].packed_key, group.to_vec()))
                        .collect::<BTreeMap<_, _>>();
                    let mut nested = Vec::new();
                    for seed in &seeds {
                        let context = crate::shared_seed::context_seed(
                            &prepared.query,
                            seed.position as usize,
                            seed.packed_key as u32,
                            seed.canonical_orientation,
                            circular,
                        )
                        .unwrap();
                        for length in [21, 31] {
                            if let Some(key) = context.key(length) {
                                let flank = u64::from((length - 15) / 2);
                                let association = QuerySeed {
                                    packed_key: key.packed().unwrap(),
                                    position: if circular {
                                        (seed.position + prepared.query_length - flank)
                                            % prepared.query_length
                                    } else {
                                        seed.position - flank
                                    },
                                    canonical_orientation: seed.canonical_orientation,
                                };
                                reference
                                    .entry(association.packed_key)
                                    .or_default()
                                    .push(association);
                                nested.push(association);
                            }
                        }
                    }
                    assert_eq!(
                        possible_context_associations(&prepared.query, circular),
                        nested.len() as u64
                    );
                    prepared
                        .positions_by_key
                        .add_nested(nested, prepared.query_length);
                    assert_eq!(prepared.positions_by_key.len(), reference.len());
                    for ((actual_key, actual), (expected_key, expected)) in
                        prepared.positions_by_key.iter().zip(&reference)
                    {
                        assert_eq!(actual_key, expected_key);
                        assert_eq!(actual.len(), expected.len());
                        assert!(
                            actual == expected,
                            "length={length} variant={variant} circular={circular} first_difference={:?}",
                            actual.iter().zip(expected).position(|(a, b)| a != b)
                        );
                    }
                    println!(
                        "flat length={length} variant={variant} circular={circular} query_capacity={} core_capacity={} nested_capacity={} directory_capacity={}",
                        prepared.query.capacity(),
                        prepared.positions_by_key.core.capacity(),
                        prepared.positions_by_key.nested.capacity(),
                        prepared.positions_by_key.directory.capacity()
                    );
                }
            }
        }
    }

    fn sequence() -> String {
        let mut state = 7u64;
        (0..128)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(state >> 62) as usize] as char
            })
            .collect()
    }

    #[test]
    fn trace_defaults_align_long_windows_with_explicit_caps_preserved() {
        let sequence = vec![b'A'; 64_000];
        let mut workspace = AlignmentWorkspace::default();
        let explicit = TraceConfig {
            alignment: AlignmentConfig::default(),
            ..TraceConfig::default()
        };
        assert!(matches!(
            workspace.align(&sequence, &sequence, explicit.alignment),
            Err(AlignmentError::MatrixTooLarge {
                cells,
                max_cells: 4_000_000,
            }) if cells > 4_000_000
        ));

        let config = TraceConfig::default();
        for (target, strand) in [
            (sequence.clone(), Strand::Forward),
            (vec![b'T'; sequence.len()], Strand::Reverse),
        ] {
            let alignment = workspace
                .align_oriented(&sequence, &target, 100, strand, config.alignment)
                .unwrap();
            alignment.validate_cigar().unwrap();
            assert_eq!(alignment.query_interval, Interval::new(0, 64_000).unwrap());
            assert_eq!(
                alignment.target_interval,
                Interval::new(100, 64_100).unwrap()
            );
            assert_eq!(alignment.strand, strand);
            assert_eq!(alignment.cigar, "64000=");
            assert_eq!(alignment.matches, 64_000);
            assert_eq!(alignment.identity(), 1.0);
        }
    }

    #[test]
    fn candidate_census_matches_scalar_across_batch_boundary() {
        let directory = tempfile::tempdir().unwrap();
        let mut state = 97u64;
        let sequence = (0..50_000)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1);
                b"ACGT"[(state >> 62) as usize] as char
            })
            .collect::<String>();
        let fasta = directory.path().join("sample.fa");
        std::fs::write(&fasta, format!(">sample\n{sequence}\n")).unwrap();
        let jam = directory.path().join("database.jam");
        build(
            &[fasta],
            &jam,
            &BuildConfig {
                kmer_size: 21,
                fscale: 1,
                singleton: true,
                memory: 1,
                ..BuildConfig::default()
            },
        )
        .unwrap();

        let bgzf_path = directory.path().join("sample.bgz");
        let fai_path = directory.path().join("sample.bgz.fai");
        let gzi_path = directory.path().join("sample.bgz.gzi");
        let mut raw = b">contig\n".to_vec();
        for line in sequence.as_bytes().chunks(80) {
            raw.extend_from_slice(line);
            raw.push(b'\n');
        }
        let mut writer = bgzf::io::Writer::new(File::create(&bgzf_path).unwrap());
        writer.write_all(&raw).unwrap();
        writer.finish().unwrap();
        std::fs::write(
            &fai_path,
            format!("contig\t{}\t8\t80\t81\n", sequence.len()),
        )
        .unwrap();
        gzi::fs::write(&gzi_path, &gzi::Index::default()).unwrap();
        let manifest = directory.path().join("manifest.json");
        std::fs::write(
            &manifest,
            format!(
                "{{\"metagenomes\":[{{\"name\":\"sample\",\"bgzf\":\"{}\",\"fai\":\"{}\",\"gzi\":\"{}\"}}]}}",
                bgzf_path.display(),
                fai_path.display(),
                gzi_path.display()
            ),
        )
        .unwrap();
        let jidx = directory.path().join("database.jidx");
        build_local_jidx(
            &jam,
            &manifest,
            &jidx,
            JidxBuildConfig {
                k: 21,
                minimizer_window: 1,
                rescue_k15: false,
            },
        )
        .unwrap();

        let engine = TraceEngine::open(&jam, &jidx, &manifest, None).unwrap();
        let config = TraceConfig {
            use_sketch: false,
            circular: false,
            ..TraceConfig::default()
        };
        let prepared = prepare_query("batch", sequence.as_bytes(), config, 21, false).unwrap();
        assert!(prepared.positions_by_key.len() > SEED_LOOKUP_BATCH_KEYS);
        let census = engine.candidate_census(&prepared, config, 0).unwrap();

        let mut scalar_frequencies = Vec::new();
        let mut scalar_hits = BTreeMap::<MetagenomeId, u64>::new();
        for (&packed_key, query_seeds) in prepared.positions_by_key.iter() {
            let Some(seed) = engine.index.find_seed(packed_key).unwrap() else {
                continue;
            };
            scalar_frequencies.push((packed_key, seed.document_frequency()));
            let query_positions = u64::try_from(query_seeds.len()).unwrap();
            for document in engine.index.seed_documents(seed).unwrap() {
                *scalar_hits.entry(document.metagenome_id()).or_default() +=
                    document.occurrence_count() * query_positions;
            }
        }

        assert_eq!(census.frequencies, scalar_frequencies);
        assert_eq!(census.candidates.len(), scalar_hits.len());
        assert!(census.candidates.iter().all(|candidate| {
            scalar_hits.get(&candidate.id).copied() == Some(candidate.exact_seed_hits)
        }));

        let packed_order = census
            .frequencies
            .iter()
            .map(|&(key, _)| key)
            .collect::<Vec<_>>();
        assert!(packed_order.len() > SEED_LOOKUP_BATCH_KEYS);
        let mut rarity_order = packed_order.clone();
        rarity_order.reverse();
        assert_ne!(rarity_order, packed_order);
        assert!(
            engine
                .trace_selected(&prepared, Vec::new(), &rarity_order, None, config)
                .unwrap()
                .is_empty()
        );
        assert!(
            engine
                .trace_selected(&prepared, Vec::new(), &packed_order, None, config)
                .unwrap()
                .is_empty()
        );

        let missing = (0..1u64 << 42)
            .find(|key| !prepared.positions_by_key.contains_key(key))
            .unwrap();
        let mut repeated_and_missing = rarity_order.clone();
        repeated_and_missing.insert(SEED_LOOKUP_BATCH_KEYS - 1, missing);
        repeated_and_missing.insert(SEED_LOOKUP_BATCH_KEYS, rarity_order[SEED_LOOKUP_BATCH_KEYS]);
        assert!(
            engine
                .trace_selected(&prepared, Vec::new(), &repeated_and_missing, None, config)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn traces_an_arbitrary_query_deterministically() {
        let directory = tempfile::tempdir().unwrap();
        let sequence = sequence();
        let fasta = directory.path().join("sample.fa");
        std::fs::write(&fasta, format!(">sample\n{sequence}\n")).unwrap();
        let fasta2 = directory.path().join("sample2.fa");
        std::fs::write(&fasta2, format!(">sample2\n{sequence}\n")).unwrap();
        let jam = directory.path().join("database.jam");
        build(
            &[fasta, fasta2],
            &jam,
            &BuildConfig {
                kmer_size: 5,
                fscale: 1,
                singleton: true,
                memory: 1,
                ..BuildConfig::default()
            },
        )
        .unwrap();

        let bgzf_path = directory.path().join("sample.bgz");
        let fai_path = directory.path().join("sample.bgz.fai");
        let gzi_path = directory.path().join("sample.bgz.gzi");
        let mut raw = b">contig\n".to_vec();
        for line in sequence.as_bytes().chunks(16) {
            raw.extend_from_slice(line);
            raw.push(b'\n');
        }
        let mut writer = bgzf::io::Writer::new(File::create(&bgzf_path).unwrap());
        writer.write_all(&raw).unwrap();
        writer.finish().unwrap();
        std::fs::write(&fai_path, b"contig\t128\t8\t16\t17\n").unwrap();
        gzi::fs::write(&gzi_path, &gzi::Index::default()).unwrap();
        let bgzf_path2 = directory.path().join("sample2.bgz");
        let fai_path2 = directory.path().join("sample2.bgz.fai");
        let gzi_path2 = directory.path().join("sample2.bgz.gzi");
        std::fs::copy(&bgzf_path, &bgzf_path2).unwrap();
        std::fs::copy(&fai_path, &fai_path2).unwrap();
        std::fs::copy(&gzi_path, &gzi_path2).unwrap();
        let manifest = directory.path().join("manifest.json");
        std::fs::write(
            &manifest,
            format!(
                "{{\"metagenomes\":[{{\"name\":\"sample\",\"bgzf\":\"{}\",\"fai\":\"{}\",\"gzi\":\"{}\"}},{{\"name\":\"sample2\",\"bgzf\":\"{}\",\"fai\":\"{}\",\"gzi\":\"{}\"}}]}}",
                bgzf_path.display(),
                fai_path.display(),
                gzi_path.display(),
                bgzf_path2.display(),
                fai_path2.display(),
                gzi_path2.display()
            ),
        )
        .unwrap();
        let jidx = directory.path().join("database.jidx");
        build_local_jidx(
            &jam,
            &manifest,
            &jidx,
            JidxBuildConfig {
                k: 5,
                minimizer_window: 4,
                rescue_k15: false,
            },
        )
        .unwrap();

        let engine = TraceEngine::open(&jam, &jidx, &manifest, None).unwrap();
        let wrong_manifest = directory.path().join("wrong-manifest.json");
        std::fs::write(&wrong_manifest, b"{}").unwrap();
        assert!(TraceEngine::open(&jam, &jidx, wrong_manifest, None).is_err());
        engine.verify_index().unwrap();
        let config = TraceConfig {
            min_seed_hits: 2,
            diagonal_bin_bases: 32,
            flank_bases: 32,
            min_aligned_bases: 32,
            endpoint_bases: 32,
            alignment: AlignmentConfig {
                band_width: 32,
                ..AlignmentConfig::default()
            },
            ..TraceConfig::default()
        };
        let direct_config = TraceConfig {
            use_sketch: false,
            ..config
        };
        let prepared =
            prepare_query("census", sequence.as_bytes(), direct_config, 5, false).unwrap();
        let census = engine
            .candidate_census(&prepared, direct_config, 0)
            .unwrap();
        let mut explicit_hits = HashMap::<MetagenomeId, u64>::new();
        for (&packed_key, query_seeds) in prepared.positions_by_key.iter() {
            let Some(seed) = engine.index.find_seed(packed_key).unwrap() else {
                continue;
            };
            for document in engine.index.seed_documents(seed).unwrap() {
                let occurrences = engine
                    .index
                    .seed_document_occurrences(seed, document)
                    .unwrap();
                *explicit_hits.entry(document.metagenome_id()).or_default() +=
                    u64::try_from(occurrences.len() * query_seeds.len()).unwrap();
            }
        }
        assert!(census.candidates.iter().all(|candidate| {
            explicit_hits.get(&candidate.id).copied() == Some(candidate.exact_seed_hits)
        }));
        let mut rare_first = census.frequencies;
        rare_first.sort_unstable_by_key(|&(key, frequency)| (frequency, key));
        let rare_first = rare_first
            .into_iter()
            .map(|(key, _)| key)
            .collect::<Vec<_>>();
        let rare_first_result = engine
            .trace_selected(
                &prepared,
                census.candidates,
                &rare_first,
                None,
                direct_config,
            )
            .unwrap();
        let reverse_census = engine
            .candidate_census(&prepared, direct_config, 0)
            .unwrap();
        let reverse_order = rare_first.iter().rev().copied().collect::<Vec<_>>();
        let reverse_result = engine
            .trace_selected(
                &prepared,
                reverse_census.candidates,
                &reverse_order,
                None,
                direct_config,
            )
            .unwrap();
        assert_eq!(rare_first_result, reverse_result);
        let cache_fixed = std::mem::size_of::<Option<CachedSeedLookups>>() + 4096;
        let cache_per_key = std::mem::size_of::<CachedDocumentGroup>()
            + usize::try_from(engine.index.document_count()).unwrap()
                * std::mem::size_of::<SeedDocument>();
        let one_group_cache_bytes = cache_fixed + cache_per_key;
        let full_cache_bytes = cache_fixed + cache_per_key * prepared.positions_by_key.len();
        for cache_bytes in [0, cache_fixed, one_group_cache_bytes, full_cache_bytes] {
            let mut cached = engine
                .candidate_census(&prepared, direct_config, cache_bytes)
                .unwrap();
            if let Some(lookups) = &cached.lookups {
                let allocated = cache_fixed
                    + lookups.groups.capacity() * std::mem::size_of::<CachedDocumentGroup>()
                    + lookups.documents.capacity() * std::mem::size_of::<SeedDocument>();
                assert!(allocated <= lookups._reservation.bytes);
                assert_eq!(
                    lookups.documents.len(),
                    lookups
                        .groups
                        .iter()
                        .map(|group| group.document_count)
                        .sum::<usize>()
                );
                if cache_bytes == one_group_cache_bytes {
                    assert_eq!(lookups.groups.len(), 1);
                    assert!(lookups.frozen);
                }
                if cache_bytes == full_cache_bytes {
                    assert!(
                        prepared
                            .positions_by_key
                            .keys()
                            .all(|&key| lookups.get(key).is_some())
                    );
                    assert!(
                        prepared
                            .positions_by_key
                            .keys()
                            .any(|&key| matches!(lookups.get(key), Some(None)))
                    );
                } else {
                    assert!(
                        prepared
                            .positions_by_key
                            .keys()
                            .any(|&key| lookups.get(key).is_none())
                    );
                }
                for &key in prepared.positions_by_key.keys() {
                    match lookups.get(key) {
                        Some(Some(cached)) => {
                            assert_eq!(Some(cached.seed), engine.index.find_seed(key).unwrap());
                            assert_eq!(
                                cached.documents,
                                engine.index.seed_documents(cached.seed).unwrap()
                            );
                        }
                        Some(None) => assert_eq!(None, engine.index.find_seed(key).unwrap()),
                        None => {}
                    }
                }
            } else {
                assert!(cache_bytes <= cache_fixed || cfg!(not(unix)));
            }
            let result = engine
                .trace_selected(
                    &prepared,
                    cached.candidates,
                    &reverse_order,
                    cached.lookups.as_ref(),
                    direct_config,
                )
                .unwrap();
            assert_eq!(rare_first_result, result);
            if let Some(lookups) = &mut cached.lookups {
                lookups.header_sha256[0] ^= 1;
                assert!(
                    engine
                        .trace_selected(
                            &prepared,
                            Vec::new(),
                            &reverse_order,
                            Some(lookups),
                            direct_config
                        )
                        .is_err()
                );
                lookups.header_sha256[0] ^= 1;
                lookups.query_identity[0] ^= 1;
                assert!(
                    engine
                        .trace_selected(
                            &prepared,
                            Vec::new(),
                            &reverse_order,
                            Some(lookups),
                            direct_config
                        )
                        .is_err()
                );
            }
        }
        let first = engine
            .search("plasmid", sequence.as_bytes(), config)
            .unwrap();
        let padded_query = [
            b"N".repeat(64),
            sequence.as_bytes().to_vec(),
            b"N".repeat(64),
        ]
        .concat();
        let padded_target = [
            b"A".repeat(64),
            sequence.as_bytes().to_vec(),
            b"T".repeat(64),
        ]
        .concat();
        let fragments = engine
            .align_tasks(
                &padded_query,
                &[AlignmentTask {
                    metagenome_id: 0,
                    contig_id: 0,
                    strand: Strand::Forward,
                    query_start: 0,
                    query_span: padded_query.len() as u64,
                    target_start: 0,
                    target_end: padded_target.len() as u64,
                    diagonal_offset: 0,
                    parent: 0,
                    island: None,
                }],
                &BTreeMap::from([(
                    (0, 0),
                    vec![LoadedRange {
                        offset: 0,
                        end: padded_target.len() as u64,
                        sequence: padded_target,
                    }],
                )]),
                TraceConfig {
                    endpoint_bases: 64,
                    circular: false,
                    ..config
                },
            )
            .unwrap();
        let fragments = fragments
            .into_iter()
            .filter_map(|outcome| outcome.fragment)
            .collect::<Vec<_>>();
        assert_eq!(fragments.len(), 1);
        assert_eq!(fragments[0].1.alignment.identity(), 1.0);
        assert_eq!(
            fragments[0].1.alignment.query_interval,
            Interval::new(64, 192).unwrap()
        );
        let second = engine
            .search("plasmid", sequence.as_bytes(), config)
            .unwrap();
        assert_eq!(first, second);
        assert_eq!(first.completion, SearchCompletion::Complete);
        assert_eq!(first.index.seed_k, 5);
        assert!(!first.index.rescue_k15);
        assert_eq!(first.index.manifest_sha256.len(), 64);
        assert_eq!(first.index.body_sha256.len(), 64);
        assert_eq!(first.metagenomes.len(), 2);
        assert_eq!(first.metagenomes[0].mosaic.covered_bases, 128);
        assert_eq!(first.metagenomes[0].contigs[0].name, "contig");
        assert!(first.metagenomes[0].bgzf_blocks_decoded > 0);
        serde_json::to_vec(&first).unwrap();

        let query = directory.path().join("query.fa");
        std::fs::write(&query, format!(">sample description\n{sequence}\n")).unwrap();
        let output = directory.path().join("trace.jsonl");
        handle_trace_command(TraceArgs {
            query: query.clone(),
            input: TraceInput::Shard {
                database: jam.clone(),
                index: jidx.clone(),
                manifest: manifest.clone(),
            },
            audit_index: false,
            output: output.clone(),
            query_id: None,
            config,
            s3: None,
            force: false,
        })
        .unwrap();
        let published = std::fs::read_to_string(output).unwrap();
        assert!(published.ends_with('\n'));
        assert_eq!(
            serde_json::from_str::<serde_json::Value>(&published).unwrap()["query_id"],
            "sample"
        );

        std::fs::write(&query, format!(">second\n{sequence}\n>first\n{sequence}\n")).unwrap();
        let batch_output = directory.path().join("batch.jsonl");
        handle_trace_command(TraceArgs {
            query,
            input: TraceInput::Shard {
                database: jam,
                index: jidx.clone(),
                manifest,
            },
            audit_index: false,
            output: batch_output.clone(),
            query_id: None,
            config,
            s3: None,
            force: false,
        })
        .unwrap();
        let mut batch = std::fs::read_to_string(batch_output)
            .unwrap()
            .lines()
            .map(|line| serde_json::from_str::<serde_json::Value>(line).unwrap())
            .collect::<Vec<_>>();
        let mut singles = ["second", "first"].map(|id| {
            serde_json::to_value(engine.search(id, sequence.as_bytes(), config).unwrap()).unwrap()
        });
        for result in batch.iter_mut().chain(singles.iter_mut()) {
            for trace in result["metagenomes"].as_array_mut().unwrap() {
                for field in [
                    "compressed_bytes_read",
                    "range_requests",
                    "bgzf_blocks_decoded",
                ] {
                    trace[field] = serde_json::json!(0);
                }
            }
        }
        assert_eq!(batch, singles);

        let direct = engine
            .search(
                "direct",
                sequence.as_bytes(),
                TraceConfig {
                    use_sketch: false,
                    ..config
                },
            )
            .unwrap();
        assert_eq!(direct.completion, SearchCompletion::Complete);
        assert_eq!(direct.metagenomes.len(), 2);
        assert_eq!(direct.metagenomes[0].shared_hashes, 0);
        assert_eq!(direct.metagenomes[0].mosaic.covered_bases, 128);

        let absent = (0..1 << 10)
            .map(|packed| {
                (0..5)
                    .rev()
                    .map(|shift| b"ACGT"[(packed >> (shift * 2)) & 3] as char)
                    .collect::<String>()
            })
            .find(|word| {
                let reverse_complement = word
                    .bytes()
                    .rev()
                    .map(|base| match base {
                        b'A' => 'T',
                        b'C' => 'G',
                        b'G' => 'C',
                        b'T' => 'A',
                        _ => unreachable!(),
                    })
                    .collect::<String>();
                !sequence.contains(word) && !sequence.contains(&reverse_complement)
            })
            .unwrap();
        let partial = format!("{}{}", &sequence[..64], absent.repeat(13));
        let sketch_miss = engine
            .search(
                "sketch-miss",
                partial.as_bytes(),
                TraceConfig {
                    min_containment: 1.0,
                    circular: false,
                    ..config
                },
            )
            .unwrap();
        assert_eq!(sketch_miss.metagenomes.len(), 2);
        assert!(
            sketch_miss
                .metagenomes
                .iter()
                .all(|metagenome| metagenome.shared_hashes == 0)
        );

        std::fs::remove_file(bgzf_path2).unwrap();
        let capped = engine
            .search(
                "capped",
                sequence.as_bytes(),
                TraceConfig {
                    max_metagenomes: 1,
                    ..config
                },
            )
            .unwrap();
        assert_eq!(capped.metagenomes.len(), 1);
        assert_eq!(capped.metagenomes[0].name, "sample");
        assert_eq!(
            capped.completion,
            SearchCompletion::CandidateBudgetExceeded {
                candidates_omitted: 1
            }
        );
        #[cfg(unix)]
        {
            let cached = engine
                .candidate_census(&prepared, direct_config, full_cache_bytes)
                .unwrap();
            let times = std::fs::FileTimes::new().set_modified(
                std::time::SystemTime::UNIX_EPOCH + std::time::Duration::from_secs(42),
            );
            std::fs::OpenOptions::new()
                .write(true)
                .open(jidx)
                .unwrap()
                .set_times(times)
                .unwrap();
            assert!(matches!(
                engine.trace_selected(
                    &prepared,
                    cached.candidates,
                    &rare_first,
                    cached.lookups.as_ref(),
                    direct_config,
                ),
                Err(TraceError::Invalid("cached seed lookup identity"))
            ));
        }
    }

    #[test]
    fn maps_fully_wrapped_query_intervals() {
        assert_eq!(
            query_segments(120, Interval::new(10, 20).unwrap(), 128, true).unwrap(),
            vec![Interval::new(2, 12).unwrap()]
        );
    }

    #[test]
    fn retains_every_repeated_query_seed_position() {
        let seeds = query_seeds(b"AAAAAA", 3, false, false).unwrap();
        assert_eq!(
            seeds.iter().map(|seed| seed.position).collect::<Vec<_>>(),
            vec![0, 1, 2, 3]
        );
        assert!(
            seeds
                .iter()
                .all(|seed| seed.packed_key == seeds[0].packed_key)
        );
    }

    #[test]
    fn region_hits_inline_first_promotes_in_order() {
        let first = SeedHit {
            query: 1,
            target: 2,
            diagonal: 1,
        };
        let second = SeedHit {
            query: 3,
            target: 5,
            diagonal: 2,
        };
        let third = SeedHit {
            query: 8,
            target: 13,
            diagonal: 5,
        };
        let mut hits = RegionHits::default();
        assert!(matches!(&hits, RegionHits::Empty));
        hits.push(first);
        assert!(matches!(&hits, RegionHits::One(hit) if *hit == first));
        hits.push(second);
        hits.push(third);
        assert!(matches!(
            &hits,
            RegionHits::Many(values) if values == &[first, second, third]
        ));
    }

    #[test]
    fn island_window_accepts_forty_base_eighty_percent_trace_without_exact_fifteen_mer() {
        let substitute = |base: u8| if base == b'A' { b'C' } else { b'A' };
        let core = window_dna(0x40b0_5eed, 40);
        let config = TraceConfig::default();
        let mut alignment_config = config.alignment;
        alignment_config.diagonal_offset = 100;
        let run = |mismatches: &[usize]| {
            let mut fragment = core.clone();
            for &at in mismatches {
                fragment[at] = substitute(fragment[at]);
            }
            let longest_exact = fragment
                .iter()
                .zip(&core)
                .fold((0, 0), |(run, best), (a, b)| {
                    let run = if a == b { run + 1 } else { 0 };
                    (run, best.max(run))
                })
                .1;
            assert!(longest_exact < 15);
            let left = window_dna(0x51, 900);
            let right = window_dna(0x52, 900);
            let mut target_left = window_dna(0x61, 1_000);
            let mut target_right = window_dna(0x62, 1_000);
            // Unrelated flanks start with forced mismatches, so the trace cannot grow by chance.
            for offset in 1..=3 {
                target_left[1_000 - offset] = substitute(left[900 - offset]);
                target_right[offset - 1] = substitute(right[offset - 1]);
            }
            let query = [left, fragment, right].concat();
            let target = [target_left, core.clone(), target_right].concat();
            let mut workspace = AlignmentWorkspace::default();
            let window = align_task_window(
                &mut workspace,
                &query,
                &target,
                0,
                Strand::Forward,
                alignment_config,
                config,
            )
            .unwrap()
            .unwrap();
            let task = AlignmentTask {
                metagenome_id: 0,
                contig_id: 0,
                strand: Strand::Forward,
                query_start: 0,
                query_span: query.len() as u64,
                target_start: 0,
                target_end: target.len() as u64,
                diagonal_offset: 100,
                parent: 0,
                island: Some(crate::trace_islands::InnerEdges {
                    query_left: true,
                    query_right: true,
                    target_left: true,
                    target_right: true,
                }),
            };
            let contact = crate::trace_islands::touches_inner_edge(
                &task,
                &window.selected,
                config.endpoint_bases as u64,
            );
            (window.selected, contact)
        };
        let eight: Vec<usize> = (2..40).step_by(5).collect();
        let (accepted, contact) = run(&eight);
        assert_eq!(accepted.query_interval, Interval::new(900, 940).unwrap());
        assert_eq!((accepted.matches, accepted.substitutions), (32, 8));
        assert!(alignment_accepted(&accepted, config) && !contact);
        let nine = [eight.as_slice(), &[20]].concat();
        let (rejected, _) = run(&nine);
        assert!(!alignment_accepted(&rejected, config));
    }

    #[test]
    fn islands_split_only_distant_anchor_groups_on_long_contigs() {
        let config = TraceConfig {
            circular: false,
            ..TraceConfig::default()
        };
        let key = envelope_key(Strand::Forward);
        let anchors = |positions: &[u64]| {
            positions
                .iter()
                .map(|&query| SeedHit {
                    query,
                    target: query + 20_000,
                    diagonal: 20_000,
                })
                .collect::<Vec<_>>()
        };
        let plan = |positions: &[u64], contig_length: u64| {
            let hits = anchors(positions);
            let mut region = RegionAccumulator::new(hits[0]);
            for &hit in &hits[1..] {
                region.add(hit);
            }
            region.support = 0..hits.len();
            let envelope = fragment_envelope(&region, key, 45_000, contig_length, config).unwrap();
            crate::trace_islands::plan_islands(
                &region,
                &envelope,
                &hits,
                key,
                45_000,
                contig_length,
                config,
            )
            .unwrap()
            .map(|islands| {
                islands
                    .iter()
                    .map(|island| (island.region.support.clone(), island.edges))
                    .collect::<Vec<_>>()
            })
        };
        assert_eq!(plan(&[5_000, 5_500, 6_000], 200_000), None);
        assert_eq!(plan(&[5_000, 35_000], 64 * 1024), None);
        let split = plan(&[5_000, 5_900, 35_000], 200_000).unwrap();
        assert_eq!(split.len(), 2);
        assert_eq!((split[0].0.clone(), split[1].0.clone()), (0..2, 2..3));
        let [first, second] = [split[0].1, split[1].1];
        assert!(first.query_left && first.target_left && first.query_right && first.target_right);
        assert!(
            second.query_left && second.target_left && second.query_right && second.target_right
        );
    }

    fn envelope_region(query: u64, target: u64, hits: u32) -> RegionAccumulator {
        RegionAccumulator {
            query_start: query,
            query_end: query + u64::from(hits.saturating_sub(1)) * 20,
            target_start: target,
            target_end: target + u64::from(hits.saturating_sub(1)) * 20,
            diagonal_min: i128::from(target) - i128::from(query),
            diagonal_max: i128::from(target) - i128::from(query),
            hits,
            support: 0..0,
        }
    }

    fn envelope_key(strand: Strand) -> RegionKey {
        RegionKey {
            metagenome_id: 0,
            contig_id: 0,
            strand,
            k: 15,
        }
    }

    #[test]
    fn short_contig_projection_keeps_prefix_before_late_anchor() {
        let config = TraceConfig {
            circular: false,
            ..TraceConfig::default()
        };
        for offset in [384, 700, 780] {
            let query = 400 + offset;
            let region = envelope_region(query, offset, 1);
            for strand in [Strand::Forward, Strand::Reverse] {
                let envelope =
                    fragment_envelope(&region, envelope_key(strand), 2_000, 800, config).unwrap();
                assert_eq!((envelope.target_start, envelope.target_end), (0, 800));
                assert!(envelope.query_start <= 400);
                assert!(envelope.query_start + envelope.query_span >= 1_200);
            }
        }
        for (offset, hits) in [(512, 14), (640, 7)] {
            let region = envelope_region(400 + offset, offset, hits);
            for strand in [Strand::Forward, Strand::Reverse] {
                let envelope =
                    fragment_envelope(&region, envelope_key(strand), 2_000, 800, config).unwrap();
                assert_eq!((envelope.target_start, envelope.target_end), (0, 800));
                assert!(envelope.query_start <= 400);
                assert!(envelope.query_start + envelope.query_span >= 1_200);
            }
        }
    }

    #[test]
    fn long_contig_envelope_is_bounded_and_keeps_displaced_fragment() {
        let config = TraceConfig {
            circular: false,
            ..TraceConfig::default()
        };
        for offset in [384, 700, 780] {
            let query = 400 + offset;
            let target = 50_000 + offset;
            let region = envelope_region(query, target, 1);
            let forward = fragment_envelope(
                &region,
                envelope_key(Strand::Forward),
                2_000,
                100_000,
                config,
            )
            .unwrap();
            assert!(forward.target_end - forward.target_start <= 2 * 1_024 + 15);
            assert!(forward.query_start <= 400);
            assert!(forward.query_start + forward.query_span >= 1_200);

            let oriented_target = 100_000 - target - 15;
            let reverse_region = envelope_region(query, oriented_target, 1);
            let reverse = fragment_envelope(
                &reverse_region,
                envelope_key(Strand::Reverse),
                2_000,
                100_000,
                config,
            )
            .unwrap();
            assert!(reverse.target_start <= target);
            assert!(reverse.target_end >= target + 15);
            assert!(reverse.target_end - reverse.target_start <= 2 * 1_024 + 15);
        }
    }

    #[test]
    fn envelope_clips_real_ends_and_wraps_only_circular_queries() {
        assert_eq!(
            projected_query_window(-400, 2_400, 0, 2_000, 2_000, true).unwrap(),
            (0, 2_000)
        );
        let linear = TraceConfig {
            circular: false,
            ..TraceConfig::default()
        };
        let at_start = envelope_region(5, 5, 3);
        let start = fragment_envelope(
            &at_start,
            envelope_key(Strand::Forward),
            2_000,
            100_000,
            linear,
        )
        .unwrap();
        assert_eq!((start.query_start, start.target_start), (0, 0));

        let crossing = envelope_region(1_950, 150, 1);
        let circular = fragment_envelope(
            &crossing,
            envelope_key(Strand::Forward),
            2_000,
            400,
            TraceConfig::default(),
        )
        .unwrap();
        assert!(circular.query_start > 1_700);
        assert!(circular.query_start + circular.query_span > 2_000);

        let absent_target_end = envelope_region(100, 390, 1);
        assert!(
            fragment_envelope(
                &absent_target_end,
                envelope_key(Strand::Forward),
                2_000,
                400,
                linear,
            )
            .is_err()
        );
        let absent_query_end = envelope_region(1_990, 100, 1);
        assert!(
            fragment_envelope(
                &absent_query_end,
                envelope_key(Strand::Forward),
                2_000,
                400,
                linear,
            )
            .is_err()
        );
    }

    #[test]
    fn envelope_respects_custom_band_and_workspace_reserve() {
        let region = envelope_region(50_000, 30_000, 1);
        let config = TraceConfig {
            circular: false,
            alignment: AlignmentConfig {
                band_width: 64,
                max_cells: 300_000,
                ..AlignmentConfig::default()
            },
            ..TraceConfig::default()
        };
        let envelope = fragment_envelope(
            &region,
            envelope_key(Strand::Forward),
            100_000,
            60_000,
            config,
        )
        .unwrap();
        assert_eq!((envelope.target_start, envelope.target_end), (0, 60_000));
        assert!(
            envelope_fits_workspace(
                envelope.query_span,
                envelope.target_end - envelope.target_start,
                envelope.diagonal_offset,
                config,
            )
            .unwrap()
        );
        // A smaller resident bound no longer rejects the task: the full contig window no longer
        // fits, so the bounded window is chosen and runs in resident chunks when needed.
        let smaller = TraceConfig {
            alignment: AlignmentConfig {
                max_cells: 250_000,
                ..config.alignment
            },
            ..config
        };
        let bounded = fragment_envelope(
            &region,
            envelope_key(Strand::Forward),
            100_000,
            60_000,
            smaller,
        )
        .unwrap();
        assert!(bounded.target_start > 0 && bounded.target_end < 60_000);
        assert!(
            !envelope_fits_workspace(
                envelope.query_span,
                60_000,
                envelope.diagonal_offset,
                smaller
            )
            .unwrap()
        );
        assert!(
            fragment_envelope(
                &region,
                envelope_key(Strand::Forward),
                100_000,
                60_000,
                TraceConfig {
                    endpoint_bases: 600,
                    alignment: AlignmentConfig {
                        max_cells: 300_000,
                        ..smaller.alignment
                    },
                    ..smaller
                },
            )
            .is_err()
        );
    }

    #[test]
    fn envelope_diagonal_places_first_anchor_on_band_center() {
        let linear = TraceConfig {
            circular: false,
            ..TraceConfig::default()
        };
        let cases = [
            (envelope_region(1_184, 700, 3), 2_000, 800, linear),
            (envelope_region(1_184, 50_784, 3), 2_000, 100_000, linear),
            (
                envelope_region(1_950, 150, 1),
                2_000,
                400,
                TraceConfig::default(),
            ),
            (
                envelope_region(40, 90_000, 2),
                2_000,
                100_000,
                TraceConfig::default(),
            ),
        ];
        for (region, query_length, contig_length, config) in cases {
            for strand in [Strand::Forward, Strand::Reverse] {
                let envelope = fragment_envelope(
                    &region,
                    envelope_key(strand),
                    query_length,
                    contig_length,
                    config,
                )
                .unwrap();
                let oriented_start = match strand {
                    Strand::Forward => envelope.target_start,
                    Strand::Reverse => contig_length - envelope.target_end,
                };
                let query_relative =
                    (region.query_start + query_length - envelope.query_start) % query_length;
                let target_relative = region.target_start - oriented_start;
                assert_eq!(
                    envelope.diagonal_offset,
                    target_relative as i64 - query_relative as i64
                );
                let query = envelope.query_span as usize;
                let target = (envelope.target_end - envelope.target_start) as usize;
                let row = crate::alignment::band_row(
                    query_relative as usize,
                    target,
                    envelope.diagonal_offset,
                    config.alignment.band_width,
                )
                .unwrap()
                .unwrap();
                assert!(row.0 <= target_relative as usize && target_relative as usize <= row.1);
                let mut workspace = AlignmentWorkspace::default();
                workspace.enable_timing();
                let sequence = vec![b'A'; query.max(target)];
                let _ = workspace.align(
                    &sequence[..query],
                    &sequence[..target],
                    AlignmentConfig {
                        diagonal_offset: envelope.diagonal_offset,
                        ..config.alignment
                    },
                );
                assert_eq!(
                    workspace.work.local_cells as usize,
                    crate::alignment::band_cells(
                        query,
                        target,
                        envelope.diagonal_offset,
                        config.alignment.band_width
                    )
                    .unwrap()
                );
            }
        }
    }

    fn exact_alignment_pairs(
        alignment: &Alignment,
        query: &[u8],
        target: &[u8],
    ) -> std::collections::BTreeSet<(u64, u64)> {
        alignment.validate_cigar().unwrap();
        let mut query_position = alignment.query_interval.start as usize;
        let physical_start = alignment.target_interval.start as usize;
        let physical_end = alignment.target_interval.end as usize;
        let mut target_sequence = target[physical_start..physical_end].to_vec();
        if alignment.strand == Strand::Reverse {
            target_sequence = window_reverse_complement(&target_sequence);
        }
        let mut target_position = 0usize;
        let mut pairs = std::collections::BTreeSet::new();
        for run in &alignment.edit_script {
            for _ in 0..run.length {
                match run.operation {
                    crate::alignment::EditOperation::Equal => {
                        assert_eq!(query[query_position], target_sequence[target_position]);
                        let physical = match alignment.strand {
                            Strand::Forward => physical_start + target_position,
                            Strand::Reverse => physical_end - 1 - target_position,
                        };
                        pairs.insert((query_position as u64, physical as u64));
                        query_position += 1;
                        target_position += 1;
                    }
                    crate::alignment::EditOperation::Substitution => {
                        assert_ne!(query[query_position], target_sequence[target_position]);
                        query_position += 1;
                        target_position += 1;
                    }
                    crate::alignment::EditOperation::Insertion => target_position += 1,
                    crate::alignment::EditOperation::Deletion => query_position += 1,
                }
            }
        }
        assert_eq!(query_position as u64, alignment.query_interval.end);
        assert_eq!(target_position, target_sequence.len());
        pairs
    }

    #[test]
    fn short_contig_query_projection_separates_neighboring_band_optima() {
        let first = window_dna(101, 400);
        let second = window_dna(103, 400);
        let mut query = window_dna(107, 2_400);
        query[300..700].copy_from_slice(&first);
        query[1_500..1_900].copy_from_slice(&second);
        for diagonal_distance in [64u64, 65] {
            let second_target_start = 1_700 + diagonal_distance as usize;
            let mut oriented_target = window_dna(109, 4_802);
            oriented_target[500..900].copy_from_slice(&first);
            oriented_target[second_target_start..second_target_start + 400]
                .copy_from_slice(&second);
            let region = envelope_region(315, 515, 18);
            for strand in [Strand::Forward, Strand::Reverse] {
                let target = if strand == Strand::Forward {
                    oriented_target.clone()
                } else {
                    window_reverse_complement(&oriented_target)
                };
                let config = TraceConfig {
                    circular: false,
                    endpoint_bases: 0,
                    ..TraceConfig::default()
                };
                let envelope = fragment_envelope(
                    &region,
                    envelope_key(strand),
                    query.len() as u64,
                    target.len() as u64,
                    config,
                )
                .unwrap();
                assert_eq!((envelope.target_start, envelope.target_end), (0, 4_802));
                assert!(envelope.query_span < query.len() as u64);
                let query_window =
                    linearize_query(&query, envelope.query_start, envelope.query_span, false)
                        .unwrap();
                let query_relative = 315 - envelope.query_start;
                let mut bounded_config = config.alignment;
                bounded_config.diagonal_offset = 515 - query_relative as i64;
                let bounded = AlignmentWorkspace::default()
                    .align_oriented(&query_window, &target, 0, strand, bounded_config)
                    .unwrap();
                let (first_start, first_end) = match strand {
                    Strand::Forward => (500, 900),
                    Strand::Reverse => (3_902, 4_302),
                };
                assert!(
                    bounded.target_interval.start <= first_start
                        && bounded.target_interval.end >= first_end
                );
                let support = exact_alignment_pairs(&bounded, &query_window, &target);
                let expected = (0..400u64)
                    .map(|offset| {
                        let query_position = 300 + offset - envelope.query_start;
                        let target_position = match strand {
                            Strand::Forward => 500 + offset,
                            Strand::Reverse => 4_301 - offset,
                        };
                        (query_position, target_position)
                    })
                    .collect::<std::collections::BTreeSet<_>>();
                assert!(expected.is_subset(&support));

                let full = AlignmentWorkspace::default()
                    .align_oriented(
                        &query,
                        &target,
                        0,
                        strand,
                        AlignmentConfig {
                            diagonal_offset: 200,
                            ..config.alignment
                        },
                    )
                    .unwrap();
                let (second_start, second_end) = match strand {
                    Strand::Forward => {
                        (second_target_start as u64, second_target_start as u64 + 400)
                    }
                    Strand::Reverse => (
                        4_802 - second_target_start as u64 - 400,
                        4_802 - second_target_start as u64,
                    ),
                };
                assert!(
                    full.target_interval.start <= second_start
                        && full.target_interval.end >= second_end
                );
            }

            let key = envelope_key(Strand::Forward);
            let grouped = form_regions(
                BTreeMap::from([(
                    key,
                    RegionHits::Many(vec![
                        SeedHit {
                            query: 315,
                            target: 515,
                            diagonal: 200,
                        },
                        SeedHit {
                            query: 335,
                            target: 535,
                            diagonal: 200,
                        },
                        SeedHit {
                            query: 1_515,
                            target: 1_715 + diagonal_distance,
                            diagonal: i128::from(200 + diagonal_distance),
                        },
                        SeedHit {
                            query: 1_535,
                            target: 1_735 + diagonal_distance,
                            diagonal: i128::from(200 + diagonal_distance),
                        },
                    ]),
                )]),
                64,
            );
            assert_eq!(grouped.len(), usize::from(diagonal_distance > 64) + 1);
        }
    }

    #[test]
    fn circular_seed_crossing_origin_matches_full_window_oracle() {
        let mut state = 0x9e3779b97f4a7c15u64;
        let query = (0..800)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize]
            })
            .collect::<Vec<_>>();
        let mut target = b"N".repeat(1_200);
        target[400..600].copy_from_slice(&query[600..800]);
        target[600..800].copy_from_slice(&query[..200]);

        let region = envelope_region(795, 595, 1);
        let config = TraceConfig::default();
        let envelope = fragment_envelope(
            &region,
            envelope_key(Strand::Forward),
            query.len() as u64,
            target.len() as u64,
            config,
        )
        .unwrap();
        assert_eq!(envelope.query_span, query.len() as u64);
        assert_eq!((envelope.target_start, envelope.target_end), (0, 1_200));
        let query_window =
            linearize_query(&query, envelope.query_start, envelope.query_span, true).unwrap();
        let query_relative = (795 + query.len() as u64 - envelope.query_start) % query.len() as u64;
        let diagonal_offset = i64::try_from(i128::from(595) - i128::from(query_relative)).unwrap();
        let mut bounded_config = config.alignment;
        bounded_config.diagonal_offset = diagonal_offset;
        let bounded = AlignmentWorkspace::default()
            .align_oriented(&query_window, &target, 0, Strand::Forward, bounded_config)
            .unwrap();
        let oracle = AlignmentWorkspace::default()
            .align_oriented(
                &query_window,
                &target,
                0,
                Strand::Forward,
                AlignmentConfig {
                    band_width: 2_048,
                    max_cells: 2_000_000,
                    ..config.alignment
                },
            )
            .unwrap();
        assert_eq!(bounded, oracle);
        assert_eq!(bounded.matches, 400);
        assert_eq!(bounded.query_interval.len(), 400);
        assert_eq!(bounded.target_interval, Interval::new(400, 800).unwrap());
    }

    type RegionRow = (RegionKey, u64, u64, u64, u64, i128, i128, u32);

    fn reference_form_regions(
        hits_by_contig: BTreeMap<RegionKey, Vec<SeedHit>>,
        max_diagonal_drift: u64,
    ) -> Vec<(RegionKey, RegionAccumulator)> {
        let mut output = Vec::new();
        for (key, mut hits) in hits_by_contig {
            hits.sort_unstable_by_key(|hit| (hit.query, hit.target));
            let mut regions = Vec::<RegionAccumulator>::new();
            for hit in hits {
                if let Some(region) = regions
                    .iter_mut()
                    .rev()
                    .find(|region| region.accepts(hit, max_diagonal_drift))
                {
                    region.add(hit);
                } else {
                    regions.push(RegionAccumulator::new(hit));
                }
            }
            output.extend(regions.into_iter().map(|region| (key, region)));
        }
        output
    }

    fn synthetic_region_events(
        groups: u32,
        hits_per_group: u32,
        extra_groups: u32,
    ) -> Vec<(RegionKey, SeedHit)> {
        assert!(groups.is_power_of_two() && extra_groups <= groups);
        let mut events = Vec::with_capacity((groups * hits_per_group + extra_groups) as usize);
        for hit in 0..hits_per_group + u32::from(extra_groups != 0) {
            for ordinal in 0..groups {
                let id = ordinal.wrapping_mul(2_654_435_761) & (groups - 1);
                if hit >= hits_per_group + u32::from(id < extra_groups) {
                    continue;
                }
                let query = u64::from(hit) * 21;
                let diagonal = i128::from(id % 97);
                events.push((
                    RegionKey {
                        metagenome_id: id / 4_096,
                        contig_id: id,
                        strand: if id.is_multiple_of(2) {
                            Strand::Forward
                        } else {
                            Strand::Reverse
                        },
                        k: if id.is_multiple_of(3) { 15 } else { 21 },
                    },
                    SeedHit {
                        query,
                        target: query + u64::try_from(diagonal).unwrap(),
                        diagonal,
                    },
                ));
            }
        }
        events
    }

    fn region_rows(regions: Vec<(RegionKey, RegionAccumulator)>) -> Vec<RegionRow> {
        regions
            .into_iter()
            .map(|(key, region)| {
                (
                    key,
                    region.query_start,
                    region.query_end,
                    region.target_start,
                    region.target_end,
                    region.diagonal_min,
                    region.diagonal_max,
                    region.hits,
                )
            })
            .collect()
    }

    fn timed_region_pipeline(
        events: &[(RegionKey, SeedHit)],
        candidate: bool,
    ) -> (u128, Vec<RegionRow>) {
        let started = std::time::Instant::now();
        let regions = if candidate {
            let mut groups = BTreeMap::<RegionKey, RegionHits>::new();
            for &(key, hit) in events {
                groups.entry(key).or_default().push(hit);
            }
            form_regions(groups, 64)
        } else {
            let mut groups = BTreeMap::<RegionKey, Vec<SeedHit>>::new();
            for &(key, hit) in events {
                groups.entry(key).or_default().push(hit);
            }
            reference_form_regions(groups, 64)
        };
        std::hint::black_box(&regions);
        let elapsed = started.elapsed().as_nanos();
        (elapsed, region_rows(regions))
    }

    fn region_checksum(rows: &[RegionRow]) -> u64 {
        rows.iter().fold(0xcbf29ce484222325u64, |mut state, row| {
            let strand = match row.0.strand {
                Strand::Forward => 0,
                Strand::Reverse => 1,
            };
            for value in [
                u64::from(row.0.metagenome_id),
                u64::from(row.0.contig_id),
                strand,
                u64::from(row.0.k),
                row.1,
                row.2,
                row.3,
                row.4,
                row.5 as u64,
                (row.5 >> 64) as u64,
                row.6 as u64,
                (row.6 >> 64) as u64,
                u64::from(row.7),
            ] {
                state = state.wrapping_mul(0x100000001b3) ^ value;
            }
            state
        })
    }

    #[test]
    #[ignore = "actual-source region storage diagnostic"]
    fn region_hits_inline_actual_map_diagnostic() {
        let cases = [
            ("all_singleton", synthetic_region_events(65_536, 1, 0)),
            (
                "mostly_singleton",
                synthetic_region_events(65_536, 1, 3_558),
            ),
            ("multihit", synthetic_region_events(8_192, 8, 0)),
        ];
        let candidate_first = std::env::var_os("JAM_REGION_HITS_CANDIDATE_FIRST").is_some();
        for (name, events) in cases {
            let (first_ns, first) = timed_region_pipeline(&events, candidate_first);
            let (second_ns, second) = timed_region_pipeline(&events, !candidate_first);
            let (reference_ns, reference, candidate_ns, candidate) = if candidate_first {
                (second_ns, second, first_ns, first)
            } else {
                (first_ns, first, second_ns, second)
            };
            assert_eq!(candidate, reference);
            println!(
                "region_hit_inline_oracle_v1 case={name} events={} rows={} reference_ns={reference_ns} candidate_ns={candidate_ns} checksum={}",
                events.len(),
                reference.len(),
                region_checksum(&reference)
            );
        }
    }

    #[test]
    fn groups_hits_across_diagonal_boundaries_and_small_indels() {
        let key = RegionKey {
            metagenome_id: 0,
            contig_id: 0,
            strand: Strand::Forward,
            k: 21,
        };
        let grouped = form_regions(
            BTreeMap::from([(
                key,
                RegionHits::Many(vec![
                    SeedHit {
                        query: 0,
                        target: 63,
                        diagonal: 63,
                    },
                    SeedHit {
                        query: 20,
                        target: 85,
                        diagonal: 65,
                    },
                ]),
            )]),
            4,
        );
        assert_eq!(grouped.len(), 1);
        assert_eq!(grouped[0].0, key);
        assert_eq!(grouped[0].1.hits, 2);
    }

    #[test]
    fn separates_hits_beyond_diagonal_drift() {
        let key = RegionKey {
            metagenome_id: 0,
            contig_id: 0,
            strand: Strand::Forward,
            k: 21,
        };
        let grouped = form_regions(
            BTreeMap::from([(
                key,
                RegionHits::Many(vec![
                    SeedHit {
                        query: 0,
                        target: 10,
                        diagonal: 10,
                    },
                    SeedHit {
                        query: 20,
                        target: 40,
                        diagonal: 20,
                    },
                ]),
            )]),
            4,
        );
        assert_eq!(grouped.len(), 2);
    }

    #[test]
    fn coalesces_only_overlapping_or_adjacent_target_spans() {
        let mut spans = vec![(100, 110), (15, 20), (0, 10), (8, 15), (200, 205)];
        coalesce_spans(&mut spans);
        assert_eq!(spans, vec![(0, 20), (100, 110), (200, 205)]);
    }

    #[test]
    fn extracts_every_tagged_k15_rescue_position() {
        let seeds = query_seeds(b"ACGTACGTACGTACGTACGTA", 21, true, false).unwrap();
        assert_eq!(
            seeds
                .iter()
                .filter(|seed| seed.packed_key & RESCUE_K15_TAG != 0)
                .map(|seed| seed.position)
                .collect::<Vec<_>>(),
            (0..7).collect::<Vec<_>>()
        );
        assert_eq!(
            seeds
                .iter()
                .filter(|seed| seed.packed_key & RESCUE_K15_TAG == 0)
                .count(),
            1
        );
        assert_eq!(minimum_region_hits(15, 2), 3);
        assert_eq!(minimum_region_hits(21, 2), 2);
    }

    fn window_dna(mut state: u64, length: usize) -> Vec<u8> {
        (0..length)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize]
            })
            .collect()
    }

    fn window_reverse_complement(sequence: &[u8]) -> Vec<u8> {
        sequence
            .iter()
            .rev()
            .map(|base| match base {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => b'N',
            })
            .collect()
    }

    fn complete_window_alignment(
        query: &[u8],
        target: &[u8],
        query_start: u64,
        query_span: u64,
        target_start: u64,
        target_end: u64,
        strand: Strand,
        diagonal_offset: i64,
        config: TraceConfig,
    ) -> Alignment {
        let query_window =
            linearize_query(query, query_start, query_span, config.circular).unwrap();
        let target_window = &target[target_start as usize..target_end as usize];
        let mut alignment_config = config.alignment;
        alignment_config.diagonal_offset = diagonal_offset;
        let mut workspace = AlignmentWorkspace::default();
        let core = workspace
            .align_oriented(
                &query_window,
                target_window,
                target_start,
                strand,
                alignment_config,
            )
            .unwrap();
        let completed = workspace
            .complete_endpoints(
                core.clone(),
                &query_window,
                target_window,
                target_start,
                config.endpoint_bases,
                alignment_config,
            )
            .unwrap()
            .alignment;
        if completed.identity() >= config.min_identity {
            completed
        } else {
            core
        }
    }

    fn compare_envelope_to_full_window(
        query: &[u8],
        target: &[u8],
        query_anchor: u64,
        oriented_target_anchor: u64,
        strand: Strand,
        config: TraceConfig,
    ) {
        let region = envelope_region(query_anchor, oriented_target_anchor, 1);
        let envelope = fragment_envelope(
            &region,
            envelope_key(strand),
            query.len() as u64,
            target.len() as u64,
            config,
        )
        .unwrap();
        let query_relative = if query_anchor >= envelope.query_start {
            query_anchor - envelope.query_start
        } else {
            query_anchor + query.len() as u64 - envelope.query_start
        };
        let oriented_start = match strand {
            Strand::Forward => envelope.target_start,
            Strand::Reverse => target.len() as u64 - envelope.target_end,
        };
        let bounded = complete_window_alignment(
            query,
            target,
            envelope.query_start,
            envelope.query_span,
            envelope.target_start,
            envelope.target_end,
            strand,
            i64::try_from(
                i128::from(oriented_target_anchor - oriented_start) - i128::from(query_relative),
            )
            .unwrap(),
            config,
        );
        let oracle = complete_window_alignment(
            query,
            target,
            0,
            query.len() as u64,
            0,
            target.len() as u64,
            strand,
            i64::try_from(i128::from(oriented_target_anchor) - i128::from(query_anchor)).unwrap(),
            TraceConfig {
                endpoint_bases: config.endpoint_bases,
                alignment: AlignmentConfig {
                    band_width: 256,
                    max_cells: 2_000_000,
                    ..config.alignment
                },
                ..config
            },
        );
        assert_eq!(bounded.score, oracle.score);
        assert_eq!(bounded.strand, oracle.strand);
        assert_eq!(bounded.target_interval, oracle.target_interval);
        assert_eq!(bounded.matches, oracle.matches);
        assert_eq!(bounded.substitutions, oracle.substitutions);
        assert_eq!(bounded.insertions, oracle.insertions);
        assert_eq!(bounded.deletions, oracle.deletions);
        assert_eq!(bounded.cigar, oracle.cigar);
        assert_eq!(bounded.edit_script, oracle.edit_script);
        assert_eq!(
            query_segments(
                envelope.query_start,
                bounded.query_interval,
                query.len() as u64,
                config.circular,
            )
            .unwrap(),
            query_segments(
                0,
                oracle.query_interval,
                query.len() as u64,
                config.circular
            )
            .unwrap(),
        );
    }

    #[test]
    fn fragment_envelopes_match_full_windows_across_bounded_error_cases() {
        let component = window_dna(0x123456789abcdef0, 800);
        let config = TraceConfig {
            circular: false,
            ..TraceConfig::default()
        };
        for anchor in [384, 700, 780] {
            let mut query = window_dna(0x5555555555555555, 1_400);
            query[300..1_100].copy_from_slice(&component);
            for strand in [Strand::Forward, Strand::Reverse] {
                let mut oriented_target = window_dna(0xaaaaaaaaaaaaaaaa, 1_800);
                oriented_target[500..1_300].copy_from_slice(&component);
                let target = if strand == Strand::Forward {
                    oriented_target
                } else {
                    window_reverse_complement(&oriented_target)
                };
                compare_envelope_to_full_window(
                    &query,
                    &target,
                    300 + anchor,
                    500 + anchor,
                    strand,
                    config,
                );
            }
        }

        for (name, query_component, target_component, anchor) in [
            (
                "insertion",
                component.clone(),
                {
                    let mut value = component.clone();
                    value.splice(400..400, window_dna(31, 30));
                    value
                },
                700,
            ),
            (
                "deletion",
                component.clone(),
                {
                    let mut value = component.clone();
                    value.drain(385..415);
                    value
                },
                700,
            ),
            (
                "ambiguity",
                {
                    let mut value = component.clone();
                    value[360..392].fill(b'N');
                    value
                },
                {
                    let mut value = component.clone();
                    value[360..392].fill(b'N');
                    value
                },
                700,
            ),
            (
                "absent_left",
                component.clone(),
                component[200..].to_vec(),
                700,
            ),
            (
                "absent_right",
                component.clone(),
                component[..600].to_vec(),
                384,
            ),
        ] {
            let mut query = window_dna(0x5555555555555555, 1_400);
            query[300..1_100].copy_from_slice(&query_component);
            let mut target = window_dna(0xaaaaaaaaaaaaaaaa, 1_800);
            target[500..500 + target_component.len()].copy_from_slice(&target_component);
            let target_anchor = match name {
                "insertion" if anchor >= 400 => 500 + anchor + 30,
                "deletion" if anchor >= 415 => 500 + anchor - 30,
                "absent_left" => 500 + anchor - 200,
                _ => 500 + anchor,
            };
            compare_envelope_to_full_window(
                &query,
                &target,
                300 + anchor,
                target_anchor,
                Strand::Forward,
                config,
            );
        }

        for target_start in [0, 1_000] {
            let mut query = window_dna(0x5555555555555555, 1_400);
            query[300..1_100].copy_from_slice(&component);
            let mut target = window_dna(0xaaaaaaaaaaaaaaaa, 1_800);
            target[target_start..target_start + 800].copy_from_slice(&component);
            compare_envelope_to_full_window(
                &query,
                &target,
                1_000,
                (target_start + 700) as u64,
                Strand::Forward,
                config,
            );
        }
    }

    #[test]
    fn long_contig_displaced_anchor_matches_larger_bounded_reference() {
        let component = window_dna(0x123456789abcdef0, 800);
        let mut query = window_dna(0x5555555555555555, 1_400);
        query[300..1_100].copy_from_slice(&component);
        let mut target = window_dna(0xaaaaaaaaaaaaaaaa, 66_000);
        target[50_000..50_800].copy_from_slice(&component);
        let query_anchor = 1_000;
        let target_anchor = 50_700;
        let config = TraceConfig {
            circular: false,
            ..TraceConfig::default()
        };
        let region = envelope_region(query_anchor, target_anchor, 1);
        let envelope = fragment_envelope(
            &region,
            envelope_key(Strand::Forward),
            query.len() as u64,
            target.len() as u64,
            config,
        )
        .unwrap();
        assert!(envelope.target_start > 0 && envelope.target_end < target.len() as u64);
        let bounded = complete_window_alignment(
            &query,
            &target,
            envelope.query_start,
            envelope.query_span,
            envelope.target_start,
            envelope.target_end,
            Strand::Forward,
            i64::try_from(
                i128::from(target_anchor - envelope.target_start)
                    - i128::from(query_anchor - envelope.query_start),
            )
            .unwrap(),
            config,
        );
        let reference_start = 49_000;
        let reference_end = 52_000;
        let reference = complete_window_alignment(
            &query,
            &target,
            0,
            query.len() as u64,
            reference_start,
            reference_end,
            Strand::Forward,
            i64::try_from(i128::from(target_anchor - reference_start) - i128::from(query_anchor))
                .unwrap(),
            TraceConfig {
                alignment: AlignmentConfig {
                    band_width: 256,
                    max_cells: 2_000_000,
                    ..config.alignment
                },
                ..config
            },
        );
        assert_eq!(bounded.score, reference.score);
        assert_eq!(bounded.target_interval, reference.target_interval);
        assert_eq!(bounded.cigar, reference.cigar);
        assert_eq!(
            query_segments(
                envelope.query_start,
                bounded.query_interval,
                query.len() as u64,
                false,
            )
            .unwrap(),
            query_segments(0, reference.query_interval, query.len() as u64, false).unwrap(),
        );
    }

    #[test]
    fn circular_retry_can_join_two_subthreshold_fragment_halves() {
        let target = window_dna(71, 40);
        let query = [
            target[20..].to_vec(),
            b"N".repeat(40),
            target[..20].to_vec(),
        ]
        .concat();
        let config = TraceConfig {
            circular: true,
            min_aligned_bases: 40,
            ..TraceConfig::default()
        };
        let task = AlignmentTask {
            metagenome_id: 0,
            contig_id: 0,
            strand: Strand::Forward,
            query_start: 0,
            query_span: query.len() as u64,
            target_start: 0,
            target_end: target.len() as u64,
            diagonal_offset: 20,
            parent: 0,
            island: None,
        };
        let mut initial_config = config.alignment;
        initial_config.diagonal_offset = task.diagonal_offset;
        let mut workspace = AlignmentWorkspace::default();
        let initial = align_task_window(
            &mut workspace,
            &query,
            &target,
            0,
            Strand::Forward,
            initial_config,
            config,
        )
        .unwrap()
        .unwrap();
        assert_eq!(initial.core.query_interval.len(), 20);
        assert_eq!(initial.selected.query_interval.len(), 20);
        assert!(!alignment_accepted(&initial.selected, config));

        let (retry_start, retry_diagonal) =
            circular_retry(&initial.core, &task, query.len() as u64)
                .unwrap()
                .unwrap();
        assert_eq!(retry_start, 20);
        let retry_query = linearize_query(&query, retry_start, query.len() as u64, true).unwrap();
        initial_config.diagonal_offset = retry_diagonal;
        let retry = align_task_window(
            &mut workspace,
            &retry_query,
            &target,
            0,
            Strand::Forward,
            initial_config,
            config,
        )
        .unwrap()
        .unwrap();
        assert_eq!(retry.selected.query_interval.len(), 40);
        assert_eq!(
            retry.selected.target_interval,
            Interval::new(0, 40).unwrap()
        );
        assert!(alignment_accepted(&retry.selected, config));
        assert!(retry_improves(&initial.selected, &retry.selected, false));
    }

    #[test]
    fn distant_collinear_anchor_task_runs_in_resident_chunks_across_workers() {
        use crate::jidx_writer::{ContigInput, JidxInput, JidxWriter, MetagenomeInput};
        // Shape of the rejected BCF tasks: two anchored segments joined by a long gap on a
        // contig longer than the short-contig window, so every run uses the same bounded window.
        let directory = tempfile::tempdir().unwrap();
        let target = window_dna(0x0f1e_2d3c_4b5a_6978, 70_000);
        let mut query = window_dna(0x1357_9bdf_2468_ace0, 12_000);
        query[1_000..1_150].copy_from_slice(&target[20_000..20_150]);
        query[11_000..11_150].copy_from_slice(&target[30_000..30_150]);
        let bgzf_path = directory.path().join("target.bgz");
        let mut raw = b">contig\n".to_vec();
        for line in target.chunks(80) {
            raw.extend_from_slice(line);
            raw.push(b'\n');
        }
        let mut writer = bgzf::io::Writer::new(File::create(&bgzf_path).unwrap());
        let mut blocks = Vec::new();
        for (ordinal, chunk) in raw.chunks(32_000).enumerate() {
            if ordinal > 0 {
                blocks.push((writer.position(), (ordinal * 32_000) as u64));
            }
            writer.write_all(chunk).unwrap();
            writer.flush().unwrap();
        }
        writer.finish().unwrap();
        let gzi_path = directory.path().join("target.gzi");
        gzi::fs::write(&gzi_path, &gzi::Index::from(blocks)).unwrap();
        let bytes = std::fs::read(&bgzf_path).unwrap();
        let reference = directory.path().join("reference.jidx");
        let mut jidx = JidxWriter::new(
            &reference,
            &JidxInput {
                k: 15,
                rescue_k15: false,
                minimizer_window: 16,
                jam_sha256: [1; 32],
                manifest_sha256: [2; 32],
            },
        )
        .unwrap();
        jidx.begin_metagenome(MetagenomeInput {
            name: "target".into(),
            bgzf_uri: bgzf_path.to_str().unwrap().to_owned(),
            bgzf_bytes: bytes.len() as u64,
            bgzf_sha256: crate::jidx::sha256(&bytes),
            gzi: std::fs::read(gzi_path).unwrap(),
        })
        .unwrap();
        jidx.begin_contig(ContigInput {
            name: "contig".into(),
            length: target.len() as u64,
            fasta_offset: 8,
            line_bases: 80,
            line_width: 81,
        })
        .unwrap();
        jidx.finish().unwrap();
        let shared = directory.path().join("target.shared");
        crate::shared_writer::build_shared_index(&reference, &shared, 16).unwrap();

        let unchunked = TraceConfig {
            use_sketch: false,
            circular: false,
            ..TraceConfig::default()
        };
        let chunked = TraceConfig {
            alignment: AlignmentConfig {
                max_cells: 500_000,
                ..unchunked.alignment
            },
            ..unchunked
        };
        let without_reads = |mut result: TraceResult| {
            for metagenome in &mut result.metagenomes {
                metagenome.compressed_bytes_read = 0;
                metagenome.range_requests = 0;
                metagenome.bgzf_blocks_decoded = 0;
            }
            result
        };
        let expected = without_reads(
            TraceEngine::open_shared(&shared, None)
                .unwrap()
                .search("gap", &query, unchunked)
                .unwrap(),
        );
        let fragments = expected.metagenomes[0]
            .mosaic
            .primary
            .iter()
            .map(|selected| selected.fragment.alignment.query_interval)
            .collect::<Vec<_>>();
        // One task yields one local alignment: the better of the two anchored segments.
        assert!(fragments.iter().any(|interval| {
            [1_000, 11_000]
                .iter()
                .any(|&start| interval.start <= start && interval.end >= start + 150)
        }));
        for workers in [1, 4] {
            let mut engine = TraceEngine::open_shared(&shared, None).unwrap();
            engine.observed = true;
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap();
            let result = pool.install(|| engine.search("gap", &query, chunked).unwrap());
            assert_eq!(without_reads(result), expected);
            let work = engine.batch_stats().alignment_work;
            assert!(work.local_chunked_passes > 0 && work.local_chunks > work.local_chunked_passes);
            assert!(work.local_recomputed_cells > 0);
        }
    }

    /// Writes one BGZF source per metagenome, each with two contigs in several blocks, and a
    /// shared index over them. Returns the index, the source paths and the contigs in ID order.
    fn write_range_sources(
        directory: &Path,
        metagenomes: usize,
    ) -> (PathBuf, Vec<PathBuf>, Vec<Vec<u8>>) {
        use crate::jidx_writer::{ContigInput, JidxInput, JidxWriter, MetagenomeInput};
        let reference = directory.join("sources.jidx");
        let mut jidx = JidxWriter::new(
            &reference,
            &JidxInput {
                k: 15,
                rescue_k15: false,
                minimizer_window: 64,
                jam_sha256: [1; 32],
                manifest_sha256: [2; 32],
            },
        )
        .unwrap();
        let (mut paths, mut contigs) = (Vec::new(), Vec::new());
        for metagenome in 0..metagenomes {
            let sequences = [0, 1].map(|contig| {
                window_dna((3 + 2 * metagenome + contig) as u64, 6_000 + 700 * contig)
            });
            let mut raw = Vec::new();
            let mut offsets = Vec::new();
            for (name, sequence) in ["a", "b"].iter().zip(&sequences) {
                raw.extend_from_slice(format!(">{name}\n").as_bytes());
                offsets.push(raw.len() as u64);
                for line in sequence.chunks(80) {
                    raw.extend_from_slice(line);
                    raw.push(b'\n');
                }
            }
            let path = directory.join(format!("source-{metagenome}.bgz"));
            let mut writer = bgzf::io::Writer::new(File::create(&path).unwrap());
            let mut blocks = Vec::new();
            for (ordinal, chunk) in raw.chunks(3_000).enumerate() {
                if ordinal > 0 {
                    blocks.push((writer.position(), (ordinal * 3_000) as u64));
                }
                writer.write_all(chunk).unwrap();
                writer.flush().unwrap();
            }
            writer.finish().unwrap();
            let mut gzi = gzi::io::Writer::new(Vec::new());
            gzi.write_index(&gzi::Index::from(blocks)).unwrap();
            let bytes = std::fs::read(&path).unwrap();
            jidx.begin_metagenome(MetagenomeInput {
                name: format!("source-{metagenome}"),
                bgzf_uri: path.to_str().unwrap().to_owned(),
                bgzf_bytes: bytes.len() as u64,
                bgzf_sha256: sha256(&bytes),
                gzi: gzi.into_inner(),
            })
            .unwrap();
            for ((name, sequence), offset) in ["a", "b"].iter().zip(&sequences).zip(offsets) {
                jidx.begin_contig(ContigInput {
                    name: (*name).into(),
                    length: sequence.len() as u64,
                    fasta_offset: offset,
                    line_bases: 80,
                    line_width: 81,
                })
                .unwrap();
            }
            paths.push(path);
            contigs.extend(sequences);
        }
        jidx.finish().unwrap();
        let shared = directory.join("sources.shared");
        crate::shared_writer::build_shared_index(&reference, &shared, 64).unwrap();
        (shared, paths, contigs)
    }

    /// Overlapping, adjacent and repeated spans on every contig, interleaved across metagenomes.
    fn range_tasks(engine: &TraceEngine, contigs: &[Vec<u8>]) -> Vec<AlignmentTask> {
        let mut tasks = Vec::new();
        for (start, end) in [
            (10, 2_500),
            (5_000, 5_990),
            (2_000, 3_100),
            (3_100, 3_300),
            (10, 2_500),
        ] {
            for contig_id in (0..contigs.len() as u32).rev() {
                let contig = engine.index.contig(contig_id).unwrap().unwrap();
                tasks.push(AlignmentTask {
                    metagenome_id: contig.metagenome_id,
                    contig_id,
                    strand: Strand::Forward,
                    query_start: 0,
                    query_span: 1,
                    target_start: start,
                    target_end: end,
                    diagonal_offset: 0,
                    parent: 0,
                    island: None,
                });
            }
        }
        tasks
    }

    type ComparableRanges = BTreeMap<(MetagenomeId, ContigId), Vec<(u64, u64, Vec<u8>)>>;

    fn comparable(
        loaded: BTreeMap<(MetagenomeId, ContigId), Vec<LoadedRange>>,
    ) -> ComparableRanges {
        loaded
            .into_iter()
            .map(|(key, ranges)| {
                let ranges = ranges
                    .into_iter()
                    .map(|range| (range.offset, range.end, range.sequence));
                (key, ranges.collect())
            })
            .collect()
    }

    /// Serial reference: one reader per metagenome in ID order, each coalesced span read once.
    fn serial_ranges(
        engine: &TraceEngine,
        tasks: &[AlignmentTask],
        cache: Option<&Arc<BgzfBlockCache>>,
    ) -> (
        ComparableRanges,
        HashMap<MetagenomeId, (crate::range_source::RangeStats, u64)>,
    ) {
        let mut spans = BTreeMap::<(MetagenomeId, ContigId), Vec<(u64, u64)>>::new();
        for task in tasks {
            spans
                .entry((task.metagenome_id, task.contig_id))
                .or_default()
                .push((task.target_start, task.target_end));
        }
        spans.values_mut().for_each(coalesce_spans);
        let (mut loaded, mut reads) = (BTreeMap::new(), HashMap::new());
        for metagenome_id in spans.keys().map(|&(id, _)| id).collect::<BTreeSet<_>>() {
            let source = engine.index.metagenome(metagenome_id).unwrap().unwrap();
            let mut reader = match cache {
                Some(cache) => BgzfReader::open_with_cache(source, None, true, Arc::clone(cache)),
                None => BgzfReader::open(source, None, true),
            }
            .unwrap();
            for (&key, contig_spans) in spans.range((metagenome_id, 0)..=(metagenome_id, u32::MAX))
            {
                let contig = engine.index.contig(key.1).unwrap().unwrap();
                let ranges = contig_spans.iter().map(|&(start, end)| {
                    (
                        start,
                        end,
                        reader.read_contig_range(contig, start, end).unwrap(),
                    )
                });
                loaded.insert(key, ranges.collect());
            }
            reads.insert(
                metagenome_id,
                (reader.range_stats(), reader.blocks_decoded()),
            );
        }
        (loaded, reads)
    }

    #[test]
    fn parallel_range_loads_match_serial_reads_and_accounting() {
        let directory = tempfile::tempdir().unwrap();
        let (shared, _, contigs) = write_range_sources(directory.path(), 4);
        let engine = TraceEngine::open_shared(&shared, None).unwrap();
        let tasks = range_tasks(&engine, &contigs);
        for cached in [false, true] {
            // A fresh default batch cache per load keeps every block resident, as in one query.
            let cache = || {
                cached
                    .then(|| Arc::new(BgzfBlockCache::new(DEFAULT_BATCH_BGZF_CACHE_BYTES).unwrap()))
            };
            let (expected, expected_reads) = serial_ranges(&engine, &tasks, cache().as_ref());
            assert_eq!(expected.len(), contigs.len());
            for (&(_, contig_id), ranges) in &expected {
                assert_eq!(ranges.len(), 2);
                for (start, end, sequence) in ranges {
                    assert_eq!(
                        sequence,
                        &contigs[contig_id as usize][*start as usize..*end as usize]
                    );
                }
            }
            assert!(
                expected_reads
                    .values()
                    .all(|(stats, blocks)| stats.read_requests > 0 && *blocks >= 3)
            );
            for threads in [1, 3, 8] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                let (loaded, reads) = pool
                    .install(|| engine.load_ranges(&tasks, true, cache().as_ref()))
                    .unwrap();
                assert_eq!(
                    comparable(loaded),
                    expected,
                    "{threads} threads, cache {cached}"
                );
                assert_eq!(reads, expected_reads, "{threads} threads, cache {cached}");
            }
        }
    }

    #[test]
    fn parallel_range_loads_respect_shared_block_admission() {
        use crate::bgzf_cache::{MAX_BGZF_BLOCK_BYTES, MAX_CONCURRENT_BGZF_DECODES};
        let directory = tempfile::tempdir().unwrap();
        let (shared, _, contigs) = write_range_sources(directory.path(), 4);
        let engine = TraceEngine::open_shared(&shared, None).unwrap();
        let tasks = range_tasks(&engine, &contigs);
        let (expected, _) = serial_ranges(&engine, &tasks, None);
        // The cache holds fewer blocks than the loads read, so parallel readers share a small
        // byte budget and evict blocks.
        let cache = Arc::new(BgzfBlockCache::new(4 * MAX_BGZF_BLOCK_BYTES).unwrap());
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(8)
            .build()
            .unwrap();
        let (loaded, _) = pool
            .install(|| engine.load_ranges(&tasks, true, Some(&cache)))
            .unwrap();
        assert_eq!(comparable(loaded), expected);
        let stats = cache.stats();
        assert!(stats.peak_loading_blocks <= MAX_CONCURRENT_BGZF_DECODES);
        assert!(stats.peak_accounted_bytes <= stats.capacity_bytes);
        assert_eq!((stats.loading_blocks, stats.reserved_bytes), (0, 0));
        assert!(stats.evictions > 0);
    }

    fn wait_until(mut ready: impl FnMut() -> bool) {
        let started = Instant::now();
        while !ready() {
            assert!(started.elapsed() < std::time::Duration::from_secs(5));
            std::thread::sleep(std::time::Duration::from_millis(1));
        }
    }

    #[cfg(unix)]
    #[test]
    fn identity_change_inside_hit_loop_fails_before_output_publication() {
        use crate::shared_format::SharedError;
        let directory = tempfile::tempdir().unwrap();
        let (shared, _, contigs) = write_range_sources(directory.path(), 1);
        let (id, sequence) = ("identity-phase-boundary", &contigs[0][..900]);
        let query = directory.path().join("identity.fa");
        std::fs::write(
            &query,
            format!(">{id}\n{}\n", String::from_utf8_lossy(sequence)),
        )
        .unwrap();
        let config = TraceConfig {
            use_sketch: false,
            circular: false,
            ..TraceConfig::default()
        };
        let trace = |output: &Path| {
            handle_trace_command(TraceArgs {
                query: query.clone(),
                input: TraceInput::Shared {
                    path: shared.clone(),
                    read_stats: None,
                    query_topology_header: false,
                },
                audit_index: false,
                output: output.to_owned(),
                query_id: None,
                config,
                s3: None,
                force: false,
            })
        };
        // Contig metadata inside the hit loop is read without identity checks; this changes the
        // file identity after those reads and before the check that ends the phase.
        let change_identity_after_hit_loop = |seconds: u64| {
            let index = shared.clone();
            *HIT_LOOP_END_HOOK.lock().unwrap() = Some((
                id.to_owned(),
                Box::new(move || {
                    let times = std::fs::FileTimes::new().set_modified(
                        std::time::SystemTime::UNIX_EPOCH + std::time::Duration::from_secs(seconds),
                    );
                    let file = File::options().write(true).open(index).unwrap();
                    file.set_times(times).unwrap();
                }),
            ));
        };

        let unchanged = directory.path().join("unchanged.jsonl");
        trace(&unchanged).unwrap();
        let result = serde_json::from_str::<serde_json::Value>(
            &std::fs::read_to_string(&unchanged).unwrap(),
        )
        .unwrap();
        assert!(
            result["metagenomes"][0]["mosaic"]["covered_bases"]
                .as_u64()
                .unwrap()
                > 0
        );
        let engine = TraceEngine::open_shared_observed(&shared, None, true).unwrap();
        engine.search(id, sequence, config).unwrap();
        assert!(engine.batch_stats().geometric_hits > 0);

        change_identity_after_hit_loop(42);
        let changed = directory.path().join("changed.jsonl");
        let error = trace(&changed).unwrap_err();
        assert!(HIT_LOOP_END_HOOK.lock().unwrap().is_none());
        assert!(
            format!("{error:#}").contains(&SharedError::SourceChanged.to_string()),
            "{error:#}"
        );
        assert!(!changed.exists());
        assert!(std::fs::read_dir(directory.path()).unwrap().all(|entry| {
            !entry
                .unwrap()
                .file_name()
                .to_string_lossy()
                .starts_with(".jam-trace-")
        }));

        // The phase check fails before regions are formed from the hits.
        let engine = TraceEngine::open_shared_observed(&shared, None, true).unwrap();
        change_identity_after_hit_loop(43);
        let error = engine.search(id, sequence, config).err().unwrap();
        assert!(HIT_LOOP_END_HOOK.lock().unwrap().is_none());
        assert_eq!(
            error.to_string(),
            TraceError::Io(io::Error::other(SharedError::SourceChanged)).to_string()
        );
        assert_eq!(engine.batch_stats().geometric_hits, 0);
    }
}
