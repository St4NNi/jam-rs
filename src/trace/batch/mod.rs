//! Reusable grouped execution for multi-query trace workers.
//!
//! Batch mode deliberately stops at a small adapter boundary.  The batch
//! planner owns query validation, witness preparation, candidate routing, and
//! request unioning.  A caller supplies the archive/alignment adapter that
//! knows how to decode one JMA page and how to turn the decoded data into its
//! normal per-query result.  This keeps batch scheduling independent from the
//! current JMA reader and from the alignment implementation.

use crate::router::postings::PostingSource;
use crate::router::search::{CandidateRouter, TieredQueryWitness};
use crate::router::witness::extract_nested_witnesses;
use crate::router::{RoutedCandidate, WitnessKey, WitnessScheme};
use crate::trace::contracts::{
    BatchCandidateWork, BatchResultKey, PreparedBatchQuery, WitnessPlanningRequest,
};
use crate::trace::model::QueryDescriptor;
use crate::trace::query_plan::{QueryPlanError, plan_witness_tier};
use std::collections::{BTreeMap, BTreeSet};
use std::error::Error;
use std::fmt::{self, Display, Formatter};

pub mod fixture;
pub mod runtime;

/// One query record supplied to [`prepare_batch_queries`].
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BatchQueryRecord {
    pub descriptor: QueryDescriptor,
    pub sequence: Vec<u8>,
}

impl BatchQueryRecord {
    #[must_use]
    pub fn new(descriptor: QueryDescriptor, sequence: Vec<u8>) -> Self {
        Self {
            descriptor,
            sequence,
        }
    }
}

/// Settings shared by all records in one batch preparation pass.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BatchPrepareOptions {
    pub witness_planning: WitnessPlanningRequest,
    pub query_window_size: u32,
}

impl BatchPrepareOptions {
    #[must_use]
    pub const fn new(witness_planning: WitnessPlanningRequest, query_window_size: u32) -> Self {
        Self {
            witness_planning,
            query_window_size,
        }
    }
}

/// Prepared records and the one witness scheme shared by the batch.
#[derive(Clone, Debug, PartialEq)]
pub struct BatchPlan {
    pub scheme: WitnessScheme,
    pub queries: Vec<PreparedBatchQuery>,
}

impl BatchPlan {
    #[must_use]
    pub fn query(&self, query_index: u32) -> Option<&PreparedBatchQuery> {
        self.queries
            .get(usize::try_from(query_index).ok()?)
            .filter(|query| query.query_index == query_index)
    }

    #[must_use]
    pub fn query_count(&self) -> usize {
        self.queries.len()
    }
}

/// Planning and input-validation failures for a complete query batch.
#[derive(Debug, Eq, PartialEq)]
pub enum BatchPlanError {
    EmptyBatch,
    EmptyQueryId,
    DuplicateQueryId(String),
    EmptySequence(String),
    QueryLengthMismatch {
        query_id: String,
        descriptor_length: u64,
        sequence_length: u64,
    },
    ZeroWindowSize,
    WitnessPlanning(String),
    WitnessExtraction(String),
}

impl Display for BatchPlanError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyBatch => formatter.write_str("batch must contain at least one query"),
            Self::EmptyQueryId => formatter.write_str("query_id must not be empty"),
            Self::DuplicateQueryId(query_id) => {
                write!(formatter, "query_id is duplicated: {query_id}")
            }
            Self::EmptySequence(query_id) => {
                write!(formatter, "query sequence is empty: {query_id}")
            }
            Self::QueryLengthMismatch {
                query_id,
                descriptor_length,
                sequence_length,
            } => write!(
                formatter,
                "query length mismatch for {query_id}: descriptor={descriptor_length}, sequence={sequence_length}"
            ),
            Self::ZeroWindowSize => {
                formatter.write_str("query window size must be greater than zero")
            }
            Self::WitnessPlanning(message) => {
                write!(formatter, "witness planning failed: {message}")
            }
            Self::WitnessExtraction(message) => {
                write!(formatter, "witness extraction failed: {message}")
            }
        }
    }
}

impl Error for BatchPlanError {}

/// Validate all records before allocating witness state for any record.
pub fn validate_batch_queries(records: &[BatchQueryRecord]) -> Result<(), BatchPlanError> {
    if records.is_empty() {
        return Err(BatchPlanError::EmptyBatch);
    }

    let mut query_ids = BTreeSet::new();
    for record in records {
        let query_id = record.descriptor.query_id.trim();
        if query_id.is_empty() {
            return Err(BatchPlanError::EmptyQueryId);
        }
        if !query_ids.insert(record.descriptor.query_id.clone()) {
            return Err(BatchPlanError::DuplicateQueryId(
                record.descriptor.query_id.clone(),
            ));
        }
        if record.sequence.is_empty() {
            return Err(BatchPlanError::EmptySequence(
                record.descriptor.query_id.clone(),
            ));
        }
        let sequence_length = u64::try_from(record.sequence.len()).unwrap_or(u64::MAX);
        if record.descriptor.query_length != sequence_length {
            return Err(BatchPlanError::QueryLengthMismatch {
                query_id: record.descriptor.query_id.clone(),
                descriptor_length: record.descriptor.query_length,
                sequence_length,
            });
        }
    }
    Ok(())
}

/// Collect records from a manifest/FASTA adapter and validate the complete
/// batch before witness preparation begins.
pub fn load_batch_queries<I>(records: I) -> Result<Vec<BatchQueryRecord>, BatchPlanError>
where
    I: IntoIterator<Item = BatchQueryRecord>,
{
    let records = records.into_iter().collect::<Vec<_>>();
    validate_batch_queries(&records)?;
    Ok(records)
}

/// Prepare each query's dense witness set exactly once.
pub fn prepare_batch_queries(
    records: &[BatchQueryRecord],
    scheme: &WitnessScheme,
    options: BatchPrepareOptions,
) -> Result<BatchPlan, BatchPlanError> {
    validate_batch_queries(records)?;
    if options.query_window_size == 0 {
        return Err(BatchPlanError::ZeroWindowSize);
    }

    let mut queries = Vec::with_capacity(records.len());
    for (query_index, record) in records.iter().enumerate() {
        let witness_plan = plan_witness_tier(options.witness_planning, scheme)
            .map_err(|error: QueryPlanError| BatchPlanError::WitnessPlanning(error.to_string()))?;
        let nested = extract_nested_witnesses(&record.sequence, scheme, options.query_window_size)
            .map_err(|error| BatchPlanError::WitnessExtraction(error.to_string()))?;
        // The base set is materialized once.  Coarser tiers are represented by
        // filtering these already-created occurrences in `tiered_witnesses`.
        let witnesses = nested
            .at_scale(scheme.base_scale)
            .map_err(|error| BatchPlanError::WitnessExtraction(error.to_string()))?;
        queries.push(PreparedBatchQuery {
            query_index: u32::try_from(query_index).unwrap_or(u32::MAX),
            descriptor: record.descriptor.clone(),
            sequence: record.sequence.clone(),
            witnesses,
            witness_plan,
        });
    }

    Ok(BatchPlan {
        scheme: scheme.clone(),
        queries,
    })
}

/// Expand one prepared dense set into the sparse-to-dense tier sequence used
/// by the collection router.  No sequence scan or k-mer hashing occurs here.
#[must_use]
pub fn tiered_witnesses(
    query: &PreparedBatchQuery,
    scheme: &WitnessScheme,
) -> Vec<TieredQueryWitness> {
    let mut scales = match query.witness_plan.selected_witness_tier {
        Some(selected) => {
            let mut scales = vec![selected];
            scales.extend(query.witness_plan.available_denser_tiers.iter().copied());
            scales
        }
        None => scheme.available_scales.clone(),
    };
    scales.sort_unstable_by(|left, right| right.cmp(left));
    scales.dedup();

    let mut expanded = Vec::new();
    for witness in &query.witnesses {
        for &scale in &scales {
            if scheme
                .includes_hash(witness.key.jamhash, scale)
                .unwrap_or(false)
            {
                expanded.push(TieredQueryWitness::new(scale, witness.clone()));
            }
        }
    }
    expanded
}

/// Adapter boundary for collection routing.  It is intentionally independent
/// of the concrete `.jam` or JWR reader.
pub trait BatchRouter {
    type Error: Error + Send + Sync + 'static;

    fn route_query(
        &mut self,
        query: &PreparedBatchQuery,
        witnesses: &[TieredQueryWitness],
    ) -> Result<Vec<RoutedCandidate>, Self::Error>;
}

impl<'a, S: PostingSource> BatchRouter for CandidateRouter<'a, S> {
    type Error = crate::router::search::SearchError;

    fn route_query(
        &mut self,
        _query: &PreparedBatchQuery,
        witnesses: &[TieredQueryWitness],
    ) -> Result<Vec<RoutedCandidate>, Self::Error> {
        self.search(witnesses)
    }
}

/// Context attached to a routing failure, including the stable query ID.
#[derive(Debug)]
pub struct BatchRoutingError<E> {
    pub query_id: String,
    pub source: E,
}

impl<E: Display> Display for BatchRoutingError<E> {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "routing query {} failed: {}",
            self.query_id, self.source
        )
    }
}

impl<E: Error + 'static> Error for BatchRoutingError<E> {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(&self.source)
    }
}

/// Route every prepared query and attach its index to each selected candidate.
pub fn route_batch<R: BatchRouter>(
    plan: &BatchPlan,
    router: &mut R,
) -> Result<Vec<BatchCandidateWork>, BatchRoutingError<R::Error>> {
    let mut works = Vec::new();
    for query in &plan.queries {
        let witnesses = tiered_witnesses(query, &plan.scheme);
        let candidates =
            router
                .route_query(query, &witnesses)
                .map_err(|source| BatchRoutingError {
                    query_id: query.descriptor.query_id.clone(),
                    source,
                })?;
        works.extend(candidates.into_iter().map(|candidate| BatchCandidateWork {
            query_index: query.query_index,
            metagenome_id: candidate.metagenome_id.clone(),
            candidate,
        }));
    }
    Ok(works)
}

/// Group selected candidates by JMA object.  The map's key order and each
/// vector's query-index order are stable, so a caller may use the result for
/// deterministic planning even before choosing deterministic output order.
#[must_use]
pub fn group_candidate_work(
    works: &[BatchCandidateWork],
) -> BTreeMap<String, Vec<BatchCandidateWork>> {
    let mut grouped = BTreeMap::<String, Vec<BatchCandidateWork>>::new();
    for work in works {
        grouped
            .entry(work.metagenome_id.clone())
            .or_default()
            .push(work.clone());
    }
    for candidates in grouped.values_mut() {
        candidates.sort_by(|left, right| {
            left.query_index
                .cmp(&right.query_index)
                .then_with(|| left.metagenome_id.cmp(&right.metagenome_id))
        });
    }
    grouped
}

/// Exact key plus the physical page needed to inspect it.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct BatchKeyPageRequest {
    pub page_id: u64,
    pub key: WitnessKey,
}

impl BatchKeyPageRequest {
    #[must_use]
    pub const fn new(page_id: u64, key: WitnessKey) -> Self {
        Self { page_id, key }
    }
}

/// One sequence-block range requested from a selected JMA object.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct BatchSequenceBlockRequest {
    pub contig_id: u32,
    pub start: u64,
    pub end: u64,
}

impl BatchSequenceBlockRequest {
    #[must_use]
    pub const fn new(contig_id: u32, start: u64, end: u64) -> Self {
        Self {
            contig_id,
            start,
            end,
        }
    }

    #[must_use]
    pub const fn is_nonempty(self) -> bool {
        self.start < self.end
    }
}

/// Read requests emitted by one candidate-specific adapter.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct BatchCandidateReadPlan {
    pub key_pages: Vec<BatchKeyPageRequest>,
    pub sequence_blocks: Vec<BatchSequenceBlockRequest>,
}

impl BatchCandidateReadPlan {
    /// Sort and deduplicate requests without changing their logical ranges.
    pub fn normalize(&mut self) -> Result<(), BatchReadPlanError> {
        if self
            .sequence_blocks
            .iter()
            .any(|request| !request.is_nonempty())
        {
            return Err(BatchReadPlanError::EmptySequenceBlock);
        }
        self.key_pages.sort_unstable();
        self.key_pages.dedup();
        self.sequence_blocks.sort_unstable();
        self.sequence_blocks.dedup();
        Ok(())
    }
}

#[derive(Debug, Eq, PartialEq)]
pub enum BatchReadPlanError {
    EmptySequenceBlock,
}

impl Display for BatchReadPlanError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySequenceBlock => {
                formatter.write_str("sequence block request must have start < end")
            }
        }
    }
}

impl Error for BatchReadPlanError {}

/// A decoded page returned by an archive adapter.  One value is shared by all
/// candidates requesting that physical page.
#[derive(Clone, Debug)]
pub struct BatchKeyPage<P> {
    pub page_id: u64,
    pub page: P,
}

/// A decoded sequence block returned by an archive adapter.
#[derive(Clone, Debug)]
pub struct BatchSequenceBlock<B> {
    pub request: BatchSequenceBlockRequest,
    pub block: B,
}

/// Result record carrying both sides of the biological comparison key.
#[derive(Clone, Debug, PartialEq)]
pub struct BatchResultRecord<R> {
    pub key: BatchResultKey,
    pub result: R,
}

/// Optional output-order policy.  Execution order follows metagenome groups;
/// deterministic order is query ID order followed by metagenome ID.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum BatchOutputOrder {
    Execution,
    #[default]
    Deterministic,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct BatchExecutionOptions {
    pub output_order: BatchOutputOrder,
}

/// Counters proving that grouped work reuses object opens and page reads.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct BatchExecutionMetrics {
    pub query_count: u64,
    pub candidate_count: u64,
    pub metagenome_open_count: u64,
    pub exact_key_requests: u64,
    pub key_page_requests: u64,
    pub key_pages_decoded: u64,
    pub sequence_block_requests: u64,
    pub sequence_blocks_decoded: u64,
    pub records_emitted: u64,
}

/// Complete grouped execution output.
#[derive(Debug, PartialEq)]
pub struct BatchExecutionOutput<R> {
    pub records: Vec<BatchResultRecord<R>>,
    pub metrics: BatchExecutionMetrics,
}

/// JMA/alignment adapter used by [`execute_batch`].
///
/// `plan_candidate` is called once per selected query/candidate pair.  The
/// executor unions those plans per metagenome, calls `open_metagenome` once,
/// reads and decodes each shared page once, then calls `evaluate_candidate`
/// independently for every pair.  The adapter can therefore wrap the normal
/// JMA lookup and alignment path without sharing mutable mosaic state.
pub trait BatchExecutionAdapter {
    type KeyPage: Clone;
    type SequenceBlock: Clone;
    type Result;
    type Error: Error + Send + Sync + 'static;

    fn plan_candidate(
        &mut self,
        query: &PreparedBatchQuery,
        work: &BatchCandidateWork,
    ) -> Result<BatchCandidateReadPlan, Self::Error>;

    fn open_metagenome(&mut self, metagenome_id: &str) -> Result<(), Self::Error>;

    fn read_key_pages(
        &mut self,
        metagenome_id: &str,
        requests: &[BatchKeyPageRequest],
    ) -> Result<Vec<BatchKeyPage<Self::KeyPage>>, Self::Error>;

    fn read_sequence_blocks(
        &mut self,
        metagenome_id: &str,
        requests: &[BatchSequenceBlockRequest],
    ) -> Result<Vec<BatchSequenceBlock<Self::SequenceBlock>>, Self::Error>;

    fn evaluate_candidate(
        &mut self,
        query: &PreparedBatchQuery,
        work: &BatchCandidateWork,
        key_pages: &BTreeMap<u64, Self::KeyPage>,
        sequence_blocks: &BTreeMap<BatchSequenceBlockRequest, Self::SequenceBlock>,
    ) -> Result<Self::Result, Self::Error>;
}

/// Errors are contextualized with the query or metagenome responsible for a
/// failed adapter call.  Missing pages fail closed instead of causing a
/// silent broader scan.
#[derive(Debug)]
pub enum BatchExecutionError<E> {
    InvalidQueryIndex(u32),
    DuplicateResultKey(BatchResultKey),
    InvalidReadPlan {
        query_index: u32,
        source: BatchReadPlanError,
    },
    Adapter {
        metagenome_id: String,
        query_index: Option<u32>,
        source: E,
    },
    DuplicateKeyPage(u64),
    MissingKeyPage(u64),
    DuplicateSequenceBlock(BatchSequenceBlockRequest),
    MissingSequenceBlock(BatchSequenceBlockRequest),
}

impl<E: Display> Display for BatchExecutionError<E> {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidQueryIndex(index) => write!(formatter, "invalid query index {index}"),
            Self::DuplicateResultKey(key) => write!(
                formatter,
                "duplicate batch result key {}/{}",
                key.query_id, key.metagenome_id
            ),
            Self::InvalidReadPlan {
                query_index,
                source,
            } => {
                write!(
                    formatter,
                    "invalid read plan for query {query_index}: {source}"
                )
            }
            Self::Adapter {
                metagenome_id,
                query_index,
                source,
            } => write!(
                formatter,
                "batch adapter failed for metagenome {metagenome_id}, query {query_index:?}: {source}"
            ),
            Self::DuplicateKeyPage(page_id) => {
                write!(formatter, "duplicate decoded key page {page_id}")
            }
            Self::MissingKeyPage(page_id) => {
                write!(formatter, "missing decoded key page {page_id}")
            }
            Self::DuplicateSequenceBlock(request) => write!(
                formatter,
                "duplicate decoded sequence block contig={} range={}:{}",
                request.contig_id, request.start, request.end
            ),
            Self::MissingSequenceBlock(request) => write!(
                formatter,
                "missing decoded sequence block contig={} range={}:{}",
                request.contig_id, request.start, request.end
            ),
        }
    }
}

impl<E: Error + 'static> Error for BatchExecutionError<E> {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::InvalidReadPlan { source, .. } => Some(source),
            Self::Adapter { source, .. } => Some(source),
            _ => None,
        }
    }
}

struct PlannedCandidate {
    work: BatchCandidateWork,
    read_plan: BatchCandidateReadPlan,
}

/// Execute selected work grouped by metagenome.
pub fn execute_batch<A: BatchExecutionAdapter>(
    plan: &BatchPlan,
    works: &[BatchCandidateWork],
    adapter: &mut A,
    options: BatchExecutionOptions,
) -> Result<BatchExecutionOutput<A::Result>, BatchExecutionError<A::Error>> {
    let mut seen_result_keys = BTreeSet::new();
    let mut grouped = BTreeMap::<String, Vec<PlannedCandidate>>::new();
    for work in works {
        let query = plan
            .query(work.query_index)
            .ok_or(BatchExecutionError::InvalidQueryIndex(work.query_index))?;
        if work.metagenome_id != work.candidate.metagenome_id {
            return Err(BatchExecutionError::InvalidQueryIndex(work.query_index));
        }
        let result_key = BatchResultKey {
            query_id: query.descriptor.query_id.clone(),
            metagenome_id: work.metagenome_id.clone(),
        };
        if !seen_result_keys.insert(result_key) {
            return Err(BatchExecutionError::DuplicateResultKey(BatchResultKey {
                query_id: query.descriptor.query_id.clone(),
                metagenome_id: work.metagenome_id.clone(),
            }));
        }
        let mut read_plan =
            adapter
                .plan_candidate(query, work)
                .map_err(|source| BatchExecutionError::Adapter {
                    metagenome_id: work.metagenome_id.clone(),
                    query_index: Some(work.query_index),
                    source,
                })?;
        read_plan
            .normalize()
            .map_err(|source| BatchExecutionError::InvalidReadPlan {
                query_index: work.query_index,
                source,
            })?;
        grouped
            .entry(work.metagenome_id.clone())
            .or_default()
            .push(PlannedCandidate {
                work: work.clone(),
                read_plan,
            });
    }

    let mut records = Vec::with_capacity(works.len());
    let mut metrics = BatchExecutionMetrics {
        query_count: u64::try_from(plan.queries.len()).unwrap_or(u64::MAX),
        candidate_count: u64::try_from(works.len()).unwrap_or(u64::MAX),
        ..BatchExecutionMetrics::default()
    };

    for (metagenome_id, candidates) in &mut grouped {
        candidates.sort_by_key(|candidate| candidate.work.query_index);
        adapter
            .open_metagenome(metagenome_id)
            .map_err(|source| BatchExecutionError::Adapter {
                metagenome_id: metagenome_id.clone(),
                query_index: None,
                source,
            })?;
        metrics.metagenome_open_count = metrics.metagenome_open_count.saturating_add(1);

        let mut key_requests = BTreeSet::new();
        let mut sequence_requests = BTreeSet::new();
        for candidate in candidates.iter() {
            key_requests.extend(candidate.read_plan.key_pages.iter().copied());
            sequence_requests.extend(candidate.read_plan.sequence_blocks.iter().copied());
        }
        let key_requests = key_requests.into_iter().collect::<Vec<_>>();
        let sequence_requests = sequence_requests.into_iter().collect::<Vec<_>>();
        metrics.exact_key_requests = metrics
            .exact_key_requests
            .saturating_add(u64::try_from(key_requests.len()).unwrap_or(u64::MAX));
        metrics.sequence_block_requests = metrics
            .sequence_block_requests
            .saturating_add(u64::try_from(sequence_requests.len()).unwrap_or(u64::MAX));

        let key_pages = if key_requests.is_empty() {
            Vec::new()
        } else {
            adapter
                .read_key_pages(metagenome_id, &key_requests)
                .map_err(|source| BatchExecutionError::Adapter {
                    metagenome_id: metagenome_id.clone(),
                    query_index: None,
                    source,
                })?
        };
        let mut key_page_map = BTreeMap::new();
        for page in key_pages {
            if key_page_map.insert(page.page_id, page.page).is_some() {
                return Err(BatchExecutionError::DuplicateKeyPage(page.page_id));
            }
        }
        let unique_key_page_count = key_requests
            .iter()
            .map(|request| request.page_id)
            .collect::<BTreeSet<_>>()
            .len();
        metrics.key_page_requests = metrics
            .key_page_requests
            .saturating_add(u64::try_from(unique_key_page_count).unwrap_or(u64::MAX));
        metrics.key_pages_decoded = metrics
            .key_pages_decoded
            .saturating_add(u64::try_from(key_page_map.len()).unwrap_or(u64::MAX));
        for page_id in key_requests
            .iter()
            .map(|request| request.page_id)
            .collect::<BTreeSet<_>>()
        {
            if !key_page_map.contains_key(&page_id) {
                return Err(BatchExecutionError::MissingKeyPage(page_id));
            }
        }

        let sequence_blocks = if sequence_requests.is_empty() {
            Vec::new()
        } else {
            adapter
                .read_sequence_blocks(metagenome_id, &sequence_requests)
                .map_err(|source| BatchExecutionError::Adapter {
                    metagenome_id: metagenome_id.clone(),
                    query_index: None,
                    source,
                })?
        };
        let mut sequence_block_map = BTreeMap::new();
        for block in sequence_blocks {
            if sequence_block_map
                .insert(block.request, block.block)
                .is_some()
            {
                return Err(BatchExecutionError::DuplicateSequenceBlock(block.request));
            }
        }
        metrics.sequence_blocks_decoded = metrics
            .sequence_blocks_decoded
            .saturating_add(u64::try_from(sequence_block_map.len()).unwrap_or(u64::MAX));
        for request in &sequence_requests {
            if !sequence_block_map.contains_key(request) {
                return Err(BatchExecutionError::MissingSequenceBlock(*request));
            }
        }

        for candidate in candidates.iter() {
            let query = plan.query(candidate.work.query_index).ok_or(
                BatchExecutionError::InvalidQueryIndex(candidate.work.query_index),
            )?;
            let result = adapter
                .evaluate_candidate(query, &candidate.work, &key_page_map, &sequence_block_map)
                .map_err(|source| BatchExecutionError::Adapter {
                    metagenome_id: metagenome_id.clone(),
                    query_index: Some(candidate.work.query_index),
                    source,
                })?;
            records.push(BatchResultRecord {
                key: BatchResultKey {
                    query_id: query.descriptor.query_id.clone(),
                    metagenome_id: metagenome_id.clone(),
                },
                result,
            });
        }
    }

    if matches!(options.output_order, BatchOutputOrder::Deterministic) {
        records.sort_by(|left, right| left.key.cmp(&right.key));
    }
    metrics.records_emitted = u64::try_from(records.len()).unwrap_or(u64::MAX);
    Ok(BatchExecutionOutput { records, metrics })
}
