//! Runtime adapter for grouped access to public [`TraceArchive`] objects.
//!
//! The batch executor owns grouping and independent result dispatch.  This
//! adapter supplies the concrete JMA side: it opens one archive per selected
//! metagenome, translates exact router witnesses to public [`SeedKey`] values,
//! performs one bounded seed lookup per unioned physical page, and coalesces
//! sequence-block reads through [`TraceArchive::read_sequences`].

use super::{
    BatchCandidateReadPlan, BatchCandidateWork, BatchExecutionAdapter, BatchKeyPage,
    BatchKeyPageRequest, BatchSequenceBlock, BatchSequenceBlockRequest,
};
use crate::archive::{
    ArchiveError, SeedKey, SeedLookupResult, SeedMatch, SeedOccurrence, SeedSchemeId,
    SequenceRequest, TraceArchive,
};
use crate::resource::ResourceOpenOptions;
use crate::router::WitnessKey;
use crate::trace::catalog::TraceCatalog;
use crate::trace::contracts::PreparedBatchQuery;
use crate::trace::raw::{AssemblyResource, RawError, open_resource};
use std::collections::{BTreeMap, BTreeSet};
use std::error::Error;
use std::fmt::{self, Display, Formatter};

/// One exact JMA key match decoded from a public archive lookup.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RuntimeKeyMatch {
    pub key: WitnessKey,
    pub occurrences: Vec<SeedOccurrence>,
}

/// Decoded seed results for one physical key page.  The page is shared by all
/// candidates in the metagenome group; exact candidate selection happens in
/// [`JmaBatchRuntime::evaluate_candidate`].
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct RuntimeKeyPage {
    pub requests: Vec<BatchKeyPageRequest>,
    pub matches: Vec<RuntimeKeyMatch>,
}

/// Sequence bytes returned for one half-open block request.
pub type RuntimeSequenceBlock = Vec<u8>;

/// Per-query evidence produced by the runtime resource stage.
///
/// This is intentionally not a biological presence call and does not replace
/// the normal alignment/mosaic record.  The trace runner can use the exact
/// matches and decoded blocks as the input to its existing per-query alignment
/// path while retaining this result key.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RuntimeTraceResult {
    pub query_id: String,
    pub metagenome_id: String,
    pub exact_key_matches: u32,
    pub occurrence_count: u64,
    pub sequence_bases: u64,
}

/// Errors crossing the runtime adapter boundary.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum RuntimeError {
    InvalidSchemeId(u32),
    MissingArchive(String),
    InvalidWitness(String),
    InvalidSequenceRequest(String),
    Factory(String),
    Archive(String),
}

impl Display for RuntimeError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidSchemeId(id) => write!(formatter, "seed scheme ID must be non-zero: {id}"),
            Self::MissingArchive(metagenome_id) => {
                write!(formatter, "metagenome archive is not open: {metagenome_id}")
            }
            Self::InvalidWitness(message) => write!(formatter, "invalid witness: {message}"),
            Self::InvalidSequenceRequest(message) => {
                write!(formatter, "invalid sequence request: {message}")
            }
            Self::Factory(message) => write!(formatter, "archive factory failed: {message}"),
            Self::Archive(message) => write!(formatter, "archive operation failed: {message}"),
        }
    }
}

impl Error for RuntimeError {}

/// Factory used by [`JmaBatchRuntime`] to open a selected metagenome archive.
pub trait BatchArchiveFactory {
    type Archive: TraceArchive;

    fn open(&mut self, metagenome_id: &str) -> Result<Self::Archive, RuntimeError>;
}

/// Native catalog-backed factory.  The expected catalog checksum is applied to
/// every open and full-object fallback remains disabled by default.
#[derive(Clone, Debug)]
pub struct NativeJmaFactory {
    catalog: TraceCatalog,
    options: ResourceOpenOptions,
}

impl NativeJmaFactory {
    #[must_use]
    pub fn new(catalog: TraceCatalog, mut options: ResourceOpenOptions) -> Self {
        options.allow_full_download_fallback = false;
        Self { catalog, options }
    }

    #[must_use]
    pub fn catalog(&self) -> &TraceCatalog {
        &self.catalog
    }
}

impl BatchArchiveFactory for NativeJmaFactory {
    type Archive = crate::archive::NativeJmaArchive<AssemblyResource>;

    fn open(&mut self, metagenome_id: &str) -> Result<Self::Archive, RuntimeError> {
        let entry = self.catalog.get(metagenome_id).ok_or_else(|| {
            RuntimeError::Factory(format!("catalog entry is missing: {metagenome_id}"))
        })?;
        let mut options = self.options.clone();
        options.expected_sha256 = Some(entry.sha256.clone());
        let resource = open_resource(entry.resource(), options)
            .map_err(|error: RawError| RuntimeError::Factory(error.to_string()))?;
        crate::archive::NativeJmaArchive::from_resource(resource)
            .map_err(|error| RuntimeError::Factory(error.to_string()))
    }
}

/// Grouped JMA runtime.  The archive map is deliberately retained for the
/// lifetime of the runtime, so repeated batch calls may reuse an already-open
/// resource while one `execute_batch` call still opens each group once.
pub struct JmaBatchRuntime<F: BatchArchiveFactory> {
    factory: F,
    scheme: SeedSchemeId,
    max_occurrences: Option<u32>,
    archives: BTreeMap<String, F::Archive>,
    candidate_plans: BTreeMap<(u32, String), BatchCandidateReadPlan>,
}

impl<F: BatchArchiveFactory> JmaBatchRuntime<F> {
    pub fn new(
        factory: F,
        scheme: SeedSchemeId,
        max_occurrences: Option<u32>,
    ) -> Result<Self, RuntimeError> {
        if scheme.0 == 0 {
            return Err(RuntimeError::InvalidSchemeId(scheme.0));
        }
        Ok(Self {
            factory,
            scheme,
            max_occurrences,
            archives: BTreeMap::new(),
            candidate_plans: BTreeMap::new(),
        })
    }

    #[must_use]
    pub fn scheme(&self) -> SeedSchemeId {
        self.scheme
    }

    /// Override the conservative default plan for one selected candidate.
    /// This is the hook for the runner to install chain-derived sequence
    /// ranges while retaining the same grouped archive reads.
    pub fn set_candidate_plan(
        &mut self,
        query_index: u32,
        metagenome_id: impl Into<String>,
        plan: BatchCandidateReadPlan,
    ) {
        self.candidate_plans
            .insert((query_index, metagenome_id.into()), plan);
    }

    #[must_use]
    pub fn open_archive_count(&self) -> usize {
        self.archives.len()
    }

    fn archive(&self, metagenome_id: &str) -> Result<&F::Archive, RuntimeError> {
        self.archives
            .get(metagenome_id)
            .ok_or_else(|| RuntimeError::MissingArchive(metagenome_id.to_string()))
    }

    fn default_plan(&self, work: &BatchCandidateWork) -> BatchCandidateReadPlan {
        let mut plan = BatchCandidateReadPlan::default();
        let mut seen_keys = BTreeSet::new();
        for shared in &work.candidate.shared_witnesses {
            if seen_keys.insert(shared.key) {
                plan.key_pages.push(BatchKeyPageRequest::new(
                    page_id_for_key(shared.key),
                    shared.key,
                ));
            }
        }
        let mut seen_blocks = BTreeSet::new();
        for occurrence in &work.candidate.positional_witnesses {
            let start = occurrence.position;
            let Some(end) = start.checked_add(u64::from(crate::router::WITNESS_K)) else {
                continue;
            };
            let request = BatchSequenceBlockRequest::new(occurrence.contig_id, start, end);
            if request.is_nonempty() && seen_blocks.insert(request) {
                plan.sequence_blocks.push(request);
            }
        }
        plan
    }

    fn plan_for(&self, work: &BatchCandidateWork) -> Result<BatchCandidateReadPlan, RuntimeError> {
        Ok(self
            .candidate_plans
            .get(&(work.query_index, work.metagenome_id.clone()))
            .cloned()
            .unwrap_or_else(|| self.default_plan(work)))
    }
}

impl<F: BatchArchiveFactory> BatchExecutionAdapter for JmaBatchRuntime<F> {
    type KeyPage = RuntimeKeyPage;
    type SequenceBlock = RuntimeSequenceBlock;
    type Result = RuntimeTraceResult;
    type Error = RuntimeError;

    fn plan_candidate(
        &mut self,
        _query: &PreparedBatchQuery,
        work: &BatchCandidateWork,
    ) -> Result<BatchCandidateReadPlan, Self::Error> {
        self.plan_for(work)
    }

    fn open_metagenome(&mut self, metagenome_id: &str) -> Result<(), Self::Error> {
        if !self.archives.contains_key(metagenome_id) {
            let archive = self.factory.open(metagenome_id)?;
            self.archives.insert(metagenome_id.to_string(), archive);
        }
        Ok(())
    }

    fn read_key_pages(
        &mut self,
        metagenome_id: &str,
        requests: &[BatchKeyPageRequest],
    ) -> Result<Vec<BatchKeyPage<Self::KeyPage>>, Self::Error> {
        let archive = self.archive(metagenome_id)?;
        let mut grouped = BTreeMap::<u64, Vec<BatchKeyPageRequest>>::new();
        for request in requests {
            grouped.entry(request.page_id).or_default().push(*request);
        }
        let mut key_pages = BTreeMap::<WitnessKey, BTreeSet<u64>>::new();
        for (page_id, page_requests) in &grouped {
            for request in page_requests {
                key_pages.entry(request.key).or_default().insert(*page_id);
            }
        }
        let unique_keys = key_pages.keys().copied().collect::<Vec<_>>();
        let keys = unique_keys
            .iter()
            .copied()
            .map(seed_key)
            .collect::<Result<Vec<_>, _>>()?;
        let lookup = archive
            .lookup_seeds_bounded(self.scheme, &keys, self.max_occurrences)
            .map_err(archive_error)?;
        let matches_by_page = runtime_matches(&unique_keys, &key_pages, &keys, lookup)?;
        let mut pages = Vec::with_capacity(grouped.len());
        for (page_id, page_requests) in grouped {
            let matches = matches_by_page.get(&page_id).cloned().unwrap_or_default();
            pages.push(BatchKeyPage {
                page_id,
                page: RuntimeKeyPage {
                    requests: page_requests,
                    matches,
                },
            });
        }
        Ok(pages)
    }

    fn read_sequence_blocks(
        &mut self,
        metagenome_id: &str,
        requests: &[BatchSequenceBlockRequest],
    ) -> Result<Vec<BatchSequenceBlock<Self::SequenceBlock>>, Self::Error> {
        let archive = self.archive(metagenome_id)?;
        let archive_requests = requests
            .iter()
            .map(|request| {
                SequenceRequest::new(request.contig_id, request.start, request.end, false)
                    .map_err(|error| RuntimeError::InvalidSequenceRequest(error.to_string()))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let slices = archive
            .read_sequences(&archive_requests)
            .map_err(archive_error)?;
        slices
            .into_iter()
            .map(|slice| {
                let request = BatchSequenceBlockRequest::new(
                    slice.request.contig_id,
                    slice.request.start,
                    slice.request.end,
                );
                Ok(BatchSequenceBlock {
                    request,
                    block: slice.bases,
                })
            })
            .collect()
    }

    fn evaluate_candidate(
        &mut self,
        query: &PreparedBatchQuery,
        work: &BatchCandidateWork,
        key_pages: &BTreeMap<u64, Self::KeyPage>,
        sequence_blocks: &BTreeMap<BatchSequenceBlockRequest, Self::SequenceBlock>,
    ) -> Result<Self::Result, Self::Error> {
        let read_plan = self.plan_for(work)?;
        let candidate_keys = work
            .candidate
            .shared_witnesses
            .iter()
            .map(|witness| witness.key)
            .collect::<BTreeSet<_>>();
        let mut exact_key_matches = 0_u32;
        let mut occurrence_count = 0_u64;
        for request in &read_plan.key_pages {
            if !candidate_keys.contains(&request.key) {
                continue;
            }
            let Some(page) = key_pages.get(&request.page_id) else {
                continue;
            };
            if let Some(matched) = page
                .matches
                .iter()
                .find(|matched| matched.key == request.key)
            {
                exact_key_matches = exact_key_matches.saturating_add(1);
                occurrence_count = occurrence_count
                    .saturating_add(u64::try_from(matched.occurrences.len()).unwrap_or(u64::MAX));
            }
        }
        let sequence_bases = read_plan
            .sequence_blocks
            .iter()
            .filter_map(|request| sequence_blocks.get(request))
            .map(|block| u64::try_from(block.len()).unwrap_or(u64::MAX))
            .fold(0_u64, u64::saturating_add);
        Ok(RuntimeTraceResult {
            query_id: query.descriptor.query_id.clone(),
            metagenome_id: work.metagenome_id.clone(),
            exact_key_matches,
            occurrence_count,
            sequence_bases,
        })
    }
}

fn page_id_for_key(key: WitnessKey) -> u64 {
    key.jamhash >> 56
}

fn seed_key(key: WitnessKey) -> Result<SeedKey, RuntimeError> {
    let bytes = key.packed.to_be_bytes();
    let verification = bytes[2..].to_vec();
    WitnessKey::checked(key.packed, key.jamhash)
        .map_err(|error| RuntimeError::InvalidWitness(error.to_string()))?;
    Ok(SeedKey {
        digest: key.jamhash,
        verification,
    })
}

fn runtime_matches(
    unique_keys: &[WitnessKey],
    key_pages: &BTreeMap<WitnessKey, BTreeSet<u64>>,
    keys: &[SeedKey],
    lookup: SeedLookupResult,
) -> Result<BTreeMap<u64, Vec<RuntimeKeyMatch>>, RuntimeError> {
    let mut matches = BTreeMap::<u64, Vec<RuntimeKeyMatch>>::new();
    for SeedMatch {
        key_index,
        occurrences,
    } in lookup.matches
    {
        let key = unique_keys.get(key_index).copied().ok_or_else(|| {
            RuntimeError::Archive("seed match key index is out of bounds".to_string())
        })?;
        let expected = keys
            .get(key_index)
            .ok_or_else(|| RuntimeError::Archive("seed query key is unavailable".to_string()))?;
        if expected.digest != key.jamhash {
            return Err(RuntimeError::Archive(
                "seed match digest mismatch".to_string(),
            ));
        }
        let pages = key_pages.get(&key).ok_or_else(|| {
            RuntimeError::Archive("seed match key has no requested page".to_string())
        })?;
        for page_id in pages {
            matches.entry(*page_id).or_default().push(RuntimeKeyMatch {
                key,
                occurrences: occurrences.clone(),
            });
        }
    }
    Ok(matches)
}

fn archive_error(error: ArchiveError) -> RuntimeError {
    RuntimeError::Archive(error.to_string())
}
