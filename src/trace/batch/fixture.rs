//! Small in-memory archive/alignment adapter for batch contract tests.
//!
//! The fixture intentionally models only the shared-resource boundary.  Its
//! result is a deterministic trace signature, not a biological presence
//! call.  It is useful for proving that grouping does not change per-query
//! results while physical JMA opens, key pages, and sequence blocks are
//! shared.

use super::{
    BatchCandidateReadPlan, BatchCandidateWork, BatchExecutionAdapter, BatchKeyPage,
    BatchSequenceBlock, BatchSequenceBlockRequest,
};
use crate::router::WitnessKey;
use crate::trace::contracts::PreparedBatchQuery;
use std::collections::{BTreeMap, BTreeSet};
use std::error::Error;
use std::fmt::{self, Display, Formatter};

/// Decoded fixture key page.  Exact packed keys remain visible so tests can
/// prove that page reuse does not turn hash equality into evidence.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct FixtureKeyPage {
    pub keys: BTreeSet<WitnessKey>,
}

/// Decoded fixture sequence block.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FixtureSequenceBlock {
    pub bases: Vec<u8>,
}

/// Deterministic stand-in for one per-query alignment/mosaic result.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FixtureTraceResult {
    pub query_id: String,
    pub metagenome_id: String,
    pub exact_key_hits: u32,
    pub sequence_bases: u64,
    pub biological_signature: u64,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum FixtureError {
    UnknownCandidate {
        query_index: u32,
        metagenome_id: String,
    },
    UnknownKeyPage {
        metagenome_id: String,
        page_id: u64,
    },
    UnknownSequenceBlock {
        metagenome_id: String,
        request: BatchSequenceBlockRequest,
    },
    WrongMetagenome {
        expected: String,
        actual: String,
    },
}

impl Display for FixtureError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownCandidate {
                query_index,
                metagenome_id,
            } => write!(
                formatter,
                "fixture has no plan for query {query_index} and metagenome {metagenome_id}"
            ),
            Self::UnknownKeyPage {
                metagenome_id,
                page_id,
            } => write!(
                formatter,
                "fixture has no key page {page_id} for {metagenome_id}"
            ),
            Self::UnknownSequenceBlock {
                metagenome_id,
                request,
            } => write!(
                formatter,
                "fixture has no sequence block {}:{}-{} for {}",
                request.contig_id, request.start, request.end, metagenome_id
            ),
            Self::WrongMetagenome { expected, actual } => {
                write!(formatter, "fixture opened {actual}, expected {expected}")
            }
        }
    }
}

impl Error for FixtureError {}

/// Observable read/open counters from [`InMemoryBatchAdapter`].
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct FixtureMetrics {
    pub opened_metagenomes: Vec<String>,
    pub key_page_reads: Vec<(String, u64)>,
    pub sequence_block_reads: Vec<(String, BatchSequenceBlockRequest)>,
}

/// In-memory adapter with explicit candidate plans and object pages.
#[derive(Clone, Debug, Default)]
pub struct InMemoryBatchAdapter {
    candidate_plans: BTreeMap<(u32, String), BatchCandidateReadPlan>,
    key_pages: BTreeMap<(String, u64), FixtureKeyPage>,
    sequence_blocks: BTreeMap<(String, BatchSequenceBlockRequest), FixtureSequenceBlock>,
    metrics: FixtureMetrics,
}

impl InMemoryBatchAdapter {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    pub fn insert_candidate_plan(
        &mut self,
        query_index: u32,
        metagenome_id: impl Into<String>,
        read_plan: BatchCandidateReadPlan,
    ) {
        self.candidate_plans
            .insert((query_index, metagenome_id.into()), read_plan);
    }

    pub fn insert_key_page(
        &mut self,
        metagenome_id: impl Into<String>,
        page_id: u64,
        keys: impl IntoIterator<Item = WitnessKey>,
    ) {
        self.key_pages.insert(
            (metagenome_id.into(), page_id),
            FixtureKeyPage {
                keys: keys.into_iter().collect(),
            },
        );
    }

    pub fn insert_sequence_block(
        &mut self,
        metagenome_id: impl Into<String>,
        request: BatchSequenceBlockRequest,
        bases: Vec<u8>,
    ) {
        self.sequence_blocks.insert(
            (metagenome_id.into(), request),
            FixtureSequenceBlock { bases },
        );
    }

    #[must_use]
    pub fn metrics(&self) -> &FixtureMetrics {
        &self.metrics
    }

    /// Reset only physical-I/O observations.  Candidate plans and object data
    /// remain reusable for a separate-query comparison.
    pub fn reset_metrics(&mut self) {
        self.metrics = FixtureMetrics::default();
    }

    fn plan_for(
        &self,
        query_index: u32,
        metagenome_id: &str,
    ) -> Result<BatchCandidateReadPlan, FixtureError> {
        self.candidate_plans
            .get(&(query_index, metagenome_id.to_string()))
            .cloned()
            .ok_or_else(|| FixtureError::UnknownCandidate {
                query_index,
                metagenome_id: metagenome_id.to_string(),
            })
    }
}

impl BatchExecutionAdapter for InMemoryBatchAdapter {
    type KeyPage = FixtureKeyPage;
    type SequenceBlock = FixtureSequenceBlock;
    type Result = FixtureTraceResult;
    type Error = FixtureError;

    fn plan_candidate(
        &mut self,
        _query: &PreparedBatchQuery,
        work: &BatchCandidateWork,
    ) -> Result<BatchCandidateReadPlan, Self::Error> {
        self.plan_for(work.query_index, &work.metagenome_id)
    }

    fn open_metagenome(&mut self, metagenome_id: &str) -> Result<(), Self::Error> {
        self.metrics
            .opened_metagenomes
            .push(metagenome_id.to_string());
        Ok(())
    }

    fn read_key_pages(
        &mut self,
        metagenome_id: &str,
        requests: &[super::BatchKeyPageRequest],
    ) -> Result<Vec<BatchKeyPage<Self::KeyPage>>, Self::Error> {
        let page_ids = requests
            .iter()
            .map(|request| request.page_id)
            .collect::<BTreeSet<_>>();
        let mut pages = Vec::with_capacity(page_ids.len());
        for page_id in page_ids {
            self.metrics
                .key_page_reads
                .push((metagenome_id.to_string(), page_id));
            let page = self
                .key_pages
                .get(&(metagenome_id.to_string(), page_id))
                .cloned()
                .ok_or_else(|| FixtureError::UnknownKeyPage {
                    metagenome_id: metagenome_id.to_string(),
                    page_id,
                })?;
            pages.push(BatchKeyPage { page_id, page });
        }
        Ok(pages)
    }

    fn read_sequence_blocks(
        &mut self,
        metagenome_id: &str,
        requests: &[BatchSequenceBlockRequest],
    ) -> Result<Vec<BatchSequenceBlock<Self::SequenceBlock>>, Self::Error> {
        let mut blocks = Vec::with_capacity(requests.len());
        for &request in requests {
            self.metrics
                .sequence_block_reads
                .push((metagenome_id.to_string(), request));
            let block = self
                .sequence_blocks
                .get(&(metagenome_id.to_string(), request))
                .cloned()
                .ok_or_else(|| FixtureError::UnknownSequenceBlock {
                    metagenome_id: metagenome_id.to_string(),
                    request,
                })?;
            blocks.push(BatchSequenceBlock { request, block });
        }
        Ok(blocks)
    }

    fn evaluate_candidate(
        &mut self,
        query: &PreparedBatchQuery,
        work: &BatchCandidateWork,
        key_pages: &BTreeMap<u64, Self::KeyPage>,
        sequence_blocks: &BTreeMap<BatchSequenceBlockRequest, Self::SequenceBlock>,
    ) -> Result<Self::Result, Self::Error> {
        let read_plan = self.plan_for(work.query_index, &work.metagenome_id)?;
        let exact_key_hits = read_plan
            .key_pages
            .iter()
            .filter(|request| {
                key_pages
                    .get(&request.page_id)
                    .is_some_and(|page| page.keys.contains(&request.key))
            })
            .count();
        let sequence_bases = read_plan
            .sequence_blocks
            .iter()
            .filter_map(|request| sequence_blocks.get(request))
            .map(|block| u64::try_from(block.bases.len()).unwrap_or(u64::MAX))
            .fold(0_u64, u64::saturating_add);
        // This signature stands in for a full alignment/mosaic result.  It is
        // intentionally derived only from the candidate's own requests, so a
        // unioned page cannot leak evidence from a different query.
        let biological_signature = u64::try_from(exact_key_hits)
            .unwrap_or(u64::MAX)
            .saturating_mul(1_000_003)
            .saturating_add(sequence_bases);
        Ok(FixtureTraceResult {
            query_id: query.descriptor.query_id.clone(),
            metagenome_id: work.metagenome_id.clone(),
            exact_key_hits: u32::try_from(exact_key_hits).unwrap_or(u32::MAX),
            sequence_bases,
            biological_signature,
        })
    }
}
