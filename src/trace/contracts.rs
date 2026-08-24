//! Shared query-planning, batch, and pre-allocation memory contracts.

use super::model::QueryDescriptor;
use crate::router::{QueryWitness, RoutedCandidate};
use serde::{Deserialize, Serialize};
use thiserror::Error;

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct WitnessPlanningRequest {
    pub min_trace_length: u64,
    pub target_identity: f64,
    pub max_zero_witness_risk: f64,
    pub strict: bool,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct WitnessTierPlan {
    pub model_assumptions: Vec<String>,
    pub requested_min_trace_length: u64,
    pub requested_identity: f64,
    pub selected_witness_tier: Option<u32>,
    pub estimated_zero_witness_risk: f64,
    pub available_denser_tiers: Vec<u32>,
    pub warning: Option<String>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PreparedBatchQuery {
    pub query_index: u32,
    pub descriptor: QueryDescriptor,
    pub sequence: Vec<u8>,
    pub witnesses: Vec<QueryWitness>,
    pub witness_plan: WitnessTierPlan,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BatchCandidateWork {
    pub query_index: u32,
    pub metagenome_id: String,
    pub candidate: RoutedCandidate,
}

#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct BatchResultKey {
    pub query_id: String,
    pub metagenome_id: String,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub struct TraceMemoryEstimate {
    pub candidate_accumulator_bytes: u64,
    pub router_key_page_bytes: u64,
    pub router_posting_bytes: u64,
    pub jma_key_page_bytes: u64,
    pub jma_occurrence_bytes: u64,
    pub anchor_bytes: u64,
    pub chain_bytes: u64,
    pub sequence_block_bytes: u64,
    pub alignment_score_bytes: u64,
    pub traceback_bytes: u64,
    pub output_bytes: u64,
}

impl TraceMemoryEstimate {
    #[must_use]
    pub fn total_bytes(self) -> u64 {
        [
            self.candidate_accumulator_bytes,
            self.router_key_page_bytes,
            self.router_posting_bytes,
            self.jma_key_page_bytes,
            self.jma_occurrence_bytes,
            self.anchor_bytes,
            self.chain_bytes,
            self.sequence_block_bytes,
            self.alignment_score_bytes,
            self.traceback_bytes,
            self.output_bytes,
        ]
        .into_iter()
        .fold(0_u64, u64::saturating_add)
    }
}

/// An acquired value must release its bytes when it is dropped.
pub trait CandidateMemoryBudget: Send + Sync {
    type Reservation: Send;

    fn acquire(
        &self,
        estimate: TraceMemoryEstimate,
    ) -> Result<Self::Reservation, MemoryBudgetError>;
}

#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum MemoryBudgetError {
    #[error("candidate requires {requested_bytes} bytes, exceeding budget {budget_bytes}")]
    CandidateTooLarge {
        requested_bytes: u64,
        budget_bytes: u64,
    },
    #[error("memory budget is closed")]
    Closed,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn memory_estimate_accounts_before_candidate_work() {
        let estimate = TraceMemoryEstimate {
            candidate_accumulator_bytes: 1,
            router_key_page_bytes: 2,
            router_posting_bytes: 4,
            jma_key_page_bytes: 8,
            jma_occurrence_bytes: 16,
            anchor_bytes: 32,
            chain_bytes: 64,
            sequence_block_bytes: 128,
            alignment_score_bytes: 256,
            traceback_bytes: 512,
            output_bytes: 1_024,
        };
        assert_eq!(estimate.total_bytes(), 2_047);
    }
}
