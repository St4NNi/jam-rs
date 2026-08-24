//! Checked planning for logical byte ranges and physical range reads.
//!
//! The planner is intentionally transport-neutral.  Callers mark ranges
//! already served by a page cache, receive deterministic physical ranges for
//! the remaining requests, and can map every returned byte back to its
//! logical request.  A physical range is formed only when all configured gap,
//! over-read, and size limits remain satisfied.

use super::ByteRange;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use thiserror::Error;

/// Limits applied while coalescing remote or local range requests.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct RangePlanLimits {
    /// Largest unrequested gap that may be fetched between two ranges.
    pub max_gap_bytes: u64,
    /// Largest unrequested portion allowed in one physical range.
    pub max_overread_bytes: u64,
    /// Largest physical range request.
    pub max_physical_range_bytes: u64,
    /// Maximum number of physical requests in one deterministic I/O batch.
    pub max_concurrency: usize,
}

impl Default for RangePlanLimits {
    fn default() -> Self {
        Self {
            max_gap_bytes: 4096,
            max_overread_bytes: 8192,
            max_physical_range_bytes: 1024 * 1024,
            max_concurrency: 4,
        }
    }
}

impl RangePlanLimits {
    fn validate(self) -> Result<Self, ReadPlanError> {
        if self.max_physical_range_bytes == 0 {
            return Err(ReadPlanError::InvalidLimit(
                "max_physical_range_bytes must be greater than zero",
            ));
        }
        if self.max_concurrency == 0 {
            return Err(ReadPlanError::InvalidLimit(
                "max_concurrency must be greater than zero",
            ));
        }
        Ok(self)
    }
}

/// A logical read request.  `cache_hit` means the caller already has the
/// exact bytes and therefore this request contributes no physical I/O.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct RangeRequest {
    pub id: u64,
    pub range: ByteRange,
    pub cache_hit: bool,
}

impl RangeRequest {
    #[must_use]
    pub const fn new(id: u64, range: ByteRange) -> Self {
        Self {
            id,
            range,
            cache_hit: false,
        }
    }

    #[must_use]
    pub const fn cached(id: u64, range: ByteRange) -> Self {
        Self {
            id,
            range,
            cache_hit: true,
        }
    }

    #[must_use]
    pub const fn with_cache_hit(mut self, cache_hit: bool) -> Self {
        self.cache_hit = cache_hit;
        self
    }
}

/// Errors found before any transport request is allowed.
#[derive(Debug, Error, Eq, PartialEq)]
pub enum ReadPlanError {
    #[error("logical range {request_id} overflows")]
    RangeOverflow { request_id: u64 },
    #[error("logical range {request_id} exceeds the physical range cap")]
    RangeExceedsPhysicalCap { request_id: u64 },
    #[error("duplicate logical request id {0}")]
    DuplicateRequestId(u64),
    #[error("invalid range planner limit: {0}")]
    InvalidLimit(&'static str),
    #[error("range accounting overflow")]
    AccountingOverflow,
    #[error("physical response {0} is not part of this plan")]
    UnknownPhysicalIndex(usize),
    #[error("physical response {0} was supplied more than once")]
    DuplicatePhysicalResponse(usize),
    #[error("physical response {0} is shorter than its useful decoded bytes")]
    ResponseShorterThanUseful(usize),
}

/// One physical range and all logical requests it covers.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct PhysicalRange {
    pub index: usize,
    pub range: ByteRange,
    pub logical_request_ids: Vec<u64>,
    pub useful_logical_bytes: u64,
    pub overread_bytes: u64,
}

/// Mapping from one logical request into a physical response.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct LogicalMapping {
    pub request_id: u64,
    pub logical_range: ByteRange,
    pub physical_index: Option<usize>,
    pub physical_offset: Option<u64>,
    pub cache_hit: bool,
}

/// A deterministic batch of physical indexes.  The indexes are in ascending
/// plan order; callers may execute each batch concurrently and batches in
/// order when deterministic I/O is required.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct ReadBatch {
    pub ordinal: usize,
    pub physical_indices: Vec<usize>,
}

/// Bytes observed after a physical request completes.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PhysicalReadObservation {
    pub physical_index: usize,
    pub returned_bytes: u64,
    pub useful_decoded_bytes: u64,
}

/// Range I/O and cache accounting.  `physical_returned_bytes` is initialized
/// to the planned range total and can be replaced with transport observations
/// using [`ReadPlan::accounting_with_returns`].
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub struct RangePlanCounters {
    pub logical_bytes: u64,
    pub physical_requested_bytes: u64,
    pub physical_returned_bytes: u64,
    pub useful_decoded_bytes: u64,
    pub cache_bytes: u64,
    pub overread_bytes: u64,
    pub physical_requests: u64,
    pub cache_hits: u64,
    pub cache_misses: u64,
}

#[derive(Clone, Debug)]
struct WorkingGroup {
    start: u64,
    end: u64,
    requests: Vec<RangeRequest>,
    intervals: Vec<(u64, u64)>,
}

impl WorkingGroup {
    fn new(request: RangeRequest, end: u64) -> Self {
        Self {
            start: request.range.offset,
            end,
            intervals: vec![(request.range.offset, end)],
            requests: vec![request],
        }
    }

    fn can_add(
        &self,
        request: RangeRequest,
        end: u64,
        limits: RangePlanLimits,
    ) -> Result<bool, ReadPlanError> {
        let gap = request.range.offset.saturating_sub(self.end);
        if gap > limits.max_gap_bytes {
            return Ok(false);
        }
        let merged_end = self.end.max(end);
        let merged_len = merged_end
            .checked_sub(self.start)
            .ok_or(ReadPlanError::AccountingOverflow)?;
        if merged_len > limits.max_physical_range_bytes {
            return Ok(false);
        }
        let mut intervals = self.intervals.clone();
        intervals.push((request.range.offset, end));
        let useful = union_bytes(&mut intervals)?;
        let overread = merged_len
            .checked_sub(useful)
            .ok_or(ReadPlanError::AccountingOverflow)?;
        Ok(overread <= limits.max_overread_bytes)
    }

    fn add(&mut self, request: RangeRequest, end: u64) {
        self.end = self.end.max(end);
        self.intervals.push((request.range.offset, end));
        self.requests.push(request);
    }
}

/// A complete deterministic physical range plan.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReadPlan {
    requests: Vec<RangeRequest>,
    mappings: Vec<LogicalMapping>,
    physical_ranges: Vec<PhysicalRange>,
    batches: Vec<ReadBatch>,
    counters: RangePlanCounters,
    limits: RangePlanLimits,
}

impl ReadPlan {
    #[must_use]
    pub fn requests(&self) -> &[RangeRequest] {
        &self.requests
    }

    #[must_use]
    pub fn mappings(&self) -> &[LogicalMapping] {
        &self.mappings
    }

    #[must_use]
    pub fn physical_ranges(&self) -> &[PhysicalRange] {
        &self.physical_ranges
    }

    #[must_use]
    pub fn batches(&self) -> &[ReadBatch] {
        &self.batches
    }

    #[must_use]
    pub fn limits(&self) -> RangePlanLimits {
        self.limits
    }

    #[must_use]
    pub fn counters(&self) -> RangePlanCounters {
        self.counters
    }

    /// Recompute returned/useful/over-read counters from physical transport
    /// observations.  Every physical index must occur at most once; omitted
    /// indexes retain the plan's requested-byte estimate.
    pub fn accounting_with_returns(
        &self,
        observations: impl IntoIterator<Item = PhysicalReadObservation>,
    ) -> Result<RangePlanCounters, ReadPlanError> {
        let mut counters = self.counters;
        counters.physical_returned_bytes = 0;
        counters.useful_decoded_bytes = counters.cache_bytes;
        let mut seen = BTreeSet::new();
        for observation in observations {
            if observation.physical_index >= self.physical_ranges.len() {
                return Err(ReadPlanError::UnknownPhysicalIndex(
                    observation.physical_index,
                ));
            }
            if !seen.insert(observation.physical_index) {
                return Err(ReadPlanError::DuplicatePhysicalResponse(
                    observation.physical_index,
                ));
            }
            if observation.returned_bytes < observation.useful_decoded_bytes {
                return Err(ReadPlanError::ResponseShorterThanUseful(
                    observation.physical_index,
                ));
            }
            counters.physical_returned_bytes = counters
                .physical_returned_bytes
                .checked_add(observation.returned_bytes)
                .ok_or(ReadPlanError::AccountingOverflow)?;
            counters.useful_decoded_bytes = counters
                .useful_decoded_bytes
                .checked_add(observation.useful_decoded_bytes)
                .ok_or(ReadPlanError::AccountingOverflow)?;
        }
        for physical in &self.physical_ranges {
            if seen.contains(&physical.index) {
                continue;
            }
            counters.physical_returned_bytes = counters
                .physical_returned_bytes
                .checked_add(physical.range.length)
                .ok_or(ReadPlanError::AccountingOverflow)?;
            counters.useful_decoded_bytes = counters
                .useful_decoded_bytes
                .checked_add(physical.useful_logical_bytes)
                .ok_or(ReadPlanError::AccountingOverflow)?;
        }
        counters.overread_bytes = counters.physical_returned_bytes.saturating_sub(
            counters
                .useful_decoded_bytes
                .saturating_sub(counters.cache_bytes),
        );
        Ok(counters)
    }
}

/// Checked range planner.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct RangePlanner {
    limits: RangePlanLimits,
}

impl RangePlanner {
    pub fn new(limits: RangePlanLimits) -> Result<Self, ReadPlanError> {
        Ok(Self {
            limits: limits.validate()?,
        })
    }

    #[must_use]
    pub const fn with_defaults() -> Self {
        Self {
            limits: RangePlanLimits {
                max_gap_bytes: 4096,
                max_overread_bytes: 8192,
                max_physical_range_bytes: 1024 * 1024,
                max_concurrency: 4,
            },
        }
    }

    #[must_use]
    pub const fn limits(self) -> RangePlanLimits {
        self.limits
    }

    /// Sort, coalesce, map, account, and batch logical ranges.
    pub fn plan(
        self,
        requests: impl IntoIterator<Item = RangeRequest>,
    ) -> Result<ReadPlan, ReadPlanError> {
        let requests = requests.into_iter().collect::<Vec<_>>();
        let mut ids = BTreeSet::new();
        let mut sorted = Vec::with_capacity(requests.len());
        for request in &requests {
            if !ids.insert(request.id) {
                return Err(ReadPlanError::DuplicateRequestId(request.id));
            }
            let end = request
                .range
                .offset
                .checked_add(request.range.length)
                .ok_or(ReadPlanError::RangeOverflow {
                    request_id: request.id,
                })?;
            if !request.cache_hit && request.range.length > self.limits.max_physical_range_bytes {
                return Err(ReadPlanError::RangeExceedsPhysicalCap {
                    request_id: request.id,
                });
            }
            sorted.push((*request, end));
        }
        sorted.sort_by_key(|(request, _)| (request.range.offset, request.id));

        let mut physical_ranges = Vec::new();
        let mut mappings_by_id = BTreeMap::new();
        let mut current: Option<WorkingGroup> = None;
        for (request, end) in sorted {
            if request.cache_hit || request.range.length == 0 {
                mappings_by_id.insert(
                    request.id,
                    LogicalMapping {
                        request_id: request.id,
                        logical_range: request.range,
                        physical_index: None,
                        physical_offset: None,
                        cache_hit: request.cache_hit,
                    },
                );
                continue;
            }
            let Some(group) = current.as_mut() else {
                current = Some(WorkingGroup::new(request, end));
                continue;
            };
            if group.can_add(request, end, self.limits)? {
                group.add(request, end);
            } else {
                flush_group(group, &mut physical_ranges, &mut mappings_by_id)?;
                current = Some(WorkingGroup::new(request, end));
            }
        }
        if let Some(group) = current {
            flush_group(&group, &mut physical_ranges, &mut mappings_by_id)?;
        }

        // Retain caller order for mappings and counters while physical ranges
        // remain sorted for deterministic transport behavior.
        let mappings = requests
            .iter()
            .map(|request| {
                mappings_by_id
                    .get(&request.id)
                    .copied()
                    .unwrap_or(LogicalMapping {
                        request_id: request.id,
                        logical_range: request.range,
                        physical_index: None,
                        physical_offset: None,
                        cache_hit: request.cache_hit,
                    })
            })
            .collect::<Vec<_>>();

        let batches = physical_ranges
            .chunks(self.limits.max_concurrency)
            .enumerate()
            .map(|(ordinal, chunk)| ReadBatch {
                ordinal,
                physical_indices: chunk.iter().map(|physical| physical.index).collect(),
            })
            .collect::<Vec<_>>();
        let counters = plan_counters(&requests, &physical_ranges)?;
        Ok(ReadPlan {
            requests,
            mappings,
            physical_ranges,
            batches,
            counters,
            limits: self.limits,
        })
    }
}

impl Default for RangePlanner {
    fn default() -> Self {
        Self::with_defaults()
    }
}

fn flush_group(
    group: &WorkingGroup,
    physical_ranges: &mut Vec<PhysicalRange>,
    mappings_by_id: &mut BTreeMap<u64, LogicalMapping>,
) -> Result<(), ReadPlanError> {
    let index = physical_ranges.len();
    let mut intervals = group.intervals.clone();
    let useful = union_bytes(&mut intervals)?;
    let physical_len = group
        .end
        .checked_sub(group.start)
        .ok_or(ReadPlanError::AccountingOverflow)?;
    let overread = physical_len
        .checked_sub(useful)
        .ok_or(ReadPlanError::AccountingOverflow)?;
    let logical_request_ids = group
        .requests
        .iter()
        .map(|request| request.id)
        .collect::<Vec<_>>();
    physical_ranges.push(PhysicalRange {
        index,
        range: ByteRange {
            offset: group.start,
            length: physical_len,
        },
        logical_request_ids,
        useful_logical_bytes: useful,
        overread_bytes: overread,
    });
    for request in &group.requests {
        mappings_by_id.insert(
            request.id,
            LogicalMapping {
                request_id: request.id,
                logical_range: request.range,
                physical_index: Some(index),
                physical_offset: Some(request.range.offset - group.start),
                cache_hit: false,
            },
        );
    }
    Ok(())
}

fn union_bytes(intervals: &mut [(u64, u64)]) -> Result<u64, ReadPlanError> {
    intervals.sort_unstable_by_key(|(start, end)| (*start, *end));
    let Some((mut current_start, mut current_end)) = intervals.first().copied() else {
        return Ok(0);
    };
    let mut total = 0u64;
    for (start, end) in intervals.iter().copied().skip(1) {
        if start <= current_end {
            current_end = current_end.max(end);
        } else {
            total = total
                .checked_add(
                    current_end
                        .checked_sub(current_start)
                        .ok_or(ReadPlanError::AccountingOverflow)?,
                )
                .ok_or(ReadPlanError::AccountingOverflow)?;
            current_start = start;
            current_end = end;
        }
    }
    total
        .checked_add(
            current_end
                .checked_sub(current_start)
                .ok_or(ReadPlanError::AccountingOverflow)?,
        )
        .ok_or(ReadPlanError::AccountingOverflow)
}

fn plan_counters(
    requests: &[RangeRequest],
    physical_ranges: &[PhysicalRange],
) -> Result<RangePlanCounters, ReadPlanError> {
    let mut counters = RangePlanCounters::default();
    for request in requests {
        counters.logical_bytes = counters
            .logical_bytes
            .checked_add(request.range.length)
            .ok_or(ReadPlanError::AccountingOverflow)?;
        if request.cache_hit {
            counters.cache_hits = counters.cache_hits.saturating_add(1);
            counters.cache_bytes = counters
                .cache_bytes
                .checked_add(request.range.length)
                .ok_or(ReadPlanError::AccountingOverflow)?;
        } else {
            counters.cache_misses = counters.cache_misses.saturating_add(1);
        }
    }
    for physical in physical_ranges {
        counters.physical_requested_bytes = counters
            .physical_requested_bytes
            .checked_add(physical.range.length)
            .ok_or(ReadPlanError::AccountingOverflow)?;
        counters.physical_returned_bytes = counters
            .physical_returned_bytes
            .checked_add(physical.range.length)
            .ok_or(ReadPlanError::AccountingOverflow)?;
        counters.useful_decoded_bytes = counters
            .useful_decoded_bytes
            .checked_add(physical.useful_logical_bytes)
            .ok_or(ReadPlanError::AccountingOverflow)?;
        counters.overread_bytes = counters
            .overread_bytes
            .checked_add(physical.overread_bytes)
            .ok_or(ReadPlanError::AccountingOverflow)?;
    }
    counters.useful_decoded_bytes = counters
        .useful_decoded_bytes
        .checked_add(counters.cache_bytes)
        .ok_or(ReadPlanError::AccountingOverflow)?;
    counters.physical_requests =
        u64::try_from(physical_ranges.len()).map_err(|_| ReadPlanError::AccountingOverflow)?;
    Ok(counters)
}
