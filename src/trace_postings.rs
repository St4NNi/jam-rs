use crate::jidx_reader::SeedOccurrence;
use crate::shared_reader::{SharedGroup, SharedMember, SharedReader};
use crate::trace::TraceError;
use crate::trace_batch::{BatchPosting, worker_cpu_ns};
use crate::trace_index::{TraceDocument as SeedDocument, TraceSeed};
use rayon::prelude::*;
use std::mem::size_of;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

const POSTING_MEMBER_ROWS: usize = 1024;
const POSTING_POSITION_ROWS: usize = 4096;
const POSTING_FRAGMENTS: usize = 128;
const POSTING_LANES: usize = 16;
const POSTING_PARALLEL_ROWS: u64 = 4096;

#[derive(Default)]
pub(crate) struct PostingExecution {
    pub(crate) retained_bytes: usize,
    pub(crate) peak_bytes: usize,
    pub(crate) complete: bool,
    pub(crate) membership_ns: u64,
    pub(crate) position_ns: u64,
    pub(crate) histogram: [u64; 16],
    pub(crate) member_tasks: usize,
    pub(crate) position_tasks: usize,
    pub(crate) task_hash: u64,
    pub(crate) admitted_member_rows: u64,
    pub(crate) admitted_position_rows: u64,
    pub(crate) peak_parallel_tasks: usize,
    pub(crate) scratch_bytes: usize,
    pub(crate) member_copies: u64,
    pub(crate) initialized_position_bytes: usize,
    pub(crate) member_worker_elapsed_ns: u64,
    pub(crate) position_worker_elapsed_ns: u64,
    pub(crate) member_worker_cpu_ns: u64,
    pub(crate) position_worker_cpu_ns: u64,
    pub(crate) cpu_unavailable: usize,
}

#[derive(Clone, Copy, Default)]
struct PostingPlan {
    documents: bool,
    positions: bool,
}

#[derive(Clone, Copy)]
struct MemberFragment {
    ordinal: usize,
    start: u32,
    count: usize,
}

#[derive(Default)]
struct PostingTiming {
    elapsed_ns: u64,
    cpu_ns: Option<u64>,
}

struct MemberLane {
    fragments: Vec<MemberFragment>,
    members: Vec<SharedMember>,
    failed_fragment: usize,
    error: Option<TraceError>,
    timing: PostingTiming,
}

struct PositionFragment<'a> {
    ordinal: usize,
    member_ordinal: usize,
    group: SharedGroup,
    member: SharedMember,
    start: u64,
    output: &'a mut [SeedOccurrence],
}

struct PositionLane<'a> {
    fragments: Vec<PositionFragment<'a>>,
    error: Option<TraceError>,
    timing: PostingTiming,
}

// This pure seam takes bytes remaining AFTER fixed work-table admission.
// It preserves member-then-position admission before the next group.
fn posting_admission(remaining: usize, members: usize, occurrences: u64) -> Option<(bool, usize)> {
    let documents = members.checked_mul(size_of::<SeedDocument>())?;
    let after_documents = remaining.checked_sub(documents)?;
    let positions = usize::try_from(occurrences)
        .ok()
        .and_then(|count| count.checked_mul(size_of::<SeedOccurrence>()))
        .and_then(|bytes| {
            bytes.checked_add(members.checked_mul(size_of::<Vec<SeedOccurrence>>())?)
        });
    match positions.filter(|&bytes| bytes <= after_documents) {
        Some(bytes) => Some((true, documents + bytes)),
        None => Some((false, documents)),
    }
}

fn posting_vec<T>(count: usize) -> Result<Vec<T>, TraceError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| crate::shared_format::SharedError::ResourceLimit)?;
    // No following task is permitted to grow a buffer. Reject an allocator's
    // larger reported capacity before any dispatch or publication.
    if result.capacity() != count {
        return Err(crate::shared_format::SharedError::ResourceLimit.into());
    }
    Ok(result)
}

fn posting_group(entries: &[(u64, Option<TraceSeed>)], ordinal: usize) -> SharedGroup {
    match entries[ordinal].1 {
        Some(TraceSeed::Shared(group)) => group,
        _ => unreachable!("posting planner only admits shared groups"),
    }
}

fn posting_shared_error(message: &'static str) -> TraceError {
    crate::shared_format::SharedError::Invalid(message).into()
}

fn posting_hash(hash: &mut u64, words: impl IntoIterator<Item = u64>) {
    for word in words {
        *hash = (*hash ^ word).wrapping_mul(0x100_0000_01b3);
    }
}

fn posting_timing(started: Option<Instant>, cpu: Option<u64>) -> PostingTiming {
    PostingTiming {
        elapsed_ns: started.map_or(0, |start| start.elapsed().as_nanos() as u64),
        cpu_ns: cpu.and_then(|start| worker_cpu_ns()?.checked_sub(start)),
    }
}

fn record_posting_timing(
    execution: &mut PostingExecution,
    timing: &PostingTiming,
    membership: bool,
    observed: bool,
) {
    if !observed {
        return;
    }
    let (elapsed, cpu) = if membership {
        (
            &mut execution.member_worker_elapsed_ns,
            &mut execution.member_worker_cpu_ns,
        )
    } else {
        (
            &mut execution.position_worker_elapsed_ns,
            &mut execution.position_worker_cpu_ns,
        )
    };
    *elapsed += timing.elapsed_ns;
    if let Some(value) = timing.cpu_ns {
        *cpu += value;
    } else {
        execution.cpu_unavailable += 1;
    }
}

fn run_member_lane(
    reader: &SharedReader,
    entries: &[(u64, Option<TraceSeed>)],
    lane: &mut MemberLane,
    observed: bool,
    active: &AtomicUsize,
    peak: &AtomicUsize,
) {
    let started = observed.then(Instant::now);
    let cpu = observed.then(worker_cpu_ns).flatten();
    if observed {
        peak.fetch_max(
            active.fetch_add(1, Ordering::Relaxed) + 1,
            Ordering::Relaxed,
        );
    }
    lane.error = (|| {
        lane.failed_fragment = usize::MAX;
        let operation = reader.posting_operation()?;
        for (fragment_index, fragment) in lane.fragments.iter().enumerate() {
            lane.failed_fragment = fragment_index;
            operation.append_member_range(
                posting_group(entries, fragment.ordinal),
                fragment.start,
                fragment.count,
                &mut lane.members,
            )?;
        }
        lane.failed_fragment = usize::MAX;
        operation.finish()?;
        Ok::<(), TraceError>(())
    })()
    .err();
    if observed {
        active.fetch_sub(1, Ordering::Relaxed);
    }
    lane.timing = posting_timing(started, cpu);
}

fn run_position_lane(
    reader: &SharedReader,
    lane: &mut PositionLane<'_>,
    observed: bool,
    active: &AtomicUsize,
    peak: &AtomicUsize,
) {
    let started = observed.then(Instant::now);
    let cpu = observed.then(worker_cpu_ns).flatten();
    if observed {
        peak.fetch_max(
            active.fetch_add(1, Ordering::Relaxed) + 1,
            Ordering::Relaxed,
        );
    }
    lane.error = (|| {
        let operation = reader.posting_operation()?;
        for fragment in &mut lane.fragments {
            operation.fill_occurrence_block(
                fragment.group,
                fragment.member,
                fragment.start,
                fragment.output,
            )?;
        }
        operation.finish()?;
        Ok::<(), TraceError>(())
    })()
    .err();
    if observed {
        active.fetch_sub(1, Ordering::Relaxed);
    }
    lane.timing = posting_timing(started, cpu);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn posting_admission_preserves_interleaving_and_exact_boundaries() {
        let document = size_of::<SeedDocument>();
        let positions = size_of::<Vec<SeedOccurrence>>() + 2 * size_of::<SeedOccurrence>();
        let exact = document + positions;
        assert_eq!(posting_admission(exact, 1, 2), Some((true, exact)));
        assert_eq!(posting_admission(exact - 1, 1, 2), Some((false, document)));
        assert_eq!(posting_admission(document - 1, 1, 2), None);
        let mut remaining = exact + 2 * document - 1;
        remaining -= posting_admission(remaining, 1, 2).unwrap().1;
        assert_eq!(posting_admission(remaining, 2, 2), None);
        let later = posting_admission(remaining, 1, 1).unwrap();
        assert!(later.0);
        assert!(later.1 <= remaining);
        assert_eq!(posting_admission(usize::MAX, usize::MAX, 1), None);
        assert_eq!(
            posting_admission(exact, 1, u64::MAX),
            Some((false, document))
        );
    }
}
