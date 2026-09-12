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

struct PostingRead<'a> {
    reader: &'a SharedReader,
    entries: &'a [(u64, Option<TraceSeed>)],
    parallel: bool,
    observed: bool,
    active: AtomicUsize,
    peak: AtomicUsize,
}

pub(crate) fn prepare_shared_postings(
    reader: &SharedReader,
    entries: &[(u64, Option<TraceSeed>)],
    postings: &mut [Option<BatchPosting>],
    available: usize,
    observed: bool,
) -> Result<Option<PostingExecution>, TraceError> {
    if entries.len() != postings.len() || postings.iter().any(Option::is_some) {
        return Err(TraceError::Invalid("posting slots"));
    }
    let result = prepare_shared_postings_inner(reader, entries, postings, available, observed);
    if result.is_err() {
        for posting in postings {
            *posting = None;
        }
    }
    result
}

// available excludes existing lookup state and ordinal posting slots. The caller
// retains its enclosing reservation until this function's work tables are gone.
fn prepare_shared_postings_inner(
    reader: &SharedReader,
    entries: &[(u64, Option<TraceSeed>)],
    postings: &mut [Option<BatchPosting>],
    available: usize,
    observed: bool,
) -> Result<Option<PostingExecution>, TraceError> {
    let mut all_members = 0u64;
    let mut work = 0u64;
    for &(_, seed) in entries {
        match seed {
            Some(TraceSeed::Shared(group)) => {
                all_members = all_members.saturating_add(u64::from(group.member_count()));
                work = work
                    .saturating_add(u64::from(group.member_count()))
                    .saturating_add(group.occurrence_count());
            }
            Some(_) => return Err(TraceError::Invalid("posting index")),
            None => {}
        }
    }
    let mut execution = PostingExecution {
        complete: true,
        task_hash: 0xcbf2_9ce4_8422_2325,
        ..PostingExecution::default()
    };
    if all_members == 0 {
        return Ok(Some(execution));
    }
    let lane_count = if work < POSTING_PARALLEL_ROWS {
        1
    } else {
        POSTING_LANES
    };
    let member_capacity = all_members.min(POSTING_MEMBER_ROWS as u64) as usize;
    // Both lane tables stay allocated through both passes. This deliberately
    // charges their sum rather than releasing a still-live reservation.
    let workspace = entries
        .len()
        .checked_mul(size_of::<PostingPlan>())
        .and_then(|bytes| {
            bytes.checked_add(lane_count.checked_mul(
                size_of::<MemberLane>()
                    + size_of::<PositionLane<'_>>()
                    + POSTING_FRAGMENTS
                        * (size_of::<MemberFragment>() + size_of::<PositionFragment<'_>>())
                    + member_capacity * size_of::<SharedMember>(),
            )?)
        })
        .and_then(|bytes| bytes.checked_add(4096));
    let Some(workspace) = workspace.filter(|&bytes| bytes <= available) else {
        return Ok(None);
    };
    let mut plans = posting_vec(entries.len())?;
    plans.resize(entries.len(), PostingPlan::default());
    let mut member_lanes = posting_vec(lane_count)?;
    let mut position_lanes = posting_vec(lane_count)?;
    for _ in 0..lane_count {
        member_lanes.push(MemberLane {
            fragments: posting_vec(POSTING_FRAGMENTS)?,
            members: posting_vec(member_capacity)?,
            failed_fragment: 0,
            error: None,
            timing: PostingTiming::default(),
        });
        position_lanes.push(PositionLane {
            fragments: posting_vec(POSTING_FRAGMENTS)?,
            error: None,
            timing: PostingTiming::default(),
        });
    }
    execution.scratch_bytes = workspace;
    for (ordinal, &(_, seed)) in entries.iter().enumerate() {
        let Some(TraceSeed::Shared(group)) = seed else {
            continue;
        };
        let count = group.member_count() as usize;
        let remaining = available - workspace - execution.retained_bytes;
        let Some((positions, bytes)) =
            posting_admission(remaining, count, group.occurrence_count())
        else {
            execution.complete = false;
            continue;
        };
        let documents = posting_vec(count)?;
        let occurrences = positions.then(|| posting_vec(count)).transpose()?;
        postings[ordinal] = Some(BatchPosting {
            documents,
            occurrences,
        });
        plans[ordinal] = PostingPlan {
            documents: true,
            positions,
        };
        execution.retained_bytes += bytes;
        execution.admitted_member_rows += u64::from(group.member_count());
        if positions {
            execution.admitted_position_rows += group.occurrence_count();
        }
    }
    execution.peak_bytes = workspace + execution.retained_bytes;
    let parallel = execution
        .admitted_member_rows
        .saturating_add(execution.admitted_position_rows)
        >= POSTING_PARALLEL_ROWS;
    let read = PostingRead {
        reader,
        entries,
        parallel,
        observed,
        active: AtomicUsize::new(0),
        peak: AtomicUsize::new(0),
    };
    let started = observed.then(Instant::now);
    let membership_error =
        execute_member_waves(&read, postings, &plans, &mut member_lanes, &mut execution)?;
    execution.membership_ns = started.map_or(0, |start| start.elapsed().as_nanos() as u64);
    let before = membership_error
        .as_ref()
        .map_or(entries.len(), |&(ordinal, _)| ordinal);
    let started = observed.then(Instant::now);
    initialize_posting_positions(postings, &plans, before, &mut execution)?;
    let position_result =
        execute_position_waves(&read, postings, before, position_lanes, &mut execution);
    execution.position_ns = started.map_or(0, |start| start.elapsed().as_nanos() as u64);
    execution.peak_parallel_tasks = read.peak.load(Ordering::Relaxed);
    position_result?;
    if let Some((_, error)) = membership_error {
        return Err(error);
    }
    Ok(Some(execution))
}

fn execute_member_waves(
    read: &PostingRead<'_>,
    postings: &mut [Option<BatchPosting>],
    plans: &[PostingPlan],
    lanes: &mut [MemberLane],
    execution: &mut PostingExecution,
) -> Result<Option<(usize, TraceError)>, TraceError> {
    let mut ordinal = 0;
    let mut member_start = 0u32;
    let mut occurrence_sum = 0u64;
    loop {
        let mut used = 0;
        for lane in lanes.iter_mut() {
            lane.fragments.clear();
            lane.members.clear();
            lane.error = None;
            let mut rows = 0;
            while rows < lane.members.capacity() && lane.fragments.len() < POSTING_FRAGMENTS {
                while ordinal < plans.len() && !plans[ordinal].documents {
                    ordinal += 1;
                }
                if ordinal == plans.len() {
                    break;
                }
                let group = posting_group(read.entries, ordinal);
                let count = ((group.member_count() - member_start) as usize)
                    .min(lane.members.capacity() - rows);
                lane.fragments.push(MemberFragment {
                    ordinal,
                    start: member_start,
                    count,
                });
                rows += count;
                member_start += count as u32;
                if member_start == group.member_count() {
                    ordinal += 1;
                    member_start = 0;
                }
            }
            if lane.fragments.is_empty() {
                break;
            }
            used += 1;
            execution.member_tasks += 1;
            if read.observed {
                posting_hash(&mut execution.task_hash, [0, lane.fragments.len() as u64]);
                for fragment in &lane.fragments {
                    posting_hash(
                        &mut execution.task_hash,
                        [
                            fragment.ordinal as u64,
                            u64::from(fragment.start),
                            fragment.count as u64,
                        ],
                    );
                }
            }
        }
        if used == 0 {
            return Ok(None);
        }
        let run = |lane: &mut MemberLane| {
            run_member_lane(
                read.reader,
                read.entries,
                lane,
                read.observed,
                &read.active,
                &read.peak,
            )
        };
        if read.parallel {
            lanes[..used].par_iter_mut().for_each(run);
        } else {
            lanes[..used].iter_mut().for_each(run);
        }
        for lane in &lanes[..used] {
            record_posting_timing(execution, &lane.timing, true, read.observed);
        }
        for lane in &mut lanes[..used] {
            if lane.failed_fragment == usize::MAX
                && let Some(error) = lane.error.take()
            {
                return Err(error);
            }
        }
        for lane in &mut lanes[..used] {
            let mut decoded = lane.members.iter().copied();
            for (fragment_index, fragment) in lane.fragments.iter().enumerate() {
                let group = posting_group(read.entries, fragment.ordinal);
                let posting = postings[fragment.ordinal].as_mut().unwrap();
                if posting.documents.len() != fragment.start as usize {
                    return Ok(Some((
                        fragment.ordinal,
                        TraceError::Invalid("posting member progress"),
                    )));
                }
                if fragment.start == 0 {
                    occurrence_sum = 0;
                }
                for member in decoded.by_ref().take(fragment.count) {
                    if posting
                        .documents
                        .last()
                        .is_some_and(|previous| previous.metagenome_id() >= member.metagenome_id)
                    {
                        return Ok(Some((
                            fragment.ordinal,
                            posting_shared_error("member order"),
                        )));
                    }
                    let Some(sum) = occurrence_sum.checked_add(member.occurrence_count()) else {
                        return Ok(Some((
                            fragment.ordinal,
                            posting_shared_error("member occurrence count"),
                        )));
                    };
                    occurrence_sum = sum;
                    posting.documents.push(SeedDocument::Shared(member));
                    execution.member_copies += 1;
                }
                if posting.documents.len() == group.member_count() as usize {
                    if occurrence_sum != group.occurrence_count() {
                        return Ok(Some((
                            fragment.ordinal,
                            posting_shared_error("member occurrence count"),
                        )));
                    }
                    if read.observed {
                        execution.histogram[occurrence_sum.ilog2().min(15) as usize] += 1;
                    }
                }
                if fragment_index == lane.failed_fragment
                    && let Some(error) = lane.error.take()
                {
                    return Ok(Some((fragment.ordinal, error)));
                }
                if posting.documents.len() != fragment.start as usize + fragment.count {
                    return Ok(Some((
                        fragment.ordinal,
                        TraceError::Invalid("posting member progress"),
                    )));
                }
            }
        }
    }
}

fn initialize_posting_positions(
    postings: &mut [Option<BatchPosting>],
    plans: &[PostingPlan],
    before: usize,
    execution: &mut PostingExecution,
) -> Result<(), TraceError> {
    for (posting, plan) in postings[..before].iter_mut().zip(plans) {
        if !plan.positions {
            continue;
        }
        let posting = posting.as_mut().unwrap();
        let occurrences = posting.occurrences.as_mut().unwrap();
        for document in &posting.documents {
            let count = usize::try_from(document.occurrence_count())
                .map_err(|_| TraceError::Invalid("posting position count"))?;
            let mut values = posting_vec(count)?;
            values.resize(
                count,
                SeedOccurrence {
                    contig_id: 0,
                    position: 0,
                    canonical_orientation: false,
                },
            );
            execution.initialized_position_bytes += values.capacity() * size_of::<SeedOccurrence>();
            occurrences.push(values);
        }
    }
    Ok(())
}

fn execute_position_waves<'a>(
    read: &PostingRead<'_>,
    postings: &'a mut [Option<BatchPosting>],
    before: usize,
    mut lanes: Vec<PositionLane<'a>>,
    execution: &mut PostingExecution,
) -> Result<(), TraceError> {
    let mut blocks = postings[..before]
        .iter_mut()
        .enumerate()
        .flat_map(|(ordinal, posting)| {
            posting.iter_mut().flat_map(move |posting| {
                let group = posting_group(read.entries, ordinal);
                let documents = &posting.documents;
                posting.occurrences.iter_mut().flat_map(move |members| {
                    members.iter_mut().zip(documents).enumerate().flat_map(
                        move |(member_ordinal, (values, document))| {
                            let SeedDocument::Shared(member) = *document else {
                                unreachable!()
                            };
                            values.chunks_mut(POSTING_POSITION_ROWS).enumerate().map(
                                move |(block, output)| PositionFragment {
                                    ordinal,
                                    member_ordinal,
                                    group,
                                    member,
                                    start: (block * POSTING_POSITION_ROWS) as u64,
                                    output,
                                },
                            )
                        },
                    )
                })
            })
        });
    let mut pending = None;
    loop {
        for lane in &mut lanes {
            lane.fragments.clear();
            lane.error = None;
        }
        let mut used = 0;
        for lane in lanes.iter_mut() {
            let mut rows = 0;
            while rows < POSTING_POSITION_ROWS && lane.fragments.len() < POSTING_FRAGMENTS {
                let Some(fragment) = pending.take().or_else(|| blocks.next()) else {
                    break;
                };
                if rows + fragment.output.len() > POSTING_POSITION_ROWS {
                    pending = Some(fragment);
                    break;
                }
                rows += fragment.output.len();
                lane.fragments.push(fragment);
            }
            if lane.fragments.is_empty() {
                break;
            }
            used += 1;
            execution.position_tasks += 1;
            if read.observed {
                posting_hash(&mut execution.task_hash, [1, lane.fragments.len() as u64]);
                for fragment in &lane.fragments {
                    posting_hash(
                        &mut execution.task_hash,
                        [
                            fragment.ordinal as u64,
                            fragment.member_ordinal as u64,
                            fragment.start,
                            fragment.output.len() as u64,
                        ],
                    );
                }
            }
        }
        if used == 0 {
            return Ok(());
        }
        let run = |lane: &mut PositionLane<'_>| {
            run_position_lane(read.reader, lane, read.observed, &read.active, &read.peak)
        };
        if read.parallel {
            lanes[..used].par_iter_mut().for_each(run);
        } else {
            lanes[..used].iter_mut().for_each(run);
        }
        for lane in &lanes[..used] {
            record_posting_timing(execution, &lane.timing, false, read.observed);
        }
        for lane in &mut lanes[..used] {
            if let Some(error) = lane.error.take() {
                return Err(error);
            }
        }
    }
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
