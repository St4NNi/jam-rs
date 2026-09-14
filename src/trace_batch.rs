use crate::bgzf_cache::BgzfBlockCache;
use crate::jidx_reader::SEED_LOOKUP_BATCH_KEYS;
use crate::trace::{CacheReservation, LOOKUP_CACHE_AVAILABLE, TraceError};
use crate::trace_index::{
    TraceCacheIdentity, TraceDocument as SeedDocument, TraceIndex, TraceSeed,
};
use rayon::prelude::*;
use std::ops::Range;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

pub(crate) fn phase_cost_enabled() -> bool {
    static ENABLED: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ENABLED.get_or_init(|| std::env::var_os("JAM_PHASE_COST").is_some_and(|v| v == "1"))
}

pub(crate) fn phase_stamp() -> Option<(Instant, u64)> {
    if !phase_cost_enabled() {
        return None;
    }
    #[cfg(target_os = "linux")]
    {
        let mut time = libc::timespec {
            tv_sec: 0,
            tv_nsec: 0,
        };
        // SAFETY: time is initialized and writable for this process CPU sample.
        if unsafe { libc::clock_gettime(libc::CLOCK_PROCESS_CPUTIME_ID, &mut time) } == 0 {
            return Some((
                Instant::now(),
                time.tv_sec as u64 * 1_000_000_000 + time.tv_nsec as u64,
            ));
        }
    }
    None
}

pub(crate) fn phase_elapsed(start: Option<(Instant, u64)>) -> [u64; 2] {
    match (start, phase_stamp()) {
        (Some((wall, cpu)), Some((end_wall, end_cpu))) => [
            end_wall.duration_since(wall).as_nanos() as u64,
            end_cpu.saturating_sub(cpu),
        ],
        _ => [0; 2],
    }
}

#[cfg(target_os = "linux")]
pub(crate) fn worker_cpu_ns() -> Option<u64> {
    let mut time = libc::timespec {
        tv_sec: 0,
        tv_nsec: 0,
    };
    // SAFETY: the output pointer refers to an initialized timespec for this call.
    if unsafe { libc::clock_gettime(libc::CLOCK_THREAD_CPUTIME_ID, &mut time) } != 0 {
        return None;
    }
    u64::try_from(time.tv_sec)
        .ok()?
        .checked_mul(1_000_000_000)?
        .checked_add(u64::try_from(time.tv_nsec).ok()?)
}

#[cfg(not(target_os = "linux"))]
pub(crate) fn worker_cpu_ns() -> Option<u64> {
    None
}

pub(crate) struct SharedSeedLookups {
    pub(crate) identity: TraceCacheIdentity,
    pub(crate) entries: Vec<(u64, Option<TraceSeed>)>,
    pub(crate) query_entries: Vec<u64>,
    pub(crate) attempted_keys: usize,
    pub(crate) query_ranges: Vec<Range<usize>>,
    pub(crate) postings: Vec<Option<BatchPosting>>,
    pub(crate) postings_complete: bool,
    pub(crate) capacity_bytes: usize,
    pub(crate) peak_capacity_bound: usize,
    pub(crate) lookup_ns: u64,
    pub(crate) phase_postings_ns: [u64; 2],
    pub(crate) membership_ns: u64,
    pub(crate) position_ns: u64,
    pub(crate) posting_execution: crate::trace_postings::PostingExecution,
    usage: Option<PostingUsage>,
    pub(crate) distinct_cores: u64,
    pub(crate) split_core_resolutions: u64,
    pub(crate) context_reuse_histogram_log2: [u64; 16],
    pub(crate) context_occurrence_histogram_log2: [u64; 16],
    pub(crate) lookup_tasks: usize,
    pub(crate) lookup_plan_hash: u64,
    pub(crate) lookup_dispatch_ns: u64,
    pub(crate) lookup_parallel_ns: u64,
    pub(crate) lookup_compute_ns: u64,
    pub(crate) lookup_reduce_ns: u64,
    pub(crate) lookup_dispatch_to_start_ns: u64,
    pub(crate) _reservation: CacheReservation<'static>,
}

pub(crate) struct BatchPosting {
    pub(crate) documents: Vec<SeedDocument>,
    pub(crate) occurrences: Option<Vec<Vec<crate::jidx_reader::SeedOccurrence>>>,
}

struct PostingUsage {
    offsets: Vec<usize>,
    used: Vec<AtomicU64>,
}

impl PostingUsage {
    fn new(postings: &[Option<BatchPosting>], available: usize) -> Option<Self> {
        let members = postings.iter().flatten().try_fold(0usize, |sum, posting| {
            sum.checked_add(posting.occurrences.as_ref().map_or(0, Vec::len))
        })?;
        let offsets = postings.len().checked_add(1)?;
        let words = members.div_ceil(64);
        let bytes = offsets
            .checked_mul(std::mem::size_of::<usize>())?
            .checked_add(words.checked_mul(std::mem::size_of::<AtomicU64>())?)?;
        if bytes > available {
            return None;
        }
        let mut result = Self {
            offsets: Vec::new(),
            used: Vec::new(),
        };
        result.offsets.try_reserve_exact(offsets).ok()?;
        result.used.try_reserve_exact(words).ok()?;
        if result.capacity_bytes() > available {
            return None;
        }
        result.used.resize_with(words, || AtomicU64::new(0));
        let mut offset = 0;
        for posting in postings {
            result.offsets.push(offset);
            offset += posting
                .as_ref()
                .and_then(|posting| posting.occurrences.as_ref())
                .map_or(0, Vec::len);
        }
        result.offsets.push(offset);
        Some(result)
    }

    fn capacity_bytes(&self) -> usize {
        self.offsets.capacity() * std::mem::size_of::<usize>()
            + self.used.capacity() * std::mem::size_of::<AtomicU64>()
    }
}

impl SharedSeedLookups {
    pub(crate) fn mark_positions_used(&self, ordinal: usize, member: usize) {
        if let Some(usage) = &self.usage {
            let start = usage.offsets[ordinal];
            let offset = start + member;
            assert!(offset < usage.offsets[ordinal + 1]);
            usage.used[offset / 64].fetch_or(1 << (offset % 64), Ordering::Relaxed);
        }
    }

    pub(crate) fn unused_positions(&self) -> Option<(u64, u64)> {
        let usage = self.usage.as_ref()?;
        let mut lists = 0;
        let mut rows = 0;
        for (ordinal, posting) in self.postings.iter().enumerate() {
            let Some(positions) = posting
                .as_ref()
                .and_then(|posting| posting.occurrences.as_ref())
            else {
                continue;
            };
            for (member, positions) in positions.iter().enumerate() {
                let offset = usage.offsets[ordinal] + member;
                if usage.used[offset / 64].load(Ordering::Relaxed) & (1 << (offset % 64)) == 0 {
                    lists += 1;
                    rows += positions.len() as u64;
                }
            }
        }
        Some((lists, rows))
    }

    pub(crate) fn posting(
        &self,
        key: u64,
        ordinal: Option<usize>,
    ) -> Result<Option<&BatchPosting>, TraceError> {
        let ordinal = if let Some(ordinal) = ordinal {
            if !self
                .entries
                .get(ordinal)
                .is_some_and(|&(found, seed)| found == key && seed.is_some())
            {
                return Err(TraceError::Invalid("batch posting ordinal"));
            }
            ordinal
        } else if let Ok(ordinal) = self.entries.binary_search_by_key(&key, |entry| entry.0) {
            ordinal
        } else {
            return Ok(None);
        };
        self.postings
            .get(ordinal)
            .map(Option::as_ref)
            .ok_or(TraceError::Invalid("batch posting ordinal"))
    }
}

pub(crate) struct TraceBatch {
    pub(crate) lookups: Option<SharedSeedLookups>,
    pub(crate) sequence: Arc<BgzfBlockCache>,
}

pub(crate) struct SharedCoreLookups {
    pub(crate) identity: TraceCacheIdentity,
    pub(crate) groups: Vec<crate::shared_reader::SharedGroup>,
    pub(crate) peak_capacity_bound: usize,
    pub(crate) tasks: usize,
    pub(crate) requests: crate::shared_reader::CoreRequestCounts,
    pub(crate) lookup_ns: u64,
    _reservation: CacheReservation<'static>,
}

impl SharedCoreLookups {
    pub(crate) fn capacity_bytes(&self) -> usize {
        self._reservation.bytes
    }
}

pub(crate) fn prepare_cores(
    index: &TraceIndex,
    requests: impl IntoIterator<Item = u32>,
    request_count: usize,
    observed: bool,
) -> Result<Option<Arc<SharedCoreLookups>>, TraceError> {
    prepare_cores_inner(index, requests, request_count, observed, true, None)
}

#[cfg(any(test, feature = "bench-internals"))]
pub(crate) fn prepare_cores_late(
    index: &TraceIndex,
    requests: impl IntoIterator<Item = u32>,
    request_count: usize,
    observed: bool,
) -> Result<Option<Arc<SharedCoreLookups>>, TraceError> {
    prepare_cores_inner(index, requests, request_count, observed, false, None)
}

pub(crate) fn prepare_screened_cores(
    index: &TraceIndex,
    requests: impl IntoIterator<Item = u32>,
    request_count: usize,
    observed: bool,
    operation: &crate::shared_reader::SharedCoreOperation<'_>,
) -> Result<Option<Arc<SharedCoreLookups>>, TraceError> {
    let TraceIndex::Shared(reader) = index else {
        return Err(TraceError::Invalid("screened core index"));
    };
    if !operation.belongs_to(reader) {
        return Err(TraceError::Invalid("screened core reader"));
    }
    prepare_cores_inner(
        index,
        requests,
        request_count,
        observed,
        true,
        Some(operation),
    )
}

fn prepare_cores_inner(
    index: &TraceIndex,
    requests: impl IntoIterator<Item = u32>,
    request_count: usize,
    observed: bool,
    early: bool,
    screened: Option<&crate::shared_reader::SharedCoreOperation<'_>>,
) -> Result<Option<Arc<SharedCoreLookups>>, TraceError> {
    let TraceIndex::Shared(reader) = index else {
        return Ok(None);
    };
    let started = observed.then(Instant::now);
    let Some(identity) = index.cache_file_identity()? else {
        return Ok(None);
    };
    let Some(key_bytes) = request_count.checked_mul(std::mem::size_of::<u32>()) else {
        return Ok(None);
    };
    let Some(key_reservation) = CacheReservation::acquire(&LOOKUP_CACHE_AVAILABLE, key_bytes)
    else {
        return Ok(None);
    };
    let mut keys = Vec::new();
    if key_bytes > lookup_budget(index)
        || keys.try_reserve_exact(request_count).is_err()
        || keys.capacity() > request_count
    {
        return Ok(None);
    }
    let owned_operation = if early && screened.is_none() {
        reader.core_operation()?
    } else {
        None
    };
    let operation = screened.or(owned_operation.as_ref());
    let mut counts = crate::shared_reader::CoreRequestCounts::default();
    for key in requests {
        if counts.attempted == request_count {
            return Err(TraceError::Invalid("core request count"));
        }
        counts.attempted += 1;
        let (covered, keep) = match &operation {
            Some(operation) if screened.is_none() => operation.screen(key)?,
            Some(_) => (false, true),
            None => (false, true),
        };
        counts.covered += usize::from(covered);
        counts.prescreened += usize::from(screened.is_some());
        counts.uncovered += usize::from(!covered && screened.is_none());
        counts.rejected += usize::from(!keep);
        if keep {
            keys.push(key);
        }
    }
    if counts.attempted != request_count {
        return Err(TraceError::Invalid("core request count"));
    }
    counts.retained = keys.len();
    keys.par_sort_unstable();
    keys.dedup();
    counts.planned = keys.len();
    if let Some(operation) = &operation {
        operation.record(counts);
    }
    let tasks = if reader.has_core_prefixes() {
        core_prefix_ranges(&keys).count()
    } else {
        keys.len().div_ceil(SEED_LOOKUP_BATCH_KEYS)
    };
    let group_bytes = std::mem::size_of::<crate::shared_reader::SharedGroup>();
    let wave_tasks = CORE_LOOKUP_TASKS;
    let overhead = 4096 + std::mem::size_of::<SharedCoreLookups>();
    let Some(base_bytes) = tasks
        .checked_mul(
            std::mem::size_of::<Range<usize>>()
                + std::mem::size_of::<Vec<crate::shared_reader::SharedGroup>>(),
        )
        .and_then(|bytes| bytes.checked_add(key_bytes))
        .and_then(|bytes| bytes.checked_add(overhead))
        .and_then(|bytes| {
            bytes.checked_add(wave_tasks.checked_mul(std::mem::size_of::<
                Result<Vec<crate::shared_reader::SharedGroup>, TraceError>,
            >())?)
        })
        .filter(|&bytes| bytes <= lookup_budget(index))
    else {
        return Ok(None);
    };
    let Some(mut reservation) =
        CacheReservation::acquire(&LOOKUP_CACHE_AVAILABLE, base_bytes - key_bytes)
    else {
        return Ok(None);
    };
    let mut ranges = Vec::new();
    let mut chunks = Vec::new();
    let mut wave = Vec::new();
    if ranges.try_reserve_exact(tasks).is_err()
        || chunks.try_reserve_exact(tasks).is_err()
        || wave.try_reserve_exact(wave_tasks).is_err()
        || ranges.capacity() > tasks
        || chunks.capacity() > tasks
        || wave.capacity() > wave_tasks
    {
        return Ok(None);
    }
    if reader.has_core_prefixes() {
        ranges.extend(core_prefix_ranges(&keys));
    } else {
        ranges.extend(
            (0..keys.len())
                .step_by(SEED_LOOKUP_BATCH_KEYS)
                .map(|start| start..(start + SEED_LOOKUP_BATCH_KEYS).min(keys.len())),
        );
    }
    let lookup = |chunk: &[u32]| {
        let mut groups = Vec::new();
        if reader.has_core_prefixes() {
            match &operation {
                Some(operation) => operation.resolve_sorted_cores_into(chunk, &mut groups)?,
                None => reader.resolve_sorted_cores_into(chunk, &mut groups)?,
            }
            return Ok::<_, TraceError>(groups);
        }
        let mut contexts = Vec::new();
        contexts
            .try_reserve_exact(chunk.len())
            .map_err(|_| TraceError::Shared(crate::shared_format::SharedError::ResourceLimit))?;
        if contexts.capacity() > chunk.len() {
            return Err(TraceError::Shared(
                crate::shared_format::SharedError::ResourceLimit,
            ));
        }
        contexts.extend(
            chunk
                .iter()
                .copied()
                .map(crate::shared_seed::SharedKey::core),
        );
        let found = reader.find_many(&contexts)?;
        groups
            .try_reserve_exact(found.iter().flatten().count())
            .map_err(|_| TraceError::Shared(crate::shared_format::SharedError::ResourceLimit))?;
        if groups.capacity() > chunk.len() {
            return Err(TraceError::Shared(
                crate::shared_format::SharedError::ResourceLimit,
            ));
        }
        groups.extend(found.into_iter().flatten());
        Ok(groups)
    };
    let resource_limit = |error: &TraceError| {
        matches!(
            error,
            TraceError::Shared(crate::shared_format::SharedError::ResourceLimit)
        )
    };
    let mut bytes = base_bytes;
    let mut retained_bytes = base_bytes;
    let mut count = 0usize;
    let mut start = 0;
    while start < ranges.len() {
        let row_bytes = 2 * group_bytes
            + if operation.is_none() {
                reader.core_filter_workspace_per_key()
            } else {
                0
            }
            + if reader.has_core_prefixes() {
                0
            } else {
                std::mem::size_of::<crate::shared_seed::SharedKey>()
                    + std::mem::size_of::<Option<crate::shared_reader::SharedGroup>>()
            };
        let mut end = start;
        let mut active_bytes = 0;
        while end < ranges.len() && end - start < wave_tasks {
            let next = active_bytes + ranges[end].len() * row_bytes;
            if retained_bytes + next > lookup_budget(index) {
                break;
            }
            active_bytes = next;
            end += 1;
        }
        if end == start {
            return Ok(None);
        }
        let Some(mut active) = CacheReservation::acquire(&LOOKUP_CACHE_AVAILABLE, active_bytes)
        else {
            return Ok(None);
        };
        reservation.bytes += active.bytes;
        active.bytes = 0;
        bytes = bytes.max(retained_bytes + active_bytes);
        ranges[start..end]
            .par_iter()
            .map(|range| lookup(&keys[range.clone()]))
            .collect_into_vec(&mut wave);
        let mut exhausted = false;
        for chunk in wave.drain(..) {
            match chunk {
                Ok(chunk) => {
                    retained_bytes += chunk.capacity() * group_bytes;
                    count += chunk.len();
                    chunks.push(chunk);
                }
                Err(error) if resource_limit(&error) => exhausted = true,
                Err(error) => return Err(error),
            }
        }
        if exhausted {
            return Ok(None);
        }
        reservation.retain(retained_bytes - key_bytes);
        start = end;
    }
    let merge_bytes = count * group_bytes;
    if retained_bytes + merge_bytes > lookup_budget(index) {
        return Ok(None);
    }
    let Some(mut merge) = CacheReservation::acquire(&LOOKUP_CACHE_AVAILABLE, merge_bytes) else {
        return Ok(None);
    };
    reservation.bytes += merge.bytes;
    merge.bytes = 0;
    bytes = bytes.max(retained_bytes + merge_bytes);
    let mut groups = Vec::new();
    if groups.try_reserve_exact(count).is_err() || groups.capacity() > count {
        return Ok(None);
    }
    for chunk in chunks {
        groups.extend(chunk);
    }
    drop(wave);
    drop(keys);
    drop(key_reservation);
    drop(ranges);
    if let Some(operation) = owned_operation {
        operation.finish()?;
    }
    if index.cache_file_identity()? != Some(identity) {
        return Err(TraceError::Invalid("resolved core identity"));
    }
    reservation.retain(
        4096 + std::mem::size_of::<SharedCoreLookups>()
            + groups.capacity() * std::mem::size_of::<crate::shared_reader::SharedGroup>(),
    );
    Ok(Some(Arc::new(SharedCoreLookups {
        identity,
        groups,
        peak_capacity_bound: bytes,
        tasks,
        requests: counts,
        lookup_ns: started.map_or(0, |started| started.elapsed().as_nanos() as u64),
        _reservation: reservation,
    })))
}

pub(crate) fn lookup_budget(index: &TraceIndex) -> usize {
    if index.is_shared() {
        128 * 1024 * 1024
    } else {
        32 * 1024 * 1024
    }
}

const QUERY_LOOKUP_ROW_BYTES: usize = std::mem::size_of::<(u64, usize)>()
    + std::mem::size_of::<u64>()
    + std::mem::size_of::<usize>()
    + std::mem::size_of::<(u64, Option<TraceSeed>)>();

fn lookup_chunk_keys() -> usize {
    SEED_LOOKUP_BATCH_KEYS
}

/// Target number of exact core lookup tasks, so a sparse survivor set still spans the pool.
const CORE_LOOKUP_TASKS: usize = 64;
/// Smallest exact core lookup task, so small request sets do not pay per-task dispatch.
const MIN_CORE_LOOKUP_TASK_KEYS: usize = 256;

#[cfg(test)]
thread_local! {
    pub(crate) static CONTEXT_ALLOCATION_FAILURE: std::cell::Cell<Option<usize>> = const { std::cell::Cell::new(None) };
}

#[cfg(test)]
pub(crate) fn core_lookup_task_count(keys: &[u32]) -> usize {
    core_prefix_ranges(keys).count()
}

/// Splits sorted cores into lookup tasks whose size depends only on the key count, never on the
/// worker count. A task keeps each core prefix group whole unless the group alone exceeds one
/// batch.
fn core_prefix_ranges(keys: &[u32]) -> impl Iterator<Item = Range<usize>> + '_ {
    let target = keys
        .len()
        .div_ceil(CORE_LOOKUP_TASKS)
        .clamp(MIN_CORE_LOOKUP_TASK_KEYS, SEED_LOOKUP_BATCH_KEYS);
    let mut start = 0;
    std::iter::from_fn(move || {
        if start == keys.len() {
            return None;
        }
        let mut end = (start + target).min(keys.len());
        if end < keys.len() && keys[end - 1] >> 14 == keys[end] >> 14 {
            let prefix = keys[end] >> 14;
            let limit = (start + SEED_LOOKUP_BATCH_KEYS).min(keys.len());
            let group_end = end + keys[end..limit].partition_point(|key| key >> 14 == prefix);
            if group_end == keys.len() || keys[group_end] >> 14 != prefix {
                end = group_end;
            } else {
                let group_start =
                    start + keys[start..end].partition_point(|key| key >> 14 < prefix);
                end = if group_start > start {
                    group_start
                } else {
                    limit
                };
            }
        }
        let range = start..end;
        start = end;
        Some(range)
    })
}

fn lookup_ranges(keys: &[u64], prefixes: bool) -> impl Iterator<Item = (Range<usize>, bool)> + '_ {
    let mut start = 0;
    let core = |key| crate::shared_seed::SharedKey::unpack(key).map(|key| key.core);
    std::iter::from_fn(move || {
        if start == keys.len() {
            return None;
        }
        let mut end = (start + lookup_chunk_keys()).min(keys.len());
        let mut split = false;
        let group = |key| core(key).map(|core| if prefixes { core >> 14 } else { core });
        if end < keys.len() && group(keys[end - 1]) == group(keys[end]) {
            let boundary = end;
            while end > start && group(keys[end - 1]) == group(keys[boundary]) {
                end -= 1;
            }
            if end == start {
                end = boundary;
                if prefixes {
                    while end > start && core(keys[end - 1]) == core(keys[boundary]) {
                        end -= 1;
                    }
                }
                if !prefixes || end == start {
                    end = boundary;
                    split = true;
                }
            }
        }
        let range = start..end;
        start = end;
        Some((range, split))
    })
}

fn lookup_workspace(index: &TraceIndex, requests: usize) -> Option<usize> {
    if index.is_shared() {
        let tasks = requests
            .div_ceil(lookup_chunk_keys())
            .checked_mul(2)?
            .checked_add(1)?;
        let workers = rayon::current_num_threads().min(tasks);
        workers
            .checked_mul(lookup_chunk_keys())?
            .checked_mul(
                std::mem::size_of::<(usize, crate::shared_seed::SharedKey)>()
                    + std::mem::size_of::<crate::shared_seed::SharedKey>()
                    + std::mem::size_of::<Option<crate::shared_reader::SharedGroup>>()
                    + std::mem::size_of::<Option<TraceSeed>>(),
            )?
            .checked_add(tasks.checked_mul(
                std::mem::size_of::<(u64, &[u64], &mut [(u64, Option<TraceSeed>)])>()
                    + std::mem::size_of::<usize>()
                    + std::mem::size_of::<Vec<(usize, crate::shared_reader::SharedMember)>>()
                    + std::mem::size_of::<Result<[u64; 3], TraceError>>(),
            )?)
    } else {
        Some(0)
    }
}

pub(crate) fn lookup_bytes(index: &TraceIndex, requests: usize, queries: usize) -> Option<usize> {
    requests
        .checked_mul(QUERY_LOOKUP_ROW_BYTES)?
        .checked_add(lookup_workspace(index, requests)?)?
        .checked_add(4096)?
        .checked_add(queries.checked_mul(std::mem::size_of::<Range<usize>>())?)
}

#[cfg(test)]
pub(crate) fn prepare_lookup(
    index: &TraceIndex,
    requests: Vec<(u64, usize)>,
    query_count: usize,
    observed: bool,
) -> Result<Option<SharedSeedLookups>, TraceError> {
    prepare_lookup_with_cores(index, requests, query_count, observed, None)
}

pub(crate) fn prepare_lookup_with_cores(
    index: &TraceIndex,
    mut requests: Vec<(u64, usize)>,
    query_count: usize,
    observed: bool,
    cores: Option<&SharedCoreLookups>,
) -> Result<Option<SharedSeedLookups>, TraceError> {
    let started = observed.then(Instant::now);
    let budget = lookup_budget(index);
    let Some(bytes) = lookup_bytes(index, requests.capacity(), query_count) else {
        return Ok(None);
    };
    if bytes > budget {
        return Ok(None);
    }
    let workspace = lookup_workspace(index, requests.len()).unwrap();
    let Some(mut reservation) = CacheReservation::acquire(&LOOKUP_CACHE_AVAILABLE, budget) else {
        return Ok(None);
    };
    let Some(identity) = index.cache_file_identity()? else {
        return Ok(None);
    };
    if index.is_shared() {
        requests.par_sort_unstable_by_key(|&(key, _)| {
            let context = crate::shared_seed::SharedKey::unpack(key).unwrap();
            (context.core, context.context_code().unwrap())
        });
    } else {
        requests.sort_unstable_by_key(|request| request.0);
    }
    let mut context_reuse_histogram_log2 = [0; 16];
    if observed {
        for same_key in requests.chunk_by(|left, right| left.0 == right.0) {
            context_reuse_histogram_log2[same_key.len().ilog2().min(15) as usize] += 1;
        }
    }
    let mut keys = Vec::new();
    keys.try_reserve_exact(requests.len())
        .map_err(|_| TraceError::Invalid("batch query key allocation"))?;
    keys.extend(requests.iter().map(|request| request.0));
    keys.dedup();
    let attempted_keys = keys.len();
    let mut entries = Vec::new();
    if entries.try_reserve_exact(keys.len()).is_err() {
        return Ok(None);
    }
    let mut peak_capacity_bound = 4096
        + std::mem::size_of::<SharedSeedLookups>()
        + requests.capacity() * std::mem::size_of::<(u64, usize)>()
        + keys.capacity() * std::mem::size_of::<u64>()
        + entries.capacity() * std::mem::size_of::<(u64, Option<TraceSeed>)>()
        + keys.len() * std::mem::size_of::<usize>()
        + query_count * std::mem::size_of::<Range<usize>>()
        + workspace;
    if peak_capacity_bound > reservation.bytes {
        return Ok(None);
    }
    let core = |key| crate::shared_seed::SharedKey::unpack(key).map(|key| key.core);
    let distinct_cores = if index.is_shared() {
        keys.iter()
            .enumerate()
            .filter(|&(i, key)| i == 0 || core(*key) != core(keys[i - 1]))
            .count() as u64
    } else {
        0
    };
    let mut split_core_resolutions = 0;
    let mut lookup_tasks = 0;
    let mut lookup_plan_hash = 0xcbf2_9ce4_8422_2325u64;
    let mut lookup_dispatch_ns = 0;
    let mut lookup_parallel_ns = 0;
    let mut lookup_compute_ns = 0;
    let mut lookup_reduce_ns = 0;
    let mut lookup_dispatch_to_start_ns = 0;
    let mut context_members: Vec<Vec<(usize, crate::shared_reader::SharedMember)>> = Vec::new();
    let mut context_member_bytes = 0;
    if let TraceIndex::Shared(reader) = index {
        let dispatch = observed.then(Instant::now);
        let limit = lookup_chunk_keys();
        entries.extend(keys.iter().map(|&key| (key, None)));
        let mut remaining = entries.as_mut_slice();
        let mut tasks = Vec::new();
        tasks
            .try_reserve_exact(keys.len().div_ceil(limit) * 2 + 1)
            .map_err(|_| TraceError::Invalid("lookup task allocation"))?;
        for (range, split) in lookup_ranges(&keys, reader.has_core_prefixes()) {
            let first_key = range.start;
            split_core_resolutions += u64::from(split);
            let chunk = &keys[range];
            if observed {
                for word in [chunk.len() as u64, u64::from(split)]
                    .into_iter()
                    .chain(chunk.iter().copied())
                {
                    lookup_plan_hash = (lookup_plan_hash ^ word).wrapping_mul(0x100_0000_01b3);
                }
            }
            let cores = 1 + chunk
                .windows(2)
                .filter(|pair| core(pair[0]) != core(pair[1]))
                .count();
            let weight = cores as u64 * u64::from(reader.core_count().max(1).ilog2() + 1)
                + chunk.len() as u64;
            let (slots, tail) = remaining.split_at_mut(chunk.len());
            remaining = tail;
            tasks.push((weight, chunk, slots, first_key));
        }
        tasks.sort_unstable_by_key(|task| std::cmp::Reverse(task.0));
        lookup_tasks = tasks.len();
        if context_members.try_reserve_exact(tasks.len()).is_err()
            || context_members.capacity() != tasks.len()
        {
            return Ok(None);
        }
        context_member_bytes = context_members.capacity()
            * std::mem::size_of::<Vec<(usize, crate::shared_reader::SharedMember)>>();
        let context_header_bytes = context_member_bytes;
        let descriptor_bytes = reservation
            .bytes
            .saturating_sub(peak_capacity_bound)
            .min(16 * 1024 * 1024);
        let descriptor_count = if cores.is_some() {
            descriptor_bytes
                / tasks.len().max(1)
                / std::mem::size_of::<(usize, crate::shared_reader::SharedMember)>()
        } else {
            0
        };
        for _task in 0..tasks.len() {
            let mut members = Vec::new();
            #[cfg(test)]
            let descriptor_count = if CONTEXT_ALLOCATION_FAILURE.get() == Some(_task) {
                CONTEXT_ALLOCATION_FAILURE.set(None);
                usize::MAX
            } else {
                descriptor_count
            };
            if members.try_reserve_exact(descriptor_count).is_err()
                || members.capacity() != descriptor_count
            {
                members = Vec::new();
            }
            context_member_bytes += members.capacity()
                * std::mem::size_of::<(usize, crate::shared_reader::SharedMember)>();
            context_members.push(members);
        }
        peak_capacity_bound += context_member_bytes - context_header_bytes;
        if peak_capacity_bound > reservation.bytes {
            return Ok(None);
        }
        lookup_dispatch_ns = dispatch.map_or(0, |start| start.elapsed().as_nanos() as u64);
        let dispatch = observed.then(Instant::now);
        let timings = tasks
            .into_par_iter()
            .zip(context_members.par_iter_mut())
            .map(|((_, keys, slots, first_key), members)| {
                let start = observed.then(Instant::now);
                let queued = start.zip(dispatch).map_or(0, |(start, dispatch)| {
                    start.duration_since(dispatch).as_nanos() as u64
                });
                let seeds = if let Some(cores) = cores {
                    index.find_seeds_in_cores_with_members(keys, cores, members)?
                } else {
                    index.find_seeds_batch(keys)?
                };
                let compute = start.map_or(0, |start| start.elapsed().as_nanos() as u64);
                let start = observed.then(Instant::now);
                for (slot, seed) in slots.iter_mut().zip(seeds) {
                    slot.1 = seed;
                }
                for member in members {
                    member.0 += first_key;
                }
                Ok::<_, TraceError>([
                    queued,
                    compute,
                    start.map_or(0, |start| start.elapsed().as_nanos() as u64),
                ])
            })
            .collect::<Vec<_>>();
        lookup_parallel_ns = dispatch.map_or(0, |start| start.elapsed().as_nanos() as u64);
        for timing in timings {
            let [queued, compute, reduce] = timing?;
            lookup_dispatch_to_start_ns += queued;
            lookup_compute_ns += compute;
            lookup_reduce_ns += reduce;
        }
    } else {
        for chunk in keys.chunks(SEED_LOOKUP_BATCH_KEYS) {
            entries.extend(chunk.iter().copied().zip(index.find_seeds_batch(chunk)?));
        }
    }
    let mut entry = 0;
    requests.retain(|&(key, _)| {
        while entries[entry].0 != key {
            entry += 1;
        }
        entries[entry].1.is_some()
    });
    if context_members.iter().any(|members| !members.is_empty()) {
        let mut ordinals = Vec::new();
        if ordinals.try_reserve_exact(keys.len()).is_err() || ordinals.capacity() != keys.len() {
            if let TraceIndex::Shared(reader) = index {
                reader
                    .record_context_posting_members(0, context_members.iter().map(Vec::len).sum());
            }
            context_members = Vec::new();
            context_member_bytes = 0;
        } else {
            ordinals.extend((0..entries.len()).filter(|&ordinal| entries[ordinal].1.is_some()));
            ordinals.sort_unstable_by_key(|&ordinal| entries[ordinal].0);
            for (new, old) in ordinals.into_iter().enumerate() {
                keys[old] = new as u64;
            }
            for members in &mut context_members {
                for member in members.iter_mut() {
                    member.0 = keys[member.0] as usize;
                }
                members.sort_unstable_by_key(|row| (row.0, row.1.metagenome_id));
            }
        }
    }
    entries.retain(|entry| entry.1.is_some());
    entries.sort_unstable_by_key(|entry| entry.0);
    for request in &mut requests {
        request.0 = entries
            .binary_search_by_key(&request.0, |entry| entry.0)
            .map_err(|_| TraceError::Invalid("successful query association"))?
            as u64;
    }
    requests.sort_unstable_by_key(|&(key, query)| (query, key));
    keys.clear();
    let mut query_ranges = Vec::with_capacity(query_count);
    let mut request = 0;
    for query in 0..query_count {
        let start = request;
        while request < requests.len() && requests[request].1 == query {
            keys.push(requests[request].0);
            request += 1;
        }
        query_ranges.push(start..request);
    }
    drop(requests);
    keys.shrink_to_fit();
    entries.shrink_to_fit();
    let lookup_ns = started.map_or(0, |started| started.elapsed().as_nanos() as u64);
    let posting_phase = phase_stamp();
    let mut membership_ns = 0;
    let mut position_ns = 0;
    let mut capacity_bytes = 4096
        + std::mem::size_of::<SharedSeedLookups>()
        + query_ranges.capacity() * std::mem::size_of::<Range<usize>>()
        + keys.capacity() * std::mem::size_of::<u64>()
        + entries.capacity() * std::mem::size_of::<(u64, Option<TraceSeed>)>();
    let slot_bytes = entries
        .len()
        .checked_mul(std::mem::size_of::<Option<BatchPosting>>());
    if slot_bytes.is_none_or(|bytes| {
        capacity_bytes
            .saturating_add(context_member_bytes)
            .saturating_add(bytes)
            > reservation.bytes
    }) {
        if let TraceIndex::Shared(reader) = index {
            reader.record_context_posting_members(0, context_members.iter().map(Vec::len).sum());
        }
        context_members = Vec::new();
        context_member_bytes = 0;
    }
    if slot_bytes.is_none_or(|bytes| {
        capacity_bytes
            .saturating_add(context_member_bytes)
            .saturating_add(bytes)
            > reservation.bytes
    }) {
        return Ok(None);
    }
    let mut postings = Vec::new();
    postings
        .try_reserve_exact(entries.len())
        .map_err(|_| TraceError::Invalid("posting slots allocation"))?;
    let slot_bytes = postings.capacity() * std::mem::size_of::<Option<BatchPosting>>();
    if capacity_bytes
        .saturating_add(context_member_bytes)
        .saturating_add(slot_bytes)
        > reservation.bytes
    {
        if let TraceIndex::Shared(reader) = index {
            reader.record_context_posting_members(0, context_members.iter().map(Vec::len).sum());
        }
        context_members = Vec::new();
        context_member_bytes = 0;
        if capacity_bytes.saturating_add(slot_bytes) > reservation.bytes {
            return Ok(None);
        }
    }
    capacity_bytes += slot_bytes;
    peak_capacity_bound = peak_capacity_bound.max(capacity_bytes);
    postings.resize_with(entries.len(), || None);
    let mut postings_complete = true;
    let mut context_occurrence_histogram_log2 = [0; 16];
    let execution = if let TraceIndex::Shared(reader) = index {
        match crate::trace_postings::prepare_shared_postings_with_members(
            reader,
            &entries,
            &mut postings,
            reservation
                .bytes
                .saturating_sub(capacity_bytes + context_member_bytes),
            observed,
            &context_members,
        ) {
            Err(TraceError::Shared(crate::shared_format::SharedError::ResourceLimit)) => None,
            result => result?,
        }
    } else {
        None
    };
    let serial_entries = if let Some(execution) = &execution {
        peak_capacity_bound =
            peak_capacity_bound.max(capacity_bytes + context_member_bytes + execution.peak_bytes);
        capacity_bytes += execution.retained_bytes;
        membership_ns = execution.membership_ns;
        position_ns = execution.position_ns;
        postings_complete = execution.complete;
        context_occurrence_histogram_log2 = execution.histogram;
        &[][..]
    } else {
        if let TraceIndex::Shared(reader) = index {
            reader.record_context_posting_members(0, context_members.iter().map(Vec::len).sum());
        }
        entries.as_slice()
    };
    drop(context_members);
    for (ordinal, &(_, seed)) in serial_entries.iter().enumerate() {
        let Some(seed) = seed else {
            continue;
        };
        let document_bytes = (seed.document_frequency() as usize)
            .saturating_mul(std::mem::size_of::<SeedDocument>());
        let document_workspace = (seed.document_frequency() as usize).saturating_mul(match index {
            TraceIndex::Shared(_) => std::mem::size_of::<crate::shared_reader::SharedMember>(),
            TraceIndex::Shard(_) => std::mem::size_of::<crate::jidx_reader::SeedDocument>(),
            TraceIndex::Owner(_) => std::mem::size_of::<crate::owner_format::OwnerDocument>(),
        });
        if capacity_bytes
            .saturating_add(document_bytes)
            .saturating_add(document_workspace)
            > reservation.bytes
        {
            postings_complete = false;
            continue;
        }
        let started = observed.then(Instant::now);
        let documents = index.seed_documents(seed)?;
        if observed {
            let occurrences = documents
                .iter()
                .map(|document| document.occurrence_count())
                .sum::<u64>();
            context_occurrence_histogram_log2[occurrences.ilog2().min(15) as usize] += 1;
        }
        membership_ns += started.map_or(0, |started| started.elapsed().as_nanos() as u64);
        let bytes = documents.capacity() * std::mem::size_of::<SeedDocument>();
        peak_capacity_bound = peak_capacity_bound.max(
            capacity_bytes
                .saturating_add(bytes)
                .saturating_add(document_workspace),
        );
        if capacity_bytes.saturating_add(bytes) > reservation.bytes {
            postings_complete = false;
            continue;
        }
        capacity_bytes += bytes;
        let position_bytes = documents
            .iter()
            .try_fold(0usize, |sum, doc| {
                usize::try_from(doc.occurrence_count())
                    .ok()?
                    .checked_mul(std::mem::size_of::<crate::jidx_reader::SeedOccurrence>())?
                    .checked_add(sum)
            })
            .and_then(|bytes| {
                bytes.checked_add(
                    documents.len()
                        * std::mem::size_of::<Vec<crate::jidx_reader::SeedOccurrence>>(),
                )
            });
        let occurrences = if position_bytes
            .is_some_and(|bytes| capacity_bytes.saturating_add(bytes) <= reservation.bytes)
        {
            let started = observed.then(Instant::now);
            let values = documents
                .iter()
                .map(|&doc| index.seed_document_occurrences(seed, doc))
                .collect::<Result<Vec<_>, _>>()?;
            position_ns += started.map_or(0, |started| started.elapsed().as_nanos() as u64);
            let charge = values.capacity()
                * std::mem::size_of::<Vec<crate::jidx_reader::SeedOccurrence>>()
                + values
                    .iter()
                    .map(|v| {
                        v.capacity() * std::mem::size_of::<crate::jidx_reader::SeedOccurrence>()
                    })
                    .sum::<usize>();
            peak_capacity_bound = peak_capacity_bound.max(capacity_bytes.saturating_add(charge));
            if capacity_bytes.saturating_add(charge) <= reservation.bytes {
                capacity_bytes += charge;
                Some(values)
            } else {
                None
            }
        } else {
            None
        };
        postings[ordinal] = Some(BatchPosting {
            documents,
            occurrences,
        });
    }
    if index.cache_file_identity()? != Some(identity) {
        return Err(TraceError::Invalid("shared seed lookup identity"));
    }
    let usage = observed
        .then(|| PostingUsage::new(&postings, reservation.bytes - capacity_bytes))
        .flatten();
    capacity_bytes += usage.as_ref().map_or(0, PostingUsage::capacity_bytes);
    peak_capacity_bound = peak_capacity_bound.max(capacity_bytes);
    reservation.retain(capacity_bytes);
    Ok(Some(SharedSeedLookups {
        identity,
        entries,
        query_entries: keys,
        attempted_keys,
        query_ranges,
        postings,
        postings_complete,
        capacity_bytes,
        peak_capacity_bound,
        lookup_ns,
        phase_postings_ns: phase_elapsed(posting_phase),
        membership_ns,
        position_ns,
        posting_execution: execution.unwrap_or_default(),
        usage,
        distinct_cores,
        split_core_resolutions,
        context_reuse_histogram_log2,
        context_occurrence_histogram_log2,
        lookup_tasks,
        lookup_plan_hash: if observed { lookup_plan_hash } else { 0 },
        lookup_dispatch_ns,
        lookup_parallel_ns,
        lookup_compute_ns,
        lookup_reduce_ns,
        lookup_dispatch_to_start_ns,
        _reservation: reservation,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn core_lookup_tasks_cover_keys_and_keep_prefix_groups_whole() {
        let reference = |keys: &[u32]| {
            // The former split: whole batches ending before a shared prefix group.
            let mut ranges = Vec::new();
            let mut start = 0;
            while start < keys.len() {
                let mut end = (start + SEED_LOOKUP_BATCH_KEYS).min(keys.len());
                if end < keys.len() {
                    let prefix = keys[end] >> 14;
                    while end > start && keys[end - 1] >> 14 == prefix {
                        end -= 1;
                    }
                    if end == start {
                        end = (start + SEED_LOOKUP_BATCH_KEYS).min(keys.len());
                    }
                }
                ranges.push(start..end);
                start = end;
            }
            ranges
        };
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 34) as u32
        };
        let sorted = |mut keys: Vec<u32>| {
            keys.sort_unstable();
            keys.dedup();
            keys
        };
        let sparse = sorted((0..11_839).map(|_| next()).collect());
        let dense = sorted((0..2_200_000).map(|_| next()).collect());
        assert!(dense.len() >= CORE_LOOKUP_TASKS * SEED_LOOKUP_BATCH_KEYS);
        let one_prefix = (0..40_000).map(|core| 7 << 14 | core).collect::<Vec<_>>();
        let straddling = sorted(
            (0..300)
                .map(|core| 1 << 14 | core)
                .chain((0..300).map(|core| 2 << 14 | core))
                .chain((0..20_000).map(|_| next()))
                .collect(),
        );
        for keys in [
            vec![],
            vec![0, (1 << 30) - 1],
            sparse,
            dense,
            one_prefix,
            straddling,
        ] {
            let ranges = core_prefix_ranges(&keys).collect::<Vec<_>>();
            assert_eq!(
                ranges.iter().map(|range| range.len()).sum::<usize>(),
                keys.len()
            );
            for pair in ranges.windows(2) {
                assert_eq!(pair[0].end, pair[1].start);
            }
            for range in &ranges {
                assert!(!range.is_empty() && range.len() <= SEED_LOOKUP_BATCH_KEYS);
                // A group spans two tasks only when it alone exceeds one batch.
                if range.end < keys.len() && keys[range.end - 1] >> 14 == keys[range.end] >> 14 {
                    let prefix = keys[range.end] >> 14;
                    assert!(
                        keys.iter().filter(|key| *key >> 14 == prefix).count()
                            > SEED_LOOKUP_BATCH_KEYS
                    );
                }
            }
            if keys.len() <= MIN_CORE_LOOKUP_TASK_KEYS {
                assert!(ranges.len() <= 1);
            } else if keys.len() <= CORE_LOOKUP_TASKS * SEED_LOOKUP_BATCH_KEYS {
                assert!(ranges.len() > 1 || keys.iter().all(|key| key >> 14 == keys[0] >> 14));
            }
            if keys.len() >= CORE_LOOKUP_TASKS * SEED_LOOKUP_BATCH_KEYS {
                assert_eq!(ranges, reference(&keys));
            }
        }
        assert_eq!(
            core_prefix_ranges(&sorted((0..11_839).map(|core| core * 90_000).collect())).count(),
            47
        );
    }

    #[test]
    fn worker_cpu_clock_is_monotone_when_available() {
        let Some(start) = worker_cpu_ns() else {
            return;
        };
        let mut value = 1u64;
        for i in 0..10_000 {
            value = std::hint::black_box(value.wrapping_mul(3) ^ i);
        }
        std::hint::black_box(value);
        assert!(worker_cpu_ns().is_some_and(|end| end > start));
    }

    #[test]
    fn lookup_task_ranges_are_identical_across_workers_and_bound_heavy_cores() {
        let keys = [(1, 20_000), (2, 40_000), (3, 20_000)]
            .into_iter()
            .flat_map(|(core, count)| {
                (0..count).map(move |context| {
                    crate::shared_seed::SharedKey {
                        core,
                        context,
                        length: 31,
                    }
                    .packed()
                    .unwrap()
                })
            })
            .collect::<Vec<_>>();
        let expected = vec![
            (0..20_000, false),
            (20_000..52_768, true),
            (52_768..80_000, false),
        ];
        for workers in [1, 4, 8, 16] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap();
            for prefixes in [false, true] {
                let actual = pool.install(|| lookup_ranges(&keys, prefixes).collect::<Vec<_>>());
                assert_eq!(actual, expected);
                assert_eq!(
                    actual.iter().map(|(range, _)| range.len()).sum::<usize>(),
                    keys.len()
                );
                assert!(
                    actual
                        .iter()
                        .all(|(range, _)| range.len() <= SEED_LOOKUP_BATCH_KEYS)
                );
            }
        }
    }

    #[test]
    fn prefix_tasks_keep_requested_prefixes_together_across_workers() {
        let cores = [(0, 16_000), (1, 12_000), (2, 14_000), (3, 13_000)]
            .into_iter()
            .flat_map(|(prefix, count)| (0..count).map(move |low| (prefix << 14) | low))
            .collect::<Vec<_>>();
        let keys = cores
            .iter()
            .map(|&core| u64::from(core))
            .collect::<Vec<_>>();
        // Core lookup tasks target 1/64 of the keys and keep each prefix whole; context lookups
        // still fill whole batches.
        let core_expected = vec![0..16_000, 16_000..28_000, 28_000..42_000, 42_000..55_000];
        let context_expected = vec![0..28_000, 28_000..55_000];
        for workers in [1, 4, 8, 16] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap();
            pool.install(|| {
                assert_eq!(
                    core_prefix_ranges(&cores).collect::<Vec<_>>(),
                    core_expected
                );
                assert_eq!(
                    lookup_ranges(&keys, true)
                        .map(|(range, split)| {
                            assert!(!split);
                            range
                        })
                        .collect::<Vec<_>>(),
                    context_expected
                );
            });
        }
    }
}
