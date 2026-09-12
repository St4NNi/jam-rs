use crate::bgzf_cache::BgzfBlockCache;
use crate::jidx_reader::SEED_LOOKUP_BATCH_KEYS;
use crate::trace::{CacheReservation, LOOKUP_CACHE_AVAILABLE, TraceError};
use crate::trace_index::{
    TraceCacheIdentity, TraceDocument as SeedDocument, TraceIndex, TraceSeed,
};
use rayon::prelude::*;
use std::ops::Range;
use std::sync::Arc;
use std::time::Instant;

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
    pub(crate) membership_ns: u64,
    pub(crate) position_ns: u64,
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

impl SharedSeedLookups {
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
    keys.try_reserve_exact(request_count)
        .map_err(|_| TraceError::Invalid("core request allocation"))?;
    if keys.capacity() > request_count {
        return Ok(None);
    }
    for key in requests {
        if keys.len() == request_count {
            return Err(TraceError::Invalid("core request count"));
        }
        keys.push(key);
    }
    if keys.len() != request_count {
        return Err(TraceError::Invalid("core request count"));
    }
    keys.par_sort_unstable();
    keys.dedup();
    let tasks = if reader.has_core_prefixes() {
        core_prefix_ranges(&keys).count()
    } else {
        keys.len().div_ceil(SEED_LOOKUP_BATCH_KEYS)
    };
    let Some(bytes) = keys
        .capacity()
        .checked_mul(std::mem::size_of::<u32>())
        .and_then(|bytes| {
            bytes.checked_add(
                keys.len()
                    .checked_mul(2 * std::mem::size_of::<crate::shared_reader::SharedGroup>())?,
            )
        })
        .and_then(|bytes| {
            bytes.checked_add(tasks.checked_mul(std::mem::size_of::<
                Result<Vec<crate::shared_reader::SharedGroup>, TraceError>,
            >())?)
        })
        .and_then(|bytes| {
            bytes.checked_add(if reader.has_core_prefixes() {
                0
            } else {
                rayon::current_num_threads()
                    .min(tasks)
                    .checked_mul(SEED_LOOKUP_BATCH_KEYS)?
                    .checked_mul(std::mem::size_of::<crate::shared_seed::SharedKey>())?
            })
        })
        .and_then(|bytes| bytes.checked_add(4096 + std::mem::size_of::<SharedCoreLookups>()))
        .and_then(|bytes| {
            bytes.checked_add(tasks.checked_mul(std::mem::size_of::<Range<usize>>())?)
        })
        .filter(|&bytes| bytes <= lookup_budget(index))
    else {
        return Ok(None);
    };
    let Some(mut reservation) =
        CacheReservation::acquire(&LOOKUP_CACHE_AVAILABLE, bytes - key_reservation.bytes)
    else {
        return Ok(None);
    };
    let ranges = if reader.has_core_prefixes() {
        let mut ranges = Vec::new();
        ranges
            .try_reserve_exact(tasks)
            .map_err(|_| TraceError::Invalid("core task allocation"))?;
        if ranges.capacity() > tasks {
            return Ok(None);
        }
        ranges.extend(core_prefix_ranges(&keys));
        Some(ranges)
    } else {
        None
    };
    let lookup = |chunk: &[u32]| {
        let mut groups = Vec::new();
        if reader.has_core_prefixes() {
            reader.resolve_sorted_cores_into(chunk, &mut groups)?;
            return Ok::<_, TraceError>(groups);
        }
        let contexts = chunk
            .iter()
            .copied()
            .map(crate::shared_seed::SharedKey::core)
            .collect::<Vec<_>>();
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
    let chunks = if let Some(ranges) = &ranges {
        ranges
            .par_iter()
            .map(|range| lookup(&keys[range.clone()]))
            .collect::<Vec<_>>()
    } else {
        keys.par_chunks(SEED_LOOKUP_BATCH_KEYS)
            .map(lookup)
            .collect::<Vec<_>>()
    };
    let resource_limit = |error: &TraceError| {
        matches!(
            error,
            TraceError::Shared(crate::shared_format::SharedError::ResourceLimit)
        )
    };
    if chunks
        .iter()
        .any(|chunk| chunk.as_ref().is_err_and(resource_limit))
    {
        for chunk in chunks {
            if let Err(error) = chunk
                && !resource_limit(&error)
            {
                return Err(error);
            }
        }
        return Ok(None);
    }
    let count = chunks.iter().try_fold(0usize, |sum, chunk| {
        chunk.as_ref().map(|chunk| sum + chunk.len())
    });
    let count = match count {
        Ok(count) => count,
        Err(_) => {
            for chunk in chunks {
                chunk?;
            }
            unreachable!()
        }
    };
    let mut groups = Vec::new();
    if groups.try_reserve_exact(count).is_err() {
        return Ok(None);
    }
    if groups.capacity() > count {
        return Ok(None);
    }
    for chunk in chunks {
        groups.extend(chunk?);
    }
    drop(keys);
    drop(key_reservation);
    drop(ranges);
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
    + std::mem::size_of::<(u64, Option<TraceSeed>)>();

fn lookup_chunk_keys() -> usize {
    SEED_LOOKUP_BATCH_KEYS
}

fn core_prefix_ranges(keys: &[u32]) -> impl Iterator<Item = Range<usize>> + '_ {
    let mut start = 0;
    std::iter::from_fn(move || {
        if start == keys.len() {
            return None;
        }
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
            tasks.push((weight, chunk, slots));
        }
        tasks.sort_unstable_by_key(|task| std::cmp::Reverse(task.0));
        lookup_tasks = tasks.len();
        lookup_dispatch_ns = dispatch.map_or(0, |start| start.elapsed().as_nanos() as u64);
        let dispatch = observed.then(Instant::now);
        let timings = tasks
            .into_par_iter()
            .map(|(_, keys, slots)| {
                let start = observed.then(Instant::now);
                let queued = start.zip(dispatch).map_or(0, |(start, dispatch)| {
                    start.duration_since(dispatch).as_nanos() as u64
                });
                let seeds = if let Some(cores) = cores {
                    index.find_seeds_in_cores(keys, cores)?
                } else {
                    index.find_seeds_batch(keys)?
                };
                let compute = start.map_or(0, |start| start.elapsed().as_nanos() as u64);
                let start = observed.then(Instant::now);
                for (slot, seed) in slots.iter_mut().zip(seeds) {
                    slot.1 = seed;
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
    if slot_bytes.is_none_or(|bytes| capacity_bytes.saturating_add(bytes) > reservation.bytes) {
        return Ok(None);
    }
    let mut postings = Vec::new();
    postings
        .try_reserve_exact(entries.len())
        .map_err(|_| TraceError::Invalid("posting slots allocation"))?;
    let slot_bytes = postings.capacity() * std::mem::size_of::<Option<BatchPosting>>();
    if capacity_bytes.saturating_add(slot_bytes) > reservation.bytes {
        return Ok(None);
    }
    capacity_bytes += slot_bytes;
    peak_capacity_bound = peak_capacity_bound.max(capacity_bytes);
    postings.resize_with(entries.len(), || None);
    let mut postings_complete = true;
    let mut context_occurrence_histogram_log2 = [0; 16];
    for (ordinal, &(_, seed)) in entries.iter().enumerate() {
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
        membership_ns,
        position_ns,
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
        let expected = vec![0..28_000, 28_000..55_000];
        for workers in [1, 4, 8, 16] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap();
            pool.install(|| {
                assert_eq!(core_prefix_ranges(&cores).collect::<Vec<_>>(), expected);
                assert_eq!(
                    lookup_ranges(&keys, true)
                        .map(|(range, split)| {
                            assert!(!split);
                            range
                        })
                        .collect::<Vec<_>>(),
                    expected
                );
            });
        }
    }
}
