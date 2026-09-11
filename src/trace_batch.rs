use crate::bgzf_cache::BgzfBlockCache;
use crate::jidx_reader::SEED_LOOKUP_BATCH_KEYS;
use crate::trace::{CacheReservation, LOOKUP_CACHE_AVAILABLE, TraceError};
use crate::trace_index::{
    TraceCacheIdentity, TraceDocument as SeedDocument, TraceIndex, TraceSeed,
};
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::ops::Range;
use std::sync::Arc;
use std::time::Instant;

pub(crate) struct SharedSeedLookups {
    pub(crate) identity: TraceCacheIdentity,
    pub(crate) entries: Vec<(u64, Option<TraceSeed>)>,
    pub(crate) query_entries: Vec<u64>,
    pub(crate) attempted_keys: usize,
    pub(crate) query_ranges: Vec<Range<usize>>,
    pub(crate) postings: BTreeMap<u64, BatchPosting>,
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
    mut keys: Vec<u32>,
    observed: bool,
) -> Result<Option<Arc<SharedCoreLookups>>, TraceError> {
    let TraceIndex::Shared(reader) = index else {
        return Ok(None);
    };
    let started = observed.then(Instant::now);
    let Some(identity) = index.cache_file_identity()? else {
        return Ok(None);
    };
    keys.par_sort_unstable();
    keys.dedup();
    let tasks = keys.len().div_ceil(SEED_LOOKUP_BATCH_KEYS);
    let Some(bytes) = keys
        .capacity()
        .checked_mul(std::mem::size_of::<u32>())
        .and_then(|bytes| {
            bytes.checked_add(keys.len().checked_mul(
                std::mem::size_of::<Option<crate::shared_reader::SharedGroup>>()
                    + std::mem::size_of::<crate::shared_reader::SharedGroup>(),
            )?)
        })
        .and_then(|bytes| {
            bytes.checked_add(tasks.checked_mul(std::mem::size_of::<
                Result<Vec<Option<crate::shared_reader::SharedGroup>>, TraceError>,
            >())?)
        })
        .and_then(|bytes| {
            bytes.checked_add(
                rayon::current_num_threads()
                    .min(tasks)
                    .checked_mul(SEED_LOOKUP_BATCH_KEYS)?
                    .checked_mul(std::mem::size_of::<crate::shared_seed::SharedKey>())?,
            )
        })
        .and_then(|bytes| bytes.checked_add(4096 + std::mem::size_of::<SharedCoreLookups>()))
        .filter(|&bytes| bytes <= lookup_budget(index))
    else {
        return Ok(None);
    };
    let Some(mut reservation) = CacheReservation::acquire(&LOOKUP_CACHE_AVAILABLE, bytes) else {
        return Ok(None);
    };
    let chunks = keys
        .par_chunks(SEED_LOOKUP_BATCH_KEYS)
        .map(|chunk| {
            let contexts = chunk
                .iter()
                .copied()
                .map(crate::shared_seed::SharedKey::core)
                .collect::<Vec<_>>();
            reader.find_many(&contexts).map_err(TraceError::from)
        })
        .collect::<Vec<_>>();
    let count = chunks.iter().try_fold(0usize, |sum, chunk| {
        chunk
            .as_ref()
            .map(|chunk| sum + chunk.iter().filter(|group| group.is_some()).count())
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
    groups
        .try_reserve_exact(count)
        .map_err(|_| TraceError::Invalid("resolved core allocation"))?;
    for chunk in chunks {
        groups.extend(chunk?.into_iter().flatten());
    }
    drop(keys);
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
    (SEED_LOOKUP_BATCH_KEYS / rayon::current_num_threads()).max(1)
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
        let mut start = 0;
        while start < keys.len() {
            let mut end = (start + limit).min(keys.len());
            if end < keys.len() && core(keys[end - 1]) == core(keys[end]) {
                let boundary = end;
                while end > start && core(keys[end - 1]) == core(keys[boundary]) {
                    end -= 1;
                }
                if end == start {
                    end = boundary;
                    split_core_resolutions += 1;
                }
            }
            let chunk = &keys[start..end];
            let cores = 1 + chunk
                .windows(2)
                .filter(|pair| core(pair[0]) != core(pair[1]))
                .count();
            let weight = cores as u64 * u64::from(reader.core_count().max(1).ilog2() + 1)
                + chunk.len() as u64;
            let (slots, tail) = remaining.split_at_mut(chunk.len());
            remaining = tail;
            tasks.push((weight, chunk, slots));
            start = end;
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
    let mut postings = BTreeMap::new();
    let mut postings_complete = true;
    let mut context_occurrence_histogram_log2 = [0; 16];
    for &(key, seed) in &entries {
        let Some(seed) = seed else {
            continue;
        };
        let document_bytes = 128usize.saturating_add(
            (seed.document_frequency() as usize)
                .saturating_mul(std::mem::size_of::<SeedDocument>()),
        );
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
        let bytes = 128 + documents.capacity() * std::mem::size_of::<SeedDocument>();
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
        postings.insert(
            key,
            BatchPosting {
                documents,
                occurrences,
            },
        );
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
        lookup_dispatch_ns,
        lookup_parallel_ns,
        lookup_compute_ns,
        lookup_reduce_ns,
        lookup_dispatch_to_start_ns,
        _reservation: reservation,
    }))
}
