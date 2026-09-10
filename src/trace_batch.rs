use crate::bgzf_cache::BgzfBlockCache;
use crate::jidx_reader::{SEED_LOOKUP_BATCH_KEYS, SeedEntry};
use crate::trace::{CacheReservation, LOOKUP_CACHE_AVAILABLE, TraceError};
use crate::trace_index::{TraceCacheIdentity, TraceDocument as SeedDocument, TraceIndex};
use std::collections::BTreeMap;
use std::sync::Arc;
use std::time::Instant;

pub(crate) struct SharedSeedLookups {
    pub(crate) identity: TraceCacheIdentity,
    pub(crate) entries: Vec<(u64, Option<SeedEntry>)>,
    pub(crate) postings: BTreeMap<u64, BatchPosting>,
    pub(crate) postings_complete: bool,
    pub(crate) capacity_bytes: usize,
    pub(crate) lookup_ns: u64,
    pub(crate) membership_ns: u64,
    pub(crate) position_ns: u64,
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

pub(crate) fn lookup_budget(index: &TraceIndex) -> usize {
    if index.is_shared() {
        128 * 1024 * 1024
    } else {
        32 * 1024 * 1024
    }
}

pub(crate) fn prepare_lookup(
    index: &TraceIndex,
    mut keys: Vec<u64>,
    observed: bool,
) -> Result<Option<SharedSeedLookups>, TraceError> {
    let started = observed.then(Instant::now);
    let budget = lookup_budget(index);
    let key_count = Some(keys.len());
    let row_bytes = std::mem::size_of::<u64>() + std::mem::size_of::<(u64, Option<SeedEntry>)>();
    let Some(bytes) = key_count.and_then(|keys| keys.checked_mul(row_bytes)?.checked_add(4096))
    else {
        return Ok(None);
    };
    if bytes > budget {
        return Ok(None);
    }
    let Some(reservation) = CacheReservation::acquire(&LOOKUP_CACHE_AVAILABLE, budget) else {
        return Ok(None);
    };
    let Some(identity) = index.cache_file_identity()? else {
        return Ok(None);
    };
    keys.sort_unstable();
    keys.dedup();
    if index.is_shared() {
        keys.sort_unstable_by_key(|key| (std::cmp::Reverse(key >> 62), *key));
    }
    let mut entries = Vec::new();
    if entries.try_reserve_exact(keys.len()).is_err() {
        return Ok(None);
    }
    for chunk in keys.chunks(SEED_LOOKUP_BATCH_KEYS) {
        entries.extend(chunk.iter().copied().zip(index.find_seeds_batch(chunk)?));
    }
    let lookup_ns = started.map_or(0, |started| started.elapsed().as_nanos() as u64);
    let mut membership_ns = 0;
    let mut position_ns = 0;
    let mut capacity_bytes = 4096
        + keys.capacity() * std::mem::size_of::<u64>()
        + entries.capacity() * std::mem::size_of::<(u64, Option<SeedEntry>)>();
    let mut postings = BTreeMap::new();
    let mut postings_complete = true;
    for &(key, seed) in &entries {
        let Some(seed) = seed else {
            continue;
        };
        let started = observed.then(Instant::now);
        let documents = index.seed_documents(seed)?;
        membership_ns += started.map_or(0, |started| started.elapsed().as_nanos() as u64);
        let bytes = 128 + documents.capacity() * std::mem::size_of::<SeedDocument>();
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
    entries.sort_unstable_by_key(|entry| entry.0);
    if index.cache_file_identity()? != Some(identity) {
        return Err(TraceError::Invalid("shared seed lookup identity"));
    }
    Ok(Some(SharedSeedLookups {
        identity,
        entries,
        postings,
        postings_complete,
        capacity_bytes,
        lookup_ns,
        membership_ns,
        position_ns,
        _reservation: reservation,
    }))
}
