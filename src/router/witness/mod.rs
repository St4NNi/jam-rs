//! Query-side witness selection for the collection router.
//!
//! The router stores one dense k=21 key set.  Coarser levels are views over
//! that same set: a witness carries the scales for which its already computed
//! `jamhash_u64_v1` passes the threshold.  In particular, this module never
//! makes one copy of a key per tier.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap};

use needletail::Sequence;
use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::router::{QueryWitness, RouterContractError, WITNESS_K, WitnessKey, WitnessScheme};

/// The fixed k-mer size of collection witnesses.
pub const WITNESS_K_USIZE: usize = WITNESS_K as usize;

/// Window sizes used by the short-trace evaluation matrix.
pub const DEFAULT_QUERY_WINDOW_SIZES: [u32; 4] = [128, 256, 512, 1024];

/// A witness occurrence retained once while its nested tier membership is
/// represented by `retained_scales`.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct NestedWitness {
    pub key: WitnessKey,
    pub query_position: u64,
    /// Whether the canonical packed k-mer is read from the reverse query
    /// orientation at `query_position`.
    pub query_reverse: bool,
    pub query_window_ids: Vec<u32>,
    pub retained_scales: Vec<u32>,
}

impl NestedWitness {
    /// Convert this occurrence to the shared handoff type.
    #[must_use]
    pub fn query_witness(&self) -> QueryWitness {
        QueryWitness {
            key: self.key,
            query_position: self.query_position,
            query_reverse: self.query_reverse,
            query_window_ids: self.query_window_ids.clone(),
        }
    }

    /// Whether this occurrence is available at `scale` without recomputing
    /// or duplicating its packed key.
    #[must_use]
    pub fn retained_at(&self, scale: u32) -> bool {
        self.retained_scales.binary_search(&scale).is_ok()
    }
}

/// All dense query witnesses and their nested tier membership.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct NestedWitnessSet {
    pub sequence_length: u64,
    pub window_size: u32,
    pub scheme: WitnessScheme,
    pub witnesses: Vec<NestedWitness>,
}

impl NestedWitnessSet {
    /// Number of logical query windows, including windows that contain only
    /// ambiguous bases.  Callers can use `eligible_window_count` when they
    /// need to exclude a final short window.
    #[must_use]
    pub fn query_window_count(&self) -> u32 {
        query_window_count(self.sequence_length, self.window_size)
    }

    /// Number of windows that are long enough to contain a complete k=21
    /// witness.  Ambiguous content remains eligible; it is evidence-free but
    /// should not silently change the requested query geometry.
    #[must_use]
    pub fn eligible_window_count(&self) -> u32 {
        eligible_window_count(self.sequence_length, self.window_size, WITNESS_K_USIZE)
    }

    /// Return the same query occurrences filtered to one available tier.
    pub fn at_scale(&self, scale: u32) -> Result<Vec<QueryWitness>, WitnessError> {
        self.scheme.includes_hash(1, scale)?;
        Ok(self
            .witnesses
            .iter()
            .filter(|witness| witness.retained_at(scale))
            .map(NestedWitness::query_witness)
            .collect())
    }

    /// Return the unique packed keys represented by this set.  Repeated
    /// occurrences in the query do not duplicate collection-router keys.
    #[must_use]
    pub fn unique_keys(&self) -> Vec<WitnessKey> {
        let mut keys: Vec<_> = self.witnesses.iter().map(|witness| witness.key).collect();
        keys.sort_unstable();
        keys.dedup();
        keys
    }
}

/// One local bottom-r occurrence.  The `query_window_ids` field is merged when
/// the same occurrence is selected by adjacent/overlapping windows.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct LocalBottomWitness {
    pub key: WitnessKey,
    pub query_position: u64,
    pub query_reverse: bool,
    pub query_window_ids: Vec<u32>,
}

impl LocalBottomWitness {
    #[must_use]
    pub fn query_witness(&self) -> QueryWitness {
        QueryWitness {
            key: self.key,
            query_position: self.query_position,
            query_reverse: self.query_reverse,
            query_window_ids: self.query_window_ids.clone(),
        }
    }
}

/// Result of local bottom-r selection, including geometry useful for candidate
/// evidence reporting.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct LocalBottomWitnessSet {
    pub sequence_length: u64,
    pub window_size: u32,
    pub r: usize,
    pub eligible_window_count: u32,
    pub witnesses: Vec<LocalBottomWitness>,
}

impl LocalBottomWitnessSet {
    #[must_use]
    pub fn query_witnesses(&self) -> Vec<QueryWitness> {
        self.witnesses
            .iter()
            .map(LocalBottomWitness::query_witness)
            .collect()
    }
}

/// One member of a fixed bottom-k sketch.  This is intentionally a bounded
/// comparison implementation rather than the default nested witness path.
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct FixedBottomWitness {
    pub key: WitnessKey,
    pub query_position: u64,
    pub query_reverse: bool,
    pub query_window_ids: Vec<u32>,
}

impl FixedBottomWitness {
    #[must_use]
    pub fn query_witness(&self) -> QueryWitness {
        QueryWitness {
            key: self.key,
            query_position: self.query_position,
            query_reverse: self.query_reverse,
            query_window_ids: self.query_window_ids.clone(),
        }
    }
}

/// Witness extraction and local-selection failures.
#[derive(Debug, Error, Eq, PartialEq)]
pub enum WitnessError {
    #[error(transparent)]
    Contract(#[from] RouterContractError),
    #[error("query window size must be greater than zero")]
    ZeroWindowSize,
    #[error("local bottom-r must be greater than zero")]
    ZeroBottomR,
    #[error("fixed bottom-k must be greater than zero")]
    ZeroBottomK,
    #[error("query window id does not fit in u32")]
    WindowIdOverflow,
}

/// Extract one dense nested k=21 witness set.
pub fn extract_nested_witnesses(
    sequence: &[u8],
    scheme: &WitnessScheme,
    window_size: u32,
) -> Result<NestedWitnessSet, WitnessError> {
    scheme.validate()?;
    validate_window_size(window_size)?;

    let sequence_length = u64::try_from(sequence.len()).unwrap_or(u64::MAX);
    let mut witnesses = Vec::new();
    for (position, key, query_reverse) in valid_observations(sequence) {
        let query_window_ids = query_window_ids(position, WITNESS_K_USIZE, window_size)?;
        let retained_scales = scheme.retained_scales(key.jamhash)?;
        if retained_scales.is_empty() {
            continue;
        }
        witnesses.push(NestedWitness {
            key,
            query_position: position,
            query_reverse,
            query_window_ids,
            retained_scales,
        });
    }

    Ok(NestedWitnessSet {
        sequence_length,
        window_size,
        scheme: scheme.clone(),
        witnesses,
    })
}

/// Label a k-mer occurrence with every logical query window it overlaps.
/// Window IDs are zero-based and the windows are half-open in query
/// coordinates.  Keeping this helper public lets candidate and alignment
/// stages use the exact same labels as extraction.
pub fn window_ids_for_position(
    position: u64,
    k: u8,
    window_size: u32,
) -> Result<Vec<u32>, WitnessError> {
    validate_window_size(window_size)?;
    query_window_ids(position, usize::from(k), window_size)
}

/// Compatibility-friendly alias for nested witness extraction.
pub fn build_nested_witnesses(
    sequence: &[u8],
    scheme: &WitnessScheme,
    window_size: u32,
) -> Result<NestedWitnessSet, WitnessError> {
    extract_nested_witnesses(sequence, scheme, window_size)
}

/// Return nested query witnesses for one tier without making a per-tier copy
/// during extraction.
pub fn nested_query_witnesses_at_scale(
    sequence: &[u8],
    scheme: &WitnessScheme,
    scale: u32,
    window_size: u32,
) -> Result<Vec<QueryWitness>, WitnessError> {
    extract_nested_witnesses(sequence, scheme, window_size)?.at_scale(scale)
}

/// Select a fixed-size bottom-k sketch from exact packed k=21 occurrences.
/// Only the current bounded set is retained; a high-ranked key is never
/// materialized merely because it appeared in the input.
pub fn fixed_bottom_k_witnesses(
    sequence: &[u8],
    bottom_k: usize,
    window_size: u32,
) -> Result<Vec<FixedBottomWitness>, WitnessError> {
    validate_window_size(window_size)?;
    if bottom_k == 0 {
        return Err(WitnessError::ZeroBottomK);
    }

    let mut selected: BTreeMap<(u64, u64), FixedBottomWitness> = BTreeMap::new();
    let mut heap = BinaryHeap::new();

    let mut extraction_error = None;
    for_each_valid_observation(sequence, |position, key, query_reverse| {
        if extraction_error.is_some() {
            return;
        }
        let rank = (key.jamhash, key.packed);
        if let Some(existing) = selected.get_mut(&rank) {
            match query_window_ids(position, WITNESS_K_USIZE, window_size) {
                Ok(windows) => merge_window_ids(&mut existing.query_window_ids, windows),
                Err(error) => extraction_error = Some(error),
            }
            return;
        }

        if selected.len() == bottom_k {
            let Some(worst) = heap.peek().copied() else {
                // The map and heap are kept in lockstep.  This branch only
                // protects the bounded comparison if the implementation is
                // changed later.
                return;
            };
            if rank >= worst {
                return;
            }
            heap.pop();
            selected.remove(&worst);
        }

        let windows = match query_window_ids(position, WITNESS_K_USIZE, window_size) {
            Ok(windows) => windows,
            Err(error) => {
                extraction_error = Some(error);
                return;
            }
        };
        selected.insert(
            rank,
            FixedBottomWitness {
                key,
                query_position: position,
                query_reverse,
                query_window_ids: windows,
            },
        );
        heap.push(rank);
    });

    if let Some(error) = extraction_error {
        return Err(error);
    }

    Ok(selected.into_values().collect())
}

/// Alias used by callers comparing the bounded fixed-bottom-k baseline.
pub fn extract_fixed_bottom_k(
    sequence: &[u8],
    bottom_k: usize,
    window_size: u32,
) -> Result<Vec<FixedBottomWitness>, WitnessError> {
    fixed_bottom_k_witnesses(sequence, bottom_k, window_size)
}

/// Short name for the bounded fixed-bottom-k comparison.
pub fn fixed_bottom_k(
    sequence: &[u8],
    bottom_k: usize,
    window_size: u32,
) -> Result<Vec<FixedBottomWitness>, WitnessError> {
    fixed_bottom_k_witnesses(sequence, bottom_k, window_size)
}

/// Select the r smallest exact witnesses independently in each logical query
/// window.  A rolling ordered set updates entering/leaving k-mers as windows
/// advance, so no window is sorted from scratch.
pub fn local_bottom_r_witnesses(
    sequence: &[u8],
    window_size: u32,
    r: usize,
) -> Result<LocalBottomWitnessSet, WitnessError> {
    validate_window_size(window_size)?;
    if r == 0 {
        return Err(WitnessError::ZeroBottomR);
    }

    let sequence_length = u64::try_from(sequence.len()).unwrap_or(u64::MAX);
    let kmer_positions = sequence.len().saturating_sub(WITNESS_K_USIZE - 1);
    let mut positions = vec![None; kmer_positions];
    for_each_valid_observation(sequence, |position, key, query_reverse| {
        if let Ok(position) = usize::try_from(position)
            && let Some(slot) = positions.get_mut(position)
        {
            *slot = Some(RankedObservation {
                position,
                key,
                query_reverse,
            });
        }
    });

    let windows = query_window_count(sequence_length, window_size);
    let eligible_window_count =
        eligible_window_count(sequence_length, window_size, WITNESS_K_USIZE);
    let mut rolling = RollingBottomR::default();
    let mut merged: BTreeMap<(u64, u64), LocalBottomWitness> = BTreeMap::new();

    for window_id in 0..windows {
        let window_start = u64::from(window_id) * u64::from(window_size);
        let window_end = window_start
            .saturating_add(u64::from(window_size))
            .min(sequence_length);
        if window_end.saturating_sub(window_start) < WITNESS_K as u64 {
            continue;
        }
        let start = usize::try_from(window_start).unwrap_or(positions.len());
        let end_base = usize::try_from(window_end).unwrap_or(sequence.len());
        let end = end_base
            .saturating_sub(WITNESS_K_USIZE - 1)
            .min(positions.len());
        rolling.set_range(&positions, start.min(end), end);
        for observation in rolling.smallest(r) {
            let id = (observation.key.packed, observation.position as u64);
            let witness = merged.entry(id).or_insert_with(|| LocalBottomWitness {
                key: observation.key,
                query_position: observation.position as u64,
                query_reverse: observation.query_reverse,
                query_window_ids: Vec::new(),
            });
            if witness.query_window_ids.last().copied() != Some(window_id) {
                witness.query_window_ids.push(window_id);
            }
        }
    }

    let mut witnesses: Vec<_> = merged.into_values().collect();
    witnesses.sort_unstable_by_key(|witness| (witness.query_position, witness.key));

    Ok(LocalBottomWitnessSet {
        sequence_length,
        window_size,
        r,
        eligible_window_count,
        witnesses,
    })
}

/// Return only local bottom-r query occurrences.
pub fn extract_local_bottom_r(
    sequence: &[u8],
    window_size: u32,
    r: usize,
) -> Result<Vec<QueryWitness>, WitnessError> {
    Ok(local_bottom_r_witnesses(sequence, window_size, r)?.query_witnesses())
}

/// Short name for local bottom-r extraction.
pub fn local_bottom_r(
    sequence: &[u8],
    window_size: u32,
    r: usize,
) -> Result<LocalBottomWitnessSet, WitnessError> {
    local_bottom_r_witnesses(sequence, window_size, r)
}

/// Compute the number of logical windows for a sequence.
#[must_use]
pub fn query_window_count(sequence_length: u64, window_size: u32) -> u32 {
    if window_size == 0 || sequence_length == 0 {
        return 0;
    }
    let count = sequence_length.saturating_add(u64::from(window_size) - 1) / u64::from(window_size);
    u32::try_from(count).unwrap_or(u32::MAX)
}

/// Count windows long enough to hold one complete k-mer.
#[must_use]
pub fn eligible_window_count(sequence_length: u64, window_size: u32, k: usize) -> u32 {
    if window_size == 0 || sequence_length < k as u64 {
        return 0;
    }
    let count = sequence_length / u64::from(window_size);
    let remainder = sequence_length % u64::from(window_size);
    let eligible = count + if remainder >= k as u64 { 1 } else { 0 };
    u32::try_from(eligible).unwrap_or(u32::MAX)
}

fn validate_window_size(window_size: u32) -> Result<(), WitnessError> {
    if window_size == 0 {
        Err(WitnessError::ZeroWindowSize)
    } else {
        Ok(())
    }
}

fn query_window_ids(position: u64, k: usize, window_size: u32) -> Result<Vec<u32>, WitnessError> {
    let first = position / u64::from(window_size);
    let last = position.saturating_add(k.saturating_sub(1) as u64) / u64::from(window_size);
    let first = u32::try_from(first).map_err(|_| WitnessError::WindowIdOverflow)?;
    let last = u32::try_from(last).map_err(|_| WitnessError::WindowIdOverflow)?;
    Ok((first..=last).collect())
}

fn valid_observations(sequence: &[u8]) -> Vec<(u64, WitnessKey, bool)> {
    let normalized = sequence.normalize(false);
    normalized
        .bit_kmers(WITNESS_K, true)
        .filter_map(|(position, kmer, reverse)| {
            WitnessKey::from_packed(kmer.0)
                .ok()
                .map(|key| (position as u64, key, reverse))
        })
        .collect()
}

fn for_each_valid_observation(sequence: &[u8], mut callback: impl FnMut(u64, WitnessKey, bool)) {
    // Normalize once so lower-case and U are handled like the existing trace
    // seed path.  The callback keeps fixed bottom-k selection streaming and
    // therefore bounded by the requested heap size.
    let normalized = sequence.normalize(false);
    for (position, kmer, reverse) in normalized.bit_kmers(WITNESS_K, true) {
        if let Ok(key) = WitnessKey::from_packed(kmer.0) {
            callback(position as u64, key, reverse);
        }
    }
}

fn merge_window_ids(target: &mut Vec<u32>, incoming: Vec<u32>) {
    for id in incoming {
        if !target.contains(&id) {
            target.push(id);
        }
    }
    target.sort_unstable();
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct RankedObservation {
    key: WitnessKey,
    position: usize,
    query_reverse: bool,
}

impl Ord for RankedObservation {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key
            .jamhash
            .cmp(&other.key.jamhash)
            .then_with(|| self.key.packed.cmp(&other.key.packed))
            .then_with(|| self.position.cmp(&other.position))
            .then_with(|| self.query_reverse.cmp(&other.query_reverse))
    }
}

impl PartialOrd for RankedObservation {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Default)]
struct RollingBottomR {
    current_start: usize,
    current_end: usize,
    ordered: BTreeSet<RankedObservation>,
}

impl RollingBottomR {
    fn set_range(&mut self, positions: &[Option<RankedObservation>], start: usize, end: usize) {
        while self.current_start < start {
            if let Some(observation) = positions.get(self.current_start).copied().flatten() {
                self.ordered.remove(&observation);
            }
            self.current_start += 1;
        }
        // K-mer starts in the gap between the previous window's final
        // complete k-mer and this window's first base overlap a boundary and
        // belong to neither non-overlapping logical window.
        while self.current_end < start {
            self.current_end += 1;
        }
        while self.current_end < end {
            if let Some(observation) = positions.get(self.current_end).copied().flatten() {
                self.ordered.insert(observation);
            }
            self.current_end += 1;
        }
        // A caller can request an empty range after a previous non-empty
        // range.  Remove any positions that now lie beyond the end while
        // preserving the monotonic rolling state.
        while self.current_end > end {
            self.current_end -= 1;
            if let Some(observation) = positions.get(self.current_end).copied().flatten() {
                self.ordered.remove(&observation);
            }
        }
    }

    fn smallest(&self, r: usize) -> impl Iterator<Item = RankedObservation> + '_ {
        self.ordered.iter().copied().take(r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::router::{HashAlgorithmId, WitnessScheme};

    fn scheme() -> WitnessScheme {
        WitnessScheme {
            scheme_id: 1,
            k: WITNESS_K,
            base_scale: 20,
            available_scales: vec![20, 50, 100, 200, 500],
            hash_id: HashAlgorithmId::JamhashU64V1,
            zero_excluded: true,
        }
    }

    fn sequence(length: usize) -> Vec<u8> {
        let alphabet = b"ACGT";
        let mut state = 0x1234_5678_9abc_def0_u64;
        (0..length)
            .map(|_| {
                state ^= state << 7;
                state ^= state >> 9;
                state ^= state << 8;
                alphabet[(state as usize) % alphabet.len()]
            })
            .collect()
    }

    fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
        sequence
            .iter()
            .rev()
            .map(|base| match base {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => b'N',
            })
            .collect()
    }

    #[test]
    fn nested_set_stores_each_occurrence_once_and_derives_tiers() {
        let set = extract_nested_witnesses(&sequence(800), &scheme(), 128).unwrap();
        assert!(!set.witnesses.is_empty());
        assert!(
            set.witnesses
                .iter()
                .all(|witness| witness.retained_scales.windows(2).all(|p| p[0] < p[1]))
        );
        assert!(
            set.witnesses
                .iter()
                .all(|witness| witness.retained_scales[0] == 20)
        );
        assert_eq!(set.at_scale(20).unwrap().len(), set.witnesses.len());
        assert!(set.at_scale(500).unwrap().len() <= set.witnesses.len());
        assert_eq!(set.unique_keys().len(), {
            let mut keys = set
                .witnesses
                .iter()
                .map(|witness| witness.key)
                .collect::<Vec<_>>();
            keys.sort_unstable();
            keys.dedup();
            keys.len()
        });
    }

    #[test]
    fn reverse_complement_has_same_exact_canonical_keys() {
        let forward = sequence(600);
        let reverse = reverse_complement(&forward);
        let left = extract_nested_witnesses(&forward, &scheme(), 256).unwrap();
        let right = extract_nested_witnesses(&reverse, &scheme(), 256).unwrap();
        let mut left_keys = left.unique_keys();
        let mut right_keys = right.unique_keys();
        left_keys.sort_unstable();
        right_keys.sort_unstable();
        assert_eq!(left_keys, right_keys);
    }

    #[test]
    fn zero_hash_and_ambiguous_windows_are_not_witnesses() {
        let mut data = vec![b'A'; 21];
        data.extend_from_slice(b"NNNNNNNNNNNNNNNNNNNNNN");
        data.extend_from_slice(&sequence(200));
        let set = extract_nested_witnesses(&data, &scheme(), 128).unwrap();
        assert!(set.witnesses.iter().all(|witness| witness.key.packed != 0));
        assert!(set.witnesses.iter().all(|witness| witness.key.jamhash != 0));
    }

    #[test]
    fn local_bottom_r_handles_requested_window_sizes() {
        let data = sequence(2_100);
        for window in DEFAULT_QUERY_WINDOW_SIZES {
            let result = local_bottom_r_witnesses(&data, window, 2).unwrap();
            assert_eq!(result.window_size, window);
            assert_eq!(result.r, 2);
            assert!(
                result
                    .witnesses
                    .iter()
                    .all(|witness| !witness.query_window_ids.is_empty())
            );
        }
    }

    #[test]
    fn local_bottom_r_uses_ranked_set_and_fixed_bottom_is_bounded() {
        let data = sequence(10_000);
        let local = local_bottom_r_witnesses(&data, 128, 4).unwrap();
        assert!(local.witnesses.len() <= local.eligible_window_count as usize * 4);
        let fixed = fixed_bottom_k_witnesses(&data, 8, 128).unwrap();
        assert_eq!(fixed.len(), 8);
        let ranks = fixed
            .iter()
            .map(|witness| (witness.key.jamhash, witness.key.packed))
            .collect::<Vec<_>>();
        assert!(ranks.windows(2).all(|pair| pair[0] < pair[1]));
    }
}
