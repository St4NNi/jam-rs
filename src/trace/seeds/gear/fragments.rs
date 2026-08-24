//! Content-defined Gear fragments and exact-run verification.

use super::tables::{GearHasher, GearTable, GearTableKind, normalized_byte};
use std::collections::BTreeMap;

/// Recovered fine/coarse Gear stream settings from `c671736`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct GearStream {
    pub min_size: u32,
    pub target_size: u32,
    pub max_size: u32,
}

impl GearStream {
    #[must_use]
    pub const fn fine() -> Self {
        Self {
            min_size: 64,
            target_size: 192,
            max_size: 384,
        }
    }

    #[must_use]
    pub const fn coarse() -> Self {
        Self {
            min_size: 256,
            target_size: 768,
            max_size: 1_536,
        }
    }
}

/// Parameters for one content-defined stream.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct GearConfig {
    pub min_size: u32,
    pub target_size: u32,
    pub max_size: u32,
    pub table_kind: GearTableKind,
    pub table_seed: u64,
}

impl GearConfig {
    #[must_use]
    pub const fn fine(table_kind: GearTableKind, table_seed: u64) -> Self {
        let stream = GearStream::fine();
        Self {
            min_size: stream.min_size,
            target_size: stream.target_size,
            max_size: stream.max_size,
            table_kind,
            table_seed,
        }
    }

    #[must_use]
    pub const fn coarse(table_kind: GearTableKind, table_seed: u64) -> Self {
        let stream = GearStream::coarse();
        Self {
            min_size: stream.min_size,
            target_size: stream.target_size,
            max_size: stream.max_size,
            table_kind,
            table_seed,
        }
    }

    fn validate(self) -> Result<(), GearError> {
        if self.min_size == 0
            || self.min_size > self.target_size
            || self.target_size > self.max_size
        {
            return Err(GearError::InvalidChunkBounds {
                min: self.min_size,
                target: self.target_size,
                max: self.max_size,
            });
        }
        Ok(())
    }

    #[must_use]
    pub fn table(self) -> GearTable {
        GearTable::new(self.table_kind, self.table_seed)
    }
}

/// Errors from invalid chunk policies or coordinates.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GearError {
    InvalidChunkBounds { min: u32, target: u32, max: u32 },
    CoordinateOverflow,
}

impl std::fmt::Display for GearError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidChunkBounds { min, target, max } => {
                write!(
                    formatter,
                    "invalid Gear bounds min={min}, target={target}, max={max}"
                )
            }
            Self::CoordinateOverflow => formatter.write_str("Gear coordinate overflow"),
        }
    }
}

impl std::error::Error for GearError {}

/// Which strand views are emitted by a stream.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum FragmentationMode {
    Forward,
    BothStrands,
    /// Use the union of forward and reverse-complement cut coordinates and a
    /// canonical orientation for each resulting fragment.  The union makes
    /// boundaries invariant under reverse complementation.
    StrandSymmetric,
}

/// Orientation of the sequence represented by a fragment.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum FragmentOrientation {
    Forward,
    Reverse,
}

/// One exact-fragment lookup record.  `start` is always a forward-coordinate
/// interval in the source contig, even for a reverse-oriented fragment.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct ExactFragment {
    pub contig_id: u32,
    pub start: u64,
    pub length: u32,
    pub orientation: FragmentOrientation,
    pub digest: u64,
}

impl ExactFragment {
    #[must_use]
    pub fn end(self) -> Option<u64> {
        self.start.checked_add(u64::from(self.length))
    }
}

/// A digest-only candidate index.  Callers must use `verify_exact_fragment`
/// before treating a candidate as evidence.
#[derive(Clone, Debug, Default)]
pub struct ExactFragmentIndex {
    by_digest: BTreeMap<u64, Vec<ExactFragment>>,
}

impl ExactFragmentIndex {
    pub fn insert(&mut self, fragment: ExactFragment) {
        self.by_digest
            .entry(fragment.digest)
            .or_default()
            .push(fragment);
    }

    #[must_use]
    pub fn candidates(&self, digest: u64) -> &[ExactFragment] {
        self.by_digest.get(&digest).map_or(&[], Vec::as_slice)
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.by_digest.values().map(Vec::len).sum()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.by_digest.is_empty()
    }
}

/// Split one sequence in its own orientation.  The returned tuples contain a
/// forward-coordinate start only after `fragment_sequence` maps reverse
/// records back to the source sequence.
fn split_oriented(
    sequence: &[u8],
    config: GearConfig,
) -> Result<Vec<(usize, usize, u64)>, GearError> {
    config.validate()?;
    if sequence.is_empty() {
        return Ok(Vec::new());
    }
    let table = config.table();
    let min = usize::try_from(config.min_size).map_err(|_| GearError::CoordinateOverflow)?;
    let target = usize::try_from(config.target_size).map_err(|_| GearError::CoordinateOverflow)?;
    let max = usize::try_from(config.max_size).map_err(|_| GearError::CoordinateOverflow)?;
    // Start the geometric boundary process at `min`, and calibrate the
    // expected post-minimum waiting length to `target - min`.  Waiting until
    // `target` before testing a 1/target mask biases ordinary DNA toward
    // `target + min` (and, for repetitive tables, toward `max`).
    let post_minimum = target
        .checked_sub(min)
        .ok_or(GearError::CoordinateOverflow)?
        .max(1);
    let cut_threshold =
        u64::MAX / u64::try_from(post_minimum).map_err(|_| GearError::CoordinateOverflow)?;
    let mut result = Vec::new();
    let mut start = 0usize;
    let mut hasher = GearHasher::default();
    for (position, byte) in sequence.iter().copied().enumerate() {
        let hash = hasher.update(&table, normalized_byte(byte));
        let length = position + 1 - start;
        let boundary = length >= max || (length >= min && hash <= cut_threshold);
        if boundary {
            let end = position + 1;
            result.push((start, end - start, digest_fragment(&sequence[start..end])));
            start = end;
            hasher.reset();
        }
    }
    if start < sequence.len() {
        result.push((
            start,
            sequence.len() - start,
            digest_fragment(&sequence[start..sequence.len()]),
        ));
    }
    Ok(result)
}

/// Reverse-complement a byte sequence while normalizing lower-case/U bases.
pub(crate) fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|byte| match normalized_byte(*byte) {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

fn mapped_reverse_start(
    sequence_length: usize,
    start: usize,
    length: usize,
) -> Result<usize, GearError> {
    let end = start
        .checked_add(length)
        .ok_or(GearError::CoordinateOverflow)?;
    sequence_length
        .checked_sub(end)
        .ok_or(GearError::CoordinateOverflow)
}

/// Select the canonical orientation from the exact normalized bytes.  The
/// digest is deliberately not consulted: a digest tie or collision must not
/// choose an orientation inconsistently.
fn select_canonical_fragment<'a>(
    forward: &'a [u8],
    reverse: &'a [u8],
) -> (&'a [u8], FragmentOrientation) {
    if reverse < forward {
        (reverse, FragmentOrientation::Reverse)
    } else {
        (forward, FragmentOrientation::Forward)
    }
}

fn records_from_boundaries(
    sequence: &[u8],
    contig_id: u32,
    boundaries: &[usize],
    orientation: FragmentOrientation,
    canonical_digest: bool,
) -> Result<Vec<ExactFragment>, GearError> {
    let mut records = Vec::new();
    for window in boundaries.windows(2) {
        let start = window[0];
        let end = window[1];
        if end <= start {
            continue;
        }
        let fragment = sequence
            .get(start..end)
            .ok_or(GearError::CoordinateOverflow)?;
        let normalized_forward = fragment
            .iter()
            .copied()
            .map(normalized_byte)
            .collect::<Vec<_>>();
        let reverse = reverse_complement(&normalized_forward);
        let (selected_bytes, canonical_orientation) = if canonical_digest {
            select_canonical_fragment(&normalized_forward, &reverse)
        } else {
            (normalized_forward.as_slice(), orientation)
        };
        records.push(ExactFragment {
            contig_id,
            start: u64::try_from(start).map_err(|_| GearError::CoordinateOverflow)?,
            length: u32::try_from(end - start).map_err(|_| GearError::CoordinateOverflow)?,
            orientation: canonical_orientation,
            digest: digest_fragment(selected_bytes),
        });
    }
    Ok(records)
}

fn boundaries_from_parts(
    parts: &[(usize, usize, u64)],
    sequence_length: usize,
) -> Result<Vec<usize>, GearError> {
    let mut boundaries = vec![0usize];
    for (start, length, _) in parts {
        boundaries.push(
            start
                .checked_add(*length)
                .ok_or(GearError::CoordinateOverflow)?,
        );
    }
    boundaries.push(sequence_length);
    boundaries.sort_unstable();
    boundaries.dedup();
    Ok(boundaries)
}

/// Return sorted half-open boundaries for one Gear policy.
pub fn fragment_boundaries(sequence: &[u8], config: GearConfig) -> Result<Vec<u64>, GearError> {
    let parts = split_oriented(sequence, config)?;
    let boundaries = boundaries_from_parts(&parts, sequence.len())?;
    boundaries
        .into_iter()
        .map(|boundary| u64::try_from(boundary).map_err(|_| GearError::CoordinateOverflow))
        .collect()
}

fn symmetric_boundaries(sequence: &[u8], config: GearConfig) -> Result<Vec<usize>, GearError> {
    let forward = split_oriented(sequence, config)?;
    let reverse = reverse_complement(sequence);
    let reverse_parts = split_oriented(&reverse, config)?;
    let length = sequence.len();
    let mut cuts = vec![0usize, length];
    for (start, part_length, _) in forward {
        cuts.push(
            start
                .checked_add(part_length)
                .ok_or(GearError::CoordinateOverflow)?,
        );
        cuts.push(start);
    }
    for (start, part_length, _) in reverse_parts {
        let mapped_start = mapped_reverse_start(length, start, part_length)?;
        let mapped_end = mapped_start
            .checked_add(part_length)
            .ok_or(GearError::CoordinateOverflow)?;
        cuts.push(mapped_start);
        cuts.push(mapped_end);
    }
    cuts.sort_unstable();
    cuts.dedup();
    Ok(cuts)
}

/// Fragment a contig and retain only coordinates/digests, never a copied
/// complete sequence.  Reverse records map to forward source coordinates.
pub fn fragment_sequence(
    sequence: &[u8],
    contig_id: u32,
    config: GearConfig,
    mode: FragmentationMode,
) -> Result<Vec<ExactFragment>, GearError> {
    config.validate()?;
    if sequence.is_empty() {
        return Ok(Vec::new());
    }
    let mut records = match mode {
        FragmentationMode::Forward => {
            let parts = split_oriented(sequence, config)?;
            let boundaries = boundaries_from_parts(&parts, sequence.len())?;
            records_from_boundaries(
                sequence,
                contig_id,
                &boundaries,
                FragmentOrientation::Forward,
                false,
            )?
        }
        FragmentationMode::BothStrands => {
            let forward_parts = split_oriented(sequence, config)?;
            let forward_boundaries = boundaries_from_parts(&forward_parts, sequence.len())?;
            let mut result = records_from_boundaries(
                sequence,
                contig_id,
                &forward_boundaries,
                FragmentOrientation::Forward,
                false,
            )?;
            let reverse = reverse_complement(sequence);
            let reverse_parts = split_oriented(&reverse, config)?;
            for (start, length, _) in reverse_parts {
                let forward_start = mapped_reverse_start(sequence.len(), start, length)?;
                let reverse_end = start
                    .checked_add(length)
                    .ok_or(GearError::CoordinateOverflow)?;
                let fragment = &reverse[start..reverse_end];
                result.push(ExactFragment {
                    contig_id,
                    start: u64::try_from(forward_start)
                        .map_err(|_| GearError::CoordinateOverflow)?,
                    length: u32::try_from(length).map_err(|_| GearError::CoordinateOverflow)?,
                    orientation: FragmentOrientation::Reverse,
                    digest: digest_fragment(fragment),
                });
            }
            result
        }
        FragmentationMode::StrandSymmetric => {
            let boundaries = symmetric_boundaries(sequence, config)?;
            records_from_boundaries(
                sequence,
                contig_id,
                &boundaries,
                FragmentOrientation::Forward,
                true,
            )?
        }
    };
    records.sort_unstable_by_key(|fragment| {
        (
            fragment.start,
            fragment.orientation,
            fragment.length,
            fragment.digest,
        )
    });
    Ok(records)
}

/// Decode one fragment in its recorded orientation.  This is deliberately a
/// range operation; callers choose the source sequence/contig and no index
/// function copies an entire archive.
pub fn fragment_bytes(sequence: &[u8], fragment: ExactFragment) -> Option<Vec<u8>> {
    let start = usize::try_from(fragment.start).ok()?;
    let length = usize::try_from(fragment.length).ok()?;
    let end = start.checked_add(length)?;
    let source = sequence.get(start..end)?;
    Some(match fragment.orientation {
        FragmentOrientation::Forward => source.iter().map(|byte| normalized_byte(*byte)).collect(),
        FragmentOrientation::Reverse => reverse_complement(source),
    })
}

/// Stable digest for exact candidate lookup.  It is not evidence without
/// sequence verification and intentionally differs from `jamhash_u64_v1`.
#[must_use]
pub fn digest_fragment(sequence: &[u8]) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    for byte in sequence.iter().copied().map(normalized_byte) {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(0x1000_0000_01b3);
    }
    // Keep zero as a valid digest value: the exact verifier, never the digest,
    // is authoritative.
    hash
}

/// Verify digest and sequence bytes for a query/target fragment pair.
pub fn verify_exact_fragment(
    query_sequence: &[u8],
    query_fragment: ExactFragment,
    target_sequence: &[u8],
    target_fragment: ExactFragment,
) -> bool {
    if query_fragment.length != target_fragment.length
        || query_fragment.digest != target_fragment.digest
    {
        return false;
    }
    let Some(query) = fragment_bytes(query_sequence, query_fragment) else {
        return false;
    };
    let Some(target) = fragment_bytes(target_sequence, target_fragment) else {
        return false;
    };
    query == target
}

/// A verified contiguous fragment match in target-axis coordinates.  Reverse
/// records must be converted by the caller using `target_axis_start` so a
/// chain can remain monotonically increasing.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct VerifiedFragment {
    pub query_start: u64,
    pub target_axis_start: u64,
    pub length: u64,
    pub contig_id: u32,
    pub orientation: FragmentOrientation,
}

impl VerifiedFragment {
    pub fn query_end(self) -> Result<u64, GearError> {
        self.query_start
            .checked_add(self.length)
            .ok_or(GearError::CoordinateOverflow)
    }

    pub fn target_axis_end(self) -> Result<u64, GearError> {
        self.target_axis_start
            .checked_add(self.length)
            .ok_or(GearError::CoordinateOverflow)
    }
}

/// A merged exact run.  `query_span`/`target_span` include only matching
/// adjacent gaps; `verified_bases` is the sum of exact fragment lengths.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct ExactRun {
    pub query_start: u64,
    pub query_end: u64,
    pub target_start: u64,
    pub target_end: u64,
    pub contig_id: u32,
    pub orientation: FragmentOrientation,
    pub fragment_count: u32,
    pub verified_bases: u64,
}

impl ExactRun {
    #[must_use]
    pub fn is_contiguous(self) -> bool {
        let Some(query_span) = self.query_end.checked_sub(self.query_start) else {
            return false;
        };
        let Some(target_span) = self.target_end.checked_sub(self.target_start) else {
            return false;
        };
        self.verified_bases == query_span && self.verified_bases == target_span
    }

    /// Return a direct equal CIGAR only when every spanned base was verified.
    #[must_use]
    pub fn direct_cigar(self) -> Option<String> {
        self.is_contiguous()
            .then(|| format!("{}=", self.verified_bases))
    }
}

/// Alias used by callers that want to emphasize query/target matching.
pub type FragmentMatch = VerifiedFragment;

/// Merge deterministic compatible Gear matches.  Equal query/target gap
/// lengths are retained as a span but are not called direct `=` evidence until
/// the caller verifies those gap bases separately.
pub fn merge_exact_runs(matches: &[VerifiedFragment]) -> Result<Vec<ExactRun>, GearError> {
    let mut ordered = matches.to_vec();
    ordered.sort_unstable_by_key(|item| {
        (
            item.contig_id,
            item.orientation,
            item.query_start,
            item.target_axis_start,
            item.length,
        )
    });
    let mut runs: Vec<ExactRun> = Vec::new();
    for item in ordered {
        let item_query_end = item.query_end()?;
        let item_target_end = item.target_axis_end()?;
        let mut merged = false;
        if let Some(previous) = runs.last_mut()
            && previous.contig_id == item.contig_id
            && previous.orientation == item.orientation
            && item.query_start >= previous.query_end
            && item.target_axis_start >= previous.target_end
        {
            let query_gap = item
                .query_start
                .checked_sub(previous.query_end)
                .ok_or(GearError::CoordinateOverflow)?;
            let target_gap = item
                .target_axis_start
                .checked_sub(previous.target_end)
                .ok_or(GearError::CoordinateOverflow)?;
            if query_gap == target_gap {
                previous.query_end = item_query_end;
                previous.target_end = item_target_end;
                previous.fragment_count = previous
                    .fragment_count
                    .checked_add(1)
                    .ok_or(GearError::CoordinateOverflow)?;
                previous.verified_bases = previous
                    .verified_bases
                    .checked_add(item.length)
                    .ok_or(GearError::CoordinateOverflow)?;
                merged = true;
            }
        }
        if !merged {
            runs.push(ExactRun {
                query_start: item.query_start,
                query_end: item_query_end,
                target_start: item.target_axis_start,
                target_end: item_target_end,
                contig_id: item.contig_id,
                orientation: item.orientation,
                fragment_count: 1,
                verified_bases: item.length,
            });
        }
    }
    Ok(runs)
}

/// Summary used by boundary-distribution experiments.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FragmentLengthStats {
    pub count: u64,
    pub min: u64,
    pub max: u64,
    pub mean: f64,
    pub p50: u64,
    pub p95: u64,
}

/// Explicit accounting for the short fragments introduced by the
/// strand-symmetric union of two independent boundary streams.  The union is
/// experimental: it may create fragments below the configured minimum, so a
/// caller must inspect this report rather than treating the minimum as a hard
/// invariant for this mode.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StrandSymmetricDistribution {
    pub stats: FragmentLengthStats,
    pub configured_min: u64,
    pub below_min_count: u64,
    pub below_min_bases: u64,
    pub below_min_fraction: f64,
}

/// Quantify the short-fragment behavior of the strand-symmetric boundary
/// union.  `below_min_count` is expected to be non-zero for many ordinary
/// sequences because the union retains cuts from both strand streams.
pub fn strand_symmetric_distribution(
    sequence: &[u8],
    config: GearConfig,
) -> Result<StrandSymmetricDistribution, GearError> {
    let fragments = fragment_sequence(sequence, 0, config, FragmentationMode::StrandSymmetric)?;
    let stats = boundary_distribution(&fragments);
    let configured_min = u64::from(config.min_size);
    let mut below_min_count = 0u64;
    let mut below_min_bases = 0u64;
    for fragment in fragments {
        let length = u64::from(fragment.length);
        if length < configured_min {
            below_min_count = below_min_count
                .checked_add(1)
                .ok_or(GearError::CoordinateOverflow)?;
            below_min_bases = below_min_bases
                .checked_add(length)
                .ok_or(GearError::CoordinateOverflow)?;
        }
    }
    let below_min_fraction = if stats.count == 0 {
        0.0
    } else {
        below_min_count as f64 / stats.count as f64
    };
    Ok(StrandSymmetricDistribution {
        stats,
        configured_min,
        below_min_count,
        below_min_bases,
        below_min_fraction,
    })
}

/// Compute bounded fragment-length distribution statistics.
#[must_use]
pub fn boundary_distribution(fragments: &[ExactFragment]) -> FragmentLengthStats {
    let mut lengths = fragments
        .iter()
        .map(|fragment| u64::from(fragment.length))
        .collect::<Vec<_>>();
    if lengths.is_empty() {
        return FragmentLengthStats {
            count: 0,
            min: 0,
            max: 0,
            mean: 0.0,
            p50: 0,
            p95: 0,
        };
    }
    lengths.sort_unstable();
    let count = lengths.len() as u64;
    let sum = lengths.iter().sum::<u64>();
    let percentile = |numerator: usize, denominator: usize| {
        let index = (lengths.len() * numerator)
            .div_ceil(denominator)
            .saturating_sub(1);
        lengths[index.min(lengths.len() - 1)]
    };
    FragmentLengthStats {
        count,
        min: lengths[0],
        max: lengths[lengths.len() - 1],
        mean: sum as f64 / count as f64,
        p50: percentile(50, 100),
        p95: percentile(95, 100),
    }
}

#[cfg(test)]
mod tests {
    use super::{FragmentOrientation, reverse_complement, select_canonical_fragment};

    #[test]
    fn canonical_selection_uses_exact_bytes_on_a_digest_tie() {
        let forward = b"ACGT";
        let reverse = reverse_complement(forward);
        assert_eq!(forward, reverse.as_slice());
        let (selected, orientation) = select_canonical_fragment(forward, &reverse);
        assert_eq!(selected, forward);
        assert_eq!(orientation, FragmentOrientation::Forward);

        let forward = b"GTCA";
        let reverse = b"ACGT";
        let (selected, orientation) = select_canonical_fragment(forward, reverse);
        assert_eq!(selected, reverse);
        assert_eq!(orientation, FragmentOrientation::Reverse);
    }
}
