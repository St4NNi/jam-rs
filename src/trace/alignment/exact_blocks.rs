//! Packed exact-anchor extension for trace alignment.
//!
//! The affine-gap kernel is intentionally kept in [`super`].  This module is
//! the small, deterministic fast path which can run before that kernel: it
//! verifies an anchor against the packed canonical key, extends the exact
//! run in both directions, and merges adjacent runs without ever joining
//! different contigs or strands.

use crate::trace::anchors::Anchor;
use crate::trace::model::{BaseInterval, EditOperation, EditRun, Strand};
use needletail::Sequence;
use std::cmp::{max, min};
use thiserror::Error;

const BASES_PER_WORD: usize = 32;

/// Minimum exact seed span accepted by the trace workflow.
pub const MIN_EXACT_K: u8 = 21;

/// A two-bit representation used by exact extension.
///
/// Only A/C/G/T (and U as T) are packed.  The original bytes and an invalid
/// bit mask are retained so ambiguous bases can still use the alignment
/// kernel's case-insensitive, non-wildcard equality semantics.  A word-wise
/// comparison is used whenever both ranges are unambiguous and word aligned;
/// the scalar path is the correctness fallback for all other ranges.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PackedTwoBit {
    bytes: Vec<u8>,
    words: Vec<u64>,
    invalid_words: Vec<u64>,
}

impl PackedTwoBit {
    /// Pack an ASCII DNA sequence.  Lower-case bases and U are normalized for
    /// packing; every other byte is marked invalid and compared as a byte.
    #[must_use]
    pub fn from_bytes(bytes: &[u8]) -> Self {
        let word_count = bytes.len().div_ceil(BASES_PER_WORD);
        let mut words = vec![0u64; word_count];
        let mut invalid_words = vec![0u64; word_count];
        for (index, &base) in bytes.iter().enumerate() {
            let word = index / BASES_PER_WORD;
            let shift = (index % BASES_PER_WORD) * 2;
            match base_bits(base) {
                Some(bits) => words[word] |= u64::from(bits) << shift,
                None => invalid_words[word] |= 1u64 << (index % BASES_PER_WORD),
            }
        }
        Self {
            bytes: bytes.to_vec(),
            words,
            invalid_words,
        }
    }

    /// Build the packed reverse complement of a sequence.
    #[must_use]
    pub fn reverse_complement(bytes: &[u8]) -> Self {
        let reversed = bytes
            .iter()
            .rev()
            .copied()
            .map(complement)
            .collect::<Vec<_>>();
        Self::from_bytes(&reversed)
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.bytes.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.bytes.is_empty()
    }

    /// Return the first mismatch in two packed ranges, or the requested
    /// length when all compared bases are equal.
    ///
    /// Inputs are expected to contain the requested ranges.  If either range
    /// is shorter, comparison stops at the available length.  Callers that
    /// require an exact-length check should use [`compare_range`].
    #[must_use]
    pub fn first_mismatch(
        &self,
        self_start: usize,
        other: &Self,
        other_start: usize,
        length: usize,
    ) -> usize {
        let available = length
            .min(self.bytes.len().saturating_sub(self_start))
            .min(other.bytes.len().saturating_sub(other_start));
        let mut offset = 0usize;

        // Word comparisons are valid only when both ranges start on the same
        // packed word boundary.  This keeps the operation branch-light and
        // makes the differing two-bit pair directly recoverable.
        if self_start.is_multiple_of(BASES_PER_WORD) && other_start.is_multiple_of(BASES_PER_WORD) {
            while available.saturating_sub(offset) >= BASES_PER_WORD {
                let self_word = self.words[(self_start + offset) / BASES_PER_WORD];
                let other_word = other.words[(other_start + offset) / BASES_PER_WORD];
                let self_invalid = self.invalid_words[(self_start + offset) / BASES_PER_WORD];
                let other_invalid = other.invalid_words[(other_start + offset) / BASES_PER_WORD];
                if self_invalid == 0 && other_invalid == 0 {
                    let difference = self_word ^ other_word;
                    if difference != 0 {
                        return offset + (difference.trailing_zeros() as usize / 2);
                    }
                    offset += BASES_PER_WORD;
                    continue;
                }
                break;
            }
        }

        while offset < available {
            if !equal_base(
                self.bytes[self_start + offset],
                other.bytes[other_start + offset],
            ) {
                return offset;
            }
            offset += 1;
        }
        available
    }

    /// Compare two ranges for exact, case-insensitive equality.
    #[must_use]
    pub fn compare_range(
        &self,
        self_start: usize,
        other: &Self,
        other_start: usize,
        length: usize,
    ) -> bool {
        self_start
            .checked_add(length)
            .is_some_and(|end| end <= self.len())
            && other_start
                .checked_add(length)
                .is_some_and(|end| end <= other.len())
            && self.first_mismatch(self_start, other, other_start, length) == length
    }
}

/// Return the first mismatch between two sequences, or the common length.
#[must_use]
pub fn first_mismatch(left: &[u8], right: &[u8]) -> usize {
    let left = PackedTwoBit::from_bytes(left);
    let right = PackedTwoBit::from_bytes(right);
    left.first_mismatch(0, &right, 0, left.len().min(right.len()))
}

/// Return the first mismatch between a query and the reverse complement of a
/// forward-stored target, or the common length.
#[must_use]
pub fn first_mismatch_reverse_complement(query: &[u8], target_forward: &[u8]) -> usize {
    let query = PackedTwoBit::from_bytes(query);
    let target = PackedTwoBit::reverse_complement(target_forward);
    query.first_mismatch(0, &target, 0, query.len().min(target.len()))
}

/// Compare one forward-oriented target range exactly.
#[must_use]
pub fn compare_range(
    query: &[u8],
    target: &[u8],
    query_start: usize,
    target_start: usize,
    length: usize,
) -> bool {
    let query = PackedTwoBit::from_bytes(query);
    let target = PackedTwoBit::from_bytes(target);
    query.compare_range(query_start, &target, target_start, length)
}

/// Compare a query range with a forward-stored target range in the requested
/// relative orientation.  `target_start` always refers to the forward
/// target's k-mer start.
#[must_use]
pub fn compare_oriented_range(
    query: &[u8],
    target_forward: &[u8],
    query_start: usize,
    target_start: usize,
    length: usize,
    strand: Strand,
) -> bool {
    if strand == Strand::Forward {
        return compare_range(query, target_forward, query_start, target_start, length);
    }
    let target_end = match target_start.checked_add(length) {
        Some(end) if end <= target_forward.len() => end,
        _ => return false,
    };
    let oriented_start = target_forward.len() - target_end;
    let query = PackedTwoBit::from_bytes(query);
    let target = PackedTwoBit::reverse_complement(target_forward);
    query.compare_range(query_start, &target, oriented_start, length)
}

/// One verified exact block.  Both target intervals are retained: the
/// forward interval is the stable output coordinate, while the oriented one
/// makes chain/block merging unambiguous on the reverse strand.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ExactBlock {
    pub query_interval: BaseInterval,
    pub target_interval: BaseInterval,
    pub oriented_target_interval: BaseInterval,
    pub contig_id: u32,
    pub strand: Strand,
    pub anchor_count: u32,
    pub minimum_anchor_k: u8,
}

impl ExactBlock {
    #[must_use]
    pub fn len(&self) -> u64 {
        self.query_interval.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.query_interval.is_empty()
    }

    /// Direct edit representation used by the alignment/mosaic stage.  No
    /// affine-gap matrix is needed for a verified exact block.
    #[must_use]
    pub fn edit_runs(&self) -> Vec<EditRun> {
        if self.is_empty() {
            Vec::new()
        } else {
            vec![EditRun {
                operation: EditOperation::Equal,
                // A single EditRun is u32-sized in the frozen JSON model.
                // The trace sequence window is bounded well below this in
                // normal operation; callers handling a larger block should
                // use `edit_runs_checked` to avoid truncation.
                length: u32::try_from(self.len()).unwrap_or(u32::MAX),
            }]
        }
    }

    /// Checked form of [`Self::edit_runs`] for callers that need to reject a
    /// block too large for one serialized edit run.
    pub fn edit_runs_checked(&self) -> Result<Vec<EditRun>, ExactBlockError> {
        if self.is_empty() {
            return Ok(Vec::new());
        }
        Ok(vec![EditRun {
            operation: EditOperation::Equal,
            length: u32::try_from(self.len()).map_err(|_| ExactBlockError::RunTooLong)?,
        }])
    }

    /// CIGAR for the direct exact block (`=` only).
    #[must_use]
    pub fn cigar(&self) -> String {
        if self.is_empty() {
            String::new()
        } else {
            format!("{}=", self.len())
        }
    }
}

/// Build an exact block around one positional anchor.
pub fn extend_anchor(
    query: &[u8],
    target_forward: &[u8],
    anchor: Anchor,
) -> Result<ExactBlock, ExactBlockError> {
    let query = PackedTwoBit::from_bytes(query);
    let target_forward = PackedTwoBit::from_bytes(target_forward);
    let target_reverse = (anchor.strand == Strand::Reverse)
        .then(|| PackedTwoBit::reverse_complement(&target_forward.bytes));
    extend_anchor_packed(
        &query,
        &target_forward,
        target_reverse.as_ref(),
        target_forward.len(),
        anchor,
    )
}

/// Extend one anchor using caller-owned packed sequences.
///
/// Passing a precomputed reverse-complement packing avoids rebuilding large
/// target windows once per anchor. A reverse anchor requires
/// `target_reverse_complement`; forward anchors ignore that argument.
pub fn extend_anchor_packed(
    query: &PackedTwoBit,
    target_forward: &PackedTwoBit,
    target_reverse_complement: Option<&PackedTwoBit>,
    target_length: usize,
    anchor: Anchor,
) -> Result<ExactBlock, ExactBlockError> {
    if target_length != target_forward.len() {
        return Err(ExactBlockError::TargetLengthMismatch {
            declared: target_length,
            packed: target_forward.len(),
        });
    }
    let k = usize::from(anchor.k);
    if k == 0 {
        return Err(ExactBlockError::ZeroAnchorSpan);
    }
    if anchor.k < MIN_EXACT_K {
        return Err(ExactBlockError::InvalidAnchorK(anchor.k));
    }
    let query_start =
        usize::try_from(anchor.query_position).map_err(|_| ExactBlockError::CoordinateOverflow)?;
    let query_end = query_start
        .checked_add(k)
        .ok_or(ExactBlockError::CoordinateOverflow)?;
    if query_end > query.len() {
        return Err(ExactBlockError::QueryOutOfBounds {
            start: anchor.query_position,
            end: anchor
                .query_position
                .checked_add(u64::from(anchor.k))
                .ok_or(ExactBlockError::CoordinateOverflow)?,
            length: query.len() as u64,
        });
    }

    let target_start =
        usize::try_from(anchor.target_position).map_err(|_| ExactBlockError::CoordinateOverflow)?;
    let target_end = target_start
        .checked_add(k)
        .ok_or(ExactBlockError::CoordinateOverflow)?;
    if target_end > target_forward.len() {
        return Err(ExactBlockError::TargetOutOfBounds {
            start: anchor.target_position,
            end: anchor
                .target_position
                .checked_add(u64::from(anchor.k))
                .ok_or(ExactBlockError::CoordinateOverflow)?,
            length: target_forward.len() as u64,
        });
    }

    let oriented_target = match anchor.strand {
        Strand::Forward => target_forward,
        Strand::Reverse => {
            target_reverse_complement.ok_or(ExactBlockError::MissingReverseComplement)?
        }
    };
    let oriented_start = if anchor.strand == Strand::Forward {
        target_start
    } else {
        target_length
            .checked_sub(target_end)
            .ok_or(ExactBlockError::CoordinateOverflow)?
    };
    let oriented_end = oriented_start
        .checked_add(k)
        .ok_or(ExactBlockError::CoordinateOverflow)?;

    if !query.compare_range(query_start, oriented_target, oriented_start, k) {
        return Err(ExactBlockError::AnchorMismatch);
    }
    if let Some(expected) = canonical_packed(&query.bytes[query_start..query_end]) {
        if expected != anchor.canonical_kmer {
            return Err(ExactBlockError::PackedKeyMismatch {
                expected,
                observed: anchor.canonical_kmer,
            });
        }
        if crate::jamhash_u64_v1(expected) != anchor.hash {
            return Err(ExactBlockError::HashMismatch {
                expected: crate::jamhash_u64_v1(expected),
                observed: anchor.hash,
            });
        }
        if expected == 0 || anchor.hash == 0 {
            return Err(ExactBlockError::ZeroHash);
        }
    } else {
        return Err(ExactBlockError::AmbiguousAnchor);
    }

    let left_available = min(query_start, oriented_start);
    let left_mismatch = first_mismatch_reversed_ranges(
        query,
        query_start,
        oriented_target,
        oriented_start,
        left_available,
    );
    let left_match = left_available.min(left_mismatch);

    let right_available = query
        .len()
        .saturating_sub(query_end)
        .min(oriented_target.len().saturating_sub(oriented_end));
    let right_match =
        query.first_mismatch(query_end, oriented_target, oriented_end, right_available);

    let block_query_start = query_start
        .checked_sub(left_match)
        .ok_or(ExactBlockError::CoordinateOverflow)?;
    let block_query_end = query_end
        .checked_add(right_match)
        .ok_or(ExactBlockError::CoordinateOverflow)?;
    let block_oriented_start = oriented_start
        .checked_sub(left_match)
        .ok_or(ExactBlockError::CoordinateOverflow)?;
    let block_oriented_end = oriented_end
        .checked_add(right_match)
        .ok_or(ExactBlockError::CoordinateOverflow)?;
    let target_interval = if anchor.strand == Strand::Forward {
        BaseInterval::new(block_oriented_start as u64, block_oriented_end as u64)
    } else {
        BaseInterval::new(
            target_length
                .checked_sub(block_oriented_end)
                .ok_or(ExactBlockError::CoordinateOverflow)? as u64,
            target_length
                .checked_sub(block_oriented_start)
                .ok_or(ExactBlockError::CoordinateOverflow)? as u64,
        )
    }
    .map_err(|_| ExactBlockError::CoordinateOverflow)?;

    Ok(ExactBlock {
        query_interval: BaseInterval::new(block_query_start as u64, block_query_end as u64)
            .map_err(|_| ExactBlockError::CoordinateOverflow)?,
        target_interval,
        oriented_target_interval: BaseInterval::new(
            block_oriented_start as u64,
            block_oriented_end as u64,
        )
        .map_err(|_| ExactBlockError::CoordinateOverflow)?,
        contig_id: anchor.contig_id,
        strand: anchor.strand,
        anchor_count: 1,
        minimum_anchor_k: anchor.k,
    })
}

/// Extend and merge exact anchors.  Input order does not affect output order.
pub fn extend_and_merge_anchors(
    query: &[u8],
    target_forward: &[u8],
    anchors: &[Anchor],
) -> Result<Vec<ExactBlock>, ExactBlockError> {
    let query_packed = PackedTwoBit::from_bytes(query);
    let target_packed = PackedTwoBit::from_bytes(target_forward);
    let target_reverse = anchors
        .iter()
        .any(|anchor| anchor.strand == Strand::Reverse)
        .then(|| PackedTwoBit::reverse_complement(&target_packed.bytes));
    let mut blocks = anchors
        .iter()
        .copied()
        .map(|anchor| {
            extend_anchor_packed(
                &query_packed,
                &target_packed,
                target_reverse.as_ref(),
                target_forward.len(),
                anchor,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(merge_exact_blocks(&mut blocks))
}

/// Merge compatible verified blocks.  Blocks with different contig IDs or
/// strands remain separate, as do runs separated by an unverified gap.
#[must_use]
pub fn merge_exact_blocks(blocks: &mut [ExactBlock]) -> Vec<ExactBlock> {
    blocks.sort_by_key(|block| {
        (
            block.contig_id,
            strand_rank(block.strand),
            block.query_interval.start,
            block.oriented_target_interval.start,
        )
    });
    let mut merged = Vec::with_capacity(blocks.len());
    for block in blocks.iter().cloned() {
        let Some(previous) = merged.last_mut() else {
            merged.push(block);
            continue;
        };
        if !compatible_for_merge(previous, &block) {
            merged.push(block);
            continue;
        }
        previous.query_interval = BaseInterval::new(
            min(previous.query_interval.start, block.query_interval.start),
            max(previous.query_interval.end, block.query_interval.end),
        )
        .expect("verified intervals are ordered");
        previous.oriented_target_interval = BaseInterval::new(
            min(
                previous.oriented_target_interval.start,
                block.oriented_target_interval.start,
            ),
            max(
                previous.oriented_target_interval.end,
                block.oriented_target_interval.end,
            ),
        )
        .expect("verified intervals are ordered");
        previous.target_interval = if previous.strand == Strand::Forward {
            previous.oriented_target_interval
        } else {
            // The target length is not stored in a block.  For reverse blocks
            // the forward interval union is ordered independently and is
            // therefore sufficient for output/coverage.  Compatible reverse
            // blocks have the same oriented-to-forward diagonal, so this
            // union cannot bridge an unrelated interval.
            BaseInterval::new(
                min(previous.target_interval.start, block.target_interval.start),
                max(previous.target_interval.end, block.target_interval.end),
            )
            .expect("verified intervals are ordered")
        };
        previous.anchor_count = previous.anchor_count.saturating_add(block.anchor_count);
        previous.minimum_anchor_k = previous.minimum_anchor_k.min(block.minimum_anchor_k);
    }
    merged
}

fn compatible_for_merge(left: &ExactBlock, right: &ExactBlock) -> bool {
    left.contig_id == right.contig_id
        && left.strand == right.strand
        && left.query_interval.start <= right.query_interval.end
        && right.query_interval.start <= left.query_interval.end
        && left.oriented_target_interval.start <= right.oriented_target_interval.end
        && right.oriented_target_interval.start <= left.oriented_target_interval.end
        && left.query_interval.start as i128 - left.oriented_target_interval.start as i128
            == right.query_interval.start as i128 - right.oriented_target_interval.start as i128
}

fn first_mismatch_reversed_ranges(
    left: &PackedTwoBit,
    left_end: usize,
    right: &PackedTwoBit,
    right_end: usize,
    length: usize,
) -> usize {
    let left_reversed = PackedTwoBit::from_bytes(
        &left.bytes[..left_end]
            .iter()
            .rev()
            .copied()
            .collect::<Vec<_>>(),
    );
    let right_reversed = PackedTwoBit::from_bytes(
        &right.bytes[..right_end]
            .iter()
            .rev()
            .copied()
            .collect::<Vec<_>>(),
    );
    left_reversed.first_mismatch(0, &right_reversed, 0, length)
}

fn canonical_packed(sequence: &[u8]) -> Option<u64> {
    if sequence.is_empty() || sequence.len() > 31 {
        return None;
    }
    let normalized = sequence.normalize(false);
    normalized
        .bit_kmers(u8::try_from(sequence.len()).ok()?, true)
        .next()
        .map(|(_, kmer, _)| kmer.0)
}

fn base_bits(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

fn equal_base(left: u8, right: u8) -> bool {
    left.eq_ignore_ascii_case(&right)
}

fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        _ => b'N',
    }
}

const fn strand_rank(strand: Strand) -> u8 {
    match strand {
        Strand::Forward => 0,
        Strand::Reverse => 1,
    }
}

/// Errors raised before an exact block can be handed to alignment.
#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum ExactBlockError {
    #[error("exact anchor span must be greater than zero")]
    ZeroAnchorSpan,
    #[error("exact anchor span {0} is below the k=21 workflow minimum")]
    InvalidAnchorK(u8),
    #[error("exact anchor coordinate overflowed")]
    CoordinateOverflow,
    #[error("query anchor range [{start}, {end}) exceeds query length {length}")]
    QueryOutOfBounds { start: u64, end: u64, length: u64 },
    #[error("target anchor range [{start}, {end}) exceeds target length {length}")]
    TargetOutOfBounds { start: u64, end: u64, length: u64 },
    #[error("reverse anchor extension requires a packed reverse-complement target")]
    MissingReverseComplement,
    #[error("declared target length {declared} differs from packed target length {packed}")]
    TargetLengthMismatch { declared: usize, packed: usize },
    #[error("anchor sequence bytes are not an exact match")]
    AnchorMismatch,
    #[error("anchor k-mer contains an ambiguous base")]
    AmbiguousAnchor,
    #[error("anchor packed key differs: expected {expected}, observed {observed}")]
    PackedKeyMismatch { expected: u64, observed: u64 },
    #[error("anchor hash differs: expected {expected}, observed {observed}")]
    HashMismatch { expected: u64, observed: u64 },
    #[error("zero packed k-mer/hash is not an exact retained anchor")]
    ZeroHash,
    #[error("exact block exceeds the u32 serialized edit-run limit")]
    RunTooLong,
}
