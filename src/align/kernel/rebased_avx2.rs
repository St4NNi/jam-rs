//! The AVX2 recurrence of the rebased i16 planes: sixteen cells of one wave per block.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::AlignmentConfig;
use super::geometry::{BestCell, START, WaveRange};
use super::rebased::{Offset16, RebaseState, zero_offset};
use super::waves::{ScoreWave, WAVE_PAD, WaveScore};
use std::arch::x86_64::*;

/// The inputs of rows `span.0..=span.1` of one wave.
#[derive(Clone, Copy)]
pub(super) struct Span<'a> {
    pub(super) query: &'a [u8],
    pub(super) contig: &'a [u8],
    pub(super) config: AlignmentConfig,
    pub(super) gap_open_score: i32,
    pub(super) wave: usize,
    pub(super) span: (usize, usize),
    pub(super) current_range: WaveRange,
    pub(super) stride: usize,
    pub(super) older: ScoreWave<'a, Offset16>,
    pub(super) previous: ScoreWave<'a, Offset16>,
    pub(super) wave_base: usize,
}

/// The offset that a block compares against to skip the best-cell search.
#[inline]
fn threshold_offset(best_score: i32, base: i32) -> i16 {
    (best_score.max(1) - 1 - base).clamp(i32::from(i16::MIN), i32::from(i16::MAX)) as i16
}

#[target_feature(enable = "avx2")]
#[inline]
unsafe fn threshold_vector(best_score: i32, base: i32) -> __m256i {
    _mm256_set1_epi16(threshold_offset(best_score, base))
}

/// The largest of sixteen i16 lanes: flipping the low 15 bits makes the unsigned minimum of the
/// flipped lanes the signed maximum.
#[target_feature(enable = "avx2")]
#[inline]
fn horizontal_max(lanes: __m256i) -> i32 {
    let halves = _mm_max_epi16(
        _mm256_castsi256_si128(lanes),
        _mm256_extracti128_si256::<1>(lanes),
    );
    let flipped = _mm_xor_si128(halves, _mm_set1_epi16(0x7FFF));
    let minimum = (_mm_cvtsi128_si32(_mm_minpos_epu16(flipped)) & 0xffff) as u16;
    i32::from((minimum ^ 0x7FFF) as i16)
}

/// The smallest of sixteen i16 lanes: flipping the sign bit makes the unsigned minimum of the
/// flipped lanes the signed minimum.
#[target_feature(enable = "avx2")]
#[inline]
fn horizontal_min(lanes: __m256i) -> i32 {
    let halves = _mm_min_epi16(
        _mm256_castsi256_si128(lanes),
        _mm256_extracti128_si256::<1>(lanes),
    );
    let flipped = _mm_xor_si128(halves, _mm_set1_epi16(i16::MIN));
    let minimum = (_mm_cvtsi128_si32(_mm_minpos_epu16(flipped)) & 0xffff) as u16;
    i32::from((minimum ^ 0x8000) as i16)
}

/// Larger score and the state that won, as 0, 1 or 2 shifted left by `SHIFT`: the second candidate
/// only when strictly larger than the first, the third only when strictly larger than both.
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn choose_state<const SHIFT: i32>(m: __m256i, d: __m256i, i: __m256i) -> (__m256i, __m256i) {
    let d_better = _mm256_cmpgt_epi16(d, m);
    let first = _mm256_max_epi16(m, d);
    let i_better = _mm256_cmpgt_epi16(i, first);
    let state = _mm256_blendv_epi8(
        _mm256_and_si256(d_better, _mm256_set1_epi16(1 << SHIFT)),
        _mm256_set1_epi16(2 << SHIFT),
        i_better,
    );
    (_mm256_max_epi16(first, i), state)
}

/// Fills rows `start..=end` of a wave in blocks of sixteen i16 offsets, equal to the i32 kernel
/// while `rebase.violation` is clear, and returns the first row it did not fill. Callers must
/// enable AVX2, pass lowercase sequences, and prove the bounds the debug assertions check.
#[inline(always)]
pub(super) unsafe fn fill_span16_rebased<const TRACEBACK: bool, const BEST: bool>(
    blocks: Span<'_>,
    current: &mut [Offset16],
    traceback: &mut [u16],
    best: &mut BestCell,
    rebase: &mut RebaseState,
) -> usize {
    let Span {
        query,
        contig,
        config,
        gap_open_score,
        wave,
        span: (start, end),
        current_range,
        stride,
        older,
        previous,
        wave_base,
    } = blocks;
    unsafe {
        if end + 1 - start < 16 {
            return start;
        }
        debug_assert!(start > 0 && end < wave && end <= query.len());
        debug_assert!(current_range.contains(start) && current_range.contains(end));
        debug_assert!(older.range.contains(start - 1) && older.range.contains(end - 1));
        // Rows just outside the previous wave read its zeroed spare cells.
        debug_assert!(previous.range.start <= start && end <= previous.range.end + 1);
        debug_assert!(current.len() >= 3 * stride);
        debug_assert!(wave - end >= 1 && wave - start <= contig.len() + 15);
        let base = rebase.base;
        let planes = |scores: ScoreWave<'_, Offset16>, shift: usize| {
            let first = scores.scores.as_ptr().cast::<i16>();
            [0, 1, 2].map(|state| {
                first
                    .wrapping_add(state * scores.stride + WAVE_PAD)
                    .wrapping_sub(scores.range.start + shift)
            })
        };
        let current_first = current.as_mut_ptr().cast::<i16>();
        // The rebased configuration bounds every scoring term to i16.
        let span = SpanPlanes {
            query: query.as_ptr().wrapping_sub(1),
            contig: contig.as_ptr().wrapping_add(wave).wrapping_sub(16),
            older: planes(older, 1),
            previous: planes(previous, 0),
            current: [0, 1, 2].map(|state| {
                current_first
                    .wrapping_add(state * stride + WAVE_PAD)
                    .wrapping_sub(current_range.start)
            }),
            trace: traceback
                .as_mut_ptr()
                .wrapping_add(wave_base)
                .wrapping_sub(current_range.start),
            wave,
            base,
            floor: _mm256_set1_epi16(zero_offset(base)),
            gap_open: _mm256_set1_epi16(gap_open_score as i16),
            gap_extend: _mm256_set1_epi16(config.gap_extend_score as i16),
            match_score: _mm256_set1_epi16(config.match_score as i16),
            mismatch_score: _mm256_set1_epi16(config.mismatch_score as i16),
            reverse: _mm_setr_epi8(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0),
        };
        let mut stats = SpanStats {
            wave_max: _mm256_set1_epi16(i16::MIN),
            wave_low: _mm256_set1_epi16(i16::MAX),
            threshold: if BEST {
                threshold_vector(best.score, base)
            } else {
                _mm256_setzero_si256()
            },
            check_low: rebase.check_low,
        };
        let last_block = end + 1 - 16;
        let mut row = start;
        while row <= last_block {
            fill_block::<TRACEBACK, BEST>(&span, row, &mut stats, best);
            row += 16;
        }
        if row <= end {
            fill_block::<TRACEBACK, BEST>(&span, last_block, &mut stats, best);
        }
        rebase.vector_max = rebase.vector_max.max(horizontal_max(stats.wave_max));
        if stats.check_low {
            rebase.vector_low = rebase.vector_low.min(horizontal_min(stats.wave_low));
        }
        end + 1
    }
}

/// What every block of one span reads: the sequences and planes, where row `r` of a plane is
/// `plane.wrapping_add(r)`, and the scores as vectors.
struct SpanPlanes {
    query: *const u8,
    contig: *const u8,
    older: [*const i16; 3],
    previous: [*const i16; 3],
    current: [*mut i16; 3],
    trace: *mut u16,
    wave: usize,
    base: i32,
    floor: __m256i,
    gap_open: __m256i,
    gap_extend: __m256i,
    match_score: __m256i,
    mismatch_score: __m256i,
    reverse: __m128i,
}

/// What the blocks of a span collect for the base policy and the best cell: the largest
/// best-state offset, the smallest winning candidate, and the offset a block must pass.
struct SpanStats {
    wave_max: __m256i,
    wave_low: __m256i,
    threshold: __m256i,
    check_low: bool,
}

/// The three states of one block: the winning candidates, the offsets kept after the floor, and
/// the traceback cells.
struct BlockStates {
    raw: [__m256i; 3],
    scores: [__m256i; 3],
    trace: __m256i,
}

/// Fills the sixteen cells from `row` on. Callers must enable AVX2 and keep the block inside
/// the span that `fill_span16_rebased` checks.
#[inline(always)]
unsafe fn fill_block<const TRACEBACK: bool, const BEST: bool>(
    span: &SpanPlanes,
    row: usize,
    stats: &mut SpanStats,
    best: &mut BestCell,
) {
    unsafe {
        let query_bytes = _mm_loadu_si128(span.query.wrapping_add(row).cast());
        let contig_bytes = _mm_loadu_si128(span.contig.wrapping_sub(row).cast());
        let contig_bytes = _mm_shuffle_epi8(contig_bytes, span.reverse);
        // The caller folded both sequences to lowercase, so a plain compare ignores case.
        let equal = _mm256_cvtepi8_epi16(_mm_cmpeq_epi8(query_bytes, contig_bytes));
        let substitution = _mm256_blendv_epi8(span.mismatch_score, span.match_score, equal);
        let load = |plane: *const i16| _mm256_loadu_si256(plane.wrapping_add(row).cast());
        let older = span.older.map(load);
        let left = span.previous.map(load);
        let above = span.previous.map(|plane| load(plane.wrapping_sub(1)));
        let states = if TRACEBACK {
            traced_states(span, substitution, older, left, above)
        } else {
            plain_states(span, substitution, older, left, above)
        };
        for (plane, scores) in span.current.iter().zip(states.scores) {
            _mm256_storeu_si256(plane.wrapping_add(row).cast(), scores);
        }
        if TRACEBACK {
            _mm256_storeu_si256(span.trace.wrapping_add(row).cast(), states.trace);
        }
        let [matched, deleted, inserted] = states.scores;
        let lane_best = _mm256_max_epi16(matched, _mm256_max_epi16(deleted, inserted));
        stats.wave_max = _mm256_max_epi16(stats.wave_max, lane_best);
        // Under a sentinel-range base the smallest winning candidate reveals values that the
        // offsets cannot hold.
        if stats.check_low {
            let [matched, deleted, inserted] = states.raw;
            let low = _mm256_min_epi16(matched, _mm256_min_epi16(deleted, inserted));
            stats.wave_low = _mm256_min_epi16(stats.wave_low, low);
        }
        // A block whose lanes all lie at or below the threshold, the best score as a clamped
        // offset, cannot change the best cell.
        if BEST && _mm256_movemask_epi8(_mm256_cmpgt_epi16(lane_best, stats.threshold)) != 0 {
            consider_block(span, row, lane_best, states.scores, stats, best);
        }
    }
}

/// The states of a block with the state each came from, packed into traceback cells.
#[inline(always)]
unsafe fn traced_states(
    span: &SpanPlanes,
    substitution: __m256i,
    [older_m, older_d, older_i]: [__m256i; 3],
    [left_m, left_d, left_i]: [__m256i; 3],
    [above_m, above_d, above_i]: [__m256i; 3],
) -> BlockStates {
    unsafe {
        let (floor, open, extend) = (span.floor, span.gap_open, span.gap_extend);
        let (diagonal, diagonal_state) = choose_state::<3>(older_m, older_d, older_i);
        let matched = _mm256_adds_epi16(diagonal, substitution);
        // The diagonal cell scored zero: a new alignment starts here.
        let match_previous = _mm256_blendv_epi8(
            diagonal_state,
            _mm256_set1_epi16(i16::from(START) << 3),
            _mm256_cmpeq_epi16(diagonal, floor),
        );
        let (deleted, deletion_previous) = choose_state::<5>(
            _mm256_adds_epi16(left_m, open),
            _mm256_adds_epi16(left_d, extend),
            _mm256_adds_epi16(left_i, open),
        );
        let (inserted, insertion_previous) = choose_state::<7>(
            _mm256_adds_epi16(above_m, open),
            _mm256_adds_epi16(above_d, open),
            _mm256_adds_epi16(above_i, extend),
        );
        // A positive state sets its own bit and the shifted state it came from.
        let bits = |raw: __m256i, previous: __m256i, bit: i16| {
            _mm256_and_si256(
                _mm256_cmpgt_epi16(raw, floor),
                _mm256_or_si256(previous, _mm256_set1_epi16(bit)),
            )
        };
        let trace = _mm256_or_si256(
            _mm256_or_si256(
                bits(matched, match_previous, 1),
                bits(deleted, deletion_previous, 2),
            ),
            bits(inserted, insertion_previous, 4),
        );
        let raw = [matched, deleted, inserted];
        BlockStates {
            raw,
            scores: raw.map(|state| _mm256_max_epi16(state, floor)),
            trace,
        }
    }
}

/// The states of a block without traceback. Adding one constant keeps the maximum, so the two
/// candidates that share a gap score meet before the add.
#[inline(always)]
unsafe fn plain_states(
    span: &SpanPlanes,
    substitution: __m256i,
    [older_m, older_d, older_i]: [__m256i; 3],
    [left_m, left_d, left_i]: [__m256i; 3],
    [above_m, above_d, above_i]: [__m256i; 3],
) -> BlockStates {
    unsafe {
        let (floor, open, extend) = (span.floor, span.gap_open, span.gap_extend);
        let diagonal = _mm256_max_epi16(_mm256_max_epi16(older_m, older_d), older_i);
        let matched = _mm256_adds_epi16(diagonal, substitution);
        let deleted = _mm256_max_epi16(
            _mm256_adds_epi16(_mm256_max_epi16(left_m, left_i), open),
            _mm256_adds_epi16(left_d, extend),
        );
        let inserted = _mm256_max_epi16(
            _mm256_adds_epi16(_mm256_max_epi16(above_m, above_d), open),
            _mm256_adds_epi16(above_i, extend),
        );
        let raw = [matched, deleted, inserted];
        BlockStates {
            raw,
            scores: raw.map(|state| _mm256_max_epi16(state, floor)),
            trace: _mm256_setzero_si256(),
        }
    }
}

/// Offers the best lane of a block that passed the threshold to the best cell, and raises the
/// threshold when the best cell changes.
#[inline(always)]
unsafe fn consider_block(
    span: &SpanPlanes,
    row: usize,
    lane_best: __m256i,
    scores: [__m256i; 3],
    stats: &mut SpanStats,
    best: &mut BestCell,
) {
    unsafe {
        let block_offset = horizontal_max(lane_best) as i16;
        let block_score = i32::from(block_offset) + span.base;
        if block_score <= 0 || block_score < best.score {
            return;
        }
        let mut lanes = [[0i16; 16]; 3];
        for (lane, state) in lanes.iter_mut().zip(scores) {
            _mm256_storeu_si256(lane.as_mut_ptr().cast(), state);
        }
        let best_lanes = _mm256_movemask_epi8(_mm256_cmpeq_epi16(
            lane_best,
            _mm256_set1_epi16(block_offset),
        )) as u32;
        // Each lane sets two mask bits.
        let lane = best_lanes.trailing_zeros() as usize / 2;
        let lane_scores = lanes.map(|state| Offset16(state[lane]).read(span.base));
        best.consider_lane(span.wave, row, lane, lane_scores);
        stats.threshold = threshold_vector(best.score, span.base);
    }
}
