//! The AVX2 recurrence of the rebased i16 planes: sixteen cells of one wave per block.

use super::AlignmentConfig;
use super::geometry::{BestCell, Cell, DELETION, INSERTION, MATCH, START, WaveRange};
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

/// AVX2 recurrence for rows `start..=end` of one wave in blocks of sixteen, with i16 offsets from
/// `rebase.base`, the floor offset as local zero and saturating adds. It equals the i32 kernel
/// while `rebase.violation` is clear. Returns the first row it did not fill. Callers must enable
/// AVX2, pass lowercase sequences, and prove the bounds that the debug assertions check.
#[allow(unused_assignments)]
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
        let check_low = rebase.check_low;
        let floor = _mm256_set1_epi16(zero_offset(base));
        let gap_open = _mm256_set1_epi16(gap_open_score as i16);
        let gap_extend = _mm256_set1_epi16(config.gap_extend_score as i16);
        // The rebased configuration bounds every scoring term to i16.
        let match_score = _mm256_set1_epi16(config.match_score as i16);
        let mismatch_score = _mm256_set1_epi16(config.mismatch_score as i16);
        let reverse = _mm_setr_epi8(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);

        // Row `r` of a plane is `plane.add(r)`; the pointers carry the range starts.
        let query_at = query.as_ptr().wrapping_sub(1);
        let contig_at = contig.as_ptr().wrapping_add(wave).wrapping_sub(16);
        let planes = |scores: ScoreWave<'_, Offset16>, shift: usize| {
            let first = scores.scores.as_ptr().cast::<i16>();
            [0, 1, 2].map(|state| {
                first
                    .wrapping_add(state * scores.stride + WAVE_PAD)
                    .wrapping_sub(scores.range.start + shift)
            })
        };
        let older_at = planes(older, 1);
        let previous_at = planes(previous, 0);
        let current_first = current.as_mut_ptr().cast::<i16>();
        let current_at = [0, 1, 2].map(|state| {
            current_first
                .wrapping_add(state * stride + WAVE_PAD)
                .wrapping_sub(current_range.start)
        });
        let trace_at = traceback
            .as_mut_ptr()
            .wrapping_add(wave_base)
            .wrapping_sub(current_range.start);

        let mut wave_max = _mm256_set1_epi16(i16::MIN);
        let mut wave_low = _mm256_set1_epi16(i16::MAX);
        let mut threshold = if BEST {
            threshold_vector(best.score, base)
        } else {
            _mm256_setzero_si256()
        };

        macro_rules! block {
            ($row:expr) => {{
                let row: usize = $row;
                let query_bytes = _mm_loadu_si128(query_at.wrapping_add(row).cast());
                let contig_bytes =
                    _mm_loadu_si128(contig_at.wrapping_sub(row).cast());
                let contig_bytes = _mm_shuffle_epi8(contig_bytes, reverse);
                // The caller folded both sequences to lowercase, so a plain compare is a case-blind one.
                let equal = _mm256_cvtepi8_epi16(_mm_cmpeq_epi8(query_bytes, contig_bytes));
                let substitution = _mm256_blendv_epi8(mismatch_score, match_score, equal);
                let load = |plane: *const i16, at: usize| {
                    _mm256_loadu_si256(plane.wrapping_add(at).cast::<__m256i>())
                };
                let older_m = load(older_at[MATCH as usize], row);
                let older_d = load(older_at[DELETION as usize], row);
                let older_i = load(older_at[INSERTION as usize], row);
                let left_m = load(previous_at[MATCH as usize], row);
                let left_d = load(previous_at[DELETION as usize], row);
                let left_i = load(previous_at[INSERTION as usize], row);
                let above_m = load(previous_at[MATCH as usize].wrapping_sub(1), row);
                let above_d = load(previous_at[DELETION as usize].wrapping_sub(1), row);
                let above_i = load(previous_at[INSERTION as usize].wrapping_sub(1), row);

                let (match_raw, match_scores, deletion_raw, deletion_scores);
                let (insertion_raw, insertion_scores, trace);
                if TRACEBACK {
                    let (diagonal, diagonal_state) = choose_state::<3>(older_m, older_d, older_i);
                    match_raw = _mm256_adds_epi16(diagonal, substitution);
                    let match_positive = _mm256_cmpgt_epi16(match_raw, floor);
                    match_scores = _mm256_max_epi16(match_raw, floor);
                    // The diagonal cell scored zero: a new alignment starts here.
                    let match_start = _mm256_cmpeq_epi16(diagonal, floor);
                    let match_previous = _mm256_blendv_epi8(
                        diagonal_state,
                        _mm256_set1_epi16(i16::from(START) << 3),
                        match_start,
                    );
                    let (raw, deletion_previous) = choose_state::<5>(
                        _mm256_adds_epi16(left_m, gap_open),
                        _mm256_adds_epi16(left_d, gap_extend),
                        _mm256_adds_epi16(left_i, gap_open),
                    );
                    deletion_raw = raw;
                    let deletion_positive = _mm256_cmpgt_epi16(deletion_raw, floor);
                    deletion_scores = _mm256_max_epi16(deletion_raw, floor);
                    let (raw, insertion_previous) = choose_state::<7>(
                        _mm256_adds_epi16(above_m, gap_open),
                        _mm256_adds_epi16(above_d, gap_open),
                        _mm256_adds_epi16(above_i, gap_extend),
                    );
                    insertion_raw = raw;
                    let insertion_positive = _mm256_cmpgt_epi16(insertion_raw, floor);
                    insertion_scores = _mm256_max_epi16(insertion_raw, floor);

                    // Each state contributes its positive bit and its shifted previous state,
                    // both only when the state is positive.
                    let match_bits = _mm256_and_si256(
                        match_positive,
                        _mm256_or_si256(match_previous, _mm256_set1_epi16(1)),
                    );
                    let deletion_bits = _mm256_and_si256(
                        deletion_positive,
                        _mm256_or_si256(deletion_previous, _mm256_set1_epi16(2)),
                    );
                    let insertion_bits = _mm256_and_si256(
                        insertion_positive,
                        _mm256_or_si256(insertion_previous, _mm256_set1_epi16(4)),
                    );
                    let bits = _mm256_or_si256(_mm256_or_si256(match_bits, deletion_bits), insertion_bits);
                    trace = bits;
                } else {
                    // Without a traceback only scores matter, and adding one constant keeps the
                    // maximum: the two candidates that share a gap penalty meet before the add.
                    let diagonal = _mm256_max_epi16(_mm256_max_epi16(older_m, older_d), older_i);
                    match_raw = _mm256_adds_epi16(diagonal, substitution);
                    match_scores = _mm256_max_epi16(match_raw, floor);
                    deletion_raw = _mm256_max_epi16(
                        _mm256_adds_epi16(_mm256_max_epi16(left_m, left_i), gap_open),
                        _mm256_adds_epi16(left_d, gap_extend),
                    );
                    deletion_scores = _mm256_max_epi16(deletion_raw, floor);
                    insertion_raw = _mm256_max_epi16(
                        _mm256_adds_epi16(_mm256_max_epi16(above_m, above_d), gap_open),
                        _mm256_adds_epi16(above_i, gap_extend),
                    );
                    insertion_scores = _mm256_max_epi16(insertion_raw, floor);
                    trace = _mm256_setzero_si256();
                }
                _mm256_storeu_si256(current_at[MATCH as usize].wrapping_add(row).cast(), match_scores);
                _mm256_storeu_si256(
                    current_at[DELETION as usize].wrapping_add(row).cast(),
                    deletion_scores,
                );
                _mm256_storeu_si256(
                    current_at[INSERTION as usize].wrapping_add(row).cast(),
                    insertion_scores,
                );
                if TRACEBACK {
                    _mm256_storeu_si256(trace_at.wrapping_add(row).cast(), trace);
                }
                let lane_best =
                    _mm256_max_epi16(match_scores, _mm256_max_epi16(deletion_scores, insertion_scores));

                // Wave statistics for the base policy: the largest best-state offset, and the
                // smallest winning candidate (which reveals unrepresentable values under a
                // sentinel-range base).
                wave_max = _mm256_max_epi16(wave_max, lane_best);
                if check_low {
                    wave_low = _mm256_min_epi16(
                        wave_low,
                        _mm256_min_epi16(match_raw, _mm256_min_epi16(deletion_raw, insertion_raw)),
                    );
                }

                // BestCell only changes for a positive score at or above the current best, so a
                // block whose lanes all lie at or below `threshold` cannot change it. The threshold
                // is that bound as an i16 offset, clamped, which can only loosen the test.
                if BEST && _mm256_movemask_epi8(_mm256_cmpgt_epi16(lane_best, threshold)) != 0 {
                    // Offsets can be negative, so flip the low 15 bits: the unsigned minimum of
                    // the flipped halves is the signed maximum.
                    let halves = _mm_max_epi16(
                        _mm256_castsi256_si128(lane_best),
                        _mm256_extracti128_si256::<1>(lane_best),
                    );
                    let flipped = _mm_xor_si128(halves, _mm_set1_epi16(0x7FFF));
                    let minimum = (_mm_cvtsi128_si32(_mm_minpos_epu16(flipped)) & 0xffff) as u16;
                    let block_offset = (minimum ^ 0x7FFF) as i16;
                    let block_score = i32::from(block_offset) + base;
                    if block_score > 0 && block_score >= best.score {
                        let mut matches = [0i16; 16];
                        let mut deletions = [0i16; 16];
                        let mut insertions = [0i16; 16];
                        _mm256_storeu_si256(matches.as_mut_ptr().cast(), match_scores);
                        _mm256_storeu_si256(deletions.as_mut_ptr().cast(), deletion_scores);
                        _mm256_storeu_si256(insertions.as_mut_ptr().cast(), insertion_scores);
                        let best_lanes = _mm256_movemask_epi8(_mm256_cmpeq_epi16(
                            lane_best,
                            _mm256_set1_epi16(block_offset),
                        )) as u32;
                        // Each lane sets two mask bits. As in `fill_wave_avx2`, the first lane
                        // with the block maximum is the lexicographically earliest candidate.
                        let block_lane = best_lanes.trailing_zeros() as usize / 2;
                        let query_index = row + block_lane;
                        let contig_index = wave - query_index;
                        best.consider(
                            query_index,
                            contig_index,
                            Cell {
                                scores: [
                                    Offset16(matches[block_lane]).read(base),
                                    Offset16(deletions[block_lane]).read(base),
                                    Offset16(insertions[block_lane]).read(base),
                                ],
                                previous: [0; 3],
                            },
                        );
                        threshold = threshold_vector(best.score, base);
                    }
                }
            }};
        }

        let last_block = end + 1 - 16;
        let mut row = start;
        while row <= last_block {
            block!(row);
            row += 16;
        }
        if row <= end {
            block!(last_block);
        }
        rebase.vector_max = rebase.vector_max.max(horizontal_max(wave_max));
        if check_low {
            rebase.vector_low = rebase.vector_low.min(horizontal_min(wave_low));
        }
        end + 1
    }
}
