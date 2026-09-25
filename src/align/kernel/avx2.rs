//! The AVX2 recurrence for eight cells of one wave with i32 scores.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::AlignmentConfig;
use super::geometry::{BestCell, DELETION, INSERTION, MATCH, START, WaveRange};
use super::waves::{ScoreWave, WAVE_PAD};
use std::arch::x86_64::*;

/// The inputs of one block of eight cells, rows `row..row + 8` of `wave`.
#[derive(Clone, Copy)]
pub(super) struct Block<'a> {
    pub(super) query: &'a [u8],
    pub(super) contig: &'a [u8],
    pub(super) config: AlignmentConfig,
    pub(super) gap_open_score: i32,
    pub(super) wave: usize,
    pub(super) row: usize,
    pub(super) current_range: WaveRange,
    pub(super) stride: usize,
    pub(super) older: ScoreWave<'a>,
    pub(super) previous: ScoreWave<'a>,
    pub(super) wave_base: usize,
}

#[target_feature(enable = "avx2")]
unsafe fn load_wave(scores: ScoreWave<'_>, state: usize, row: usize) -> __m256i {
    debug_assert!(scores.range.contains(row));
    debug_assert!(scores.range.contains(row + 7));
    let offset = state * scores.stride + WAVE_PAD + row - scores.range.start;
    // SAFETY: the caller proves rows row..row+7 are inside this wave and each state plane
    // has stride elements in a buffer of exactly three strides.
    unsafe { _mm256_loadu_si256(scores.scores.as_ptr().add(offset).cast()) }
}

/// The larger candidate and its state 0, 1 or 2: the second only when strictly larger than the
/// first, the third only when strictly larger than both.
#[target_feature(enable = "avx2")]
fn choose_avx2(m: __m256i, d: __m256i, i: __m256i) -> (__m256i, __m256i) {
    let zero = _mm256_setzero_si256();
    let one = _mm256_set1_epi32(1);
    let two = _mm256_set1_epi32(2);
    let d_better = _mm256_cmpgt_epi32(d, m);
    let mut score = _mm256_blendv_epi8(m, d, d_better);
    let mut state = _mm256_blendv_epi8(zero, one, d_better);
    let i_better = _mm256_cmpgt_epi32(i, score);
    score = _mm256_blendv_epi8(score, i, i_better);
    state = _mm256_blendv_epi8(state, two, i_better);
    (score, state)
}

#[target_feature(enable = "avx2")]
fn saturating_add_avx2(left: __m256i, right: __m256i) -> __m256i {
    let sum = _mm256_add_epi32(left, right);
    let overflow = _mm256_srai_epi32::<31>(_mm256_and_si256(
        _mm256_xor_si256(left, sum),
        _mm256_xor_si256(right, sum),
    ));
    let saturation = _mm256_xor_si256(_mm256_srai_epi32::<31>(left), _mm256_set1_epi32(i32::MAX));
    _mm256_blendv_epi8(sum, saturation, overflow)
}

/// AVX2 recurrence for eight cells. Inlined into its single AVX2-enabled caller, the vector
/// instantiation of `fill_waves`. Callers must enable AVX2 and prove that rows `row..row + 8`
/// lie in the current wave, that `row > 0`, and that `row + 7 < wave`.
#[inline(always)]
pub(super) unsafe fn fill_wave_avx2<const TRACEBACK: bool>(
    block: Block<'_>,
    current: &mut [i32],
    traceback: &mut [u16],
    best: &mut BestCell,
) {
    let Block {
        query,
        contig,
        config,
        wave,
        row,
        current_range,
        stride,
        wave_base,
        ..
    } = block;
    unsafe {
        debug_assert!(row > 0);
        debug_assert!(row + 7 < wave);
        debug_assert!(row + 7 <= query.len());
        debug_assert!(current_range.contains(row));
        debug_assert!(current_range.contains(row + 7));
        debug_assert!(current.len() >= 3 * stride);
        let first_contig = wave - row;
        debug_assert!(first_contig >= 8 && first_contig <= contig.len());

        let query_bytes = _mm_loadl_epi64(query.as_ptr().add(row - 1).cast());
        let contig_bytes = _mm_loadl_epi64(contig.as_ptr().add(first_contig - 8).cast());
        let reverse = _mm_setr_epi8(7, 6, 5, 4, 3, 2, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1);
        let contig_bytes = _mm_shuffle_epi8(contig_bytes, reverse);
        let equal = _mm256_cvtepi8_epi32(_mm_cmpeq_epi8(query_bytes, contig_bytes));
        let substitution = _mm256_blendv_epi8(
            _mm256_set1_epi32(config.mismatch_score),
            _mm256_set1_epi32(config.match_score),
            equal,
        );
        let states = block_states(&block, substitution);

        let wave_offset = row - current_range.start;
        let cell_offset = WAVE_PAD + wave_offset;
        let [matched, deleted, inserted] = states.scores;
        let cells = current.as_mut_ptr().add(cell_offset);
        _mm256_storeu_si256(cells.cast(), matched);
        _mm256_storeu_si256(cells.add(stride).cast(), deleted);
        _mm256_storeu_si256(cells.add(2 * stride).cast(), inserted);
        if TRACEBACK {
            debug_assert!(wave_base + wave_offset + 8 <= traceback.len());
            // SAFETY: the wave layout stores this wave's cells contiguously, and rows row..=row+7
            // lie inside current_range, so the eight cells are inside the laid-out traceback.
            let cells = traceback.as_mut_ptr().add(wave_base + wave_offset);
            _mm_storeu_si128(cells.cast(), pack_traceback(&states));
        }
        consider_block(&block, states.scores, best);
    }
}

/// The three states of eight cells: the scores kept after the zero floor, which states are
/// positive, and the state each came from.
struct BlockStates {
    scores: [__m256i; 3],
    positive: [__m256i; 3],
    previous: [__m256i; 3],
}

/// The states of the eight cells of `block` with the substitution scores of their bases.
#[inline(always)]
unsafe fn block_states(block: &Block<'_>, substitution: __m256i) -> BlockStates {
    let Block {
        config,
        gap_open_score,
        row,
        older,
        previous,
        ..
    } = *block;
    unsafe {
        let zero = _mm256_setzero_si256();
        let (open, extend) = (
            _mm256_set1_epi32(gap_open_score),
            _mm256_set1_epi32(config.gap_extend_score),
        );
        let (diagonal, diagonal_state) = choose_avx2(
            load_wave(older, MATCH as usize, row - 1),
            load_wave(older, DELETION as usize, row - 1),
            load_wave(older, INSERTION as usize, row - 1),
        );
        let match_raw = saturating_add_avx2(diagonal, substitution);
        let match_start = _mm256_cmpeq_epi32(match_raw, substitution);
        let match_previous = _mm256_blendv_epi8(
            diagonal_state,
            _mm256_set1_epi32(i32::from(START)),
            match_start,
        );
        // The cell to the left (one contig base back, same query row) opens or extends a
        // deletion.
        let (deletion_raw, deletion_previous) = choose_avx2(
            _mm256_add_epi32(load_wave(previous, MATCH as usize, row), open),
            _mm256_add_epi32(load_wave(previous, DELETION as usize, row), extend),
            _mm256_add_epi32(load_wave(previous, INSERTION as usize, row), open),
        );
        // The cell above (one query base back, same contig column) opens or extends an
        // insertion.
        let (insertion_raw, insertion_previous) = choose_avx2(
            _mm256_add_epi32(load_wave(previous, MATCH as usize, row - 1), open),
            _mm256_add_epi32(load_wave(previous, DELETION as usize, row - 1), open),
            _mm256_add_epi32(load_wave(previous, INSERTION as usize, row - 1), extend),
        );
        let match_positive = _mm256_cmpgt_epi32(match_raw, zero);
        let deletion_positive = _mm256_cmpgt_epi32(deletion_raw, zero);
        let insertion_positive = _mm256_cmpgt_epi32(insertion_raw, zero);
        BlockStates {
            scores: [
                _mm256_and_si256(match_raw, match_positive),
                _mm256_and_si256(deletion_raw, deletion_positive),
                _mm256_and_si256(insertion_raw, insertion_positive),
            ],
            positive: [match_positive, deletion_positive, insertion_positive],
            previous: [match_previous, deletion_previous, insertion_previous],
        }
    }
}

/// The traceback cells of eight states: bit `s` marks a positive state `s`, and bits `3 + 2s`
/// and `4 + 2s` hold the state it came from.
#[inline(always)]
unsafe fn pack_traceback(states: &BlockStates) -> __m128i {
    let [match_positive, deletion_positive, insertion_positive] = states.positive;
    let [match_previous, deletion_previous, insertion_previous] = states.previous;
    unsafe {
        let mut trace = _mm256_and_si256(match_positive, _mm256_set1_epi32(1));
        trace = _mm256_or_si256(
            trace,
            _mm256_slli_epi32::<3>(_mm256_and_si256(match_previous, match_positive)),
        );
        trace = _mm256_or_si256(
            trace,
            _mm256_and_si256(deletion_positive, _mm256_set1_epi32(2)),
        );
        trace = _mm256_or_si256(
            trace,
            _mm256_slli_epi32::<5>(_mm256_and_si256(deletion_previous, deletion_positive)),
        );
        trace = _mm256_or_si256(
            trace,
            _mm256_and_si256(insertion_positive, _mm256_set1_epi32(4)),
        );
        trace = _mm256_or_si256(
            trace,
            _mm256_slli_epi32::<7>(_mm256_and_si256(insertion_previous, insertion_positive)),
        );
        // Traces use nine bits, so unsigned saturation to 16 bits is exact. Packing places
        // lanes 0..3 and 4..7 in separate 128-bit halves; the permutation joins them in order.
        let zero = _mm256_setzero_si256();
        let packed = _mm256_permute4x64_epi64::<0xD8>(_mm256_packus_epi32(trace, zero));
        _mm256_castsi256_si128(packed)
    }
}

/// Offers the best lane of a block to the best cell. A block below the current best score
/// cannot change it, as the best cell only changes for a higher score or an equal, earlier cell.
#[inline(always)]
unsafe fn consider_block(block: &Block<'_>, scores: [__m256i; 3], best: &mut BestCell) {
    unsafe {
        let [matched, deleted, inserted] = scores;
        let (lane_best, _) = choose_avx2(matched, deleted, inserted);
        let mut block_best = _mm256_max_epi32(
            lane_best,
            _mm256_permute2x128_si256::<0x01>(lane_best, lane_best),
        );
        block_best = _mm256_max_epi32(block_best, _mm256_shuffle_epi32::<0x4e>(block_best));
        block_best = _mm256_max_epi32(block_best, _mm256_shuffle_epi32::<0xb1>(block_best));
        let block_score = _mm256_extract_epi32::<0>(block_best);
        if block_score > 0 && block_score >= best.score {
            let mut lanes = [[0i32; 8]; 3];
            _mm256_storeu_si256(lanes[0].as_mut_ptr().cast(), matched);
            _mm256_storeu_si256(lanes[1].as_mut_ptr().cast(), deleted);
            _mm256_storeu_si256(lanes[2].as_mut_ptr().cast(), inserted);
            let best_lanes = _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpeq_epi32(
                lane_best, block_best,
            )));
            let lane = best_lanes.trailing_zeros() as usize;
            let scores = [lanes[0][lane], lanes[1][lane], lanes[2][lane]];
            best.consider_lane(block.wave, block.row, lane, scores);
        }
    }
}
