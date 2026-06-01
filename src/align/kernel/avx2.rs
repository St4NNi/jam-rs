//! The AVX2 recurrence for eight cells of one wave with i32 scores. It serves passes the rebased
//! i16 kernel cannot represent.

use super::AlignmentConfig;
use super::geometry::{BestCell, Cell, DELETION, INSERTION, MATCH, START, WaveRange};
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

/// Lanes where two bytes are equal, ignoring the case of letters.
#[target_feature(enable = "avx2")]
fn equal_ascii_case(left: __m128i, right: __m128i) -> __m128i {
    let case_bit = _mm_set1_epi8(0x20);
    let lower = _mm_or_si128(left, case_bit);
    let letter = _mm_and_si128(
        _mm_cmpgt_epi8(lower, _mm_set1_epi8((b'a' - 1) as i8)),
        _mm_cmpgt_epi8(_mm_set1_epi8((b'z' + 1) as i8), lower),
    );
    let case_pair = _mm_cmpeq_epi8(_mm_xor_si128(left, right), case_bit);
    _mm_or_si128(
        _mm_cmpeq_epi8(left, right),
        _mm_and_si128(case_pair, letter),
    )
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
        gap_open_score,
        wave,
        row,
        current_range,
        stride,
        older,
        previous,
        wave_base,
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
        let equal = equal_ascii_case(query_bytes, contig_bytes);
        let equal = _mm256_cvtepi8_epi32(equal);
        let substitution = _mm256_blendv_epi8(
            _mm256_set1_epi32(config.mismatch_score),
            _mm256_set1_epi32(config.match_score),
            equal,
        );

        let (diagonal, diagonal_state) = choose_avx2(
            load_wave(older, MATCH as usize, row - 1),
            load_wave(older, DELETION as usize, row - 1),
            load_wave(older, INSERTION as usize, row - 1),
        );
        let match_raw = saturating_add_avx2(diagonal, substitution);
        let zero = _mm256_setzero_si256();
        let match_positive = _mm256_cmpgt_epi32(match_raw, zero);
        let match_scores = _mm256_and_si256(match_raw, match_positive);
        let match_start = _mm256_cmpeq_epi32(match_raw, substitution);
        let match_previous = _mm256_blendv_epi8(
            diagonal_state,
            _mm256_set1_epi32(i32::from(START)),
            match_start,
        );

        // The cell to the left (one contig base back, same query row) opens or extends a
        // deletion.
        let left_m = _mm256_add_epi32(
            load_wave(previous, MATCH as usize, row),
            _mm256_set1_epi32(gap_open_score),
        );
        let left_d = _mm256_add_epi32(
            load_wave(previous, DELETION as usize, row),
            _mm256_set1_epi32(config.gap_extend_score),
        );
        let left_i = _mm256_add_epi32(
            load_wave(previous, INSERTION as usize, row),
            _mm256_set1_epi32(gap_open_score),
        );
        let (deletion_raw, deletion_previous) = choose_avx2(left_m, left_d, left_i);
        let deletion_positive = _mm256_cmpgt_epi32(deletion_raw, zero);
        let deletion_scores = _mm256_and_si256(deletion_raw, deletion_positive);

        // The cell above (one query base back, same contig column) opens or extends an
        // insertion.
        let above_m = _mm256_add_epi32(
            load_wave(previous, MATCH as usize, row - 1),
            _mm256_set1_epi32(gap_open_score),
        );
        let above_d = _mm256_add_epi32(
            load_wave(previous, DELETION as usize, row - 1),
            _mm256_set1_epi32(gap_open_score),
        );
        let above_i = _mm256_add_epi32(
            load_wave(previous, INSERTION as usize, row - 1),
            _mm256_set1_epi32(config.gap_extend_score),
        );
        let (insertion_raw, insertion_previous) = choose_avx2(above_m, above_d, above_i);
        let insertion_positive = _mm256_cmpgt_epi32(insertion_raw, zero);
        let insertion_scores = _mm256_and_si256(insertion_raw, insertion_positive);

        let wave_offset = row - current_range.start;
        let cell_offset = WAVE_PAD + wave_offset;
        _mm256_storeu_si256(current.as_mut_ptr().add(cell_offset).cast(), match_scores);
        _mm256_storeu_si256(
            current.as_mut_ptr().add(stride + cell_offset).cast(),
            deletion_scores,
        );
        _mm256_storeu_si256(
            current.as_mut_ptr().add(2 * stride + cell_offset).cast(),
            insertion_scores,
        );

        if TRACEBACK {
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
            // lanes 0..3 and 4..7 in separate 128-bit halves; the permutation joins them in lane
            // order.
            let packed = _mm256_permute4x64_epi64::<0xD8>(_mm256_packus_epi32(trace, zero));
            debug_assert!(wave_base + wave_offset + 8 <= traceback.len());
            // SAFETY: the wave layout stores this wave's cells contiguously, and rows row..=row+7
            // lie inside current_range, so the eight cells are inside the laid-out traceback.
            _mm_storeu_si128(
                traceback.as_mut_ptr().add(wave_base + wave_offset).cast(),
                _mm256_castsi256_si128(packed),
            );
        }
        let (lane_best_scores, _) = choose_avx2(match_scores, deletion_scores, insertion_scores);
        let mut block_best = _mm256_max_epi32(
            lane_best_scores,
            _mm256_permute2x128_si256::<0x01>(lane_best_scores, lane_best_scores),
        );
        block_best = _mm256_max_epi32(block_best, _mm256_shuffle_epi32::<0x4e>(block_best));
        block_best = _mm256_max_epi32(block_best, _mm256_shuffle_epi32::<0xb1>(block_best));
        let block_score = _mm256_extract_epi32::<0>(block_best);
        // BestCell only changes for a higher score or an equal, earlier cell, so a block below
        // the current best cannot change it.
        if block_score > 0 && block_score >= best.score {
            let mut matches = [0i32; 8];
            let mut deletions = [0i32; 8];
            let mut insertions = [0i32; 8];
            _mm256_storeu_si256(matches.as_mut_ptr().cast(), match_scores);
            _mm256_storeu_si256(deletions.as_mut_ptr().cast(), deletion_scores);
            _mm256_storeu_si256(insertions.as_mut_ptr().cast(), insertion_scores);
            let best_lanes = _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpeq_epi32(
                lane_best_scores,
                block_best,
            )));
            let block_lane = best_lanes.trailing_zeros() as usize;
            // Lanes are in ascending query-index order. The first lane with the block maximum
            // is therefore lexicographically earliest; every lesser or later equal lane loses
            // after this unchanged BestCell comparison.
            let query_index = row + block_lane;
            let contig_index = wave - query_index;
            best.consider(
                query_index,
                contig_index,
                Cell {
                    scores: [
                        matches[block_lane],
                        deletions[block_lane],
                        insertions[block_lane],
                    ],
                    previous: [0; 3],
                },
            );
        }
    }
}
