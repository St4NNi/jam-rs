//! Score planes of the wave kernel and the scalar recurrence of one cell.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

// Each wave keeps three planes of `stride` elements, with `WAVE_PAD` spare elements on each side
// so vector blocks can read one neighbour outside the previous wave.

use super::AlignmentConfig;
use super::geometry::{
    BestCell, Cell, DELETION, INSERTION, MATCH, START, WaveRange, best_score, choose,
    encode_traceback,
};
use super::rebased::RebaseState;

/// Spare elements before and after the rows of each plane. 32 keeps the rows of every plane on
/// a 64-byte boundary.
pub(super) const WAVE_PAD: usize = 32;

/// Sixteen i32 on a 64-byte boundary. The wave buffers are chunk vectors, so the planes start on
/// a boundary and a vector load or store of one block never straddles a cache line.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, bytemuck::Pod, bytemuck::Zeroable)]
#[repr(C, align(64))]
pub(super) struct WaveChunk(pub(super) [i32; 16]);

impl WaveChunk {
    pub(super) const LEN: usize = 16;
}

/// Grows one wave buffer to `length` i32 elements without shrinking its capacity.
pub(super) fn prepare_wave(scores: &mut Vec<WaveChunk>, length: usize) {
    let chunks = length.div_ceil(WaveChunk::LEN);
    if scores.capacity() < chunks {
        scores.reserve_exact(chunks - scores.len());
    }
    scores.resize(chunks, WaveChunk::default());
}

/// One finished wave as the next waves read it.
#[derive(Clone, Copy)]
pub(super) struct ScoreWave<'a, S = i32> {
    pub(super) range: WaveRange,
    pub(super) scores: &'a [S],
    pub(super) stride: usize,
}

impl<S: WaveScore> ScoreWave<'_, S> {
    /// The three state scores of `row`, if the wave holds it.
    pub(super) fn cell(self, row: usize, base: i32) -> Option<[i32; 3]> {
        self.range.contains(row).then(|| {
            let offset = WAVE_PAD + row - self.range.start;
            [
                self.scores[offset].read(base),
                self.scores[self.stride + offset].read(base),
                self.scores[2 * self.stride + offset].read(base),
            ]
        })
    }
}

/// Score plane element of the wave recurrence: plain i32 scores, or i16 offsets from a base.
pub(super) trait WaveScore: bytemuck::Pod {
    /// Cells of one vector block.
    const LANES: usize;
    /// True for the rebased i16 planes, whose scores are offsets from `WaveCarry::base`.
    const REBASED: bool = false;

    /// The score of a stored element under `base`.
    fn read(self, base: i32) -> i32;

    /// The element that stores `score` under `base`.
    fn write(score: i32, base: i32) -> Self;

    /// Evaluates rows `start..=end` of one wave in vector blocks of LANES cells, and returns the
    /// first row it did not fill (`start` when the span is shorter than a block). Callers must
    /// enable AVX2 and prove the block bounds that `fill_wave_avx2` documents.
    #[cfg(target_arch = "x86_64")]
    #[allow(clippy::too_many_arguments)]
    unsafe fn fill_blocks<const TRACEBACK: bool, const BEST: bool>(
        query: &[u8],
        contig: &[u8],
        config: AlignmentConfig,
        gap_open_score: i32,
        wave: usize,
        span: (usize, usize),
        current_range: WaveRange,
        stride: usize,
        older: ScoreWave<'_, Self>,
        previous: ScoreWave<'_, Self>,
        current: &mut [Self],
        traceback: &mut [u16],
        wave_base: usize,
        best: &mut BestCell,
        rebase: &mut RebaseState,
    ) -> usize;
}

impl WaveScore for i32 {
    const LANES: usize = 8;

    fn read(self, _base: i32) -> i32 {
        self
    }

    fn write(score: i32, _base: i32) -> Self {
        score
    }

    /// Blocks of eight, the last one overlapping finished cells (a repeat writes the same
    /// values).
    #[cfg(target_arch = "x86_64")]
    #[inline(always)]
    unsafe fn fill_blocks<const TRACEBACK: bool, const BEST: bool>(
        query: &[u8],
        contig: &[u8],
        config: AlignmentConfig,
        gap_open_score: i32,
        wave: usize,
        (start, end): (usize, usize),
        current_range: WaveRange,
        stride: usize,
        older: ScoreWave<'_, Self>,
        previous: ScoreWave<'_, Self>,
        current: &mut [Self],
        traceback: &mut [u16],
        wave_base: usize,
        best: &mut BestCell,
        _rebase: &mut RebaseState,
    ) -> usize {
        let block = |row| super::avx2::Block {
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
        };
        let mut row = start;
        // SAFETY: the caller proves every block is inside the span, as `fill_wave_avx2` requires.
        unsafe {
            while row
                .checked_add(Self::LANES - 1)
                .is_some_and(|last| last <= end)
            {
                super::avx2::fill_wave_avx2::<TRACEBACK>(block(row), current, traceback, best);
                row += Self::LANES;
            }
            if row <= end && end + 1 - start >= Self::LANES {
                let last = end + 1 - Self::LANES;
                super::avx2::fill_wave_avx2::<TRACEBACK>(block(last), current, traceback, best);
                row = end + 1;
            }
        }
        row
    }
}

/// Fills one cell of a wave with the scalar recurrence.
#[allow(clippy::too_many_arguments)]
pub(super) fn fill_wave_scalar<S: WaveScore, const TRACEBACK: bool>(
    query: &[u8],
    contig: &[u8],
    config: AlignmentConfig,
    gap_open_score: i32,
    wave: usize,
    row: usize,
    current_range: WaveRange,
    stride: usize,
    older: Option<ScoreWave<'_, S>>,
    previous: Option<ScoreWave<'_, S>>,
    current: &mut [S],
    traceback: &mut [u16],
    wave_base: usize,
    best: &mut BestCell,
    rebase: &mut RebaseState,
) {
    let base = rebase.base;
    let contig_index = wave - row;
    let wave_offset = WAVE_PAD + row - current_range.start;
    let trace_index = wave_base + row - current_range.start;
    if row == 0 && contig_index == 0 {
        current[wave_offset] = S::write(0, base);
        current[stride + wave_offset] = S::write(0, base);
        current[2 * stride + wave_offset] = S::write(0, base);
        if TRACEBACK {
            traceback[trace_index] = 0;
        }
        return;
    }
    let mut cell = Cell::default();
    if row > 0
        && contig_index > 0
        && let Some(previous) = older.and_then(|wave| wave.cell(row - 1, base))
    {
        let (score, state) = best_score(previous);
        let substitution = if query[row - 1].eq_ignore_ascii_case(&contig[contig_index - 1]) {
            config.match_score
        } else {
            config.mismatch_score
        };
        let score = score.saturating_add(substitution);
        if score > 0 {
            cell.scores[MATCH as usize] = score;
            cell.previous[MATCH as usize] = if score == substitution { START } else { state };
        }
    }
    if let Some(previous) = previous.and_then(|wave| wave.cell(row, base)) {
        let (score, state) = choose([
            (
                previous[DELETION as usize] + config.gap_extend_score,
                DELETION,
            ),
            (previous[MATCH as usize] + gap_open_score, MATCH),
            (previous[INSERTION as usize] + gap_open_score, INSERTION),
        ]);
        if score > 0 {
            cell.scores[DELETION as usize] = score;
            cell.previous[DELETION as usize] = state;
        }
    }
    if row > 0
        && let Some(previous) = previous.and_then(|wave| wave.cell(row - 1, base))
    {
        let (score, state) = choose([
            (
                previous[INSERTION as usize] + config.gap_extend_score,
                INSERTION,
            ),
            (previous[MATCH as usize] + gap_open_score, MATCH),
            (previous[DELETION as usize] + gap_open_score, DELETION),
        ]);
        if score > 0 {
            cell.scores[INSERTION as usize] = score;
            cell.previous[INSERTION as usize] = state;
        }
    }
    current[wave_offset] = S::write(cell.scores[MATCH as usize], base);
    current[stride + wave_offset] = S::write(cell.scores[DELETION as usize], base);
    current[2 * stride + wave_offset] = S::write(cell.scores[INSERTION as usize], base);
    if S::REBASED {
        rebase.scalar_cell(&cell);
    }
    if TRACEBACK {
        traceback[trace_index] = encode_traceback(cell);
    }
    best.consider(row, contig_index, cell);
}
