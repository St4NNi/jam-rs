//! Filling the band wave by wave: the scalar and vector recurrences over one wave, the rotation
//! of the three wave buffers, the traceback layout of a chunk and its checkpoints.

use super::geometry::{BestCell, WaveCarry, WaveRange, wave_bounds, wave_range};
#[cfg(test)]
use super::rebased::TEST_HOOKS;
use super::rebased::{Offset16, RebaseState, SENTINEL_BASE};
use super::waves::{ScoreWave, WaveScore, fill_wave_scalar};
use super::{AlignmentConfig, AlignmentWorkspace};
use crate::align::AlignmentError;

/// The score element of a pass: i16 offsets from a moving base, which only vector passes use
/// and only when the configuration allows them, or plain i32 scores.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum ScoreWidth {
    I32,
    I16,
}

/// A pass with i16 offsets met a value they cannot hold; the pass reruns with i32 scores.
#[derive(Debug)]
pub(super) struct RebaseViolation;

/// The constants of one pass over the band.
#[derive(Clone, Copy)]
pub(super) struct Pass<'a> {
    pub(super) query: &'a [u8],
    pub(super) contig: &'a [u8],
    pub(super) config: AlignmentConfig,
    pub(super) stride: usize,
}

impl AlignmentWorkspace {
    /// Runs `fill_waves` with i32 planes and eight AVX2 lanes in one out-of-line body, so the
    /// drivers share one inlined copy of the vector kernel.
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    #[inline(never)]
    unsafe fn fill_i32_avx2<const TRACEBACK: bool, const BEST: bool>(
        &mut self,
        pass: Pass<'_>,
        waves: (usize, usize),
        carry: &mut WaveCarry,
        best: &mut BestCell,
    ) -> Result<(), RebaseViolation> {
        // SAFETY: this function enables AVX2 for the vector instantiation.
        unsafe { self.fill_waves::<i32, true, TRACEBACK, BEST>(pass, waves, carry, best) }
    }

    /// Runs `fill_waves` with rebased i16 planes and sixteen AVX2 lanes. Callers must check the
    /// configuration with `rebased_local_scores`.
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    #[inline(never)]
    unsafe fn fill_i16_avx2<const TRACEBACK: bool, const BEST: bool>(
        &mut self,
        pass: Pass<'_>,
        waves: (usize, usize),
        carry: &mut WaveCarry,
        best: &mut BestCell,
    ) -> Result<(), RebaseViolation> {
        // SAFETY: this function enables AVX2 and the caller checks the configuration.
        unsafe { self.fill_waves::<Offset16, true, TRACEBACK, BEST>(pass, waves, carry, best) }
    }

    /// Fills waves with scores of `width`; BEST false lets the i16 kernel skip the best cell,
    /// which the others always track. Callers must enable AVX2 when VECTOR is true, and pass
    /// `I16` only for a vector pass whose configuration `rebased_local_scores` accepts.
    #[inline(always)]
    pub(super) unsafe fn fill_waves_selected<
        const VECTOR: bool,
        const TRACEBACK: bool,
        const BEST: bool,
    >(
        &mut self,
        pass: Pass<'_>,
        width: ScoreWidth,
        waves: (usize, usize),
        carry: &mut WaveCarry,
        best: &mut BestCell,
    ) -> Result<(), RebaseViolation> {
        #[cfg(target_arch = "x86_64")]
        if VECTOR {
            // SAFETY: forwarded from this function's contract.
            return unsafe {
                match width {
                    ScoreWidth::I16 => {
                        self.fill_i16_avx2::<TRACEBACK, BEST>(pass, waves, carry, best)
                    }
                    ScoreWidth::I32 => {
                        self.fill_i32_avx2::<TRACEBACK, BEST>(pass, waves, carry, best)
                    }
                }
            };
        }
        debug_assert_eq!(width, ScoreWidth::I32);
        // SAFETY: the scalar instantiation never calls AVX2 code.
        unsafe { self.fill_waves::<i32, false, TRACEBACK, BEST>(pass, waves, carry, best) }
    }

    /// Fills waves `first..=last`. The wave buffers and the two predecessor ranges rotate in
    /// registers: reloading them from memory that the previous wave just wrote stalls on store
    /// forwarding. Callers must enable AVX2 when VECTOR is true.
    #[inline(always)]
    pub(super) unsafe fn fill_waves<
        S: WaveScore,
        const VECTOR: bool,
        const TRACEBACK: bool,
        const BEST: bool,
    >(
        &mut self,
        pass: Pass<'_>,
        (first_wave, last_wave): (usize, usize),
        carry: &mut WaveCarry,
        best: &mut BestCell,
    ) -> Result<(), RebaseViolation> {
        if first_wave > last_wave {
            return Ok(());
        }
        let (query_len, contig_len, stride) = (pass.query.len(), pass.contig.len(), pass.stride);
        let active = wave_bounds(query_len, contig_len, pass.config).and_then(|(first, last)| {
            let first = first.max(first_wave);
            let last = last.min(last_wave);
            (first <= last).then_some((first, last))
        });
        let Some((active_first, active_last)) = active else {
            self.skip_empty_waves(last_wave - first_wave + 1, carry);
            return Ok(());
        };
        self.skip_empty_waves(active_first - first_wave, carry);
        let mut rebase = RebaseState::new(carry.base);
        #[cfg(test)]
        let forced_violation = TEST_HOOKS.get().violation_wave;
        let [first_buffer, second_buffer, third_buffer] = &mut self.waves;
        let (mut older_buffer, mut previous_buffer, mut current_buffer) =
            (first_buffer, second_buffer, third_buffer);
        let mut older_range = carry.older;
        let mut previous_range = carry.previous;
        for wave in active_first..=active_last {
            let current_range = wave_range(query_len, contig_len, pass.config, wave);
            if let Some(current_range) = current_range {
                if S::REBASED {
                    rebase.start_wave();
                    // The next wave loads these cells with vector loads, so they are written a
                    // whole wave early: a load that overlaps a store still in flight stalls.
                    rebase.write_pads(
                        bytemuck::cast_slice_mut(&mut current_buffer[..]),
                        current_range,
                        stride,
                    );
                }
                let wave_base = if TRACEBACK {
                    self.wave_offsets[wave - first_wave]
                } else {
                    0
                };
                // An i16 view of an i32 wave buffer holds the same three planes of stride cells.
                let older_scores: &[S] = bytemuck::cast_slice(&older_buffer[..]);
                let previous_scores: &[S] = bytemuck::cast_slice(&previous_buffer[..]);
                let older = older_range.map(|range| ScoreWave {
                    range,
                    scores: older_scores,
                    stride,
                });
                let previous = previous_range.map(|range| ScoreWave {
                    range,
                    scores: previous_scores,
                    stride,
                });
                let cells = WaveCells {
                    wave,
                    range: current_range,
                    wave_base,
                    older,
                    previous,
                };
                // SAFETY: forwarded from this function's contract.
                unsafe {
                    fill_wave::<S, VECTOR, TRACEBACK, BEST>(
                        pass,
                        cells,
                        bytemuck::cast_slice_mut(&mut current_buffer[..]),
                        &mut self.compact_cells,
                        best,
                        &mut rebase,
                    );
                }
                if S::REBASED {
                    rebase.finish_wave(
                        bytemuck::cast_slice_mut(&mut current_buffer[..]),
                        current_range,
                        bytemuck::cast_slice_mut(&mut previous_buffer[..]),
                        previous_range,
                        stride,
                    );
                    #[cfg(test)]
                    if forced_violation.is_some_and(|after| wave >= after) {
                        rebase.violation = true;
                    }
                    if rebase.violation {
                        carry.base = rebase.base;
                        return Err(RebaseViolation);
                    }
                }
            }
            (older_buffer, previous_buffer, current_buffer) =
                (previous_buffer, current_buffer, older_buffer);
            older_range = previous_range;
            previous_range = current_range;
        }
        carry.older = older_range;
        carry.previous = previous_range;
        self.waves.rotate_left((active_last - active_first + 1) % 3);
        carry.base = rebase.base;
        self.skip_empty_waves(last_wave - active_last, carry);
        Ok(())
    }

    /// Empty waves still rotate score buffers and age the two predecessor ranges.
    pub(super) fn skip_empty_waves(&mut self, count: usize, carry: &mut WaveCarry) {
        if count > 0 {
            self.waves.rotate_left(count % 3);
            carry.older = if count == 1 { carry.previous } else { None };
            carry.previous = None;
        }
    }

    /// Wave layout of the traceback cells in waves first..=last: each wave's cells are
    /// contiguous, so one vector block writes its traces together. Returns the cell count.
    pub(super) fn layout_chunk(
        &mut self,
        query_len: usize,
        contig_len: usize,
        config: AlignmentConfig,
        first_wave: usize,
        last_wave: usize,
    ) -> Result<usize, AlignmentError> {
        let mut total = 0usize;
        self.wave_offsets.clear();
        let active = wave_bounds(query_len, contig_len, config).and_then(|(first, last)| {
            let first = first.max(first_wave);
            let last = last.min(last_wave);
            (first <= last).then_some((first, last))
        });
        if let Some((first, last)) = active {
            self.wave_offsets.resize(first - first_wave, 0);
            for wave in first..=last {
                self.wave_offsets.push(total);
                let width =
                    wave_range(query_len, contig_len, config, wave).map_or(0, WaveRange::width);
                total = total
                    .checked_add(width)
                    .ok_or(AlignmentError::LengthOverflow)?;
            }
        }
        let offsets = last_wave
            .checked_sub(first_wave)
            .map_or(Some(1), |length| length.checked_add(2))
            .ok_or(AlignmentError::LengthOverflow)?;
        self.wave_offsets.resize(offsets, total);
        self.compact_cells.resize(total, 0);
        Ok(total)
    }

    pub(super) fn save_checkpoint(&mut self, carry: WaveCarry, stride: usize) {
        self.chunk_carries.push(carry);
        for wave in &self.waves[..2] {
            let scores: &[i32] = bytemuck::cast_slice(&wave[..]);
            self.checkpoints.extend_from_slice(&scores[..3 * stride]);
        }
    }
}

/// One wave to fill: its rows, where its traceback cells start, and the two waves before it.
struct WaveCells<'a, S> {
    wave: usize,
    range: WaveRange,
    wave_base: usize,
    older: Option<ScoreWave<'a, S>>,
    previous: Option<ScoreWave<'a, S>>,
}

/// Fills the cells of one wave: vector blocks over the rows whose neighbours all lie in the two
/// waves before, scalar cells at the edges. Callers must enable AVX2 when VECTOR is true.
#[inline(always)]
unsafe fn fill_wave<S: WaveScore, const VECTOR: bool, const TRACEBACK: bool, const BEST: bool>(
    pass: Pass<'_>,
    cells: WaveCells<'_, S>,
    current: &mut [S],
    traceback: &mut [u16],
    best: &mut BestCell,
    rebase: &mut RebaseState,
) {
    let Pass {
        query,
        contig,
        config,
        stride,
    } = pass;
    let WaveCells {
        wave,
        range: current_range,
        wave_base,
        older,
        previous,
    } = cells;
    let gap_open_score = config.gap_open();
    let scalar = |row: usize,
                  current: &mut [S],
                  traceback: &mut [u16],
                  best: &mut BestCell,
                  rebase: &mut RebaseState| {
        fill_wave_scalar::<S, TRACEBACK>(
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
            current,
            traceback,
            wave_base,
            best,
            rebase,
        );
    };
    let mut row = current_range.start;
    #[cfg(target_arch = "x86_64")]
    if let Some((start, end, older, previous)) =
        vector_span::<S, VECTOR>(wave, current_range, older, previous, rebase)
    {
        while row < start {
            scalar(row, current, traceback, best, rebase);
            row += 1;
        }
        // SAFETY: vector_span proves LANES current, diagonal, left, and above cells of every
        // block are in their wave slices. start>=1 and end<=wave-1 prove the query and
        // reverse-loaded contig bytes are also in bounds.
        row = unsafe {
            S::fill_blocks::<TRACEBACK, BEST>(
                query,
                contig,
                config,
                gap_open_score,
                wave,
                (start, end),
                current_range,
                stride,
                older,
                previous,
                current,
                traceback,
                wave_base,
                best,
                rebase,
            )
        };
    }
    while row <= current_range.end {
        scalar(row, current, traceback, best, rebase);
        row += 1;
    }
}

/// The rows of a wave whose three neighbours all lie in the two waves before, for vector blocks.
#[cfg(target_arch = "x86_64")]
#[inline(always)]
#[allow(clippy::type_complexity)]
fn vector_span<'a, S: WaveScore, const VECTOR: bool>(
    wave: usize,
    current_range: WaveRange,
    older: Option<ScoreWave<'a, S>>,
    previous: Option<ScoreWave<'a, S>>,
    rebase: &RebaseState,
) -> Option<(usize, usize, ScoreWave<'a, S>, ScoreWave<'a, S>)> {
    older
        .zip(previous)
        .filter(|_| VECTOR)
        .and_then(|(older, previous)| {
            // Rebased planes read a zeroed spare cell one row outside the previous wave. A
            // sentinel-range base would take it for a low candidate, so those waves keep their
            // edge cells scalar.
            let (previous_low, previous_high) = if S::REBASED && rebase.base <= SENTINEL_BASE {
                (previous.range.start, previous.range.end.saturating_add(1))
            } else {
                (previous.range.start.saturating_add(1), previous.range.end)
            };
            let start = current_range
                .start
                .max(older.range.start.saturating_add(1))
                .max(previous_low)
                .max(1);
            let end = current_range
                .end
                .min(older.range.end.saturating_add(1))
                .min(previous_high)
                .min(wave.saturating_sub(1));
            (start <= end).then_some((start, end, older, previous))
        })
}
