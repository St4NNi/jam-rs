//! Filling the band wave by wave: the scalar and vector recurrences over one wave, the rotation
//! of the three wave buffers, the traceback layout of a chunk and its checkpoints.

use super::geometry::{BestCell, WaveCarry, WaveRange, wave_bounds, wave_range};
#[cfg(test)]
use super::rebased::REBASE_FORCED_VIOLATION;
use super::rebased::RebaseState;
use super::waves::{ScoreWave, WaveScore, fill_wave_scalar};
use super::{AlignmentConfig, AlignmentWorkspace};
use crate::align::AlignmentError;

/// The element type of a vector pass: rebased i16 offsets when the configuration allows them,
/// i32 scores otherwise.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum LocalKernel {
    Wide,
    Rebased,
}

/// The constants of one pass over the band.
#[derive(Clone, Copy)]
pub(super) struct Pass<'a> {
    pub(super) query: &'a [u8],
    pub(super) contig: &'a [u8],
    pub(super) config: AlignmentConfig,
    pub(super) stride: usize,
}

impl AlignmentWorkspace {
    /// Fills waves with the selected kernel. Callers must enable AVX2 when VECTOR is true, and
    /// pass the rebased kernel only when `select_kernel` chose it.
    #[inline(always)]
    pub(super) unsafe fn fill_waves_selected<
        const VECTOR: bool,
        const TRACE: bool,
        const BEST: bool,
    >(
        &mut self,
        pass: Pass<'_>,
        kernel: LocalKernel,
        waves: (usize, usize),
        carry: &mut WaveCarry,
        best: &mut BestCell,
    ) {
        debug_assert_eq!(kernel, LocalKernel::Wide);
        // SAFETY: the scalar instantiation never calls AVX2 code.
        unsafe { self.fill_waves::<i32, false, TRACE, BEST>(pass, waves, carry, best) }
    }

    /// Fills waves `first..=last`. The wave buffers and the two predecessor ranges rotate in
    /// registers: reloading them from memory that the previous wave just wrote stalls on store
    /// forwarding. Callers must enable AVX2 when VECTOR is true.
    #[inline(always)]
    pub(super) unsafe fn fill_waves<
        S: WaveScore,
        const VECTOR: bool,
        const TRACE: bool,
        const BEST: bool,
    >(
        &mut self,
        pass: Pass<'_>,
        (first_wave, last_wave): (usize, usize),
        carry: &mut WaveCarry,
        best: &mut BestCell,
    ) {
        if first_wave > last_wave {
            return;
        }
        let (query_len, contig_len, stride) = (pass.query.len(), pass.contig.len(), pass.stride);
        let active = wave_bounds(query_len, contig_len, pass.config).and_then(|(first, last)| {
            let first = first.max(first_wave);
            let last = last.min(last_wave);
            (first <= last).then_some((first, last))
        });
        let Some((active_first, active_last)) = active else {
            self.skip_empty_waves(last_wave - first_wave + 1, carry);
            return;
        };
        self.skip_empty_waves(active_first - first_wave, carry);
        let mut rebase = RebaseState::new(carry.base);
        #[cfg(test)]
        let forced_violation = REBASE_FORCED_VIOLATION.get();
        let interrupt = &self.interrupt;
        let [first_buffer, second_buffer, third_buffer] = &mut self.waves;
        let (mut older_buffer, mut previous_buffer, mut current_buffer) =
            (first_buffer, second_buffer, third_buffer);
        let mut older_range = carry.older;
        let mut previous_range = carry.previous;
        for wave in active_first..=active_last {
            if wave % 64 == 0
                && interrupt
                    .as_ref()
                    .is_some_and(|flag| flag.load(std::sync::atomic::Ordering::Relaxed))
            {
                return;
            }
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
                let wave_base = if TRACE {
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
                    fill_wave::<S, VECTOR, TRACE, BEST>(
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
                        self.rebase_violation = true;
                        carry.base = rebase.base;
                        return;
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
unsafe fn fill_wave<S: WaveScore, const VECTOR: bool, const TRACE: bool, const BEST: bool>(
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
        fill_wave_scalar::<S, TRACE>(
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
    while row <= current_range.end {
        scalar(row, current, traceback, best, rebase);
        row += 1;
    }
}
