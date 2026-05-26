//! The rebased 16-bit score representation. See `Offset16` for the representation and its
//! exactness argument; `rebased_avx2` holds its vector recurrence.

use super::AlignmentConfig;
use super::geometry::{Cell, WaveRange};
use super::waves::WAVE_PAD;

#[cfg(test)]
thread_local! {
    /// Replaces the rebase thresholds (high, settle, low), so short inputs rebase often.
    pub(super) static REBASE_THRESHOLDS: std::cell::Cell<Option<(i32, i32, i32)>> =
        const { std::cell::Cell::new(None) };
    /// Makes the rebased pass report a violation after this many waves.
    pub(super) static REBASE_FORCED_VIOLATION: std::cell::Cell<Option<usize>> =
        const { std::cell::Cell::new(None) };
    /// Largest base, shifts up, shifts down and i32 reruns of the rebased passes on this thread.
    pub(super) static REBASE_REACHED: std::cell::Cell<(i32, usize, usize, usize)> =
        const { std::cell::Cell::new((0, 0, 0, 0)) };
}

/// Largest per-wave score spread the rebased kernel accepts a configuration for.
const REBASE_MAX_SPREAD: i64 = 8_000;

/// Sentinel offset of a stored score that is zero, or too far below the base to represent.
pub(super) const OFFSET_FLOOR: i16 = i16::MIN;
/// Bases above this cannot store a zero score exactly, so zero states use `OFFSET_FLOOR`.
pub(super) const SENTINEL_BASE: i32 = 32_768;
/// A winning candidate below this offset, under a sentinel-range base, is not representable.
const LOW_VIOLATION: i32 = -30_000;
/// Wave maxima above `REBASE_HIGH` move the base up so the maximum lands on `REBASE_SETTLE`;
/// with a nonzero base, maxima below `REBASE_LOW` move it down.
const REBASE_HIGH: i32 = 24_000;
const REBASE_SETTLE: i32 = 8_000;
const REBASE_LOW: i32 = 4_000;

/// i16 offsets from a per-pass base; a wave's best scores lie within the spread that `rebased_local_scores` bounds.
/// A stored zero is `-base` or `OFFSET_FLOOR`; only band-edge gap states are zero next to large scores,
/// and another candidate wins there. An unrepresentable winner sets `RebaseState::violation` (rerun on i32).
#[derive(Clone, Copy, Debug, bytemuck::Pod, bytemuck::Zeroable)]
#[repr(transparent)]
pub(super) struct Offset16(pub(super) i16);

impl super::waves::WaveScore for Offset16 {
    const LANES: usize = 16;
    const REBASED: bool = true;

    fn read(self, base: i32) -> i32 {
        if self.0 == OFFSET_FLOOR {
            0
        } else {
            i32::from(self.0) + base
        }
    }

    fn write(score: i32, base: i32) -> Self {
        if score <= 0 {
            Self(zero_offset(base))
        } else {
            Self((score - base).clamp(-32_767, 32_767) as i16)
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[inline(always)]
    unsafe fn fill_blocks<const TRACE: bool, const BEST: bool>(
        query: &[u8],
        contig: &[u8],
        config: AlignmentConfig,
        gap_open_score: i32,
        wave: usize,
        span: (usize, usize),
        current_range: WaveRange,
        stride: usize,
        older: super::waves::ScoreWave<'_, Self>,
        previous: super::waves::ScoreWave<'_, Self>,
        current: &mut [Self],
        traceback: &mut [u16],
        wave_base: usize,
        best: &mut super::geometry::BestCell,
        rebase: &mut RebaseState,
    ) -> usize {
        let blocks = super::rebased_avx2::Span {
            query,
            contig,
            config,
            gap_open_score,
            wave,
            span,
            current_range,
            stride,
            older,
            previous,
            wave_base,
        };
        // SAFETY: forwarded from this function's contract.
        unsafe {
            super::rebased_avx2::fill_span16_rebased::<TRACE, BEST>(
                blocks, current, traceback, best, rebase,
            )
        }
    }
}

/// Whether the rebased i16 kernel is exact for this configuration on any window length: the
/// scoring terms fit i16 and one wave's score spread stays below the bound. Two cells of a wave
/// join through a gap of at most `2 * band + 1` positions, costing at most
/// `open + extend * (2 * band + 1)`, plus `match` per position.
pub(super) fn rebased_local_scores(config: AlignmentConfig) -> bool {
    let positions = 2 * i64::from(config.band_width) + 2;
    let spread = positions
        .checked_mul(i64::from(config.match_score) + i64::from(config.gap_extend_score).abs())
        .and_then(|spread| {
            spread.checked_add(
                i64::from(config.gap_open_score).abs()
                    + i64::from(config.gap_extend_score).abs()
                    + i64::from(config.mismatch_score).abs(),
            )
        });
    config.validate().is_ok()
        && spread.is_some_and(|spread| spread <= REBASE_MAX_SPREAD)
        && [
            config.match_score,
            config.mismatch_score,
            config.gap_open(),
            config.gap_extend_score,
        ]
        .iter()
        .all(|&score| i16::try_from(score).is_ok())
}

/// The offset that stores a zero score under `base`.
pub(super) fn zero_offset(base: i32) -> i16 {
    (-base).max(i32::from(i16::MIN)) as i16
}

/// Per-pass state of the rebased kernel: the base, and what the current wave shows about it.
pub(super) struct RebaseState {
    pub(super) base: i32,
    /// Largest stored best-state offset and smallest winning raw offset of the vector blocks of
    /// the current wave, as `i16::MIN` and `i16::MAX` when there were none.
    pub(super) vector_max: i32,
    pub(super) vector_low: i32,
    /// The same for scalar cells: largest best-state offset (`i32::MIN` if none).
    scalar_max: i32,
    /// Set when a value is not representable; the pass result must be discarded.
    pub(super) violation: bool,
    /// Whether the current wave tracks its smallest winning candidate: only a base in the
    /// sentinel range can meet an unrepresentable one.
    pub(super) check_low: bool,
    high: i32,
    settle: i32,
    low: i32,
}

impl RebaseState {
    pub(super) fn new(base: i32) -> Self {
        #[cfg(test)]
        let (high, settle, low) =
            REBASE_THRESHOLDS
                .get()
                .unwrap_or((REBASE_HIGH, REBASE_SETTLE, REBASE_LOW));
        #[cfg(not(test))]
        let (high, settle, low) = (REBASE_HIGH, REBASE_SETTLE, REBASE_LOW);
        Self {
            base,
            vector_max: i32::from(i16::MIN),
            vector_low: i32::from(i16::MAX),
            scalar_max: i32::MIN,
            violation: false,
            check_low: false,
            high,
            settle,
            low,
        }
    }

    pub(super) fn start_wave(&mut self) {
        self.vector_max = i32::from(i16::MIN);
        self.vector_low = i32::from(i16::MAX);
        self.scalar_max = i32::MIN;
        self.check_low = self.base > SENTINEL_BASE;
    }

    /// Records a scalar cell, whose scores are `cell.scores` in the score domain.
    pub(super) fn scalar_cell(&mut self, cell: &Cell) {
        let best = cell.scores.into_iter().max().unwrap_or(0);
        if best > 0 {
            self.scalar_max = self.scalar_max.max(best - self.base);
        }
        for score in cell.scores {
            if score > 0 {
                let offset = score - self.base;
                if offset > i32::from(i16::MAX)
                    || (self.base > SENTINEL_BASE && offset < LOW_VIOLATION)
                {
                    self.violation = true;
                }
            }
        }
    }

    /// Ends a wave: checks representability and moves the base to keep the wave maximum in range.
    /// `current` and `previous` are the two waves the next wave reads.
    #[inline(always)]
    pub(super) fn finish_wave(
        &mut self,
        current: &mut [Offset16],
        current_range: WaveRange,
        previous: &mut [Offset16],
        previous_range: Option<WaveRange>,
        stride: usize,
    ) {
        let wave_max = if self.vector_max > i32::from(i16::MIN) {
            self.vector_max.max(self.scalar_max)
        } else {
            self.scalar_max
        };
        if self.base > SENTINEL_BASE && self.vector_low < LOW_VIOLATION {
            self.violation = true;
        }
        if self.violation || wave_max == i32::MIN {
            return;
        }
        if wave_max > self.high {
            let delta = wave_max - self.settle;
            self.shift(
                current,
                current_range,
                previous,
                previous_range,
                stride,
                delta,
                true,
            );
            self.base += delta;
            self.write_pads(current, current_range, stride);
            #[cfg(test)]
            REBASE_REACHED.set({
                let (base, up, down, fallbacks) = REBASE_REACHED.get();
                (base.max(self.base), up + 1, down, fallbacks)
            });
        } else if self.base > 0 && wave_max < self.low {
            let delta = (self.settle - wave_max).min(self.base);
            self.shift(
                current,
                current_range,
                previous,
                previous_range,
                stride,
                delta,
                false,
            );
            self.base -= delta;
            self.write_pads(current, current_range, stride);
            #[cfg(test)]
            REBASE_REACHED.set({
                let (base, up, down, fallbacks) = REBASE_REACHED.get();
                (base, up, down + 1, fallbacks)
            });
        }
    }

    /// Sets the spare cells around each plane of `wave` to zero under the current base. The next
    /// wave reads them as its neighbours outside this wave, so they act as absent cells.
    #[inline(always)]
    pub(super) fn write_pads(&self, wave: &mut [Offset16], range: WaveRange, stride: usize) {
        let zero = Offset16(zero_offset(self.base));
        for state in 0..3 {
            wave[state * stride + WAVE_PAD - 1] = zero;
            wave[state * stride + WAVE_PAD + range.width()] = zero;
        }
    }

    /// Re-expresses both waves under a base moved by `delta` (up: scores stay, offsets fall).
    #[allow(clippy::too_many_arguments)]
    fn shift(
        &self,
        current: &mut [Offset16],
        current_range: WaveRange,
        previous: &mut [Offset16],
        previous_range: Option<WaveRange>,
        stride: usize,
        delta: i32,
        up: bool,
    ) {
        let new_base = if up {
            self.base + delta
        } else {
            self.base - delta
        };
        for (wave, range) in [(current, Some(current_range)), (previous, previous_range)] {
            let Some(range) = range else { continue };
            for state in 0..3 {
                let first = state * stride + WAVE_PAD;
                for element in &mut wave[first..first + range.width()] {
                    element.0 = if up {
                        if element.0 == OFFSET_FLOOR {
                            OFFSET_FLOOR
                        } else {
                            (i32::from(element.0) - delta).max(i32::from(i16::MIN)) as i16
                        }
                    } else if element.0 == OFFSET_FLOOR {
                        zero_offset(new_base)
                    } else {
                        (i32::from(element.0) + delta).min(i32::from(i16::MAX)) as i16
                    };
                }
            }
        }
    }
}
