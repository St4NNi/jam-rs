//! The shape of the band in anti-diagonal waves, and the cell states of the local recurrence.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

// Wave `w` holds the cells with `row + column = w` (row: query, column: contig); a cell lies in
// the band when `|column - row - diagonal_offset| <= band_width`.

use super::AlignmentConfig;

/// The state of a cell that pairs a query base with a contig base.
pub(super) const MATCH: u8 = 0;
/// The state of a cell that ends in contig bases only (CIGAR `D`).
pub(super) const DELETION: u8 = 1;
/// The state of a cell that ends in query bases only (CIGAR `I`).
pub(super) const INSERTION: u8 = 2;
/// The traceback mark of a cell where a local alignment starts.
pub(super) const START: u8 = 3;
/// The state of no candidate; it loses every tie.
pub(super) const UNREACHABLE: u8 = 4;

/// Score wave ranges carried into the next wave, saved at each resident chunk start.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct WaveCarry {
    pub(super) older: Option<WaveRange>,
    pub(super) previous: Option<WaveRange>,
    /// Base of the rebased i16 planes: a stored offset `o` is the score `o + base`. Zero for the
    /// other element types.
    pub(super) base: i32,
}

/// The inclusive rows of one wave inside the band.
#[derive(Clone, Copy, Debug)]
pub(super) struct WaveRange {
    pub(super) start: usize,
    pub(super) end: usize,
}

impl WaveRange {
    /// The number of rows.
    pub(super) fn width(self) -> usize {
        self.end - self.start + 1
    }

    pub(super) fn contains(self, row: usize) -> bool {
        row >= self.start && row <= self.end
    }
}

/// The first and last wave with a cell in the band, if any.
pub(super) fn wave_bounds(
    query_len: usize,
    contig_len: usize,
    config: AlignmentConfig,
) -> Option<(usize, usize)> {
    // For 0<=i<=m, 0<=j<=n and d-b<=j-i<=d+b, extrema of i+j are these
    // clipped bounds. Interior waves can still be empty, notably for a zero-width band.
    let m = i128::try_from(query_len).ok()?;
    let n = i128::try_from(contig_len).ok()?;
    let d = i128::from(config.diagonal_offset);
    let b = i128::from(config.band_width);
    let low = d.checked_sub(b)?;
    let high = d.checked_add(b)?;
    let first = 0.max(low).max(high.checked_neg()?);
    let last = m
        .checked_add(n)?
        .min(m.checked_mul(2)?.checked_add(high)?)
        .min(n.checked_mul(2)?.checked_sub(low)?);
    (first <= last).then_some((usize::try_from(first).ok()?, usize::try_from(last).ok()?))
}

/// The rows of `wave` inside the band, if any.
#[inline]
pub(super) fn wave_range(
    query_len: usize,
    contig_len: usize,
    config: AlignmentConfig,
    wave: usize,
) -> Option<WaveRange> {
    // Below this every intermediate value fits i64. No slice in memory comes near it, and a
    // diagonal this far off holds no cell of a slice that does fit.
    const SMALL: u64 = 1 << 60;
    let diagonal = config.diagonal_offset;
    if (query_len as u64) >= SMALL
        || (contig_len as u64) >= SMALL
        || (wave as u64) >= SMALL
        || diagonal.unsigned_abs() >= SMALL
    {
        return None;
    }
    let (query_len, contig_len, wave) = (query_len as i64, contig_len as i64, wave as i64);
    let band = i64::from(config.band_width);
    // ceil(x / 2) is -floor(-x / 2), and an arithmetic shift floors.
    let lower = 0
        .max(wave - contig_len)
        .max(-((diagonal + band - wave) >> 1));
    let upper = query_len.min(wave).min((wave - (diagonal - band)) >> 1);
    (lower <= upper).then_some(WaveRange {
        start: lower as usize,
        end: upper as usize,
    })
}

/// The three state scores of one cell and the state each came from.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct Cell {
    pub(super) scores: [i32; 3],
    pub(super) previous: [u8; 3],
}

impl Cell {
    pub(super) fn best_score(self) -> (i32, u8) {
        best_score(self.scores)
    }
}

/// The best cell so far: the highest score, then the smallest (row, column, state), so the
/// result does not depend on the order the cells are visited in.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct BestCell {
    pub(super) score: i32,
    pub(super) query_index: usize,
    pub(super) contig_index: usize,
    pub(super) state: u8,
}

impl BestCell {
    /// Offers lane `lane` of a vector block that starts at `row` of `wave`. Lanes run in query
    /// order, so the first lane with the block maximum is the earliest cell; later ones lose.
    #[inline(always)]
    pub(super) fn consider_lane(&mut self, wave: usize, row: usize, lane: usize, scores: [i32; 3]) {
        let query_index = row + lane;
        let cell = Cell {
            scores,
            previous: [0; 3],
        };
        self.consider(query_index, wave - query_index, cell);
    }

    pub(super) fn consider(&mut self, query_index: usize, contig_index: usize, cell: Cell) {
        let (score, state) = cell.best_score();
        if score > self.score
            || (score == self.score
                && score > 0
                && (query_index, contig_index, state)
                    < (self.query_index, self.contig_index, self.state))
        {
            *self = Self {
                score,
                query_index,
                contig_index,
                state,
            };
        }
    }
}

/// The highest positive candidate, ties going to the smaller state; `(0, UNREACHABLE)` if none
/// is positive.
pub(super) fn choose<const N: usize>(values: [(i32, u8); N]) -> (i32, u8) {
    values.into_iter().fold((0, UNREACHABLE), |best, value| {
        if value.0 > best.0 || (value.0 == best.0 && value.1 < best.1) {
            value
        } else {
            best
        }
    })
}

/// The best of the three state scores of a cell.
pub(super) fn best_score(scores: [i32; 3]) -> (i32, u8) {
    choose([
        (scores[MATCH as usize], MATCH),
        (scores[DELETION as usize], DELETION),
        (scores[INSERTION as usize], INSERTION),
    ])
}

/// The traceback of one cell in nine bits: bit `s` marks a positive state `s`, and bits
/// `3 + 2s` and `4 + 2s` hold the state it came from.
pub(super) fn encode_traceback(cell: Cell) -> u16 {
    let mut traceback = 0u16;
    for state in 0..=INSERTION {
        let index = state as usize;
        if cell.scores[index] > 0 {
            debug_assert!(cell.previous[index] <= START);
            traceback |= 1 << state;
            traceback |= u16::from(cell.previous[index]) << (3 + 2 * state);
        }
    }
    traceback
}

/// Whether `state` of a traceback cell scored positive.
pub(super) fn traceback_positive(traceback: u16, state: u8) -> bool {
    traceback & (1 << state) != 0
}

/// The state that `state` of a traceback cell came from.
pub(super) fn traceback_previous(traceback: u16, state: u8) -> u8 {
    ((traceback >> (3 + 2 * state)) & 3) as u8
}
