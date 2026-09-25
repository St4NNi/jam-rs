//! The pass driver: picks the score width, plans how a pass keeps its traceback, and runs it.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

// A pass keeps its traceback resident, or in chunks of at most `max_cells` cells that are
// filled again from checkpoints.

use super::fill::{Pass, RebaseViolation, ScoreWidth};
use super::geometry::{BestCell, WaveCarry, WaveRange, wave_bounds, wave_range};
use super::rebased::rebased_local_scores;
#[cfg(test)]
use super::rebased::{TEST_HOOKS, update_hooks};
use super::waves::{WAVE_PAD, prepare_wave};
use super::{AlignmentConfig, AlignmentWorkspace, RawAlignment};
use crate::align::AlignmentError;

/// The folded byte of an unknown query base and of an unknown contig base. They differ, so an
/// unknown base never matches, not even the same code on the other side.
const QUERY_UNKNOWN: u8 = 0;
const CONTIG_UNKNOWN: u8 = 1;

/// Writes the query then the contig into `output`: A, C, G and T in lowercase, every other byte
/// as the unknown code of its side.
fn fold_bases(query: &[u8], contig: &[u8], output: &mut Vec<u8>) {
    let fold = |unknown: u8| {
        move |&base: &u8| match base.to_ascii_lowercase() {
            known @ (b'a' | b'c' | b'g' | b't') => known,
            _ => unknown,
        }
    };
    output.clear();
    output.extend(query.iter().map(fold(QUERY_UNKNOWN)));
    output.extend(contig.iter().map(fold(CONTIG_UNKNOWN)));
}

/// Vector passes above `max_cells` fill all waves once without traceback, then fill each chunk
/// of this many cells again for the traceback. That is faster than filling the chunks with
/// traceback twice; passes within `max_cells` keep their traceback resident in one pass.
const FINE_CHECKPOINT_CELLS: usize = 65_536;

/// Why a pass ended without an alignment: its i16 offsets met a value they cannot hold, or it
/// failed.
enum PassStop {
    Rebase,
    Failed(AlignmentError),
}

impl From<AlignmentError> for PassStop {
    fn from(error: AlignmentError) -> Self {
        Self::Failed(error)
    }
}

impl From<RebaseViolation> for PassStop {
    fn from(_: RebaseViolation) -> Self {
        Self::Rebase
    }
}

/// How a pass keeps its traceback: the waves of every chunk, and whether the chunks are filled
/// first without traceback.
struct PassLayout {
    stride: usize,
    checkpointed: bool,
    score_first: bool,
    checkpoint_limit: usize,
}

impl AlignmentWorkspace {
    /// The best local alignment of `query` and `contig` (already on the right strand), in the
    /// positions of both slices.
    pub(super) fn align_raw(
        &mut self,
        query: &[u8],
        contig: &[u8],
        config: AlignmentConfig,
    ) -> Result<RawAlignment, AlignmentError> {
        #[cfg(target_arch = "x86_64")]
        if self.backend.vector() {
            // SAFETY: `vector` is true only on hosts with AVX2.
            return unsafe { self.align_raw_avx2(query, contig, config) };
        }
        // SAFETY: the scalar instantiation never calls AVX2 code.
        unsafe { self.align_waves::<false>(query, contig, config) }
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    unsafe fn align_raw_avx2(
        &mut self,
        query: &[u8],
        contig: &[u8],
        config: AlignmentConfig,
    ) -> Result<RawAlignment, AlignmentError> {
        // SAFETY: this function enables AVX2 for the inlined wave driver.
        unsafe { self.align_waves::<true>(query, contig, config) }
    }

    /// Aligns with i16 offsets when a vector pass allows them. A pass that meets a value they
    /// cannot hold reruns with i32 scores, so the result is always that of the i32 kernel.
    /// Callers must enable AVX2 when VECTOR is true.
    #[inline(always)]
    unsafe fn align_waves<const VECTOR: bool>(
        &mut self,
        query: &[u8],
        contig: &[u8],
        config: AlignmentConfig,
    ) -> Result<RawAlignment, AlignmentError> {
        #[cfg(test)]
        let force_i32 = TEST_HOOKS.get().force_i32;
        #[cfg(not(test))]
        let force_i32 = false;
        // Every pass compares bytes directly, so it reads folded copies of both.
        let mut folded = std::mem::take(&mut self.folded);
        fold_bases(query, contig, &mut folded);
        let (query, contig) = folded.split_at(query.len());
        let mut result = Err(PassStop::Rebase);
        if VECTOR && !force_i32 && rebased_local_scores(config) {
            // SAFETY: forwarded from this function's contract.
            result = unsafe { self.run_pass::<VECTOR>(query, contig, config, ScoreWidth::I16) };
            #[cfg(test)]
            if matches!(result, Err(PassStop::Rebase)) {
                update_hooks(|hooks| hooks.reruns += 1);
            }
        }
        if matches!(result, Err(PassStop::Rebase)) {
            // SAFETY: forwarded from this function's contract.
            result = unsafe { self.run_pass::<VECTOR>(query, contig, config, ScoreWidth::I32) };
        }
        self.folded = folded;
        match result {
            Ok(raw) => Ok(raw),
            Err(PassStop::Failed(error)) => Err(error),
            // i32 scores hold every value, so this pass cannot stop for a rebase.
            Err(PassStop::Rebase) => Err(AlignmentError::InvalidTraceback),
        }
    }

    /// Decides how a pass keeps its traceback and sizes the wave buffers.
    fn plan_pass<const VECTOR: bool>(
        &mut self,
        query: &[u8],
        contig: &[u8],
        config: AlignmentConfig,
    ) -> Result<PassLayout, AlignmentError> {
        config.validate()?;
        if query.is_empty() {
            return Err(AlignmentError::EmptyQuery);
        }
        if contig.is_empty() {
            return Err(AlignmentError::EmptyContig);
        }
        let (active_first, active_last) = wave_bounds(query.len(), contig.len(), config)
            .ok_or(AlignmentError::BandExcludesInput)?;
        let mut max_wave_width = 0usize;
        let mut total_cells = 0usize;
        for wave in active_first..=active_last {
            if let Some(range) = wave_range(query.len(), contig.len(), config, wave) {
                max_wave_width = max_wave_width.max(range.width());
                total_cells = total_cells
                    .checked_add(range.width())
                    .ok_or(AlignmentError::LengthOverflow)?;
            }
        }
        if total_cells == 0 {
            return Err(AlignmentError::BandExcludesInput);
        }
        // A pass above max_cells keeps at most max_cells traceback cells resident and recomputes
        // earlier waves from checkpoints; every wave must fit on its own.
        let chunked = total_cells > config.max_cells;
        if chunked && max_wave_width > config.max_cells {
            return Err(AlignmentError::MatrixTooLarge {
                cells: max_wave_width,
                max_cells: config.max_cells,
            });
        }
        let score_first = VECTOR
            && chunked
            && config.max_cells > FINE_CHECKPOINT_CELLS
            && max_wave_width <= FINE_CHECKPOINT_CELLS;
        let stride = max_wave_width
            .checked_add(2 * WAVE_PAD)
            .ok_or(AlignmentError::LengthOverflow)?
            .next_multiple_of(32);
        let wave_scores = stride
            .checked_mul(3)
            .ok_or(AlignmentError::LengthOverflow)?;
        for wave in &mut self.waves {
            prepare_wave(wave, wave_scores);
        }
        Ok(PassLayout {
            stride,
            checkpointed: chunked || score_first,
            score_first,
            checkpoint_limit: if score_first {
                FINE_CHECKPOINT_CELLS
            } else {
                config.max_cells
            },
        })
    }

    /// Runs one pass with scores of `width` and traces the best cell back. Callers must enable
    /// AVX2 when VECTOR is true.
    #[inline(always)]
    unsafe fn run_pass<const VECTOR: bool>(
        &mut self,
        query: &[u8],
        contig: &[u8],
        config: AlignmentConfig,
        width: ScoreWidth,
    ) -> Result<RawAlignment, PassStop> {
        let layout = self.plan_pass::<VECTOR>(query, contig, config)?;
        let last_wave = query.len() + contig.len();
        let pass = Pass {
            query,
            contig,
            config,
            stride: layout.stride,
        };
        let mut best = BestCell::default();
        let mut carry = WaveCarry::default();
        if layout.checkpointed {
            // SAFETY: forwarded from this function's contract.
            unsafe { self.fill_checkpointed::<VECTOR>(pass, width, &layout, &mut best)? };
        } else {
            self.layout_chunk(query.len(), contig.len(), config, 0, last_wave)?;
            // SAFETY: forwarded from this function's contract.
            unsafe {
                self.fill_waves_selected::<VECTOR, true, true>(
                    pass,
                    width,
                    (0, last_wave),
                    &mut carry,
                    &mut best,
                )?;
            }
        }
        if best.score <= 0 {
            return Err(AlignmentError::NoAlignment.into());
        }
        let (query_start, contig_start) = if layout.checkpointed {
            // SAFETY: forwarded from this function's contract.
            unsafe { self.traceback_chunked::<VECTOR>(pass, width, best)? }
        } else {
            self.traceback_compact(query, contig, config, best)?
        };
        Ok(RawAlignment {
            score: best.score,
            query_interval: (query_start, best.query_index),
            contig_interval: (contig_start, best.contig_index),
            runs: crate::align::cigar::runs_of(&self.operations)?,
        })
    }

    /// Fills the waves in chunks of at most `checkpoint_limit` cells, saving the two score
    /// waves before each chunk.
    #[inline(always)]
    unsafe fn fill_checkpointed<const VECTOR: bool>(
        &mut self,
        pass: Pass<'_>,
        width: ScoreWidth,
        layout: &PassLayout,
        best: &mut BestCell,
    ) -> Result<(), PassStop> {
        let (query_len, contig_len) = (pass.query.len(), pass.contig.len());
        let config = pass.config;
        let last_wave = query_len + contig_len;
        let (active_first, active_last) =
            wave_bounds(query_len, contig_len, config).ok_or(AlignmentError::BandExcludesInput)?;
        self.chunk_waves.clear();
        self.chunk_carries.clear();
        self.checkpoints.clear();
        let wave_cells =
            |wave| wave_range(query_len, contig_len, config, wave).map_or(0, WaveRange::width);
        let mut carry = WaveCarry::default();
        let mut first = 0;
        while first <= last_wave {
            let mut last = first;
            let mut cells = wave_cells(first);
            if last < active_first {
                last = active_first - 1;
            }
            while last < last_wave && cells + wave_cells(last + 1) <= layout.checkpoint_limit {
                last += 1;
                cells += wave_cells(last);
                if last >= active_last {
                    last = last_wave;
                    break;
                }
            }
            self.save_checkpoint(carry, layout.stride);
            self.chunk_waves.push((first, last));
            // SAFETY: forwarded from this function's contract.
            unsafe {
                if layout.score_first {
                    self.fill_waves_selected::<VECTOR, false, true>(
                        pass,
                        width,
                        (first, last),
                        &mut carry,
                        best,
                    )?;
                } else {
                    self.layout_chunk(query_len, contig_len, config, first, last)?;
                    self.fill_waves_selected::<VECTOR, true, true>(
                        pass,
                        width,
                        (first, last),
                        &mut carry,
                        best,
                    )?;
                }
            }
            first = last + 1;
        }
        Ok(())
    }
}
