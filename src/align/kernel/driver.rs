//! The pass driver: picks the kernel, plans how a pass keeps its traceback (all resident, or in
//! chunks of at most `max_cells` cells refilled from checkpoints), and runs the pass.

use super::fill::{LocalKernel, Pass};
use super::geometry::{BestCell, WaveCarry, WaveRange, wave_bounds, wave_range};
#[cfg(test)]
use super::rebased::REBASE_REACHED;
use super::rebased::rebased_local_scores;
use super::waves::{WAVE_PAD, prepare_wave};
use super::{AlignmentConfig, AlignmentWorkspace, RawAlignment};
use crate::align::AlignmentError;

/// Vector passes above `max_cells` fill all waves once without traceback, then fill each chunk
/// of this many cells again for the traceback. That is faster than filling the chunks with
/// traceback twice; passes within `max_cells` keep their traceback resident in one pass.
const FINE_CHECKPOINT_CELLS: usize = 65_536;

#[cfg(test)]
thread_local! {
    /// Keeps the rebased kernel off, so tests can compare it with the i32 kernel.
    pub(super) static WIDE_LOCAL_KERNEL: std::cell::Cell<bool> =
        const { std::cell::Cell::new(false) };
}

/// Rebased i16 offsets for every vector pass the configuration allows, i32 scores otherwise.
fn select_kernel<const VECTOR: bool>(config: AlignmentConfig, allow_rebased: bool) -> LocalKernel {
    #[cfg(test)]
    let wide = WIDE_LOCAL_KERNEL.get();
    #[cfg(not(test))]
    let wide = false;
    if VECTOR && !wide && allow_rebased && rebased_local_scores(config) {
        LocalKernel::Rebased
    } else {
        LocalKernel::Wide
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
        self.check_interrupt()?;
        // SAFETY: the scalar instantiation never calls AVX2 code.
        unsafe { self.align_raw_waves::<false>(query, contig, config) }
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
        unsafe { self.align_raw_waves::<true>(query, contig, config) }
    }

    /// Aligns with the wave kernels. A rebased pass that meets an unrepresentable value stops, and
    /// the whole window reruns on the i32 kernel, so the result is always that of the i32 kernel.
    /// Callers must enable AVX2 when VECTOR is true.
    #[inline(always)]
    unsafe fn align_raw_waves<const VECTOR: bool>(
        &mut self,
        query: &[u8],
        contig: &[u8],
        config: AlignmentConfig,
    ) -> Result<RawAlignment, AlignmentError> {
        self.rebase_violation = false;
        // SAFETY: forwarded from this function's contract.
        let result = unsafe { self.align_waves_with::<VECTOR>(query, contig, config, true) };
        if !self.rebase_violation {
            return result;
        }
        self.rebase_violation = false;
        #[cfg(test)]
        REBASE_REACHED.set({
            let (base, up, down, fallbacks) = REBASE_REACHED.get();
            (base, up, down, fallbacks + 1)
        });
        // SAFETY: forwarded from this function's contract.
        unsafe { self.align_waves_with::<VECTOR>(query, contig, config, false) }
    }

    #[inline(always)]
    unsafe fn align_waves_with<const VECTOR: bool>(
        &mut self,
        query: &[u8],
        contig: &[u8],
        config: AlignmentConfig,
        allow_rebased: bool,
    ) -> Result<RawAlignment, AlignmentError> {
        let kernel = select_kernel::<VECTOR>(config, allow_rebased);
        if kernel == LocalKernel::Wide {
            // SAFETY: forwarded from this function's contract.
            return unsafe { self.align_waves_kernel::<VECTOR>(query, contig, config, kernel) };
        }
        // The rebased kernel compares bytes directly, so it reads lowercase copies of both.
        let mut folded = std::mem::take(&mut self.folded);
        folded.clear();
        folded.extend(query.iter().chain(contig).map(u8::to_ascii_lowercase));
        let (query, contig) = folded.split_at(query.len());
        // SAFETY: forwarded from this function's contract.
        let result = unsafe { self.align_waves_kernel::<VECTOR>(query, contig, config, kernel) };
        self.folded = folded;
        result
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

    #[inline(always)]
    unsafe fn align_waves_kernel<const VECTOR: bool>(
        &mut self,
        query: &[u8],
        contig: &[u8],
        config: AlignmentConfig,
        kernel: LocalKernel,
    ) -> Result<RawAlignment, AlignmentError> {
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
            unsafe { self.fill_checkpointed::<VECTOR>(pass, kernel, &layout, &mut best)? };
        } else {
            self.layout_chunk(query.len(), contig.len(), config, 0, last_wave)?;
            // SAFETY: forwarded from this function's contract.
            unsafe {
                self.fill_waves_selected::<VECTOR, true, true>(
                    pass,
                    kernel,
                    (0, last_wave),
                    &mut carry,
                    &mut best,
                );
            }
        }
        self.check_interrupt()?;
        if self.rebase_violation || best.score <= 0 {
            return Err(AlignmentError::NoAlignment);
        }
        let (query_start, contig_start) = if layout.checkpointed {
            // SAFETY: forwarded from this function's contract.
            unsafe { self.traceback_chunked::<VECTOR>(pass, kernel, best)? }
        } else {
            self.traceback_compact(query, contig, config, best)?
        };
        Ok(RawAlignment {
            score: best.score,
            query_interval: (query_start, best.query_index),
            contig_interval: (contig_start, best.contig_index),
            runs: super::traceback::runs_from_operations(&self.operations)?,
        })
    }

    /// Fills the waves in chunks of at most `checkpoint_limit` cells, saving the two score
    /// waves before each chunk.
    #[inline(always)]
    unsafe fn fill_checkpointed<const VECTOR: bool>(
        &mut self,
        pass: Pass<'_>,
        kernel: LocalKernel,
        layout: &PassLayout,
        best: &mut BestCell,
    ) -> Result<(), AlignmentError> {
        let (query_len, contig_len) = (pass.query.len(), pass.contig.len());
        let config = pass.config;
        let last_wave = query_len + contig_len;
        let (active_first, active_last) =
            wave_bounds(query_len, contig_len, config).ok_or(AlignmentError::BandExcludesInput)?;
        self.chunk_waves.clear();
        self.chunk_carries.clear();
        self.checkpoints.clear();
        let width =
            |wave| wave_range(query_len, contig_len, config, wave).map_or(0, WaveRange::width);
        let mut carry = WaveCarry::default();
        let mut first = 0;
        while first <= last_wave {
            self.check_interrupt()?;
            let mut last = first;
            let mut cells = width(first);
            if last < active_first {
                last = active_first - 1;
            }
            while last < last_wave && cells + width(last + 1) <= layout.checkpoint_limit {
                last += 1;
                cells += width(last);
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
                        kernel,
                        (first, last),
                        &mut carry,
                        best,
                    );
                } else {
                    self.layout_chunk(query_len, contig_len, config, first, last)?;
                    self.fill_waves_selected::<VECTOR, true, true>(
                        pass,
                        kernel,
                        (first, last),
                        &mut carry,
                        best,
                    );
                }
            }
            if self.rebase_violation {
                return Ok(());
            }
            first = last + 1;
        }
        Ok(())
    }
}
