//! The traceback from the best cell to the start of the local alignment.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::fill::{Pass, ScoreWidth};
use super::geometry::{
    BestCell, DELETION, INSERTION, MATCH, START, WaveCarry, traceback_positive, traceback_previous,
    wave_range,
};
use super::{AlignmentConfig, AlignmentWorkspace};
use crate::align::{AlignmentError, EditOperation};

impl AlignmentWorkspace {
    /// Traces back over a pass whose traceback cells are all resident. Returns the start cell.
    pub(super) fn traceback_compact(
        &mut self,
        query: &[u8],
        contig: &[u8],
        config: AlignmentConfig,
        best: BestCell,
    ) -> Result<(usize, usize), AlignmentError> {
        self.operations.clear();
        let mut position = (best.query_index, best.contig_index, best.state);
        self.trace_operations(query, contig, config, &mut position, 0)?;
        self.operations.reverse();
        Ok((position.0, position.1))
    }

    /// Traces back over a checkpointed pass: fills each chunk again from its checkpoint, from
    /// the chunk of the best cell down to the start. Callers must enable AVX2 when VECTOR is true.
    #[inline(always)]
    pub(super) unsafe fn traceback_chunked<const VECTOR: bool>(
        &mut self,
        pass: Pass<'_>,
        width: ScoreWidth,
        best: BestCell,
    ) -> Result<(usize, usize), AlignmentError> {
        let (query, contig, config, stride) = (pass.query, pass.contig, pass.config, pass.stride);
        self.operations.clear();
        let mut position = (best.query_index, best.contig_index, best.state);
        let mut chunk = self
            .chunk_waves
            .partition_point(|&(first, _)| first <= position.0 + position.1)
            .checked_sub(1)
            .ok_or(AlignmentError::InvalidTraceback)?;
        loop {
            let (first, last) = self.chunk_waves[chunk];
            let stop = last.min(position.0 + position.1);
            let mut carry: WaveCarry = self.chunk_carries[chunk];
            let saved = &self.checkpoints[chunk * 6 * stride..(chunk + 1) * 6 * stride];
            for (wave, saved) in self.waves[..2]
                .iter_mut()
                .zip(saved.chunks_exact(3 * stride))
            {
                let scores: &mut [i32] = bytemuck::cast_slice_mut(&mut wave[..]);
                scores[..3 * stride].copy_from_slice(saved);
            }
            self.layout_chunk(query.len(), contig.len(), config, first, stop)?;
            let mut ignored = BestCell::default();
            // The first fill of these waves held every value, so a refill cannot stop.
            // SAFETY: forwarded from this function's contract.
            unsafe {
                self.fill_waves_selected::<VECTOR, true, false>(
                    pass,
                    width,
                    (first, stop),
                    &mut carry,
                    &mut ignored,
                )
            }
            .map_err(|_| AlignmentError::InvalidTraceback)?;
            if self.trace_operations(query, contig, config, &mut position, first)? {
                break;
            }
            chunk = chunk
                .checked_sub(1)
                .ok_or(AlignmentError::InvalidTraceback)?;
        }
        self.operations.reverse();
        Ok((position.0, position.1))
    }

    /// Appends reversed operations until the local start, or returns false before a wave
    /// below first_wave whose traceback cells are not resident.
    fn trace_operations(
        &mut self,
        query: &[u8],
        contig: &[u8],
        config: AlignmentConfig,
        position: &mut (usize, usize, u8),
        first_wave: usize,
    ) -> Result<bool, AlignmentError> {
        let (mut query_index, mut contig_index, mut state) = *position;
        let complete = loop {
            if query_index == 0 && contig_index == 0 {
                break true;
            }
            if query_index + contig_index < first_wave {
                break false;
            }
            let wave = query_index + contig_index;
            let index = wave_range(query.len(), contig.len(), config, wave)
                .filter(|range| range.contains(query_index))
                .and_then(|range| {
                    let base = *self.wave_offsets.get(wave - first_wave)?;
                    (base + query_index - range.start < *self.wave_offsets.last()?)
                        .then_some(base + query_index - range.start)
                })
                .ok_or(AlignmentError::TracebackOutsideBand)?;
            let traceback = self.compact_cells[index];
            if state > INSERTION || !traceback_positive(traceback, state) {
                break true;
            }
            let previous = traceback_previous(traceback, state);
            match state {
                MATCH if query_index > 0 && contig_index > 0 => {
                    self.operations
                        .push(if query[query_index - 1] == contig[contig_index - 1] {
                            EditOperation::Match
                        } else {
                            EditOperation::Substitution
                        });
                    query_index -= 1;
                    contig_index -= 1;
                }
                DELETION if contig_index > 0 => {
                    self.operations.push(EditOperation::Deletion);
                    contig_index -= 1;
                }
                INSERTION if query_index > 0 => {
                    self.operations.push(EditOperation::Insertion);
                    query_index -= 1;
                }
                _ => return Err(AlignmentError::InvalidTraceback),
            }
            if previous == START {
                break true;
            }
            state = previous;
        };
        *position = (query_index, contig_index, state);
        Ok(complete)
    }
}
