//! The banded affine-gap local alignment kernel.
//!
//! The band is filled in anti-diagonal waves. On hosts with AVX2, waves run in vector blocks of
//! sixteen i16 offsets from a moving base (`rebased`), or of eight i32 scores when the offsets
//! cannot represent a pass; other hosts run the same waves cell by cell. All paths give the same
//! alignment. Traceback cells take two bytes; long passes keep only a chunk of them and fill
//! earlier chunks again from checkpoints during the traceback.

#[cfg(target_arch = "x86_64")]
mod avx2;
mod driver;
mod endpoints;
mod fill;
mod geometry;
mod rebased;
#[cfg(target_arch = "x86_64")]
mod rebased_avx2;
mod traceback;
mod waves;

use super::{Alignment, AlignmentError, EditOperation, EditRun, Interval, Strand, complement};
use geometry::WaveCarry;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use waves::WaveChunk;

/// The default limit of resident traceback cells of one pass.
pub const DEFAULT_MAX_CELLS: usize = 4_000_000;

/// Scores and band of one alignment. The band holds the cells whose diagonal (contig position
/// minus query position) lies within `band_width` of `diagonal_offset`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AlignmentConfig {
    pub match_score: i32,
    pub mismatch_score: i32,
    /// Added once for each gap, on top of one `gap_extend_score` per gap base.
    pub gap_open_score: i32,
    pub gap_extend_score: i32,
    pub band_width: u32,
    pub diagonal_offset: i64,
    /// Traceback cells a pass keeps resident; longer passes recompute chunks from checkpoints.
    pub max_cells: usize,
}

impl Default for AlignmentConfig {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_score: -3,
            gap_open_score: -5,
            gap_extend_score: -1,
            band_width: 128,
            diagonal_offset: 0,
            max_cells: DEFAULT_MAX_CELLS,
        }
    }
}

impl AlignmentConfig {
    /// Checks that a match scores positive, the other terms at most zero, and cells are allowed.
    pub fn validate(self) -> Result<(), AlignmentError> {
        if self.match_score <= 0
            || self.mismatch_score > 0
            || self.gap_open_score > 0
            || self.gap_extend_score > 0
            || self.max_cells == 0
        {
            return Err(AlignmentError::InvalidConfig);
        }
        Ok(())
    }

    /// The score of the first base of a gap.
    pub fn gap_open(self) -> i32 {
        self.gap_open_score.saturating_add(self.gap_extend_score)
    }
}

/// A query window and a contig window, with the positions where they start. Results are in
/// those positions: `query_offset + i` for query window byte `i`, likewise for the contig.
#[derive(Clone, Copy, Debug)]
pub struct AlignmentInput<'a> {
    pub query: &'a [u8],
    pub query_offset: u64,
    /// The contig window on the forward strand; the reverse strand reads its reverse complement.
    pub contig: &'a [u8],
    pub contig_offset: u64,
}

/// The kernels a workspace runs.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
enum Backend {
    /// AVX2 when the host has it, scalar otherwise.
    #[default]
    Detected,
    Scalar,
}

impl Backend {
    /// Whether calls run the AVX2 kernels.
    fn vector(self) -> bool {
        match self {
            #[cfg(target_arch = "x86_64")]
            Self::Detected => is_x86_feature_detected!("avx2"),
            #[cfg(not(target_arch = "x86_64"))]
            Self::Detected => false,
            Self::Scalar => false,
        }
    }
}

/// Reusable buffers of the kernel. One workspace serves one thread; buffers keep their capacity
/// between alignments.
#[derive(Debug, Default)]
pub struct AlignmentWorkspace {
    backend: Backend,
    interrupt: Option<Arc<AtomicBool>>,
    /// Traceback cells of the resident waves, two bytes each.
    compact_cells: Vec<u16>,
    /// The score planes of three consecutive waves.
    waves: [Vec<WaveChunk>; 3],
    /// First and last wave of each chunk of a checkpointed pass.
    chunk_waves: Vec<(usize, usize)>,
    chunk_carries: Vec<WaveCarry>,
    /// The two score waves before each chunk.
    checkpoints: Vec<i32>,
    /// Set by a rebased pass that met a value it cannot represent; the pass reruns on i32.
    rebase_violation: bool,
    /// Traceback offset of each laid-out wave, relative to the layout's first wave, plus the
    /// total; one wave's cells are stored contiguously in query-row order.
    wave_offsets: Vec<usize>,
    operations: Vec<EditOperation>,
    /// The reverse complement of the contig window for the reverse strand.
    reverse: Vec<u8>,
    /// The query then the contig in ASCII lowercase, for the rebased kernel's byte compare.
    folded: Vec<u8>,
    endpoint: endpoints::EndpointWorkspace,
}

/// A local alignment in the positions of the oriented query and contig slices.
struct RawAlignment {
    score: i32,
    query_interval: (usize, usize),
    contig_interval: (usize, usize),
    runs: Vec<EditRun>,
}

impl AlignmentWorkspace {
    /// A workspace that runs the scalar kernels even on hosts with AVX2.
    pub fn scalar() -> Self {
        Self {
            backend: Backend::Scalar,
            ..Self::default()
        }
    }

    /// Stops running and later alignments with `NoAlignment` once `flag` is set.
    pub fn set_interrupt(&mut self, flag: Arc<AtomicBool>) {
        self.interrupt = Some(flag);
    }

    fn check_interrupt(&self) -> Result<(), AlignmentError> {
        if self
            .interrupt
            .as_ref()
            .is_some_and(|flag| flag.load(Ordering::Relaxed))
        {
            return Err(AlignmentError::NoAlignment);
        }
        Ok(())
    }

    /// The best local alignment of the query window to the contig window on `strand`, within the
    /// band of `config`. On the reverse strand the band's diagonals count positions in the
    /// reverse complement of the contig window.
    pub fn align(
        &mut self,
        input: AlignmentInput<'_>,
        strand: Strand,
        config: AlignmentConfig,
    ) -> Result<Alignment, AlignmentError> {
        let raw = match strand {
            Strand::Forward => self.align_raw(input.query, input.contig, config)?,
            Strand::Reverse => {
                let mut reverse = std::mem::take(&mut self.reverse);
                reverse_complement_into(input.contig, &mut reverse);
                let result = self.align_raw(input.query, &reverse, config);
                self.reverse = reverse;
                result?
            }
        };
        finish(raw, input, strand)
    }
}

/// Writes the reverse complement of `sequence` into `output`.
fn reverse_complement_into(sequence: &[u8], output: &mut Vec<u8>) {
    output.clear();
    output.reserve(sequence.len());
    output.extend(sequence.iter().rev().map(|base| complement(*base)));
}

/// The alignment of `raw` in the positions of `input`: contig positions on the forward strand.
fn finish(
    raw: RawAlignment,
    input: AlignmentInput<'_>,
    strand: Strand,
) -> Result<Alignment, AlignmentError> {
    let (start, end) = raw.contig_interval;
    let (start, end) = match strand {
        Strand::Forward => (start, end),
        Strand::Reverse => (input.contig.len() - end, input.contig.len() - start),
    };
    let shift = |offset: u64, position: usize| {
        u64::try_from(position)
            .ok()
            .and_then(|position| offset.checked_add(position))
            .ok_or(AlignmentError::LengthOverflow)
    };
    let alignment = Alignment {
        score: raw.score,
        strand,
        query_interval: Interval::new(
            shift(input.query_offset, raw.query_interval.0)?,
            shift(input.query_offset, raw.query_interval.1)?,
        )?,
        contig_interval: Interval::new(
            shift(input.contig_offset, start)?,
            shift(input.contig_offset, end)?,
        )?,
        runs: raw.runs,
    };
    alignment.validate()?;
    Ok(alignment)
}
