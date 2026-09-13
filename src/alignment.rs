use serde::ser::SerializeStruct;
use serde::{Deserialize, Serialize};
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
use std::fmt::Write as _;
use std::sync::{Condvar, Mutex};
use thiserror::Error;

const MATCH: u8 = 0;
const INSERTION: u8 = 1;
const DELETION: u8 = 2;
const START: u8 = 3;
const UNREACHABLE: u8 = 4;

pub const DEFAULT_MAX_CELLS: usize = 4_000_000;

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct Interval {
    pub start: u64,
    pub end: u64,
}

impl Interval {
    pub fn new(start: u64, end: u64) -> Result<Self, AlignmentError> {
        if start > end {
            return Err(AlignmentError::ReversedInterval { start, end });
        }
        Ok(Self { start, end })
    }

    pub fn len(self) -> u64 {
        self.end - self.start
    }

    pub fn is_empty(self) -> bool {
        self.start == self.end
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Strand {
    Forward,
    Reverse,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum EditOperation {
    #[serde(rename = "=")]
    Equal,
    #[serde(rename = "X")]
    Substitution,
    #[serde(rename = "I")]
    Insertion,
    #[serde(rename = "D")]
    Deletion,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct EditRun {
    pub operation: EditOperation,
    pub length: u32,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AlignmentConfig {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open_score: i32,
    pub gap_extend_score: i32,
    pub band_width: u32,
    pub diagonal_offset: i64,
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
    fn validate(self) -> Result<(), AlignmentError> {
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
}

#[derive(Clone, Debug, Eq, PartialEq, Deserialize)]
pub struct Alignment {
    pub score: i32,
    pub strand: Strand,
    pub query_interval: Interval,
    pub target_interval: Interval,
    pub matches: u64,
    pub substitutions: u64,
    pub insertions: u64,
    pub deletions: u64,
    pub cigar: String,
    pub edit_script: Vec<EditRun>,
}

impl Serialize for Alignment {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut record = serializer.serialize_struct("Alignment", 11)?;
        record.serialize_field("score", &self.score)?;
        record.serialize_field("strand", &self.strand)?;
        record.serialize_field("query_interval", &self.query_interval)?;
        record.serialize_field("target_interval", &self.target_interval)?;
        record.serialize_field("matches", &self.matches)?;
        record.serialize_field("substitutions", &self.substitutions)?;
        record.serialize_field("insertions", &self.insertions)?;
        record.serialize_field("deletions", &self.deletions)?;
        record.serialize_field("cigar", &self.cigar)?;
        record.serialize_field("edit_script", &self.edit_script)?;
        record.serialize_field("identity", &self.identity())?;
        record.end()
    }
}

impl Alignment {
    pub fn identity(&self) -> f64 {
        let total = self
            .matches
            .saturating_add(self.substitutions)
            .saturating_add(self.insertions)
            .saturating_add(self.deletions);
        if total == 0 {
            0.0
        } else {
            self.matches as f64 / total as f64
        }
    }

    pub fn validate_cigar(&self) -> Result<(), AlignmentError> {
        let runs = parse_cigar(&self.cigar)?;
        if runs != self.edit_script {
            return Err(AlignmentError::CigarMismatch);
        }
        let summary = summarize_runs(&runs)?;
        if summary.query_bases != self.query_interval.len()
            || summary.target_bases != self.target_interval.len()
        {
            return Err(AlignmentError::CigarSpanMismatch);
        }
        if (
            summary.matches,
            summary.substitutions,
            summary.insertions,
            summary.deletions,
        ) != (
            self.matches,
            self.substitutions,
            self.insertions,
            self.deletions,
        ) {
            return Err(AlignmentError::CigarCountMismatch);
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct AlignmentWork {
    pub local_init_cpu_ns: u64,
    pub local_matrix_cpu_ns: u64,
    pub local_trace_cpu_ns: u64,
    pub endpoint_init_cpu_ns: u64,
    pub endpoint_matrix_cpu_ns: u64,
    pub endpoint_trace_cpu_ns: u64,
    pub endpoint_total_cpu_ns: u64,
    pub local_cells: u64,
    pub local_vector8_cells: u64,
    pub local_vector16_cells: u64,
    pub local_scalar_cells: u64,
    pub local_boundary_cells: u64,
    pub local_narrow_eligible_passes: u64,
    pub endpoint_cells: u64,
    pub endpoint_recurrence_cells: u64,
    pub local_passes: u64,
    pub local_chunked_passes: u64,
    pub local_chunks: u64,
    pub local_recomputed_cells: u64,
    pub endpoint_passes: u64,
    pub reverse_bytes: u64,
    pub endpoint_scratch_bytes: u64,
    pub growth_operations: u64,
    pub capacity_bytes: u64,
}

impl AlignmentWork {
    pub(crate) fn add(&mut self, other: Self) {
        self.local_init_cpu_ns += other.local_init_cpu_ns;
        self.local_matrix_cpu_ns += other.local_matrix_cpu_ns;
        self.local_trace_cpu_ns += other.local_trace_cpu_ns;
        self.endpoint_init_cpu_ns += other.endpoint_init_cpu_ns;
        self.endpoint_matrix_cpu_ns += other.endpoint_matrix_cpu_ns;
        self.endpoint_trace_cpu_ns += other.endpoint_trace_cpu_ns;
        self.endpoint_total_cpu_ns += other.endpoint_total_cpu_ns;
        self.local_cells += other.local_cells;
        self.local_vector8_cells += other.local_vector8_cells;
        self.local_vector16_cells += other.local_vector16_cells;
        self.local_scalar_cells += other.local_scalar_cells;
        self.local_boundary_cells += other.local_boundary_cells;
        self.local_narrow_eligible_passes += other.local_narrow_eligible_passes;
        self.endpoint_cells += other.endpoint_cells;
        self.endpoint_recurrence_cells += other.endpoint_recurrence_cells;
        self.local_passes += other.local_passes;
        self.local_chunked_passes += other.local_chunked_passes;
        self.local_chunks += other.local_chunks;
        self.local_recomputed_cells += other.local_recomputed_cells;
        self.endpoint_passes += other.endpoint_passes;
        self.reverse_bytes += other.reverse_bytes;
        self.endpoint_scratch_bytes += other.endpoint_scratch_bytes;
        self.growth_operations += other.growth_operations;
        self.capacity_bytes += other.capacity_bytes;
    }
}

pub(crate) fn observed_cpu(observed: bool) -> Option<u64> {
    observed.then(crate::trace_batch::worker_cpu_ns).flatten()
}

pub(crate) fn elapsed_cpu(start: Option<u64>) -> u64 {
    start
        .and_then(|start| crate::trace_batch::worker_cpu_ns().map(|end| end.saturating_sub(start)))
        .unwrap_or(0)
}

#[derive(Debug, Default)]
pub struct AlignmentWorkspace {
    observed: bool,
    resident_chunks: bool,
    pub(crate) work: AlignmentWork,
    traceback_nanoseconds: u64,
    endpoint_nanoseconds: u64,
    cells: Vec<Cell>,
    #[cfg(target_arch = "x86_64")]
    compact_cells: Vec<u16>,
    #[cfg(target_arch = "x86_64")]
    waves: [Vec<i32>; 3],
    #[cfg(target_arch = "x86_64")]
    chunk_waves: Vec<(usize, usize)>,
    #[cfg(target_arch = "x86_64")]
    chunk_carries: Vec<WaveCarry>,
    #[cfg(target_arch = "x86_64")]
    checkpoints: Vec<i32>,
    /// Traceback offset of each laid-out wave, relative to the layout's first wave, plus the
    /// total; one wave's cells are stored contiguously in query-row order.
    #[cfg(target_arch = "x86_64")]
    wave_offsets: Vec<usize>,
    endpoint_cells: EndpointWorkspace,
    row_offsets: Vec<usize>,
    row_starts: Vec<usize>,
    row_widths: Vec<usize>,
    operations: Vec<EditOperation>,
    reverse: Vec<u8>,
}

pub(crate) const TRACE_ALIGNMENT_BUDGET_BYTES: usize = 2 * 1024 * 1024 * 1024;
// Covers retained workspace buffers and conservative simultaneous task temporaries.
const TRACE_ALIGNMENT_ALLOCATION_COUNT: usize = 26;

#[derive(Clone, Copy, Debug, Error, Eq, PartialEq)]
pub(crate) enum AlignmentAdmissionError {
    #[error("alignment workspace byte bound overflows")]
    ByteOverflow,
    #[error("alignment workspace requires {requested} bytes, exceeding workspace budget {budget}")]
    RequestExceedsBudget { requested: usize, budget: usize },
}

struct AlignmentBytePool {
    budget: usize,
    available: Mutex<usize>,
    wake: Condvar,
}

impl AlignmentBytePool {
    const fn new(budget: usize) -> Self {
        Self {
            budget,
            available: Mutex::new(budget),
            wake: Condvar::new(),
        }
    }

    fn acquire(&self, bytes: usize) -> Result<AlignmentBytePermit<'_>, AlignmentAdmissionError> {
        if bytes > self.budget {
            return Err(AlignmentAdmissionError::RequestExceedsBudget {
                requested: bytes,
                budget: self.budget,
            });
        }
        let mut available = self
            .available
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        while *available < bytes {
            available = self
                .wake
                .wait(available)
                .unwrap_or_else(|poisoned| poisoned.into_inner());
        }
        *available -= bytes;
        Ok(AlignmentBytePermit { pool: self, bytes })
    }

    #[cfg(test)]
    fn available(&self) -> usize {
        *self
            .available
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
    }
}

struct AlignmentBytePermit<'a> {
    pool: &'a AlignmentBytePool,
    bytes: usize,
}

impl Drop for AlignmentBytePermit<'_> {
    fn drop(&mut self) {
        let mut available = self
            .pool
            .available
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        *available += self.bytes;
        debug_assert!(*available <= self.pool.budget);
        self.pool.wake.notify_all();
    }
}

static TRACE_ALIGNMENT_POOL: AlignmentBytePool =
    AlignmentBytePool::new(TRACE_ALIGNMENT_BUDGET_BYTES);

pub(crate) struct TraceAlignmentWorkspace {
    // Accumulated result vectors returned by completed tasks are caller-owned and excluded.
    workspace: AlignmentWorkspace,
    _permit: AlignmentBytePermit<'static>,
}

impl TraceAlignmentWorkspace {
    pub(crate) fn acquire(bytes: usize) -> Result<Self, AlignmentAdmissionError> {
        let permit = TRACE_ALIGNMENT_POOL.acquire(bytes)?;
        Ok(Self {
            workspace: AlignmentWorkspace {
                resident_chunks: true,
                ..AlignmentWorkspace::default()
            },
            _permit: permit,
        })
    }

    pub(crate) fn workspace_mut(&mut self) -> &mut AlignmentWorkspace {
        &mut self.workspace
    }
}

/// Bytes for one reused trace workspace. Each maximum is taken over the tasks it may run, so
/// retained capacities are covered without combining the spans of different tasks.
pub(crate) fn trace_alignment_bytes(
    max_query_bases: usize,
    max_target_bases: usize,
    max_path_bases: usize,
    max_local_cells: usize,
    endpoint_bases: usize,
    config: AlignmentConfig,
) -> Result<usize, AlignmentAdmissionError> {
    let query_rows = max_query_bases
        .checked_add(1)
        .ok_or(AlignmentAdmissionError::ByteOverflow)?;
    let band_columns = usize::try_from(config.band_width)
        .map_err(|_| AlignmentAdmissionError::ByteOverflow)?
        .checked_mul(2)
        .and_then(|value| value.checked_add(1))
        .ok_or(AlignmentAdmissionError::ByteOverflow)?;
    let local_cells = max_local_cells;
    let core_cells = local_cells.min(config.max_cells);
    // Longer traced tasks keep one resident chunk plus two saved score waves per chunk.
    // Greedy whole-wave chunks satisfy chunks <= 2 * ceil(cells / max_cells).
    let chunks = if local_cells > config.max_cells {
        local_cells
            .div_ceil(config.max_cells)
            .checked_mul(2)
            .and_then(|value| value.checked_add(1))
            .ok_or(AlignmentAdmissionError::ByteOverflow)?
    } else {
        0
    };
    let wave_width = band_columns.div_ceil(2).min(query_rows);
    let endpoint_query = max_query_bases.min(endpoint_bases);
    let endpoint_target = max_target_bases.min(endpoint_bases);
    let endpoint_cells = endpoint_query
        .checked_add(1)
        .and_then(|rows| {
            endpoint_target
                .checked_add(1)
                .and_then(|columns| rows.checked_mul(columns))
        })
        .ok_or(AlignmentAdmissionError::ByteOverflow)?
        .min(config.max_cells);
    let path_bases = max_path_bases;
    let endpoint_path_bases = endpoint_query
        .checked_add(endpoint_target)
        .ok_or(AlignmentAdmissionError::ByteOverflow)?;

    let retained = checked_sum(&[
        doubled_vec_bytes::<Cell>(core_cells)?,
        doubled_vec_bytes::<u8>(endpoint_cells)?,
        doubled_vec_bytes::<[i32; 3]>(
            endpoint_target
                .checked_add(1)
                .ok_or(AlignmentAdmissionError::ByteOverflow)?,
        )?
        .checked_mul(2)
        .ok_or(AlignmentAdmissionError::ByteOverflow)?,
        doubled_vec_bytes::<i32>(
            endpoint_query
                .checked_add(9)
                .and_then(|stride| stride.checked_mul(9))
                .ok_or(AlignmentAdmissionError::ByteOverflow)?,
        )?,
        doubled_vec_bytes::<usize>(query_rows)?
            .checked_mul(3)
            .ok_or(AlignmentAdmissionError::ByteOverflow)?,
        doubled_vec_bytes::<usize>(
            query_rows
                .checked_add(max_target_bases)
                .and_then(|waves| waves.checked_add(1))
                .ok_or(AlignmentAdmissionError::ByteOverflow)?,
        )?,
        doubled_vec_bytes::<EditOperation>(path_bases)?,
        doubled_vec_bytes::<u8>(max_target_bases)?,
        doubled_vec_bytes::<i32>(
            wave_width
                .checked_mul(9)
                .ok_or(AlignmentAdmissionError::ByteOverflow)?,
        )?,
        doubled_vec_bytes::<i32>(
            chunks
                .checked_mul(6)
                .and_then(|value| value.checked_mul(wave_width))
                .ok_or(AlignmentAdmissionError::ByteOverflow)?,
        )?,
        doubled_vec_bytes::<[Option<(usize, usize)>; 3]>(chunks)?,
    ])?;
    let run_entries = path_bases
        .checked_mul(3)
        .and_then(|value| {
            endpoint_path_bases
                .checked_mul(4)
                .and_then(|extra| value.checked_add(extra))
        })
        .ok_or(AlignmentAdmissionError::ByteOverflow)?;
    let cigar_bytes = path_bases
        .checked_mul(12)
        .and_then(|value| {
            endpoint_path_bases
                .checked_mul(8)
                .and_then(|extra| value.checked_add(extra))
        })
        .ok_or(AlignmentAdmissionError::ByteOverflow)?;
    let temporary = checked_sum(&[
        doubled_vec_bytes::<u8>(max_query_bases)?,
        doubled_vec_bytes::<u8>(endpoint_query)?,
        doubled_vec_bytes::<u8>(endpoint_target)?,
        doubled_vec_bytes::<EditOperation>(endpoint_path_bases)?,
        doubled_vec_bytes::<EditRun>(run_entries)?,
        cigar_bytes,
    ])?;
    let allocator_margin = TRACE_ALIGNMENT_ALLOCATION_COUNT
        .checked_mul(4096)
        .ok_or(AlignmentAdmissionError::ByteOverflow)?;
    retained
        .checked_add(temporary)
        .and_then(|value| value.checked_add(allocator_margin))
        .ok_or(AlignmentAdmissionError::ByteOverflow)
}

fn doubled_vec_bytes<T>(length: usize) -> Result<usize, AlignmentAdmissionError> {
    length
        .checked_mul(2)
        .and_then(|capacity| capacity.checked_mul(std::mem::size_of::<T>()))
        .ok_or(AlignmentAdmissionError::ByteOverflow)
}

fn checked_sum(values: &[usize]) -> Result<usize, AlignmentAdmissionError> {
    values.iter().try_fold(0usize, |total, value| {
        total
            .checked_add(*value)
            .ok_or(AlignmentAdmissionError::ByteOverflow)
    })
}

impl AlignmentWorkspace {
    pub(crate) fn retained_bytes(&self) -> usize {
        let bytes = self.cells.capacity() * std::mem::size_of::<Cell>()
            + self.endpoint_cells.capacity_bytes()
            + (self.row_offsets.capacity()
                + self.row_starts.capacity()
                + self.row_widths.capacity())
                * std::mem::size_of::<usize>()
            + self.operations.capacity() * std::mem::size_of::<EditOperation>()
            + self.reverse.capacity();
        #[cfg(target_arch = "x86_64")]
        {
            bytes
                + self.compact_cells.capacity() * 2
                + self.waves.iter().map(|w| w.capacity() * 4).sum::<usize>()
                + self.chunk_waves.capacity() * std::mem::size_of::<(usize, usize)>()
                + self.chunk_carries.capacity() * std::mem::size_of::<WaveCarry>()
                + self.checkpoints.capacity() * 4
                + self.wave_offsets.capacity() * std::mem::size_of::<usize>()
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            bytes
        }
    }

    pub(crate) fn enable_timing(&mut self) {
        self.observed = true;
    }

    pub(crate) fn traceback_nanoseconds(&self) -> u64 {
        self.traceback_nanoseconds
    }

    pub(crate) fn endpoint_nanoseconds(&self) -> u64 {
        self.endpoint_nanoseconds
    }
    pub fn capacity_cells(&self) -> usize {
        #[cfg(target_arch = "x86_64")]
        {
            self.cells.capacity().max(self.compact_cells.capacity())
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            self.cells.capacity()
        }
    }

    pub fn align(
        &mut self,
        query: &[u8],
        target: &[u8],
        config: AlignmentConfig,
    ) -> Result<Alignment, AlignmentError> {
        let raw = self.align_raw(query, target, config)?;
        let cpu = observed_cpu(self.observed);
        let started = self.observed.then(std::time::Instant::now);
        let result = finish(raw, Strand::Forward, 0, target.len());
        if let Some(started) = started {
            self.traceback_nanoseconds += started.elapsed().as_nanos() as u64;
            self.work.local_trace_cpu_ns += elapsed_cpu(cpu);
        }
        result
    }

    pub fn align_oriented(
        &mut self,
        query: &[u8],
        target: &[u8],
        target_offset: u64,
        strand: Strand,
        config: AlignmentConfig,
    ) -> Result<Alignment, AlignmentError> {
        let raw = match strand {
            Strand::Forward => self.align_raw(query, target, config)?,
            Strand::Reverse => {
                let mut reverse = std::mem::take(&mut self.reverse);
                reverse.clear();
                reverse.reserve(target.len());
                reverse.extend(target.iter().rev().map(|base| complement(*base)));
                if self.observed {
                    self.work.reverse_bytes += target.len() as u64;
                }
                let result = self.align_raw(query, &reverse, config);
                self.reverse = reverse;
                result?
            }
        };
        let cpu = observed_cpu(self.observed);
        let started = self.observed.then(std::time::Instant::now);
        let result = finish(raw, strand, target_offset, target.len());
        if let Some(started) = started {
            self.traceback_nanoseconds += started.elapsed().as_nanos() as u64;
            self.work.local_trace_cpu_ns += elapsed_cpu(cpu);
        }
        result
    }

    pub fn complete_endpoints(
        &mut self,
        core: Alignment,
        query: &[u8],
        target: &[u8],
        target_offset: u64,
        max_extension: usize,
        config: AlignmentConfig,
    ) -> Result<EndpointCompletion, AlignmentError> {
        let cpu = observed_cpu(self.observed);
        let started = self.observed.then(std::time::Instant::now);
        config.validate()?;
        let result = match core.strand {
            Strand::Forward => complete_endpoints(
                &mut self.endpoint_cells,
                core,
                query,
                target,
                target_offset,
                max_extension,
                config,
                self.observed.then_some(&mut self.work),
            ),
            Strand::Reverse => {
                let mut reverse = std::mem::take(&mut self.reverse);
                reverse.clear();
                reverse.reserve(target.len());
                reverse.extend(target.iter().rev().map(|base| complement(*base)));
                if self.observed {
                    self.work.reverse_bytes += target.len() as u64;
                }
                let result = complete_endpoints(
                    &mut self.endpoint_cells,
                    core,
                    query,
                    &reverse,
                    target_offset,
                    max_extension,
                    config,
                    self.observed.then_some(&mut self.work),
                );
                self.reverse = reverse;
                result
            }
        };
        if let Some(started) = started {
            self.work.endpoint_total_cpu_ns += elapsed_cpu(cpu);
            self.endpoint_nanoseconds += started.elapsed().as_nanos() as u64;
        }
        result
    }

    fn align_raw(
        &mut self,
        query: &[u8],
        target: &[u8],
        config: AlignmentConfig,
    ) -> Result<RawAlignment, AlignmentError> {
        #[cfg(target_arch = "x86_64")]
        if is_x86_feature_detected!("avx2") {
            // SAFETY: the runtime feature check guards every AVX2 instruction in this path.
            return unsafe { self.align_raw_avx2(query, target, config) };
        }
        let result = self.align_raw_scalar(query, target, config);
        #[cfg(target_arch = "x86_64")]
        if self.resident_chunks && matches!(result, Err(AlignmentError::MatrixTooLarge { .. })) {
            return self.align_raw_waves_scalar(query, target, config);
        }
        result
    }

    fn align_raw_scalar(
        &mut self,
        query: &[u8],
        target: &[u8],
        config: AlignmentConfig,
    ) -> Result<RawAlignment, AlignmentError> {
        let init_cpu = observed_cpu(self.observed);
        config.validate()?;
        if query.is_empty() {
            return Err(AlignmentError::EmptyQuery);
        }
        if target.is_empty() {
            return Err(AlignmentError::EmptyTarget);
        }

        self.prepare_rows(query.len(), target.len(), config)?;
        let total_cells = self
            .row_offsets
            .last()
            .copied()
            .unwrap_or(0)
            .checked_add(self.row_widths.last().copied().unwrap_or(0))
            .ok_or(AlignmentError::LengthOverflow)?;
        if total_cells == 0 {
            return Err(AlignmentError::BandExcludesInput);
        }
        if total_cells > config.max_cells {
            return Err(AlignmentError::MatrixTooLarge {
                cells: total_cells,
                max_cells: config.max_cells,
            });
        }
        if self.observed {
            self.work.growth_operations += u64::from(self.cells.capacity() < total_cells);
        }
        if self.cells.len() < total_cells {
            self.cells.resize(total_cells, Cell::default());
        } else {
            self.cells[..total_cells].fill(Cell::default());
            self.cells.truncate(total_cells);
        }

        if self.observed {
            self.work.local_passes += 1;
            self.work.local_cells += total_cells as u64;
            self.work.local_scalar_cells += total_cells as u64;
            self.work.local_boundary_cells += self.row_widths[0] as u64
                + self
                    .row_starts
                    .iter()
                    .zip(&self.row_widths)
                    .skip(1)
                    .filter(|&(start, width)| *start == 0 && *width > 0)
                    .count() as u64;
            self.work.local_init_cpu_ns += elapsed_cpu(init_cpu);
        }
        let matrix_cpu = observed_cpu(self.observed);
        let gap_open_score = gap_open(config);
        let mut best = BestCell::default();
        for query_index in 0..=query.len() {
            let row_offset = self.row_offsets[query_index];
            let row_start = self.row_starts[query_index];
            let row_width = self.row_widths[query_index];
            let previous_metadata = query_index.checked_sub(1).map(|previous_index| {
                (
                    self.row_offsets[previous_index],
                    self.row_starts[previous_index],
                    self.row_widths[previous_index],
                )
            });
            let (completed_rows, current_and_following) = self.cells.split_at_mut(row_offset);
            let current_row = &mut current_and_following[..row_width];
            let previous_row = previous_metadata
                .map(|(offset, start, width)| (start, &completed_rows[offset..offset + width]));
            for local_index in 0..row_width {
                let target_index = row_start + local_index;
                if query_index == 0 && target_index == 0 {
                    continue;
                }
                let mut cell = Cell::default();
                if let Some((previous_start, previous_cells)) = previous_row
                    && target_index > 0
                    && let Some(previous_index) = (target_index - 1).checked_sub(previous_start)
                    && let Some(previous) = previous_cells.get(previous_index).copied()
                {
                    let (score, state) = previous.best_score();
                    let substitution =
                        if query[query_index - 1].eq_ignore_ascii_case(&target[target_index - 1]) {
                            config.match_score
                        } else {
                            config.mismatch_score
                        };
                    let score = score.saturating_add(substitution);
                    if score > 0 {
                        cell.scores[MATCH as usize] = score;
                        cell.previous[MATCH as usize] =
                            if score == substitution { START } else { state };
                    }
                }
                if local_index > 0 {
                    let previous = current_row[local_index - 1];
                    let (score, state) = choose([
                        (
                            previous.scores[INSERTION as usize] + config.gap_extend_score,
                            INSERTION,
                        ),
                        (previous.scores[MATCH as usize] + gap_open_score, MATCH),
                        (
                            previous.scores[DELETION as usize] + gap_open_score,
                            DELETION,
                        ),
                    ]);
                    if score > 0 {
                        cell.scores[INSERTION as usize] = score;
                        cell.previous[INSERTION as usize] = state;
                    }
                }
                if let Some((previous_start, previous_cells)) = previous_row
                    && let Some(previous_index) = target_index.checked_sub(previous_start)
                    && let Some(previous) = previous_cells.get(previous_index).copied()
                {
                    let (score, state) = choose([
                        (
                            previous.scores[DELETION as usize] + config.gap_extend_score,
                            DELETION,
                        ),
                        (previous.scores[MATCH as usize] + gap_open_score, MATCH),
                        (
                            previous.scores[INSERTION as usize] + gap_open_score,
                            INSERTION,
                        ),
                    ]);
                    if score > 0 {
                        cell.scores[DELETION as usize] = score;
                        cell.previous[DELETION as usize] = state;
                    }
                }
                current_row[local_index] = cell;
                best.consider(query_index, target_index, cell);
            }
        }
        self.work.local_matrix_cpu_ns += elapsed_cpu(matrix_cpu);
        if best.score <= 0 {
            return Err(AlignmentError::NoAlignment);
        }

        let (query_start, target_start) = self.traceback_scalar(query, target, best)?;
        let edit_script = runs_from_operations(&self.operations)?;
        let summary = summarize_runs(&edit_script)?;
        Ok(RawAlignment {
            score: best.score,
            query_interval: Interval::new(query_start as u64, best.query_index as u64)?,
            target_interval: Interval::new(target_start as u64, best.target_index as u64)?,
            cigar: cigar_from_runs(&edit_script)?,
            edit_script,
            summary,
        })
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    unsafe fn align_raw_avx2(
        &mut self,
        query: &[u8],
        target: &[u8],
        config: AlignmentConfig,
    ) -> Result<RawAlignment, AlignmentError> {
        // SAFETY: this function enables AVX2 for the inlined wave driver.
        unsafe { self.align_raw_waves::<true>(query, target, config) }
    }

    #[cfg(target_arch = "x86_64")]
    fn align_raw_waves_scalar(
        &mut self,
        query: &[u8],
        target: &[u8],
        config: AlignmentConfig,
    ) -> Result<RawAlignment, AlignmentError> {
        // SAFETY: the scalar instantiation never calls AVX2 code.
        unsafe { self.align_raw_waves::<false>(query, target, config) }
    }

    #[cfg(target_arch = "x86_64")]
    #[inline(always)]
    unsafe fn align_raw_waves<const VECTOR: bool>(
        &mut self,
        query: &[u8],
        target: &[u8],
        config: AlignmentConfig,
    ) -> Result<RawAlignment, AlignmentError> {
        let init_cpu = observed_cpu(self.observed);
        config.validate()?;
        if query.is_empty() {
            return Err(AlignmentError::EmptyQuery);
        }
        if target.is_empty() {
            return Err(AlignmentError::EmptyTarget);
        }

        self.prepare_rows(query.len(), target.len(), config)?;
        let total_cells = self
            .row_offsets
            .last()
            .copied()
            .unwrap_or(0)
            .checked_add(self.row_widths.last().copied().unwrap_or(0))
            .ok_or(AlignmentError::LengthOverflow)?;
        if total_cells == 0 {
            return Err(AlignmentError::BandExcludesInput);
        }
        let last_wave = query
            .len()
            .checked_add(target.len())
            .ok_or(AlignmentError::LengthOverflow)?;
        let mut max_wave_width = 0usize;
        let mut wave_cells = 0usize;
        for wave in 0..=last_wave {
            if let Some(range) = wave_range(query.len(), target.len(), config, wave) {
                max_wave_width = max_wave_width.max(range.width());
                wave_cells = wave_cells.saturating_add(range.width());
            }
        }
        debug_assert_eq!(wave_cells, total_cells);
        // Traced tasks may exceed max_cells: they keep at most max_cells traceback cells
        // resident and recompute earlier waves from checkpoints. Direct callers keep the limit.
        let chunked = total_cells > config.max_cells;
        if chunked && (!self.resident_chunks || max_wave_width > config.max_cells) {
            return Err(AlignmentError::MatrixTooLarge {
                cells: total_cells,
                max_cells: config.max_cells,
            });
        }
        if !chunked {
            self.layout_chunk(query.len(), target.len(), config, 0, last_wave)?;
        }
        let wave_scores = max_wave_width
            .checked_mul(3)
            .ok_or(AlignmentError::LengthOverflow)?;
        for wave in &mut self.waves {
            if self.observed {
                self.work.growth_operations += u64::from(wave.capacity() < wave_scores);
            }
            prepare_wave(wave, wave_scores);
        }

        if self.observed {
            self.work.local_narrow_eligible_passes +=
                u64::from(narrow_local_scores(query.len(), target.len(), config));
            self.work.local_passes += 1;
            self.work.local_cells += total_cells as u64;
            self.work.local_chunked_passes += u64::from(chunked);
            self.work.local_init_cpu_ns += elapsed_cpu(init_cpu);
        }
        let matrix_cpu = observed_cpu(self.observed);
        let mut best = BestCell::default();
        let mut carry = WaveCarry::default();
        if chunked {
            self.chunk_waves.clear();
            self.chunk_carries.clear();
            self.checkpoints.clear();
            let width = |wave| {
                wave_range(query.len(), target.len(), config, wave).map_or(0, WaveRange::width)
            };
            let mut first = 0;
            while first <= last_wave {
                let mut last = first;
                let mut cells = width(first);
                while last < last_wave && cells + width(last + 1) <= config.max_cells {
                    last += 1;
                    cells += width(last);
                }
                self.save_checkpoint(carry, max_wave_width);
                self.chunk_waves.push((first, last));
                self.layout_chunk(query.len(), target.len(), config, first, last)?;
                // SAFETY: forwarded from this function's contract.
                unsafe {
                    self.fill_waves_selected::<VECTOR>(
                        query,
                        target,
                        config,
                        (first, last, max_wave_width),
                        &mut carry,
                        &mut best,
                        self.observed,
                    );
                }
                self.work.local_chunks += u64::from(self.observed);
                first = last + 1;
            }
        } else {
            // SAFETY: forwarded from this function's contract.
            unsafe {
                self.fill_waves_selected::<VECTOR>(
                    query,
                    target,
                    config,
                    (0, last_wave, max_wave_width),
                    &mut carry,
                    &mut best,
                    self.observed,
                );
            }
        }
        self.work.local_matrix_cpu_ns += elapsed_cpu(matrix_cpu);
        if best.score <= 0 {
            return Err(AlignmentError::NoAlignment);
        }

        let (query_start, target_start) = if chunked {
            // SAFETY: forwarded from this function's contract.
            unsafe {
                self.traceback_chunked::<VECTOR>(query, target, config, max_wave_width, best)?
            }
        } else {
            self.traceback_compact(query, target, config, best)?
        };
        let edit_script = runs_from_operations(&self.operations)?;
        let summary = summarize_runs(&edit_script)?;
        Ok(RawAlignment {
            score: best.score,
            query_interval: Interval::new(query_start as u64, best.query_index as u64)?,
            target_interval: Interval::new(target_start as u64, best.target_index as u64)?,
            cigar: cigar_from_runs(&edit_script)?,
            edit_script,
            summary,
        })
    }

    /// Runs `fill_waves` with its AVX2 recurrence in one out-of-line body, so the three wave
    /// drivers share one inlined copy of the vector kernel. Callers must enable AVX2.
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    #[inline(never)]
    #[allow(clippy::too_many_arguments)]
    unsafe fn fill_waves_avx2(
        &mut self,
        query: &[u8],
        target: &[u8],
        config: AlignmentConfig,
        waves: (usize, usize, usize),
        carry: &mut WaveCarry,
        best: &mut BestCell,
        count: bool,
    ) {
        // SAFETY: this function enables AVX2 for the vector instantiation.
        unsafe { self.fill_waves::<i32, true>(query, target, config, waves, carry, best, count) }
    }

    /// Fills waves with the vector or scalar recurrence. Callers must enable AVX2 when VECTOR is
    /// true.
    #[cfg(target_arch = "x86_64")]
    #[inline(always)]
    #[allow(clippy::too_many_arguments)]
    unsafe fn fill_waves_selected<const VECTOR: bool>(
        &mut self,
        query: &[u8],
        target: &[u8],
        config: AlignmentConfig,
        waves: (usize, usize, usize),
        carry: &mut WaveCarry,
        best: &mut BestCell,
        count: bool,
    ) {
        // SAFETY: forwarded from this function's contract.
        unsafe {
            if VECTOR {
                self.fill_waves_avx2(query, target, config, waves, carry, best, count)
            } else {
                self.fill_waves::<i32, false>(query, target, config, waves, carry, best, count)
            }
        }
    }

    /// Fills waves first..=last with the retained recurrence. Callers must enable AVX2 when
    /// VECTOR is true.
    #[cfg(target_arch = "x86_64")]
    #[inline(always)]
    #[allow(clippy::too_many_arguments)]
    unsafe fn fill_waves<S: WaveScore, const VECTOR: bool>(
        &mut self,
        query: &[u8],
        target: &[u8],
        config: AlignmentConfig,
        (first_wave, last_wave, stride): (usize, usize, usize),
        carry: &mut WaveCarry,
        best: &mut BestCell,
        count: bool,
    ) {
        let gap_open_score = gap_open(config);
        for wave in first_wave..=last_wave {
            let current_range = wave_range(query.len(), target.len(), config, wave);
            if let Some(current_range) = current_range {
                let wave_base = self.wave_offsets[wave - first_wave];
                let (older_previous, current) = self.waves.split_at_mut(2);
                let older_scores: &[S] = bytemuck::cast_slice(&older_previous[0][..]);
                let previous_scores: &[S] = bytemuck::cast_slice(&older_previous[1][..]);
                let older = carry.older.map(|range| ScoreWave {
                    range,
                    scores: older_scores,
                    stride,
                });
                let previous = carry.previous.map(|range| ScoreWave {
                    range,
                    scores: previous_scores,
                    stride,
                });
                let current: &mut [S] = bytemuck::cast_slice_mut(&mut current[0][..]);
                let mut row = current_range.start;
                let vector_range =
                    older
                        .zip(previous)
                        .filter(|_| VECTOR)
                        .and_then(|(older, previous)| {
                            let start = current_range
                                .start
                                .max(older.range.start.saturating_add(1))
                                .max(previous.range.start.saturating_add(1))
                                .max(1);
                            let end = current_range
                                .end
                                .min(older.range.end.saturating_add(1))
                                .min(previous.range.end)
                                .min(wave.saturating_sub(1));
                            (start <= end).then_some((start, end, older, previous))
                        });
                if count {
                    let vectors = vector_range.map_or(0, |(start, end, _, _)| {
                        Some(end - start + 1)
                            .filter(|&width| width >= S::LANES)
                            .unwrap_or(0)
                    });
                    self.work.local_vector8_cells += vectors as u64;
                    self.work.local_scalar_cells += (current_range.width() - vectors) as u64;
                    self.work.local_boundary_cells += u64::from(current_range.contains(0))
                        + u64::from(current_range.contains(wave))
                        - u64::from(wave == 0 && current_range.contains(0));
                }
                if let Some((start, end, older, previous)) = vector_range {
                    while row < start {
                        fill_wave_scalar(
                            query,
                            target,
                            config,
                            gap_open_score,
                            wave,
                            row,
                            current_range,
                            stride,
                            Some(older),
                            Some(previous),
                            current,
                            &mut self.compact_cells,
                            wave_base,
                            best,
                        );
                        row += 1;
                    }
                    while row
                        .checked_add(S::LANES - 1)
                        .is_some_and(|last| last <= end)
                    {
                        // SAFETY: vector_range proves LANES current, diagonal, left, and above
                        // cells are in their wave slices. row>=1 and row+LANES-1<=wave-1 prove
                        // the query and reverse-loaded target bytes are also in bounds.
                        unsafe {
                            S::fill_block(
                                query,
                                target,
                                config,
                                gap_open_score,
                                wave,
                                row,
                                current_range,
                                stride,
                                older,
                                previous,
                                current,
                                &mut self.compact_cells,
                                wave_base,
                                best,
                            );
                        }
                        row += S::LANES;
                    }
                    // Cells of one wave depend only on the two earlier waves, so a final block
                    // overlapping finished cells rewrites identical scores and traces, and
                    // BestCell ignores a repeated identical cell.
                    if row <= end && end + 1 - start >= S::LANES {
                        // SAFETY: the last LANES rows through end lie inside vector_range, as for
                        // the blocks above.
                        unsafe {
                            S::fill_block(
                                query,
                                target,
                                config,
                                gap_open_score,
                                wave,
                                end + 1 - S::LANES,
                                current_range,
                                stride,
                                older,
                                previous,
                                current,
                                &mut self.compact_cells,
                                wave_base,
                                best,
                            );
                        }
                        row = end + 1;
                    }
                }
                while row <= current_range.end {
                    fill_wave_scalar(
                        query,
                        target,
                        config,
                        gap_open_score,
                        wave,
                        row,
                        current_range,
                        stride,
                        older,
                        previous,
                        current,
                        &mut self.compact_cells,
                        wave_base,
                        best,
                    );
                    row += 1;
                }
            }
            self.waves.rotate_left(1);
            carry.older = carry.previous;
            carry.previous = current_range;
        }
    }

    /// Wave layout of the traceback cells in waves first..=last: each wave's cells are
    /// contiguous, so one vector block writes its eight traces together.
    #[cfg(target_arch = "x86_64")]
    fn layout_chunk(
        &mut self,
        query_len: usize,
        target_len: usize,
        config: AlignmentConfig,
        first_wave: usize,
        last_wave: usize,
    ) -> Result<usize, AlignmentError> {
        let mut total = 0usize;
        self.wave_offsets.clear();
        for wave in first_wave..=last_wave {
            self.wave_offsets.push(total);
            let width = wave_range(query_len, target_len, config, wave).map_or(0, WaveRange::width);
            total = total
                .checked_add(width)
                .ok_or(AlignmentError::LengthOverflow)?;
        }
        self.wave_offsets.push(total);
        if self.observed {
            self.work.growth_operations += u64::from(self.compact_cells.capacity() < total);
        }
        self.compact_cells.resize(total, 0);
        Ok(total)
    }

    #[cfg(target_arch = "x86_64")]
    fn save_checkpoint(&mut self, carry: WaveCarry, stride: usize) {
        self.chunk_carries.push(carry);
        self.checkpoints
            .extend_from_slice(&self.waves[0][..3 * stride]);
        self.checkpoints
            .extend_from_slice(&self.waves[1][..3 * stride]);
    }

    #[cfg(target_arch = "x86_64")]
    #[inline(always)]
    unsafe fn traceback_chunked<const VECTOR: bool>(
        &mut self,
        query: &[u8],
        target: &[u8],
        config: AlignmentConfig,
        stride: usize,
        best: BestCell,
    ) -> Result<(usize, usize), AlignmentError> {
        let cpu = observed_cpu(self.observed);
        let started = self.observed.then(std::time::Instant::now);
        self.operations.clear();
        let mut position = (best.query_index, best.target_index, best.state);
        let mut chunk = self
            .chunk_waves
            .partition_point(|&(first, _)| first <= position.0 + position.1)
            .checked_sub(1)
            .ok_or(AlignmentError::InvalidTraceback)?;
        loop {
            let (first, last) = self.chunk_waves[chunk];
            let stop = last.min(position.0 + position.1);
            let mut carry = self.chunk_carries[chunk];
            let saved = &self.checkpoints[chunk * 6 * stride..(chunk + 1) * 6 * stride];
            self.waves[0][..3 * stride].copy_from_slice(&saved[..3 * stride]);
            self.waves[1][..3 * stride].copy_from_slice(&saved[3 * stride..]);
            let cells = self.layout_chunk(query.len(), target.len(), config, first, stop)?;
            if self.observed {
                self.work.local_recomputed_cells += cells as u64;
            }
            let mut ignored = BestCell::default();
            // SAFETY: forwarded from this function's contract.
            unsafe {
                self.fill_waves_selected::<VECTOR>(
                    query,
                    target,
                    config,
                    (first, stop, stride),
                    &mut carry,
                    &mut ignored,
                    false,
                );
            }
            if self.trace_operations(query, target, config, &mut position, first)? {
                break;
            }
            chunk = chunk
                .checked_sub(1)
                .ok_or(AlignmentError::InvalidTraceback)?;
        }
        self.operations.reverse();
        if let Some(started) = started {
            self.traceback_nanoseconds += started.elapsed().as_nanos() as u64;
            self.work.local_trace_cpu_ns += elapsed_cpu(cpu);
        }
        Ok((position.0, position.1))
    }

    fn prepare_rows(
        &mut self,
        query_len: usize,
        target_len: usize,
        config: AlignmentConfig,
    ) -> Result<(), AlignmentError> {
        self.row_offsets.clear();
        self.row_starts.clear();
        self.row_widths.clear();
        let mut total = 0usize;
        for query_index in 0..=query_len {
            let (start, width) = band_row(
                query_index,
                target_len,
                config.diagonal_offset,
                config.band_width,
            )?
            .map_or((0, 0), |(start, end)| (start, end - start + 1));
            self.row_offsets.push(total);
            self.row_starts.push(start);
            self.row_widths.push(width);
            total = total
                .checked_add(width)
                .ok_or(AlignmentError::LengthOverflow)?;
        }
        Ok(())
    }

    fn cell_index_checked(&self, query_index: usize, target_index: usize) -> Option<usize> {
        let start = *self.row_starts.get(query_index)?;
        let width = *self.row_widths.get(query_index)?;
        (target_index >= start && target_index < start + width)
            .then(|| self.row_offsets[query_index] + target_index - start)
    }

    fn traceback_scalar(
        &mut self,
        query: &[u8],
        target: &[u8],
        best: BestCell,
    ) -> Result<(usize, usize), AlignmentError> {
        let cpu = observed_cpu(self.observed);
        let started = self.observed.then(std::time::Instant::now);
        let mut query_index = best.query_index;
        let mut target_index = best.target_index;
        let mut state = best.state;
        self.operations.clear();
        while query_index > 0 || target_index > 0 {
            let index = self
                .cell_index_checked(query_index, target_index)
                .ok_or(AlignmentError::TracebackOutsideBand)?;
            let cell = self.cells[index];
            if state > DELETION || cell.scores[state as usize] <= 0 {
                break;
            }
            let previous = cell.previous[state as usize];
            match state {
                MATCH if query_index > 0 && target_index > 0 => {
                    self.operations.push(
                        if query[query_index - 1].eq_ignore_ascii_case(&target[target_index - 1]) {
                            EditOperation::Equal
                        } else {
                            EditOperation::Substitution
                        },
                    );
                    query_index -= 1;
                    target_index -= 1;
                }
                INSERTION if target_index > 0 => {
                    self.operations.push(EditOperation::Insertion);
                    target_index -= 1;
                }
                DELETION if query_index > 0 => {
                    self.operations.push(EditOperation::Deletion);
                    query_index -= 1;
                }
                _ => return Err(AlignmentError::InvalidTraceback),
            }
            if previous == START {
                break;
            }
            state = previous;
        }
        self.operations.reverse();
        if let Some(started) = started {
            self.traceback_nanoseconds += started.elapsed().as_nanos() as u64;
            self.work.local_trace_cpu_ns += elapsed_cpu(cpu);
        }
        Ok((query_index, target_index))
    }

    #[cfg(target_arch = "x86_64")]
    fn traceback_compact(
        &mut self,
        query: &[u8],
        target: &[u8],
        config: AlignmentConfig,
        best: BestCell,
    ) -> Result<(usize, usize), AlignmentError> {
        let cpu = observed_cpu(self.observed);
        let started = self.observed.then(std::time::Instant::now);
        self.operations.clear();
        let mut position = (best.query_index, best.target_index, best.state);
        self.trace_operations(query, target, config, &mut position, 0)?;
        self.operations.reverse();
        if let Some(started) = started {
            self.traceback_nanoseconds += started.elapsed().as_nanos() as u64;
            self.work.local_trace_cpu_ns += elapsed_cpu(cpu);
        }
        Ok((position.0, position.1))
    }

    /// Appends reversed operations until the local start, or returns false before a wave
    /// below first_wave whose traceback cells are not resident.
    #[cfg(target_arch = "x86_64")]
    fn trace_operations(
        &mut self,
        query: &[u8],
        target: &[u8],
        config: AlignmentConfig,
        position: &mut (usize, usize, u8),
        first_wave: usize,
    ) -> Result<bool, AlignmentError> {
        let (mut query_index, mut target_index, mut state) = *position;
        let complete = loop {
            if query_index == 0 && target_index == 0 {
                break true;
            }
            if query_index + target_index < first_wave {
                break false;
            }
            let wave = query_index + target_index;
            let index = wave_range(query.len(), target.len(), config, wave)
                .filter(|range| range.contains(query_index))
                .and_then(|range| {
                    let base = *self.wave_offsets.get(wave - first_wave)?;
                    (base + query_index - range.start < *self.wave_offsets.last()?)
                        .then_some(base + query_index - range.start)
                })
                .ok_or(AlignmentError::TracebackOutsideBand)?;
            let traceback = self.compact_cells[index];
            if state > DELETION || !traceback_positive(traceback, state) {
                break true;
            }
            let previous = traceback_previous(traceback, state);
            match state {
                MATCH if query_index > 0 && target_index > 0 => {
                    self.operations.push(
                        if query[query_index - 1].eq_ignore_ascii_case(&target[target_index - 1]) {
                            EditOperation::Equal
                        } else {
                            EditOperation::Substitution
                        },
                    );
                    query_index -= 1;
                    target_index -= 1;
                }
                INSERTION if target_index > 0 => {
                    self.operations.push(EditOperation::Insertion);
                    target_index -= 1;
                }
                DELETION if query_index > 0 => {
                    self.operations.push(EditOperation::Deletion);
                    query_index -= 1;
                }
                _ => return Err(AlignmentError::InvalidTraceback),
            }
            if previous == START {
                break true;
            }
            state = previous;
        };
        *position = (query_index, target_index, state);
        Ok(complete)
    }
}

/// Inclusive target columns stored for one query row of a clipped diagonal band.
pub(crate) fn band_row(
    query_index: usize,
    target_len: usize,
    diagonal_offset: i64,
    band_width: u32,
) -> Result<Option<(usize, usize)>, AlignmentError> {
    let target_len = i128::try_from(target_len).map_err(|_| AlignmentError::LengthOverflow)?;
    let center = i128::try_from(query_index).map_err(|_| AlignmentError::LengthOverflow)?
        + i128::from(diagonal_offset);
    let low = center - i128::from(band_width);
    let high = center + i128::from(band_width);
    if high < 0 || low > target_len {
        return Ok(None);
    }
    let start = usize::try_from(low.max(0)).map_err(|_| AlignmentError::LengthOverflow)?;
    let end = usize::try_from(high.min(target_len)).map_err(|_| AlignmentError::LengthOverflow)?;
    Ok(Some((start, end)))
}

/// Exact local cells stored for query rows 0..=query_len, counted without row metadata.
pub(crate) fn band_cells(
    query_len: usize,
    target_len: usize,
    diagonal_offset: i64,
    band_width: u32,
) -> Result<usize, AlignmentError> {
    (0..=query_len).try_fold(0usize, |total, query_index| {
        let width = band_row(query_index, target_len, diagonal_offset, band_width)?
            .map_or(0, |(start, end)| end - start + 1);
        total
            .checked_add(width)
            .ok_or(AlignmentError::LengthOverflow)
    })
}

#[cfg(target_arch = "x86_64")]
fn prepare_wave(scores: &mut Vec<i32>, length: usize) {
    if scores.capacity() < length {
        scores.reserve_exact(length - scores.len());
    }
    scores.resize(length, 0);
}

/// Score wave ranges carried into the next wave, saved at each resident chunk start.
#[cfg(target_arch = "x86_64")]
#[derive(Clone, Copy, Debug, Default)]
struct WaveCarry {
    older: Option<WaveRange>,
    previous: Option<WaveRange>,
}

#[cfg(target_arch = "x86_64")]
#[derive(Clone, Copy, Debug)]
struct WaveRange {
    start: usize,
    end: usize,
}

#[cfg(target_arch = "x86_64")]
impl WaveRange {
    fn width(self) -> usize {
        self.end - self.start + 1
    }

    fn contains(self, row: usize) -> bool {
        row >= self.start && row <= self.end
    }
}

#[cfg(target_arch = "x86_64")]
fn wave_range(
    query_len: usize,
    target_len: usize,
    config: AlignmentConfig,
    wave: usize,
) -> Option<WaveRange> {
    let query_len = i128::try_from(query_len).ok()?;
    let target_len = i128::try_from(target_len).ok()?;
    let wave = i128::try_from(wave).ok()?;
    let diagonal = i128::from(config.diagonal_offset);
    let band = i128::from(config.band_width);
    let lower = 0
        .max(wave - target_len)
        .max(ceil_div2(wave - (diagonal + band)));
    let upper = query_len
        .min(wave)
        .min((wave - (diagonal - band)).div_euclid(2));
    if lower > upper {
        return None;
    }
    Some(WaveRange {
        start: usize::try_from(lower).ok()?,
        end: usize::try_from(upper).ok()?,
    })
}

#[cfg(target_arch = "x86_64")]
fn ceil_div2(value: i128) -> i128 {
    -(-value).div_euclid(2)
}

#[cfg(target_arch = "x86_64")]
#[derive(Clone, Copy)]
struct ScoreWave<'a, S = i32> {
    range: WaveRange,
    scores: &'a [S],
    stride: usize,
}

#[cfg(target_arch = "x86_64")]
impl<S: WaveScore> ScoreWave<'_, S> {
    fn cell(self, row: usize) -> Option<[i32; 3]> {
        self.range.contains(row).then(|| {
            let offset = row - self.range.start;
            [
                self.scores[offset].get(),
                self.scores[self.stride + offset].get(),
                self.scores[2 * self.stride + offset].get(),
            ]
        })
    }
}

/// Score plane element of the local wave recurrence.
#[cfg(target_arch = "x86_64")]
trait WaveScore: bytemuck::Pod {
    const LANES: usize;

    fn get(self) -> i32;

    fn put(score: i32) -> Self;

    /// Evaluates LANES cells of one wave. Callers must enable AVX2 and prove the block bounds
    /// that `fill_wave_avx2` documents.
    #[allow(clippy::too_many_arguments)]
    unsafe fn fill_block(
        query: &[u8],
        target: &[u8],
        config: AlignmentConfig,
        gap_open_score: i32,
        wave: usize,
        row: usize,
        current_range: WaveRange,
        stride: usize,
        older: ScoreWave<'_, Self>,
        previous: ScoreWave<'_, Self>,
        current: &mut [Self],
        traceback: &mut [u16],
        wave_base: usize,
        best: &mut BestCell,
    );
}

#[cfg(target_arch = "x86_64")]
impl WaveScore for i32 {
    const LANES: usize = 8;

    fn get(self) -> i32 {
        self
    }

    fn put(score: i32) -> Self {
        score
    }

    #[inline(always)]
    unsafe fn fill_block(
        query: &[u8],
        target: &[u8],
        config: AlignmentConfig,
        gap_open_score: i32,
        wave: usize,
        row: usize,
        current_range: WaveRange,
        stride: usize,
        older: ScoreWave<'_, Self>,
        previous: ScoreWave<'_, Self>,
        current: &mut [Self],
        traceback: &mut [u16],
        wave_base: usize,
        best: &mut BestCell,
    ) {
        // SAFETY: forwarded from this function's contract.
        unsafe {
            fill_wave_avx2(
                query,
                target,
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
            )
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[allow(clippy::too_many_arguments)]
fn fill_wave_scalar<S: WaveScore>(
    query: &[u8],
    target: &[u8],
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
) {
    let target_index = wave - row;
    let wave_offset = row - current_range.start;
    let trace_index = wave_base + row - current_range.start;
    if row == 0 && target_index == 0 {
        current[wave_offset] = S::put(0);
        current[stride + wave_offset] = S::put(0);
        current[2 * stride + wave_offset] = S::put(0);
        traceback[trace_index] = 0;
        return;
    }
    let mut cell = Cell::default();
    if row > 0
        && target_index > 0
        && let Some(previous) = older.and_then(|wave| wave.cell(row - 1))
    {
        let (score, state) = best_score(previous);
        let substitution = if query[row - 1].eq_ignore_ascii_case(&target[target_index - 1]) {
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
    if let Some(previous) = previous.and_then(|wave| wave.cell(row)) {
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
    if row > 0
        && let Some(previous) = previous.and_then(|wave| wave.cell(row - 1))
    {
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
    current[wave_offset] = S::put(cell.scores[MATCH as usize]);
    current[stride + wave_offset] = S::put(cell.scores[INSERTION as usize]);
    current[2 * stride + wave_offset] = S::put(cell.scores[DELETION as usize]);
    traceback[trace_index] = encode_traceback(cell);
    best.consider(row, target_index, cell);
}

#[cfg(target_arch = "x86_64")]
fn narrow_local_scores(query_len: usize, target_len: usize, config: AlignmentConfig) -> bool {
    // Stored local states are in [0, min(lengths)*match]: gaps and mismatches cannot
    // increase them. Every raw candidate is a stored state plus one scoring term.
    let upper = i64::try_from(query_len.min(target_len))
        .ok()
        .and_then(|len| len.checked_mul(i64::from(config.match_score)))
        .and_then(|score| score.checked_add(i64::from(config.match_score)));
    config.validate().is_ok()
        && upper.is_some_and(|score| score <= i64::from(i16::MAX))
        && [
            config.match_score,
            config.mismatch_score,
            gap_open(config),
            config.gap_extend_score,
        ]
        .iter()
        .all(|&score| i16::try_from(score).is_ok())
}

#[cfg(target_arch = "x86_64")]
fn best_score(scores: [i32; 3]) -> (i32, u8) {
    choose([
        (scores[MATCH as usize], MATCH),
        (scores[INSERTION as usize], INSERTION),
        (scores[DELETION as usize], DELETION),
    ])
}

#[cfg(target_arch = "x86_64")]
fn encode_traceback(cell: Cell) -> u16 {
    let mut traceback = 0u16;
    for state in 0..=DELETION {
        let index = state as usize;
        if cell.scores[index] > 0 {
            debug_assert!(cell.previous[index] <= START);
            traceback |= 1 << state;
            traceback |= u16::from(cell.previous[index]) << (3 + 2 * state);
        }
    }
    traceback
}

#[cfg(target_arch = "x86_64")]
fn traceback_positive(traceback: u16, state: u8) -> bool {
    traceback & (1 << state) != 0
}

#[cfg(target_arch = "x86_64")]
fn traceback_previous(traceback: u16, state: u8) -> u8 {
    ((traceback >> (3 + 2 * state)) & 3) as u8
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn load_wave(scores: ScoreWave<'_>, state: usize, row: usize) -> __m256i {
    debug_assert!(scores.range.contains(row));
    debug_assert!(scores.range.contains(row + 7));
    let offset = state * scores.stride + row - scores.range.start;
    // SAFETY: the caller proves rows row..row+7 are inside this wave and each state plane
    // has stride elements in a buffer of exactly three strides.
    unsafe { _mm256_loadu_si256(scores.scores.as_ptr().add(offset).cast()) }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn choose_avx2(m: __m256i, i: __m256i, d: __m256i) -> (__m256i, __m256i) {
    let zero = _mm256_setzero_si256();
    let one = _mm256_set1_epi32(1);
    let two = _mm256_set1_epi32(2);
    let i_better = _mm256_cmpgt_epi32(i, m);
    let mut score = _mm256_blendv_epi8(m, i, i_better);
    let mut state = _mm256_blendv_epi8(zero, one, i_better);
    let d_better = _mm256_cmpgt_epi32(d, score);
    score = _mm256_blendv_epi8(score, d, d_better);
    state = _mm256_blendv_epi8(state, two, d_better);
    (score, state)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn saturating_add_avx2(left: __m256i, right: __m256i) -> __m256i {
    let sum = _mm256_add_epi32(left, right);
    let overflow = _mm256_srai_epi32::<31>(_mm256_and_si256(
        _mm256_xor_si256(left, sum),
        _mm256_xor_si256(right, sum),
    ));
    let saturation = _mm256_xor_si256(_mm256_srai_epi32::<31>(left), _mm256_set1_epi32(i32::MAX));
    _mm256_blendv_epi8(sum, saturation, overflow)
}

#[cfg(target_arch = "x86_64")]
unsafe fn lowercase_ascii_8(bytes: __m128i) -> __m128i {
    unsafe {
        let upper = _mm_and_si128(
            _mm_cmpgt_epi8(bytes, _mm_set1_epi8((b'A' - 1) as i8)),
            _mm_cmpgt_epi8(_mm_set1_epi8((b'Z' + 1) as i8), bytes),
        );
        _mm_or_si128(bytes, _mm_and_si128(upper, _mm_set1_epi8(0x20)))
    }
}

/// AVX2 recurrence for eight cells. Inlined into its single AVX2-enabled caller, the vector
/// instantiation of `fill_waves`; callers must enable AVX2.
#[cfg(target_arch = "x86_64")]
#[allow(clippy::too_many_arguments)]
#[inline(always)]
unsafe fn fill_wave_avx2(
    query: &[u8],
    target: &[u8],
    config: AlignmentConfig,
    gap_open_score: i32,
    wave: usize,
    row: usize,
    current_range: WaveRange,
    stride: usize,
    older: ScoreWave<'_>,
    previous: ScoreWave<'_>,
    current: &mut [i32],
    traceback: &mut [u16],
    wave_base: usize,
    best: &mut BestCell,
) {
    unsafe {
        debug_assert!(row > 0);
        debug_assert!(row + 7 < wave);
        debug_assert!(row + 7 <= query.len());
        debug_assert!(current_range.contains(row));
        debug_assert!(current_range.contains(row + 7));
        debug_assert!(current.len() >= 3 * stride);
        let first_target = wave - row;
        debug_assert!(first_target >= 8 && first_target <= target.len());

        let query_bytes = _mm_loadl_epi64(query.as_ptr().add(row - 1).cast());
        let target_bytes = _mm_loadl_epi64(target.as_ptr().add(first_target - 8).cast());
        let reverse = _mm_setr_epi8(7, 6, 5, 4, 3, 2, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1);
        let target_bytes = _mm_shuffle_epi8(target_bytes, reverse);
        let equal = _mm_cmpeq_epi8(
            lowercase_ascii_8(query_bytes),
            lowercase_ascii_8(target_bytes),
        );
        let equal = _mm256_cvtepi8_epi32(equal);
        let substitution = _mm256_blendv_epi8(
            _mm256_set1_epi32(config.mismatch_score),
            _mm256_set1_epi32(config.match_score),
            equal,
        );

        let (diagonal, diagonal_state) = choose_avx2(
            load_wave(older, MATCH as usize, row - 1),
            load_wave(older, INSERTION as usize, row - 1),
            load_wave(older, DELETION as usize, row - 1),
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

        let left_m = _mm256_add_epi32(
            load_wave(previous, MATCH as usize, row),
            _mm256_set1_epi32(gap_open_score),
        );
        let left_i = _mm256_add_epi32(
            load_wave(previous, INSERTION as usize, row),
            _mm256_set1_epi32(config.gap_extend_score),
        );
        let left_d = _mm256_add_epi32(
            load_wave(previous, DELETION as usize, row),
            _mm256_set1_epi32(gap_open_score),
        );
        let (insertion_raw, insertion_previous) = choose_avx2(left_m, left_i, left_d);
        let insertion_positive = _mm256_cmpgt_epi32(insertion_raw, zero);
        let insertion_scores = _mm256_and_si256(insertion_raw, insertion_positive);

        let above_m = _mm256_add_epi32(
            load_wave(previous, MATCH as usize, row - 1),
            _mm256_set1_epi32(gap_open_score),
        );
        let above_i = _mm256_add_epi32(
            load_wave(previous, INSERTION as usize, row - 1),
            _mm256_set1_epi32(gap_open_score),
        );
        let above_d = _mm256_add_epi32(
            load_wave(previous, DELETION as usize, row - 1),
            _mm256_set1_epi32(config.gap_extend_score),
        );
        let (deletion_raw, deletion_previous) = choose_avx2(above_m, above_i, above_d);
        let deletion_positive = _mm256_cmpgt_epi32(deletion_raw, zero);
        let deletion_scores = _mm256_and_si256(deletion_raw, deletion_positive);

        let wave_offset = row - current_range.start;
        _mm256_storeu_si256(current.as_mut_ptr().add(wave_offset).cast(), match_scores);
        _mm256_storeu_si256(
            current.as_mut_ptr().add(stride + wave_offset).cast(),
            insertion_scores,
        );
        _mm256_storeu_si256(
            current.as_mut_ptr().add(2 * stride + wave_offset).cast(),
            deletion_scores,
        );

        let mut trace = _mm256_and_si256(match_positive, _mm256_set1_epi32(1));
        trace = _mm256_or_si256(
            trace,
            _mm256_slli_epi32::<3>(_mm256_and_si256(match_previous, match_positive)),
        );
        trace = _mm256_or_si256(
            trace,
            _mm256_and_si256(insertion_positive, _mm256_set1_epi32(2)),
        );
        trace = _mm256_or_si256(
            trace,
            _mm256_slli_epi32::<5>(_mm256_and_si256(insertion_previous, insertion_positive)),
        );
        trace = _mm256_or_si256(
            trace,
            _mm256_and_si256(deletion_positive, _mm256_set1_epi32(4)),
        );
        trace = _mm256_or_si256(
            trace,
            _mm256_slli_epi32::<7>(_mm256_and_si256(deletion_previous, deletion_positive)),
        );

        // Traces use nine bits, so unsigned saturation to 16 bits is exact. Packing places lanes
        // 0..3 and 4..7 in separate 128-bit halves; the permutation joins them in lane order.
        let packed = _mm256_permute4x64_epi64::<0xD8>(_mm256_packus_epi32(trace, zero));
        debug_assert!(wave_base + wave_offset + 8 <= traceback.len());
        // SAFETY: the wave layout stores this wave's cells contiguously, and rows row..=row+7
        // lie inside current_range, so the eight cells are inside the laid-out traceback.
        _mm_storeu_si128(
            traceback.as_mut_ptr().add(wave_base + wave_offset).cast(),
            _mm256_castsi256_si128(packed),
        );
        let (lane_best_scores, _) = choose_avx2(match_scores, insertion_scores, deletion_scores);
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
            let mut insertions = [0i32; 8];
            let mut deletions = [0i32; 8];
            _mm256_storeu_si256(matches.as_mut_ptr().cast(), match_scores);
            _mm256_storeu_si256(insertions.as_mut_ptr().cast(), insertion_scores);
            _mm256_storeu_si256(deletions.as_mut_ptr().cast(), deletion_scores);
            let best_lanes = _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpeq_epi32(
                lane_best_scores,
                block_best,
            )));
            let block_lane = best_lanes.trailing_zeros() as usize;
            // Lanes are in ascending query-index order. The first lane with the block maximum
            // is therefore lexicographically earliest; every lesser or later equal lane loses
            // after this unchanged BestCell comparison.
            let query_index = row + block_lane;
            let target_index = wave - query_index;
            best.consider(
                query_index,
                target_index,
                Cell {
                    scores: [
                        matches[block_lane],
                        insertions[block_lane],
                        deletions[block_lane],
                    ],
                    previous: [0; 3],
                },
            );
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct Cell {
    scores: [i32; 3],
    previous: [u8; 3],
}

impl Cell {
    fn best_score(self) -> (i32, u8) {
        choose([
            (self.scores[MATCH as usize], MATCH),
            (self.scores[INSERTION as usize], INSERTION),
            (self.scores[DELETION as usize], DELETION),
        ])
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct BestCell {
    score: i32,
    query_index: usize,
    target_index: usize,
    state: u8,
}

impl BestCell {
    fn consider(&mut self, query_index: usize, target_index: usize, cell: Cell) {
        let (score, state) = cell.best_score();
        if score > self.score
            || (score == self.score
                && score > 0
                && (query_index, target_index, state)
                    < (self.query_index, self.target_index, self.state))
        {
            *self = Self {
                score,
                query_index,
                target_index,
                state,
            };
        }
    }
}

struct RawAlignment {
    score: i32,
    query_interval: Interval,
    target_interval: Interval,
    cigar: String,
    edit_script: Vec<EditRun>,
    summary: RunSummary,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct EndpointMetrics {
    pub left_query_bases: usize,
    pub left_target_bases: usize,
    pub right_query_bases: usize,
    pub right_target_bases: usize,
    pub matrix_cells: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct EndpointCompletion {
    pub alignment: Alignment,
    pub metrics: EndpointMetrics,
}

const NEGATIVE: i32 = i32::MIN / 4;

#[derive(Debug, Default)]
struct EndpointMatrix {
    rows: [Vec<[i32; 3]>; 2],
    previous: Vec<u8>,
    /// Three anti-diagonal score waves of the AVX2 kernel: wave, state, then query row.
    waves: Vec<i32>,
    operations: Vec<EditOperation>,
}

#[derive(Debug, Default)]
struct EndpointWorkspace {
    matrix: EndpointMatrix,
    left_query: Vec<u8>,
    left_target: Vec<u8>,
}

impl EndpointWorkspace {
    fn capacity_bytes(&self) -> usize {
        self.matrix
            .rows
            .iter()
            .map(|row| row.capacity() * std::mem::size_of::<[i32; 3]>())
            .sum::<usize>()
            + self.matrix.previous.capacity()
            + self.matrix.waves.capacity() * std::mem::size_of::<i32>()
            + self.matrix.operations.capacity() * std::mem::size_of::<EditOperation>()
            + self.left_query.capacity()
            + self.left_target.capacity()
    }
}

#[cfg(all(test, target_arch = "x86_64"))]
thread_local! {
    /// Selects the scalar endpoint kernel on AVX2 hosts, so tests can compare both kernels.
    static SCALAR_ENDPOINT_KERNEL: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

#[cfg(test)]
#[derive(Clone, Copy, Debug)]
struct EndpointCell {
    scores: [i32; 3],
    previous: [u8; 3],
}

#[cfg(test)]
impl Default for EndpointCell {
    fn default() -> Self {
        Self {
            scores: [NEGATIVE; 3],
            previous: [START; 3],
        }
    }
}

struct EndpointResult {
    runs: Vec<EditRun>,
    query_bases: usize,
    target_bases: usize,
    matrix_cells: usize,
}

#[allow(clippy::too_many_arguments)]
fn complete_endpoints(
    workspace: &mut EndpointWorkspace,
    mut core: Alignment,
    query: &[u8],
    target: &[u8],
    target_offset: u64,
    max_extension: usize,
    config: AlignmentConfig,
    mut work: Option<&mut AlignmentWork>,
) -> Result<EndpointCompletion, AlignmentError> {
    let forward_target_len = target.len();
    core.validate_cigar()?;
    let query_start =
        usize::try_from(core.query_interval.start).map_err(|_| AlignmentError::LengthOverflow)?;
    let query_end =
        usize::try_from(core.query_interval.end).map_err(|_| AlignmentError::LengthOverflow)?;
    let local_start = core
        .target_interval
        .start
        .checked_sub(target_offset)
        .and_then(|value| usize::try_from(value).ok())
        .ok_or(AlignmentError::EndpointOutsideWindow)?;
    let local_end = core
        .target_interval
        .end
        .checked_sub(target_offset)
        .and_then(|value| usize::try_from(value).ok())
        .ok_or(AlignmentError::EndpointOutsideWindow)?;
    if query_start > query_end
        || query_end > query.len()
        || local_start > local_end
        || local_end > forward_target_len
    {
        return Err(AlignmentError::EndpointOutsideWindow);
    }
    let (target_start, target_end) = match core.strand {
        Strand::Forward => (local_start, local_end),
        Strand::Reverse => (
            forward_target_len - local_end,
            forward_target_len - local_start,
        ),
    };

    let left_query_len = query_start.min(max_extension);
    let left_target_len = target_start.min(max_extension);
    let EndpointWorkspace {
        matrix,
        left_query,
        left_target,
    } = workspace;
    left_query.clear();
    left_query.extend(
        query[query_start - left_query_len..query_start]
            .iter()
            .rev()
            .copied(),
    );
    left_target.clear();
    left_target.extend(
        target[target_start - left_target_len..target_start]
            .iter()
            .rev()
            .copied(),
    );
    if let Some(work) = work.as_deref_mut() {
        work.endpoint_scratch_bytes += (left_query.len() + left_target.len()) as u64;
    }
    let mut left =
        anchored_semiglobal(matrix, left_query, left_target, config, work.as_deref_mut())?;
    left.runs.reverse();

    let query_limit = query_end.saturating_add(max_extension).min(query.len());
    let target_limit = target_end.saturating_add(max_extension).min(target.len());
    let right = anchored_semiglobal(
        matrix,
        &query[query_end..query_limit],
        &target[target_end..target_limit],
        config,
        work,
    )?;

    let mut runs = Vec::with_capacity(left.runs.len() + core.edit_script.len() + right.runs.len());
    append_runs(&mut runs, left.runs)?;
    append_runs(&mut runs, core.edit_script)?;
    append_runs(&mut runs, right.runs)?;
    let summary = summarize_runs(&runs)?;
    core.score = score_runs(&runs, config);
    core.query_interval = Interval::new(
        u64::try_from(query_start - left.query_bases)
            .map_err(|_| AlignmentError::LengthOverflow)?,
        u64::try_from(query_end + right.query_bases).map_err(|_| AlignmentError::LengthOverflow)?,
    )?;
    let oriented_start = target_start - left.target_bases;
    let oriented_end = target_end + right.target_bases;
    let (forward_start, forward_end) = match core.strand {
        Strand::Forward => (oriented_start, oriented_end),
        Strand::Reverse => (
            forward_target_len - oriented_end,
            forward_target_len - oriented_start,
        ),
    };
    core.target_interval = Interval::new(
        target_offset
            .checked_add(u64::try_from(forward_start).map_err(|_| AlignmentError::LengthOverflow)?)
            .ok_or(AlignmentError::LengthOverflow)?,
        target_offset
            .checked_add(u64::try_from(forward_end).map_err(|_| AlignmentError::LengthOverflow)?)
            .ok_or(AlignmentError::LengthOverflow)?,
    )?;
    core.matches = summary.matches;
    core.substitutions = summary.substitutions;
    core.insertions = summary.insertions;
    core.deletions = summary.deletions;
    core.cigar = cigar_from_runs(&runs)?;
    core.edit_script = runs;
    core.validate_cigar()?;

    Ok(EndpointCompletion {
        alignment: core,
        metrics: EndpointMetrics {
            left_query_bases: left.query_bases,
            left_target_bases: left.target_bases,
            right_query_bases: right.query_bases,
            right_target_bases: right.target_bases,
            matrix_cells: left.matrix_cells.saturating_add(right.matrix_cells),
        },
    })
}

/// Keeps the best endpoint cell: higher score, then the later (row + column, row, column) and
/// the lower state. The comparison is a total order, so the visiting order does not matter.
fn endpoint_consider(
    best: &mut (i32, usize, usize, u8),
    row: usize,
    column: usize,
    scores: [i32; 3],
) {
    let (score, state) = maximum(scores);
    let candidate = (
        row.saturating_add(column),
        row,
        column,
        std::cmp::Reverse(state),
    );
    let current = (
        best.1.saturating_add(best.2),
        best.1,
        best.2,
        std::cmp::Reverse(best.3),
    );
    if score > best.0 || (score == best.0 && candidate > current) {
        *best = (score, row, column, state);
    }
}

/// Anchored semiglobal recurrence evaluated along anti-diagonals. A cell reads only the two
/// previous anti-diagonals, the scores, saturation, tie order and traceback bits equal the row
/// kernel, and the best cell comes from the last row and column through `endpoint_consider`.
/// Callers must enable AVX2 and set the border tracebacks.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn endpoint_waves_avx2(
    matrix: &mut EndpointMatrix,
    query: &[u8],
    target: &[u8],
    config: AlignmentConfig,
    best: &mut (i32, usize, usize, u8),
) {
    let (query_len, target_len) = (query.len(), target.len());
    let columns = target_len + 1;
    let stride = query_len + 9;
    matrix.waves.clear();
    matrix.waves.resize(9 * stride, NEGATIVE);
    let open = gap_open(config);
    let border = |length: usize| {
        config.gap_open_score.saturating_add(
            config
                .gap_extend_score
                .saturating_mul(i32::try_from(length).unwrap_or(i32::MAX)),
        )
    };
    let index = |slot: usize, state: u8, row: usize| (slot * 3 + state as usize) * stride + row;
    let mut slots = [0usize, 1, 2];
    for wave in 0..=query_len + target_len {
        let [older, previous, current] = slots;
        let waves = &mut matrix.waves;
        let read = |waves: &Vec<i32>, slot, row| {
            [
                waves[index(slot, MATCH, row)],
                waves[index(slot, INSERTION, row)],
                waves[index(slot, DELETION, row)],
            ]
        };
        let low = wave.saturating_sub(target_len);
        let high = query_len.min(wave);
        let mut row = low;
        while row <= high {
            let column = wave - row;
            if row == 0 || column == 0 {
                let cell = match (row, column) {
                    (0, 0) => [0, NEGATIVE, NEGATIVE],
                    (0, _) => [NEGATIVE, border(column), NEGATIVE],
                    _ => [NEGATIVE, NEGATIVE, border(row)],
                };
                for state in [MATCH, INSERTION, DELETION] {
                    waves[index(current, state, row)] = cell[state as usize];
                }
                row += 1;
                continue;
            }
            let last = high.min(wave - 1);
            if last + 1 - row < 8 {
                let diagonal = read(waves, older, row - 1);
                let above = read(waves, previous, row - 1);
                let left = read(waves, previous, row);
                let (score, match_state) = maximum(diagonal);
                let substitution = if query[row - 1].eq_ignore_ascii_case(&target[column - 1]) {
                    config.match_score
                } else {
                    config.mismatch_score
                };
                let (deletion, deletion_state) = maximum([
                    above[MATCH as usize].saturating_add(open),
                    above[INSERTION as usize].saturating_add(open),
                    above[DELETION as usize].saturating_add(config.gap_extend_score),
                ]);
                let (insertion, insertion_state) = maximum([
                    left[MATCH as usize].saturating_add(open),
                    left[INSERTION as usize].saturating_add(config.gap_extend_score),
                    left[DELETION as usize].saturating_add(open),
                ]);
                waves[index(current, MATCH, row)] = score.saturating_add(substitution);
                waves[index(current, INSERTION, row)] = insertion;
                waves[index(current, DELETION, row)] = deletion;
                matrix.previous[row * columns + column] =
                    match_state | (insertion_state << 2) | (deletion_state << 4);
                row += 1;
                continue;
            }
            // Full blocks, then one block ending at `last`; overlapping cells are rewritten with
            // identical values because they read only earlier anti-diagonals.
            let mut block = row;
            loop {
                let start = block.min(last - 7);
                // SAFETY: rows start..start+8 are interior cells of this wave, so their query
                // bytes, reversed target bytes and wave rows start-1..start+8 are in bounds.
                unsafe {
                    let load = |slot: usize, state: u8, row: usize| {
                        _mm256_loadu_si256(waves.as_ptr().add(index(slot, state, row)).cast())
                    };
                    let query_bytes = _mm_loadl_epi64(query.as_ptr().add(start - 1).cast());
                    let target_bytes =
                        _mm_loadl_epi64(target.as_ptr().add(wave - start - 8).cast());
                    let reverse =
                        _mm_setr_epi8(7, 6, 5, 4, 3, 2, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1);
                    let equal = _mm256_cvtepi8_epi32(_mm_cmpeq_epi8(
                        lowercase_ascii_8(query_bytes),
                        lowercase_ascii_8(_mm_shuffle_epi8(target_bytes, reverse)),
                    ));
                    let substitution = _mm256_blendv_epi8(
                        _mm256_set1_epi32(config.mismatch_score),
                        _mm256_set1_epi32(config.match_score),
                        equal,
                    );
                    let (diagonal, match_state) = choose_avx2(
                        load(older, MATCH, start - 1),
                        load(older, INSERTION, start - 1),
                        load(older, DELETION, start - 1),
                    );
                    let open = _mm256_set1_epi32(open);
                    let extend = _mm256_set1_epi32(config.gap_extend_score);
                    let (deletion, deletion_state) = choose_avx2(
                        saturating_add_avx2(load(previous, MATCH, start - 1), open),
                        saturating_add_avx2(load(previous, INSERTION, start - 1), open),
                        saturating_add_avx2(load(previous, DELETION, start - 1), extend),
                    );
                    let (insertion, insertion_state) = choose_avx2(
                        saturating_add_avx2(load(previous, MATCH, start), open),
                        saturating_add_avx2(load(previous, INSERTION, start), extend),
                        saturating_add_avx2(load(previous, DELETION, start), open),
                    );
                    let scores = waves.as_mut_ptr();
                    _mm256_storeu_si256(
                        scores.add(index(current, MATCH, start)).cast(),
                        saturating_add_avx2(diagonal, substitution),
                    );
                    _mm256_storeu_si256(
                        scores.add(index(current, INSERTION, start)).cast(),
                        insertion,
                    );
                    _mm256_storeu_si256(
                        scores.add(index(current, DELETION, start)).cast(),
                        deletion,
                    );
                    let packed = _mm256_or_si256(
                        match_state,
                        _mm256_or_si256(
                            _mm256_slli_epi32::<2>(insertion_state),
                            _mm256_slli_epi32::<4>(deletion_state),
                        ),
                    );
                    let mut lanes = [0i32; 8];
                    _mm256_storeu_si256(lanes.as_mut_ptr().cast(), packed);
                    let first = start * columns + wave - start;
                    for (lane, &value) in lanes.iter().enumerate() {
                        matrix.previous[first + lane * (columns - 1)] = value as u8;
                    }
                }
                if start + 7 == last {
                    break;
                }
                block += 8;
            }
            row = last + 1;
        }
        let waves = &matrix.waves;
        if let Some(row) = wave.checked_sub(target_len).filter(|&row| row <= query_len) {
            endpoint_consider(best, row, target_len, read(waves, current, row));
        }
        if let Some(column) = wave
            .checked_sub(query_len)
            .filter(|&column| column <= target_len)
        {
            endpoint_consider(best, query_len, column, read(waves, current, query_len));
        }
        slots = [previous, current, older];
    }
}

fn anchored_semiglobal(
    matrix: &mut EndpointMatrix,
    query: &[u8],
    target: &[u8],
    config: AlignmentConfig,
    mut work: Option<&mut AlignmentWork>,
) -> Result<EndpointResult, AlignmentError> {
    if query.is_empty() || target.is_empty() {
        return Ok(EndpointResult {
            runs: Vec::new(),
            query_bases: 0,
            target_bases: 0,
            matrix_cells: 0,
        });
    }
    let rows = query
        .len()
        .checked_add(1)
        .ok_or(AlignmentError::LengthOverflow)?;
    let columns = target
        .len()
        .checked_add(1)
        .ok_or(AlignmentError::LengthOverflow)?;
    let matrix_cells = rows
        .checked_mul(columns)
        .ok_or(AlignmentError::LengthOverflow)?;
    if matrix_cells > config.max_cells {
        return Err(AlignmentError::MatrixTooLarge {
            cells: matrix_cells,
            max_cells: config.max_cells,
        });
    }
    let init_cpu = observed_cpu(work.is_some());
    if let Some(work) = work.as_deref_mut() {
        work.endpoint_passes += 1;
        work.endpoint_cells += matrix_cells as u64;
        work.endpoint_recurrence_cells += (query.len() * target.len()) as u64;
        work.growth_operations += u64::from(matrix.previous.capacity() < matrix_cells)
            + matrix
                .rows
                .iter()
                .filter(|row| row.capacity() < columns)
                .count() as u64;
    }
    // Two bits per predecessor; the high two bits are reserved zero.
    const INITIAL_PREVIOUS: u8 = START | (START << 2) | (START << 4);
    matrix.previous.resize(matrix_cells, INITIAL_PREVIOUS);
    matrix.previous.fill(INITIAL_PREVIOUS);
    for row in &mut matrix.rows {
        row.resize(columns, [NEGATIVE; 3]);
    }
    matrix.rows[0][0] = [0, NEGATIVE, NEGATIVE];
    for column in 1..columns {
        matrix.rows[0][column] = [NEGATIVE; 3];
        matrix.rows[0][column][INSERTION as usize] = config.gap_open_score.saturating_add(
            config
                .gap_extend_score
                .saturating_mul(i32::try_from(column).unwrap_or(i32::MAX)),
        );
        let state = if column == 1 { MATCH } else { INSERTION };
        matrix.previous[column] = (INITIAL_PREVIOUS & !(3 << 2)) | (state << 2);
    }
    if let Some(work) = work.as_deref_mut() {
        work.endpoint_init_cpu_ns += elapsed_cpu(init_cpu);
    }
    let matrix_cpu = observed_cpu(work.is_some());
    let mut best = (NEGATIVE, 0usize, 0usize, START);
    #[cfg(target_arch = "x86_64")]
    let vector = is_x86_feature_detected!("avx2");
    #[cfg(all(test, target_arch = "x86_64"))]
    let vector = vector && !SCALAR_ENDPOINT_KERNEL.get();
    #[cfg(not(target_arch = "x86_64"))]
    let vector = false;
    #[cfg(target_arch = "x86_64")]
    if vector {
        for row in 1..rows {
            let state = if row == 1 { MATCH } else { DELETION };
            matrix.previous[row * columns] = (INITIAL_PREVIOUS & !(3 << 4)) | (state << 4);
        }
        // SAFETY: the runtime feature check guards every AVX2 instruction in this call.
        unsafe { endpoint_waves_avx2(matrix, query, target, config, &mut best) };
    }
    if !vector {
        endpoint_consider(&mut best, 0, columns - 1, matrix.rows[0][columns - 1]);
    }
    for row in (1..rows).filter(|_| !vector) {
        let first = row * columns;
        let [previous, current] = &mut matrix.rows;
        let previous = &previous[..columns];
        let current = &mut current[..columns];
        let traceback = &mut matrix.previous[first..first + columns];
        current[0] = [NEGATIVE; 3];
        current[0][DELETION as usize] = config.gap_open_score.saturating_add(
            config
                .gap_extend_score
                .saturating_mul(i32::try_from(row).unwrap_or(i32::MAX)),
        );
        let state = if row == 1 { MATCH } else { DELETION };
        traceback[0] = (INITIAL_PREVIOUS & !(3 << 4)) | (state << 4);
        for column in 1..columns {
            let (score, state) = maximum(previous[column - 1]);
            let left = current[column - 1];
            let cell = &mut current[column];
            cell[MATCH as usize] = score.saturating_add(
                if query[row - 1].eq_ignore_ascii_case(&target[column - 1]) {
                    config.match_score
                } else {
                    config.mismatch_score
                },
            );
            let mut packed = state;

            let above = previous[column];
            let (score, state) = maximum([
                above[MATCH as usize].saturating_add(gap_open(config)),
                above[INSERTION as usize].saturating_add(gap_open(config)),
                above[DELETION as usize].saturating_add(config.gap_extend_score),
            ]);
            cell[DELETION as usize] = score;
            packed |= state << 4;

            let (score, state) = maximum([
                left[MATCH as usize].saturating_add(gap_open(config)),
                left[INSERTION as usize].saturating_add(config.gap_extend_score),
                left[DELETION as usize].saturating_add(gap_open(config)),
            ]);
            cell[INSERTION as usize] = score;
            traceback[column] = packed | (state << 2);
        }
        if row + 1 == rows {
            for (column, &scores) in current.iter().enumerate() {
                endpoint_consider(&mut best, row, column, scores);
            }
        } else {
            endpoint_consider(&mut best, row, columns - 1, current[columns - 1]);
        }
        matrix.rows.swap(0, 1);
    }

    if best.3 > DELETION {
        return Err(AlignmentError::NoAlignment);
    }

    if let Some(work) = work.as_deref_mut() {
        work.endpoint_matrix_cpu_ns += elapsed_cpu(matrix_cpu);
    }
    let trace_cpu = observed_cpu(work.is_some());
    let mut row = best.1;
    let mut column = best.2;
    let mut state = best.3;
    let operations = &mut matrix.operations;
    operations.clear();
    operations.reserve(row.saturating_add(column));
    while row > 0 || column > 0 {
        let index = row * columns + column;
        let packed = matrix.previous[index];
        debug_assert_eq!(packed & 0xc0, 0);
        let previous = (packed >> (2 * state)) & 3;
        match state {
            MATCH if row > 0 && column > 0 => {
                operations.push(
                    if query[row - 1].eq_ignore_ascii_case(&target[column - 1]) {
                        EditOperation::Equal
                    } else {
                        EditOperation::Substitution
                    },
                );
                row -= 1;
                column -= 1;
            }
            INSERTION if column > 0 => {
                operations.push(EditOperation::Insertion);
                column -= 1;
            }
            DELETION if row > 0 => {
                operations.push(EditOperation::Deletion);
                row -= 1;
            }
            _ => return Err(AlignmentError::InvalidTraceback),
        }
        state = previous;
    }
    operations.reverse();
    let runs = runs_from_operations(operations)?;
    if let Some(work) = work {
        work.endpoint_trace_cpu_ns += elapsed_cpu(trace_cpu);
        work.endpoint_scratch_bytes += operations.len() as u64;
    }
    Ok(EndpointResult {
        runs,
        query_bases: best.1,
        target_bases: best.2,
        matrix_cells,
    })
}

#[cfg(test)]
fn anchored_semiglobal_reference(
    cells: &mut Vec<EndpointCell>,
    query: &[u8],
    target: &[u8],
    config: AlignmentConfig,
) -> Result<EndpointResult, AlignmentError> {
    if query.is_empty() || target.is_empty() {
        return Ok(EndpointResult {
            runs: Vec::new(),
            query_bases: 0,
            target_bases: 0,
            matrix_cells: 0,
        });
    }
    let rows = query
        .len()
        .checked_add(1)
        .ok_or(AlignmentError::LengthOverflow)?;
    let columns = target
        .len()
        .checked_add(1)
        .ok_or(AlignmentError::LengthOverflow)?;
    let matrix_cells = rows
        .checked_mul(columns)
        .ok_or(AlignmentError::LengthOverflow)?;
    if matrix_cells > config.max_cells {
        return Err(AlignmentError::MatrixTooLarge {
            cells: matrix_cells,
            max_cells: config.max_cells,
        });
    }
    cells.resize(matrix_cells, EndpointCell::default());
    cells.fill(EndpointCell::default());
    cells[0].scores[MATCH as usize] = 0;
    for (column, cell) in cells.iter_mut().enumerate().take(columns).skip(1) {
        cell.scores[INSERTION as usize] = config.gap_open_score.saturating_add(
            config
                .gap_extend_score
                .saturating_mul(i32::try_from(column).unwrap_or(i32::MAX)),
        );
        cell.previous[INSERTION as usize] = if column == 1 { MATCH } else { INSERTION };
    }
    for row in 1..rows {
        let first = row * columns;
        let (previous, current) = cells.split_at_mut(first);
        let previous = &previous[first - columns..first];
        let current = &mut current[..columns];
        current[0].scores[DELETION as usize] = config.gap_open_score.saturating_add(
            config
                .gap_extend_score
                .saturating_mul(i32::try_from(row).unwrap_or(i32::MAX)),
        );
        current[0].previous[DELETION as usize] = if row == 1 { MATCH } else { DELETION };
        for column in 1..columns {
            let (score, state) = maximum(previous[column - 1].scores);
            let left = current[column - 1].scores;
            let cell = &mut current[column];
            cell.scores[MATCH as usize] = score.saturating_add(
                if query[row - 1].eq_ignore_ascii_case(&target[column - 1]) {
                    config.match_score
                } else {
                    config.mismatch_score
                },
            );
            cell.previous[MATCH as usize] = state;

            let above = previous[column].scores;
            let (score, state) = maximum([
                above[MATCH as usize].saturating_add(gap_open(config)),
                above[INSERTION as usize].saturating_add(gap_open(config)),
                above[DELETION as usize].saturating_add(config.gap_extend_score),
            ]);
            cell.scores[DELETION as usize] = score;
            cell.previous[DELETION as usize] = state;

            let (score, state) = maximum([
                left[MATCH as usize].saturating_add(gap_open(config)),
                left[INSERTION as usize].saturating_add(config.gap_extend_score),
                left[DELETION as usize].saturating_add(gap_open(config)),
            ]);
            cell.scores[INSERTION as usize] = score;
            cell.previous[INSERTION as usize] = state;
        }
    }

    let mut best = (NEGATIVE, 0usize, 0usize, START);
    for row in 0..rows {
        let first_column = if row + 1 == rows { 0 } else { columns - 1 };
        for column in first_column..columns {
            if (row + 1 != rows && column + 1 != columns) || (row == 0 && column == 0) {
                continue;
            }
            let (score, state) = maximum(cells[row * columns + column].scores);
            let candidate = (
                row.saturating_add(column),
                row,
                column,
                std::cmp::Reverse(state),
            );
            let current = (
                best.1.saturating_add(best.2),
                best.1,
                best.2,
                std::cmp::Reverse(best.3),
            );
            if score > best.0 || (score == best.0 && candidate > current) {
                best = (score, row, column, state);
            }
        }
    }
    if best.3 > DELETION {
        return Err(AlignmentError::NoAlignment);
    }

    let mut row = best.1;
    let mut column = best.2;
    let mut state = best.3;
    let mut operations = Vec::with_capacity(row.saturating_add(column));
    while row > 0 || column > 0 {
        let index = row * columns + column;
        let previous = cells[index].previous[state as usize];
        match state {
            MATCH if row > 0 && column > 0 => {
                operations.push(
                    if query[row - 1].eq_ignore_ascii_case(&target[column - 1]) {
                        EditOperation::Equal
                    } else {
                        EditOperation::Substitution
                    },
                );
                row -= 1;
                column -= 1;
            }
            INSERTION if column > 0 => {
                operations.push(EditOperation::Insertion);
                column -= 1;
            }
            DELETION if row > 0 => {
                operations.push(EditOperation::Deletion);
                row -= 1;
            }
            _ => return Err(AlignmentError::InvalidTraceback),
        }
        state = previous;
    }
    operations.reverse();
    Ok(EndpointResult {
        runs: runs_from_operations(&operations)?,
        query_bases: best.1,
        target_bases: best.2,
        matrix_cells,
    })
}

fn maximum(values: [i32; 3]) -> (i32, u8) {
    values
        .into_iter()
        .enumerate()
        .skip(1)
        .fold((values[0], MATCH), |best, (state, score)| {
            if score > best.0 {
                (score, state as u8)
            } else {
                best
            }
        })
}

fn append_runs(output: &mut Vec<EditRun>, runs: Vec<EditRun>) -> Result<(), AlignmentError> {
    for run in runs {
        if let Some(last) = output.last_mut()
            && last.operation == run.operation
        {
            last.length = last
                .length
                .checked_add(run.length)
                .ok_or(AlignmentError::RunTooLong)?;
        } else {
            output.push(run);
        }
    }
    Ok(())
}

fn score_runs(runs: &[EditRun], config: AlignmentConfig) -> i32 {
    runs.iter().fold(0, |score, run| {
        let length = i32::try_from(run.length).unwrap_or(i32::MAX);
        let delta = match run.operation {
            EditOperation::Equal => config.match_score.saturating_mul(length),
            EditOperation::Substitution => config.mismatch_score.saturating_mul(length),
            EditOperation::Insertion | EditOperation::Deletion => config
                .gap_open_score
                .saturating_add(config.gap_extend_score.saturating_mul(length)),
        };
        score.saturating_add(delta)
    })
}

fn finish(
    raw: RawAlignment,
    strand: Strand,
    target_offset: u64,
    target_len: usize,
) -> Result<Alignment, AlignmentError> {
    let target_len = u64::try_from(target_len).map_err(|_| AlignmentError::LengthOverflow)?;
    let (start, end) = if strand == Strand::Reverse {
        (
            target_len
                .checked_sub(raw.target_interval.end)
                .ok_or(AlignmentError::LengthOverflow)?,
            target_len
                .checked_sub(raw.target_interval.start)
                .ok_or(AlignmentError::LengthOverflow)?,
        )
    } else {
        (raw.target_interval.start, raw.target_interval.end)
    };
    let alignment = Alignment {
        score: raw.score,
        strand,
        query_interval: raw.query_interval,
        target_interval: Interval::new(
            target_offset
                .checked_add(start)
                .ok_or(AlignmentError::LengthOverflow)?,
            target_offset
                .checked_add(end)
                .ok_or(AlignmentError::LengthOverflow)?,
        )?,
        matches: raw.summary.matches,
        substitutions: raw.summary.substitutions,
        insertions: raw.summary.insertions,
        deletions: raw.summary.deletions,
        cigar: raw.cigar,
        edit_script: raw.edit_script,
    };
    alignment.validate_cigar()?;
    Ok(alignment)
}

fn choose<const N: usize>(values: [(i32, u8); N]) -> (i32, u8) {
    values.into_iter().fold((0, UNREACHABLE), |best, value| {
        if value.0 > best.0 || (value.0 == best.0 && value.1 < best.1) {
            value
        } else {
            best
        }
    })
}

fn gap_open(config: AlignmentConfig) -> i32 {
    config
        .gap_open_score
        .saturating_add(config.gap_extend_score)
}

fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        _ => b'N',
    }
}

fn runs_from_operations(operations: &[EditOperation]) -> Result<Vec<EditRun>, AlignmentError> {
    let mut runs: Vec<EditRun> = Vec::new();
    for &operation in operations {
        if let Some(run) = runs.last_mut()
            && run.operation == operation
        {
            run.length = run
                .length
                .checked_add(1)
                .ok_or(AlignmentError::RunTooLong)?;
        } else {
            runs.push(EditRun {
                operation,
                length: 1,
            });
        }
    }
    Ok(runs)
}

pub fn parse_cigar(cigar: &str) -> Result<Vec<EditRun>, AlignmentError> {
    let mut runs = Vec::new();
    let mut length = 0u32;
    for byte in cigar.bytes() {
        if byte.is_ascii_digit() {
            length = length
                .checked_mul(10)
                .and_then(|value| value.checked_add(u32::from(byte - b'0')))
                .ok_or(AlignmentError::RunTooLong)?;
            continue;
        }
        if length == 0 {
            return Err(AlignmentError::InvalidCigar);
        }
        let operation = match byte {
            b'=' => EditOperation::Equal,
            b'X' => EditOperation::Substitution,
            b'I' => EditOperation::Insertion,
            b'D' => EditOperation::Deletion,
            _ => return Err(AlignmentError::InvalidCigar),
        };
        if runs
            .last()
            .is_some_and(|run: &EditRun| run.operation == operation)
        {
            return Err(AlignmentError::InvalidCigar);
        }
        runs.push(EditRun { operation, length });
        length = 0;
    }
    if length != 0 || runs.is_empty() {
        return Err(AlignmentError::InvalidCigar);
    }
    Ok(runs)
}

pub(crate) fn cigar_from_runs(runs: &[EditRun]) -> Result<String, AlignmentError> {
    let mut cigar = String::new();
    for run in runs {
        let operation = match run.operation {
            EditOperation::Equal => '=',
            EditOperation::Substitution => 'X',
            EditOperation::Insertion => 'I',
            EditOperation::Deletion => 'D',
        };
        write!(&mut cigar, "{}{operation}", run.length).map_err(|_| AlignmentError::CigarWrite)?;
    }
    Ok(cigar)
}

#[derive(Clone, Copy, Debug, Default)]
struct RunSummary {
    query_bases: u64,
    target_bases: u64,
    matches: u64,
    substitutions: u64,
    insertions: u64,
    deletions: u64,
}

fn summarize_runs(runs: &[EditRun]) -> Result<RunSummary, AlignmentError> {
    let mut summary = RunSummary::default();
    for run in runs {
        if run.length == 0 {
            return Err(AlignmentError::InvalidCigar);
        }
        let length = u64::from(run.length);
        match run.operation {
            EditOperation::Equal => {
                summary.query_bases += length;
                summary.target_bases += length;
                summary.matches += length;
            }
            EditOperation::Substitution => {
                summary.query_bases += length;
                summary.target_bases += length;
                summary.substitutions += length;
            }
            EditOperation::Insertion => {
                summary.target_bases += length;
                summary.insertions += length;
            }
            EditOperation::Deletion => {
                summary.query_bases += length;
                summary.deletions += length;
            }
        }
    }
    Ok(summary)
}

#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum AlignmentError {
    #[error("alignment query is empty")]
    EmptyQuery,
    #[error("alignment target is empty")]
    EmptyTarget,
    #[error("alignment configuration is invalid")]
    InvalidConfig,
    #[error("alignment length overflows")]
    LengthOverflow,
    #[error("alignment band excludes the input")]
    BandExcludesInput,
    #[error("alignment matrix requires {cells} cells, exceeding {max_cells}")]
    MatrixTooLarge { cells: usize, max_cells: usize },
    #[error("no positive-scoring local alignment was found")]
    NoAlignment,
    #[error("alignment traceback left the band")]
    TracebackOutsideBand,
    #[error("alignment traceback is inconsistent")]
    InvalidTraceback,
    #[error("alignment interval is reversed: {start}..{end}")]
    ReversedInterval { start: u64, end: u64 },
    #[error("alignment edit run exceeds u32")]
    RunTooLong,
    #[error("CIGAR text is invalid")]
    InvalidCigar,
    #[error("CIGAR text differs from the edit script")]
    CigarMismatch,
    #[error("CIGAR span differs from alignment coordinates")]
    CigarSpanMismatch,
    #[error("CIGAR counts differ from alignment counts")]
    CigarCountMismatch,
    #[error("failed to write CIGAR text")]
    CigarWrite,
    #[error("alignment core is outside the supplied endpoint window")]
    EndpointOutsideWindow,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config() -> AlignmentConfig {
        AlignmentConfig {
            band_width: 8,
            ..AlignmentConfig::default()
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn local_lane_counts_and_narrow_bounds() {
        assert!(narrow_local_scores(6_247, 41_651, config()));
        assert!(narrow_local_scores(16_382, 20_000, config()));
        assert!(!narrow_local_scores(16_383, 20_000, config()));
        for changed in [
            AlignmentConfig {
                mismatch_score: i32::MIN,
                ..config()
            },
            AlignmentConfig {
                match_score: i32::MAX,
                ..config()
            },
            AlignmentConfig {
                gap_open_score: i32::MIN,
                ..config()
            },
        ] {
            assert!(!narrow_local_scores(32, 32, changed));
        }
        if !is_x86_feature_detected!("avx2") {
            return;
        }
        for length in [7, 8, 15, 16, 17, 63] {
            for diagonal_offset in [-17, 0, 17] {
                let cfg = AlignmentConfig {
                    band_width: 16,
                    diagonal_offset,
                    ..config()
                };
                let mut workspace = AlignmentWorkspace::default();
                workspace.enable_timing();
                let _ = workspace.align(&vec![b'A'; length], &vec![b'A'; length + 5], cfg);
                let work = workspace.work;
                let mut slots = 0;
                let mut boundary = 0;
                for q in 0..=length {
                    for t in 0..=length + 5 {
                        if ((t as i64 - q as i64) - diagonal_offset).abs() <= 16 {
                            slots += 1;
                            boundary += u64::from(q == 0 || t == 0);
                        }
                    }
                }
                assert_eq!(work.local_cells, slots);
                assert_eq!(work.local_boundary_cells, boundary);
                assert_eq!(
                    work.local_vector8_cells + work.local_vector16_cells + work.local_scalar_cells,
                    slots
                );
            }
        }
    }

    #[test]
    fn band_cells_match_kernel_row_and_wave_layouts() {
        let brute = |m: usize, n: usize, d: i64, w: u32| {
            (0..=m)
                .flat_map(|i| (0..=n).map(move |j| (i, j)))
                .filter(|&(i, j)| (j as i64 - i as i64 - d).abs() <= i64::from(w))
                .count()
        };
        let cases: [(usize, usize, i64, u32); 16] = [
            (0, 0, 0, 0),
            (1, 1, 0, 0),
            (40, 3, 0, 4),
            (3, 40, 0, 4),
            (60, 60, -45, 8),
            (60, 60, 45, 8),
            (60, 20, 30, 8),
            (20, 20, 29, 8),
            (20, 20, -29, 8),
            (20, 20, 28, 8),
            (20, 20, -28, 8),
            (500, 7, -3, 128),
            (7, 500, 3, 128),
            (300, 280, -17, 16),
            (900, 120, -700, 64),
            (120, 900, 700, 64),
        ];
        for (m, n, d, w) in cases {
            let expected = brute(m, n, d, w);
            assert_eq!(band_cells(m, n, d, w).unwrap(), expected, "{m} {n} {d} {w}");
            let cfg = AlignmentConfig {
                band_width: w,
                diagonal_offset: d,
                ..config()
            };
            let mut workspace = AlignmentWorkspace::default();
            workspace.prepare_rows(m, n, cfg).unwrap();
            for i in 0..=m {
                let row = band_row(i, n, d, w).unwrap();
                assert_eq!(
                    (workspace.row_starts[i], workspace.row_widths[i]),
                    row.map_or((0, 0), |(start, end)| (start, end - start + 1))
                );
            }
            assert_eq!(workspace.row_offsets[m] + workspace.row_widths[m], expected);
            if m > 0 && n > 0 && expected > 0 {
                let (query, target) = (vec![b'A'; m], vec![b'C'; n]);
                let mut scalar = AlignmentWorkspace::default();
                scalar.enable_timing();
                let _ = scalar.align_raw_scalar(&query, &target, cfg);
                assert_eq!(scalar.work.local_cells as usize, expected);
                #[cfg(target_arch = "x86_64")]
                if is_x86_feature_detected!("avx2") {
                    let mut avx2 = AlignmentWorkspace::default();
                    avx2.enable_timing();
                    // SAFETY: guarded by the runtime AVX2 check.
                    let _ = unsafe { avx2.align_raw_avx2(&query, &target, cfg) };
                    assert_eq!(avx2.work.local_cells as usize, expected);
                }
            }
            #[cfg(target_arch = "x86_64")]
            assert_eq!(
                (0..=m + n)
                    .filter_map(|wave| wave_range(m, n, cfg, wave))
                    .map(WaveRange::width)
                    .sum::<usize>(),
                expected
            );
        }
        // Rows outside a short target leave the exact count well below rows times band.
        let (m, n, d, w) = (61_200, 50_000, -5_000, 128);
        let conservative = (m + 1) * (2 * w as usize + 1);
        let exact = band_cells(m, n, d, w).unwrap();
        assert!(conservative > 15_728_640 && exact <= 15_728_640, "{exact}");
        assert_eq!(
            exact,
            (0..=m)
                .map(|i| {
                    let low = (i as i64 + d - i64::from(w)).max(0);
                    let high = (i as i64 + d + i64::from(w)).min(n as i64);
                    (high - low + 1).max(0) as usize
                })
                .sum::<usize>()
        );
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn runtime_dispatch_matches_scalar_and_forced_avx2() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }
        let cases: [(&[u8], &[u8], i64); 4] = [
            (b"ACGTACG", b"ACGTTCG", 0),
            (b"ACGTACGT", b"ACGTTGGT", 0),
            (b"ACGTACGTA", b"ACGTTACGTA", 2),
            (b"TTACGTACG", b"ACGTACG", -2),
        ];
        for (query, target, diagonal_offset) in cases {
            let config = AlignmentConfig {
                band_width: 16,
                diagonal_offset,
                max_cells: 10_000,
                ..AlignmentConfig::default()
            };
            for strand in [Strand::Forward, Strand::Reverse] {
                let oriented_target = match strand {
                    Strand::Forward => target.to_vec(),
                    Strand::Reverse => target.iter().rev().map(|base| complement(*base)).collect(),
                };
                let mut dispatched = AlignmentWorkspace::default();
                let actual = dispatched.align_oriented(query, target, 17, strand, config);
                let mut scalar = AlignmentWorkspace::default();
                let expected = scalar
                    .align_raw_scalar(query, &oriented_target, config)
                    .and_then(|raw| finish(raw, strand, 17, target.len()));
                let mut avx2 = AlignmentWorkspace::default();
                // SAFETY: this test returns early unless AVX2 is available.
                let forced = unsafe { avx2.align_raw_avx2(query, &oriented_target, config) }
                    .and_then(|raw| finish(raw, strand, 17, target.len()));
                assert_eq!(actual, expected);
                assert_eq!(forced, expected);
            }
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn resident_chunks_match_unbounded_alignment() {
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = move |bound: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % bound
        };
        let mut chunked_cases = 0;
        for case in 0..128 {
            let alphabet = if case % 3 == 0 { 2 } else { 4 };
            let mut query: Vec<u8> = (0..24 + next(200))
                .map(|_| b"ACGT"[next(alphabet) as usize])
                .collect();
            let mut target: Vec<u8> = (0..next(24)).map(|_| b"ACGT"[next(4) as usize]).collect();
            for &base in &query {
                match next(24) {
                    0 => {}
                    1 => {
                        target.push(base);
                        target.extend((0..1 + next(9)).map(|_| b"ACGT"[next(4) as usize]));
                    }
                    2 => target.push(b"ACGT"[next(4) as usize]),
                    _ => target.push(base),
                }
            }
            if case % 7 == 0 {
                query.reverse();
            }
            let band_width = 2 + next(15) as u32;
            let base = AlignmentConfig {
                band_width,
                diagonal_offset: next(11) as i64 - 5,
                max_cells: usize::MAX,
                ..AlignmentConfig::default()
            };
            let expected = AlignmentWorkspace::default().align(&query, &target, base);
            let mut reference = AlignmentWorkspace::default();
            let scalar = reference
                .align_raw_scalar(&query, &target, base)
                .and_then(|raw| finish(raw, Strand::Forward, 0, target.len()));
            assert_eq!(scalar, expected);
            let total =
                band_cells(query.len(), target.len(), base.diagonal_offset, band_width).unwrap();
            let smallest = band_width as usize + 1;
            for max_cells in [smallest, smallest + 7, total / 3 + 1, total.max(2) - 1] {
                let config = AlignmentConfig {
                    max_cells: max_cells.max(smallest),
                    ..base
                };
                let mut chunked = AlignmentWorkspace {
                    resident_chunks: true,
                    ..AlignmentWorkspace::default()
                };
                chunked.enable_timing();
                assert_eq!(
                    chunked.align(&query, &target, config),
                    expected,
                    "case {case} max_cells {}",
                    config.max_cells
                );
                let mut portable = AlignmentWorkspace {
                    resident_chunks: true,
                    ..AlignmentWorkspace::default()
                };
                let portable_result = portable
                    .align_raw_waves_scalar(&query, &target, config)
                    .and_then(|raw| finish(raw, Strand::Forward, 0, target.len()));
                assert_eq!(portable_result, expected);
                if total > config.max_cells {
                    assert!(chunked.work.local_chunks >= 2);
                    assert_eq!(chunked.work.local_chunked_passes, 1);
                    assert!(chunked.work.local_recomputed_cells <= total as u64);
                    assert!(chunked.compact_cells.capacity() <= (2 * config.max_cells).max(8));
                    assert!(matches!(
                        AlignmentWorkspace::default().align(&query, &target, config),
                        Err(AlignmentError::MatrixTooLarge { .. })
                    ));
                    chunked_cases += 1;
                }
            }
        }
        assert!(chunked_cases > 250, "{chunked_cases}");
    }

    /// Rows per wave that the AVX2 driver may fill in eight-lane blocks: interior cells whose
    /// diagonal, left and above cells are inside the band. Found by brute force over cells.
    #[cfg(target_arch = "x86_64")]
    fn vector_wave_widths(
        query_len: usize,
        target_len: usize,
        config: AlignmentConfig,
    ) -> Vec<usize> {
        let band = |query: usize, target: usize| {
            query <= query_len
                && target <= target_len
                && (target as i64 - query as i64 - config.diagonal_offset).abs()
                    <= i64::from(config.band_width)
        };
        let (low, high) = (
            config.diagonal_offset - i64::from(config.band_width),
            config.diagonal_offset + i64::from(config.band_width),
        );
        (0..=query_len + target_len)
            .map(|wave| {
                // Only rows with low <= wave - 2 * row <= high can be inside the band.
                let first = ((wave as i64 - high) / 2 - 1).max(1) as usize;
                let last = ((wave as i64 - low) / 2 + 1).clamp(0, query_len as i64) as usize;
                (first..=last.min(wave.saturating_sub(1)))
                    .filter(|&query| {
                        let target = wave - query;
                        band(query, target)
                            && band(query - 1, target - 1)
                            && band(query, target - 1)
                            && band(query - 1, target)
                    })
                    .count()
            })
            .collect()
    }

    /// Checks the scalar and AVX2 wave kernels against the unbounded scalar row kernel: equal
    /// alignments (best-cell tie order included), equal traceback cells, and vector counters
    /// that count each cell once although a final block may evaluate some lanes twice.
    #[cfg(target_arch = "x86_64")]
    fn assert_wave_kernels_match(
        query: &[u8],
        target: &[u8],
        config: AlignmentConfig,
    ) -> (Result<Alignment, AlignmentError>, Vec<usize>) {
        let forward = |raw: Result<RawAlignment, AlignmentError>| {
            raw.and_then(|raw| finish(raw, Strand::Forward, 0, target.len()))
        };
        let unbounded = AlignmentConfig {
            max_cells: usize::MAX,
            ..config
        };
        let expected =
            forward(AlignmentWorkspace::default().align_raw_scalar(query, target, unbounded));
        let [(scalar_result, scalar), (vector_result, vector)] = [false, true].map(|avx2| {
            let mut workspace = AlignmentWorkspace {
                resident_chunks: true,
                ..AlignmentWorkspace::default()
            };
            workspace.enable_timing();
            let result = if avx2 {
                // SAFETY: callers return early unless AVX2 is available.
                unsafe { workspace.align_raw_avx2(query, target, config) }
            } else {
                workspace.align_raw_waves_scalar(query, target, config)
            };
            (forward(result), workspace)
        });
        let case = format!("{config:?} query={query:?} target={target:?}");
        assert_eq!(scalar_result, expected, "{case}");
        assert_eq!(vector_result, expected, "{case}");
        assert_eq!(
            vector.compact_cells.len(),
            scalar.compact_cells.len(),
            "{case}"
        );
        let mismatch = vector
            .compact_cells
            .iter()
            .zip(&scalar.compact_cells)
            .position(|(vector, scalar)| vector != scalar);
        assert_eq!(mismatch, None, "{case}");

        let widths = vector_wave_widths(query.len(), target.len(), config);
        let blocked = widths.iter().filter(|&&width| width >= 8).sum::<usize>() as u64;
        let total = band_cells(
            query.len(),
            target.len(),
            config.diagonal_offset,
            config.band_width,
        )
        .unwrap() as u64;
        for work in [scalar.work, vector.work] {
            assert_eq!(work.local_cells, total, "{case}");
        }
        assert_eq!(vector.work.local_vector8_cells, blocked, "{case}");
        assert_eq!(
            vector.work.local_scalar_cells,
            vector.work.local_cells - blocked,
            "{case}"
        );
        assert_eq!(scalar.work.local_vector8_cells, 0);
        assert_eq!(scalar.work.local_scalar_cells, scalar.work.local_cells);
        assert_eq!(vector.work.local_chunks, scalar.work.local_chunks, "{case}");
        assert_eq!(
            vector.work.local_recomputed_cells, scalar.work.local_recomputed_cells,
            "{case}"
        );
        (expected, widths)
    }

    #[cfg(target_arch = "x86_64")]
    fn tie_and_saturation_configs() -> [AlignmentConfig; 3] {
        [
            AlignmentConfig::default(),
            AlignmentConfig {
                match_score: 1,
                mismatch_score: 0,
                gap_open_score: 0,
                gap_extend_score: 0,
                ..AlignmentConfig::default()
            },
            AlignmentConfig {
                match_score: i32::MAX,
                mismatch_score: i32::MIN,
                gap_open_score: i32::MIN,
                gap_extend_score: i32::MIN,
                ..AlignmentConfig::default()
            },
        ]
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn wave_tail_blocks_match_scalar_at_lane_boundaries_with_ties_and_saturation() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }
        let mut state = 0x3c6e_f372_fe94_f82bu64;
        let mut next = move |bound: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state as usize % bound
        };
        let mut seen = std::collections::BTreeSet::new();
        for config in tie_and_saturation_configs() {
            for band_width in [0, 3, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 23] {
                // Offsets clip waves at the matrix edges, leave early waves empty, or exclude the
                // whole input.
                for diagonal_offset in [-6, 0, 11, 80] {
                    for (query_len, target_len) in [(40, 48), (13, 61), (57, 9)] {
                        let random = (0..query_len + target_len)
                            .map(|_| b"ACGTacgt"[next(8)])
                            .collect::<Vec<_>>();
                        let uniform = [b'A'; 128];
                        for bases in [&random[..], &uniform[..]] {
                            let (query, target) = bases.split_at(query_len);
                            let config = AlignmentConfig {
                                band_width,
                                diagonal_offset,
                                ..config
                            };
                            let (_, widths) =
                                assert_wave_kernels_match(query, &target[..target_len], config);
                            seen.extend(widths);
                        }
                    }
                }
            }
        }
        for width in [0, 7, 8, 9, 15, 16, 17] {
            assert!(seen.contains(&width), "wave width {width} not covered");
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn chunked_wave_tails_match_scalar_across_carries_with_ties_and_saturation() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }
        let mut state = 0xa54f_f53a_5f1d_36f1u64;
        let mut next = move |bound: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state as usize % bound
        };
        for config in tie_and_saturation_configs() {
            for band_width in [7, 8, 9, 15, 16, 17] {
                for diagonal_offset in [-3, 0, 4] {
                    let query = (0..64).map(|_| b"ACGt"[next(4)]).collect::<Vec<_>>();
                    let mut target = Vec::new();
                    for &base in &query {
                        match next(12) {
                            0 => target.extend([base, b"ACGT"[next(4)]]),
                            1 => {}
                            2 => target.push(b"ACGT"[next(4)]),
                            _ => target.push(base),
                        }
                    }
                    let base = AlignmentConfig {
                        band_width,
                        diagonal_offset,
                        ..config
                    };
                    let total =
                        band_cells(query.len(), target.len(), diagonal_offset, band_width).unwrap();
                    let widest = (0..=query.len() + target.len())
                        .filter_map(|wave| wave_range(query.len(), target.len(), base, wave))
                        .map(WaveRange::width)
                        .max()
                        .unwrap();
                    // The smallest resident chunk, one wider by a partial block, and chunks of
                    // about half the band.
                    for max_cells in [widest, widest + 9, total / 2] {
                        assert!(widest <= max_cells && max_cells < total);
                        let config = AlignmentConfig { max_cells, ..base };
                        assert!(assert_wave_kernels_match(&query, &target, config).0.is_ok());
                    }
                }
            }
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn known_band_limited_case_matches_scalar_kernels() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }
        // Inputs of default_band_retains_known_gap256_bounded_counterexample.
        let dna = |mut state: u64, length: usize| {
            (0..length)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    b"ACGT"[(state & 3) as usize]
                })
                .collect::<Vec<_>>()
        };
        let query = dna(0x9e37_79b9_7f4a_7c15 ^ 1_000, 1_000);
        let mut homolog = query.clone();
        homolog.splice(500..500, dna(0xabcd_dcba ^ 256, 256));
        let mut target = dna(0x1234_5678_9abc_def0 ^ 4_096, 4_096);
        target[1_420..2_676].copy_from_slice(&homolog);
        // The default band stops at the 256-base insertion; a wider band in that test finds the
        // longer alignment with the dispatched kernel.
        let config = AlignmentConfig {
            diagonal_offset: 1_420,
            band_width: 128,
            max_cells: 20_000_000,
            ..AlignmentConfig::default()
        };
        let (alignment, _) = assert_wave_kernels_match(&query, &target, config);
        let alignment = alignment.unwrap();
        assert_eq!((alignment.score, alignment.cigar.as_str()), (1_000, "500="));
    }

    #[test]
    fn trace_alignment_permits_exclude_wait_and_release() {
        let pool = AlignmentBytePool::new(10);
        let first = pool.acquire(8).unwrap();
        std::thread::scope(|scope| {
            let (started_tx, started_rx) = std::sync::mpsc::channel();
            let (acquired_tx, acquired_rx) = std::sync::mpsc::channel();
            let pool_ref = &pool;
            let waiter = scope.spawn(move || {
                started_tx.send(()).unwrap();
                let second = pool_ref.acquire(3).unwrap();
                acquired_tx.send(pool_ref.available()).unwrap();
                drop(second);
            });
            started_rx.recv().unwrap();
            assert_eq!(pool.available(), 2);
            assert!(acquired_rx.try_recv().is_err());
            drop(first);
            assert_eq!(
                acquired_rx
                    .recv_timeout(std::time::Duration::from_secs(1))
                    .unwrap(),
                7
            );
            waiter.join().unwrap();
        });
        assert_eq!(pool.available(), 10);
        assert!(matches!(
            pool.acquire(11),
            Err(AlignmentAdmissionError::RequestExceedsBudget {
                requested: 11,
                budget: 10
            })
        ));
    }

    #[test]
    fn admitted_workspace_preserves_alignment_and_bounds_default_trace() {
        let mut plain = AlignmentWorkspace::default();
        let expected = plain.align(b"ACGTACGT", b"ACGTGACGT", config()).unwrap();
        let bytes = trace_alignment_bytes(8, 9, 17, 9 * 17, 0, config()).unwrap();
        let mut admitted = TraceAlignmentWorkspace::acquire(bytes).unwrap();
        let actual = admitted
            .workspace_mut()
            .align(b"ACGTACGT", b"ACGTGACGT", config())
            .unwrap();
        assert_eq!(actual, expected);
        let trace_config = AlignmentConfig {
            max_cells: 1 << 24,
            ..AlignmentConfig::default()
        };
        let bytes = trace_alignment_bytes(
            64 * 1024,
            64 * 1024,
            128 * 1024,
            (64 * 1024 + 1) * 257,
            256,
            trace_config,
        )
        .unwrap();
        #[cfg(target_pointer_width = "64")]
        assert_eq!(bytes, 550_852_370);
    }

    #[test]
    fn affine_alignment_validates_its_cigar() {
        let mut workspace = AlignmentWorkspace::default();
        let alignment = workspace
            .align(b"ACGTACGT", b"ACGTGACGT", config())
            .unwrap();
        assert_eq!(alignment.cigar, "4=1I4=");
        assert_eq!(alignment.insertions, 1);
        assert_eq!(alignment.identity(), 8.0 / 9.0);
        let json = serde_json::to_value(&alignment).unwrap();
        assert_eq!(json["identity"].as_f64().unwrap(), alignment.identity());
        assert_eq!(
            serde_json::from_value::<Alignment>(json).unwrap(),
            alignment
        );
        alignment.validate_cigar().unwrap();
    }

    #[test]
    fn reverse_interval_uses_forward_coordinates() {
        let mut workspace = AlignmentWorkspace::default();
        let alignment = workspace
            .align_oriented(b"AGGACTT", b"AAGTCCT", 100, Strand::Reverse, config())
            .unwrap();
        assert_eq!(alignment.target_interval, Interval::new(100, 107).unwrap());
        assert_eq!(alignment.strand, Strand::Reverse);
    }

    #[test]
    fn complements_every_accepted_iupac_symbol() {
        for (&base, &expected) in b"ACGTURYSWKMBDHVN".iter().zip(b"TGCAAYRSWMKVHDBN") {
            assert_eq!(complement(base), expected);
            assert_eq!(complement(base.to_ascii_lowercase()), expected);
        }
    }

    #[test]
    fn workspace_reuses_matrix_capacity() {
        let mut workspace = AlignmentWorkspace::default();
        workspace.align(b"ACGT", b"ACGT", config()).unwrap();
        let capacity = workspace.capacity_cells();
        workspace.align(b"ACG", b"ACG", config()).unwrap();
        assert_eq!(workspace.capacity_cells(), capacity);
    }

    #[test]
    fn malformed_cigar_is_rejected() {
        assert_eq!(
            parse_cigar("4=1I4").unwrap_err(),
            AlignmentError::InvalidCigar
        );
        assert_eq!(
            parse_cigar("4=1I1I").unwrap_err(),
            AlignmentError::InvalidCigar
        );
    }

    fn core(strand: Strand, target_interval: Interval) -> Alignment {
        Alignment {
            score: 8,
            strand,
            query_interval: Interval::new(2, 6).unwrap(),
            target_interval,
            matches: 4,
            substitutions: 0,
            insertions: 0,
            deletions: 0,
            cigar: "4=".into(),
            edit_script: vec![EditRun {
                operation: EditOperation::Equal,
                length: 4,
            }],
        }
    }

    #[test]
    fn completes_both_endpoints_within_bound() {
        let mut workspace = AlignmentWorkspace::default();
        let completed = workspace
            .complete_endpoints(
                core(Strand::Forward, Interval::new(102, 106).unwrap()),
                b"TTACGTAA",
                b"TTACGTAA",
                100,
                2,
                config(),
            )
            .unwrap();
        assert_eq!(
            completed.alignment.query_interval,
            Interval::new(0, 8).unwrap()
        );
        assert_eq!(
            completed.alignment.target_interval,
            Interval::new(100, 108).unwrap()
        );
        assert_eq!(completed.alignment.cigar, "8=");
        assert_eq!(completed.metrics.left_query_bases, 2);
        assert_eq!(completed.metrics.right_query_bases, 2);
    }

    #[test]
    fn completion_preserves_reverse_coordinates() {
        let mut workspace = AlignmentWorkspace::default();
        let completed = workspace
            .complete_endpoints(
                core(Strand::Reverse, Interval::new(102, 106).unwrap()),
                b"TTACGCAA",
                b"TTGCGTAA",
                100,
                2,
                config(),
            )
            .unwrap();
        assert_eq!(
            completed.alignment.target_interval,
            Interval::new(100, 108).unwrap()
        );
        assert_eq!(completed.alignment.cigar, "8=");
    }

    fn assert_endpoint_case(
        matrix: &mut EndpointMatrix,
        reference: &mut Vec<EndpointCell>,
        query: &[u8],
        target: &[u8],
        config: AlignmentConfig,
    ) {
        let actual = anchored_semiglobal(matrix, query, target, config, None);
        let expected = anchored_semiglobal_reference(reference, query, target, config);
        match (actual, expected) {
            (Ok(actual), Ok(expected)) => {
                assert_eq!(
                    actual.runs, expected.runs,
                    "query={query:?} target={target:?}"
                );
                assert_eq!(actual.query_bases, expected.query_bases);
                assert_eq!(actual.target_bases, expected.target_bases);
                assert_eq!(actual.matrix_cells, expected.matrix_cells);
                assert!(matrix.previous.iter().all(|packed| packed & 0xc0 == 0));
            }
            (Err(actual), Err(expected)) => assert_eq!(actual, expected),
            (Ok(_), Err(expected)) => panic!(
                "candidate succeeded for query={query:?} target={target:?}, reference failed: {expected}"
            ),
            (Err(actual), Ok(_)) => panic!(
                "candidate failed for query={query:?} target={target:?}, reference succeeded: {actual}"
            ),
        }
    }

    fn ac_strings(max_len: usize) -> Vec<Vec<u8>> {
        let mut strings = Vec::new();
        for len in 0..=max_len {
            for bits in 0..(1usize << len) {
                strings.push(
                    (0..len)
                        .map(|index| if bits & (1 << index) == 0 { b'A' } else { b'C' })
                        .collect(),
                );
            }
        }
        strings
    }

    #[test]
    fn endpoint_kernel_matches_reference_exhaustively() {
        let strings = ac_strings(4);
        let configs = [
            config(),
            AlignmentConfig {
                match_score: 1,
                mismatch_score: 0,
                gap_open_score: 0,
                gap_extend_score: 0,
                ..config()
            },
            AlignmentConfig {
                match_score: 7,
                mismatch_score: -11,
                gap_open_score: -13,
                gap_extend_score: -17,
                ..config()
            },
            AlignmentConfig {
                match_score: i32::MAX,
                mismatch_score: i32::MIN,
                gap_open_score: i32::MIN,
                gap_extend_score: i32::MIN,
                ..config()
            },
        ];
        let mut matrix = EndpointMatrix::default();
        let mut reference = Vec::new();
        for config in configs {
            for query in &strings {
                for target in &strings {
                    assert_endpoint_case(&mut matrix, &mut reference, query, target, config);
                }
            }
        }
    }

    #[test]
    fn endpoint_kernel_matches_reference_for_bounded_mixed_inputs() {
        const BASES: &[u8] = b"ACGTNRYKMacgtnryk";
        let mut state = 0x6a09_e667_f3bc_c909u64;
        let mut matrix = EndpointMatrix::default();
        let mut reference = Vec::new();
        for case in 0..256usize {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let query_len = if case % 17 == 0 {
                257
            } else {
                (state as usize) % 96
            };
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let target_len = if case % 19 == 0 {
                257
            } else {
                (state as usize) % 96
            };
            let mut make_sequence = |len| {
                (0..len)
                    .map(|_| {
                        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                        BASES[(state >> 59) as usize % BASES.len()]
                    })
                    .collect::<Vec<_>>()
            };
            let query = make_sequence(query_len);
            let target = make_sequence(target_len);
            let config = AlignmentConfig {
                match_score: 3,
                mismatch_score: -7,
                gap_open_score: -9,
                gap_extend_score: -2,
                max_cells: 70_000,
                ..config()
            };
            assert_endpoint_case(&mut matrix, &mut reference, &query, &target, config);
        }
        for (query_len, target_len) in [(1, 257), (257, 1), (256, 257), (257, 256), (257, 257)] {
            let query = vec![b'A'; query_len];
            let target = vec![b'C'; target_len];
            assert_endpoint_case(
                &mut matrix,
                &mut reference,
                &query,
                &target,
                AlignmentConfig {
                    max_cells: 70_000,
                    ..config()
                },
            );
        }
    }

    #[test]
    fn endpoint_kernel_matches_reference_for_wide_saturating_and_tied_inputs() {
        const BASES: &[u8] = b"ACgtN";
        let configs = [
            AlignmentConfig {
                match_score: i32::MAX,
                mismatch_score: i32::MIN,
                gap_open_score: i32::MIN,
                gap_extend_score: i32::MIN,
                ..config()
            },
            AlignmentConfig {
                match_score: 1,
                mismatch_score: 0,
                gap_open_score: 0,
                gap_extend_score: 0,
                ..config()
            },
        ];
        let mut state = 0xbb67_ae85_84ca_a73bu64;
        let mut matrix = EndpointMatrix::default();
        let mut reference = Vec::new();
        for case in 0..96usize {
            let mut next = |modulus: usize| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                (state >> 33) as usize % modulus
            };
            let (query_len, target_len) = (8 + next(40), 8 + next(40));
            let query = (0..query_len)
                .map(|_| BASES[next(BASES.len())])
                .collect::<Vec<_>>();
            let target = (0..target_len)
                .map(|_| BASES[next(BASES.len())])
                .collect::<Vec<_>>();
            let config = configs[case % configs.len()];
            assert_endpoint_case(&mut matrix, &mut reference, &query, &target, config);
        }
    }

    #[test]
    fn endpoint_kernel_reuses_buffers_and_preserves_errors() {
        let mut matrix = EndpointMatrix::default();
        let mut reference = Vec::new();
        for length in [257, 1, 96, 3] {
            assert_endpoint_case(
                &mut matrix,
                &mut reference,
                &vec![b'A'; length],
                &vec![b'A'; length],
                AlignmentConfig {
                    max_cells: 70_000,
                    ..config()
                },
            );
        }
        let limited = AlignmentConfig {
            max_cells: 89,
            ..config()
        };
        assert_endpoint_case(
            &mut matrix,
            &mut reference,
            b"AAAAAAAA",
            b"AAAAAAAAA",
            limited,
        );
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn endpoint_waves_match_scalar_kernel_bits_at_lane_boundaries() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }
        let lengths = [1, 7, 8, 9, 15, 16, 17, 24, 25];
        let mut state = 0x510e_527f_ade6_82d1u64;
        let mut next = move |bound: usize| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 33) as usize % bound
        };
        let mut vector = EndpointMatrix::default();
        let mut scalar = EndpointMatrix::default();
        for (case, config) in tie_and_saturation_configs().into_iter().enumerate() {
            for query_len in lengths {
                for target_len in lengths {
                    let mut sequence = |length| {
                        (0..length)
                            .map(|_| if case == 1 { b'A' } else { b"ACgtN"[next(5)] })
                            .collect::<Vec<_>>()
                    };
                    let (query, target) = (sequence(query_len), sequence(target_len));
                    let actual = anchored_semiglobal(&mut vector, &query, &target, config, None);
                    SCALAR_ENDPOINT_KERNEL.set(true);
                    let expected = anchored_semiglobal(&mut scalar, &query, &target, config, None);
                    SCALAR_ENDPOINT_KERNEL.set(false);
                    let summary = |result: Result<EndpointResult, AlignmentError>| {
                        result.map(|result| {
                            (
                                result.runs,
                                result.query_bases,
                                result.target_bases,
                                result.matrix_cells,
                            )
                        })
                    };
                    let inputs = format!("{config:?} query={query:?} target={target:?}");
                    assert_eq!(summary(actual), summary(expected), "{inputs}");
                    assert_eq!(vector.previous, scalar.previous, "{inputs}");
                }
            }
        }
        assert!(!vector.waves.is_empty());
        assert!(scalar.waves.is_empty());
    }

    #[derive(serde::Deserialize)]
    struct RetainedEndpointFixture {
        query: Vec<u8>,
        target: Vec<u8>,
        target_start: u64,
        scoring: [i32; 4],
        band_width: u32,
        diagonal_offset: i64,
        max_cells: usize,
        endpoint_bases: usize,
        circular: bool,
        min_identity: f64,
        min_aligned_bases: u64,
        core: Alignment,
        completed: Alignment,
        metrics: [usize; 5],
    }

    #[test]
    #[ignore = "requires retained BCF fixtures in JAM_ALIGNMENT_FIXTURES"]
    fn retained_bcf_endpoint_fixtures_match() {
        let root = std::env::var_os("JAM_ALIGNMENT_FIXTURES").expect("fixture directory");
        let mut paths = std::fs::read_dir(root)
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .filter(|path| {
                path.extension()
                    .is_some_and(|extension| extension == "json")
                    && path
                        .file_name()
                        .unwrap()
                        .to_string_lossy()
                        .starts_with("task-")
            })
            .collect::<Vec<_>>();
        paths.sort();
        assert!(!paths.is_empty());
        assert!(paths.len() <= 16, "fixture set must remain bounded");
        for path in paths {
            let fixture: RetainedEndpointFixture =
                serde_json::from_reader(std::fs::File::open(&path).unwrap()).unwrap();
            let config = AlignmentConfig {
                match_score: fixture.scoring[0],
                mismatch_score: fixture.scoring[1],
                gap_open_score: fixture.scoring[2],
                gap_extend_score: fixture.scoring[3],
                band_width: fixture.band_width,
                diagonal_offset: fixture.diagonal_offset,
                max_cells: fixture.max_cells,
            };
            let mut workspace = AlignmentWorkspace::default();
            let core = workspace
                .align_oriented(
                    &fixture.query,
                    &fixture.target,
                    fixture.target_start,
                    fixture.core.strand,
                    config,
                )
                .unwrap();
            assert_eq!(core, fixture.core, "{}", path.display());
            let completion = workspace
                .complete_endpoints(
                    core.clone(),
                    &fixture.query,
                    &fixture.target,
                    fixture.target_start,
                    fixture.endpoint_bases,
                    config,
                )
                .unwrap();
            assert_eq!(
                completion.alignment,
                fixture.completed,
                "{}",
                path.display()
            );
            assert_eq!(
                [
                    completion.metrics.left_query_bases,
                    completion.metrics.left_target_bases,
                    completion.metrics.right_query_bases,
                    completion.metrics.right_target_bases,
                    completion.metrics.matrix_cells
                ],
                fixture.metrics,
                "{}",
                path.display()
            );
            let selected = if completion.alignment.identity() >= fixture.min_identity {
                &completion.alignment
            } else {
                &core
            };
            let expected = if fixture.completed.identity() >= fixture.min_identity {
                &fixture.completed
            } else {
                &fixture.core
            };
            assert_eq!(selected, expected, "{}", path.display());
            assert_eq!(
                selected.identity() >= fixture.min_identity
                    && selected.query_interval.len() >= fixture.min_aligned_bases,
                expected.identity() >= fixture.min_identity
                    && expected.query_interval.len() >= fixture.min_aligned_bases
            );
            let _ = fixture.circular;
        }
    }

    #[test]
    fn semiglobal_completion_leaves_outer_overhang_free() {
        let mut cells = EndpointMatrix::default();
        let result =
            anchored_semiglobal(&mut cells, b"ACGTACGT", b"ACGTACGTCCCC", config(), None).unwrap();
        assert_eq!((result.query_bases, result.target_bases), (8, 8));
        assert_eq!(cigar_from_runs(&result.runs).unwrap(), "8=");
    }

    #[test]
    fn semiglobal_workspace_growth_matches_fresh_workspace() {
        let mut reused = EndpointMatrix::default();
        anchored_semiglobal(&mut reused, b"A", b"A", config(), None).unwrap();
        let reused_result = anchored_semiglobal(&mut reused, b"C", b"AAA", config(), None).unwrap();
        let fresh_result =
            anchored_semiglobal(&mut EndpointMatrix::default(), b"C", b"AAA", config(), None)
                .unwrap();

        assert_eq!(reused_result.runs, fresh_result.runs);
        assert_eq!(reused_result.query_bases, fresh_result.query_bases);
        assert_eq!(reused_result.target_bases, fresh_result.target_bases);
    }
}
