//! Reusable worker-local buffers for trace candidates.

use std::collections::VecDeque;
use std::ops::{Deref, DerefMut};
use std::sync::{Arc, Condvar, Mutex, MutexGuard, PoisonError};

/// Initial capacities for the buffers used by one trace worker.
///
/// Values are element counts for the corresponding typed buffer, not bytes.
/// The byte estimate is available through [`Self::estimated_bytes`] for
/// admission accounting.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct ScratchCapacities {
    pub key_decode: usize,
    pub posting_decode: usize,
    pub anchors: usize,
    pub chains: usize,
    pub sequence: usize,
    pub score_rows: usize,
    pub traceback: usize,
    pub cigar: usize,
    pub edit_script: usize,
    pub output: usize,
}

impl ScratchCapacities {
    #[must_use]
    pub const fn estimated_bytes(self) -> u64 {
        (self.key_decode as u64)
            .saturating_add(self.posting_decode as u64)
            .saturating_add(self.sequence as u64)
            .saturating_add(self.cigar as u64)
            .saturating_add(self.edit_script as u64)
            .saturating_add(self.output as u64)
            .saturating_add((self.anchors as u64).saturating_mul(std::mem::size_of::<u64>() as u64))
            .saturating_add((self.chains as u64).saturating_mul(std::mem::size_of::<u64>() as u64))
            .saturating_add(
                (self.score_rows as u64).saturating_mul(std::mem::size_of::<i32>() as u64),
            )
            .saturating_add(
                (self.traceback as u64).saturating_mul(std::mem::size_of::<i32>() as u64),
            )
    }
}

/// One worker's reusable typed buffers.
///
/// The byte-oriented buffers intentionally remain independent from the
/// concrete archive/model types.  Callers can decode into them directly or
/// use the typed `u64`/`i32` buffers for compact anchor, chain, and score
/// state.  [`TraceWorkerScratch::reset`] only clears lengths; capacities stay
/// allocated for the next candidate.
#[derive(Debug, Default)]
pub struct TraceWorkerScratch {
    key_decode: Vec<u8>,
    posting_decode: Vec<u8>,
    anchors: Vec<u64>,
    chains: Vec<u64>,
    sequence: Vec<u8>,
    score_rows: Vec<i32>,
    traceback: Vec<i32>,
    cigar: Vec<u8>,
    edit_script: Vec<u8>,
    output: Vec<u8>,
}

impl TraceWorkerScratch {
    #[must_use]
    pub fn with_capacities(capacities: ScratchCapacities) -> Self {
        Self {
            key_decode: Vec::with_capacity(capacities.key_decode),
            posting_decode: Vec::with_capacity(capacities.posting_decode),
            anchors: Vec::with_capacity(capacities.anchors),
            chains: Vec::with_capacity(capacities.chains),
            sequence: Vec::with_capacity(capacities.sequence),
            score_rows: Vec::with_capacity(capacities.score_rows),
            traceback: Vec::with_capacity(capacities.traceback),
            cigar: Vec::with_capacity(capacities.cigar),
            edit_script: Vec::with_capacity(capacities.edit_script),
            output: Vec::with_capacity(capacities.output),
        }
    }

    #[must_use]
    pub fn capacities(&self) -> ScratchCapacities {
        ScratchCapacities {
            key_decode: self.key_decode.capacity(),
            posting_decode: self.posting_decode.capacity(),
            anchors: self.anchors.capacity(),
            chains: self.chains.capacity(),
            sequence: self.sequence.capacity(),
            score_rows: self.score_rows.capacity(),
            traceback: self.traceback.capacity(),
            cigar: self.cigar.capacity(),
            edit_script: self.edit_script.capacity(),
            output: self.output.capacity(),
        }
    }

    /// Ensures each buffer has at least the requested capacity.
    pub fn reserve(&mut self, capacities: ScratchCapacities) {
        reserve_to(&mut self.key_decode, capacities.key_decode);
        reserve_to(&mut self.posting_decode, capacities.posting_decode);
        reserve_to(&mut self.anchors, capacities.anchors);
        reserve_to(&mut self.chains, capacities.chains);
        reserve_to(&mut self.sequence, capacities.sequence);
        reserve_to(&mut self.score_rows, capacities.score_rows);
        reserve_to(&mut self.traceback, capacities.traceback);
        reserve_to(&mut self.cigar, capacities.cigar);
        reserve_to(&mut self.edit_script, capacities.edit_script);
        reserve_to(&mut self.output, capacities.output);
    }

    /// Clears all logical contents while retaining allocations.
    pub fn reset(&mut self) {
        self.key_decode.clear();
        self.posting_decode.clear();
        self.anchors.clear();
        self.chains.clear();
        self.sequence.clear();
        self.score_rows.clear();
        self.traceback.clear();
        self.cigar.clear();
        self.edit_script.clear();
        self.output.clear();
    }

    #[must_use]
    pub fn key_decode(&mut self) -> &mut Vec<u8> {
        &mut self.key_decode
    }

    #[must_use]
    pub fn posting_decode(&mut self) -> &mut Vec<u8> {
        &mut self.posting_decode
    }

    #[must_use]
    pub fn anchors(&mut self) -> &mut Vec<u64> {
        &mut self.anchors
    }

    #[must_use]
    pub fn chains(&mut self) -> &mut Vec<u64> {
        &mut self.chains
    }

    #[must_use]
    pub fn sequence(&mut self) -> &mut Vec<u8> {
        &mut self.sequence
    }

    #[must_use]
    pub fn score_rows(&mut self) -> &mut Vec<i32> {
        &mut self.score_rows
    }

    #[must_use]
    pub fn traceback(&mut self) -> &mut Vec<i32> {
        &mut self.traceback
    }

    #[must_use]
    pub fn cigar(&mut self) -> &mut Vec<u8> {
        &mut self.cigar
    }

    #[must_use]
    pub fn edit_script(&mut self) -> &mut Vec<u8> {
        &mut self.edit_script
    }

    #[must_use]
    pub fn output(&mut self) -> &mut Vec<u8> {
        &mut self.output
    }
}

fn reserve_to<T>(buffer: &mut Vec<T>, requested: usize) {
    let additional = requested.saturating_sub(buffer.capacity());
    if additional != 0 {
        buffer.reserve_exact(additional);
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ScratchPoolError {
    InvalidWorkerCount,
    Closed,
}

impl std::fmt::Display for ScratchPoolError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidWorkerCount => formatter.write_str("worker count must be positive"),
            Self::Closed => formatter.write_str("scratch pool is closed"),
        }
    }
}

impl std::error::Error for ScratchPoolError {}

#[derive(Debug)]
struct PoolState {
    available: VecDeque<TraceWorkerScratch>,
    closed: bool,
}

#[derive(Debug)]
struct PoolInner {
    state: Mutex<PoolState>,
    changed: Condvar,
}

/// A fixed-size blocking pool of worker scratch buffers.
#[derive(Clone, Debug)]
pub struct TraceWorkerScratchPool {
    inner: Arc<PoolInner>,
}

/// Short alias for [`TraceWorkerScratchPool`].
pub type ScratchPool = TraceWorkerScratchPool;

/// A checked-out worker scratch buffer returned to the pool on drop.
pub struct ScratchLease {
    inner: Arc<PoolInner>,
    scratch: Option<TraceWorkerScratch>,
}

impl TraceWorkerScratchPool {
    /// Creates and preallocates one scratch buffer per worker.
    pub fn new(
        worker_count: usize,
        capacities: ScratchCapacities,
    ) -> Result<Self, ScratchPoolError> {
        if worker_count == 0 {
            return Err(ScratchPoolError::InvalidWorkerCount);
        }
        let mut available = VecDeque::with_capacity(worker_count);
        for _ in 0..worker_count {
            available.push_back(TraceWorkerScratch::with_capacities(capacities));
        }
        Ok(Self {
            inner: Arc::new(PoolInner {
                state: Mutex::new(PoolState {
                    available,
                    closed: false,
                }),
                changed: Condvar::new(),
            }),
        })
    }

    /// Checks out one scratch buffer in FIFO worker order.
    pub fn acquire(&self) -> Result<ScratchLease, ScratchPoolError> {
        let mut state = lock_pool(&self.inner);
        loop {
            if state.closed {
                return Err(ScratchPoolError::Closed);
            }
            if let Some(scratch) = state.available.pop_front() {
                drop(state);
                return Ok(ScratchLease {
                    inner: Arc::clone(&self.inner),
                    scratch: Some(scratch),
                });
            }
            state = self
                .inner
                .changed
                .wait(state)
                .unwrap_or_else(PoisonError::into_inner);
        }
    }

    /// Checks out a buffer without waiting.
    pub fn try_acquire(&self) -> Result<Option<ScratchLease>, ScratchPoolError> {
        let mut state = lock_pool(&self.inner);
        if state.closed {
            return Err(ScratchPoolError::Closed);
        }
        if let Some(scratch) = state.available.pop_front() {
            drop(state);
            return Ok(Some(ScratchLease {
                inner: Arc::clone(&self.inner),
                scratch: Some(scratch),
            }));
        }
        Ok(None)
    }

    /// Closes the pool after checked-out leases return.
    pub fn close(&self) {
        let mut state = lock_pool(&self.inner);
        state.closed = true;
        self.inner.changed.notify_all();
    }

    #[must_use]
    pub fn available(&self) -> usize {
        lock_pool(&self.inner).available.len()
    }

    #[must_use]
    pub fn is_closed(&self) -> bool {
        lock_pool(&self.inner).closed
    }
}

impl Deref for ScratchLease {
    type Target = TraceWorkerScratch;

    fn deref(&self) -> &Self::Target {
        self.scratch.as_ref().expect("scratch lease is active")
    }
}

impl DerefMut for ScratchLease {
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.scratch.as_mut().expect("scratch lease is active")
    }
}

impl ScratchLease {
    /// Returns the buffer to the caller and prevents it from being returned
    /// to the pool.  This is useful when transferring ownership at shutdown.
    pub fn into_inner(mut self) -> TraceWorkerScratch {
        self.scratch.take().expect("scratch lease is active")
    }
}

impl Drop for ScratchLease {
    fn drop(&mut self) {
        let Some(mut scratch) = self.scratch.take() else {
            return;
        };
        scratch.reset();
        let mut state = lock_pool(&self.inner);
        if !state.closed {
            state.available.push_back(scratch);
        }
        drop(state);
        self.inner.changed.notify_all();
    }
}

fn lock_pool(inner: &PoolInner) -> MutexGuard<'_, PoolState> {
    inner.state.lock().unwrap_or_else(PoisonError::into_inner)
}
