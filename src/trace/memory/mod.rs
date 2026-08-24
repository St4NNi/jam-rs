//! Bounded memory admission and reusable per-worker trace scratch space.
//!
//! A trace candidate can touch several independently sized resources before it
//! is serialized.  [`TraceMemorySemaphore`] admits the complete estimate as
//! one FIFO reservation, so a worker cannot start a candidate and discover
//! that only part of its working set fits.  The reservation releases its
//! bytes on drop, including when candidate processing returns an error or
//! panics.

mod early_cap;
mod scratch;
mod semaphore;

pub use early_cap::{FixedTopK, HeapRank, PerContigTopK, retain_per_contig_top_k, retain_top_k};
pub use scratch::{
    ScratchCapacities, ScratchLease, ScratchPool, ScratchPoolError, TraceWorkerScratch,
    TraceWorkerScratchPool,
};
pub use semaphore::{
    BlockingMemorySemaphore, GlobalMemoryBudget, MemoryReservation, MemorySemaphore,
    TraceMemoryReservation, TraceMemorySemaphore,
};
