//! FIFO, blocking admission for candidate working sets.

use crate::trace::contracts::{CandidateMemoryBudget, MemoryBudgetError, TraceMemoryEstimate};
use std::collections::VecDeque;
use std::sync::{Arc, Condvar, Mutex, MutexGuard, PoisonError};

#[derive(Debug)]
struct Waiter {
    ticket: u64,
    bytes: u64,
}

#[derive(Debug)]
struct State {
    capacity: u64,
    in_use: u64,
    next_ticket: u64,
    closed: bool,
    waiters: VecDeque<Waiter>,
}

#[derive(Debug)]
struct Inner {
    state: Mutex<State>,
    changed: Condvar,
}

/// A process-local global memory budget for candidate processing.
///
/// Admission is FIFO.  A later, smaller request cannot bypass an earlier
/// request, which keeps progress deterministic under contention and prevents
/// a stream of small candidates from starving a large candidate.  Capacity is
/// measured in bytes and is never exceeded by admitted reservations.
#[derive(Clone, Debug)]
pub struct TraceMemorySemaphore {
    inner: Arc<Inner>,
}

/// Compatibility-oriented name for [`TraceMemorySemaphore`].
pub type MemorySemaphore = TraceMemorySemaphore;

/// Compatibility-oriented name for [`TraceMemorySemaphore`].
pub type GlobalMemoryBudget = TraceMemorySemaphore;

/// Explicit name for callers that want to document blocking admission.
pub type BlockingMemorySemaphore = TraceMemorySemaphore;

/// A released-on-drop reservation held by one candidate.
#[derive(Debug)]
pub struct TraceMemoryReservation {
    inner: Arc<Inner>,
    bytes: u64,
}

/// Compatibility-oriented name for [`TraceMemoryReservation`].
pub type MemoryReservation = TraceMemoryReservation;

impl TraceMemorySemaphore {
    /// Creates a semaphore with the supplied byte capacity.
    #[must_use]
    pub fn new(capacity_bytes: u64) -> Self {
        Self {
            inner: Arc::new(Inner {
                state: Mutex::new(State {
                    capacity: capacity_bytes,
                    in_use: 0,
                    next_ticket: 0,
                    closed: false,
                    waiters: VecDeque::new(),
                }),
                changed: Condvar::new(),
            }),
        }
    }

    /// Creates a semaphore from a platform-sized byte count.
    #[must_use]
    pub fn from_usize(capacity_bytes: usize) -> Self {
        Self::new(capacity_bytes as u64)
    }

    /// The configured byte capacity.
    #[must_use]
    pub fn capacity_bytes(&self) -> u64 {
        lock_state(&self.inner).capacity
    }

    /// Bytes held by currently admitted candidates.
    #[must_use]
    pub fn in_use_bytes(&self) -> u64 {
        lock_state(&self.inner).in_use
    }

    /// Bytes not currently held by admitted candidates.
    #[must_use]
    pub fn available_bytes(&self) -> u64 {
        let state = lock_state(&self.inner);
        state.capacity.saturating_sub(state.in_use)
    }

    /// Number of candidates waiting for admission.
    #[must_use]
    pub fn waiting_candidates(&self) -> usize {
        lock_state(&self.inner).waiters.len()
    }

    /// Returns whether no future candidate can be admitted.
    #[must_use]
    pub fn is_closed(&self) -> bool {
        lock_state(&self.inner).closed
    }

    /// Closes admission and wakes every blocked waiter.
    ///
    /// Existing reservations remain valid and release normally.  Calling
    /// `close` more than once is harmless.
    pub fn close(&self) {
        let mut state = lock_state(&self.inner);
        state.closed = true;
        self.inner.changed.notify_all();
    }

    /// Admits a candidate without waiting when enough capacity is available.
    ///
    /// `Ok(None)` means that the request would fit in the budget but another
    /// reservation or an earlier waiter currently owns the available bytes.
    pub fn try_acquire(
        &self,
        estimate: TraceMemoryEstimate,
    ) -> Result<Option<TraceMemoryReservation>, MemoryBudgetError> {
        let bytes = estimate.total_bytes();
        let mut state = lock_state(&self.inner);
        validate_request(&state, bytes)?;
        if !state.waiters.is_empty() || state.capacity.saturating_sub(state.in_use) < bytes {
            return Ok(None);
        }
        state.in_use = state
            .in_use
            .checked_add(bytes)
            .expect("validated capacity prevents memory-use overflow");
        drop(state);
        Ok(Some(TraceMemoryReservation {
            inner: Arc::clone(&self.inner),
            bytes,
        }))
    }

    /// Acquires a complete candidate reservation, blocking behind earlier
    /// requests when the currently available bytes are insufficient.
    pub fn acquire(
        &self,
        estimate: TraceMemoryEstimate,
    ) -> Result<TraceMemoryReservation, MemoryBudgetError> {
        self.acquire_bytes(estimate.total_bytes())
    }

    /// Acquires a candidate reservation and runs `work` while it is held.
    ///
    /// The reservation is dropped before this function returns, which makes
    /// the release-after-serialization boundary explicit to callers.
    pub fn with_candidate<T, E, F>(
        &self,
        estimate: TraceMemoryEstimate,
        work: F,
    ) -> Result<T, MemoryBudgetError>
    where
        F: FnOnce(&TraceMemoryReservation) -> Result<T, E>,
        E: Into<MemoryBudgetError>,
    {
        let reservation = self.acquire(estimate)?;
        let result = work(&reservation);
        drop(reservation);
        result.map_err(Into::into)
    }

    /// Runs a candidate with an infallible closure and always releases its
    /// reservation before returning.
    pub fn with_candidate_result<T, F>(
        &self,
        estimate: TraceMemoryEstimate,
        work: F,
    ) -> Result<T, MemoryBudgetError>
    where
        F: FnOnce(&TraceMemoryReservation) -> T,
    {
        let reservation = self.acquire(estimate)?;
        let result = work(&reservation);
        drop(reservation);
        Ok(result)
    }

    fn acquire_bytes(&self, bytes: u64) -> Result<TraceMemoryReservation, MemoryBudgetError> {
        let mut state = lock_state(&self.inner);
        validate_request(&state, bytes)?;

        let ticket = state.next_ticket;
        state.next_ticket = state.next_ticket.wrapping_add(1);
        state.waiters.push_back(Waiter { ticket, bytes });

        loop {
            if state.closed {
                remove_waiter(&mut state.waiters, ticket);
                drop(state);
                self.inner.changed.notify_all();
                return Err(MemoryBudgetError::Closed);
            }

            let at_front = state
                .waiters
                .front()
                .is_some_and(|waiter| waiter.ticket == ticket);
            if at_front && state.capacity.saturating_sub(state.in_use) >= bytes {
                let waiter = state
                    .waiters
                    .pop_front()
                    .expect("front waiter was checked above");
                debug_assert_eq!(waiter.ticket, ticket);
                debug_assert_eq!(waiter.bytes, bytes);
                state.in_use = state
                    .in_use
                    .checked_add(bytes)
                    .expect("validated capacity prevents memory-use overflow");
                drop(state);
                self.inner.changed.notify_all();
                return Ok(TraceMemoryReservation {
                    inner: Arc::clone(&self.inner),
                    bytes,
                });
            }

            state = wait_state(&self.inner.changed, state);
        }
    }
}

impl CandidateMemoryBudget for TraceMemorySemaphore {
    type Reservation = TraceMemoryReservation;

    fn acquire(
        &self,
        estimate: TraceMemoryEstimate,
    ) -> Result<Self::Reservation, MemoryBudgetError> {
        TraceMemorySemaphore::acquire(self, estimate)
    }
}

impl TraceMemoryReservation {
    /// Bytes held by this reservation.
    #[must_use]
    pub const fn bytes(&self) -> u64 {
        self.bytes
    }

    /// Returns whether this reservation still owns bytes.
    #[must_use]
    pub const fn is_active(&self) -> bool {
        self.bytes != 0
    }
}

impl Drop for TraceMemoryReservation {
    fn drop(&mut self) {
        if self.bytes == 0 {
            return;
        }
        let mut state = lock_state(&self.inner);
        state.in_use = state
            .in_use
            .checked_sub(self.bytes)
            .expect("reservation bytes must be accounted by the semaphore");
        self.bytes = 0;
        drop(state);
        self.inner.changed.notify_all();
    }
}

fn validate_request(state: &State, bytes: u64) -> Result<(), MemoryBudgetError> {
    if bytes > state.capacity {
        return Err(MemoryBudgetError::CandidateTooLarge {
            requested_bytes: bytes,
            budget_bytes: state.capacity,
        });
    }
    if state.closed {
        return Err(MemoryBudgetError::Closed);
    }
    Ok(())
}

fn remove_waiter(waiters: &mut VecDeque<Waiter>, ticket: u64) {
    if let Some(index) = waiters.iter().position(|waiter| waiter.ticket == ticket) {
        waiters.remove(index);
    }
}

fn lock_state(inner: &Inner) -> MutexGuard<'_, State> {
    inner.state.lock().unwrap_or_else(PoisonError::into_inner)
}

fn wait_state<'a>(changed: &Condvar, state: MutexGuard<'a, State>) -> MutexGuard<'a, State> {
    changed.wait(state).unwrap_or_else(PoisonError::into_inner)
}
