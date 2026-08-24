//! FIFO bounded producer/consumer queue for candidate work.

use std::collections::VecDeque;
use std::sync::{Arc, Condvar, Mutex, MutexGuard, PoisonError};

/// Errors produced by queue operations.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum QueueError {
    /// A producer attempted to enqueue after closure.
    Closed,
    /// A nonblocking producer found the queue full.
    Full,
    /// A queue with no slots cannot provide bounded progress.
    InvalidCapacity,
}

impl std::fmt::Display for QueueError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Closed => formatter.write_str("candidate work queue is closed"),
            Self::Full => formatter.write_str("candidate work queue is full"),
            Self::InvalidCapacity => formatter.write_str("candidate work queue capacity is zero"),
        }
    }
}

impl std::error::Error for QueueError {}

/// Snapshot of queue state for diagnostics and tests.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct QueueState {
    pub capacity: usize,
    pub len: usize,
    pub closed: bool,
}

#[derive(Debug)]
struct Inner<T> {
    state: Mutex<State<T>>,
    changed: Condvar,
}

#[derive(Debug)]
struct State<T> {
    capacity: usize,
    closed: bool,
    items: VecDeque<T>,
}

/// A FIFO queue with a hard item bound and blocking push/pop operations.
///
/// Both producer and consumer waiters are woken together and re-check the
/// queue under the mutex.  Items remain FIFO, close drains already-enqueued
/// work, and all waiters are released on close.  The queue therefore cannot
/// strand a producer behind a consumer that stopped waiting.
#[derive(Clone, Debug)]
pub struct CandidateWorkQueue<T> {
    inner: Arc<Inner<T>>,
}

/// More descriptive alias used by trace callers.
pub type BoundedCandidateQueue<T> = CandidateWorkQueue<T>;

/// A small scheduler façade for candidate producers and workers.
#[derive(Clone, Debug)]
pub struct CandidateScheduler<T> {
    queue: CandidateWorkQueue<T>,
}

impl<T> CandidateWorkQueue<T> {
    /// Creates a queue with at least one bounded slot.
    pub fn new(capacity: usize) -> Result<Self, QueueError> {
        if capacity == 0 {
            return Err(QueueError::InvalidCapacity);
        }
        Ok(Self {
            inner: Arc::new(Inner {
                state: Mutex::new(State {
                    capacity,
                    closed: false,
                    items: VecDeque::with_capacity(capacity),
                }),
                changed: Condvar::new(),
            }),
        })
    }

    /// Enqueues an item, waiting for space if necessary.
    pub fn push(&self, item: T) -> Result<(), QueueError> {
        let mut state = lock_state(&self.inner);
        loop {
            if state.closed {
                return Err(QueueError::Closed);
            }
            if state.items.len() < state.capacity {
                state.items.push_back(item);
                drop(state);
                self.inner.changed.notify_all();
                return Ok(());
            }
            state = self
                .inner
                .changed
                .wait(state)
                .unwrap_or_else(PoisonError::into_inner);
        }
    }

    /// Enqueues an item without waiting for space.
    pub fn try_push(&self, item: T) -> Result<(), QueueError> {
        let mut state = lock_state(&self.inner);
        if state.closed {
            return Err(QueueError::Closed);
        }
        if state.items.len() == state.capacity {
            return Err(QueueError::Full);
        }
        state.items.push_back(item);
        drop(state);
        self.inner.changed.notify_all();
        Ok(())
    }

    /// Removes an item, waiting until work arrives or the queue closes.
    ///
    /// `Ok(None)` means the queue is closed and drained.
    pub fn pop(&self) -> Result<Option<T>, QueueError> {
        let mut state = lock_state(&self.inner);
        loop {
            if let Some(item) = state.items.pop_front() {
                drop(state);
                self.inner.changed.notify_all();
                return Ok(Some(item));
            }
            if state.closed {
                return Ok(None);
            }
            state = self
                .inner
                .changed
                .wait(state)
                .unwrap_or_else(PoisonError::into_inner);
        }
    }

    /// Removes one item without waiting.
    pub fn try_pop(&self) -> Option<T> {
        let mut state = lock_state(&self.inner);
        let item = state.items.pop_front();
        if item.is_some() {
            drop(state);
            self.inner.changed.notify_all();
        }
        item
    }

    /// Closes producers while allowing workers to drain existing items.
    pub fn close(&self) {
        let mut state = lock_state(&self.inner);
        state.closed = true;
        drop(state);
        self.inner.changed.notify_all();
    }

    #[must_use]
    pub fn is_closed(&self) -> bool {
        lock_state(&self.inner).closed
    }

    #[must_use]
    pub fn len(&self) -> usize {
        lock_state(&self.inner).items.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    #[must_use]
    pub fn capacity(&self) -> usize {
        lock_state(&self.inner).capacity
    }

    #[must_use]
    pub fn state(&self) -> QueueState {
        let state = lock_state(&self.inner);
        QueueState {
            capacity: state.capacity,
            len: state.items.len(),
            closed: state.closed,
        }
    }
}

impl<T> CandidateScheduler<T> {
    pub fn new(capacity: usize) -> Result<Self, QueueError> {
        Ok(Self {
            queue: CandidateWorkQueue::new(capacity)?,
        })
    }

    #[must_use]
    pub fn queue(&self) -> &CandidateWorkQueue<T> {
        &self.queue
    }

    #[must_use]
    pub fn into_queue(self) -> CandidateWorkQueue<T> {
        self.queue
    }

    pub fn submit(&self, candidate: T) -> Result<(), QueueError> {
        self.queue.push(candidate)
    }

    pub fn next(&self) -> Result<Option<T>, QueueError> {
        self.queue.pop()
    }

    pub fn finish(&self) {
        self.queue.close();
    }
}

fn lock_state<T>(inner: &Inner<T>) -> MutexGuard<'_, State<T>> {
    inner.state.lock().unwrap_or_else(PoisonError::into_inner)
}
