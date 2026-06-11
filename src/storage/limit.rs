//! A limit on how many requests are in flight at once.

use std::sync::{Condvar, Mutex, PoisonError};

/// A counting limit shared by many threads. A thread that finds the limit reached waits until
/// another one is done, so the number of requests to a remote store stays bounded.
pub struct InFlight {
    limit: usize,
    active: Mutex<usize>,
    freed: Condvar,
}

/// One request in flight; dropping it frees its place.
pub struct Permit<'a>(&'a InFlight);

impl InFlight {
    /// A limit of `limit` requests at once (at least one).
    pub fn new(limit: usize) -> Self {
        Self {
            limit: limit.max(1),
            active: Mutex::new(0),
            freed: Condvar::new(),
        }
    }

    /// Waits for a free place and takes it.
    pub fn acquire(&self) -> Permit<'_> {
        let mut active = self.active.lock().unwrap_or_else(PoisonError::into_inner);
        while *active >= self.limit {
            active = self
                .freed
                .wait(active)
                .unwrap_or_else(PoisonError::into_inner);
        }
        *active += 1;
        Permit(self)
    }
}

impl Drop for Permit<'_> {
    fn drop(&mut self) {
        *self.0.active.lock().unwrap_or_else(PoisonError::into_inner) -= 1;
        self.0.freed.notify_one();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    #[test]
    fn never_more_than_the_limit_run_at_once() {
        let limit = InFlight::new(3);
        let (active, peak) = (AtomicUsize::new(0), AtomicUsize::new(0));
        std::thread::scope(|scope| {
            for _ in 0..12 {
                scope.spawn(|| {
                    let _permit = limit.acquire();
                    let now = active.fetch_add(1, Ordering::SeqCst) + 1;
                    peak.fetch_max(now, Ordering::SeqCst);
                    std::thread::sleep(std::time::Duration::from_millis(10));
                    active.fetch_sub(1, Ordering::SeqCst);
                });
            }
        });
        assert_eq!(peak.load(Ordering::SeqCst), 3);
    }
}
