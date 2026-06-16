//! The limits of a server: connections at once, searches in flight, and searches that wait.

use crate::search::SearchControl;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Condvar, Mutex, PoisonError};

use super::connection::POLL;

/// How much work a server takes at once.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ServerLimits {
    /// Open connections; a connection over the limit is closed at once.
    pub connections: usize,
    /// Searches that run at once.
    pub searches: usize,
    /// Searches that wait for a place; a search over the limit is refused.
    pub queue: usize,
}

impl Default for ServerLimits {
    fn default() -> Self {
        Self {
            connections: 64,
            searches: 4,
            queue: 32,
        }
    }
}

/// Why a search did not get a place.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum NoPlace {
    /// The queue is full.
    Busy,
    /// The search passed its deadline or was asked to stop while it waited.
    Stopped,
}

/// The counters behind `ServerLimits`, shared by all connections of a server.
pub struct Admission {
    limits: ServerLimits,
    connections: AtomicUsize,
    /// Searches running and searches waiting.
    searches: Mutex<(usize, usize)>,
    freed: Condvar,
}

/// One open connection; dropping it frees its place.
pub struct ConnectionPlace<'a>(&'a Admission);

/// One running search; dropping it frees its place.
pub struct SearchPlace<'a>(&'a Admission);

impl Admission {
    /// Counters for `limits`.
    pub fn new(limits: ServerLimits) -> Self {
        Self {
            limits,
            connections: AtomicUsize::new(0),
            searches: Mutex::new((0, 0)),
            freed: Condvar::new(),
        }
    }

    /// The limits.
    pub fn limits(&self) -> ServerLimits {
        self.limits
    }

    /// A place for one more connection, if the limit allows it.
    pub fn connection(&self) -> Option<ConnectionPlace<'_>> {
        let taken = self.connections.fetch_add(1, Ordering::SeqCst);
        if taken >= self.limits.connections {
            self.connections.fetch_sub(1, Ordering::SeqCst);
            return None;
        }
        Some(ConnectionPlace(self))
    }

    /// A place to run one search. It waits in the queue while all places are taken, and gives up
    /// when the queue is full or the search stops.
    pub fn search(&self, control: SearchControl<'_>) -> Result<SearchPlace<'_>, NoPlace> {
        let mut searches = self.searches.lock().unwrap_or_else(PoisonError::into_inner);
        if searches.0 < self.limits.searches.max(1) {
            searches.0 += 1;
            return Ok(SearchPlace(self));
        }
        if searches.1 >= self.limits.queue {
            return Err(NoPlace::Busy);
        }
        searches.1 += 1;
        loop {
            if control.check().is_err() {
                searches.1 -= 1;
                return Err(NoPlace::Stopped);
            }
            if searches.0 < self.limits.searches.max(1) {
                searches.1 -= 1;
                searches.0 += 1;
                return Ok(SearchPlace(self));
            }
            searches = self
                .freed
                .wait_timeout(searches, POLL)
                .unwrap_or_else(PoisonError::into_inner)
                .0;
        }
    }

    /// Searches running and waiting now.
    pub fn searches(&self) -> (usize, usize) {
        *self.searches.lock().unwrap_or_else(PoisonError::into_inner)
    }
}

impl Drop for ConnectionPlace<'_> {
    fn drop(&mut self) {
        self.0.connections.fetch_sub(1, Ordering::SeqCst);
    }
}

impl Drop for SearchPlace<'_> {
    fn drop(&mut self) {
        let mut searches = self
            .0
            .searches
            .lock()
            .unwrap_or_else(PoisonError::into_inner);
        searches.0 -= 1;
        self.0.freed.notify_one();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::AtomicBool;
    use std::time::{Duration, Instant};

    fn control(stop: &AtomicBool) -> SearchControl<'_> {
        SearchControl {
            deadline: None,
            stop,
        }
    }

    #[test]
    fn connections_over_the_limit_get_no_place() {
        let admission = Admission::new(ServerLimits {
            connections: 2,
            ..ServerLimits::default()
        });
        let first = admission.connection().unwrap();
        let _second = admission.connection().unwrap();
        assert!(admission.connection().is_none());
        drop(first);
        assert!(admission.connection().is_some());
    }

    #[test]
    fn searches_wait_in_a_bounded_queue_until_a_place_is_free() {
        let admission = Admission::new(ServerLimits {
            connections: 8,
            searches: 1,
            queue: 1,
        });
        let stop = AtomicBool::new(false);
        let running = admission.search(control(&stop)).unwrap();
        std::thread::scope(|scope| {
            let waiting = scope.spawn(|| admission.search(control(&stop)).map(|_| ()));
            while admission.searches() != (1, 1) {
                std::thread::sleep(Duration::from_millis(5));
            }
            assert_eq!(admission.search(control(&stop)).err(), Some(NoPlace::Busy));
            drop(running);
            assert_eq!(waiting.join().unwrap(), Ok(()));
        });
        assert_eq!(admission.searches(), (0, 0));
    }

    #[test]
    fn a_waiting_search_leaves_at_its_deadline_or_stop() {
        let admission = Admission::new(ServerLimits {
            connections: 8,
            searches: 1,
            queue: 4,
        });
        let stop = AtomicBool::new(false);
        let _running = admission.search(control(&stop)).unwrap();
        let deadline = SearchControl {
            deadline: Some(Instant::now() + Duration::from_millis(100)),
            stop: &stop,
        };
        assert_eq!(admission.search(deadline).err(), Some(NoPlace::Stopped));
        let asked = AtomicBool::new(true);
        assert_eq!(
            admission.search(control(&asked)).err(),
            Some(NoPlace::Stopped)
        );
        assert_eq!(admission.searches(), (1, 0));
    }
}
