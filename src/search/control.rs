//! How a running search is stopped: a deadline and a stop flag, checked between steps.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;
use thiserror::Error;

/// Why a search stopped before it finished.
#[derive(Clone, Copy, Debug, Error, Eq, PartialEq)]
pub enum Stopped {
    #[error("the search passed its deadline")]
    Deadline,
    #[error("the search was asked to stop")]
    Asked,
}

/// The deadline and the stop flag of one search. Every stage calls `check` between steps; a read
/// or an alignment that has started runs to its end.
#[derive(Clone, Copy, Debug)]
pub struct SearchControl<'a> {
    pub deadline: Option<Instant>,
    pub stop: &'a AtomicBool,
}

impl SearchControl<'_> {
    pub fn check(&self) -> Result<(), Stopped> {
        if self.stop.load(Ordering::Relaxed) {
            return Err(Stopped::Asked);
        }
        if self
            .deadline
            .is_some_and(|deadline| Instant::now() >= deadline)
        {
            return Err(Stopped::Deadline);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    #[test]
    fn the_stop_flag_and_the_deadline_stop_a_search() {
        let stop = AtomicBool::new(false);
        let mut control = SearchControl {
            deadline: None,
            stop: &stop,
        };
        assert_eq!(control.check(), Ok(()));
        control.deadline = Some(Instant::now() - Duration::from_secs(1));
        assert_eq!(control.check(), Err(Stopped::Deadline));
        stop.store(true, Ordering::Relaxed);
        assert_eq!(control.check(), Err(Stopped::Asked));
    }
}
