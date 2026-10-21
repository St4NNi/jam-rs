//! Searches of several shards at the same time, whose results are written in shard order.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use std::collections::VecDeque;
use std::thread::{self, ScopedJoinHandle};

/// Runs `run(index)` for each index of `0..count` on a thread of its own, at most `at_once` at
/// the same time, and calls `write(index, result)` in index order. A run that panicked gives
/// `write` an error. The first error of `write` is returned once the running ones end.
pub fn in_order<T: Send, E>(
    count: usize,
    at_once: usize,
    run: impl Fn(usize) -> T + Sync,
    mut write: impl FnMut(usize, thread::Result<T>) -> Result<(), E>,
) -> Result<(), E> {
    let run = &run;
    thread::scope(|scope| {
        let mut running: VecDeque<(usize, ScopedJoinHandle<'_, T>)> = VecDeque::new();
        for index in 0..count {
            if running.len() >= at_once.max(1) {
                let (first, handle) = running.pop_front().expect("one runs");
                write(first, handle.join())?;
            }
            running.push_back((index, scope.spawn(move || run(index))));
        }
        for (index, handle) in running {
            write(index, handle.join())?;
        }
        Ok(())
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::time::Duration;

    /// Later runs end first, yet the results come in order, and at most `at_once` run at once.
    #[test]
    fn results_in_order() {
        let (active, peak) = (AtomicUsize::new(0), AtomicUsize::new(0));
        let mut seen = Vec::new();
        let run = |index: usize| {
            let now = active.fetch_add(1, Ordering::SeqCst) + 1;
            peak.fetch_max(now, Ordering::SeqCst);
            thread::sleep(Duration::from_millis(20 - 2 * index as u64));
            active.fetch_sub(1, Ordering::SeqCst);
            index * 10
        };
        let write = |index, result: thread::Result<usize>| {
            seen.push((index, result.unwrap()));
            Ok::<_, ()>(())
        };
        in_order(8, 3, run, write).unwrap();
        let expected: Vec<(usize, usize)> = (0..8).map(|index| (index, index * 10)).collect();
        assert_eq!(seen, expected);
        assert_eq!(peak.load(Ordering::SeqCst), 3);
        let failed = in_order(4, 2, |index| index, |index, _| Err(index));
        assert_eq!(failed, Err(0));
    }
}
