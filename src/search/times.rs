//! Wall seconds, CPU seconds, major page faults and storage counts of the parts of one search,
//! so time spent waiting for storage shows up next to the work.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use crate::storage::StorageCounts;
use serde::{Deserialize, Serialize};
use std::time::Instant;

/// What one part of a search cost. CPU seconds and faults count the whole process, all threads,
/// while the part ran; wall seconds well above CPU seconds over threads mean waiting.
#[derive(Clone, Copy, Debug, Default, Deserialize, PartialEq, Serialize)]
pub struct PartTime {
    pub wall_seconds: f64,
    pub cpu_seconds: f64,
    /// Page faults that had to read from storage.
    pub major_faults: u64,
    /// Requests to S3, block cache lookups and pinned reads of the process while the part ran.
    pub storage: StorageCounts,
}

/// The cost of the parts of one search, in the order they run.
#[derive(Clone, Copy, Debug, Default, Deserialize, PartialEq, Serialize)]
pub struct SearchTimes {
    /// Opening the shard, once per process before its first search.
    pub open_shard_seconds: f64,
    /// Hashing the query k-mers and building the query tables of stage 3.
    pub prepare_query: PartTime,
    /// Stage 1.
    pub sketch: PartTime,
    /// Stages 2 and 3: the seed locations of all candidate entries are read first, then each
    /// worker takes an entry through the rest of both.
    pub entries: PartTime,
}

/// The start of a part: the clock and the process usage at that moment.
pub(super) struct PartStart {
    wall: Instant,
    cpu_seconds: f64,
    major_faults: u64,
    storage: StorageCounts,
}

impl PartStart {
    /// Starts a part now.
    pub(super) fn now() -> Self {
        let (cpu_seconds, major_faults) = process_usage();
        Self {
            wall: Instant::now(),
            cpu_seconds,
            major_faults,
            storage: StorageCounts::now(),
        }
    }

    /// The cost of the part from its start until now.
    pub(super) fn finish(&self) -> PartTime {
        let (cpu_seconds, major_faults) = process_usage();
        PartTime {
            wall_seconds: self.wall.elapsed().as_secs_f64(),
            cpu_seconds: (cpu_seconds - self.cpu_seconds).max(0.0),
            major_faults: major_faults.saturating_sub(self.major_faults),
            storage: StorageCounts::now().since(&self.storage),
        }
    }
}

/// CPU seconds (user and system) and major page faults of the process so far.
fn process_usage() -> (f64, u64) {
    let mut usage = std::mem::MaybeUninit::<libc::rusage>::zeroed();
    // SAFETY: getrusage fills the struct it is given and only reads the process counters.
    let usage = unsafe {
        if libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) != 0 {
            return (0.0, 0);
        }
        usage.assume_init()
    };
    let seconds = |time: libc::timeval| time.tv_sec as f64 + time.tv_usec as f64 / 1_000_000.0;
    let faults = u64::try_from(usage.ru_majflt).unwrap_or_default();
    (seconds(usage.ru_utime) + seconds(usage.ru_stime), faults)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn part_counts_work() {
        let start = PartStart::now();
        let mut sum = 0u64;
        for value in 0..2_000_000u64 {
            sum = sum.wrapping_add(value.wrapping_mul(value));
        }
        std::hint::black_box(sum);
        let part = start.finish();
        assert!(part.wall_seconds > 0.0);
        assert!(part.cpu_seconds >= 0.0);
    }
}
