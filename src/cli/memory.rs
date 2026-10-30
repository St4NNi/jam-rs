//! The RAM budget of a command when `--memory` does not give it: a share of the RAM that its
//! cgroup and the machine leave to the process.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use std::path::{Path, PathBuf};

/// RAM that one search takes besides the pins and the block cache, at most: 0.3 to 1.5 GB were
/// measured with 8 threads.
const SEARCH_BYTES: u64 = 3 << 29;
/// RAM that the open index files and the process take besides the searches, the pins and the
/// block cache: the decoded bucket tables of the ten shards of the AllTheBacteria index hold
/// 5.1 GiB alone.
const PROCESS_BYTES: u64 = 8 << 30;
/// The least default budget.
const LEAST: u64 = 1 << 30;

/// The default budget for pins, the block cache and the buffers of reads: half the RAM available
/// to the process, less what the process and `searches` searches at the same time take, but at
/// least 1 GiB.
pub(super) fn default_budget(searches: usize) -> u64 {
    budget_of(available(), searches)
}

fn budget_of(available: Option<u64>, searches: usize) -> u64 {
    let Some(available) = available else {
        return LEAST;
    };
    let taken = SEARCH_BYTES
        .saturating_mul(searches.max(1) as u64)
        .saturating_add(PROCESS_BYTES);
    (available / 2)
        .min(available.saturating_sub(taken))
        .max(LEAST)
}

/// The RAM available to the process: the lowest limit of its cgroup and the cgroups above it,
/// or the RAM the machine has available, whichever is less.
pub(super) fn available() -> Option<u64> {
    let machine = std::fs::read_to_string("/proc/meminfo").ok();
    let machine = machine.as_deref().and_then(mem_available);
    let cgroup = std::fs::read_to_string("/proc/self/cgroup").ok();
    let cgroup = cgroup
        .as_deref()
        .and_then(|text| cgroup_limit(text, Path::new(ROOT)));
    match (machine, cgroup) {
        (Some(machine), Some(cgroup)) => Some(machine.min(cgroup)),
        (machine, cgroup) => machine.or(cgroup),
    }
}

/// Where the cgroup file system is mounted.
const ROOT: &str = "/sys/fs/cgroup";

/// The `MemAvailable` line of `/proc/meminfo`, in bytes.
fn mem_available(meminfo: &str) -> Option<u64> {
    let line = meminfo
        .lines()
        .find(|line| line.starts_with("MemAvailable:"))?;
    let kib: u64 = line.split_whitespace().nth(1)?.parse().ok()?;
    kib.checked_mul(1024)
}

/// The lowest memory limit of the cgroup in `/proc/self/cgroup` (`own`) and those above it,
/// below `root`: `memory.max` of cgroup v2, else `memory.limit_in_bytes` of v1.
fn cgroup_limit(own: &str, root: &Path) -> Option<u64> {
    let (base, file, path) = match own.lines().find_map(|line| line.strip_prefix("0::")) {
        Some(path) => (root.to_path_buf(), "memory.max", path),
        None => {
            let line = own
                .lines()
                .find(|line| line.split(':').nth(1) == Some("memory"))?;
            let path = line.splitn(3, ':').nth(2)?;
            (root.join("memory"), "memory.limit_in_bytes", path)
        }
    };
    let mut directory: PathBuf = base.join(path.trim_start_matches('/'));
    let mut least: Option<u64> = None;
    loop {
        let limit = std::fs::read_to_string(directory.join(file)).ok();
        // "max" means no limit; v1 writes a huge number instead.
        if let Some(limit) = limit.and_then(|text| text.trim().parse::<u64>().ok()) {
            least = Some(least.map_or(limit, |least| least.min(limit)));
        }
        if directory == base || !directory.pop() || !directory.starts_with(&base) {
            return least;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn budget_by_ram() {
        let gib = 1u64 << 30;
        // 16 GiB: the process takes 8 GiB, so one search leaves 6.5 GiB, four leave 2 GiB and
        // eight leave the least.
        assert_eq!(budget_of(Some(16 * gib), 1), 13 * gib / 2);
        assert_eq!(budget_of(Some(16 * gib), 4), 2 * gib);
        assert_eq!(budget_of(Some(16 * gib), 8), gib);
        // 96 GiB: half, as the process and four searches take only 14 GiB.
        assert_eq!(budget_of(Some(96 * gib), 4), 48 * gib);
        // Little RAM, or none known: 1 GiB.
        assert_eq!(budget_of(Some(2 * gib), 4), gib);
        assert_eq!(budget_of(None, 1), gib);
    }

    #[test]
    fn limits_read() {
        let meminfo = "MemTotal:       32000000 kB\nMemAvailable:   20000000 kB\n";
        assert_eq!(mem_available(meminfo), Some(20_000_000 * 1024));
        assert_eq!(mem_available("MemTotal: 1 kB\n"), None);
        let root = tempfile::tempdir().unwrap();
        let job = root.path().join("slice/job");
        std::fs::create_dir_all(&job).unwrap();
        std::fs::write(job.join("memory.max"), "max\n").unwrap();
        std::fs::write(root.path().join("slice/memory.max"), "17179869184\n").unwrap();
        assert_eq!(cgroup_limit("0::/slice/job\n", root.path()), Some(16 << 30));
        assert_eq!(cgroup_limit("0::/other\n", root.path()), None);
        let v1 = root.path().join("memory/job");
        std::fs::create_dir_all(&v1).unwrap();
        std::fs::write(v1.join("memory.limit_in_bytes"), "4294967296\n").unwrap();
        assert_eq!(
            cgroup_limit("4:memory:/job\n1:cpu:/\n", root.path()),
            Some(4 << 30)
        );
    }
}
