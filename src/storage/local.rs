//! A local file read by offset.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::{ReadAt, StorageError, checked_end, readahead};
use std::fs::File;
use std::io;
use std::path::Path;

/// A local file with the size it had when it was opened. Reads are positional, so threads never
/// share a file position. The file is advised for random access, so reads are exact and
/// read-ahead hints do not wait for the data (measured on CephFS).
pub struct LocalFile {
    file: File,
    size: u64,
    /// Whether the file lies on CephFS, where each read costs much system time.
    on_ceph: bool,
}

/// The `f_type` of CephFS in `statfs`.
#[cfg(target_os = "linux")]
const CEPH_SUPER_MAGIC: i64 = 0x00c3_6400;

impl LocalFile {
    /// Opens the file at `path` for reading.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, StorageError> {
        let path = path.as_ref();
        let failed = |source| StorageError::File {
            path: path.to_path_buf(),
            source,
        };
        let file = File::open(path).map_err(failed)?;
        let size = file.metadata().map_err(failed)?.len();
        #[cfg(target_os = "linux")]
        let on_ceph = {
            use std::os::fd::AsRawFd;
            // SAFETY: the descriptor is open, and the call only sets advice.
            unsafe { libc::posix_fadvise(file.as_raw_fd(), 0, 0, libc::POSIX_FADV_RANDOM) };
            // SAFETY: the descriptor is open, and `stats` is written by the call before use.
            let mut stats: libc::statfs = unsafe { std::mem::zeroed() };
            let known = unsafe { libc::fstatfs(file.as_raw_fd(), &mut stats) } == 0;
            known && stats.f_type as i64 == CEPH_SUPER_MAGIC
        };
        #[cfg(not(target_os = "linux"))]
        let on_ceph = false;
        Ok(Self {
            file,
            size,
            on_ceph,
        })
    }
}

impl ReadAt for LocalFile {
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> io::Result<()> {
        checked_end(offset, buffer.len(), self.size)?;
        std::os::unix::fs::FileExt::read_exact_at(&self.file, buffer, offset)
    }

    fn size(&self) -> u64 {
        self.size
    }

    fn will_need(&self, offset: u64, length: u64) {
        readahead::will_need(&self.file, offset, length);
    }

    fn costly_reads(&self) -> bool {
        self.on_ceph
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_stop_end() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("data");
        let data: Vec<u8> = (0..=255).cycle().take(10_000).collect();
        std::fs::write(&path, &data).unwrap();
        let file = LocalFile::open(&path).unwrap();
        assert_eq!(file.size(), 10_000);
        let mut buffer = vec![0; 300];
        file.read_exact_at(9_700, &mut buffer).unwrap();
        assert_eq!(buffer, data[9_700..]);
        assert!(file.read_exact_at(9_701, &mut buffer).is_err());
        file.will_need(0, 4096);
        // A temporary directory is not on CephFS.
        assert!(!file.costly_reads());
        assert!(LocalFile::open(directory.path().join("missing")).is_err());
    }

    #[test]
    fn threads_share_file() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("data");
        let data: Vec<u8> = (0..=250).cycle().take(100_000).collect();
        std::fs::write(&path, &data).unwrap();
        let file = LocalFile::open(&path).unwrap();
        std::thread::scope(|scope| {
            for thread in 0..8usize {
                let (file, data) = (&file, &data);
                scope.spawn(move || {
                    for step in 0..200usize {
                        let start = (thread * 7919 + step * 331) % 99_000;
                        let mut out = vec![0; 500];
                        file.read_exact_at(start as u64, &mut out).unwrap();
                        assert_eq!(out, data[start..start + 500]);
                    }
                });
            }
        });
    }
}
