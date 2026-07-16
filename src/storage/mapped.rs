//! A local file mapped into memory and read by offset.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::{ReadAt, StorageError, read_slice};
use std::io;
use std::path::Path;

/// A local file mapped into memory and read by offset. It suits files that are read whole or
/// often, such as the sketch.
pub struct MappedFile {
    map: Option<memmap2::Mmap>,
}

impl MappedFile {
    /// Maps the file at `path`; `populate` reads it into RAM at once.
    pub fn open(path: &Path, populate: bool) -> Result<Self, StorageError> {
        let failed = |source| StorageError::File {
            path: path.to_path_buf(),
            source,
        };
        let file = std::fs::File::open(path).map_err(failed)?;
        if file.metadata().map_err(failed)?.len() == 0 {
            return Ok(Self { map: None });
        }
        let mut options = memmap2::MmapOptions::new();
        if populate {
            options.populate();
        }
        // SAFETY: the files are written once and never changed while a shard is open.
        let map = unsafe { options.map(&file).map_err(failed)? };
        Ok(Self { map: Some(map) })
    }

    /// All bytes of the file.
    pub fn bytes(&self) -> &[u8] {
        self.map.as_deref().unwrap_or_default()
    }
}

impl ReadAt for MappedFile {
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> io::Result<()> {
        read_slice(self.bytes(), offset, buffer)
    }

    fn size(&self) -> u64 {
        self.bytes().len() as u64
    }
}
