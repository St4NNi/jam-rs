//! A local file mapped into memory, or a copy of a source in anonymous memory, read by offset.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::{ReadAt, StorageError, read_slice};
use rayon::prelude::*;
use std::io;
use std::path::Path;

/// Reads of at least this many bytes ask the kernel for all their pages at once, so they are read
/// in large requests instead of one page fault after another.
const LARGE_READ: usize = 256 << 10;

/// Bytes per read when a source is copied into memory; the reads run on all threads.
const COPY_CHUNK: usize = 8 << 20;

/// A local file mapped into memory, or a copy of a source in anonymous memory, read by offset.
/// It suits files that are read whole or often, such as the sketch.
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

    /// Copies all bytes of `source` into anonymous memory, in large reads on all threads.
    pub fn read_from(source: &dyn ReadAt) -> io::Result<Self> {
        let size = usize::try_from(source.size()).map_err(io::Error::other)?;
        if size == 0 {
            return Ok(Self { map: None });
        }
        let mut map = memmap2::MmapMut::map_anon(size)?;
        map.par_chunks_mut(COPY_CHUNK)
            .enumerate()
            .try_for_each(|(index, chunk)| {
                source.read_exact_at((index * COPY_CHUNK) as u64, chunk)
            })?;
        Ok(Self {
            map: Some(map.make_read_only()?),
        })
    }

    /// Reads every page of the file into RAM now, so later reads never wait for storage.
    pub fn load(&self) -> io::Result<()> {
        match &self.map {
            Some(map) => map.advise(memmap2::Advice::PopulateRead),
            None => Ok(()),
        }
    }

    /// All bytes of the file.
    pub fn bytes(&self) -> &[u8] {
        self.map.as_deref().unwrap_or_default()
    }
}

impl ReadAt for MappedFile {
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> io::Result<()> {
        if buffer.len() >= LARGE_READ {
            self.will_need(offset, buffer.len() as u64);
        }
        read_slice(self.bytes(), offset, buffer)
    }

    fn size(&self) -> u64 {
        self.bytes().len() as u64
    }

    fn will_need(&self, offset: u64, length: u64) {
        let Some(map) = &self.map else {
            return;
        };
        let (Ok(offset), Ok(length)) = (usize::try_from(offset), usize::try_from(length)) else {
            return;
        };
        let length = length.min(map.len().saturating_sub(offset));
        if length > 0 {
            // A hint only: a failed one leaves the pages to be read when they are touched.
            let _ = map.advise_range(memmap2::Advice::WillNeed, offset, length);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hints_keep_bytes() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("data");
        let data: Vec<u8> = (0..=255).cycle().take(LARGE_READ + 5_000).collect();
        std::fs::write(&path, &data).unwrap();
        let file = MappedFile::open(&path, false).unwrap();
        file.will_need(0, 100);
        file.will_need(u64::MAX, u64::MAX);
        file.will_need(4_000, u64::MAX);
        let mut buffer = vec![0; LARGE_READ + 1_000];
        file.read_exact_at(3_000, &mut buffer).unwrap();
        assert_eq!(buffer, data[3_000..LARGE_READ + 4_000]);
        assert!(file.read_exact_at(5_000, &mut buffer).is_err());
        let copy = MappedFile::read_from(&data).unwrap();
        assert_eq!(copy.bytes(), data);
        assert!(
            MappedFile::read_from(&Vec::new())
                .unwrap()
                .bytes()
                .is_empty()
        );
    }
}
