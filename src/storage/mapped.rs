//! A local file mapped into memory and read by offset.

use super::{ReadAt, StorageError};
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
        let file = std::fs::File::open(path)?;
        if file.metadata()?.len() == 0 {
            return Ok(Self { map: None });
        }
        let mut options = memmap2::MmapOptions::new();
        if populate {
            options.populate();
        }
        // SAFETY: the files are written once and never changed while a shard is open.
        let map = unsafe { options.map(&file)? };
        Ok(Self { map: Some(map) })
    }

    /// All bytes of the file.
    pub fn bytes(&self) -> &[u8] {
        self.map.as_deref().unwrap_or_default()
    }
}

impl ReadAt for MappedFile {
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> io::Result<()> {
        let bytes = usize::try_from(offset)
            .ok()
            .and_then(|start| self.bytes().get(start..start.checked_add(buffer.len())?))
            .ok_or_else(|| io::Error::new(io::ErrorKind::UnexpectedEof, "read past the end"))?;
        buffer.copy_from_slice(bytes);
        Ok(())
    }

    fn size(&self) -> u64 {
        self.bytes().len() as u64
    }
}
