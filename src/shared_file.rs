use crate::jidx::sha256;
use crate::owner_format::{ChecksumLevel, checksum_layout};
use crate::shared_format::{HEADER_BYTES, PAGE_BYTES, Section, SharedError, SharedHeader};
use memmap2::{Mmap, MmapOptions};
use serde::Serialize;
use std::fs::File;
use std::path::Path;
use std::sync::Mutex;
use std::sync::atomic::{AtomicU64, Ordering};

#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct FileReadStats {
    pub observed: bool,
    pub requested_bytes: u64,
    pub requested_pages: u64,
    pub authenticated_pages: u64,
    pub authenticated_bytes: u64,
    pub resident_integrity_bytes: usize,
}

pub(crate) struct SharedFile {
    observed: bool,
    file: File,
    mmap: Mmap,
    identity: Option<[u64; 7]>,
    pub(crate) header: SharedHeader,
    levels: Vec<ChecksumLevel>,
    verified: Box<[AtomicU64]>,
    authentication: Mutex<()>,
    requested_bytes: AtomicU64,
    requested_pages: AtomicU64,
    authenticated_pages: AtomicU64,
}

impl SharedFile {
    pub(crate) fn open(path: impl AsRef<Path>, observed: bool) -> Result<Self, SharedError> {
        let file = File::open(path)?;
        let identity = file_identity(&file)?;
        let file_bytes = file.metadata()?.len();
        if file_bytes < HEADER_BYTES as u64 {
            return Err(SharedError::Invalid("file size"));
        }
        // SAFETY: shared generations are immutable; each operation checks retained file identity.
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        let header = SharedHeader::decode(&mmap[..HEADER_BYTES], file_bytes)?;
        let checksum_start = header.section(Section::Checksums).offset;
        let levels = checksum_layout(checksum_start / PAGE_BYTES - 1)
            .map_err(|_| SharedError::Invalid("checksum layout"))?;
        let words = usize::try_from(file_bytes.div_ceil(PAGE_BYTES).div_ceil(64))
            .map_err(|_| SharedError::ResourceLimit)?;
        if words > 1024 * 1024 {
            return Err(SharedError::ResourceLimit);
        }
        let reader = Self {
            observed,
            file,
            mmap,
            identity,
            header,
            levels,
            verified: (0..words).map(|_| AtomicU64::new(0)).collect(),
            authentication: Mutex::new(()),
            requested_bytes: AtomicU64::new(HEADER_BYTES as u64),
            requested_pages: AtomicU64::new(1),
            authenticated_pages: AtomicU64::new(1),
        };
        reader.verify_unchanged()?;
        let top = reader
            .levels
            .last()
            .ok_or(SharedError::Invalid("checksum root"))?;
        reader.authenticate((checksum_start + top.offset) / PAGE_BYTES)?;
        reader.verify_unchanged()?;
        Ok(reader)
    }

    pub(crate) fn identity(&self) -> Option<[u64; 7]> {
        self.identity
    }

    pub(crate) fn verify_unchanged(&self) -> Result<(), SharedError> {
        if self.identity.is_some() && file_identity(&self.file)? != self.identity {
            return Err(SharedError::SourceChanged);
        }
        Ok(())
    }

    pub(crate) fn section(
        &self,
        kind: Section,
        offset: u64,
        length: u64,
    ) -> Result<&[u8], SharedError> {
        let section = self.header.section(kind);
        let end = offset
            .checked_add(length)
            .filter(|&end| end <= section.length)
            .ok_or(SharedError::Invalid("section request"))?;
        let start = section.offset + offset;
        let end = section.offset + end;
        if self.observed {
            self.requested_bytes.fetch_add(length, Ordering::Relaxed);
        }
        if length != 0 {
            let first = start / PAGE_BYTES;
            let last = (end - 1) / PAGE_BYTES;
            if self.observed {
                self.requested_pages
                    .fetch_add(last - first + 1, Ordering::Relaxed);
            }
            for page in first..=last {
                self.authenticate(page)?;
            }
        }
        self.mmap
            .get(start as usize..end as usize)
            .ok_or(SharedError::Invalid("mapped extent"))
    }

    pub(crate) fn record(
        &self,
        kind: Section,
        ordinal: u64,
        size: u64,
    ) -> Result<&[u8], SharedError> {
        if size != kind.row_bytes() {
            return Err(SharedError::Invalid("record size"));
        }
        self.section(
            kind,
            ordinal
                .checked_mul(size)
                .ok_or(SharedError::Invalid("record ordinal"))?,
            size,
        )
    }

    fn known(&self, page: u64) -> bool {
        self.identity.is_some()
            && self.verified[page as usize / 64].load(Ordering::Acquire) & (1u64 << (page % 64))
                != 0
    }

    fn authenticate(&self, page: u64) -> Result<(), SharedError> {
        if page >= self.mmap.len() as u64 / PAGE_BYTES || page == 0 {
            return Err(SharedError::Invalid("authentication page"));
        }
        if self.known(page) {
            return Ok(());
        }
        let _guard = self
            .authentication
            .lock()
            .map_err(|_| SharedError::Invalid("authentication state"))?;
        self.authenticate_locked(page)
    }

    fn authenticate_locked(&self, page: u64) -> Result<(), SharedError> {
        if self.known(page) {
            return Ok(());
        }
        let checksums = self.header.section(Section::Checksums);
        let byte = page * PAGE_BYTES;
        let expected: [u8; 32] = if byte < checksums.offset {
            let hash_offset = checksums.offset + (page - 1) * 32;
            self.authenticate_locked(hash_offset / PAGE_BYTES)?;
            self.mmap[hash_offset as usize..hash_offset as usize + 32]
                .try_into()
                .unwrap()
        } else {
            let (index, level) = self
                .levels
                .iter()
                .enumerate()
                .find(|(_, level)| {
                    byte >= checksums.offset + level.offset
                        && byte < checksums.offset + level.offset + level.page_count * PAGE_BYTES
                })
                .ok_or(SharedError::Invalid("checksum page"))?;
            if let Some(parent) = self.levels.get(index + 1) {
                let page_index = (byte - checksums.offset - level.offset) / PAGE_BYTES;
                let hash_offset = checksums.offset + parent.offset + page_index * 32;
                self.authenticate_locked(hash_offset / PAGE_BYTES)?;
                self.mmap[hash_offset as usize..hash_offset as usize + 32]
                    .try_into()
                    .unwrap()
            } else {
                self.header.checksum_root_sha256
            }
        };
        let digest: [u8; 32] = sha256(&self.mmap[byte as usize..(byte + PAGE_BYTES) as usize]);
        if digest != expected {
            return Err(SharedError::ChecksumMismatch);
        }
        self.authenticated_pages.fetch_add(1, Ordering::Relaxed);
        if self.identity.is_some() {
            self.verified[page as usize / 64].fetch_or(1u64 << (page % 64), Ordering::Release);
        }
        Ok(())
    }

    pub(crate) fn stats(&self) -> FileReadStats {
        let pages = self.authenticated_pages.load(Ordering::Relaxed);
        FileReadStats {
            observed: self.observed,
            requested_bytes: self.requested_bytes.load(Ordering::Relaxed),
            requested_pages: self.requested_pages.load(Ordering::Relaxed),
            authenticated_pages: pages,
            authenticated_bytes: pages * PAGE_BYTES,
            resident_integrity_bytes: self.verified.len() * 8
                + self.levels.capacity() * std::mem::size_of::<ChecksumLevel>(),
        }
    }

    pub(crate) fn verify_checksum(&self) -> Result<(), SharedError> {
        self.verify_unchanged()?;
        if sha256(&self.mmap[HEADER_BYTES..]) != self.header.body_sha256 {
            return Err(SharedError::ChecksumMismatch);
        }
        self.verify_unchanged()
    }
}

fn file_identity(file: &File) -> std::io::Result<Option<[u64; 7]>> {
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        let m = file.metadata()?;
        Ok(Some([
            m.dev(),
            m.ino(),
            m.len(),
            m.mtime() as u64,
            m.mtime_nsec() as u64,
            m.ctime() as u64,
            m.ctime_nsec() as u64,
        ]))
    }
    #[cfg(not(unix))]
    {
        let _ = file;
        Ok(None)
    }
}
