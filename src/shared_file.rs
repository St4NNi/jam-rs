use crate::jidx::sha256;
use crate::owner_format::{ChecksumLevel, checksum_layout};
use crate::shared_format::{HEADER_BYTES, PAGE_BYTES, Section, SharedError, SharedHeader};
use memmap2::{Mmap, MmapOptions};
use serde::Serialize;
use std::fs::File;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct FileReadStats {
    pub observed: bool,
    pub identity_checks: u64,
    pub requested_bytes: u64,
    pub requested_pages: u64,
    pub authenticated_pages: u64,
    pub authenticated_bytes: u64,
    /// Page hashes computed, including failed and repeated verifications of one page.
    pub hash_attempts: u64,
    pub resident_integrity_bytes: usize,
}

pub(crate) struct SharedFile {
    observed: bool,
    identity_checks: AtomicU64,
    file: File,
    mmap: Mmap,
    identity: Option<[u64; 7]>,
    pub(crate) header: SharedHeader,
    levels: Vec<ChecksumLevel>,
    verified: Box<[AtomicU64]>,
    requested_bytes: AtomicU64,
    requested_pages: AtomicU64,
    authenticated_pages: AtomicU64,
    hash_attempts: AtomicU64,
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
            identity_checks: AtomicU64::new(1),
            file,
            mmap,
            identity,
            header,
            levels,
            verified: (0..words).map(|_| AtomicU64::new(0)).collect(),
            requested_bytes: AtomicU64::new(HEADER_BYTES as u64),
            requested_pages: AtomicU64::new(1),
            authenticated_pages: AtomicU64::new(1),
            hash_attempts: AtomicU64::new(1),
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
        if self.observed {
            self.identity_checks.fetch_add(1, Ordering::Relaxed);
        }
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
        // Errors are built only when returned; this runs for every index read.
        let Some(end) = offset
            .checked_add(length)
            .filter(|&end| end <= section.length)
        else {
            return Err(SharedError::Invalid("section request"));
        };
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
        match self.mmap.get(start as usize..end as usize) {
            Some(bytes) => Ok(bytes),
            None => Err(SharedError::Invalid("mapped extent")),
        }
    }

    pub(crate) fn record(
        &self,
        kind: Section,
        ordinal: u64,
        size: u64,
    ) -> Result<&[u8], SharedError> {
        if size != self.header.row_bytes(kind) {
            return Err(SharedError::Invalid("record size"));
        }
        let Some(offset) = ordinal.checked_mul(size) else {
            return Err(SharedError::Invalid("record ordinal"));
        };
        self.section(kind, offset, size)
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
        self.authenticate_page(page)
    }

    /// Verifies a page against its parent hash without a lock: the mapped generation is
    /// immutable, verification is idempotent and the verified bit is set atomically. Each
    /// caller that finds the bit clear hashes the page once, so N racing callers may hash it
    /// up to N times; callers that start after the bit is published do not hash it again. A
    /// failed verification never sets the bit, so every later call hashes the page again.
    /// Without a file identity no bit is kept, and every call hashes the page.
    fn authenticate_page(&self, page: u64) -> Result<(), SharedError> {
        if self.known(page) {
            return Ok(());
        }
        let checksums = self.header.section(Section::Checksums);
        let byte = page * PAGE_BYTES;
        let expected: [u8; 32] = if byte < checksums.offset {
            let hash_offset = checksums.offset + (page - 1) * 32;
            self.authenticate_page(hash_offset / PAGE_BYTES)?;
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
                self.authenticate_page(hash_offset / PAGE_BYTES)?;
                self.mmap[hash_offset as usize..hash_offset as usize + 32]
                    .try_into()
                    .unwrap()
            } else {
                self.header.checksum_root_sha256
            }
        };
        self.hash_attempts.fetch_add(1, Ordering::Relaxed);
        let digest: [u8; 32] = sha256(&self.mmap[byte as usize..(byte + PAGE_BYTES) as usize]);
        if digest != expected {
            return Err(SharedError::ChecksumMismatch);
        }
        let bit = 1u64 << (page % 64);
        // A page verified concurrently by several readers is counted once.
        if self.identity.is_none()
            || self.verified[page as usize / 64].fetch_or(bit, Ordering::AcqRel) & bit == 0
        {
            self.authenticated_pages.fetch_add(1, Ordering::Relaxed);
        }
        Ok(())
    }

    pub(crate) fn stats(&self) -> FileReadStats {
        let pages = self.authenticated_pages.load(Ordering::Relaxed);
        FileReadStats {
            observed: self.observed,
            identity_checks: self.identity_checks.load(Ordering::Relaxed),
            requested_bytes: self.requested_bytes.load(Ordering::Relaxed),
            requested_pages: self.requested_pages.load(Ordering::Relaxed),
            authenticated_pages: pages,
            authenticated_bytes: pages * PAGE_BYTES,
            hash_attempts: self.hash_attempts.load(Ordering::Relaxed),
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::jidx_reader::JidxReader;
    use crate::jidx_writer::{ContigInput, JidxInput, JidxWriter, MetagenomeInput};
    use crate::shared_seed::SharedSeed;
    use crate::shared_writer::{IndexedSeed, write_shared_index};
    use std::sync::Barrier;

    const WORKERS: usize = 8;

    /// Writes a shared index with two checksum levels, so data pages have a non-root parent.
    fn write_index(directory: &Path, corrupt: impl FnOnce(&mut [u8], u64)) -> std::path::PathBuf {
        let jidx = directory.join("metadata.jidx");
        let mut writer = JidxWriter::new(
            &jidx,
            &JidxInput {
                k: 15,
                rescue_k15: false,
                minimizer_window: 64,
                jam_sha256: [1; 32],
                manifest_sha256: [2; 32],
            },
        )
        .unwrap();
        writer
            .begin_metagenome(MetagenomeInput {
                name: "a".into(),
                bgzf_uri: "a.bgz".into(),
                bgzf_bytes: 100,
                bgzf_sha256: [3; 32],
                gzi: vec![0; 8],
            })
            .unwrap();
        writer
            .begin_contig(ContigInput {
                name: "a-contig".into(),
                length: 100_000,
                fasta_offset: 4,
                line_bases: 80,
                line_width: 81,
            })
            .unwrap();
        writer.finish().unwrap();
        let reference = JidxReader::open(&jidx).unwrap();
        let mut seeds = (0..40_000)
            .map(|core| IndexedSeed {
                member: 0,
                contig: 0,
                seed: SharedSeed {
                    core,
                    context: 0,
                    flags: 0,
                    position: u64::from(core % 90_000),
                },
            })
            .collect::<Vec<_>>();
        let path = directory.join("pages.shared");
        write_shared_index(&reference, &path, 64, &mut seeds).unwrap();
        let mut bytes = std::fs::read(&path).unwrap();
        let header = SharedHeader::decode(&bytes[..HEADER_BYTES], bytes.len() as u64).unwrap();
        let checksums = header.section(Section::Checksums).offset;
        assert!(checksum_layout(checksums / PAGE_BYTES - 1).unwrap().len() >= 2);
        corrupt(&mut bytes, checksums);
        std::fs::write(&path, bytes).unwrap();
        path
    }

    /// Runs one authentication of `page` on each worker, all released together.
    fn race(file: &SharedFile, page: u64) -> Vec<Result<(), SharedError>> {
        let barrier = Barrier::new(WORKERS);
        std::thread::scope(|scope| {
            let workers = (0..WORKERS)
                .map(|_| {
                    scope.spawn(|| {
                        barrier.wait();
                        file.authenticate(page)
                    })
                })
                .collect::<Vec<_>>();
            workers
                .into_iter()
                .map(|worker| worker.join().unwrap())
                .collect()
        })
    }

    #[cfg(unix)]
    #[test]
    fn racing_first_verification_publishes_one_page_for_other_workers() {
        let directory = tempfile::tempdir().unwrap();
        let path = write_index(directory.path(), |_, _| {});
        let file = SharedFile::open(&path, false).unwrap();
        // Page 1 verifies the checksum parent that page 2 shares.
        file.authenticate(1).unwrap();
        let before = file.stats();
        assert!(race(&file, 2).iter().all(Result::is_ok));
        let raced = file.stats();
        assert!(file.known(2));
        assert_eq!(raced.authenticated_pages - before.authenticated_pages, 1);
        let hashes = raced.hash_attempts - before.hash_attempts;
        assert!((1..=WORKERS as u64).contains(&hashes), "{hashes}");

        std::thread::scope(|scope| {
            scope.spawn(|| file.authenticate(2).unwrap());
        });
        assert_eq!(file.stats().hash_attempts, raced.hash_attempts);
        assert_eq!(file.stats().authenticated_pages, raced.authenticated_pages);
    }

    #[test]
    fn corrupted_data_page_is_hashed_again_and_never_published() {
        let directory = tempfile::tempdir().unwrap();
        let path = write_index(directory.path(), |bytes, _| {
            bytes[2 * PAGE_BYTES as usize + 17] ^= 1;
        });
        let file = SharedFile::open(&path, false).unwrap();
        file.authenticate(1).unwrap();
        let before = file.stats();
        assert!(
            race(&file, 2)
                .iter()
                .all(|result| matches!(result, Err(SharedError::ChecksumMismatch)))
        );
        let failed = file.stats();
        assert!(!file.known(2));
        assert_eq!(failed.authenticated_pages, before.authenticated_pages);
        // Without a published bit, every racing caller hashes the page.
        assert_eq!(failed.hash_attempts - before.hash_attempts, WORKERS as u64);

        let (kind, section) = file
            .header
            .section_order()
            .iter()
            .map(|&kind| (kind, file.header.section(kind)))
            .find(|(_, section)| {
                (section.offset..section.offset + section.length).contains(&(2 * PAGE_BYTES))
            })
            .unwrap();
        assert!(matches!(
            file.section(kind, 2 * PAGE_BYTES - section.offset, 1),
            Err(SharedError::ChecksumMismatch)
        ));
        assert!(!file.known(2));
        assert_eq!(file.stats().hash_attempts, failed.hash_attempts + 1);
        assert_eq!(file.stats().authenticated_pages, before.authenticated_pages);
    }

    #[test]
    fn corrupted_checksum_parent_fails_its_data_pages_on_every_read() {
        let directory = tempfile::tempdir().unwrap();
        let path = write_index(directory.path(), |bytes, checksums| {
            bytes[checksums as usize + 5] ^= 1;
        });
        // Opening verifies only the root, which the corrupted first-level page does not change.
        let file = SharedFile::open(&path, false).unwrap();
        let parent = file.header.section(Section::Checksums).offset / PAGE_BYTES;
        let sibling = 1 + 128;
        assert!(sibling < parent);
        let before = file.stats();
        assert!(
            race(&file, 1)
                .iter()
                .all(|result| matches!(result, Err(SharedError::ChecksumMismatch)))
        );
        let failed = file.stats();
        assert!(!file.known(1) && !file.known(parent));
        assert_eq!(failed.authenticated_pages, before.authenticated_pages);
        // Each caller hashes only the parent, which fails before the data page is hashed.
        assert_eq!(failed.hash_attempts - before.hash_attempts, WORKERS as u64);
        assert!(matches!(
            file.authenticate(2),
            Err(SharedError::ChecksumMismatch)
        ));
        assert_eq!(file.stats().hash_attempts, failed.hash_attempts + 1);

        // A data page under an intact first-level page still verifies.
        file.authenticate(sibling).unwrap();
        assert_eq!(file.known(sibling), file.identity().is_some());
        assert!(!file.known(parent));
    }
}
