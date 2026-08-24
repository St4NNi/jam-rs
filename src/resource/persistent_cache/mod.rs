//! Optional process-persistent cache for immutable resource pages.
//!
//! Cache entries are addressed only by the immutable object identity, section,
//! page, exact byte range, and page checksum.  The on-disk name is a digest of
//! those fields; locators and other caller metadata are deliberately not
//! persisted.  Entries use a fixed binary header so a partial write can be
//! discarded without interpreting arbitrary metadata.

use super::ByteRange;
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
use std::fs::{self, File, OpenOptions};
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicU64, Ordering};
use thiserror::Error;

const MAGIC: &[u8; 8] = b"JAMPCG01";
const VERSION: u16 = 1;
const HEADER_BYTES: usize = 120;
const ENTRY_SUFFIX: &str = ".jpc";
static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

/// Errors returned while opening or accessing the persistent page cache.
#[derive(Debug, Error)]
pub enum PersistentCacheError {
    #[error("persistent cache {operation} failed: {source}")]
    Io {
        operation: &'static str,
        #[source]
        source: io::Error,
    },
    #[error("persistent cache key is invalid: {0}")]
    InvalidKey(&'static str),
    #[error("persistent cache entry is invalid")]
    InvalidEntry,
    #[error("persistent cache page checksum does not match the key")]
    PageChecksumMismatch,
    #[error("persistent cache accounting overflow")]
    AccountingOverflow,
}

fn io_error(operation: &'static str, source: io::Error) -> PersistentCacheError {
    PersistentCacheError::Io { operation, source }
}

/// Exact identity of one immutable page.
///
/// The object checksum prevents bytes from another object version being used;
/// the page checksum protects against a stale or corrupted page.  `range` is
/// retained in the key even when callers happen to use fixed-size pages so
/// that different physical layouts cannot alias.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct PageKey {
    pub object_checksum: [u8; 32],
    pub section_id: u32,
    pub page_id: u64,
    pub range: ByteRange,
    pub page_checksum: [u8; 32],
}

impl PageKey {
    /// Construct and validate a page key from its exact range.
    pub fn new(
        object_checksum: [u8; 32],
        section_id: u32,
        page_id: u64,
        range: ByteRange,
        page_checksum: [u8; 32],
    ) -> Result<Self, PersistentCacheError> {
        range
            .offset
            .checked_add(range.length)
            .ok_or(PersistentCacheError::InvalidKey("range overflows"))?;
        Ok(Self {
            object_checksum,
            section_id,
            page_id,
            range,
            page_checksum,
        })
    }

    fn logical_eq(self, other: Self) -> bool {
        self.section_id == other.section_id
            && self.page_id == other.page_id
            && self.range == other.range
    }

    fn digest(self) -> [u8; 32] {
        let mut hasher = Sha256::new();
        hasher.update(b"jam-rs-persistent-page-cache-v1");
        hasher.update(self.object_checksum);
        hasher.update(self.section_id.to_le_bytes());
        hasher.update(self.page_id.to_le_bytes());
        hasher.update(self.range.offset.to_le_bytes());
        hasher.update(self.range.length.to_le_bytes());
        hasher.update(self.page_checksum);
        hasher.finalize().into()
    }

    /// Return the non-secret cache file path for this key.
    pub fn path(&self, cache_dir: impl AsRef<Path>) -> PathBuf {
        cache_dir
            .as_ref()
            .join(format!("{}{}", hex(&self.digest()), ENTRY_SUFFIX))
    }
}

/// A cache miss is explicit so callers can distinguish an ordinary cold miss
/// from an entry discarded because an object identity or page was stale.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CacheMissReason {
    Disabled,
    NotFound,
    StaleIdentity,
    Corrupt,
}

/// Result of looking up a page.
#[derive(Debug, Eq, PartialEq)]
pub enum CacheLookup {
    Hit(Vec<u8>),
    Miss(CacheMissReason),
}

impl CacheLookup {
    #[must_use]
    pub fn hit(&self) -> bool {
        matches!(self, Self::Hit(_))
    }

    #[must_use]
    pub fn into_option(self) -> Option<Vec<u8>> {
        match self {
            Self::Hit(bytes) => Some(bytes),
            Self::Miss(_) => None,
        }
    }
}

/// Configuration for a persistent page cache.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PersistentCacheOptions {
    pub enabled: bool,
    /// Maximum complete on-disk bytes, including entry headers.
    pub max_bytes: u64,
}

impl PersistentCacheOptions {
    #[must_use]
    pub const fn disabled() -> Self {
        Self {
            enabled: false,
            max_bytes: 0,
        }
    }
}

impl Default for PersistentCacheOptions {
    fn default() -> Self {
        Self {
            enabled: true,
            max_bytes: 4 * 1024 * 1024 * 1024,
        }
    }
}

/// Result of inserting a page.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CacheStoreOutcome {
    Stored { evicted: u64 },
    Disabled,
    Oversized,
}

/// Persistent cache counters. `bytes` is the complete on-disk size of valid
/// entries, including each fixed header. These are snapshots and do not
/// expose paths or object locators.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct PersistentCacheStats {
    pub entries: u64,
    pub bytes: u64,
    pub hits: u64,
    pub misses: u64,
    pub stale_identity_rejections: u64,
    pub corrupt_rejections: u64,
    pub writes: u64,
    pub evictions: u64,
    pub disabled_operations: u64,
}

#[derive(Clone, Debug)]
struct Entry {
    key: PageKey,
    file_name: String,
    sequence: u64,
    /// Complete on-disk entry size, including the fixed header.
    bytes: u64,
}

#[derive(Debug, Default)]
struct CacheState {
    entries: BTreeMap<String, Entry>,
    bytes: u64,
    next_sequence: u64,
    stats: PersistentCacheStats,
}

/// An optional persistent cache for exact immutable pages.
pub struct PersistentPageCache {
    cache_dir: PathBuf,
    options: PersistentCacheOptions,
    state: Mutex<CacheState>,
}

impl std::fmt::Debug for PersistentPageCache {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let state = self.state.lock().map_err(|_| std::fmt::Error)?;
        formatter
            .debug_struct("PersistentPageCache")
            .field("enabled", &self.options.enabled)
            .field("max_bytes", &self.options.max_bytes)
            .field("entries", &state.entries.len())
            .field("bytes", &state.bytes)
            .finish()
    }
}

impl PersistentPageCache {
    /// Open an enabled cache with a byte bound.  A zero bound is equivalent to
    /// disabled mode and does not create the directory.
    pub fn open(cache_dir: impl AsRef<Path>, max_bytes: u64) -> Result<Self, PersistentCacheError> {
        Self::open_with_options(
            cache_dir,
            PersistentCacheOptions {
                enabled: max_bytes != 0,
                max_bytes,
            },
        )
    }

    /// Open a cache with explicit enabled/disabled behavior.
    pub fn open_with_options(
        cache_dir: impl AsRef<Path>,
        options: PersistentCacheOptions,
    ) -> Result<Self, PersistentCacheError> {
        let cache_dir = cache_dir.as_ref().to_path_buf();
        let enabled = options.enabled && options.max_bytes != 0;
        let mut state = CacheState::default();
        if enabled {
            fs::create_dir_all(&cache_dir).map_err(|error| io_error("create directory", error))?;
            scan_entries(&cache_dir, &mut state)?;
        }
        let cache = Self {
            cache_dir,
            options: PersistentCacheOptions { enabled, ..options },
            state: Mutex::new(state),
        };
        if enabled {
            let mut state = cache.lock_state()?;
            cache.evict_locked(&mut state)?;
        }
        Ok(cache)
    }

    /// Construct a cache which never reads or writes the filesystem.
    pub fn disabled(cache_dir: impl AsRef<Path>) -> Self {
        Self {
            cache_dir: cache_dir.as_ref().to_path_buf(),
            options: PersistentCacheOptions::disabled(),
            state: Mutex::new(CacheState::default()),
        }
    }

    #[must_use]
    pub fn enabled(&self) -> bool {
        self.options.enabled
    }

    #[must_use]
    pub fn max_bytes(&self) -> u64 {
        self.options.max_bytes
    }

    #[must_use]
    pub fn cache_dir(&self) -> &Path {
        &self.cache_dir
    }

    /// Return the digest-based path used by a page key.  The path contains no
    /// locator, token, or caller-provided metadata.
    pub fn entry_path(&self, key: &PageKey) -> PathBuf {
        key.path(&self.cache_dir)
    }

    /// Look up a page and validate its complete header, range, length, and
    /// checksum before returning bytes.
    pub fn get(&self, key: &PageKey) -> Result<CacheLookup, PersistentCacheError> {
        validate_key(key)?;
        let mut state = self.lock_state()?;
        if !self.options.enabled {
            state.stats.disabled_operations = state.stats.disabled_operations.saturating_add(1);
            return Ok(CacheLookup::Miss(CacheMissReason::Disabled));
        }

        let file_name = file_name(key);
        if let Some(entry) = state.entries.get(&file_name).cloned() {
            match read_page_file(&self.cache_dir.join(&file_name), key) {
                Ok(bytes) => {
                    state.stats.hits = state.stats.hits.saturating_add(1);
                    return Ok(CacheLookup::Hit(bytes));
                }
                Err(PageReadError::Invalid) => {
                    self.remove_entry_locked(&mut state, &entry.file_name)?;
                    state.stats.corrupt_rejections =
                        state.stats.corrupt_rejections.saturating_add(1);
                    state.stats.misses = state.stats.misses.saturating_add(1);
                    return Ok(CacheLookup::Miss(CacheMissReason::Corrupt));
                }
                Err(PageReadError::Io(error)) => return Err(error),
            }
        }

        // A cache instance may have been reconstructed while another process
        // populated a page.  Check an exact digest path once before deciding
        // that it is a cold miss; this also keeps cross-process writes safe.
        let path = self.cache_dir.join(&file_name);
        if path.is_file() {
            match read_page_file(&path, key) {
                Ok(bytes) => {
                    let sequence = read_sequence(&path)?;
                    let payload_len = u64::try_from(bytes.len())
                        .map_err(|_| PersistentCacheError::AccountingOverflow)?;
                    let entry_len = entry_size(payload_len)?;
                    state.entries.insert(
                        file_name.clone(),
                        Entry {
                            key: *key,
                            file_name: file_name.clone(),
                            sequence,
                            bytes: entry_len,
                        },
                    );
                    state.bytes = state.bytes.saturating_add(entry_len);
                    state.next_sequence = state.next_sequence.max(sequence.saturating_add(1));
                    state.stats.entries = state.entries.len() as u64;
                    state.stats.bytes = state.bytes;
                    state.stats.hits = state.stats.hits.saturating_add(1);
                    return Ok(CacheLookup::Hit(bytes));
                }
                Err(PageReadError::Invalid) => {
                    remove_file_if_present(&path)?;
                    state.stats.corrupt_rejections =
                        state.stats.corrupt_rejections.saturating_add(1);
                    state.stats.misses = state.stats.misses.saturating_add(1);
                    return Ok(CacheLookup::Miss(CacheMissReason::Corrupt));
                }
                Err(PageReadError::Io(error)) => return Err(error),
            }
        }

        let stale_names = state
            .entries
            .values()
            .filter(|entry| entry.key.logical_eq(*key) && entry.key != *key)
            .map(|entry| entry.file_name.clone())
            .collect::<Vec<_>>();
        if !stale_names.is_empty() {
            for stale_name in stale_names {
                self.remove_entry_locked(&mut state, &stale_name)?;
            }
            state.stats.stale_identity_rejections =
                state.stats.stale_identity_rejections.saturating_add(1);
            state.stats.misses = state.stats.misses.saturating_add(1);
            return Ok(CacheLookup::Miss(CacheMissReason::StaleIdentity));
        }

        state.stats.misses = state.stats.misses.saturating_add(1);
        Ok(CacheLookup::Miss(CacheMissReason::NotFound))
    }

    /// Insert a complete page using a temporary file and atomic rename.
    pub fn put(
        &self,
        key: &PageKey,
        bytes: &[u8],
    ) -> Result<CacheStoreOutcome, PersistentCacheError> {
        validate_key(key)?;
        let expected_len = usize::try_from(key.range.length)
            .map_err(|_| PersistentCacheError::InvalidKey("range is too large"))?;
        if bytes.len() != expected_len {
            return Err(PersistentCacheError::InvalidKey(
                "page length does not match the exact range",
            ));
        }
        if sha256(bytes) != key.page_checksum {
            return Err(PersistentCacheError::PageChecksumMismatch);
        }

        let mut state = self.lock_state()?;
        if !self.options.enabled {
            state.stats.disabled_operations = state.stats.disabled_operations.saturating_add(1);
            return Ok(CacheStoreOutcome::Disabled);
        }
        let payload_len =
            u64::try_from(bytes.len()).map_err(|_| PersistentCacheError::AccountingOverflow)?;
        let entry_len = entry_size(payload_len)?;
        if entry_len > self.options.max_bytes {
            return Ok(CacheStoreOutcome::Oversized);
        }

        let sequence = state.next_sequence;
        state.next_sequence = state.next_sequence.saturating_add(1);
        let file_name = file_name(key);
        let destination = self.cache_dir.join(&file_name);
        write_atomic(
            &self.cache_dir,
            &destination,
            &encode_header(*key, sequence, payload_len),
            bytes,
        )?;

        if let Some(previous) = state.entries.remove(&file_name) {
            state.bytes = state.bytes.saturating_sub(previous.bytes);
        }
        state.bytes = state.bytes.saturating_add(entry_len);
        state.entries.insert(
            file_name.clone(),
            Entry {
                key: *key,
                file_name,
                sequence,
                bytes: entry_len,
            },
        );
        state.stats.writes = state.stats.writes.saturating_add(1);
        let evicted = self.evict_locked(&mut state)?;
        state.stats.entries = state.entries.len() as u64;
        state.stats.bytes = state.bytes;
        Ok(CacheStoreOutcome::Stored { evicted })
    }

    /// Remove all pages belonging to one object checksum.
    pub fn invalidate_object(
        &self,
        object_checksum: &[u8; 32],
    ) -> Result<u64, PersistentCacheError> {
        let mut state = self.lock_state()?;
        if !self.options.enabled {
            state.stats.disabled_operations = state.stats.disabled_operations.saturating_add(1);
            return Ok(0);
        }
        let names = state
            .entries
            .values()
            .filter(|entry| &entry.key.object_checksum == object_checksum)
            .map(|entry| entry.file_name.clone())
            .collect::<Vec<_>>();
        let count =
            u64::try_from(names.len()).map_err(|_| PersistentCacheError::AccountingOverflow)?;
        for name in names {
            self.remove_entry_locked(&mut state, &name)?;
        }
        Ok(count)
    }

    /// Current cache statistics.
    pub fn stats(&self) -> Result<PersistentCacheStats, PersistentCacheError> {
        let state = self.lock_state()?;
        Ok(PersistentCacheStats {
            entries: state.entries.len() as u64,
            bytes: state.bytes,
            ..state.stats
        })
    }

    /// Remove all valid entries.  Temporary files from an interrupted write
    /// are also ignored and left for the next explicit cache cleanup.
    pub fn clear(&self) -> Result<u64, PersistentCacheError> {
        let mut state = self.lock_state()?;
        let names = state.entries.keys().cloned().collect::<Vec<_>>();
        let count =
            u64::try_from(names.len()).map_err(|_| PersistentCacheError::AccountingOverflow)?;
        for name in names {
            self.remove_entry_locked(&mut state, &name)?;
        }
        Ok(count)
    }

    fn lock_state(&self) -> Result<std::sync::MutexGuard<'_, CacheState>, PersistentCacheError> {
        self.state.lock().map_err(|_| PersistentCacheError::Io {
            operation: "lock cache state",
            source: io::Error::other("cache state lock poisoned"),
        })
    }

    fn remove_entry_locked(
        &self,
        state: &mut CacheState,
        file_name: &str,
    ) -> Result<(), PersistentCacheError> {
        if let Some(entry) = state.entries.remove(file_name) {
            state.bytes = state.bytes.saturating_sub(entry.bytes);
            remove_file_if_present(&self.cache_dir.join(file_name))?;
        }
        Ok(())
    }

    fn evict_locked(&self, state: &mut CacheState) -> Result<u64, PersistentCacheError> {
        if !self.options.enabled {
            return Ok(0);
        }
        let mut evicted = 0u64;
        while state.bytes > self.options.max_bytes {
            let Some(name) = state
                .entries
                .values()
                .min_by_key(|entry| (entry.sequence, entry.file_name.clone()))
                .map(|entry| entry.file_name.clone())
            else {
                break;
            };
            self.remove_entry_locked(state, &name)?;
            evicted = evicted.saturating_add(1);
        }
        state.stats.evictions = state.stats.evictions.saturating_add(evicted);
        Ok(evicted)
    }
}

/// Compatibility aliases for callers that prefer the shorter page-cache name.
pub type PageCache = PersistentPageCache;
pub type PageCacheKey = PageKey;

fn validate_key(key: &PageKey) -> Result<(), PersistentCacheError> {
    key.range
        .offset
        .checked_add(key.range.length)
        .ok_or(PersistentCacheError::InvalidKey("range overflows"))?;
    Ok(())
}

fn file_name(key: &PageKey) -> String {
    format!("{}{}", hex(&key.digest()), ENTRY_SUFFIX)
}

fn hex(bytes: &[u8]) -> String {
    const DIGITS: &[u8; 16] = b"0123456789abcdef";
    let mut output = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        output.push(DIGITS[(byte >> 4) as usize] as char);
        output.push(DIGITS[(byte & 0x0f) as usize] as char);
    }
    output
}

fn sha256(bytes: &[u8]) -> [u8; 32] {
    Sha256::digest(bytes).into()
}

fn entry_size(payload_len: u64) -> Result<u64, PersistentCacheError> {
    u64::try_from(HEADER_BYTES)
        .map_err(|_| PersistentCacheError::AccountingOverflow)?
        .checked_add(payload_len)
        .ok_or(PersistentCacheError::AccountingOverflow)
}

fn encode_header(key: PageKey, sequence: u64, payload_len: u64) -> [u8; HEADER_BYTES] {
    let mut bytes = [0u8; HEADER_BYTES];
    bytes[0..8].copy_from_slice(MAGIC);
    bytes[8..10].copy_from_slice(&VERSION.to_le_bytes());
    bytes[10..12].copy_from_slice(&0u16.to_le_bytes());
    bytes[12..20].copy_from_slice(&sequence.to_le_bytes());
    bytes[20..52].copy_from_slice(&key.object_checksum);
    bytes[52..56].copy_from_slice(&key.section_id.to_le_bytes());
    bytes[56..64].copy_from_slice(&key.page_id.to_le_bytes());
    bytes[64..72].copy_from_slice(&key.range.offset.to_le_bytes());
    bytes[72..80].copy_from_slice(&key.range.length.to_le_bytes());
    bytes[80..112].copy_from_slice(&key.page_checksum);
    bytes[112..120].copy_from_slice(&payload_len.to_le_bytes());
    bytes
}

fn decode_header(bytes: &[u8; HEADER_BYTES]) -> Result<(PageKey, u64, u64), PersistentCacheError> {
    if &bytes[0..8] != MAGIC || u16::from_le_bytes(bytes[8..10].try_into().unwrap()) != VERSION {
        return Err(PersistentCacheError::InvalidEntry);
    }
    if u16::from_le_bytes(bytes[10..12].try_into().unwrap()) != 0 {
        return Err(PersistentCacheError::InvalidEntry);
    }
    let sequence = u64::from_le_bytes(bytes[12..20].try_into().unwrap());
    let mut object_checksum = [0u8; 32];
    object_checksum.copy_from_slice(&bytes[20..52]);
    let section_id = u32::from_le_bytes(bytes[52..56].try_into().unwrap());
    let page_id = u64::from_le_bytes(bytes[56..64].try_into().unwrap());
    let offset = u64::from_le_bytes(bytes[64..72].try_into().unwrap());
    let length = u64::from_le_bytes(bytes[72..80].try_into().unwrap());
    let mut page_checksum = [0u8; 32];
    page_checksum.copy_from_slice(&bytes[80..112]);
    let payload_len = u64::from_le_bytes(bytes[112..120].try_into().unwrap());
    let key = PageKey::new(
        object_checksum,
        section_id,
        page_id,
        ByteRange { offset, length },
        page_checksum,
    )?;
    if payload_len != length {
        return Err(PersistentCacheError::InvalidEntry);
    }
    Ok((key, sequence, payload_len))
}

enum PageReadError {
    Invalid,
    Io(PersistentCacheError),
}

fn read_page_file(path: &Path, expected: &PageKey) -> Result<Vec<u8>, PageReadError> {
    let mut file =
        File::open(path).map_err(|error| PageReadError::Io(io_error("open entry", error)))?;
    let total_len = file
        .metadata()
        .map_err(|error| PageReadError::Io(io_error("stat entry", error)))?
        .len();
    let mut header = [0u8; HEADER_BYTES];
    file.read_exact(&mut header).map_err(|error| {
        if error.kind() == io::ErrorKind::UnexpectedEof {
            PageReadError::Invalid
        } else {
            PageReadError::Io(io_error("read entry header", error))
        }
    })?;
    let (key, _sequence, payload_len) =
        decode_header(&header).map_err(|_| PageReadError::Invalid)?;
    let expected_total = entry_size(payload_len).map_err(|_| PageReadError::Invalid)?;
    if key != *expected || total_len != expected_total {
        return Err(PageReadError::Invalid);
    }
    let payload_size = usize::try_from(payload_len).map_err(|_| PageReadError::Invalid)?;
    let mut payload = vec![0u8; payload_size];
    file.read_exact(&mut payload).map_err(|error| {
        if error.kind() == io::ErrorKind::UnexpectedEof {
            PageReadError::Invalid
        } else {
            PageReadError::Io(io_error("read entry payload", error))
        }
    })?;
    if sha256(&payload) != expected.page_checksum {
        return Err(PageReadError::Invalid);
    }
    Ok(payload)
}

fn read_sequence(path: &Path) -> Result<u64, PersistentCacheError> {
    let mut file = File::open(path).map_err(|error| io_error("open entry", error))?;
    let mut header = [0u8; HEADER_BYTES];
    file.read_exact(&mut header)
        .map_err(|error| io_error("read entry header", error))?;
    decode_header(&header).map(|(_, sequence, _)| sequence)
}

fn scan_entries(cache_dir: &Path, state: &mut CacheState) -> Result<(), PersistentCacheError> {
    let directory = fs::read_dir(cache_dir).map_err(|error| io_error("read directory", error))?;
    for directory_entry in directory {
        let directory_entry =
            directory_entry.map_err(|error| io_error("read directory entry", error))?;
        let file_type = directory_entry
            .file_type()
            .map_err(|error| io_error("stat directory entry", error))?;
        if !file_type.is_file() {
            continue;
        }
        let file_name = directory_entry.file_name().to_string_lossy().into_owned();
        if file_name.starts_with(".cache-") && file_name.ends_with(".tmp") {
            remove_file_if_present(&directory_entry.path())?;
            continue;
        }
        if !file_name.ends_with(ENTRY_SUFFIX) {
            continue;
        }
        let path = directory_entry.path();
        let mut file = File::open(&path).map_err(|error| io_error("open cache entry", error))?;
        let total_len = file
            .metadata()
            .map_err(|error| io_error("stat cache entry", error))?
            .len();
        let mut header = [0u8; HEADER_BYTES];
        let valid = file.read_exact(&mut header).is_ok();
        let decoded = valid.then(|| decode_header(&header));
        drop(file);
        let Some(Ok((key, sequence, payload_len))) = decoded else {
            remove_file_if_present(&path)?;
            continue;
        };
        let Some(entry_len) = entry_size(payload_len).ok() else {
            remove_file_if_present(&path)?;
            continue;
        };
        if entry_len != total_len || file_name != file_name_for_key(&key) {
            remove_file_if_present(&path)?;
            continue;
        }
        state.bytes = state.bytes.saturating_add(entry_len);
        state.next_sequence = state.next_sequence.max(sequence.saturating_add(1));
        state.entries.insert(
            file_name.clone(),
            Entry {
                key,
                file_name,
                sequence,
                bytes: entry_len,
            },
        );
    }
    state.stats.entries = state.entries.len() as u64;
    state.stats.bytes = state.bytes;
    Ok(())
}

fn file_name_for_key(key: &PageKey) -> String {
    format!("{}{}", hex(&key.digest()), ENTRY_SUFFIX)
}

fn write_atomic(
    cache_dir: &Path,
    destination: &Path,
    header: &[u8; HEADER_BYTES],
    payload: &[u8],
) -> Result<(), PersistentCacheError> {
    let counter = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
    let temp_name = format!(".cache-{}-{counter:016x}.tmp", std::process::id());
    let temporary = cache_dir.join(temp_name);
    let result = (|| {
        let mut file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)
            .map_err(|error| io_error("create temporary entry", error))?;
        file.write_all(header)
            .map_err(|error| io_error("write entry header", error))?;
        file.write_all(payload)
            .map_err(|error| io_error("write entry payload", error))?;
        file.sync_all()
            .map_err(|error| io_error("sync temporary entry", error))?;
        fs::rename(&temporary, destination)
            .map_err(|error| io_error("publish cache entry", error))?;
        // A directory sync is supported by the Unix filesystems used by the
        // worker.  If a platform does not expose it, the renamed entry is
        // still complete and will be validated on the next open.
        if let Ok(directory) = File::open(cache_dir) {
            let _ = directory.sync_all();
        }
        Ok(())
    })();
    if result.is_err() {
        let _ = fs::remove_file(&temporary);
    }
    result
}

fn remove_file_if_present(path: &Path) -> Result<(), PersistentCacheError> {
    match fs::remove_file(path) {
        Ok(()) => Ok(()),
        Err(error) if error.kind() == io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(io_error("remove cache entry", error)),
    }
}
