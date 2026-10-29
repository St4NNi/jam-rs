//! A filter of the stored seed keys of a shard, so that most k-mers of an assembly are known to
//! be missing without reading a sketch bucket. It is built from the shard, or mapped from the
//! filter file `NAME.jamfilter` that `jam index build --filter` writes next to `NAME.jam`.
//!
//! The filter is split into independent filters of about 64 MiB, one for each range of sketch
//! buckets, so they are built and checked in parallel and a key is routed by its bucket. In the
//! file, each filter starts on its own page.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use crate::index::{IndexError, SEED_KEY_BITS, SeedFile};
use rayon::prelude::*;
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU32, Ordering};

/// Bits of the filter per stored seed: about one false hit in 185 missing keys.
const BITS_PER_SEED: u64 = 12;

/// Bytes of one filter, about: the filter count is the next power of two that keeps them at most
/// this large.
const FILTER_BYTES: u64 = 64 << 20;

/// Filter files above this size draw a warning: filters are meant for indexes of small sequences.
/// The dense index of 853,462 plasmids (1.21 G seeds) needs 1.7 GiB.
pub const FILTER_WARNING_BYTES: u64 = 8 << 30;

/// Buckets read at once while a filter is built.
const CHUNK_BUCKETS: u32 = 4096;

/// Odd numbers that pick one bit in each word of a block from the same key.
const SALTS: [u32; 8] = [
    0x47b6_137b,
    0x4497_4d91,
    0x8824_ad5b,
    0xa2b7_289d,
    0x7054_95c7,
    0x2df1_424b,
    0x9efc_4947,
    0x5c6b_fb31,
];

const EXTENSION: &str = "jamfilter";
const MAGIC: [u8; 8] = *b"JAMFILT1";
const VERSION: u32 = 2;
/// Bytes of the file header; the table of filters follows (offset and CRC-32 of each, 16 bytes),
/// then the filters, each on its own page, as blocks of eight little-endian words.
const HEADER_BYTES: usize = 64;
const TABLE_ENTRY_BYTES: usize = 16;
const BLOCK_BYTES: usize = 32;
const PAGE: u64 = 4096;

/// Split block Bloom filters: each key sets one bit in each of the eight words of one block of
/// the filter of its bucket range. They never miss a stored key.
pub struct KeyFilter {
    blocks: Blocks,
    layout: Layout,
}

enum Blocks {
    /// The blocks of all filters, one filter after another.
    Built(Vec<[u32; 8]>),
    /// The filter file, mapped read-only, so processes share its pages.
    Mapped(memmap2::Mmap),
}

/// How the filters of a shard are cut and where they lie in the file.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Layout {
    filters: u64,
    /// Blocks of each filter.
    blocks: u64,
    /// Right shift of a seed key that gives its filter.
    shift: u32,
}

impl Layout {
    fn of(seeds: u64, bucket_bits: u32) -> Self {
        let total = (seeds * BITS_PER_SEED).div_ceil(256).max(1);
        let wanted = (total * BLOCK_BYTES as u64).div_ceil(FILTER_BYTES);
        let filters = wanted.next_power_of_two().min(1 << bucket_bits);
        Self {
            filters,
            blocks: total.div_ceil(filters),
            shift: SEED_KEY_BITS - filters.trailing_zeros(),
        }
    }

    /// Bytes from the start of one filter to the next in the file.
    fn stride(&self) -> u64 {
        (self.blocks * BLOCK_BYTES as u64).next_multiple_of(PAGE)
    }

    /// Where the first filter starts in the file.
    fn data_start(&self) -> u64 {
        (HEADER_BYTES as u64 + self.filters * TABLE_ENTRY_BYTES as u64).next_multiple_of(PAGE)
    }

    fn file_bytes(&self) -> u64 {
        self.data_start() + self.filters * self.stride()
    }
}

/// What a filter file must match: the shard it was built from.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Identity {
    bucket_bits: u32,
    seeds: u64,
    index_check: u32,
}

impl KeyFilter {
    /// Reads every sketch bucket of `seed_file` and adds the key of each stored seed.
    pub fn new(seed_file: &SeedFile) -> Result<Self, IndexError> {
        let layout = Layout::of(seed_file.seed_count(), seed_file.bucket_bits());
        let filters = (0..layout.filters)
            .into_par_iter()
            .map(|filter| build_filter(seed_file, &layout, filter))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            blocks: Blocks::Built(filters.concat()),
            layout,
        })
    }

    /// The filter file of the index file `index` (`NAME.jam`): `NAME.jamfilter`. The filter is bound
    /// to the header of the seeds file as it was read when the shard was opened, never to the
    /// file at `index` now.
    pub fn path_of(index: &Path) -> PathBuf {
        index.with_extension(EXTENSION)
    }

    /// The size of the filter file of an index with the seeds `seed_file`.
    pub fn file_bytes(seed_file: &SeedFile) -> u64 {
        Layout::of(seed_file.seed_count(), seed_file.bucket_bits()).file_bytes()
    }

    /// Builds the filters of the local index file `index` on all threads and writes the filter
    /// file, through a temporary file that is renamed at the end. Calls `warn` first when the
    /// file will be larger than `limit` bytes. Returns the path of the filter file.
    pub fn write(
        index: &Path,
        seed_file: &SeedFile,
        limit: u64,
        warn: impl FnOnce(String),
    ) -> Result<PathBuf, IndexError> {
        let layout = Layout::of(seed_file.seed_count(), seed_file.bucket_bits());
        if let Some(warning) = size_warning(layout.file_bytes(), limit) {
            warn(warning);
        }
        let identity = identity(seed_file);
        let path = Self::path_of(index);
        let partial = path.with_extension(format!("{EXTENSION}.partial"));
        let file = std::fs::File::create(&partial)?;
        file.set_len(layout.file_bytes())?;
        let crcs = (0..layout.filters)
            .into_par_iter()
            .map(|filter| {
                let blocks = build_filter(seed_file, &layout, filter)?;
                let bytes: Vec<u8> = (blocks.iter().flatten())
                    .flat_map(|word| word.to_le_bytes())
                    .collect();
                let offset = layout.data_start() + filter * layout.stride();
                file.write_all_at(&bytes, offset)?;
                Ok(crc32fast::hash(&bytes))
            })
            .collect::<Result<Vec<u32>, IndexError>>()?;
        file.write_all_at(&header(identity, &layout, &crcs), 0)?;
        file.sync_all()?;
        std::fs::rename(&partial, &path)?;
        Ok(path)
    }

    /// Maps the filter file of the local index file `index`, whose seeds are `seed_file`, and
    /// checks every filter on all threads. Calls `warn` when the file is larger than `limit`
    /// bytes. `Ok(None)` when there is no filter file; an error names why it does not fit.
    pub fn open(
        index: &Path,
        seed_file: &SeedFile,
        limit: u64,
        warn: impl FnOnce(String),
    ) -> Result<Option<Self>, String> {
        let path = Self::path_of(index);
        let file = match std::fs::File::open(&path) {
            Ok(file) => file,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(None),
            Err(error) => return Err(format!("{}: {error}", path.display())),
        };
        let failed = |why: &str| format!("{} {why}", path.display());
        let wanted = identity(seed_file);
        // The file must not change while it is mapped; it is written under another name and
        // renamed, so an old mapping keeps the old pages.
        let map = unsafe { memmap2::Mmap::map(&file) }.map_err(|e| failed(&e.to_string()))?;
        let (found, layout, crcs) = read_header(&map).ok_or_else(|| failed("is damaged"))?;
        if found != wanted || layout != Layout::of(wanted.seeds, wanted.bucket_bits) {
            return Err(failed("was built for another index"));
        }
        if map.len() as u64 != layout.file_bytes() {
            return Err(failed("is cut short"));
        }
        if let Some(warning) = size_warning(layout.file_bytes(), limit) {
            warn(warning);
        }
        let filter_bytes = (layout.blocks as usize) * BLOCK_BYTES;
        let intact = crcs.par_iter().enumerate().all(|(filter, &crc)| {
            let start = (layout.data_start() + filter as u64 * layout.stride()) as usize;
            crc32fast::hash(&map[start..start + filter_bytes]) == crc
        });
        if !intact {
            return Err(failed("is damaged"));
        }
        Ok(Some(Self {
            blocks: Blocks::Mapped(map),
            layout,
        }))
    }

    /// Whether the shard may store `key`; `false` only for keys it does not store.
    pub fn may_hold(&self, key: u64) -> bool {
        let layout = &self.layout;
        let filter = key >> layout.shift;
        let (block, bits) = place_of(key, layout.blocks);
        match &self.blocks {
            Blocks::Built(blocks) => {
                let words = &blocks[(filter * layout.blocks + block) as usize];
                (words.iter().zip(bits)).all(|(word, bit)| word & bit != 0)
            }
            Blocks::Mapped(map) => {
                let start = layout.data_start() + filter * layout.stride();
                let at = (start + block * BLOCK_BYTES as u64) as usize;
                let (words, _) = map[at..at + BLOCK_BYTES].as_chunks::<4>();
                (words.iter().zip(bits)).all(|(word, bit)| u32::from_le_bytes(*word) & bit != 0)
            }
        }
    }

    /// Whether the filter is mapped from a filter file.
    pub fn is_mapped(&self) -> bool {
        matches!(self.blocks, Blocks::Mapped(_))
    }
}

/// The warning for a filter file of `bytes` bytes above `limit`.
fn size_warning(bytes: u64, limit: u64) -> Option<String> {
    let gib = bytes as f64 / f64::from(1u32 << 30);
    (bytes > limit).then(|| {
        format!(
            "the filter will be about {gib:.1} GiB, filters are meant for indexes of small \
             sequences, not for assemblies"
        )
    })
}

/// The blocks of filter `filter`: the keys of the stored seeds of its buckets.
fn build_filter(
    seed_file: &SeedFile,
    layout: &Layout,
    filter: u64,
) -> Result<Vec<[u32; 8]>, IndexError> {
    let blocks: Vec<[AtomicU32; 8]> = (0..layout.blocks).map(|_| Default::default()).collect();
    let bucket_bits = seed_file.bucket_bits();
    let check_bits = SEED_KEY_BITS - bucket_bits;
    let per_filter = ((1u64 << bucket_bits) / layout.filters) as u32;
    let first = filter as u32 * per_filter;
    (0..per_filter.div_ceil(CHUNK_BUCKETS))
        .into_par_iter()
        .try_for_each(|chunk| {
            let start = first + chunk * CHUNK_BUCKETS;
            let end = (first + per_filter).min(start + CHUNK_BUCKETS);
            let list: Vec<u32> = (start..end).collect();
            seed_file.read_sketch_buckets(&list, |bucket, part| {
                for place in 0..part.len() {
                    let key = (u64::from(bucket) << check_bits) | part.check(place);
                    let (block, bits) = place_of(key, layout.blocks);
                    for (word, bit) in blocks[block as usize].iter().zip(bits) {
                        word.fetch_or(bit, Ordering::Relaxed);
                    }
                }
                Ok(())
            })
        })?;
    Ok((blocks.into_iter())
        .map(|block| block.map(AtomicU32::into_inner))
        .collect())
}

/// The block of a key among `blocks` blocks and the bit to test in each of its words.
fn place_of(key: u64, blocks: u64) -> (u64, [u32; 8]) {
    let mixed = key.wrapping_mul(0x9e37_79b9_7f4a_7c15);
    let block = ((u128::from(mixed) * u128::from(blocks)) >> 64) as u64;
    let low = (mixed as u32) ^ ((mixed >> 32) as u32);
    (
        block,
        SALTS.map(|salt| 1u32 << (low.wrapping_mul(salt) >> 27)),
    )
}

/// The identity of the shard with the seeds `seed_file`, from the header it read when it was
/// opened, so a file replaced since then cannot lend its identity.
fn identity(seed_file: &SeedFile) -> Identity {
    Identity {
        bucket_bits: seed_file.bucket_bits(),
        seeds: seed_file.seed_count(),
        index_check: crc32fast::hash(seed_file.header()),
    }
}

/// The file header and the table of filters: magic, version, bucket bits, seeds, filters, bits
/// per seed, blocks per filter, the CRC-32 of the seeds file header, the CRC-32 of the table and of the header before it;
/// then the offset and CRC-32 of each filter.
fn header(identity: Identity, layout: &Layout, crcs: &[u32]) -> Vec<u8> {
    let mut table = Vec::with_capacity(crcs.len() * TABLE_ENTRY_BYTES);
    for (filter, crc) in crcs.iter().enumerate() {
        let offset = layout.data_start() + filter as u64 * layout.stride();
        table.extend_from_slice(&offset.to_le_bytes());
        table.extend_from_slice(&crc.to_le_bytes());
        table.extend_from_slice(&[0; 4]);
    }
    let mut head = vec![0u8; HEADER_BYTES];
    head[..8].copy_from_slice(&MAGIC);
    head[8..12].copy_from_slice(&VERSION.to_le_bytes());
    head[12..16].copy_from_slice(&identity.bucket_bits.to_le_bytes());
    head[16..24].copy_from_slice(&identity.seeds.to_le_bytes());
    head[24..32].copy_from_slice(&layout.filters.to_le_bytes());
    head[32..36].copy_from_slice(&(BITS_PER_SEED as u32).to_le_bytes());
    head[36..44].copy_from_slice(&layout.blocks.to_le_bytes());
    head[44..48].copy_from_slice(&identity.index_check.to_le_bytes());
    head[48..52].copy_from_slice(&crc32fast::hash(&table).to_le_bytes());
    let crc = crc32fast::hash(&head[..60]);
    head[60..64].copy_from_slice(&crc.to_le_bytes());
    head.extend(table);
    head
}

/// The identity, layout and filter CRC-32 of a filter file, if its header and table are whole,
/// of this version, and put each filter where the layout says.
fn read_header(bytes: &[u8]) -> Option<(Identity, Layout, Vec<u32>)> {
    let head = bytes.get(..HEADER_BYTES)?;
    let u32_at = |bytes: &[u8], at: usize| -> Option<u32> {
        Some(u32::from_le_bytes(bytes.get(at..at + 4)?.try_into().ok()?))
    };
    let u64_at = |bytes: &[u8], at: usize| -> Option<u64> {
        Some(u64::from_le_bytes(bytes.get(at..at + 8)?.try_into().ok()?))
    };
    let valid = head[..8] == MAGIC
        && u32_at(head, 8)? == VERSION
        && u32_at(head, 32)? == BITS_PER_SEED as u32
        && u32_at(head, 60)? == crc32fast::hash(&head[..60]);
    let identity = Identity {
        bucket_bits: u32_at(head, 12)?,
        seeds: u64_at(head, 16)?,
        index_check: u32_at(head, 44)?,
    };
    let filters = u64_at(head, 24)?;
    if !valid || filters == 0 || identity.bucket_bits > SEED_KEY_BITS {
        return None;
    }
    let layout = Layout::of(identity.seeds, identity.bucket_bits);
    let table_bytes = usize::try_from(filters).ok()? * TABLE_ENTRY_BYTES;
    let table = bytes.get(HEADER_BYTES..HEADER_BYTES + table_bytes)?;
    if filters != layout.filters
        || u64_at(head, 36)? != layout.blocks
        || u32_at(head, 48)? != crc32fast::hash(table)
    {
        return None;
    }
    let mut crcs = Vec::with_capacity(table_bytes / TABLE_ENTRY_BYTES);
    for (filter, entry) in table.chunks(TABLE_ENTRY_BYTES).enumerate() {
        if u64_at(entry, 0)? != layout.data_start() + filter as u64 * layout.stride() {
            return None;
        }
        crcs.push(u32_at(entry, 8)?);
    }
    Some((identity, layout, crcs))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::{random_bases, rng, test_shard};

    /// A filter file answers as the filter built from its shard, and a filter file of another
    /// shard or a damaged one is not used.
    #[test]
    fn filter_file_matches_its_index() {
        let mut rng = rng(17);
        let (one, other) = (tempfile::tempdir().unwrap(), tempfile::tempdir().unwrap());
        let entries = |rng: &mut _| {
            vec![
                vec![random_bases(rng, 20_000)],
                vec![random_bases(rng, 9_000)],
            ]
        };
        let shard = test_shard(one.path(), &entries(&mut rng));
        let index = one.path().join("shard.jam");
        let seed_file = shard.seed_file();
        let open = |index: &Path, seed_file| KeyFilter::open(index, seed_file, u64::MAX, drop);
        assert!(open(&index, seed_file).unwrap().is_none());
        let path = KeyFilter::write(&index, seed_file, u64::MAX, drop).unwrap();
        assert_eq!(path, one.path().join("shard.jamfilter"));
        let mapped = open(&index, seed_file).unwrap().unwrap();
        assert!(mapped.is_mapped());
        let built = KeyFilter::new(seed_file).unwrap();
        for key in (0..20_000u64).map(|at| at.wrapping_mul(0x2545_f491_4f6c_dd1d) >> 20) {
            assert_eq!(mapped.may_hold(key), built.may_hold(key));
        }
        // The filter file of another shard does not fit this one.
        let shard_two = test_shard(other.path(), &entries(&mut rng));
        let index_two = other.path().join("shard.jam");
        std::fs::copy(&path, KeyFilter::path_of(&index_two)).unwrap();
        assert!(open(&index_two, shard_two.seed_file()).is_err());
        // A changed bit in the header or in a filter is found.
        let length = std::fs::metadata(&path).unwrap().len() as usize;
        let layout = Layout::of(seed_file.seed_count(), seed_file.bucket_bits());
        for at in [20, layout.data_start() as usize + 5] {
            assert!(at < length);
            let mut bytes = std::fs::read(&path).unwrap();
            bytes[at] ^= 1;
            std::fs::write(&path, &bytes).unwrap();
            assert!(open(&index, seed_file).is_err(), "{at}");
            bytes[at] ^= 1;
            std::fs::write(&path, &bytes).unwrap();
        }
        assert!(open(&index, seed_file).unwrap().is_some());
        // An index file replaced at the path does not change which filter fits the open shard:
        // the filter of shard one still fits it, and does not fit the shard now at the path.
        let replacement = one.path().join("replacement.jam");
        std::fs::copy(&index_two, &replacement).unwrap();
        std::fs::rename(&replacement, &index).unwrap();
        assert!(open(&index, seed_file).unwrap().is_some());
        let now = crate::index::Shard::open(&index).unwrap();
        assert!(open(&index, now.seed_file()).is_err());
    }

    /// Writing and mapping a filter file above the limit warn, and below it do not.
    #[test]
    fn large_filter_warns() {
        let directory = tempfile::tempdir().unwrap();
        let shard = test_shard(directory.path(), &[vec![random_bases(&mut rng(3), 5_000)]]);
        let (index, seed_file) = (directory.path().join("shard.jam"), shard.seed_file());
        let bytes = KeyFilter::file_bytes(seed_file);
        for (limit, warns) in [(bytes, false), (bytes - 1, true)] {
            let mut warned = None;
            KeyFilter::write(&index, seed_file, limit, |text| warned = Some(text)).unwrap();
            assert_eq!(warned.is_some(), warns);
            let mut warned = None;
            let opened = KeyFilter::open(&index, seed_file, limit, |text| warned = Some(text));
            assert!(opened.unwrap().is_some());
            assert_eq!(warned.is_some(), warns);
            if let Some(text) = warned {
                assert!(text.contains("GiB") && text.contains("not for assemblies"));
            }
        }
    }
}
