//! A filter of the stored seed keys of a shard, so that most k-mers of an assembly are known to
//! be missing without reading a sketch bucket. It is built from the shard, or mapped from the key
//! file `NAME.jamkeys` that `jam index build --singleton` writes next to `NAME.jam`.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use crate::index::{IndexError, SEED_KEY_BITS, SeedFile};
use rayon::prelude::*;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU32, Ordering};

/// Bits of the filter per stored seed: about one false hit in 185 missing keys.
const BITS_PER_SEED: u64 = 12;

/// Buckets read at once while the filter is built.
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

/// The extension of the key file.
const KEY_FILE_EXTENSION: &str = "jamkeys";
const MAGIC: [u8; 8] = *b"JAMKEYS1";
const VERSION: u32 = 1;
/// Bytes of the key file header; the blocks follow, 32 bytes each, as little-endian words.
const HEADER_BYTES: usize = 64;
const BLOCK_BYTES: usize = 32;
/// Where the seeds file header of `NAME.jam` keeps the CRC-32 of its first 252 bytes, which covers
/// the sizes and CRC-32 of every section.
const INDEX_CHECK_AT: u64 = 252;

/// A split block Bloom filter: each key sets one bit in each of the eight words of one block.
/// It never misses a stored key.
pub struct KeyFilter {
    blocks: Blocks,
    count: usize,
}

enum Blocks {
    Built(Vec<[u32; 8]>),
    /// The key file, mapped read-only, so processes share its pages.
    Mapped(memmap2::Mmap),
}

/// What a key file must match: the shard it was built from.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Identity {
    bucket_bits: u32,
    seeds: u64,
    index_check: u32,
}

impl KeyFilter {
    /// Reads every sketch bucket of `seed_file` and adds the key of each stored seed.
    pub fn new(seed_file: &SeedFile) -> Result<Self, IndexError> {
        let count = (seed_file.seed_count() * BITS_PER_SEED)
            .div_ceil(256)
            .max(1) as usize;
        let blocks: Vec<[AtomicU32; 8]> = (0..count).map(|_| Default::default()).collect();
        let bucket_bits = seed_file.bucket_bits();
        let check_bits = SEED_KEY_BITS - bucket_bits;
        let buckets = 1u32 << bucket_bits;
        (0..buckets.div_ceil(CHUNK_BUCKETS))
            .into_par_iter()
            .try_for_each(|chunk| {
                let first = chunk * CHUNK_BUCKETS;
                let list: Vec<u32> = (first..buckets.min(first + CHUNK_BUCKETS)).collect();
                seed_file.read_sketch_buckets(&list, |bucket, part| {
                    for place in 0..part.len() {
                        let key = (u64::from(bucket) << check_bits) | part.check(place);
                        let (block, bits) = place_of(key, count);
                        for (word, bit) in blocks[block].iter().zip(bits) {
                            word.fetch_or(bit, Ordering::Relaxed);
                        }
                    }
                    Ok(())
                })
            })?;
        let blocks = (blocks.into_iter())
            .map(|block| block.map(AtomicU32::into_inner))
            .collect();
        Ok(Self {
            blocks: Blocks::Built(blocks),
            count,
        })
    }

    /// The key file of the index file `index` (`NAME.jam`): `NAME.jamkeys`.
    pub fn path_of(index: &Path) -> PathBuf {
        index.with_extension(KEY_FILE_EXTENSION)
    }

    /// Maps the key file of the local index file `index`, whose seeds are `seed_file`. `Ok(None)`
    /// when there is no key file; an error names why a key file does not fit the index.
    pub fn open(index: &Path, seed_file: &SeedFile) -> Result<Option<Self>, String> {
        let path = Self::path_of(index);
        let file = match std::fs::File::open(&path) {
            Ok(file) => file,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(None),
            Err(error) => return Err(format!("{}: {error}", path.display())),
        };
        let failed = |why: &str| format!("{} {why}", path.display());
        let wanted = identity(index, seed_file).map_err(|error| failed(&error.to_string()))?;
        // The file must not change while it is mapped; `jam index build` writes a new file and
        // renames it, so an old mapping keeps the old pages.
        let map = unsafe { memmap2::Mmap::map(&file) }.map_err(|e| failed(&e.to_string()))?;
        let (found, count) = read_header(&map).ok_or_else(|| failed("is not a key file"))?;
        let expected = (wanted.seeds * BITS_PER_SEED).div_ceil(256).max(1);
        if found != wanted || count as u64 != expected {
            return Err(failed("was built for another index"));
        }
        if map.len() != HEADER_BYTES + count * BLOCK_BYTES {
            return Err(failed("is cut short"));
        }
        Ok(Some(Self {
            blocks: Blocks::Mapped(map),
            count,
        }))
    }

    /// Builds the filter of the local index file `index` and writes its key file, through a
    /// temporary file that is renamed at the end. Returns the path of the key file.
    pub fn write(index: &Path, seed_file: &SeedFile) -> Result<PathBuf, IndexError> {
        let filter = Self::new(seed_file)?;
        let Blocks::Built(blocks) = &filter.blocks else {
            unreachable!("a new filter is built");
        };
        let identity = identity(index, seed_file)?;
        let path = Self::path_of(index);
        let partial = path.with_extension(format!("{KEY_FILE_EXTENSION}.partial"));
        let mut out = std::io::BufWriter::new(std::fs::File::create(&partial)?);
        out.write_all(&header(identity, filter.count))?;
        for block in blocks {
            for word in block {
                out.write_all(&word.to_le_bytes())?;
            }
        }
        out.into_inner()
            .map_err(|error| error.into_error())?
            .sync_all()?;
        std::fs::rename(&partial, &path)?;
        Ok(path)
    }

    /// Whether the shard may store `key`; `false` only for keys it does not store.
    pub fn may_hold(&self, key: u64) -> bool {
        let (block, bits) = place_of(key, self.count);
        match &self.blocks {
            Blocks::Built(blocks) => (blocks[block].iter().zip(bits)).all(|(w, bit)| w & bit != 0),
            Blocks::Mapped(map) => {
                let at = HEADER_BYTES + block * BLOCK_BYTES;
                let (words, _) = map[at..at + BLOCK_BYTES].as_chunks::<4>();
                (words.iter().zip(bits)).all(|(word, bit)| u32::from_le_bytes(*word) & bit != 0)
            }
        }
    }

    /// Whether the filter is mapped from a key file.
    pub fn is_mapped(&self) -> bool {
        matches!(self.blocks, Blocks::Mapped(_))
    }
}

/// The block of a key among `count` blocks and the bit to test in each of its words.
fn place_of(key: u64, count: usize) -> (usize, [u32; 8]) {
    let mixed = key.wrapping_mul(0x9e37_79b9_7f4a_7c15);
    let block = ((u128::from(mixed) * count as u128) >> 64) as usize;
    let low = (mixed as u32) ^ ((mixed >> 32) as u32);
    (
        block,
        SALTS.map(|salt| 1u32 << (low.wrapping_mul(salt) >> 27)),
    )
}

/// The identity of the local index file `index` with seeds `seed_file`.
fn identity(index: &Path, seed_file: &SeedFile) -> Result<Identity, IndexError> {
    let mut file = std::fs::File::open(index)?;
    let mut head = [0u8; INDEX_CHECK_AT as usize + 4];
    file.read_exact(&mut head)?;
    let at = INDEX_CHECK_AT as usize;
    Ok(Identity {
        bucket_bits: seed_file.bucket_bits(),
        seeds: seed_file.seed_count(),
        index_check: u32::from_le_bytes([head[at], head[at + 1], head[at + 2], head[at + 3]]),
    })
}

/// The key file header: magic, version, bucket bits, seeds, blocks, index check, bits per seed,
/// and the CRC-32 of the bytes before it.
fn header(identity: Identity, count: usize) -> [u8; HEADER_BYTES] {
    let mut head = [0u8; HEADER_BYTES];
    head[..8].copy_from_slice(&MAGIC);
    head[8..12].copy_from_slice(&VERSION.to_le_bytes());
    head[12..16].copy_from_slice(&identity.bucket_bits.to_le_bytes());
    head[16..24].copy_from_slice(&identity.seeds.to_le_bytes());
    head[24..32].copy_from_slice(&(count as u64).to_le_bytes());
    head[32..36].copy_from_slice(&identity.index_check.to_le_bytes());
    head[36..40].copy_from_slice(&(BITS_PER_SEED as u32).to_le_bytes());
    let crc = crc32fast::hash(&head[..60]);
    head[60..64].copy_from_slice(&crc.to_le_bytes());
    head
}

/// The identity and block count of a key file, if its header is whole and of this version.
fn read_header(bytes: &[u8]) -> Option<(Identity, usize)> {
    let head = bytes.get(..HEADER_BYTES)?;
    let u32_at = |at: usize| u32::from_le_bytes(head[at..at + 4].try_into().unwrap_or_default());
    let u64_at = |at: usize| u64::from_le_bytes(head[at..at + 8].try_into().unwrap_or_default());
    let valid = head[..8] == MAGIC
        && u32_at(8) == VERSION
        && u32_at(36) == BITS_PER_SEED as u32
        && u32_at(60) == crc32fast::hash(&head[..60]);
    let identity = Identity {
        bucket_bits: u32_at(12),
        seeds: u64_at(16),
        index_check: u32_at(32),
    };
    let count = usize::try_from(u64_at(24))
        .ok()
        .filter(|&count| count > 0)?;
    valid.then_some((identity, count))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::{random_bases, rng, test_shard};

    /// A key file answers as the filter built from its shard, and a key file of another shard,
    /// a damaged one or none is not used.
    #[test]
    fn key_file_matches_its_index() {
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
        assert!(KeyFilter::open(&index, seed_file).unwrap().is_none());
        let path = KeyFilter::write(&index, seed_file).unwrap();
        assert_eq!(path, one.path().join("shard.jamkeys"));
        let mapped = KeyFilter::open(&index, seed_file).unwrap().unwrap();
        assert!(mapped.is_mapped());
        let built = KeyFilter::new(seed_file).unwrap();
        for key in (0..20_000u64).map(|at| at.wrapping_mul(0x2545_f491_4f6c_dd1d) >> 20) {
            assert_eq!(mapped.may_hold(key), built.may_hold(key));
        }
        // The key file of another shard does not fit this one.
        let shard_two = test_shard(other.path(), &entries(&mut rng));
        let index_two = other.path().join("shard.jam");
        std::fs::copy(&path, KeyFilter::path_of(&index_two)).unwrap();
        assert!(KeyFilter::open(&index_two, shard_two.seed_file()).is_err());
        // A damaged header is refused.
        let mut bytes = std::fs::read(&path).unwrap();
        bytes[20] ^= 1;
        std::fs::write(&path, bytes).unwrap();
        assert!(KeyFilter::open(&index, seed_file).is_err());
    }
}
