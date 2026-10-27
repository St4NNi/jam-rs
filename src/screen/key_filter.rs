//! A filter of the stored seed keys of a shard, so that most k-mers of an assembly are known to
//! be missing without reading a sketch bucket.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use crate::index::{IndexError, SEED_KEY_BITS, SeedFile};
use rayon::prelude::*;
use std::sync::atomic::{AtomicU32, Ordering};

/// Bits of the filter per stored seed: about one false hit in 400 missing keys.
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

/// A split block Bloom filter: each key sets one bit in each of the eight words of one block.
/// It never misses a stored key.
pub(super) struct KeyFilter {
    blocks: Vec<[AtomicU32; 8]>,
}

impl KeyFilter {
    /// Reads every sketch bucket of `seed_file` and adds the key of each stored seed.
    pub fn new(seed_file: &SeedFile) -> Result<Self, IndexError> {
        let count = (seed_file.seed_count() * BITS_PER_SEED)
            .div_ceil(256)
            .max(1);
        let blocks = (0..count).map(|_| Default::default()).collect();
        let filter = Self { blocks };
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
                        filter.add((u64::from(bucket) << check_bits) | part.check(place));
                    }
                    Ok(())
                })
            })?;
        Ok(filter)
    }

    /// The block of a key and the bits to set in each of its words.
    fn place(&self, key: u64) -> (&[AtomicU32; 8], [u32; 8]) {
        let mixed = key.wrapping_mul(0x9e37_79b9_7f4a_7c15);
        let block = ((u128::from(mixed) * self.blocks.len() as u128) >> 64) as usize;
        let low = (mixed as u32) ^ ((mixed >> 32) as u32);
        let bits = SALTS.map(|salt| 1u32 << (low.wrapping_mul(salt) >> 27));
        (&self.blocks[block], bits)
    }

    fn add(&self, key: u64) {
        let (block, bits) = self.place(key);
        for (word, bit) in block.iter().zip(bits) {
            word.fetch_or(bit, Ordering::Relaxed);
        }
    }

    /// Whether the shard may store `key`; `false` only for keys it does not store.
    pub fn may_hold(&self, key: u64) -> bool {
        let (block, bits) = self.place(key);
        (block.iter().zip(bits)).all(|(word, bit)| word.load(Ordering::Relaxed) & bit != 0)
    }
}
