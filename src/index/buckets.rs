//! The two bucket encodings of the seeds file: a sketch bucket holds the check values and color set
//! ids of its stored seeds; a locator bucket holds, in the same order, their seed locations.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::IndexError;
use super::bits::{BitWriter, bits_needed, read_bits, u32_at};

const SKETCH_HEADER_BYTES: usize = 8;
const LOCATOR_HEADER_BYTES: usize = 12;
const MAX_ESCAPES: usize = (1 << 24) - 1;
/// Bits that one read takes when flags are counted.
const SCAN_BITS: u32 = 56;

/// Appends the sketch bucket of `seeds`, as (check value, color set id) with check values of
/// `check_bits` bits, ascending and distinct. An empty bucket has no bytes. Check values are
/// stored as Elias-Fano lists: the high part of each in unary, the low part at a fixed width.
pub(crate) fn encode_sketch_bucket(seeds: &[(u64, u32)], check_bits: u32, out: &mut Vec<u8>) {
    if seeds.is_empty() {
        return;
    }
    let ids: Vec<u32> = seeds.iter().map(|(_, id)| *id).collect();
    let (width, escapes) = choose_width(&ids);
    out.extend_from_slice(&(seeds.len() as u32).to_le_bytes());
    out.extend_from_slice(&((width << 24) | escapes as u32).to_le_bytes());
    let direct = (1u64 << width) - escapes as u64;
    let mut escaped = 0u64;
    let mut values = Vec::with_capacity(ids.len());
    for &id in &ids {
        if u64::from(id) < direct {
            values.push(u64::from(id));
        } else {
            values.push(direct + escaped);
            escaped += 1;
            out.extend_from_slice(&id.to_le_bytes());
        }
    }
    let (high_bits, low_bits) = check_split(seeds.len() as u64, check_bits);
    let mut writer = BitWriter::new(out);
    let mut group = 0u64;
    for (check, _) in seeds {
        while group < check >> low_bits {
            writer.push(0, 1);
            group += 1;
        }
        writer.push(1, 1);
    }
    for _ in group..1 << high_bits {
        writer.push(0, 1);
    }
    for (check, _) in seeds {
        writer.push(*check, low_bits);
    }
    for value in values {
        writer.push(value, width);
    }
    writer.finish();
}

/// A decoded view of one sketch bucket.
#[derive(Clone, Copy, Debug)]
pub struct SketchBucket<'a> {
    bytes: &'a [u8],
    seeds: u32,
    width: u32,
    escapes: u32,
    /// Bits of the low part of a check value.
    low_bits: u32,
    /// Bits of the unary high parts: one per seed and one per possible high part.
    upper_bits: u64,
}

impl<'a> SketchBucket<'a> {
    /// Checks the header and size of a bucket; empty bytes are an empty bucket.
    pub(crate) fn decode(bytes: &'a [u8], check_bits: u32) -> Result<Self, IndexError> {
        let damaged = || IndexError::format("seeds sketch bucket");
        let seeds = if bytes.is_empty() {
            0
        } else {
            u32_at(bytes, 0).ok_or_else(damaged)?
        };
        let packed = u32_at(bytes, 4).unwrap_or_default();
        let (width, escapes) = (packed >> 24, packed & 0x00ff_ffff);
        let (high_bits, low_bits) = check_split(u64::from(seeds), check_bits);
        let upper_bits = u64::from(seeds) + (1 << high_bits);
        let bucket = Self {
            bytes,
            seeds,
            width,
            escapes,
            low_bits,
            upper_bits,
        };
        if bytes.is_empty() {
            return Ok(bucket);
        }
        let bits = upper_bits + u64::from(seeds) * u64::from(low_bits + width);
        let expected = (SKETCH_HEADER_BYTES as u64 + 4 * u64::from(escapes)) + bits.div_ceil(8);
        if seeds == 0
            || width > 32
            || u64::from(escapes) > 1u64 << width
            || expected != bytes.len() as u64
        {
            return Err(damaged());
        }
        Ok(bucket)
    }

    /// Stored seeds in the bucket.
    pub fn len(&self) -> u32 {
        self.seeds
    }

    /// Whether the bucket holds no seed.
    #[cfg(test)]
    pub fn is_empty(&self) -> bool {
        self.seeds == 0
    }

    /// The place of the seed with this check value among the seeds of the bucket, and its color
    /// set id.
    pub fn find(&self, check: u64) -> Option<(u32, u32)> {
        let bits = self.bits();
        let group = check >> self.low_bits;
        // The seeds of a high part follow the zero that ends the group before it.
        let mut position = match group.checked_sub(1) {
            None => 0,
            Some(before) => select(bits, 0, self.upper_bits, before, false)? + 1,
        };
        let mut place = position.checked_sub(group)?;
        let wanted = check & ((1 << self.low_bits) - 1);
        while position < self.upper_bits && read_bits(bits, position, 1) == 1 {
            let low = self.low_part(place);
            if low >= wanted {
                let place = place as u32;
                return (low == wanted).then(|| (place, self.color_set(place)));
            }
            (place, position) = (place + 1, position + 1);
        }
        None
    }

    /// The check value of the seed at `place`.
    pub fn check(&self, place: u32) -> u64 {
        let place = u64::from(place);
        let position = select(self.bits(), 0, self.upper_bits, place, true).unwrap_or(place);
        ((position - place) << self.low_bits) | self.low_part(place)
    }

    fn low_part(&self, place: u64) -> u64 {
        let bit = self.upper_bits + place * u64::from(self.low_bits);
        read_bits(self.bits(), bit, self.low_bits)
    }

    /// The color set id of the seed at `place`.
    pub fn color_set(&self, place: u32) -> u32 {
        let start = self.upper_bits + u64::from(self.seeds) * u64::from(self.low_bits);
        let value = read_bits(
            self.bits(),
            start + u64::from(place) * u64::from(self.width),
            self.width,
        );
        let direct = (1u64 << self.width) - u64::from(self.escapes);
        if value < direct {
            return value as u32;
        }
        // The size check of `decode` keeps every escaped id inside the bucket.
        let at = SKETCH_HEADER_BYTES + 4 * (value - direct) as usize;
        u32_at(self.bytes, at).unwrap_or_default()
    }

    fn bits(&self) -> &'a [u8] {
        let start = SKETCH_HEADER_BYTES + 4 * self.escapes as usize;
        self.bytes.get(start..).unwrap_or_default()
    }
}

/// The high and low bits of the check values of a bucket of `seeds` seeds: about log2(seeds)
/// high bits, so the unary high parts cost at most three bits per seed.
fn check_split(seeds: u64, check_bits: u32) -> (u32, u32) {
    let high_bits = bits_needed(seeds).min(check_bits);
    (high_bits, check_bits - high_bits)
}

/// The id width and escape count with the fewest bits for these ids.
fn choose_width(ids: &[u32]) -> (u32, usize) {
    let mut sorted = ids.to_vec();
    sorted.sort_unstable();
    let count = sorted.len();
    let max = u64::from(sorted[count - 1]);
    let mut best = (bits_needed(max), 0);
    let mut best_cost = count as u64 * u64::from(best.0);
    for width in 0..best.0 {
        let slots = 1u64 << width;
        for escapes in 0..=count.min(MAX_ESCAPES) {
            let cost = count as u64 * u64::from(width) + 32 * escapes as u64;
            if cost >= best_cost || escapes as u64 > slots {
                break;
            }
            let direct = slots - escapes as u64;
            if escapes == count || u64::from(sorted[count - 1 - escapes]) < direct {
                (best, best_cost) = ((width, escapes), cost);
                break;
            }
        }
    }
    best
}

/// Appends the locator bucket of the seeds of one sketch bucket, in its order. `seed_pairs`
/// holds the entries of each seed's color set, `pair_locations` the seed locations of each such
/// pair of seed and entry (at least one), and `locations` the values of all, in order.
pub(crate) fn encode_locator_bucket(
    seed_pairs: &[u32],
    pair_locations: &[u32],
    locations: &[u32],
    out: &mut Vec<u8>,
) {
    if seed_pairs.is_empty() {
        return;
    }
    let location_bits = bits_needed(u64::from(locations.iter().copied().max().unwrap_or(0)));
    out.extend_from_slice(&(pair_locations.len() as u32).to_le_bytes());
    out.extend_from_slice(&(locations.len() as u32).to_le_bytes());
    out.extend_from_slice(&[location_bits as u8, 0, 0, 0]);
    let mut writer = BitWriter::new(out);
    for counts in [seed_pairs, pair_locations] {
        for &count in counts {
            for flag in 1..=count {
                writer.push(u64::from(flag == count), 1);
            }
        }
    }
    for &location in locations {
        writer.push(u64::from(location), location_bits);
    }
    writer.finish();
}

/// A decoded view of one locator bucket.
#[derive(Clone, Copy, Debug)]
pub struct LocatorBucket<'a> {
    bits: &'a [u8],
    pairs: u64,
    locations: u64,
    location_bits: u32,
}

impl<'a> LocatorBucket<'a> {
    /// Checks the header and size of a bucket; empty bytes are an empty bucket.
    pub(crate) fn decode(bytes: &'a [u8]) -> Result<Self, IndexError> {
        if bytes.is_empty() {
            return Ok(Self {
                bits: bytes,
                pairs: 0,
                locations: 0,
                location_bits: 0,
            });
        }
        let damaged = || IndexError::format("seeds locator bucket");
        let pairs = u64::from(u32_at(bytes, 0).ok_or_else(damaged)?);
        let locations = u64::from(u32_at(bytes, 4).ok_or_else(damaged)?);
        let location_bits = u32::from(*bytes.get(8).ok_or_else(damaged)?);
        let bits = pairs + locations * (1 + u64::from(location_bits));
        let expected = LOCATOR_HEADER_BYTES as u64 + bits.div_ceil(8);
        if location_bits > 32 || pairs > locations || expected != bytes.len() as u64 {
            return Err(damaged());
        }
        Ok(Self {
            bits: &bytes[LOCATOR_HEADER_BYTES..],
            pairs,
            locations,
            location_bits,
        })
    }

    /// Appends the seed location values of the seed at `place` in the entry at `rank` of its
    /// color set, or fails when the bucket holds no such pair.
    pub fn locations(&self, place: u32, rank: u32, out: &mut Vec<u32>) -> Result<(), IndexError> {
        let missing = || IndexError::format("seeds locator bucket misses a seed location");
        let first_pair = match place.checked_sub(1) {
            None => 0,
            Some(before) => {
                select(self.bits, 0, self.pairs, u64::from(before), true).ok_or_else(missing)? + 1
            }
        };
        let pair = first_pair + u64::from(rank);
        let seed_end =
            select(self.bits, 0, self.pairs, u64::from(place), true).ok_or_else(missing)?;
        if pair > seed_end {
            return Err(missing());
        }
        let flags = self.pairs;
        let first = match pair.checked_sub(1) {
            None => 0,
            Some(before) => {
                select(self.bits, flags, self.locations, before, true).ok_or_else(missing)? + 1
            }
        };
        let last = select(self.bits, flags, self.locations, pair, true).ok_or_else(missing)?;
        let values = flags + self.locations;
        for index in first..=last {
            let bit = values + index * u64::from(self.location_bits);
            out.push(read_bits(self.bits, bit, self.location_bits) as u32);
        }
        Ok(())
    }

    /// Calls `visit(place, rank, location value)` for every seed location in stored order, and
    /// stops at the first error.
    pub fn for_each(
        &self,
        mut visit: impl FnMut(u32, u32, u32) -> Result<(), IndexError>,
    ) -> Result<(), IndexError> {
        let (mut place, mut rank) = (0u32, 0u32);
        let mut pair = 0u64;
        let values = self.pairs + self.locations;
        let width = self.location_bits;
        for index in 0..self.locations {
            let value = read_bits(self.bits, values + index * u64::from(width), width);
            visit(place, rank, value as u32)?;
            if read_bits(self.bits, self.pairs + index, 1) == 1 {
                if read_bits(self.bits, pair, 1) == 1 {
                    place += 1;
                    rank = 0;
                } else {
                    rank += 1;
                }
                pair += 1;
            }
        }
        Ok(())
    }
}

/// The offset (from `start`) of the bit with index `index` among the bits equal to `one` in the
/// `count` bits that start at bit `start`.
fn select(bits: &[u8], start: u64, count: u64, index: u64, one: bool) -> Option<u64> {
    let mut left = index;
    let mut offset = 0u64;
    while offset < count {
        let width = SCAN_BITS.min((count - offset) as u32);
        let mut word = read_bits(bits, start + offset, width);
        if !one {
            word = !word & ((1 << width) - 1);
        }
        let ones = u64::from(word.count_ones());
        if left < ones {
            for _ in 0..left {
                word &= word - 1;
            }
            return Some(offset + u64::from(word.trailing_zeros()));
        }
        left -= ones;
        offset += u64::from(width);
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::rng;
    use rand::RngExt;

    #[test]
    fn sketch_buckets_decode() {
        let mut draws = rng(5);
        for (count, check_bits) in [(1usize, 20u32), (300, 23), (2_000, 36), (50, 16)] {
            let mut checks: Vec<u64> = (0..count)
                .map(|_| draws.random::<u64>() & ((1 << check_bits) - 1))
                .collect();
            checks.sort_unstable();
            checks.dedup();
            let seeds: Vec<(u64, u32)> = checks
                .iter()
                .map(|&check| {
                    let id = match draws.random_range(0..20) {
                        0 => draws.random_range(1_000_000..u32::MAX),
                        _ => draws.random_range(0..30),
                    };
                    (check, id)
                })
                .collect();
            let mut bytes = vec![7];
            encode_sketch_bucket(&seeds, check_bits, &mut bytes);
            let bucket = SketchBucket::decode(&bytes[1..], check_bits).unwrap();
            assert_eq!(bucket.len() as usize, seeds.len());
            for (place, &(check, id)) in seeds.iter().enumerate() {
                assert_eq!(bucket.find(check), Some((place as u32, id)));
                assert_eq!(bucket.check(place as u32), check);
                if !checks.contains(&(check + 1)) {
                    assert_eq!(bucket.find(check + 1), None);
                }
            }
            assert!(SketchBucket::decode(&bytes[2..], check_bits).is_err());
        }
        let empty = SketchBucket::decode(&[], 20).unwrap();
        assert!(empty.is_empty());
        assert_eq!(empty.find(0), None);
    }

    #[test]
    fn locator_buckets_decode() {
        let mut draws = rng(6);
        let seed_pairs: Vec<u32> = (0..200).map(|_| draws.random_range(1..80)).collect();
        let pairs: u32 = seed_pairs.iter().sum();
        let pair_locations: Vec<u32> = (0..pairs)
            .map(|_| match draws.random_range(0..10) {
                0 => draws.random_range(2..70),
                _ => 1,
            })
            .collect();
        let total: u32 = pair_locations.iter().sum();
        let locations: Vec<u32> = (0..total).map(|_| draws.random_range(0..1 << 21)).collect();
        let mut bytes = vec![9; 3];
        encode_locator_bucket(&seed_pairs, &pair_locations, &locations, &mut bytes);
        let bucket = LocatorBucket::decode(&bytes[3..]).unwrap();
        let mut visited = Vec::new();
        bucket
            .for_each(|place, rank, value| {
                visited.push((place, rank, value));
                Ok(())
            })
            .unwrap();
        let mut expected = Vec::new();
        let (mut pair, mut at) = (0usize, 0usize);
        for (place, &count) in seed_pairs.iter().enumerate() {
            for rank in 0..count {
                let length = pair_locations[pair] as usize;
                let mut out = Vec::new();
                bucket.locations(place as u32, rank, &mut out).unwrap();
                assert_eq!(out, locations[at..at + length]);
                expected.extend(out.iter().map(|&value| (place as u32, rank, value)));
                (pair, at) = (pair + 1, at + length);
            }
            assert!(
                bucket
                    .locations(place as u32, count, &mut Vec::new())
                    .is_err()
            );
        }
        assert_eq!(visited, expected);
        assert!(bucket.locations(200, 0, &mut Vec::new()).is_err());
        assert!(LocatorBucket::decode(&bytes[4..]).is_err());
    }
}
