//! The seeds file (stages 1 and 2): a bucket table kept in RAM, the sketch part with the color set
//! of every stored seed, and the optional locator part with its seed locations per entry.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::IndexError;
use super::bits::{u32_at, u64_at};
use super::buckets::{LocatorBucket, SketchBucket};
use super::pages::PartTable;
use super::{MAX_BUCKET_BITS, MIN_BUCKET_BITS};
use crate::storage::{ReadAt, will_need_ranges};
use integer_encoding::VarInt;
use std::ops::Range;
use std::sync::Arc;

/// The file name of the seeds file in a shard directory.
pub const SEEDS_FILE: &str = "seeds";

/// Bits of a seed's hash that the seeds file keeps: the bucket bits, then the check value.
pub const SEED_KEY_BITS: u32 = 44;

/// The key of a seed or query k-mer: the low `SEED_KEY_BITS` bits of its hash. Seeds are the
/// smallest hashes of their minimizer windows, so their high bits crowd near zero while the
/// low bits spread evenly over the buckets.
pub fn seed_key(hash: u64) -> u64 {
    hash & ((1 << SEED_KEY_BITS) - 1)
}

pub(crate) const MAGIC: [u8; 16] = *b"jam seeds\0\0\0\0\0\0\0";
pub(crate) const VERSION: u32 = 1;
pub(crate) const HEADER_BYTES: u64 = 256;
/// Pages from 64 bytes (for tests) to 64 KiB, so a start within a page fits 16 bits.
pub(crate) const PAGE_BITS: Range<u32> = 6..17;
/// The sections in file order: bucket table, sketch buckets, color set offsets, color sets,
/// locator buckets. The sketch part is sections 1 to 3, the locator part section 4.
pub(crate) const SECTIONS: usize = 5;
/// Ranges closer than this are read as one.
const JOIN_GAP: u64 = 8 << 10;
/// Joined reads stop growing at this size.
const JOIN_LIMIT: u64 = 1 << 20;

/// The id of a distinct color set; small ids belong to frequent sets.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, PartialOrd, Ord)]
pub struct ColorSetId(pub u32);

/// One stored bin and strand of a seed occurrence.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, PartialOrd, Ord)]
pub struct SeedLocation {
    /// Bin within the entry: the seed starts in bases `bin << bin_bits` onward.
    pub bin: u32,
    /// The canonical k-mer is the reverse complement of the contig strand.
    pub reverse: bool,
}

impl SeedLocation {
    /// The stored value: bin shifted left by one, strand in the low bit.
    pub fn value(self) -> u32 {
        (self.bin << 1) | u32::from(self.reverse)
    }

    /// The seed location of a stored value.
    pub fn from_value(value: u32) -> Self {
        Self {
            bin: value >> 1,
            reverse: value & 1 == 1,
        }
    }
}

/// A stored seed of the sketch part: where its seed locations lie, and its color set.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct StoredSeed {
    /// Its bucket.
    pub bucket: u32,
    /// Its place among the stored seeds of the bucket.
    pub place: u32,
    /// Its color set.
    pub color_set: ColorSetId,
}

/// The color sets that one lookup needs, read together.
#[derive(Clone, Debug, Default)]
pub struct ColorSets {
    ids: Vec<ColorSetId>,
    /// Where each set starts in `entries`, and one more for the end.
    starts: Vec<usize>,
    entries: Vec<u32>,
}

impl ColorSets {
    /// Joins color sets read in parts, whose ids ascend from part to part.
    pub fn concat(parts: Vec<ColorSets>) -> Self {
        let mut joined = Self {
            starts: vec![0],
            ..Self::default()
        };
        for part in parts {
            let offset = joined.entries.len();
            joined.ids.extend(part.ids);
            joined
                .starts
                .extend(part.starts[1..].iter().map(|start| offset + start));
            joined.entries.extend(part.entries);
        }
        joined
    }

    /// The entries of a color set, ascending, or `None` when it was not read.
    pub fn entries(&self, id: ColorSetId) -> Option<&[u32]> {
        let index = self.ids.binary_search(&id).ok()?;
        self.entries.get(self.starts[index]..self.starts[index + 1])
    }
}

/// An open seeds file, read by offset from a local file or S3.
pub struct SeedFile {
    source: Arc<dyn ReadAt>,
    bucket_bits: u32,
    entries: u64,
    seeds: u64,
    color_sets: u64,
    locations: u64,
    sections: [Range<u64>; SECTIONS],
    crcs: [u32; SECTIONS],
    sketch: PartTable,
    locators: Option<PartTable>,
}

impl SeedFile {
    /// Reads and checks the header, and reads the bucket table into RAM.
    pub fn open(source: Arc<dyn ReadAt>) -> Result<Self, IndexError> {
        let mut header = [0u8; HEADER_BYTES as usize];
        source.read_exact_at(0, &mut header)?;
        if header[..16] != MAGIC || u32_at(&header, 16) != Some(VERSION) {
            return Err(IndexError::format("not a seeds file of this version"));
        }
        if u32_at(&header, 252) != Some(crc32fast::hash(&header[..252])) {
            return Err(IndexError::Checksum("seeds header".into()));
        }
        let field = |at: usize| u64_at(&header, at).unwrap_or_default();
        let (bucket_bits, page_bits) = (u32::from(header[20]), u32::from(header[22]));
        let bucket_range = u32::from(MIN_BUCKET_BITS)..=u32::from(MAX_BUCKET_BITS);
        if !bucket_range.contains(&bucket_bits)
            || u32::from(header[21]) != SEED_KEY_BITS
            || !PAGE_BITS.contains(&page_bits)
        {
            return Err(IndexError::format("seeds file parameters"));
        }
        let sections: [Range<u64>; SECTIONS] = std::array::from_fn(|index| {
            let start = field(72 + 16 * index);
            start..start.saturating_add(field(80 + 16 * index))
        });
        let crcs = std::array::from_fn(|index| u32_at(&header, 152 + 4 * index).unwrap_or(0));
        check_sections(&sections, page_bits, source.size())?;
        let mut table = vec![0u8; (sections[0].end - sections[0].start) as usize];
        source.read_exact_at(sections[0].start, &mut table)?;
        if crc32fast::hash(&table) != crcs[0] {
            return Err(IndexError::Checksum("seeds bucket table".into()));
        }
        let buckets = 1usize << bucket_bits;
        let (sketch_pages, locator_pages) = (field(56) as usize, field(64) as usize);
        let sketch_bytes = sections[1].end - sections[1].start;
        let sketch = PartTable::decode(&table, page_bits, buckets, sketch_pages, sketch_bytes)?;
        let locators = match locator_pages {
            0 => None,
            pages => {
                let at = PartTable::encoded_bytes(buckets, sketch_pages);
                let bytes = sections[4].end - sections[4].start;
                let rest = table.get(at..).unwrap_or_default();
                Some(PartTable::decode(rest, page_bits, buckets, pages, bytes)?)
            }
        };
        Ok(Self {
            source,
            bucket_bits,
            entries: field(24),
            seeds: field(32),
            color_sets: field(40),
            locations: field(48),
            sections,
            crcs,
            sketch,
            locators,
        })
    }

    /// Bits of the seed key that select a bucket.
    pub fn bucket_bits(&self) -> u32 {
        self.bucket_bits
    }

    /// The bucket of a seed key: its top bits.
    pub fn bucket(&self, key: u64) -> u32 {
        (key >> (SEED_KEY_BITS - self.bucket_bits)) as u32
    }

    /// The check value of a seed key: its bits after the bucket bits.
    fn check(&self, key: u64) -> u64 {
        key & ((1u64 << (SEED_KEY_BITS - self.bucket_bits)) - 1)
    }

    /// Entries of the shard.
    pub fn entry_count(&self) -> u64 {
        self.entries
    }

    /// Stored seeds: distinct seed keys.
    pub fn seed_count(&self) -> u64 {
        self.seeds
    }

    /// Distinct color sets.
    pub fn color_set_count(&self) -> u64 {
        self.color_sets
    }

    /// Seed locations of the locator part, 0 without it.
    pub fn location_count(&self) -> u64 {
        self.locations
    }

    /// Whether the file has the locator part, which stage 2 needs.
    pub fn has_locators(&self) -> bool {
        self.locators.is_some()
    }

    /// The bytes of the sketch part in the file.
    pub fn sketch_part(&self) -> Range<u64> {
        self.sections[1].start..self.sections[3].end
    }

    /// The bytes of the locator part in the file, if it has one.
    pub fn locator_part(&self) -> Option<Range<u64>> {
        self.locators.as_ref().map(|_| self.sections[4].clone())
    }

    /// Finds the stored seed of each seed key; `keys` must ascend. Every bucket is read once.
    pub fn find_seeds(&self, keys: &[u64]) -> Result<Vec<Option<StoredSeed>>, IndexError> {
        let mut buckets: Vec<u32> = keys.iter().map(|&key| self.bucket(key)).collect();
        buckets.dedup();
        let mut found = Vec::with_capacity(keys.len());
        let mut next = 0;
        self.read_sketch_buckets(&buckets, |bucket, seeds| {
            while let Some(&key) = keys.get(next).filter(|&&key| self.bucket(key) == bucket) {
                found.push(seeds.find(self.check(key)).map(|(place, id)| StoredSeed {
                    bucket,
                    place,
                    color_set: ColorSetId(id),
                }));
                next += 1;
            }
            Ok(())
        })?;
        if found.len() != keys.len() {
            return Err(IndexError::Input("seed keys to look up must ascend".into()));
        }
        Ok(found)
    }

    /// Reads the sketch buckets `buckets` (ascending and distinct) and calls `visit` for each.
    pub fn read_sketch_buckets(
        &self,
        buckets: &[u32],
        mut visit: impl FnMut(u32, SketchBucket<'_>) -> Result<(), IndexError>,
    ) -> Result<(), IndexError> {
        let check_bits = SEED_KEY_BITS - self.bucket_bits;
        let ranges = self.bucket_ranges(&self.sketch, self.sections[1].start, buckets)?;
        self.read_ranges(&ranges, |index, bytes| {
            visit(buckets[index], SketchBucket::decode(bytes, check_bits)?)
        })
    }

    /// Reads the locator buckets `buckets` (ascending and distinct) and calls `visit` for each.
    /// Fails when the file has no locator part.
    pub fn read_locator_buckets(
        &self,
        buckets: &[u32],
        mut visit: impl FnMut(u32, LocatorBucket<'_>) -> Result<(), IndexError>,
    ) -> Result<(), IndexError> {
        let table = self.locators.as_ref().ok_or_else(|| {
            IndexError::Input(
                "the shard has no locator part: it was built with the sketch only".into(),
            )
        })?;
        let ranges = self.bucket_ranges(table, self.sections[4].start, buckets)?;
        self.read_ranges(&ranges, |index, bytes| {
            visit(buckets[index], LocatorBucket::decode(bytes)?)
        })
    }

    /// Reads the color sets with these ids.
    pub fn read_color_sets(
        &self,
        ids: impl IntoIterator<Item = ColorSetId>,
    ) -> Result<ColorSets, IndexError> {
        let mut ids: Vec<ColorSetId> = ids.into_iter().collect();
        ids.sort_unstable();
        ids.dedup();
        if ids
            .last()
            .is_some_and(|id| u64::from(id.0) >= self.color_sets)
        {
            return Err(IndexError::format("color set id past the last set"));
        }
        let offsets = self.sections[2].start;
        let ranges: Vec<Range<u64>> = (ids.iter())
            .map(|id| offsets + 8 * u64::from(id.0)..offsets + 8 * u64::from(id.0) + 16)
            .collect();
        let mut set_ranges = Vec::with_capacity(ids.len());
        let sets = self.sections[3].clone();
        self.read_ranges(&ranges, |_, bytes| {
            let (start, end) = (u64_at(bytes, 0), u64_at(bytes, 8));
            let range = start
                .zip(end)
                .map(|(start, end)| sets.start + start..sets.start + end);
            let range = range.filter(|range| range.start <= range.end && range.end <= sets.end);
            set_ranges.push(range.ok_or_else(|| IndexError::format("color set offsets"))?);
            Ok(())
        })?;
        let mut read = ColorSets {
            ids,
            starts: vec![0],
            entries: Vec::new(),
        };
        self.read_ranges(&set_ranges, |_, bytes| {
            decode_color_set(bytes, self.entries, &mut read.entries)
                .ok_or_else(|| IndexError::format("color set"))?;
            read.starts.push(read.entries.len());
            Ok(())
        })?;
        Ok(read)
    }

    /// Checks the CRC-32 of every section.
    pub fn verify_checksums(&self) -> Result<(), IndexError> {
        let mut buffer = vec![0u8; 1 << 20];
        for (section, crc) in self.sections.iter().zip(self.crcs) {
            let mut hasher = crc32fast::Hasher::new();
            let mut offset = section.start;
            while offset < section.end {
                let length = (section.end - offset).min(buffer.len() as u64) as usize;
                self.source.read_exact_at(offset, &mut buffer[..length])?;
                hasher.update(&buffer[..length]);
                offset += length as u64;
            }
            if hasher.finalize() != crc {
                return Err(IndexError::Checksum("seeds section".into()));
            }
        }
        Ok(())
    }

    /// The file ranges of `buckets` in a part that starts at `start`.
    fn bucket_ranges(
        &self,
        table: &PartTable,
        start: u64,
        buckets: &[u32],
    ) -> Result<Vec<Range<u64>>, IndexError> {
        if !buckets.is_sorted_by(|a, b| a < b) {
            return Err(IndexError::Input("buckets to read must ascend".into()));
        }
        (buckets.iter())
            .map(|&bucket| {
                let range = table.range(bucket as usize).ok_or_else(|| {
                    IndexError::Input(format!("bucket {bucket} is not in this file"))
                })?;
                Ok(start + range.start..start + range.end)
            })
            .collect()
    }

    /// Reads `ranges` (ascending by start) and calls `visit(index, bytes)` for each. Close ranges
    /// are read as one, and the source is told about all reads before the first one.
    fn read_ranges(
        &self,
        ranges: &[Range<u64>],
        mut visit: impl FnMut(usize, &[u8]) -> Result<(), IndexError>,
    ) -> Result<(), IndexError> {
        let mut groups: Vec<(Range<u64>, Range<usize>)> = Vec::new();
        for (index, range) in ranges.iter().enumerate() {
            match groups.last_mut() {
                Some((joined, members))
                    if range.start <= joined.end.saturating_add(JOIN_GAP)
                        && range.end.max(joined.end) - joined.start <= JOIN_LIMIT =>
                {
                    joined.end = joined.end.max(range.end);
                    members.end = index + 1;
                }
                _ => groups.push((range.clone(), index..index + 1)),
            }
        }
        let mut hints: Vec<(u64, u64)> = groups.iter().map(|(r, _)| (r.start, r.end)).collect();
        will_need_ranges(self.source.as_ref(), &mut hints);
        let mut buffer = Vec::new();
        for (joined, members) in groups {
            buffer.resize((joined.end - joined.start) as usize, 0);
            self.source.read_exact_at(joined.start, &mut buffer)?;
            for index in members {
                let range = &ranges[index];
                let at = (range.start - joined.start) as usize;
                visit(index, &buffer[at..at + (range.end - range.start) as usize])?;
            }
        }
        Ok(())
    }
}

/// Checks that the sections follow each other in order within `size` bytes, with the table
/// first and each part on a page boundary.
fn check_sections(
    sections: &[Range<u64>; SECTIONS],
    page_bits: u32,
    size: u64,
) -> Result<(), IndexError> {
    let on_page = |offset: u64| offset.is_multiple_of(1 << page_bits);
    let mut end = HEADER_BYTES;
    for section in sections {
        if section.start < end || section.end < section.start {
            return Err(IndexError::format("seeds sections"));
        }
        end = section.end;
    }
    let locators = &sections[4];
    let ordered = sections[0].start == HEADER_BYTES
        && on_page(sections[1].start)
        && (locators.is_empty() || on_page(locators.start))
        && end == size;
    if !ordered {
        return Err(IndexError::format("seeds sections"));
    }
    Ok(())
}

/// Appends a color set (entries ascending and distinct) in its file encoding: the count, then
/// the gap before each entry, all varints.
pub(crate) fn encode_color_set(entries: &[u32], out: &mut Vec<u8>) {
    let mut buffer = [0u8; 10];
    let mut put = |value: u64| {
        let used = value.encode_var(&mut buffer);
        out.extend_from_slice(&buffer[..used]);
    };
    put(entries.len() as u64);
    let mut next = 0u64;
    for &entry in entries {
        put(u64::from(entry) - next);
        next = u64::from(entry) + 1;
    }
}

/// Decodes one color set (count, then gaps) of a shard of `entries` entries into `out`.
fn decode_color_set(bytes: &[u8], entries: u64, out: &mut Vec<u32>) -> Option<()> {
    let (count, mut at) = u64::decode_var(bytes)?;
    let mut next = 0u64;
    for _ in 0..count {
        let (gap, used) = u64::decode_var(bytes.get(at..)?)?;
        at += used;
        let entry = next.checked_add(gap).filter(|entry| *entry < entries)?;
        out.push(entry as u32);
        next = entry + 1;
    }
    (at == bytes.len() && count > 0).then_some(())
}
