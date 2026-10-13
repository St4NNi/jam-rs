//! The seeds of an index file (stages 1 and 2), before its manifest: a bucket table in RAM, the
//! sketch part with color sets, and the optional locator part with locations and contig tables.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::IndexError;
use super::bits::{u32_at, u64_at};
use super::buckets::{LocatorBucket, SketchBucket};
use super::contig_table::ContigTable;
use super::pages::PartTable;
use super::{MAX_BUCKET_BITS, MIN_BUCKET_BITS};
use crate::storage::{ReadAt, will_need_ranges};
use integer_encoding::VarInt;
use std::ops::Range;
use std::sync::Arc;

/// Bits of a seed's hash that the seeds file keeps: the bucket bits, then the check value.
pub const SEED_KEY_BITS: u32 = 44;

/// Bits of the top `SEED_KEY_BITS` bits of a hash that move to the end of the key.
const KEY_TAIL_BITS: u32 = SEED_KEY_BITS - MAX_BUCKET_BITS as u32;

/// The key of a seed or query k-mer: the top `SEED_KEY_BITS` bits of its hash, with their first
/// `KEY_TAIL_BITS` bits moved to the end. Seeds are the smallest hashes of their minimizer
/// windows, so the first bits crowd near zero; the bits after them select even buckets.
pub fn seed_key(hash: u64) -> u64 {
    let top = hash >> (64 - SEED_KEY_BITS);
    let tail = top >> (SEED_KEY_BITS - KEY_TAIL_BITS);
    ((top << KEY_TAIL_BITS) & ((1 << SEED_KEY_BITS) - 1)) | tail
}

pub(crate) const MAGIC: [u8; 16] = *b"jam seeds\0\0\0\0\0\0\0";
pub(crate) const VERSION: u32 = 1;
pub(crate) const HEADER_BYTES: u64 = 256;
/// Pages from 64 bytes (for tests) to 64 KiB, so a start within a page fits 16 bits.
pub(crate) const PAGE_BITS: Range<u32> = 6..17;
/// The sections in file order: bucket table, sketch buckets, color set offsets, color sets,
/// locator buckets, contig tables. The sketch part is sections 1 to 3, the locator part 4 and 5.
pub(crate) const SECTIONS: usize = 6;
/// Where the CRC-32 of each section lies in the header.
pub(crate) const CRCS_AT: usize = 72 + 16 * SECTIONS;
/// Where the offset and bytes of the manifest (u64 each) and its CRC-32 lie in the header.
pub(crate) const MANIFEST_AT: usize = CRCS_AT + 4 * SECTIONS;
/// Ranges closer than this are read as one. On CephFS one more read costs about as much system
/// time as copying 64 KiB more (measured on worm), and on S3 each read is one request.
const JOIN_GAP: u64 = 64 << 10;
/// Joined reads stop growing at this size, so a thread holds at most this much at once.
const JOIN_LIMIT: u64 = 4 << 20;
/// Bytes of the CRC-32 that ends every non-empty bucket, color set and contig table, so each
/// read checks what it uses. It covers the unit's index, then its bytes.
pub(crate) const UNIT_CHECK_BYTES: usize = 4;

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
    sections: [Range<u64>; SECTIONS],
    #[cfg(test)]
    crcs: [u32; SECTIONS],
    /// The bytes of the manifest and their CRC-32.
    manifest: (Range<u64>, u32),
    sketch: PartTable,
    locators: Option<PartTable>,
    /// Where the contig table of each entry starts in the file, then where the last one ends;
    /// empty without the locator part.
    table_starts: Vec<u64>,
}

impl SeedFile {
    /// Reads and checks the header, and reads the bucket table into RAM.
    pub fn open(source: Arc<dyn ReadAt>) -> Result<Self, IndexError> {
        let mut header = [0u8; HEADER_BYTES as usize];
        source.read_exact_at(0, &mut header)?;
        let (sections, crcs) = header_sections(&header)?;
        let manifest = header_manifest(&header);
        let field = |at: usize| u64_at(&header, at).unwrap_or_default();
        let (bucket_bits, page_bits) = (u32::from(header[20]), u32::from(header[22]));
        let bucket_range = u32::from(MIN_BUCKET_BITS)..=u32::from(MAX_BUCKET_BITS);
        if !bucket_range.contains(&bucket_bits)
            || u32::from(header[21]) != SEED_KEY_BITS
            || !PAGE_BITS.contains(&page_bits)
        {
            return Err(IndexError::format("seeds file parameters"));
        }
        check_sections(&sections, page_bits, &manifest.0, source.size())?;
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
        let table_starts = match locators {
            Some(_) => read_table_starts(source.as_ref(), &sections[5], field(24))?,
            None => Vec::new(),
        };
        Ok(Self {
            source,
            bucket_bits,
            entries: field(24),
            seeds: field(32),
            color_sets: field(40),
            sections,
            #[cfg(test)]
            crcs,
            manifest,
            sketch,
            locators,
            table_starts,
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
        (self.locators.as_ref()).map(|_| self.sections[4].start..self.sections[5].end)
    }

    /// The contig table of an entry; an error without the locator part.
    pub fn contig_table(&self, entry: u32) -> Result<ContigTable, IndexError> {
        ContigTable::decode(&self.contig_table_bytes(entry)?)
    }

    /// The encoded contig table of an entry, as a reshard copies it.
    pub(crate) fn contig_table_bytes(&self, entry: u32) -> Result<Vec<u8>, IndexError> {
        let range = self.contig_table_range(entry)?;
        let mut bytes = vec![0u8; (range.end - range.start) as usize];
        self.source.read_exact_at(range.start, &mut bytes)?;
        let length = checked_unit(&bytes, entry, "contig table")?.len();
        bytes.truncate(length);
        Ok(bytes)
    }

    /// Tells the source that the contig table of an entry will be read soon.
    pub fn hint_contig_table(&self, entry: u32) {
        if let Ok(range) = self.contig_table_range(entry) {
            self.source.will_need(range.start, range.end - range.start);
        }
    }

    fn contig_table_range(&self, entry: u32) -> Result<Range<u64>, IndexError> {
        let starts = &self.table_starts;
        match (starts.get(entry as usize), starts.get(entry as usize + 1)) {
            (Some(&start), Some(&end)) => Ok(start..end),
            _ if starts.is_empty() => Err(IndexError::Input(
                "the shard was built with the sketch part only".into(),
            )),
            _ => Err(IndexError::Input(format!(
                "entry {entry} is not in this shard"
            ))),
        }
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
            let bytes = checked_unit(bytes, buckets[index], "sketch bucket")?;
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
            IndexError::Input("the shard was built with the sketch part only".into())
        })?;
        let ranges = self.bucket_ranges(table, self.sections[4].start, buckets)?;
        self.read_ranges(&ranges, |index, bytes| {
            let bytes = checked_unit(bytes, buckets[index], "locator bucket")?;
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
        self.read_ranges(&set_ranges, |index, bytes| {
            let bytes = checked_unit(bytes, read.ids[index].0, "color set")?;
            decode_color_set(bytes, self.entries, &mut read.entries)
                .ok_or_else(|| IndexError::format("color set"))?;
            read.starts.push(read.entries.len());
            Ok(())
        })?;
        Ok(read)
    }

    /// The bytes of the manifest in the file.
    pub fn manifest_range(&self) -> Range<u64> {
        self.manifest.0.clone()
    }

    /// Checks the CRC-32 of every section and of the manifest.
    #[cfg(test)]
    pub fn verify_checksums(&self) -> Result<(), IndexError> {
        let sections = self.sections.iter().cloned().zip(self.crcs);
        for (section, crc) in sections.chain([self.manifest.clone()]) {
            if crc_of(self.source.as_ref(), section)? != crc {
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
        // The buffer only grows, so its bytes are zeroed once, not before every read.
        let mut buffer = Vec::new();
        for (joined, members) in groups {
            let length = (joined.end - joined.start) as usize;
            if buffer.len() < length {
                buffer.resize(length, 0);
            }
            let read = &mut buffer[..length];
            self.source.read_exact_at(joined.start, read)?;
            for index in members {
                let range = &ranges[index];
                let at = (range.start - joined.start) as usize;
                visit(index, &read[at..at + (range.end - range.start) as usize])?;
            }
        }
        Ok(())
    }
}

/// The byte range and the CRC-32 of each section, from a seeds file header. Fails on another
/// format or a damaged header.
pub(crate) fn header_sections(
    header: &[u8; HEADER_BYTES as usize],
) -> Result<([Range<u64>; SECTIONS], [u32; SECTIONS]), IndexError> {
    if header[..16] != MAGIC || u32_at(header, 16) != Some(VERSION) {
        return Err(IndexError::format("not a seeds file of this version"));
    }
    if u32_at(header, 252) != Some(crc32fast::hash(&header[..252])) {
        return Err(IndexError::Checksum("seeds header".into()));
    }
    let field = |at: usize| u64_at(header, at).unwrap_or_default();
    let sections = std::array::from_fn(|index| {
        let start = field(72 + 16 * index);
        start..start.saturating_add(field(80 + 16 * index))
    });
    let crcs = std::array::from_fn(|index| u32_at(header, CRCS_AT + 4 * index).unwrap_or(0));
    Ok((sections, crcs))
}

/// The bytes of the manifest and their CRC-32, from a header checked by `header_sections`.
pub(crate) fn header_manifest(header: &[u8; HEADER_BYTES as usize]) -> (Range<u64>, u32) {
    let field = |at: usize| u64_at(header, at).unwrap_or_default();
    let start = field(MANIFEST_AT);
    let crc = u32_at(header, MANIFEST_AT + 16).unwrap_or_default();
    (start..start.saturating_add(field(MANIFEST_AT + 8)), crc)
}

/// Checks the CRC-32 of every section of a seeds file with `header` that lies within `range`.
/// `part` holds the bytes of `range`, from its offset 0, as a pinned range does.
pub(crate) fn verify_part(
    header: &[u8; HEADER_BYTES as usize],
    part: &dyn ReadAt,
    range: Range<u64>,
) -> Result<(), IndexError> {
    let (sections, crcs) = header_sections(header)?;
    let inside = |section: &Range<u64>| range.start <= section.start && section.end <= range.end;
    for (section, crc) in sections.iter().zip(crcs).filter(|(s, _)| inside(s)) {
        let local = section.start - range.start..section.end - range.start;
        if crc_of(part, local)? != crc {
            return Err(IndexError::Checksum("seeds section".into()));
        }
    }
    Ok(())
}

/// The CRC-32 of the bytes `range` of `source`, read in parts of 1 MiB.
fn crc_of(source: &dyn ReadAt, range: Range<u64>) -> Result<u32, IndexError> {
    let mut buffer = vec![0u8; (range.end - range.start).min(1 << 20) as usize];
    let mut hasher = crc32fast::Hasher::new();
    let mut offset = range.start;
    while offset < range.end {
        let length = (range.end - offset).min(buffer.len() as u64) as usize;
        source.read_exact_at(offset, &mut buffer[..length])?;
        hasher.update(&buffer[..length]);
        offset += length as u64;
    }
    Ok(hasher.finalize())
}

/// The CRC-32 that follows the stored unit `bytes` with index `index`.
pub(crate) fn unit_check(index: u32, bytes: &[u8]) -> [u8; UNIT_CHECK_BYTES] {
    let mut hasher = crc32fast::Hasher::new();
    hasher.update(&index.to_le_bytes());
    hasher.update(bytes);
    hasher.finalize().to_le_bytes()
}

/// The bytes of the stored unit `what` with index `index`, without its CRC-32, once that
/// matches. An empty unit has no CRC-32.
fn checked_unit<'a>(stored: &'a [u8], index: u32, what: &str) -> Result<&'a [u8], IndexError> {
    if stored.is_empty() {
        return Ok(stored);
    }
    let damaged = || IndexError::Checksum(format!("seeds {what} {index}"));
    let split = (stored.len().checked_sub(UNIT_CHECK_BYTES)).ok_or_else(damaged)?;
    let (bytes, check) = stored.split_at(split);
    if unit_check(index, bytes) != check {
        return Err(damaged());
    }
    Ok(bytes)
}

/// Reads where the contig table of each of `entries` entries starts, from the front of the
/// contig tables section, and checks that they ascend to the section end.
fn read_table_starts(
    source: &dyn ReadAt,
    section: &Range<u64>,
    entries: u64,
) -> Result<Vec<u64>, IndexError> {
    let count = usize::try_from(entries + 1).map_err(|_| IndexError::format("seeds entries"))?;
    let mut bytes = vec![0u8; 8 * count];
    if 8 * count as u64 > section.end - section.start {
        return Err(IndexError::format("seeds contig tables"));
    }
    source.read_exact_at(section.start, &mut bytes)?;
    let starts: Vec<u64> = (bytes.as_chunks::<8>().0.iter())
        .map(|word| section.start.saturating_add(u64::from_le_bytes(*word)))
        .collect();
    let valid = starts.first() == Some(&(section.start + 8 * count as u64))
        && starts.last() == Some(&section.end)
        && starts.is_sorted();
    if !valid {
        return Err(IndexError::format("seeds contig tables"));
    }
    Ok(starts)
}

/// Checks that the sections follow each other in order, with the table first and each part on a
/// page boundary, and that the manifest follows them up to the end of the file of `size` bytes.
fn check_sections(
    sections: &[Range<u64>; SECTIONS],
    page_bits: u32,
    manifest: &Range<u64>,
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
        && manifest.start == end
        && manifest.end > manifest.start
        && manifest.end == size;
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
