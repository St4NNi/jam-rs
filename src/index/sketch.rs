//! The sketch file (stage 1): the color set of every stored seed, found from its hash, and every
//! distinct color set once.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::bits::{BIT_PADDING, BitWriter, bits_needed, read_bits, u16_at, u32_at, u64_at};
use super::{AtPath, FINGERPRINT_BITS, IndexError, MAX_BUCKET_BITS, MIN_BUCKET_BITS, seed_key};
use crate::storage::MappedFile;
use integer_encoding::VarInt;
use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

/// The file name of the sketch in a shard directory.
pub const SKETCH_FILE: &str = "sketch";

const MAGIC: [u8; 16] = *b"jam sketch\0\0\0\0\0\0";
const VERSION: u32 = 1;
const HEADER_BYTES: usize = 128;
const BLOCK_HEADER_BYTES: usize = 8;
const MAX_ESCAPES: usize = (1 << 24) - 1;

/// The id of a distinct color set; small ids belong to frequent sets.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, PartialOrd, Ord)]
pub struct ColorSetId(pub u32);

/// An open sketch file.
pub struct Sketch {
    file: MappedFile,
    bucket_bits: u32,
    entries: u64,
    seeds: u64,
    color_sets: u64,
    crcs: [u32; 4],
    layout: Layout,
}

impl Sketch {
    /// Reads the whole sketch into RAM, so lookups never wait for storage.
    pub fn load(&self) -> Result<(), IndexError> {
        self.file.load().map_err(IndexError::from)
    }

    /// Maps the sketch at `path` and checks its header and size. Pages are read when a lookup
    /// first needs them; `load` reads them all at once.
    pub fn open(path: &Path) -> Result<Self, IndexError> {
        let file = MappedFile::open(path, false)?;
        let bytes = file.bytes();
        let header = bytes
            .get(..HEADER_BYTES)
            .ok_or_else(|| IndexError::format("sketch header is cut off"))?;
        let field = |offset| u64_at(header, offset).unwrap_or_default();
        if header[..16] != MAGIC || u32_at(header, 16) != Some(VERSION) {
            return Err(IndexError::format("not a sketch file of this version"));
        }
        if u32_at(header, 124) != Some(crc32fast::hash(&header[..124])) {
            return Err(IndexError::Checksum("sketch header".into()));
        }
        let bucket_bits = u32::from(header[20]);
        let bucket_range = u32::from(MIN_BUCKET_BITS)..=u32::from(MAX_BUCKET_BITS);
        if !bucket_range.contains(&bucket_bits) || header[21] != FINGERPRINT_BITS {
            return Err(IndexError::format("sketch parameters"));
        }
        let (entries, seeds, color_sets) = (field(24), field(32), field(40));
        let layout = Layout::new(bucket_bits, field(48), color_sets, field(56))
            .filter(|layout| layout.end == bytes.len())
            .ok_or_else(|| IndexError::format("sketch size does not match its header"))?;
        let crcs = [64, 68, 72, 76].map(|offset| u32_at(header, offset).unwrap_or_default());
        Ok(Self {
            file,
            bucket_bits,
            entries,
            seeds,
            color_sets,
            crcs,
            layout,
        })
    }

    /// The color set of the stored seed with this hash's bucket and fingerprint, if any.
    pub fn color_set(&self, hash: u64) -> Option<ColorSetId> {
        let key = seed_key(hash, self.bucket_bits + u32::from(FINGERPRINT_BITS));
        let bucket = (key >> FINGERPRINT_BITS) as usize;
        let bytes = self.file.bytes();
        let at = self.layout.directory + 8 * bucket;
        let start = usize::try_from(u64_at(bytes, at)?).ok()?;
        let end = usize::try_from(u64_at(bytes, at + 8)?).ok()?;
        if start >= end {
            return None;
        }
        let block = bytes.get(self.layout.buckets + start..self.layout.buckets_end)?;
        find_in_block(block, key as u16)
    }

    /// The distinct color sets.
    pub fn color_sets(&self) -> ColorSets<'_> {
        let bytes = self.file.bytes();
        ColorSets {
            offsets: &bytes[self.layout.set_offsets..self.layout.sets],
            sets: &bytes[self.layout.sets..],
        }
    }

    /// Bits of the hash that select a bucket.
    pub fn bucket_bits(&self) -> u32 {
        self.bucket_bits
    }

    /// Entries of the shard.
    pub fn entry_count(&self) -> u64 {
        self.entries
    }

    /// Stored seeds; seeds that share bucket and fingerprint count once.
    pub fn seed_count(&self) -> u64 {
        self.seeds
    }

    /// Distinct color sets.
    pub fn color_set_count(&self) -> u64 {
        self.color_sets
    }

    /// Checks the CRC-32 of every section.
    pub fn verify_checksums(&self) -> Result<(), IndexError> {
        let bytes = self.file.bytes();
        for (range, crc) in self.layout.sections().into_iter().zip(self.crcs) {
            if crc32fast::hash(&bytes[range]) != crc {
                return Err(IndexError::Checksum("sketch section".into()));
            }
        }
        Ok(())
    }
}

/// The distinct color sets of a sketch.
#[derive(Clone, Copy)]
pub struct ColorSets<'a> {
    offsets: &'a [u8],
    sets: &'a [u8],
}

impl<'a> ColorSets<'a> {
    /// The entries of a color set, ascending, or `None` for an unknown id.
    pub fn entries(&self, id: ColorSetId) -> Option<ColorSetEntries<'a>> {
        let at = 8 * id.0 as usize;
        let start = usize::try_from(u64_at(self.offsets, at)?).ok()?;
        let end = usize::try_from(u64_at(self.offsets, at + 8)?).ok()?;
        let bytes = self.sets.get(start..end)?;
        let (count, at) = u64::decode_var(bytes)?;
        Some(ColorSetEntries {
            bytes,
            at,
            left: count,
            count,
            next: 0,
        })
    }
}

/// The entries of one color set, ascending.
pub struct ColorSetEntries<'a> {
    bytes: &'a [u8],
    at: usize,
    left: u64,
    count: u64,
    next: u64,
}

impl ColorSetEntries<'_> {
    /// Entries in the set.
    pub fn len(&self) -> u64 {
        self.count
    }

    /// Whether the set has no entries (never true in a valid sketch).
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }
}

impl Iterator for ColorSetEntries<'_> {
    type Item = u32;

    fn next(&mut self) -> Option<u32> {
        if self.left == 0 {
            return None;
        }
        let (gap, used) = u64::decode_var(self.bytes.get(self.at..)?)?;
        self.at += used;
        let value = self.next.checked_add(gap)?;
        self.left -= 1;
        self.next = value + 1;
        u32::try_from(value).ok()
    }
}

/// Appends a color set (entries ascending and distinct) in its file encoding.
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

/// Appends the block of one bucket: its seeds as (fingerprint, color set id), fingerprints
/// ascending and distinct. `out` may already hold other blocks.
pub(crate) fn encode_bucket(seeds: &[(u16, u32)], out: &mut Vec<u8>) {
    if seeds.is_empty() {
        return;
    }
    let escaped = out.len() + BLOCK_HEADER_BYTES + 2 * seeds.len();
    let ids: Vec<u32> = seeds.iter().map(|(_, id)| *id).collect();
    let (width, escapes) = choose_width(&ids);
    out.extend_from_slice(&(seeds.len() as u32).to_le_bytes());
    out.extend_from_slice(&((width << 24) | escapes as u32).to_le_bytes());
    for (fingerprint, _) in seeds {
        out.extend_from_slice(&fingerprint.to_le_bytes());
    }
    let direct = (1u64 << width) - escapes as u64;
    let mut values = Vec::with_capacity(ids.len());
    for &id in &ids {
        if u64::from(id) < direct {
            values.push(u64::from(id));
        } else {
            values.push(direct + (out.len() - escaped) as u64 / 4);
            out.extend_from_slice(&id.to_le_bytes());
        }
    }
    let mut writer = BitWriter::new(out);
    for value in values {
        writer.push(value, width);
    }
    writer.finish();
}

/// Writes a sketch file bucket by bucket.
pub(crate) struct SketchWriter {
    file: BufWriter<File>,
    path: PathBuf,
    bucket_bits: u32,
    entries: u64,
    directory: Vec<u64>,
    written: u64,
    seeds: u64,
    crc: crc32fast::Hasher,
}

impl SketchWriter {
    /// Creates the file at `path` for a shard of `entries` entries.
    pub(crate) fn create(path: &Path, bucket_bits: u32, entries: u64) -> Result<Self, IndexError> {
        let mut file = File::create_new(path).at(path)?;
        let buckets = 1usize << bucket_bits;
        let start = (HEADER_BYTES + 8 * (buckets + 1)) as u64;
        file.seek(SeekFrom::Start(start)).at(path)?;
        Ok(Self {
            file: BufWriter::with_capacity(1 << 20, file),
            path: path.to_path_buf(),
            bucket_bits,
            entries,
            directory: Vec::with_capacity(buckets + 1),
            written: 0,
            seeds: 0,
            crc: crc32fast::Hasher::new(),
        })
    }

    /// Appends the block of `bucket` made by `encode_bucket`. Buckets come in ascending order.
    pub(crate) fn push_bucket(&mut self, bucket: usize, block: &[u8]) -> Result<(), IndexError> {
        if bucket < self.directory.len() || bucket >= 1 << self.bucket_bits {
            return Err(IndexError::format("sketch buckets out of order"));
        }
        self.directory.resize(bucket + 1, self.written);
        self.seeds += u64::from(u32_at(block, 0).unwrap_or_default());
        self.file.write_all(block).at(&self.path)?;
        self.crc.update(block);
        self.written += block.len() as u64;
        Ok(())
    }

    /// Writes the color sets (encoded by `encode_color_set`, in id order) and the header.
    pub(crate) fn finish(mut self, sets: &[impl AsRef<[u8]>]) -> Result<(), IndexError> {
        self.directory
            .resize((1 << self.bucket_bits) + 1, self.written);
        let padding = vec![
            0u8;
            (self.written as usize + BIT_PADDING).next_multiple_of(8)
                - self.written as usize
        ];
        let path = self.path.clone();
        self.file.write_all(&padding).at(&path)?;
        self.crc.update(&padding);
        let mut offsets = Vec::with_capacity(8 * (sets.len() + 1));
        let mut set_bytes = 0u64;
        offsets.extend_from_slice(&0u64.to_le_bytes());
        for set in sets {
            set_bytes += set.as_ref().len() as u64;
            offsets.extend_from_slice(&set_bytes.to_le_bytes());
        }
        self.file.write_all(&offsets).at(&path)?;
        let mut set_crc = crc32fast::Hasher::new();
        for set in sets {
            self.file.write_all(set.as_ref()).at(&path)?;
            set_crc.update(set.as_ref());
        }
        let directory: Vec<u8> = self
            .directory
            .iter()
            .flat_map(|o| o.to_le_bytes())
            .collect();
        let mut header = vec![0u8; HEADER_BYTES];
        header[..16].copy_from_slice(&MAGIC);
        header[16..20].copy_from_slice(&VERSION.to_le_bytes());
        header[20] = self.bucket_bits as u8;
        header[21] = FINGERPRINT_BITS;
        let fields = [
            self.entries,
            self.seeds,
            sets.len() as u64,
            self.written,
            set_bytes,
        ];
        for (index, value) in fields.into_iter().enumerate() {
            header[24 + 8 * index..32 + 8 * index].copy_from_slice(&value.to_le_bytes());
        }
        let crcs = [
            crc32fast::hash(&directory),
            self.crc.finalize(),
            crc32fast::hash(&offsets),
            set_crc.finalize(),
        ];
        for (index, crc) in crcs.into_iter().enumerate() {
            header[64 + 4 * index..68 + 4 * index].copy_from_slice(&crc.to_le_bytes());
        }
        let header_crc = crc32fast::hash(&header[..124]);
        header[124..].copy_from_slice(&header_crc.to_le_bytes());
        let mut file = self
            .file
            .into_inner()
            .map_err(|error| error.into_error())
            .at(&path)?;
        file.seek(SeekFrom::Start(0)).at(&path)?;
        file.write_all(&header).at(&path)?;
        file.write_all(&directory).at(&path)?;
        file.sync_all().at(&path)?;
        Ok(())
    }
}

/// Where the sections of a sketch file lie.
#[derive(Clone, Copy, Debug)]
struct Layout {
    directory: usize,
    buckets: usize,
    buckets_end: usize,
    set_offsets: usize,
    sets: usize,
    end: usize,
}

impl Layout {
    fn new(bucket_bits: u32, bucket_bytes: u64, sets: u64, set_bytes: u64) -> Option<Self> {
        let directory = HEADER_BYTES;
        let buckets = directory.checked_add(8usize.checked_mul((1usize << bucket_bits) + 1)?)?;
        let padded = usize::try_from(bucket_bytes)
            .ok()?
            .checked_add(BIT_PADDING)?
            .next_multiple_of(8);
        let buckets_end = buckets.checked_add(padded)?;
        let offsets = usize::try_from(sets).ok()?.checked_add(1)?.checked_mul(8)?;
        let set_offsets = buckets_end;
        let sets = set_offsets.checked_add(offsets)?;
        let end = sets.checked_add(usize::try_from(set_bytes).ok()?)?;
        Some(Self {
            directory,
            buckets,
            buckets_end,
            set_offsets,
            sets,
            end,
        })
    }

    fn sections(&self) -> [std::ops::Range<usize>; 4] {
        [
            self.directory..self.buckets,
            self.buckets..self.buckets_end,
            self.set_offsets..self.sets,
            self.sets..self.end,
        ]
    }
}

/// Finds a fingerprint in a bucket block and decodes its color set id.
fn find_in_block(block: &[u8], fingerprint: u16) -> Option<ColorSetId> {
    let count = u32_at(block, 0)? as usize;
    let packed = u32_at(block, 4)?;
    let (width, escapes) = (packed >> 24, (packed & 0x00ff_ffff) as usize);
    let fingerprints = block.get(BLOCK_HEADER_BYTES..BLOCK_HEADER_BYTES + 2 * count)?;
    let (mut low, mut high) = (0, count);
    while low < high {
        let middle = (low + high) / 2;
        if u16_at(fingerprints, 2 * middle)? < fingerprint {
            low = middle + 1;
        } else {
            high = middle;
        }
    }
    if u16_at(fingerprints, 2 * low)? != fingerprint || width > 32 {
        return None;
    }
    let escaped = BLOCK_HEADER_BYTES + 2 * count;
    let values = block.get(escaped + 4 * escapes..)?;
    let value = read_bits(values, low as u64 * u64::from(width), width);
    let direct = (1u64 << width).checked_sub(escapes as u64)?;
    if value < direct {
        return Some(ColorSetId(value as u32));
    }
    let index = usize::try_from(value - direct).ok()?;
    u32_at(block, escaped + 4 * index).map(ColorSetId)
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

#[cfg(test)]
#[path = "sketch_tests.rs"]
mod tests;
