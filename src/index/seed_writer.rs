//! Writes the seeds file bucket by bucket: each part goes to a temporary file, and both are joined
//! behind the header and the bucket table at the end.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::bits::u32_at;
use super::pages::{PartLayout, PartTable};
use super::seed_file::{HEADER_BYTES, MAGIC, PAGE_BITS, SECTIONS, SEED_KEY_BITS, VERSION};
use super::{AtPath, IndexError};
use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

/// The shape of a seeds file, fixed before its first bucket.
#[derive(Clone, Copy, Debug)]
pub(crate) struct SeedFileShape {
    pub(crate) bucket_bits: u32,
    /// Bytes per page, as a power of two.
    pub(crate) page_bits: u32,
    pub(crate) entries: u64,
    /// Whether the locator part is written.
    pub(crate) locators: bool,
}

/// One part, written to its own temporary file.
struct PartFile {
    file: BufWriter<File>,
    path: PathBuf,
    layout: PartLayout,
    crc: crc32fast::Hasher,
}

impl PartFile {
    fn create(path: PathBuf, page_bits: u32) -> Result<Self, IndexError> {
        let file = File::create_new(&path).at(&path)?;
        Ok(Self {
            file: BufWriter::with_capacity(1 << 20, file),
            path,
            layout: PartLayout::new(page_bits),
            crc: crc32fast::Hasher::new(),
        })
    }

    /// Appends one bucket, after the padding that keeps it within a page.
    fn push(&mut self, bucket: &[u8]) -> Result<(), IndexError> {
        let padding = vec![0u8; self.layout.place(bucket.len() as u64) as usize];
        for bytes in [&padding[..], bucket] {
            self.file.write_all(bytes).at(&self.path)?;
            self.crc.update(bytes);
        }
        Ok(())
    }

    fn finish(self) -> Result<FinishedPart, IndexError> {
        let path = self.path;
        let mut file = self
            .file
            .into_inner()
            .map_err(|error| error.into_error())
            .at(&path)?;
        file.seek(SeekFrom::Start(0)).at(&path)?;
        Ok(FinishedPart {
            bytes: self.layout.written(),
            table: self.layout.finish(),
            crc: self.crc.finalize(),
            file,
        })
    }
}

/// Writes a seeds file; buckets come in ascending order.
pub(crate) struct SeedWriter {
    path: PathBuf,
    shape: SeedFileShape,
    sketch: PartFile,
    locators: Option<PartFile>,
    next_bucket: usize,
    seeds: u64,
    locations: u64,
    _temporary: tempfile::TempDir,
}

impl SeedWriter {
    /// Creates the file at `path`, with temporary files in `temporary_directory`.
    pub(crate) fn create(
        path: &Path,
        temporary_directory: &Path,
        shape: SeedFileShape,
    ) -> Result<Self, IndexError> {
        if !PAGE_BITS.contains(&shape.page_bits) {
            return Err(IndexError::Input("pages must be 64 bytes to 64 KiB".into()));
        }
        let temporary = tempfile::Builder::new()
            .prefix("jam-seeds")
            .tempdir_in(temporary_directory)
            .at(temporary_directory)?;
        let part = |name: &str| PartFile::create(temporary.path().join(name), shape.page_bits);
        let sketch = part("sketch")?;
        let locators = shape.locators.then(|| part("locators")).transpose()?;
        File::create_new(path).at(path)?;
        Ok(Self {
            path: path.to_path_buf(),
            shape,
            sketch,
            locators,
            next_bucket: 0,
            seeds: 0,
            locations: 0,
            _temporary: temporary,
        })
    }

    /// Appends the sketch bucket and the locator bucket of `bucket`; buckets skipped since the
    /// last one are empty. The locator bucket is ignored without a locator part.
    pub(crate) fn push_bucket(
        &mut self,
        bucket: usize,
        sketch: &[u8],
        locators: &[u8],
    ) -> Result<(), IndexError> {
        if bucket < self.next_bucket || bucket >= 1 << self.shape.bucket_bits {
            return Err(IndexError::format("seeds buckets out of order"));
        }
        while self.next_bucket < bucket {
            self.push_one(&[], &[])?;
        }
        self.push_one(sketch, locators)
    }

    fn push_one(&mut self, sketch: &[u8], locators: &[u8]) -> Result<(), IndexError> {
        self.seeds += u64::from(u32_at(sketch, 0).unwrap_or_default());
        self.sketch.push(sketch)?;
        if let Some(part) = &mut self.locators {
            self.locations += u64::from(u32_at(locators, 4).unwrap_or_default());
            part.push(locators)?;
        }
        self.next_bucket += 1;
        Ok(())
    }

    /// Writes the color sets (encoded, in id order), the header and the table, and joins the
    /// parts into the file.
    pub(crate) fn finish(mut self, color_sets: &[impl AsRef<[u8]>]) -> Result<(), IndexError> {
        while self.next_bucket < 1 << self.shape.bucket_bits {
            self.push_one(&[], &[])?;
        }
        let header = HeaderFields {
            shape: self.shape,
            seeds: self.seeds,
            locations: self.locations,
            color_sets: color_sets.len() as u64,
        };
        let sketch = self.sketch.finish()?;
        let locators = self.locators.map(PartFile::finish).transpose()?;
        let mut table = Vec::new();
        sketch.table.encode(&mut table);
        if let Some(part) = &locators {
            part.table.encode(&mut table);
        }
        let mut offsets = Vec::with_capacity(8 * (color_sets.len() + 1));
        let mut set_bytes = 0u64;
        offsets.extend_from_slice(&0u64.to_le_bytes());
        let mut set_crc = crc32fast::Hasher::new();
        for set in color_sets {
            set_bytes += set.as_ref().len() as u64;
            set_crc.update(set.as_ref());
            offsets.extend_from_slice(&set_bytes.to_le_bytes());
        }
        let page = 1u64 << self.shape.page_bits;
        let mut sections = [(0u64, 0u64); SECTIONS];
        sections[0] = (HEADER_BYTES, table.len() as u64);
        let sketch_at = (HEADER_BYTES + table.len() as u64).next_multiple_of(page);
        sections[1] = (sketch_at, sketch.bytes);
        sections[2] = (
            (sketch_at + sketch.bytes).next_multiple_of(8),
            offsets.len() as u64,
        );
        sections[3] = (sections[2].0 + offsets.len() as u64, set_bytes);
        let sets_end = sections[3].0 + set_bytes;
        sections[4] = match &locators {
            Some(part) => (sets_end.next_multiple_of(page), part.bytes),
            None => (sets_end, 0),
        };
        let crcs = [
            crc32fast::hash(&table),
            sketch.crc,
            crc32fast::hash(&offsets),
            set_crc.finalize(),
            locators.as_ref().map_or(0, |part| part.crc),
        ];
        let pages = [
            sketch.table.page_count() as u64,
            locators
                .as_ref()
                .map_or(0, |part| part.table.page_count() as u64),
        ];
        let header = header.encode(&sections, &crcs, pages);
        let path = &self.path;
        let mut out = Output {
            file: BufWriter::with_capacity(1 << 20, File::create(path).at(path)?),
            written: 0,
        };
        out.put(0, &header).at(path)?;
        out.put(sections[0].0, &table).at(path)?;
        out.copy(sections[1].0, sketch).at(path)?;
        out.put(sections[2].0, &offsets).at(path)?;
        for set in color_sets {
            out.put(out.written, set.as_ref()).at(path)?;
        }
        if let Some(part) = locators {
            out.copy(sections[4].0, part).at(path)?;
        }
        let file = (out.file.into_inner())
            .map_err(|error| error.into_error())
            .at(path)?;
        file.sync_all().at(path)
    }
}

/// A finished part: its table, bytes and CRC-32, and its file ready to be read from the start.
struct FinishedPart {
    table: PartTable,
    bytes: u64,
    crc: u32,
    file: File,
}

/// The file being joined, with the bytes written so far.
struct Output {
    file: BufWriter<File>,
    written: u64,
}

impl Output {
    /// Writes zero bytes up to `at`, then `bytes`.
    fn put(&mut self, at: u64, bytes: &[u8]) -> std::io::Result<()> {
        let padding = vec![0u8; at.saturating_sub(self.written) as usize];
        self.file.write_all(&padding)?;
        self.file.write_all(bytes)?;
        self.written = at.max(self.written) + bytes.len() as u64;
        Ok(())
    }

    /// Writes zero bytes up to `at`, then the bytes of a finished part.
    fn copy(&mut self, at: u64, mut part: FinishedPart) -> std::io::Result<()> {
        self.put(at, &[])?;
        let copied = std::io::copy(&mut part.file, &mut self.file)?;
        self.written += copied;
        Ok(())
    }
}

/// The counts that the header holds besides the sections.
struct HeaderFields {
    shape: SeedFileShape,
    seeds: u64,
    locations: u64,
    color_sets: u64,
}

impl HeaderFields {
    /// The header of the file with these sections, CRC-32s and page rows per part.
    fn encode(
        &self,
        sections: &[(u64, u64); SECTIONS],
        crcs: &[u32; SECTIONS],
        pages: [u64; 2],
    ) -> Vec<u8> {
        let mut header = vec![0u8; HEADER_BYTES as usize];
        header[..16].copy_from_slice(&MAGIC);
        header[16..20].copy_from_slice(&VERSION.to_le_bytes());
        header[20] = self.shape.bucket_bits as u8;
        header[21] = SEED_KEY_BITS as u8;
        header[22] = self.shape.page_bits as u8;
        let counts = [
            self.shape.entries,
            self.seeds,
            self.color_sets,
            self.locations,
            pages[0],
            pages[1],
        ];
        let fields = counts
            .into_iter()
            .chain(sections.iter().flat_map(|&(offset, bytes)| [offset, bytes]));
        for (index, value) in fields.enumerate() {
            header[24 + 8 * index..32 + 8 * index].copy_from_slice(&value.to_le_bytes());
        }
        for (index, crc) in crcs.iter().enumerate() {
            header[152 + 4 * index..156 + 4 * index].copy_from_slice(&crc.to_le_bytes());
        }
        let header_crc = crc32fast::hash(&header[..252]);
        header[252..].copy_from_slice(&header_crc.to_le_bytes());
        header
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::buckets::{encode_locator_bucket, encode_sketch_bucket};
    use crate::index::seed_file::SeedFile;
    use crate::index::seed_file::{ColorSetId, encode_color_set};
    use crate::storage::ReadAt;
    use crate::test_support::rng;
    use rand::RngExt;
    use std::collections::BTreeMap;
    use std::sync::Arc;

    /// Random stored seeds by key: color set id and, per entry of the set, seed locations.
    type Seeds = BTreeMap<u64, (u32, Vec<Vec<u32>>)>;

    fn random_seeds(bucket_bits: u32) -> (Seeds, Vec<Vec<u32>>) {
        let mut draws = rng(0x5eed + u64::from(bucket_bits));
        let sets: Vec<Vec<u32>> = (0..300)
            .map(|_| {
                let mut set: Vec<u32> = (0..draws.random_range(1..21))
                    .map(|_| draws.random_range(0..5000))
                    .collect();
                set.sort_unstable();
                set.dedup();
                set
            })
            .collect();
        let mut seeds = Seeds::new();
        for _ in 0..20_000 {
            let key = draws.random::<u64>() >> (64 - SEED_KEY_BITS);
            let id = match draws.random_range(0..10) {
                0 => draws.random_range(0..300),
                _ => draws.random_range(0..4),
            };
            let locations = (sets[id as usize].iter())
                .map(|_| {
                    let mut values: Vec<u32> = (0..draws.random_range(1..4))
                        .map(|_| draws.random_range(0..1 << 20))
                        .collect();
                    values.sort_unstable();
                    values.dedup();
                    values
                })
                .collect();
            seeds.insert(key, (id as u32, locations));
        }
        (seeds, sets)
    }

    /// Writes the seeds into a file in `directory` and opens it from RAM.
    fn write(directory: &Path, seeds: &Seeds, sets: &[Vec<u32>], shape: SeedFileShape) -> SeedFile {
        let path = directory.join(format!("seeds-{}", shape.page_bits));
        let mut writer = SeedWriter::create(&path, directory, shape).unwrap();
        let check_bits = SEED_KEY_BITS - shape.bucket_bits;
        let bucket_of = |key: &u64| key >> check_bits;
        let all: Vec<_> = seeds.iter().collect();
        for group in all.chunk_by(|a, b| bucket_of(a.0) == bucket_of(b.0)) {
            let mut sketch = Vec::new();
            let checks: Vec<(u64, u32)> = (group.iter())
                .map(|(key, (id, _))| (**key & ((1 << check_bits) - 1), *id))
                .collect();
            encode_sketch_bucket(&checks, check_bits, &mut sketch);
            let seed_pairs: Vec<u32> = group.iter().map(|(_, (_, l))| l.len() as u32).collect();
            let lists = group.iter().flat_map(|(_, (_, lists))| lists);
            let pair_locations: Vec<u32> = lists.clone().map(|list| list.len() as u32).collect();
            let values: Vec<u32> = lists.flatten().copied().collect();
            let mut locators = Vec::new();
            encode_locator_bucket(&seed_pairs, &pair_locations, &values, &mut locators);
            let bucket = bucket_of(group[0].0) as usize;
            writer.push_bucket(bucket, &sketch, &locators).unwrap();
        }
        let encoded: Vec<Vec<u8>> = (sets.iter())
            .map(|set| {
                let mut bytes = Vec::new();
                encode_color_set(set, &mut bytes);
                bytes
            })
            .collect();
        writer.finish(&encoded).unwrap();
        SeedFile::open(Arc::new(std::fs::read(&path).unwrap())).unwrap()
    }

    #[test]
    fn seeds_read_back() {
        let directory = tempfile::tempdir().unwrap();
        for (bucket_bits, page_bits, locators) in [(10, 16, true), (12, 6, true), (8, 12, false)] {
            let (seeds, sets) = random_seeds(bucket_bits);
            let shape = SeedFileShape {
                bucket_bits,
                page_bits,
                entries: 5000,
                locators,
            };
            let file = write(directory.path(), &seeds, &sets, shape);
            file.verify_checksums().unwrap();
            assert_eq!(file.seed_count(), seeds.len() as u64);
            assert_eq!(file.has_locators(), locators);
            let mut hashes: Vec<u64> = (seeds.keys())
                .flat_map(|key| [key << 20 | 0xabcde, (key ^ 1) << 20])
                .collect();
            hashes.sort_unstable();
            let found = file.find_seeds(&hashes).unwrap();
            let ids = found.iter().flatten().map(|seed| seed.color_set);
            let read = file.read_color_sets(ids).unwrap();
            let mut wanted = Vec::new();
            for (hash, seed) in hashes.iter().zip(&found) {
                let Some((id, lists)) = seeds.get(&(hash >> 20)) else {
                    assert_eq!(*seed, None);
                    continue;
                };
                let seed = seed.unwrap();
                assert_eq!(seed.color_set, ColorSetId(*id));
                assert_eq!(read.entries(seed.color_set), Some(&sets[*id as usize][..]));
                wanted.push((seed, lists));
            }
            if !locators {
                assert!(file.read_locator_buckets(&[0], |_, _| Ok(())).is_err());
                continue;
            }
            let mut buckets: Vec<u32> = wanted.iter().map(|(seed, _)| seed.bucket).collect();
            buckets.dedup();
            let mut next = 0;
            file.read_locator_buckets(&buckets, |bucket, part| {
                while let Some((seed, lists)) = wanted.get(next).filter(|w| w.0.bucket == bucket) {
                    for (rank, list) in lists.iter().enumerate() {
                        let mut out = Vec::new();
                        part.locations(seed.place, rank as u32, &mut out).unwrap();
                        assert_eq!(&out, list);
                    }
                    next += 1;
                }
                Ok(())
            })
            .unwrap();
            assert_eq!(next, wanted.len());
        }
    }

    #[test]
    fn damage_found() {
        let directory = tempfile::tempdir().unwrap();
        let (seeds, sets) = random_seeds(9);
        let shape = SeedFileShape {
            bucket_bits: 9,
            page_bits: 10,
            entries: 5000,
            locators: true,
        };
        write(directory.path(), &seeds, &sets, shape);
        let bytes = std::fs::read(directory.path().join("seeds-10")).unwrap();
        for at in [20, 300, bytes.len() / 2, bytes.len() - 3] {
            let mut damaged = bytes.clone();
            damaged[at] ^= 0x10;
            let opened = SeedFile::open(Arc::new(damaged));
            assert!(opened.is_err() || opened.unwrap().verify_checksums().is_err());
        }
        let cut = bytes[..bytes.len() - 1].to_vec();
        assert!(SeedFile::open(Arc::new(cut)).is_err());
        assert!(Arc::new(bytes).size() > 0);
    }
}
