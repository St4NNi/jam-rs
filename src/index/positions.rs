//! The positions file (stage 2): per entry its seed locations and contig starts, read memory-mapped
//! or as byte ranges.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::bits::{BIT_PADDING, BitWriter, bits_needed, read_bits, u32_at, u64_at};
use super::{AtPath, IndexError, seed_key};
use crate::storage::ReadAt;
use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// The file name of the positions in a shard directory.
pub const POSITIONS_FILE: &str = "positions";

/// Bits of a hash that the positions store per seed: enough for the largest bucket bits and the
/// fingerprint, so shards of any size can be joined without changing their sections.
pub const POSITION_KEY_BITS: u32 = super::MAX_BUCKET_BITS as u32 + super::FINGERPRINT_BITS as u32;

const MAGIC: [u8; 16] = *b"jam positions\0\0\0";
const VERSION: u32 = 2;
const HEADER_BYTES: u64 = 128;
const ROW_BYTES: usize = 40;
/// Seed locations per slot that the slot bits aim for.
const SLOT_LOCATIONS: u64 = 16;

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

/// Where one entry's section lies and how its records are packed.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct Section {
    /// Offset of the section in the file.
    pub(crate) offset: u64,
    /// Bytes of the section.
    pub(crate) length: u64,
    /// Seed locations in the section.
    pub(crate) locations: u64,
    /// CRC-32 of the section bytes.
    pub(crate) crc: u32,
    /// Contigs of the entry.
    pub(crate) contigs: u64,
    slot_bits: u8,
    suffix_bits: u8,
    location_bits: u8,
}

impl Section {
    fn record_bits(&self) -> u64 {
        u64::from(self.suffix_bits) + u64::from(self.location_bits)
    }

    fn directory_bytes(&self) -> u64 {
        4 * ((1u64 << self.slot_bits) + 1)
    }

    /// Bytes of the contig starts and the entry's bases, before the slot directory.
    fn contig_bytes(&self) -> u64 {
        8 * (self.contigs + 1)
    }

    fn to_row(self) -> [u8; ROW_BYTES] {
        let mut row = [0u8; ROW_BYTES];
        row[..8].copy_from_slice(&self.offset.to_le_bytes());
        row[8..16].copy_from_slice(&self.length.to_le_bytes());
        row[16..24].copy_from_slice(&self.locations.to_le_bytes());
        row[24..28].copy_from_slice(&self.crc.to_le_bytes());
        row[28..31].copy_from_slice(&[self.slot_bits, self.suffix_bits, self.location_bits]);
        row[32..40].copy_from_slice(&self.contigs.to_le_bytes());
        row
    }

    fn from_row(row: &[u8]) -> Option<Self> {
        let section = Self {
            offset: u64_at(row, 0)?,
            length: u64_at(row, 8)?,
            locations: u64_at(row, 16)?,
            crc: u32_at(row, 24)?,
            slot_bits: row[28],
            suffix_bits: row[29],
            location_bits: row[30],
            contigs: u64_at(row, 32)?,
        };
        let valid = u32::from(section.slot_bits) + u32::from(section.suffix_bits)
            == POSITION_KEY_BITS
            && section.location_bits <= 32
            && section.slot_bits <= 32
            && row[31] == 0
            && section.contigs < 1 << 40
            && section.length
                == section.contig_bytes()
                    + section.directory_bytes()
                    + (section.locations * section.record_bits()).div_ceil(8)
                    + BIT_PADDING as u64;
        valid.then_some(section)
    }
}

/// Encodes one entry's section from its (key, seed location value) pairs, sorted and distinct,
/// and its contig starts followed by its bases. The offset of the returned section is 0 until a
/// writer places it.
pub(crate) fn encode_section(
    records: &[(u64, u32)],
    contig_bounds: &[u64],
) -> Result<(Section, Vec<u8>), IndexError> {
    if contig_bounds.first() != Some(&0) || !contig_bounds.is_sorted() {
        return Err(IndexError::Input("contig starts must ascend from 0".into()));
    }
    let count = records.len() as u64;
    if u32::try_from(count).is_err() {
        return Err(IndexError::Input(
            "an entry has 2^32 seed locations or more".into(),
        ));
    }
    let slot_bits = if count >= 2 * SLOT_LOCATIONS {
        (count / SLOT_LOCATIONS).ilog2().min(POSITION_KEY_BITS)
    } else {
        0
    };
    let suffix_bits = POSITION_KEY_BITS - slot_bits;
    let location_bits = bits_needed(u64::from(records.iter().map(|r| r.1).max().unwrap_or(0)));
    let mut bytes: Vec<u8> = contig_bounds.iter().flat_map(|v| v.to_le_bytes()).collect();
    let mut next = 0usize;
    for slot in 0..=1u64 << slot_bits {
        while next < records.len() && records[next].0 >> suffix_bits < slot {
            next += 1;
        }
        bytes.extend_from_slice(&(next as u32).to_le_bytes());
    }
    let mut writer = BitWriter::new(&mut bytes);
    let suffix_mask = (1u64 << suffix_bits) - 1;
    for &(key, location) in records {
        writer.push(key & suffix_mask, suffix_bits);
        writer.push(u64::from(location), location_bits);
    }
    writer.finish();
    bytes.resize(bytes.len() + BIT_PADDING, 0);
    let section = Section {
        offset: 0,
        length: bytes.len() as u64,
        locations: count,
        crc: crc32fast::hash(&bytes),
        slot_bits: slot_bits as u8,
        suffix_bits: suffix_bits as u8,
        location_bits: location_bits as u8,
        contigs: contig_bounds.len() as u64 - 1,
    };
    Ok((section, bytes))
}

/// Decodes all (key, seed location value) pairs of a section from its bytes, in stored order.
pub(crate) fn decode_section(
    section: &Section,
    bytes: &[u8],
    mut visit: impl FnMut(u64, u32),
) -> Result<(), IndexError> {
    let damaged = || IndexError::format("positions section");
    let directory = bytes
        .get(section.contig_bytes() as usize..)
        .ok_or_else(damaged)?;
    let records = directory
        .get(section.directory_bytes() as usize..)
        .ok_or_else(damaged)?;
    let (suffix_bits, location_bits) = (
        u32::from(section.suffix_bits),
        u32::from(section.location_bits),
    );
    let mut index = 0u64;
    for slot in 0..1u64 << section.slot_bits {
        let end = u64::from(u32_at(directory, 4 * (slot as usize + 1)).ok_or_else(damaged)?);
        if end < index || end > section.locations {
            return Err(damaged());
        }
        while index < end {
            let bit = index * section.record_bits();
            let suffix = read_bits(records, bit, suffix_bits);
            let location = read_bits(records, bit + u64::from(suffix_bits), location_bits);
            visit((slot << suffix_bits) | suffix, location as u32);
            index += 1;
        }
    }
    Ok(())
}

/// The positions of a shard, read by offset from a local file or S3.
pub struct Positions {
    source: Arc<dyn ReadAt>,
    bin_bits: u32,
    sections: Vec<Section>,
}

impl Positions {
    /// Reads and checks the header and the section table.
    pub fn open(source: Arc<dyn ReadAt>) -> Result<Self, IndexError> {
        let mut header = [0u8; HEADER_BYTES as usize];
        source.read_exact_at(0, &mut header)?;
        if header[..16] != MAGIC || u32_at(&header, 16) != Some(VERSION) {
            return Err(IndexError::format("not a positions file of this version"));
        }
        if u32_at(&header, 124) != Some(crc32fast::hash(&header[..124])) {
            return Err(IndexError::Checksum("positions header".into()));
        }
        if u32::from(header[20]) != POSITION_KEY_BITS {
            return Err(IndexError::format("positions key bits"));
        }
        let bin_bits = u32::from(header[21]);
        let entries = u64_at(&header, 24).unwrap_or_default();
        let table_bytes = entries
            .checked_mul(ROW_BYTES as u64)
            .filter(|bytes| HEADER_BYTES + bytes <= source.size())
            .ok_or_else(|| IndexError::format("positions table is cut off"))?;
        let mut table = vec![0u8; table_bytes as usize];
        source.read_exact_at(HEADER_BYTES, &mut table)?;
        if u32_at(&header, 32) != Some(crc32fast::hash(&table)) {
            return Err(IndexError::Checksum("positions table".into()));
        }
        let mut sections = Vec::with_capacity(entries as usize);
        let mut end = HEADER_BYTES + table_bytes;
        for row in table.as_chunks::<ROW_BYTES>().0 {
            let section = Section::from_row(row)
                .filter(|section| section.offset >= end)
                .ok_or_else(|| IndexError::format("positions section table"))?;
            end = section.offset + section.length;
            sections.push(section);
        }
        if end > source.size() {
            return Err(IndexError::format("positions file is cut off"));
        }
        Ok(Self {
            source,
            bin_bits,
            sections,
        })
    }

    /// Appends to `out` the seed locations of the stored seed with this hash's key in one entry.
    pub fn locations(
        &self,
        entry: u32,
        hash: u64,
        out: &mut Vec<SeedLocation>,
    ) -> Result<(), IndexError> {
        let section = self.section(entry)?;
        let key = seed_key(hash, POSITION_KEY_BITS);
        let suffix_bits = u32::from(section.suffix_bits);
        let slot = key >> suffix_bits;
        let directory = section.offset + section.contig_bytes();
        let mut slot_range = [0u8; 8];
        self.source
            .read_exact_at(directory + 4 * slot, &mut slot_range)?;
        let first = u64::from(u32_at(&slot_range, 0).unwrap_or_default());
        let end = u64::from(u32_at(&slot_range, 4).unwrap_or_default());
        if first >= end {
            return Ok(());
        }
        if end > section.locations {
            return Err(IndexError::format("positions slot directory"));
        }
        let record_bits = section.record_bits();
        let start_byte = first * record_bits / 8;
        let end_byte = (end * record_bits).div_ceil(8) + BIT_PADDING as u64;
        let mut records = vec![0u8; (end_byte - start_byte) as usize];
        self.source.read_exact_at(
            directory + section.directory_bytes() + start_byte,
            &mut records,
        )?;
        let suffix = key & ((1u64 << suffix_bits) - 1);
        let location_bits = u32::from(section.location_bits);
        for index in first..end {
            let bit = index * record_bits - start_byte * 8;
            if read_bits(&records, bit, suffix_bits) == suffix {
                let value = read_bits(&records, bit + u64::from(suffix_bits), location_bits);
                out.push(SeedLocation::from_value(value as u32));
            }
        }
        Ok(())
    }

    /// The start of each contig of an entry in entry coordinates, followed by the entry's bases,
    /// so contig `i` spans `bounds[i]..bounds[i + 1]`.
    pub fn contig_starts(&self, entry: u32) -> Result<Vec<u64>, IndexError> {
        let section = self.section(entry)?;
        let mut bytes = vec![0u8; section.contig_bytes() as usize];
        self.source.read_exact_at(section.offset, &mut bytes)?;
        let bounds: Vec<u64> = bytes
            .as_chunks::<8>()
            .0
            .iter()
            .map(|v| u64::from_le_bytes(*v))
            .collect();
        if bounds.first() != Some(&0) || !bounds.is_sorted() {
            return Err(IndexError::format("positions contig starts"));
        }
        Ok(bounds)
    }

    /// The section of an entry.
    pub(crate) fn section(&self, entry: u32) -> Result<&Section, IndexError> {
        self.sections
            .get(entry as usize)
            .ok_or_else(|| IndexError::Input(format!("entry {entry} is not in this shard")))
    }

    /// The bytes of an entry's section, checked against its CRC-32.
    pub(crate) fn section_bytes(&self, entry: u32) -> Result<Vec<u8>, IndexError> {
        let section = self.section(entry)?;
        let mut bytes = vec![0u8; section.length as usize];
        self.source.read_exact_at(section.offset, &mut bytes)?;
        if crc32fast::hash(&bytes) != section.crc {
            return Err(IndexError::Checksum(format!("positions of entry {entry}")));
        }
        Ok(bytes)
    }

    /// Entries in the file.
    pub fn entry_count(&self) -> u32 {
        self.sections.len() as u32
    }

    /// Bases per bin, as a power of two.
    pub fn bin_bits(&self) -> u32 {
        self.bin_bits
    }

    /// Seed locations of one entry.
    pub fn entry_locations(&self, entry: u32) -> Result<u64, IndexError> {
        Ok(self.section(entry)?.locations)
    }

    /// Checks the CRC-32 of every section.
    pub fn verify_checksums(&self) -> Result<(), IndexError> {
        (0..self.entry_count()).try_for_each(|entry| self.section_bytes(entry).map(drop))
    }
}

/// Writes a positions file section by section, in entry order.
pub(crate) struct PositionsWriter {
    file: BufWriter<File>,
    path: PathBuf,
    bin_bits: u32,
    entries: usize,
    sections: Vec<Section>,
    written: u64,
}

impl PositionsWriter {
    /// Creates the file at `path` for `entries` sections.
    pub(crate) fn create(path: &Path, bin_bits: u32, entries: usize) -> Result<Self, IndexError> {
        let mut file = File::create_new(path).at(path)?;
        let start = HEADER_BYTES + (entries * ROW_BYTES) as u64;
        file.seek(SeekFrom::Start(start)).at(path)?;
        Ok(Self {
            file: BufWriter::with_capacity(1 << 20, file),
            path: path.to_path_buf(),
            bin_bits,
            entries,
            sections: Vec::with_capacity(entries),
            written: start,
        })
    }

    /// Appends the next entry's section, as made by `encode_section` or read from another file.
    pub(crate) fn push_section(
        &mut self,
        section: &Section,
        bytes: &[u8],
    ) -> Result<(), IndexError> {
        if self.sections.len() == self.entries
            || bytes.len() as u64 != section.length
            || u32::from(section.slot_bits) + u32::from(section.suffix_bits) != POSITION_KEY_BITS
        {
            return Err(IndexError::format(
                "positions section does not fit the file",
            ));
        }
        let padding = self.written.next_multiple_of(8) - self.written;
        let path = &self.path;
        self.file.write_all(&vec![0u8; padding as usize]).at(path)?;
        self.written += padding;
        self.sections.push(Section {
            offset: self.written,
            ..*section
        });
        self.file.write_all(bytes).at(path)?;
        self.written += section.length;
        Ok(())
    }

    /// Writes the header and the section table.
    pub(crate) fn finish(self) -> Result<(), IndexError> {
        if self.sections.len() != self.entries {
            return Err(IndexError::format("positions file is missing sections"));
        }
        let table: Vec<u8> = self.sections.iter().flat_map(|s| s.to_row()).collect();
        let mut header = [0u8; HEADER_BYTES as usize];
        header[..16].copy_from_slice(&MAGIC);
        header[16..20].copy_from_slice(&VERSION.to_le_bytes());
        header[20] = POSITION_KEY_BITS as u8;
        header[21] = self.bin_bits as u8;
        header[24..32].copy_from_slice(&(self.entries as u64).to_le_bytes());
        header[32..36].copy_from_slice(&crc32fast::hash(&table).to_le_bytes());
        let header_crc = crc32fast::hash(&header[..124]);
        header[124..].copy_from_slice(&header_crc.to_le_bytes());
        let path = &self.path;
        let mut file = self
            .file
            .into_inner()
            .map_err(|error| error.into_error())
            .at(path)?;
        file.seek(SeekFrom::Start(0)).at(path)?;
        file.write_all(&header).at(path)?;
        file.write_all(&table).at(path)?;
        file.sync_all().at(path)?;
        Ok(())
    }
}

#[cfg(test)]
#[path = "positions_tests.rs"]
mod tests;
