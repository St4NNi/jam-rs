//! The positions file (stage 2): for every entry a section with the postings of its seeds, so a
//! seed of one candidate entry is found with two small reads. It works memory-mapped from local
//! disk and as byte ranges from S3.
//!
//! A seed is stored by its key: the top `key_bits` bits of its hash (bucket bits and fingerprint,
//! as in the sketch). A posting is the bin within the entry, shifted left by one, with the low bit
//! set when the canonical k-mer is the reverse complement of the contig strand.
//!
//! Layout, little endian:
//!
//! | Bytes | Field |
//! | --- | --- |
//! | 0..16 | magic `jam positions` padded with zero bytes |
//! | 16..20 | version (u32, 1) |
//! | 20 | key bits (u8) |
//! | 21 | bin bits (u8) |
//! | 22..24 | zero |
//! | 24..32 | entries e (u64) |
//! | 32..36 | CRC-32 of the section table (u32) |
//! | 36..124 | zero |
//! | 124..128 | CRC-32 of bytes 0..124 (u32) |
//! | 128.. | section table: e rows of 32 bytes, see below |
//! | sections | one per entry, in entry order, each at a multiple of 8 |
//!
//! A section table row:
//!
//! | Bytes | Field |
//! | --- | --- |
//! | 0..8 | section offset in the file (u64) |
//! | 8..16 | section length (u64) |
//! | 16..24 | postings (u64) |
//! | 24..28 | CRC-32 of the section (u32) |
//! | 28 | slot bits d (u8) |
//! | 29 | suffix bits s = key bits - d (u8) |
//! | 30 | posting bits p (u8) |
//! | 31 | zero |
//!
//! A section with r postings, sorted by key and then posting:
//!
//! | Bytes | Field |
//! | --- | --- |
//! | 4 (2^d + 1) | slot directory: index of the first posting of each slot (u32); a key's slot is its top d bits |
//! | ceil(r (s + p) / 8) | one record per posting, bit-packed, first record in the low bits: the low s bits of the key, then the posting in p bits |
//! | 8 | zero |
//!
//! d is chosen so a slot holds about 16 postings; a lookup reads 8 bytes of the directory and the
//! records of one slot. A seed that occurs at several bins has one record per bin.

use super::bits::{BIT_PADDING, BitWriter, bits_needed, crc32, read_bits, u32_at, u64_at};
use super::{IndexError, seed_key};
use crate::storage::ReadAt;
use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::Arc;

/// The file name of the positions in a shard directory.
pub const POSITIONS_FILE: &str = "positions";

const MAGIC: [u8; 16] = *b"jam positions\0\0\0";
const VERSION: u32 = 1;
const HEADER_BYTES: u64 = 128;
const ROW_BYTES: usize = 32;
/// Postings per slot that the slot bits aim for.
const SLOT_POSTINGS: u64 = 16;

/// One stored bin and strand of a seed occurrence.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, PartialOrd, Ord)]
pub struct Posting {
    /// Bin within the entry: the seed starts in bases `bin << bin_bits` onward.
    pub bin: u32,
    /// The canonical k-mer is the reverse complement of the contig strand.
    pub reverse: bool,
}

impl Posting {
    /// The stored value: bin shifted left by one, strand in the low bit.
    pub fn value(self) -> u32 {
        (self.bin << 1) | u32::from(self.reverse)
    }

    /// The posting of a stored value.
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
    /// Postings in the section.
    pub(crate) postings: u64,
    /// CRC-32 of the section bytes.
    pub(crate) crc: u32,
    slot_bits: u8,
    suffix_bits: u8,
    posting_bits: u8,
}

impl Section {
    fn record_bits(&self) -> u64 {
        u64::from(self.suffix_bits) + u64::from(self.posting_bits)
    }

    fn directory_bytes(&self) -> u64 {
        4 * ((1u64 << self.slot_bits) + 1)
    }

    fn to_row(self) -> [u8; ROW_BYTES] {
        let mut row = [0u8; ROW_BYTES];
        row[..8].copy_from_slice(&self.offset.to_le_bytes());
        row[8..16].copy_from_slice(&self.length.to_le_bytes());
        row[16..24].copy_from_slice(&self.postings.to_le_bytes());
        row[24..28].copy_from_slice(&self.crc.to_le_bytes());
        row[28..31].copy_from_slice(&[self.slot_bits, self.suffix_bits, self.posting_bits]);
        row
    }

    fn from_row(row: &[u8], key_bits: u32) -> Option<Self> {
        let section = Self {
            offset: u64_at(row, 0)?,
            length: u64_at(row, 8)?,
            postings: u64_at(row, 16)?,
            crc: u32_at(row, 24)?,
            slot_bits: row[28],
            suffix_bits: row[29],
            posting_bits: row[30],
        };
        let valid = u32::from(section.slot_bits) + u32::from(section.suffix_bits) == key_bits
            && section.posting_bits <= 32
            && section.slot_bits <= 32
            && row[31] == 0
            && section.length
                == section.directory_bytes()
                    + (section.postings * section.record_bits()).div_ceil(8)
                    + BIT_PADDING as u64;
        valid.then_some(section)
    }
}

/// Encodes one entry's section from its (key, posting value) pairs, sorted and distinct. The
/// offset of the returned section is 0 until a writer places it.
pub(crate) fn encode_section(
    records: &[(u64, u32)],
    key_bits: u32,
) -> Result<(Section, Vec<u8>), IndexError> {
    let count = records.len() as u64;
    if u32::try_from(count).is_err() {
        return Err(IndexError::Input(
            "an entry has 2^32 postings or more".into(),
        ));
    }
    let slot_bits = if count >= 2 * SLOT_POSTINGS {
        (count / SLOT_POSTINGS).ilog2().min(key_bits)
    } else {
        0
    };
    let suffix_bits = key_bits - slot_bits;
    let posting_bits = bits_needed(u64::from(records.iter().map(|r| r.1).max().unwrap_or(0)));
    let mut bytes = Vec::new();
    let mut next = 0usize;
    for slot in 0..=1u64 << slot_bits {
        while next < records.len() && records[next].0 >> suffix_bits < slot {
            next += 1;
        }
        bytes.extend_from_slice(&(next as u32).to_le_bytes());
    }
    let mut writer = BitWriter::new(&mut bytes);
    let suffix_mask = (1u64 << suffix_bits) - 1;
    for &(key, posting) in records {
        writer.push(key & suffix_mask, suffix_bits);
        writer.push(u64::from(posting), posting_bits);
    }
    writer.finish();
    bytes.resize(bytes.len() + BIT_PADDING, 0);
    let section = Section {
        offset: 0,
        length: bytes.len() as u64,
        postings: count,
        crc: crc32(&bytes),
        slot_bits: slot_bits as u8,
        suffix_bits: suffix_bits as u8,
        posting_bits: posting_bits as u8,
    };
    Ok((section, bytes))
}

/// Decodes all (key, posting value) pairs of a section from its bytes, in stored order.
pub(crate) fn decode_section(
    section: &Section,
    bytes: &[u8],
    mut visit: impl FnMut(u64, u32),
) -> Result<(), IndexError> {
    let damaged = || IndexError::format("positions section");
    let records = bytes
        .get(section.directory_bytes() as usize..)
        .ok_or_else(damaged)?;
    let (suffix_bits, posting_bits) = (
        u32::from(section.suffix_bits),
        u32::from(section.posting_bits),
    );
    let mut index = 0u64;
    for slot in 0..1u64 << section.slot_bits {
        let end = u64::from(u32_at(bytes, 4 * (slot as usize + 1)).ok_or_else(damaged)?);
        if end < index || end > section.postings {
            return Err(damaged());
        }
        while index < end {
            let bit = index * section.record_bits();
            let suffix = read_bits(records, bit, suffix_bits);
            let posting = read_bits(records, bit + u64::from(suffix_bits), posting_bits);
            visit((slot << suffix_bits) | suffix, posting as u32);
            index += 1;
        }
    }
    Ok(())
}

/// The positions of a shard, read by offset from a local file or S3.
pub struct Positions {
    source: Arc<dyn ReadAt>,
    key_bits: u32,
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
        if u32_at(&header, 124) != Some(crc32(&header[..124])) {
            return Err(IndexError::Checksum("positions header".into()));
        }
        let (key_bits, bin_bits) = (u32::from(header[20]), u32::from(header[21]));
        let entries = u64_at(&header, 24).unwrap_or_default();
        let table_bytes = entries
            .checked_mul(ROW_BYTES as u64)
            .filter(|bytes| HEADER_BYTES + bytes <= source.size())
            .ok_or_else(|| IndexError::format("positions table is cut off"))?;
        let mut table = vec![0u8; table_bytes as usize];
        source.read_exact_at(HEADER_BYTES, &mut table)?;
        if u32_at(&header, 32) != Some(crc32(&table)) {
            return Err(IndexError::Checksum("positions table".into()));
        }
        let mut sections = Vec::with_capacity(entries as usize);
        let mut end = HEADER_BYTES + table_bytes;
        for row in table.as_chunks::<ROW_BYTES>().0 {
            let section = Section::from_row(row, key_bits)
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
            key_bits,
            bin_bits,
            sections,
        })
    }

    /// Appends to `out` the postings of the stored seed with this hash's key in one entry.
    pub fn postings(
        &self,
        entry: u32,
        hash: u64,
        out: &mut Vec<Posting>,
    ) -> Result<(), IndexError> {
        let section = self.section(entry)?;
        let key = seed_key(hash, self.key_bits);
        let suffix_bits = u32::from(section.suffix_bits);
        let slot = key >> suffix_bits;
        let mut slot_range = [0u8; 8];
        self.source
            .read_exact_at(section.offset + 4 * slot, &mut slot_range)?;
        let first = u64::from(u32_at(&slot_range, 0).unwrap_or_default());
        let end = u64::from(u32_at(&slot_range, 4).unwrap_or_default());
        if first >= end {
            return Ok(());
        }
        if end > section.postings {
            return Err(IndexError::format("positions slot directory"));
        }
        let record_bits = section.record_bits();
        let start_byte = first * record_bits / 8;
        let end_byte = (end * record_bits).div_ceil(8) + BIT_PADDING as u64;
        let mut records = vec![0u8; (end_byte - start_byte) as usize];
        self.source.read_exact_at(
            section.offset + section.directory_bytes() + start_byte,
            &mut records,
        )?;
        let suffix = key & ((1u64 << suffix_bits) - 1);
        let posting_bits = u32::from(section.posting_bits);
        for index in first..end {
            let bit = index * record_bits - start_byte * 8;
            if read_bits(&records, bit, suffix_bits) == suffix {
                let value = read_bits(&records, bit + u64::from(suffix_bits), posting_bits);
                out.push(Posting::from_value(value as u32));
            }
        }
        Ok(())
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
        if crc32(&bytes) != section.crc {
            return Err(IndexError::Checksum(format!("positions of entry {entry}")));
        }
        Ok(bytes)
    }

    /// Entries in the file.
    pub fn entries(&self) -> u32 {
        self.sections.len() as u32
    }

    /// Bits of a hash stored as a seed key.
    pub fn key_bits(&self) -> u32 {
        self.key_bits
    }

    /// Bases per bin, as a power of two.
    pub fn bin_bits(&self) -> u32 {
        self.bin_bits
    }

    /// Postings of one entry.
    pub fn entry_postings(&self, entry: u32) -> Result<u64, IndexError> {
        Ok(self.section(entry)?.postings)
    }

    /// Checks the CRC-32 of every section.
    pub fn verify_checksums(&self) -> Result<(), IndexError> {
        (0..self.entries()).try_for_each(|entry| self.section_bytes(entry).map(drop))
    }
}

/// Writes a positions file section by section, in entry order.
pub(crate) struct PositionsWriter {
    file: BufWriter<File>,
    key_bits: u32,
    bin_bits: u32,
    entries: usize,
    sections: Vec<Section>,
    written: u64,
}

impl PositionsWriter {
    /// Creates the file at `path` for `entries` sections.
    pub(crate) fn create(
        path: &Path,
        key_bits: u32,
        bin_bits: u32,
        entries: usize,
    ) -> Result<Self, IndexError> {
        let mut file = File::create_new(path)?;
        let start = HEADER_BYTES + (entries * ROW_BYTES) as u64;
        file.seek(SeekFrom::Start(start))?;
        Ok(Self {
            file: BufWriter::with_capacity(1 << 20, file),
            key_bits,
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
            || u32::from(section.slot_bits) + u32::from(section.suffix_bits) != self.key_bits
        {
            return Err(IndexError::format(
                "positions section does not fit the file",
            ));
        }
        let padding = self.written.next_multiple_of(8) - self.written;
        self.file.write_all(&vec![0u8; padding as usize])?;
        self.written += padding;
        self.sections.push(Section {
            offset: self.written,
            ..*section
        });
        self.file.write_all(bytes)?;
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
        header[20] = self.key_bits as u8;
        header[21] = self.bin_bits as u8;
        header[24..32].copy_from_slice(&(self.entries as u64).to_le_bytes());
        header[32..36].copy_from_slice(&crc32(&table).to_le_bytes());
        let header_crc = crc32(&header[..124]);
        header[124..].copy_from_slice(&header_crc.to_le_bytes());
        let mut file = self.file.into_inner().map_err(|error| error.into_error())?;
        file.seek(SeekFrom::Start(0))?;
        file.write_all(&header)?;
        file.write_all(&table)?;
        file.sync_all()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests;
