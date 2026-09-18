//! The packed sequence of one entry (stage 3): contigs at 2 bits per base in checked blocks, and
//! the contig names in checked zstd blocks, so a read fetches only the blocks it needs.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::IndexError;
use super::bits::{u32_at, u64_at};
use super::contig_table::{ContigTable, Run};
use crate::seeds::base_code;
use crate::storage::{ReadAt, will_need_ranges};
use std::collections::BTreeMap;
use std::ops::Range;
use std::sync::{Arc, OnceLock};

const MAGIC: [u8; 16] = *b"jam-rs 2bit seq\0";
const VERSION: u32 = 2;
/// Bytes of the header fields at the front of the header page.
const HEADER_BYTES: u64 = 64;
/// The packed bases start after a page that holds the header, so blocks line up with pages.
pub(crate) const BASES_AT: u64 = 4096;
/// Bytes of one stored block: packed bases, then their CRC-32.
pub(crate) const BLOCK_BYTES: u64 = 4096;
/// Packed bases of one block, so 16,368 bases.
const BLOCK_DATA: u64 = BLOCK_BYTES - 4;
/// Fewest contigs whose names form one name block.
const NAME_BLOCK: u64 = 256;
/// Name blocks whose ends fit the header page after the header fields.
const INDEX_SLOTS: u64 = (BASES_AT - HEADER_BYTES) / 8;
/// The zstd level of name blocks: they are written once and read rarely.
const NAME_LEVEL: i32 = 19;

/// Collects the contigs of one entry and encodes its packed sequence file.
#[derive(Default)]
pub struct PackedWriter {
    starts: Vec<u64>,
    /// Every contig name followed by a line feed.
    names: Vec<u8>,
    runs: Vec<Run>,
    packed: Vec<u8>,
    bases: u64,
}

impl PackedWriter {
    /// An empty entry.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds the next contig; `name` must be one word.
    pub fn add_contig(&mut self, name: &[u8], sequence: &[u8]) -> Result<(), IndexError> {
        if !valid_name(name) {
            return Err(IndexError::Input(format!(
                "contig name {:?} is empty or holds white space",
                String::from_utf8_lossy(name)
            )));
        }
        self.starts.push(self.bases);
        self.names.extend_from_slice(name);
        self.names.push(b'\n');
        self.packed.reserve((sequence.len() / 4).saturating_sub(1));
        for &byte in sequence {
            let slot = self.bases & 3;
            if slot == 0 {
                self.packed.push(0);
            }
            match base_code(byte) {
                None => self.push_run(byte.to_ascii_uppercase()),
                Some(code) => {
                    if let Some(last) = self.packed.last_mut() {
                        *last |= code << (2 * slot);
                    }
                }
            }
            self.bases += 1;
        }
        Ok(())
    }

    /// Bases added so far.
    pub fn bases(&self) -> u64 {
        self.bases
    }

    /// Contigs added so far.
    pub fn contig_count(&self) -> u64 {
        self.starts.len() as u64
    }

    /// The contig table of the contigs added so far.
    pub(crate) fn contig_table(&self) -> ContigTable {
        let ends = self.starts.iter().skip(1).chain([&self.bases]);
        let lengths = self.starts.iter().zip(ends).map(|(start, end)| end - start);
        ContigTable::new(lengths, self.runs.clone())
    }

    /// The file in parts, so a caller can write them without copying the bases.
    pub fn finish(self) -> Result<PackedFile, IndexError> {
        let mut index = Vec::new();
        let mut names = Vec::new();
        let lines: Vec<&[u8]> = self.names.split_inclusive(|b| *b == b'\n').collect();
        let block_contigs = name_block(lines.len() as u64) as usize;
        for block in lines.chunks(block_contigs) {
            let start = names.len();
            names.extend(zstd::bulk::compress(&block.concat(), NAME_LEVEL)?);
            let crc = crc32fast::hash(&names[start..]);
            names.extend_from_slice(&crc.to_le_bytes());
            index.extend_from_slice(&(names.len() as u64).to_le_bytes());
        }
        let mut header = vec![0u8; BASES_AT as usize];
        header[HEADER_BYTES as usize..][..index.len()].copy_from_slice(&index);
        header[..16].copy_from_slice(&MAGIC);
        header[16..20].copy_from_slice(&VERSION.to_le_bytes());
        let fields = [self.starts.len() as u64, self.bases, names.len() as u64];
        for (at, value) in fields.into_iter().enumerate() {
            header[24 + 8 * at..32 + 8 * at].copy_from_slice(&value.to_le_bytes());
        }
        let header_crc = page_crc(&header);
        header[48..52].copy_from_slice(&header_crc.to_le_bytes());
        let block_crcs = (self.packed.chunks(BLOCK_DATA as usize))
            .map(|block| crc32fast::hash(block).to_le_bytes())
            .collect();
        Ok(PackedFile {
            header,
            packed: self.packed,
            block_crcs,
            names,
        })
    }

    /// Adds one byte other than A, C, G or T, extending the last run when it continues it.
    fn push_run(&mut self, byte: u8) {
        if let Some(run) = self.runs.last_mut()
            && run.byte == byte
            && run.start + u64::from(run.length) == self.bases
            && run.length < u32::MAX
        {
            run.length += 1;
            return;
        }
        self.runs.push(Run {
            start: self.bases,
            length: 1,
            byte,
        });
    }
}

/// A packed sequence file in parts, written in the order of [`PackedFile::parts`].
pub struct PackedFile {
    header: Vec<u8>,
    packed: Vec<u8>,
    block_crcs: Vec<[u8; 4]>,
    /// The name blocks.
    names: Vec<u8>,
}

impl PackedFile {
    /// The bytes of the file as consecutive slices.
    pub fn parts(&self) -> Vec<&[u8]> {
        let blocks = self
            .packed
            .chunks(BLOCK_DATA as usize)
            .zip(&self.block_crcs);
        let mut parts = vec![&self.header[..]];
        parts.extend(blocks.flat_map(|(block, crc)| [block, &crc[..]]));
        parts.push(&self.names);
        parts
    }
}

/// The name and length of one contig of an entry.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Contig {
    pub name: String,
    pub length: u64,
}

/// An open packed sequence with the contig table of its entry. Opening only hints the header,
/// which the first read checks; names are read only for the contigs asked for.
pub struct PackedSequence {
    source: Arc<dyn ReadAt>,
    table: ContigTable,
    /// The end of each name block, set once the header page has been read and checked.
    name_ends: OnceLock<Vec<u64>>,
}

impl PackedSequence {
    /// Opens a packed sequence with the contig table of its entry, and hints its header.
    pub fn open(source: Arc<dyn ReadAt>, table: ContigTable) -> Self {
        source.will_need(0, BASES_AT);
        Self {
            source,
            table,
            name_ends: OnceLock::new(),
        }
    }

    /// Bases of all contigs.
    pub fn bases(&self) -> u64 {
        self.table.bases()
    }

    /// Contigs of the entry.
    pub fn contig_count(&self) -> usize {
        self.table.contig_count()
    }

    /// Entry coordinates of a contig.
    pub fn contig_range(&self, contig: usize) -> Result<Range<u64>, IndexError> {
        self.table.contig_range(contig)
    }

    /// The name of a contig, read with the names of its name block.
    pub fn contig_name(&self, contig: usize) -> Result<Vec<u8>, IndexError> {
        let index = u32::try_from(contig).map_err(|_| no_contig(contig))?;
        let mut names = self.read_names(vec![index])?;
        names.remove(&index).ok_or_else(|| no_contig(contig))
    }

    /// The names and lengths of the given contigs, by contig index.
    pub fn contigs(
        &self,
        contigs: impl IntoIterator<Item = u32>,
    ) -> Result<BTreeMap<u32, Contig>, IndexError> {
        let names = self.read_names(contigs.into_iter().collect())?;
        (names.into_iter())
            .map(|(contig, name)| {
                let range = self.contig_range(contig as usize)?;
                let name = String::from_utf8_lossy(&name).into_owned();
                let length = range.end - range.start;
                Ok((contig, Contig { name, length }))
            })
            .collect()
    }

    /// Asks the source for the stored blocks that hold these ranges of contigs (contig
    /// coordinates) at once, joining close ones, so later reads find them.
    pub fn will_need<'a>(&self, ranges: impl IntoIterator<Item = (usize, &'a Range<u64>)>) {
        let mut stored: Vec<(u64, u64)> = Vec::new();
        for (contig, range) in ranges {
            let Ok(contig_range) = self.contig_range(contig) else {
                continue;
            };
            let bases = contig_range.start + range.start..contig_range.start + range.end;
            stored.extend(self.stored_range(bases).map(|r| (r.start, r.end)));
        }
        will_need_ranges(self.source.as_ref(), &mut stored);
    }

    /// Replaces `out` with bases `range` of a contig (contig coordinates) as upper case ASCII.
    pub fn read(
        &self,
        contig: usize,
        range: Range<u64>,
        out: &mut Vec<u8>,
    ) -> Result<(), IndexError> {
        let contig_range = self.contig_range(contig)?;
        if range.start > range.end || range.end > contig_range.end - contig_range.start {
            return Err(IndexError::Input(format!(
                "bases {range:?} lie outside contig {contig}"
            )));
        }
        self.read_bases(
            contig_range.start + range.start..contig_range.start + range.end,
            out,
        )
    }

    /// Replaces `out` with bases `range` in entry coordinates. Every block read is checked
    /// against its CRC-32.
    pub fn read_bases(&self, range: Range<u64>, out: &mut Vec<u8>) -> Result<(), IndexError> {
        out.clear();
        if range.start > range.end || range.end > self.bases() {
            return Err(IndexError::Input(format!(
                "bases {range:?} lie outside the entry"
            )));
        }
        let Some(stored) = self.stored_range(range.clone()) else {
            return Ok(());
        };
        self.check_header()?;
        let mut bytes = vec![0u8; (stored.end - stored.start) as usize];
        self.source.read_exact_at(stored.start, &mut bytes)?;
        let mut packed = Vec::with_capacity(bytes.len());
        for block in bytes.chunks(BLOCK_BYTES as usize) {
            let (data, crc) = block.split_at(block.len().saturating_sub(4));
            if u32_at(crc, 0) != Some(crc32fast::hash(data)) {
                return Err(IndexError::Checksum("packed bases".into()));
            }
            packed.extend_from_slice(data);
        }
        let first_byte = (stored.start - BASES_AT) / BLOCK_BYTES * BLOCK_DATA;
        out.reserve((range.end - range.start) as usize);
        for position in range.clone() {
            let byte = packed[(position / 4 - first_byte) as usize];
            out.push(b"ACGT"[usize::from((byte >> (2 * (position % 4))) & 3)]);
        }
        for run in self.table.runs_in(range.clone()) {
            let from = run.start.max(range.start) - range.start;
            let to = (run.start + u64::from(run.length)).min(range.end) - range.start;
            out[from as usize..to as usize].fill(run.byte);
        }
        Ok(())
    }

    /// The file range of the stored blocks that hold bases `range` (entry coordinates), or
    /// `None` when it is empty.
    fn stored_range(&self, range: Range<u64>) -> Option<Range<u64>> {
        if range.is_empty() {
            return None;
        }
        let first = range.start / 4 / BLOCK_DATA;
        let last = (range.end - 1) / 4 / BLOCK_DATA;
        let end = (BASES_AT + (last + 1) * BLOCK_BYTES).min(self.names_at());
        Some(BASES_AT + first * BLOCK_BYTES..end)
    }

    /// Where the name blocks start, after the bases, from the contig table.
    fn names_at(&self) -> u64 {
        let data = self.bases().div_ceil(4);
        BASES_AT + data + 4 * data.div_ceil(BLOCK_DATA)
    }

    /// Reads and checks the header page once: its CRC-32, and that it fits the contig table and
    /// the file size. Returns the end of each name block within the name blocks.
    fn check_header(&self) -> Result<&[u64], IndexError> {
        if let Some(ends) = self.name_ends.get() {
            return Ok(ends);
        }
        let mut header = vec![0u8; BASES_AT as usize];
        self.source.read_exact_at(0, &mut header)?;
        if header[..16] != MAGIC || u32_at(&header, 16) != Some(VERSION) {
            return Err(IndexError::format("not a packed sequence of this version"));
        }
        if u32_at(&header, 48) != Some(page_crc(&header)) || header[20..24] != [0; 4] {
            return Err(IndexError::Checksum("packed sequence header".into()));
        }
        let [contigs, bases, names] = [24, 32, 40].map(|at| u64_at(&header, at).unwrap_or(0));
        let blocks = contigs.div_ceil(name_block(contigs)) as usize;
        let ends: Vec<u64> = (0..blocks)
            .map(|block| u64_at(&header, HEADER_BYTES as usize + 8 * block).unwrap_or(0))
            .collect();
        if contigs != self.contig_count() as u64
            || bases != self.bases()
            || self.names_at().checked_add(names) != Some(self.source.size())
            || !ends.is_sorted()
            || ends.last().is_some_and(|&end| end != names)
        {
            return Err(IndexError::format("packed sequence layout"));
        }
        Ok(self.name_ends.get_or_init(|| ends))
    }

    /// The names of the given contigs. Only the name blocks that hold them are read, at once;
    /// their bounds come from the header page.
    fn read_names(&self, mut wanted: Vec<u32>) -> Result<BTreeMap<u32, Vec<u8>>, IndexError> {
        wanted.sort_unstable();
        wanted.dedup();
        if let Some(&last) = wanted.last() {
            self.contig_range(last as usize)?;
        }
        let ends = self.check_header()?;
        let size = name_block(self.contig_count() as u64);
        let mut blocks: Vec<u64> = wanted.iter().map(|&c| u64::from(c) / size).collect();
        blocks.dedup();
        let names_at = self.names_at();
        let ranges = blocks.iter().map(|&block| {
            let start = (block.checked_sub(1)).map_or(0, |before| ends[before as usize]);
            names_at + start..names_at + ends[block as usize]
        });
        let read = self.read_ranges(ranges)?;
        let mut found = BTreeMap::new();
        for (&block, bytes) in blocks.iter().zip(&read) {
            let names = self.block_names(block, size, bytes)?;
            let inside = wanted.iter().filter(|&&c| u64::from(c) / size == block);
            for &contig in inside {
                let name = &names[(u64::from(contig) - block * size) as usize];
                found.insert(contig, name.clone());
            }
        }
        Ok(found)
    }

    /// Reads byte ranges, hinting all of them before the first read.
    fn read_ranges(
        &self,
        ranges: impl Iterator<Item = Range<u64>>,
    ) -> Result<Vec<Vec<u8>>, IndexError> {
        let ranges: Vec<Range<u64>> = ranges.collect();
        let mut hints: Vec<(u64, u64)> = ranges.iter().map(|r| (r.start, r.end)).collect();
        will_need_ranges(self.source.as_ref(), &mut hints);
        ranges
            .into_iter()
            .map(|range| {
                let mut bytes = vec![0u8; (range.end - range.start) as usize];
                self.source.read_exact_at(range.start, &mut bytes)?;
                Ok(bytes)
            })
            .collect()
    }

    /// Checks a name block of `size` contigs against its CRC-32, unpacks it and splits it into
    /// the names of its contigs.
    fn block_names(&self, block: u64, size: u64, bytes: &[u8]) -> Result<Vec<Vec<u8>>, IndexError> {
        let damaged = || IndexError::Checksum("packed sequence names".into());
        let (packed, crc) = bytes.split_at(bytes.len().checked_sub(4).ok_or_else(damaged)?);
        if u32_at(crc, 0) != Some(crc32fast::hash(packed)) {
            return Err(damaged());
        }
        let wrong = || IndexError::format("packed sequence names");
        let lines = zstd::stream::decode_all(packed).map_err(|_| wrong())?;
        let lines = lines.strip_suffix(b"\n").ok_or_else(wrong)?;
        let names: Vec<Vec<u8>> = lines.split(|b| *b == b'\n').map(<[u8]>::to_vec).collect();
        let count = (self.contig_count() as u64 - block * size).min(size);
        if names.len() as u64 != count || !names.iter().all(|name| valid_name(name)) {
            return Err(wrong());
        }
        Ok(names)
    }
}

/// Contigs per name block for an entry of `contigs` contigs: a power of two of at least
/// `NAME_BLOCK`, so the end of every block fits the header page.
fn name_block(contigs: u64) -> u64 {
    contigs
        .div_ceil(INDEX_SLOTS)
        .next_power_of_two()
        .max(NAME_BLOCK)
}

/// The CRC-32 of a header page without its CRC field.
fn page_crc(page: &[u8]) -> u32 {
    let mut hasher = crc32fast::Hasher::new();
    hasher.update(&page[..48]);
    hasher.update(&page[52..]);
    hasher.finalize()
}

/// A contig name is one word.
fn valid_name(name: &[u8]) -> bool {
    !name.is_empty() && !name.iter().any(u8::is_ascii_whitespace)
}

fn no_contig(contig: usize) -> IndexError {
    IndexError::Input(format!("contig {contig} is not in this entry"))
}

#[cfg(test)]
#[path = "packed_tests.rs"]
mod tests;
