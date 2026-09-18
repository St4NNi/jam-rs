//! The packed sequence of one entry (stage 3): contigs at 2 bits per base with their names. The
//! contig starts and runs of other bytes lie in the contig table of the seeds file.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::IndexError;
use super::bits::{u32_at, u64_at};
use super::contig_table::{ContigTable, Run};
use crate::seeds::base_code;
use crate::storage::ReadAt;
use std::collections::BTreeMap;
use std::ops::Range;
use std::sync::Arc;

const MAGIC: [u8; 16] = *b"jam-rs 2bit seq\0";
const VERSION: u32 = 2;
const HEADER_BYTES: u64 = 128;
/// Packed bytes per checksum block, so 16,384 bases.
const BLOCK_BYTES: u64 = 4096;

/// Collects the contigs of one entry and encodes its packed sequence file.
#[derive(Default)]
pub struct PackedWriter {
    starts: Vec<u64>,
    name_ends: Vec<u64>,
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
        self.name_ends.push(self.names.len() as u64);
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
        self.name_ends.len() as u64
    }

    /// The contig table of the contigs added so far.
    pub(crate) fn contig_table(&self) -> ContigTable {
        let ends = self.starts.iter().skip(1).chain([&self.bases]);
        let lengths = self.starts.iter().zip(ends).map(|(start, end)| end - start);
        ContigTable::new(lengths, self.runs.clone())
    }

    /// The file in two parts, so a caller can write both without copying the bases.
    pub fn finish(self) -> PackedFile {
        let contigs = self.name_ends.len() as u64;
        let mut file = vec![0u8; HEADER_BYTES as usize];
        for value in &self.name_ends {
            file.extend_from_slice(&value.to_le_bytes());
        }
        file.extend_from_slice(&self.names);
        file.resize(file.len().next_multiple_of(8), 0);
        for block in self.packed.chunks(BLOCK_BYTES as usize) {
            file.extend_from_slice(&crc32fast::hash(block).to_le_bytes());
        }
        let packed_offset = file.len() as u64;
        let tables_crc = crc32fast::hash(&file[HEADER_BYTES as usize..]);
        file[..16].copy_from_slice(&MAGIC);
        file[16..20].copy_from_slice(&VERSION.to_le_bytes());
        let fields = [contigs, self.bases, self.names.len() as u64, packed_offset];
        for (index, value) in fields.into_iter().enumerate() {
            file[24 + 8 * index..32 + 8 * index].copy_from_slice(&value.to_le_bytes());
        }
        file[64..68].copy_from_slice(&tables_crc.to_le_bytes());
        let header_crc = crc32fast::hash(&file[..68]);
        file[68..72].copy_from_slice(&header_crc.to_le_bytes());
        PackedFile {
            header_and_tables: file,
            packed_bases: self.packed,
        }
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

/// A packed sequence file in two parts, one written after the other.
pub struct PackedFile {
    /// The header and the tables.
    pub header_and_tables: Vec<u8>,
    /// The packed bases.
    pub packed_bases: Vec<u8>,
}

/// The name and length of one contig of an entry.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Contig {
    pub name: String,
    pub length: u64,
}

/// An open packed sequence with the contig table of its entry. Block checksums stay in memory;
/// names and bases are read on request, so an entry that is never reported never reads names.
pub struct PackedSequence {
    source: Arc<dyn ReadAt>,
    table: ContigTable,
    block_crcs: Vec<u32>,
    packed_offset: u64,
    /// The byte count of the names without padding.
    names_bytes: u64,
    /// The CRC-32 of all tables and the running CRC-32 of the block checksums after the names,
    /// so the names complete the check of the tables when they are read.
    tables_crc: u32,
    after_names: crc32fast::Hasher,
}

/// The contig names of a packed sequence, read and checked together.
struct ContigNames {
    ends: Vec<u64>,
    bytes: Vec<u8>,
}

impl ContigNames {
    fn name(&self, contig: usize) -> Result<&[u8], IndexError> {
        let end = *self.ends.get(contig).ok_or_else(|| no_contig(contig))?;
        let start = contig
            .checked_sub(1)
            .map_or(0, |previous| self.ends[previous]);
        self.bytes
            .get(start as usize..end as usize)
            .filter(|name| valid_name(name))
            .ok_or_else(|| IndexError::format("contig name"))
    }
}

impl PackedSequence {
    /// Reads and checks the header and the block checksums, against the contig table of the
    /// entry.
    pub fn open(source: Arc<dyn ReadAt>, table: ContigTable) -> Result<Self, IndexError> {
        let mut header = [0u8; HEADER_BYTES as usize];
        source.read_exact_at(0, &mut header)?;
        let Header {
            contigs,
            bases,
            names,
            packed_offset,
        } = Header::decode(&header)?;
        let blocks = bases.div_ceil(4).div_ceil(BLOCK_BYTES);
        let tail_at = [contigs.checked_mul(8), names.checked_next_multiple_of(8)]
            .into_iter()
            .try_fold(HEADER_BYTES, |total, part| total.checked_add(part?));
        let tables_end = tail_at.and_then(|at| at.checked_add(blocks.checked_mul(4)?));
        if tables_end != Some(packed_offset)
            || packed_offset.checked_add(bases.div_ceil(4)) != Some(source.size())
            || contigs != table.contig_count() as u64
            || bases != table.bases()
        {
            return Err(IndexError::format("packed sequence layout"));
        }
        let tail_at = tail_at.unwrap_or_default();
        let mut tail = vec![0u8; (packed_offset - tail_at) as usize];
        source.read_exact_at(tail_at, &mut tail)?;
        let block_crcs = (0..tail.len() / 4)
            .map(|i| u32_at(&tail, 4 * i).unwrap_or_default())
            .collect();
        Ok(Self {
            source,
            table,
            block_crcs,
            packed_offset,
            names_bytes: names,
            tables_crc: u32_at(&header, 64).unwrap_or_default(),
            after_names: hasher_of(&tail),
        })
    }

    /// Bases of all contigs.
    pub fn bases(&self) -> u64 {
        self.table.bases()
    }

    /// Contigs of the entry.
    pub fn contig_count(&self) -> usize {
        self.table.contig_count()
    }

    /// The name of a contig, read with all names of the entry.
    pub fn contig_name(&self, contig: usize) -> Result<Vec<u8>, IndexError> {
        Ok(self.names()?.name(contig)?.to_vec())
    }

    /// The names and lengths of the given contigs, by contig index; all names are read once.
    pub fn contigs(
        &self,
        contigs: impl IntoIterator<Item = u32>,
    ) -> Result<BTreeMap<u32, Contig>, IndexError> {
        let names = self.names()?;
        contigs
            .into_iter()
            .map(|contig| {
                let name = names.name(contig as usize)?;
                let range = self.contig_range(contig as usize)?;
                let contig_record = Contig {
                    name: String::from_utf8_lossy(name).into_owned(),
                    length: range.end - range.start,
                };
                Ok((contig, contig_record))
            })
            .collect()
    }

    /// Reads the name ends and the names, and checks the CRC-32 of all tables with them.
    fn names(&self) -> Result<ContigNames, IndexError> {
        let contigs = self.contig_count();
        let mut bytes = vec![0u8; 8 * contigs + self.names_bytes.next_multiple_of(8) as usize];
        self.source.read_exact_at(HEADER_BYTES, &mut bytes)?;
        let mut crc = crc32fast::Hasher::new();
        crc.update(&bytes);
        crc.combine(&self.after_names);
        if crc.finalize() != self.tables_crc {
            return Err(IndexError::Checksum("packed sequence tables".into()));
        }
        let ends: Vec<u64> = (0..contigs)
            .map(|i| u64_at(&bytes, 8 * i).unwrap_or_default())
            .collect();
        if !ends.is_sorted() || ends.last().copied().unwrap_or(0) != self.names_bytes {
            return Err(IndexError::format("packed sequence tables"));
        }
        let names = 8 * contigs..8 * contigs + self.names_bytes as usize;
        Ok(ContigNames {
            ends,
            bytes: bytes[names].to_vec(),
        })
    }

    /// Asks the source for the packed bytes of bases `range` of a contig (contig coordinates),
    /// whole checksum blocks, so that a later read finds them. It never changes what is read.
    pub fn will_need(&self, contig: usize, range: Range<u64>) {
        let Ok(contig_range) = self.contig_range(contig) else {
            return;
        };
        let start = (contig_range.start + range.start) / 4 / BLOCK_BYTES * BLOCK_BYTES;
        let end = (contig_range.start + range.end)
            .div_ceil(4)
            .next_multiple_of(BLOCK_BYTES);
        let end = end.min(self.bases().div_ceil(4));
        if start < end {
            self.source
                .will_need(self.packed_offset + start, end - start);
        }
    }

    /// Entry coordinates of a contig.
    pub fn contig_range(&self, contig: usize) -> Result<Range<u64>, IndexError> {
        self.table.contig_range(contig)
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

    /// Replaces `out` with bases `range` in entry coordinates. Every packed block read is checked
    /// against its CRC-32.
    pub fn read_bases(&self, range: Range<u64>, out: &mut Vec<u8>) -> Result<(), IndexError> {
        out.clear();
        if range.start > range.end || range.end > self.bases() {
            return Err(IndexError::Input(format!(
                "bases {range:?} lie outside the entry"
            )));
        }
        if range.is_empty() {
            return Ok(());
        }
        let first_block = range.start / 4 / BLOCK_BYTES;
        let read_start = first_block * BLOCK_BYTES;
        let read_end =
            (((range.end - 1) / 4 / BLOCK_BYTES + 1) * BLOCK_BYTES).min(self.bases().div_ceil(4));
        let mut packed = vec![0u8; (read_end - read_start) as usize];
        self.source
            .read_exact_at(self.packed_offset + read_start, &mut packed)?;
        for (index, block) in packed.chunks(BLOCK_BYTES as usize).enumerate() {
            if self.block_crcs.get(first_block as usize + index) != Some(&crc32fast::hash(block)) {
                return Err(IndexError::Checksum("packed bases".into()));
            }
        }
        out.reserve((range.end - range.start) as usize);
        for position in range.clone() {
            let byte = packed[(position / 4 - read_start) as usize];
            out.push(b"ACGT"[usize::from((byte >> (2 * (position % 4))) & 3)]);
        }
        for run in self.table.runs_in(range.clone()) {
            let from = run.start.max(range.start) - range.start;
            let to = (run.start + u64::from(run.length)).min(range.end) - range.start;
            out[from as usize..to as usize].fill(run.byte);
        }
        Ok(())
    }
}

fn hasher_of(bytes: &[u8]) -> crc32fast::Hasher {
    let mut hasher = crc32fast::Hasher::new();
    hasher.update(bytes);
    hasher
}

/// The counts and the packed offset that the header of a packed sequence holds.
#[derive(Clone, Copy)]
struct Header {
    contigs: u64,
    bases: u64,
    names: u64,
    packed_offset: u64,
}

impl Header {
    /// Checks the magic, version and CRC-32 of a header and reads its fields.
    fn decode(header: &[u8]) -> Result<Self, IndexError> {
        if header.len() != HEADER_BYTES as usize
            || header[..16] != MAGIC
            || u32_at(header, 16) != Some(VERSION)
            || header[20..24] != [0; 4]
        {
            return Err(IndexError::format("not a packed sequence of this version"));
        }
        if u32_at(header, 68) != Some(crc32fast::hash(&header[..68]))
            || header[72..].iter().any(|b| *b != 0)
        {
            return Err(IndexError::Checksum("packed sequence header".into()));
        }
        let [contigs, bases, names, packed_offset] =
            [24, 32, 40, 48].map(|offset| u64_at(header, offset).unwrap_or_default());
        Ok(Self {
            contigs,
            bases,
            names,
            packed_offset,
        })
    }
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
