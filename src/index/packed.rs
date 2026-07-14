//! The packed sequence of one entry (stage 3): its contigs at 2 bits per base, with contig names
//! and the runs of bytes other than A, C, G and T. The format is the one of jam 0.10, so packed
//! files built by it are read unchanged. Lower case a, c, g and t are stored as upper case; other
//! bytes keep their value, in upper case.
//!
//! Layout, little endian ("entry coordinates" count bases over all contigs in order):
//!
//! | Bytes | Field |
//! | --- | --- |
//! | 0..16 | magic `jam-rs 2bit seq` and a zero byte |
//! | 16..20 | version (u32, 1) |
//! | 20..24 | zero |
//! | 24..32 | contigs c (u64) |
//! | 32..40 | bases n (u64) |
//! | 40..48 | runs r of other bytes (u64) |
//! | 48..56 | name bytes m (u64) |
//! | 56..64 | offset of the packed bases (u64) |
//! | 64..68 | CRC-32 of the tables, from byte 128 to the packed bases (u32) |
//! | 68..72 | CRC-32 of bytes 0..68 (u32) |
//! | 72..128 | zero |
//! | 8 (c + 1) | contig starts in entry coordinates, then n (u64 each) |
//! | 8 c | end of each contig name in the name bytes (u64 each) |
//! | m, padded to 8 | contig names, the first word of each FASTA header |
//! | 16 r | runs: start (u64), length (u32), byte (u8), three zero bytes |
//! | 4 ceil(n / 16384) | CRC-32 of every block of 4,096 packed bytes |
//! | ceil(n / 4) | the bases, four per byte, the first base in the low bits |

use super::IndexError;
use super::bits::{u32_at, u64_at};
use crate::seeds::base_code;
use crate::storage::ReadAt;
use std::ops::Range;
use std::sync::Arc;

const MAGIC: [u8; 16] = *b"jam-rs 2bit seq\0";
const VERSION: u32 = 1;
const HEADER_BYTES: u64 = 128;
/// Packed bytes per checksum block, so 16,384 bases.
const BLOCK_BYTES: u64 = 4096;
const RUN_BYTES: u64 = 16;

/// A run of one byte other than A, C, G or T, at entry coordinates `start..start + length`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Run {
    start: u64,
    length: u32,
    byte: u8,
}

fn valid_name(name: &[u8]) -> bool {
    !name.is_empty() && !name.iter().any(u8::is_ascii_whitespace)
}

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

    /// Bases added so far.
    pub fn bases(&self) -> u64 {
        self.bases
    }

    /// Contigs added so far.
    pub fn contig_count(&self) -> u64 {
        self.name_ends.len() as u64
    }

    /// The file in two parts: header and tables, then the packed bases, so a caller can write
    /// both without copying the bases.
    pub fn finish(mut self) -> (Vec<u8>, Vec<u8>) {
        let contigs = self.name_ends.len() as u64;
        self.starts.push(self.bases);
        let mut file = vec![0u8; HEADER_BYTES as usize];
        for value in self.starts.iter().chain(&self.name_ends) {
            file.extend_from_slice(&value.to_le_bytes());
        }
        file.extend_from_slice(&self.names);
        file.resize(file.len().next_multiple_of(8), 0);
        for run in &self.runs {
            file.extend_from_slice(&run.start.to_le_bytes());
            file.extend_from_slice(&run.length.to_le_bytes());
            file.extend_from_slice(&[run.byte, 0, 0, 0]);
        }
        for block in self.packed.chunks(BLOCK_BYTES as usize) {
            file.extend_from_slice(&crc32fast::hash(block).to_le_bytes());
        }
        let packed_offset = file.len() as u64;
        let tables_crc = crc32fast::hash(&file[HEADER_BYTES as usize..]);
        file[..16].copy_from_slice(&MAGIC);
        file[16..20].copy_from_slice(&VERSION.to_le_bytes());
        let fields = [
            contigs,
            self.bases,
            self.runs.len() as u64,
            self.names.len() as u64,
            packed_offset,
        ];
        for (index, value) in fields.into_iter().enumerate() {
            file[24 + 8 * index..32 + 8 * index].copy_from_slice(&value.to_le_bytes());
        }
        file[64..68].copy_from_slice(&tables_crc.to_le_bytes());
        let header_crc = crc32fast::hash(&file[..68]);
        file[68..72].copy_from_slice(&header_crc.to_le_bytes());
        (file, self.packed)
    }
}

/// An open packed sequence. Contig starts, runs and block checksums stay in memory; names and
/// bases are read on request.
pub struct PackedSequence {
    source: Arc<dyn ReadAt>,
    starts: Vec<u64>,
    names_offset: u64,
    name_ends: Vec<u64>,
    runs: Vec<Run>,
    block_crcs: Vec<u32>,
    packed_offset: u64,
}

impl PackedSequence {
    /// Reads and checks the header and tables.
    pub fn open(source: Arc<dyn ReadAt>) -> Result<Self, IndexError> {
        let mut header = [0u8; HEADER_BYTES as usize];
        source.read_exact_at(0, &mut header)?;
        if header[..16] != MAGIC || u32_at(&header, 16) != Some(VERSION) || header[20..24] != [0; 4]
        {
            return Err(IndexError::format("not a packed sequence of this version"));
        }
        if u32_at(&header, 68) != Some(crc32fast::hash(&header[..68]))
            || header[72..].iter().any(|b| *b != 0)
        {
            return Err(IndexError::Checksum("packed sequence header".into()));
        }
        let [contigs, bases, runs, names, packed_offset] =
            [24, 32, 40, 48, 56].map(|offset| u64_at(&header, offset).unwrap_or_default());
        let blocks = bases.div_ceil(4).div_ceil(BLOCK_BYTES);
        let tables_end = [
            contigs.checked_add(1).and_then(|n| n.checked_mul(8)),
            contigs.checked_mul(8),
            names.checked_next_multiple_of(8),
            runs.checked_mul(RUN_BYTES),
            blocks.checked_mul(4),
        ]
        .into_iter()
        .try_fold(HEADER_BYTES, |total, part| total.checked_add(part?));
        if tables_end != Some(packed_offset)
            || packed_offset.checked_add(bases.div_ceil(4)) != Some(source.size())
        {
            return Err(IndexError::format("packed sequence layout"));
        }
        let mut tables = vec![0u8; (packed_offset - HEADER_BYTES) as usize];
        source.read_exact_at(HEADER_BYTES, &mut tables)?;
        if u32_at(&header, 64) != Some(crc32fast::hash(&tables)) {
            return Err(IndexError::Checksum("packed sequence tables".into()));
        }
        Self::from_tables(
            source,
            &tables,
            [contigs, bases, runs, names, packed_offset],
        )
    }

    fn from_tables(
        source: Arc<dyn ReadAt>,
        tables: &[u8],
        [contigs, bases, runs, names, packed_offset]: [u64; 5],
    ) -> Result<Self, IndexError> {
        let words = |at: usize, count: u64| -> Vec<u64> {
            (0..count as usize)
                .map(|i| u64_at(tables, at + 8 * i).unwrap_or_default())
                .collect()
        };
        let starts = words(0, contigs + 1);
        let name_ends = words(8 * (contigs as usize + 1), contigs);
        let names_at = 8 * (2 * contigs as usize + 1);
        let runs_at = names_at + names.next_multiple_of(8) as usize;
        let runs: Vec<Run> = (0..runs as usize)
            .map(|i| {
                let at = runs_at + i * RUN_BYTES as usize;
                Run {
                    start: u64_at(tables, at).unwrap_or_default(),
                    length: u32_at(tables, at + 8).unwrap_or_default(),
                    byte: tables[at + 12],
                }
            })
            .collect();
        let crcs_at = runs_at + runs.len() * RUN_BYTES as usize;
        let block_crcs = (0..(tables.len() - crcs_at) / 4)
            .map(|i| u32_at(tables, crcs_at + 4 * i).unwrap_or_default())
            .collect();
        let ordered = starts.first() == Some(&0)
            && starts.last() == Some(&bases)
            && starts.is_sorted()
            && name_ends.is_sorted()
            && name_ends.last().copied().unwrap_or(0) == names;
        let mut previous_end = 0;
        let runs_valid = runs.iter().all(|run| {
            let end = run.start + u64::from(run.length);
            let valid = run.length > 0
                && run.start >= previous_end
                && end <= bases
                && base_code(run.byte).is_none()
                && !run.byte.is_ascii_lowercase();
            previous_end = end;
            valid
        });
        if !ordered || !runs_valid {
            return Err(IndexError::format("packed sequence tables"));
        }
        Ok(Self {
            source,
            starts,
            names_offset: HEADER_BYTES + names_at as u64,
            name_ends,
            runs,
            block_crcs,
            packed_offset,
        })
    }

    /// Bases of all contigs.
    pub fn bases(&self) -> u64 {
        self.starts.last().copied().unwrap_or(0)
    }

    /// Contigs of the entry.
    pub fn contig_count(&self) -> usize {
        self.name_ends.len()
    }

    /// The name of a contig, read from the source.
    pub fn contig_name(&self, contig: usize) -> Result<Vec<u8>, IndexError> {
        let end = *self
            .name_ends
            .get(contig)
            .ok_or_else(|| no_contig(contig))?;
        let start = contig
            .checked_sub(1)
            .map_or(0, |previous| self.name_ends[previous]);
        let mut name = vec![0; (end - start) as usize];
        self.source
            .read_exact_at(self.names_offset + start, &mut name)?;
        if !valid_name(&name) {
            return Err(IndexError::format("contig name"));
        }
        Ok(name)
    }

    /// Entry coordinates of a contig.
    pub fn contig_range(&self, contig: usize) -> Result<Range<u64>, IndexError> {
        match (self.starts.get(contig), self.starts.get(contig + 1)) {
            (Some(start), Some(end)) => Ok(*start..*end),
            _ => Err(no_contig(contig)),
        }
    }

    /// The contig that holds the base at an entry coordinate.
    pub fn contig_at(&self, position: u64) -> Option<usize> {
        if position >= self.bases() {
            return None;
        }
        Some(self.starts.partition_point(|start| *start <= position) - 1)
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
        let first = self
            .runs
            .partition_point(|run| run.start + u64::from(run.length) <= range.start);
        for run in self.runs[first..]
            .iter()
            .take_while(|run| run.start < range.end)
        {
            let from = run.start.max(range.start) - range.start;
            let to = (run.start + u64::from(run.length)).min(range.end) - range.start;
            out[from as usize..to as usize].fill(run.byte);
        }
        Ok(())
    }
}

fn no_contig(contig: usize) -> IndexError {
    IndexError::Input(format!("contig {contig} is not in this entry"))
}

#[cfg(test)]
mod tests;
