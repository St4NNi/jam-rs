//! The contig table of one entry: where each contig starts and the runs of bytes other than A, C,
//! G and T. It lies in the locator part of the seeds file, so stage 2 reads no packed sequence.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::IndexError;
use crate::seeds::base_code;
use integer_encoding::VarInt;
use std::ops::Range;

/// A run of one byte other than A, C, G or T, at entry coordinates `start..start + length`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct Run {
    pub(crate) start: u64,
    pub(crate) length: u32,
    pub(crate) byte: u8,
}

/// The contig starts and the runs of other bytes of one entry.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ContigTable {
    /// The start of each contig in entry coordinates, then the bases of the entry.
    starts: Vec<u64>,
    runs: Vec<Run>,
}

impl ContigTable {
    /// The table of contigs with these lengths and these runs (ordered, inside the entry).
    pub(crate) fn new(lengths: impl IntoIterator<Item = u64>, runs: Vec<Run>) -> Self {
        let mut starts = vec![0];
        for length in lengths {
            starts.push(starts[starts.len() - 1] + length);
        }
        Self { starts, runs }
    }

    /// Bases of all contigs.
    pub fn bases(&self) -> u64 {
        self.starts[self.starts.len() - 1]
    }

    /// Contigs of the entry.
    pub fn contig_count(&self) -> usize {
        self.starts.len() - 1
    }

    /// The start of each contig in entry coordinates, then the bases of the entry, so contig `i`
    /// spans `starts[i]..starts[i + 1]`.
    pub fn starts(&self) -> &[u64] {
        &self.starts
    }

    /// Entry coordinates of a contig.
    pub fn contig_range(&self, contig: usize) -> Result<Range<u64>, IndexError> {
        match (self.starts.get(contig), self.starts.get(contig + 1)) {
            (Some(start), Some(end)) => Ok(*start..*end),
            _ => Err(IndexError::Input(format!(
                "contig {contig} is not in this entry"
            ))),
        }
    }

    /// The runs that overlap bases `range` in entry coordinates.
    pub(crate) fn runs_in(&self, range: Range<u64>) -> impl Iterator<Item = &Run> {
        let first =
            (self.runs).partition_point(|run| run.start + u64::from(run.length) <= range.start);
        self.runs[first..]
            .iter()
            .take_while(move |run| run.start < range.end)
    }

    /// Appends the table in its file encoding, all varints: the contig count and each contig
    /// length, then the run count and per run its gap after the last run, length and byte.
    pub(crate) fn encode(&self, out: &mut Vec<u8>) {
        let mut buffer = [0u8; 10];
        let mut put = |value: u64, out: &mut Vec<u8>| {
            let used = value.encode_var(&mut buffer);
            out.extend_from_slice(&buffer[..used]);
        };
        put(self.contig_count() as u64, out);
        for pair in self.starts.windows(2) {
            put(pair[1] - pair[0], out);
        }
        put(self.runs.len() as u64, out);
        let mut end = 0;
        for run in &self.runs {
            put(run.start - end, out);
            put(u64::from(run.length), out);
            out.push(run.byte);
            end = run.start + u64::from(run.length);
        }
    }

    /// Decodes a table and checks that its runs are ordered, lie inside the entry and hold only
    /// upper case bytes other than A, C, G and T.
    pub(crate) fn decode(bytes: &[u8]) -> Result<Self, IndexError> {
        Self::try_decode(bytes).ok_or_else(|| IndexError::format("contig table"))
    }

    fn try_decode(bytes: &[u8]) -> Option<Self> {
        let mut at = 0;
        let take = |at: &mut usize| {
            let (value, used) = u64::decode_var(bytes.get(*at..)?)?;
            *at += used;
            Some(value)
        };
        let contigs = take(&mut at)?;
        // Every contig and run takes at least one byte, which bounds the counts.
        let mut starts = Vec::with_capacity(usize::try_from(contigs).ok()?.min(bytes.len()) + 1);
        starts.push(0u64);
        for _ in 0..contigs {
            let start = starts[starts.len() - 1].checked_add(take(&mut at)?)?;
            starts.push(start);
        }
        let bases = starts[starts.len() - 1];
        let count = take(&mut at)?;
        let mut runs = Vec::with_capacity(usize::try_from(count).ok()?.min(bytes.len()));
        let mut end = 0u64;
        for _ in 0..count {
            let start = end.checked_add(take(&mut at)?)?;
            let length = u32::try_from(take(&mut at)?).ok()?;
            let byte = *bytes.get(at)?;
            at += 1;
            end = start.checked_add(u64::from(length))?;
            let valid = length > 0 && base_code(byte).is_none() && !byte.is_ascii_lowercase();
            if !valid || end > bases {
                return None;
            }
            runs.push(Run {
                start,
                length,
                byte,
            });
        }
        (at == bytes.len()).then_some(Self { starts, runs })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table() -> ContigTable {
        let runs = vec![
            Run {
                start: 3,
                length: 2,
                byte: b'N',
            },
            Run {
                start: 300,
                length: 1,
                byte: b'R',
            },
        ];
        ContigTable::new([10, 0, 1_000_000], runs)
    }

    #[test]
    fn round_trip() {
        let table = table();
        let mut bytes = Vec::new();
        table.encode(&mut bytes);
        assert_eq!(ContigTable::decode(&bytes).unwrap(), table);
        assert_eq!(table.starts(), &[0, 10, 10, 1_000_010]);
        assert_eq!(table.contig_range(2).unwrap(), 10..1_000_010);
        assert!(table.contig_range(3).is_err());
        let runs: Vec<u64> = table.runs_in(4..301).map(|run| run.start).collect();
        assert_eq!(runs, [3, 300]);
        assert_eq!(table.runs_in(5..300).count(), 0);
    }

    #[test]
    fn damaged_tables_fail() {
        let mut bytes = Vec::new();
        table().encode(&mut bytes);
        assert!(ContigTable::decode(&bytes[..bytes.len() - 1]).is_err());
        let mut longer = bytes.clone();
        longer.push(0);
        assert!(ContigTable::decode(&longer).is_err());
        let mut lower = bytes.clone();
        *lower.last_mut().unwrap() = b'n';
        assert!(ContigTable::decode(&lower).is_err());
        // A run past the entry end.
        let mut past = Vec::new();
        ContigTable::new(
            [4],
            vec![Run {
                start: 3,
                length: 2,
                byte: b'N',
            }],
        )
        .encode(&mut past);
        assert!(ContigTable::decode(&past).is_err());
    }
}
