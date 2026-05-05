//! Base-level alignment: the alignment type and CIGAR, the banded alignment kernel, anchors and
//! chains.
//!
//! CIGAR strings follow SAM with the entry contig as the reference: `=` match, `X` substitution,
//! `I` bases only in the query, `D` bases only in the contig.

mod cigar;

pub use cigar::{EditOperation, EditRun, format_cigar, parse_cigar};

use serde::ser::SerializeStruct;
use serde::{Serialize, Serializer};
use thiserror::Error;

/// A half-open range of positions, `start..end`.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd, Serialize)]
pub struct Interval {
    pub start: u64,
    pub end: u64,
}

impl Interval {
    pub fn new(start: u64, end: u64) -> Result<Self, AlignmentError> {
        if start > end {
            return Err(AlignmentError::ReversedInterval { start, end });
        }
        Ok(Self { start, end })
    }

    pub fn len(self) -> u64 {
        self.end - self.start
    }

    pub fn is_empty(self) -> bool {
        self.start == self.end
    }
}

/// The contig strand that the query aligns to.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Strand {
    Forward,
    Reverse,
}

/// A gapped alignment of a query interval to a contig interval.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Alignment {
    pub score: i32,
    pub strand: Strand,
    pub query_interval: Interval,
    pub contig_interval: Interval,
    /// The edit runs from the start of both intervals on; on the reverse strand the contig runs
    /// backwards from `contig_interval.end`.
    pub runs: Vec<EditRun>,
}

/// Base counts of an alignment, one per CIGAR operation.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct EditCounts {
    pub matches: u64,
    pub substitutions: u64,
    pub insertions: u64,
    pub deletions: u64,
}

impl EditCounts {
    /// All alignment columns: every match, substitution, insertion and deletion.
    pub fn columns(self) -> u64 {
        self.matches + self.substitutions + self.insertions + self.deletions
    }

    /// Matches divided by columns; 0 for an empty alignment.
    pub fn identity(self) -> f64 {
        match self.columns() {
            0 => 0.0,
            columns => self.matches as f64 / columns as f64,
        }
    }
}

impl Alignment {
    pub fn counts(&self) -> EditCounts {
        cigar::counts(&self.runs)
    }

    pub fn identity(&self) -> f64 {
        self.counts().identity()
    }

    /// Checks that the runs are valid and span exactly both intervals.
    pub fn validate(&self) -> Result<(), AlignmentError> {
        let counts = cigar::checked_counts(&self.runs)?;
        let query_bases = counts.matches + counts.substitutions + counts.insertions;
        let contig_bases = counts.matches + counts.substitutions + counts.deletions;
        if query_bases != self.query_interval.len() || contig_bases != self.contig_interval.len() {
            return Err(AlignmentError::CigarSpanMismatch);
        }
        Ok(())
    }
}

impl Serialize for Alignment {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let counts = self.counts();
        let mut record = serializer.serialize_struct("Alignment", 10)?;
        record.serialize_field("score", &self.score)?;
        record.serialize_field("strand", &self.strand)?;
        record.serialize_field("query_interval", &self.query_interval)?;
        record.serialize_field("contig_interval", &self.contig_interval)?;
        record.serialize_field("matches", &counts.matches)?;
        record.serialize_field("substitutions", &counts.substitutions)?;
        record.serialize_field("insertions", &counts.insertions)?;
        record.serialize_field("deletions", &counts.deletions)?;
        record.serialize_field("identity", &counts.identity())?;
        record.serialize_field("cigar", &format_cigar(&self.runs))?;
        record.end()
    }
}

/// The complement of an IUPAC base; other bytes stay as they are.
pub fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        other => other,
    }
}

#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum AlignmentError {
    #[error("interval is reversed: {start}..{end}")]
    ReversedInterval { start: u64, end: u64 },
    #[error("edit run is longer than 2^32 - 1 bases")]
    RunTooLong,
    #[error("CIGAR text is invalid")]
    InvalidCigar,
    #[error("CIGAR does not span the alignment intervals")]
    CigarSpanMismatch,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn alignment(cigar: &str, query: (u64, u64), contig: (u64, u64)) -> Alignment {
        Alignment {
            score: 1,
            strand: Strand::Forward,
            query_interval: Interval::new(query.0, query.1).unwrap(),
            contig_interval: Interval::new(contig.0, contig.1).unwrap(),
            runs: parse_cigar(cigar).unwrap(),
        }
    }

    #[test]
    fn insertions_take_query_bases_and_deletions_take_contig_bases() {
        alignment("5=2I3=", (0, 10), (0, 8)).validate().unwrap();
        alignment("5=2D3=", (0, 8), (0, 10)).validate().unwrap();
        assert_eq!(
            alignment("5=2I3=", (0, 10), (0, 10)).validate(),
            Err(AlignmentError::CigarSpanMismatch)
        );
    }

    #[test]
    fn identity_counts_every_column() {
        let counts = alignment("6=1X1I2D", (0, 8), (0, 9)).counts();
        assert_eq!(counts.columns(), 10);
        assert_eq!(counts.identity(), 0.6);
        assert_eq!(EditCounts::default().identity(), 0.0);
    }

    #[test]
    fn serialized_alignments_carry_cigar_and_identity() {
        let json = serde_json::to_value(alignment("4=1X", (2, 7), (10, 15))).unwrap();
        assert_eq!(json["cigar"], "4=1X");
        assert_eq!(json["identity"], 0.8);
        assert_eq!(json["contig_interval"]["start"], 10);
    }

    #[test]
    fn complement_keeps_unknown_bytes() {
        assert_eq!(complement(b'a'), b'T');
        assert_eq!(complement(b'N'), b'N');
    }
}
