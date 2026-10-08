//! Base-level alignment: the alignment type, CIGAR, the banded kernel, anchors, chains and walks.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

// CIGAR strings follow SAM with the entry contig as the reference: `=` match, `X` substitution,
// `I` bases only in the query, `D` bases only in the contig.

mod anchors;
mod chain_alignment;
mod chains;
mod cigar;
mod kernel;
mod walk;

pub use anchors::{Anchor, AnchorConfig, QueryTable};
pub use chain_alignment::{BandConfig, align_chain};
pub use chains::{
    CHAINING_BYTES, Chain, ChainConfig, WindowChains, forward_extent, split_by_drift, window_chains,
};
#[cfg(test)]
pub use cigar::parse_cigar;
pub use cigar::{EditOperation, EditRun, format_cigar};
pub use kernel::{AlignmentConfig, AlignmentWorkspace};
pub use walk::{
    DenseResults, DenseWalk, DenseWalks, JoinContext, JoinLimits, JoinScratch, Walk,
    align_dense_walk, dense_walks, gap_continuation, holds_walk, long_enough, walk_alignment,
};

use serde::ser::SerializeStruct;
use serde::{Deserialize, Serialize, Serializer};
use thiserror::Error;

/// Whether a query base and a contig base match: the same A, C, G or T in either case. Unknown
/// bases (N and IUPAC codes) never match, not even themselves.
pub(crate) fn bases_match(query: u8, contig: u8) -> bool {
    query.eq_ignore_ascii_case(&contig)
        && matches!(query.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

/// Whether every byte is an upper-case A, C, G or T, so equal blocks of such bytes match base
/// by base. Other blocks are compared base by base with `bases_match`.
pub(crate) fn known_bases(bases: &[u8]) -> bool {
    bases
        .iter()
        .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
}

/// A half-open range of positions, `start..end`.
#[derive(Clone, Copy, Debug, Deserialize, Eq, Ord, PartialEq, PartialOrd, Serialize)]
pub struct Interval {
    /// The first position.
    pub start: u64,
    /// The position after the last.
    pub end: u64,
}

impl Interval {
    /// The interval `start..end`; fails when `start` lies after `end`.
    pub fn new(start: u64, end: u64) -> Result<Self, AlignmentError> {
        if start > end {
            return Err(AlignmentError::ReversedInterval { start, end });
        }
        Ok(Self { start, end })
    }

    /// The number of positions.
    pub fn len(self) -> u64 {
        self.end - self.start
    }

    /// Whether the interval holds no position.
    pub fn is_empty(self) -> bool {
        self.start == self.end
    }

    /// The union of `intervals`: sorted, with overlapping and touching intervals joined.
    pub fn union(mut intervals: Vec<Interval>) -> Vec<Interval> {
        intervals.sort_unstable();
        let mut union: Vec<Interval> = Vec::new();
        for interval in intervals {
            match union.last_mut() {
                Some(last) if interval.start <= last.end => last.end = last.end.max(interval.end),
                _ => union.push(interval),
            }
        }
        union
    }

    /// The summed length of `intervals`.
    pub fn total_length(intervals: &[Interval]) -> u64 {
        intervals.iter().map(|interval| interval.len()).sum()
    }

    /// This interval of positions on a circular sequence of `length` bases, which may run once
    /// past its end, as one or two intervals inside `0..length`.
    pub fn wrap(self, length: u64) -> impl Iterator<Item = Interval> {
        let (first, second) = if self.end <= length {
            (self, None)
        } else if self.start >= length {
            let (start, end) = (self.start - length, self.end - length);
            (Interval { start, end }, None)
        } else {
            let end = self.end - length;
            let second = Interval { start: 0, end };
            (
                Interval {
                    end: length,
                    ..self
                },
                Some(second),
            )
        };
        std::iter::once(first).chain(second)
    }
}

/// The contig strand that the query aligns to.
#[derive(Clone, Copy, Debug, Deserialize, Eq, Ord, PartialEq, PartialOrd, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    /// Whether this is the reverse strand.
    pub fn is_reverse(self) -> bool {
        self == Self::Reverse
    }
}

/// A gapped alignment of a query interval to a contig interval.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Alignment {
    /// The affine score of the runs.
    pub score: i32,
    /// The contig strand the query aligns to.
    pub strand: Strand,
    /// The aligned query bases; over the origin of a circular query they run past its length.
    pub query_interval: Interval,
    /// The aligned contig bases, in forward contig positions.
    pub contig_interval: Interval,
    /// The edit runs from the start of both intervals on; on the reverse strand the contig runs
    /// backwards from `contig_interval.end`.
    pub runs: Vec<EditRun>,
}

impl Alignment {
    /// The base counts of the runs.
    pub fn counts(&self) -> EditCounts {
        cigar::counts(&self.runs)
    }

    /// Matches divided by alignment columns.
    #[cfg(test)]
    pub fn identity(&self) -> f64 {
        self.counts().identity()
    }

    /// Checks that the runs are valid and span exactly both intervals.
    pub fn validate(&self) -> Result<(), AlignmentError> {
        let counts = cigar::checked_counts(&self.runs)?;
        if counts.query_bases() != self.query_interval.len()
            || counts.contig_bases() != self.contig_interval.len()
        {
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

/// Base counts of an alignment, one per CIGAR operation.
#[derive(Clone, Copy, Debug, Default, Deserialize, Eq, PartialEq, Serialize)]
pub struct EditCounts {
    /// Columns of equal bases (`=`).
    pub matches: u64,
    /// Columns of different bases (`X`).
    pub substitutions: u64,
    /// Bases only in the query (`I`).
    pub insertions: u64,
    /// Bases only in the contig (`D`).
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

    /// Whether the identity reaches `limit`, tested as `matches >= limit * columns` so that
    /// every caller agrees at the boundary.
    pub fn reaches_identity(self, limit: f64) -> bool {
        self.matches as f64 >= limit * self.columns() as f64
    }

    /// The query bases the alignment spans: matches, substitutions and insertions.
    pub fn query_bases(self) -> u64 {
        self.matches + self.substitutions + self.insertions
    }

    /// The contig bases the alignment spans: matches, substitutions and deletions.
    pub fn contig_bases(self) -> u64 {
        self.matches + self.substitutions + self.deletions
    }
}

impl std::ops::AddAssign for EditCounts {
    fn add_assign(&mut self, other: Self) {
        self.matches += other.matches;
        self.substitutions += other.substitutions;
        self.insertions += other.insertions;
        self.deletions += other.deletions;
    }
}

/// What can go wrong while finding anchors, chains and walks and while aligning.
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
    #[error("alignment query is empty")]
    EmptyQuery,
    #[error("alignment contig window is empty")]
    EmptyContig,
    #[error("alignment scores or limits are invalid")]
    InvalidConfig,
    #[error("alignment length overflows")]
    LengthOverflow,
    #[error("alignment band holds no cell of the input")]
    BandExcludesInput,
    #[error("alignment needs {cells} resident cells, more than {max_cells}")]
    MatrixTooLarge { cells: usize, max_cells: usize },
    #[error("no positive-scoring local alignment was found")]
    NoAlignment,
    #[error("alignment traceback left the band")]
    TracebackOutsideBand,
    #[error("alignment traceback is inconsistent")]
    InvalidTraceback,
    #[error("alignment to complete lies outside the given windows")]
    EndpointOutsideWindow,
    #[error("anchor configuration is invalid")]
    InvalidAnchorConfig,
    #[error("chain configuration is invalid")]
    InvalidChainConfig,
    #[error("query is too long for the query table")]
    QueryTooLong,
    #[error("{0} exceed the chaining budget")]
    ChainingBudget(Budgeted),
    #[error("alignment work was asked to stop")]
    Stopped,
}

/// What exceeded the chaining budget of an alignment window.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Budgeted {
    PeriodicChains,
    MatchStarts,
    PeriodicWalks,
    ChainingWork,
}

impl std::fmt::Display for Budgeted {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(match self {
            Self::PeriodicChains => "periodic chains",
            Self::MatchStarts => "dense match starts",
            Self::PeriodicWalks => "circular periodic walks",
            Self::ChainingWork => "predecessors examined in chaining",
        })
    }
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
    fn indels_take_bases() {
        alignment("5=2I3=", (0, 10), (0, 8)).validate().unwrap();
        alignment("5=2D3=", (0, 8), (0, 10)).validate().unwrap();
        assert_eq!(
            alignment("5=2I3=", (0, 10), (0, 10)).validate(),
            Err(AlignmentError::CigarSpanMismatch)
        );
    }

    #[test]
    fn identity_counts_columns() {
        let counts = alignment("6=1X1I2D", (0, 8), (0, 9)).counts();
        assert_eq!(counts.columns(), 10);
        assert_eq!(counts.identity(), 0.6);
        assert_eq!(EditCounts::default().identity(), 0.0);
    }

    #[test]
    fn alignment_json_fields() {
        let json = serde_json::to_value(alignment("4=1X", (2, 7), (10, 15))).unwrap();
        assert_eq!(json["cigar"], "4=1X");
        assert_eq!(json["identity"], 0.8);
        assert_eq!(json["contig_interval"]["start"], 10);
    }

    #[test]
    fn records_round_trip() {
        let counts = alignment("4=1X2I", (2, 9), (10, 15)).counts();
        let json = serde_json::to_string(&(Strand::Reverse, Interval::new(2, 9).unwrap(), counts));
        let json = json.unwrap();
        assert!(
            json.starts_with(r#"["reverse",{"start":2,"end":9},{"matches":4,"#),
            "{json}"
        );
        let back: (Strand, Interval, EditCounts) = serde_json::from_str(&json).unwrap();
        assert_eq!(
            back,
            (Strand::Reverse, Interval::new(2, 9).unwrap(), counts)
        );
        let runs = parse_cigar("3=1D").unwrap();
        let json = serde_json::to_string(&runs).unwrap();
        assert_eq!(
            json,
            r#"[{"operation":"=","length":3},{"operation":"D","length":1}]"#
        );
        assert_eq!(serde_json::from_str::<Vec<EditRun>>(&json).unwrap(), runs);
    }
}
