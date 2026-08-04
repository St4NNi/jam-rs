//! Edit runs, their CIGAR text and their base counts.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::{AlignmentError, EditCounts};
use serde::{Deserialize, Serialize};

/// One CIGAR operation; see the module documentation of `align` for its meaning.
#[derive(Clone, Copy, Debug, Deserialize, Eq, Ord, PartialEq, PartialOrd, Serialize)]
pub enum EditOperation {
    #[serde(rename = "=")]
    Match,
    #[serde(rename = "X")]
    Substitution,
    #[serde(rename = "I")]
    Insertion,
    #[serde(rename = "D")]
    Deletion,
}

impl EditOperation {
    fn symbol(self) -> char {
        match self {
            Self::Match => '=',
            Self::Substitution => 'X',
            Self::Insertion => 'I',
            Self::Deletion => 'D',
        }
    }
}

/// `length` columns of one operation.
#[derive(Clone, Copy, Debug, Deserialize, Eq, Ord, PartialEq, PartialOrd, Serialize)]
pub struct EditRun {
    /// The operation of every column of the run.
    pub operation: EditOperation,
    /// The number of columns.
    pub length: u32,
}

/// The CIGAR text of `runs`, such as `12=1X3I40=`.
pub fn format_cigar(runs: &[EditRun]) -> String {
    let mut cigar = String::with_capacity(runs.len() * 4);
    for run in runs {
        cigar.push_str(&run.length.to_string());
        cigar.push(run.operation.symbol());
    }
    cigar
}

/// Parses CIGAR text with the operations `=`, `X`, `I` and `D`. Empty runs and two neighbouring
/// runs of one operation are invalid.
pub fn parse_cigar(cigar: &str) -> Result<Vec<EditRun>, AlignmentError> {
    let mut runs: Vec<EditRun> = Vec::new();
    let mut length = 0u32;
    for byte in cigar.bytes() {
        if byte.is_ascii_digit() {
            length = length
                .checked_mul(10)
                .and_then(|value| value.checked_add(u32::from(byte - b'0')))
                .ok_or(AlignmentError::RunTooLong)?;
            continue;
        }
        let operation = match byte {
            b'=' => EditOperation::Match,
            b'X' => EditOperation::Substitution,
            b'I' => EditOperation::Insertion,
            b'D' => EditOperation::Deletion,
            _ => return Err(AlignmentError::InvalidCigar),
        };
        if length == 0 || runs.last().is_some_and(|run| run.operation == operation) {
            return Err(AlignmentError::InvalidCigar);
        }
        runs.push(EditRun { operation, length });
        length = 0;
    }
    if length != 0 || runs.is_empty() {
        return Err(AlignmentError::InvalidCigar);
    }
    Ok(runs)
}

/// Appends `more` to `runs`, joining neighbouring runs of one operation.
pub(super) fn append_runs(
    runs: &mut Vec<EditRun>,
    more: impl IntoIterator<Item = EditRun>,
) -> Result<(), AlignmentError> {
    for run in more {
        match runs.last_mut() {
            Some(last) if last.operation == run.operation => {
                last.length = last
                    .length
                    .checked_add(run.length)
                    .ok_or(AlignmentError::RunTooLong)?;
            }
            _ => runs.push(run),
        }
    }
    Ok(())
}

/// The runs of single operations in order.
pub(super) fn runs_of(operations: &[EditOperation]) -> Result<Vec<EditRun>, AlignmentError> {
    let mut runs = Vec::new();
    let single = |&operation| EditRun {
        operation,
        length: 1,
    };
    append_runs(&mut runs, operations.iter().map(single))?;
    Ok(runs)
}

/// The base counts of `runs`.
pub(super) fn counts(runs: &[EditRun]) -> EditCounts {
    let mut counts = EditCounts::default();
    for run in runs {
        let length = u64::from(run.length);
        match run.operation {
            EditOperation::Match => counts.matches += length,
            EditOperation::Substitution => counts.substitutions += length,
            EditOperation::Insertion => counts.insertions += length,
            EditOperation::Deletion => counts.deletions += length,
        }
    }
    counts
}

/// The base counts of `runs`, which must be non-empty and hold no empty run.
pub(super) fn checked_counts(runs: &[EditRun]) -> Result<EditCounts, AlignmentError> {
    if runs.is_empty() || runs.iter().any(|run| run.length == 0) {
        return Err(AlignmentError::InvalidCigar);
    }
    Ok(counts(runs))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cigar_text_round_trips() {
        let runs = parse_cigar("12=1X3I40=2D").unwrap();
        assert_eq!(runs.len(), 5);
        assert_eq!(format_cigar(&runs), "12=1X3I40=2D");
    }

    #[test]
    fn invalid_cigar_text_is_rejected() {
        for text in ["", "5", "=", "0=", "5M", "3=4=", "99999999999="] {
            assert!(parse_cigar(text).is_err(), "{text}");
        }
    }
}
