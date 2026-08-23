//! Projection of alignment operations onto query-coordinate models.
//!
//! Equal and substitution operations support query bases.  Deletions consume
//! query bases but do not support them, while insertions consume only target
//! bases.  This distinction is deliberately kept in the projection so a
//! deletion cannot silently become coverage in a downstream summary.

use crate::trace::intervals::{
    IntervalError, circular_gap_complement, circular_union, covered_length, linear_gap_complement,
};
use crate::trace::model::{
    BaseAlignment, BaseInterval, CoordinateModel, CoverageSummary, EditOperation, EditRun,
    GapRecord,
};
use std::fmt::Write as _;
use thiserror::Error;

/// A single operation projected into query and target coordinates.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ProjectedSpan {
    pub operation: EditOperation,
    /// One or two ordered non-wrapping query pieces.  Insertions have no
    /// query pieces; deletions may have two when they cross the origin.
    pub query_segments: Vec<BaseInterval>,
    /// Target coordinates are linear and relative to the supplied alignment
    /// window.  Empty target intervals are used for query deletions.
    pub target_interval: BaseInterval,
}

/// CIGAR/edit-script evidence projected onto query coordinates.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CoverageProjection {
    /// Query intervals supported by equal or substitution operations only.
    pub supported_intervals: Vec<BaseInterval>,
    /// Query intervals consumed by aligned operations, including deletions.
    pub aligned_intervals: Vec<BaseInterval>,
    /// Query intervals consumed by deletion operations.  These are not
    /// included in `supported_intervals`.
    pub deletion_intervals: Vec<BaseInterval>,
    pub spans: Vec<ProjectedSpan>,
    pub query_consumed: u64,
    pub target_consumed: u64,
    pub supported_bases: u64,
    pub aligned_query_bases: u64,
    pub matches: u64,
    pub substitutions: u64,
    pub insertions: u64,
    pub deletions: u64,
    pub soft_clipped_bases: u64,
    pub cigar: String,
}

impl CoverageProjection {
    /// Supported query intervals before a caller-specific union operation.
    #[must_use]
    pub fn supported(&self) -> &[BaseInterval] {
        &self.supported_intervals
    }

    /// Query intervals consumed by equal, substitution, or deletion runs.
    #[must_use]
    pub fn aligned_spans(&self) -> &[BaseInterval] {
        &self.aligned_intervals
    }

    /// Query intervals missing from the target because of deletions.
    #[must_use]
    pub fn internal_deletions(&self) -> &[BaseInterval] {
        &self.deletion_intervals
    }
}

/// CIGAR parsing and projection errors.
#[derive(Clone, Debug, Eq, PartialEq, Error)]
pub enum CoverageError {
    #[error("CIGAR is empty")]
    EmptyCigar,
    #[error("CIGAR operation at byte {offset} has no length")]
    MissingLength { offset: usize },
    #[error("CIGAR contains a zero-length run at byte {offset}")]
    ZeroLengthRun { offset: usize },
    #[error("CIGAR has an invalid operation {operation:?} at byte {offset}")]
    InvalidOperation { operation: char, offset: usize },
    #[error("CIGAR length overflows a supported edit run at byte {offset}")]
    LengthOverflow { offset: usize },
    #[error("query segment [{start}, {end}) is reversed")]
    ReversedQuerySegment { start: u64, end: u64 },
    #[error("query segments [{left_start}, {left_end}) and [{right_start}, {right_end}) overlap")]
    OverlappingQuerySegments {
        left_start: u64,
        left_end: u64,
        right_start: u64,
        right_end: u64,
    },
    #[error("edit script consumes {consumed} query bases but segments contain {available}")]
    QueryLengthMismatch { consumed: u64, available: u64 },
    #[error(
        "edit script consumes {consumed} target bases but alignment target length is {expected}"
    )]
    TargetLengthMismatch { consumed: u64, expected: u64 },
    #[error("alignment target interval [{start}, {end}) is reversed")]
    ReversedTargetInterval { start: u64, end: u64 },
    #[error(
        "alignment target interval length {interval_length} differs from declared target length {declared_length}"
    )]
    TargetIntervalLengthMismatch {
        interval_length: u64,
        declared_length: u64,
    },
    #[error("edit run has zero length")]
    ZeroLengthEditRun,
    #[error("coordinate arithmetic overflow while projecting an edit script")]
    CoordinateOverflow,
    #[error("interval operation failed: {0}")]
    Interval(#[from] IntervalError),
}

/// Parse the common `=`, `X`, `M`, `I`, `D`, and `S` CIGAR operations into the
/// shared edit-run representation.  `M` is treated as an aligned operation;
/// when base identities are available, callers should provide `=`/`X` or an
/// authoritative edit script instead.
pub fn parse_cigar(cigar: &str) -> Result<Vec<EditRun>, CoverageError> {
    if cigar.is_empty() {
        return Err(CoverageError::EmptyCigar);
    }
    let bytes = cigar.as_bytes();
    let mut runs = Vec::new();
    let mut number: u64 = 0;
    let mut have_number = false;
    let mut saw_operation = false;
    for (offset, byte) in bytes.iter().copied().enumerate() {
        if byte.is_ascii_digit() {
            have_number = true;
            number = number
                .checked_mul(10)
                .and_then(|value| value.checked_add(u64::from(byte - b'0')))
                .ok_or(CoverageError::LengthOverflow { offset })?;
            continue;
        }
        if !have_number {
            return Err(CoverageError::MissingLength { offset });
        }
        if number == 0 {
            return Err(CoverageError::ZeroLengthRun { offset });
        }
        let operation = match byte {
            b'=' | b'M' => EditOperation::Equal,
            b'X' => EditOperation::Substitution,
            b'I' => EditOperation::Insertion,
            b'D' => EditOperation::Deletion,
            b'S' => EditOperation::SoftClip,
            other => {
                return Err(CoverageError::InvalidOperation {
                    operation: char::from(other),
                    offset,
                });
            }
        };
        let length = u32::try_from(number).map_err(|_| CoverageError::LengthOverflow { offset })?;
        push_run(&mut runs, operation, length);
        saw_operation = true;
        number = 0;
        have_number = false;
    }
    if !have_number {
        return if !saw_operation {
            Err(CoverageError::EmptyCigar)
        } else {
            Ok(runs)
        };
    }
    Err(CoverageError::MissingLength {
        offset: bytes.len().saturating_sub(1),
    })
}

/// Convert edit runs back to a canonical CIGAR string.
#[must_use]
pub fn cigar_from_edit_script(edit_script: &[EditRun]) -> String {
    let mut runs = Vec::new();
    for run in edit_script {
        if run.length == 0 {
            continue;
        }
        push_run(&mut runs, run.operation, run.length);
    }
    let mut cigar = String::new();
    for run in runs {
        let code = cigar_code(run.operation);
        let _ = write!(&mut cigar, "{}{}", run.length, code);
    }
    cigar
}

/// Project an edit script whose query sequence is represented by one or more
/// non-overlapping ordinary circular segments.  Touching segments are valid;
/// they are useful for explicit origin-crossing representations.  Target
/// coordinates begin at zero because the projection is about query coverage;
/// the resulting spans are relative to the supplied target alignment window.
pub fn project_edit_script(
    query_segments: &[BaseInterval],
    edit_script: &[EditRun],
) -> Result<CoverageProjection, CoverageError> {
    for segment in query_segments {
        if segment.start > segment.end {
            return Err(CoverageError::ReversedQuerySegment {
                start: segment.start,
                end: segment.end,
            });
        }
    }
    let mut sorted_segments = query_segments.to_vec();
    sorted_segments.sort_unstable_by_key(|segment| (segment.start, segment.end));
    for pair in sorted_segments.windows(2) {
        let [left, right] = pair else {
            continue;
        };
        if left.end > right.start {
            return Err(CoverageError::OverlappingQuerySegments {
                left_start: left.start,
                left_end: left.end,
                right_start: right.start,
                right_end: right.end,
            });
        }
    }
    let available = query_segments.iter().try_fold(0u64, |total, segment| {
        total
            .checked_add(segment.len())
            .ok_or(CoverageError::CoordinateOverflow)
    })?;
    let mut projection = CoverageProjection {
        supported_intervals: Vec::new(),
        aligned_intervals: Vec::new(),
        deletion_intervals: Vec::new(),
        spans: Vec::new(),
        query_consumed: 0,
        target_consumed: 0,
        supported_bases: 0,
        aligned_query_bases: 0,
        matches: 0,
        substitutions: 0,
        insertions: 0,
        deletions: 0,
        soft_clipped_bases: 0,
        cigar: cigar_from_edit_script(edit_script),
    };
    let mut query_cursor = 0u64;
    let mut target_cursor = 0u64;

    for run in edit_script {
        if run.length == 0 {
            return Err(CoverageError::ZeroLengthEditRun);
        }
        let length = u64::from(run.length);
        let consumes_query = matches!(
            run.operation,
            EditOperation::Equal
                | EditOperation::Substitution
                | EditOperation::Deletion
                | EditOperation::SoftClip
        );
        let consumes_target = matches!(
            run.operation,
            EditOperation::Equal | EditOperation::Substitution | EditOperation::Insertion
        );
        let query_piece = if consumes_query {
            let end = query_cursor
                .checked_add(length)
                .ok_or(CoverageError::CoordinateOverflow)?;
            let pieces = map_query_span(query_segments, query_cursor, length)?;
            query_cursor = end;
            pieces
        } else {
            Vec::new()
        };
        let target_interval = if consumes_target {
            let end = target_cursor
                .checked_add(length)
                .ok_or(CoverageError::CoordinateOverflow)?;
            let interval = BaseInterval {
                start: target_cursor,
                end,
            };
            target_cursor = end;
            interval
        } else {
            BaseInterval {
                start: target_cursor,
                end: target_cursor,
            }
        };

        match run.operation {
            EditOperation::Equal => {
                projection.matches = projection
                    .matches
                    .checked_add(length)
                    .ok_or(CoverageError::CoordinateOverflow)?;
                projection.supported_bases = projection
                    .supported_bases
                    .checked_add(length)
                    .ok_or(CoverageError::CoordinateOverflow)?;
                projection.aligned_query_bases = projection
                    .aligned_query_bases
                    .checked_add(length)
                    .ok_or(CoverageError::CoordinateOverflow)?;
                projection
                    .supported_intervals
                    .extend(query_piece.iter().copied());
                projection
                    .aligned_intervals
                    .extend(query_piece.iter().copied());
            }
            EditOperation::Substitution => {
                projection.substitutions = projection
                    .substitutions
                    .checked_add(length)
                    .ok_or(CoverageError::CoordinateOverflow)?;
                projection.supported_bases = projection
                    .supported_bases
                    .checked_add(length)
                    .ok_or(CoverageError::CoordinateOverflow)?;
                projection.aligned_query_bases = projection
                    .aligned_query_bases
                    .checked_add(length)
                    .ok_or(CoverageError::CoordinateOverflow)?;
                projection
                    .supported_intervals
                    .extend(query_piece.iter().copied());
                projection
                    .aligned_intervals
                    .extend(query_piece.iter().copied());
            }
            EditOperation::Insertion => {
                projection.insertions = projection
                    .insertions
                    .checked_add(length)
                    .ok_or(CoverageError::CoordinateOverflow)?;
            }
            EditOperation::Deletion => {
                projection.deletions = projection
                    .deletions
                    .checked_add(length)
                    .ok_or(CoverageError::CoordinateOverflow)?;
                projection.aligned_query_bases = projection
                    .aligned_query_bases
                    .checked_add(length)
                    .ok_or(CoverageError::CoordinateOverflow)?;
                projection
                    .aligned_intervals
                    .extend(query_piece.iter().copied());
                projection
                    .deletion_intervals
                    .extend(query_piece.iter().copied());
            }
            EditOperation::SoftClip => {
                projection.soft_clipped_bases = projection
                    .soft_clipped_bases
                    .checked_add(length)
                    .ok_or(CoverageError::CoordinateOverflow)?;
            }
        }
        projection.spans.push(ProjectedSpan {
            operation: run.operation,
            query_segments: query_piece,
            target_interval,
        });
    }

    if query_cursor != available {
        return Err(CoverageError::QueryLengthMismatch {
            consumed: query_cursor,
            available,
        });
    }
    projection.query_consumed = query_cursor;
    projection.target_consumed = target_cursor;
    Ok(projection)
}

/// Parse and project a CIGAR string onto query segments.
pub fn project_cigar(
    query_segments: &[BaseInterval],
    cigar: &str,
) -> Result<CoverageProjection, CoverageError> {
    let edit_script = parse_cigar(cigar)?;
    project_edit_script(query_segments, &edit_script)
}

/// Project a complete frozen alignment and validate its query/target lengths.
pub fn project_alignment(alignment: &BaseAlignment) -> Result<CoverageProjection, CoverageError> {
    if alignment.target_interval.start > alignment.target_interval.end {
        return Err(CoverageError::ReversedTargetInterval {
            start: alignment.target_interval.start,
            end: alignment.target_interval.end,
        });
    }
    let target_interval_length = alignment.target_interval.end - alignment.target_interval.start;
    if target_interval_length != alignment.target_length {
        return Err(CoverageError::TargetIntervalLengthMismatch {
            interval_length: target_interval_length,
            declared_length: alignment.target_length,
        });
    }
    let projection = project_edit_script(&alignment.query_segments, &alignment.edit_script)?;
    if projection.query_consumed != alignment.query_length {
        return Err(CoverageError::QueryLengthMismatch {
            consumed: projection.query_consumed,
            available: alignment.query_length,
        });
    }
    if projection.target_consumed != alignment.target_length {
        return Err(CoverageError::TargetLengthMismatch {
            consumed: projection.target_consumed,
            expected: alignment.target_length,
        });
    }
    Ok(projection)
}

/// Summarise already classified primary and secondary alignments.  The
/// primary set determines `supported_bases`; secondary intervals are retained
/// for review but never increase that count.
pub fn summarize_coverage(
    plasmid_length: u64,
    primary_alignments: &[BaseAlignment],
    secondary_alignments: &[BaseAlignment],
) -> Result<CoverageSummary, CoverageError> {
    if plasmid_length == 0 {
        return Err(CoverageError::Interval(IntervalError::ZeroLength));
    }
    let mut primary_support = Vec::new();
    for alignment in primary_alignments {
        primary_support.extend(project_alignment(alignment)?.supported_intervals);
    }
    let mut secondary_support = Vec::new();
    for alignment in secondary_alignments {
        secondary_support.extend(project_alignment(alignment)?.supported_intervals);
    }
    summary_from_intervals(plasmid_length, &primary_support, &secondary_support)
}

/// Summarise one list using each alignment's existing `primary` flag.
pub fn summarize_alignments(
    plasmid_length: u64,
    alignments: &[BaseAlignment],
) -> Result<CoverageSummary, CoverageError> {
    let (primary, secondary): (Vec<_>, Vec<_>) = alignments
        .iter()
        .cloned()
        .partition(|alignment| alignment.primary);
    summarize_coverage(plasmid_length, &primary, &secondary)
}

/// Build a summary from query support intervals.  This is useful for mosaic
/// selection, where projections have already been computed and no alignment
/// objects need to be cloned again.
pub fn summary_from_intervals(
    plasmid_length: u64,
    primary_intervals: &[BaseInterval],
    secondary_intervals: &[BaseInterval],
) -> Result<CoverageSummary, CoverageError> {
    if plasmid_length == 0 {
        return Err(CoverageError::Interval(IntervalError::ZeroLength));
    }
    let primary = circular_union(primary_intervals, plasmid_length)?;
    let secondary = circular_union(secondary_intervals, plasmid_length)?;
    let gaps = circular_gap_complement(&primary, plasmid_length)?
        .into_iter()
        .map(|interval| GapRecord {
            segments: vec![interval],
            wraps_origin: false,
            length: interval.len(),
            interval,
        })
        .collect::<Vec<_>>();
    let supported_bases = covered_length(&primary);
    let largest_gap = largest_circular_gap(&gaps, plasmid_length)?;
    Ok(CoverageSummary {
        plasmid_length,
        supported_bases,
        supported_fraction: supported_bases as f64 / plasmid_length as f64,
        primary_intervals: primary,
        secondary_intervals: secondary,
        gaps,
        largest_gap,
    })
}

/// Build a coverage summary using an explicit query coordinate model.  The
/// legacy [`summary_from_intervals`] API remains circular for compatibility;
/// new query-centered callers should use this function when linear and wrap
/// coordinates must be distinguished.
pub fn summary_from_intervals_model(
    query_length: u64,
    coordinate_model: CoordinateModel,
    primary_intervals: &[BaseInterval],
    secondary_intervals: &[BaseInterval],
) -> Result<CoverageSummary, CoverageError> {
    if query_length == 0 {
        return Err(CoverageError::Interval(IntervalError::ZeroLength));
    }
    // circular_union provides shared bounds/reversal validation for both
    // coordinate models.  Its result is also the correct ordinary union.
    let primary = circular_union(primary_intervals, query_length)?;
    let secondary = circular_union(secondary_intervals, query_length)?;
    let gaps = match coordinate_model {
        CoordinateModel::Wrap => circular_gap_complement(&primary, query_length)?,
        CoordinateModel::Linear | CoordinateModel::Undetermined => {
            linear_gap_complement(&primary, query_length)
        }
    }
    .into_iter()
    .map(|interval| GapRecord {
        segments: vec![interval],
        wraps_origin: false,
        length: interval.len(),
        interval,
    })
    .collect::<Vec<_>>();
    let supported_bases = covered_length(&primary);
    let largest_gap = if matches!(coordinate_model, CoordinateModel::Wrap) {
        largest_circular_gap(&gaps, query_length)?
    } else {
        gaps.iter().map(|gap| gap.length).max().unwrap_or(0)
    };
    Ok(CoverageSummary {
        plasmid_length: query_length,
        supported_bases,
        supported_fraction: supported_bases as f64 / query_length as f64,
        primary_intervals: primary,
        secondary_intervals: secondary,
        gaps,
        largest_gap,
    })
}

/// Alias for callers that describe the operation as model-aware coverage
/// summarisation.
pub fn summarize_coverage_model(
    query_length: u64,
    coordinate_model: CoordinateModel,
    primary_intervals: &[BaseInterval],
    secondary_intervals: &[BaseInterval],
) -> Result<CoverageSummary, CoverageError> {
    summary_from_intervals_model(
        query_length,
        coordinate_model,
        primary_intervals,
        secondary_intervals,
    )
}

/// Return the largest unsupported run on a circular plasmid while preserving
/// the serialized gap representation as ordinary, non-wrapping pieces.
fn largest_circular_gap(gaps: &[GapRecord], plasmid_length: u64) -> Result<u64, CoverageError> {
    let largest_linear = gaps.iter().map(|gap| gap.length).max().unwrap_or(0);
    if gaps.len() < 2 {
        return Ok(largest_linear);
    }
    let (Some(first), Some(last)) = (gaps.first(), gaps.last()) else {
        return Ok(largest_linear);
    };
    if first.interval.start == 0 && last.interval.end == plasmid_length {
        let wrapped = first
            .length
            .checked_add(last.length)
            .ok_or(CoverageError::CoordinateOverflow)?;
        Ok(largest_linear.max(wrapped))
    } else {
        Ok(largest_linear)
    }
}

/// Return the query pieces covered by a linear query offset and length.
fn map_query_span(
    query_segments: &[BaseInterval],
    offset: u64,
    length: u64,
) -> Result<Vec<BaseInterval>, CoverageError> {
    if length == 0 {
        return Ok(Vec::new());
    }
    let available = query_segments.iter().try_fold(0u64, |total, segment| {
        total
            .checked_add(segment.len())
            .ok_or(CoverageError::CoordinateOverflow)
    })?;
    let end = offset
        .checked_add(length)
        .ok_or(CoverageError::CoordinateOverflow)?;
    if end > available {
        return Err(CoverageError::QueryLengthMismatch {
            consumed: end,
            available,
        });
    }
    let mut result = Vec::new();
    let mut cursor = 0u64;
    let mut remaining_start = offset;
    let mut remaining = length;
    for segment in query_segments {
        let segment_length = segment.len();
        let segment_end = cursor
            .checked_add(segment_length)
            .ok_or(CoverageError::CoordinateOverflow)?;
        if remaining_start >= segment_end {
            cursor = segment_end;
            continue;
        }
        let local_start = remaining_start.saturating_sub(cursor);
        let piece_length = remaining.min(segment_length - local_start);
        let piece_start = segment
            .start
            .checked_add(local_start)
            .ok_or(CoverageError::CoordinateOverflow)?;
        result.push(BaseInterval {
            start: piece_start,
            end: piece_start
                .checked_add(piece_length)
                .ok_or(CoverageError::CoordinateOverflow)?,
        });
        remaining -= piece_length;
        remaining_start = segment_end;
        cursor = segment_end;
        if remaining == 0 {
            break;
        }
    }
    Ok(result)
}

fn push_run(runs: &mut Vec<EditRun>, operation: EditOperation, length: u32) {
    if let Some(last) = runs.last_mut()
        && last.operation == operation
    {
        if let Some(total) = last.length.checked_add(length) {
            last.length = total;
        } else {
            runs.push(EditRun { operation, length });
        }
    } else {
        runs.push(EditRun { operation, length });
    }
}

const fn cigar_code(operation: EditOperation) -> char {
    match operation {
        EditOperation::Equal => '=',
        EditOperation::Substitution => 'X',
        EditOperation::Insertion => 'I',
        EditOperation::Deletion => 'D',
        EditOperation::SoftClip => 'S',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cigar_parser_and_projection_keep_deletions_out_of_support() {
        let projection = project_cigar(&[BaseInterval { start: 10, end: 19 }], "3=2I2D4X").unwrap();
        assert_eq!(projection.supported_bases, 7);
        assert_eq!(projection.insertions, 2);
        assert_eq!(projection.deletions, 2);
        assert_eq!(
            projection.deletion_intervals,
            vec![BaseInterval { start: 13, end: 15 }]
        );
        assert_eq!(
            projection.supported_intervals,
            vec![
                BaseInterval { start: 10, end: 13 },
                BaseInterval { start: 15, end: 19 }
            ]
        );
    }

    #[test]
    fn projection_maps_runs_across_origin_segments() {
        let projection = project_edit_script(
            &[
                BaseInterval { start: 8, end: 10 },
                BaseInterval { start: 0, end: 4 },
            ],
            &[
                EditRun {
                    operation: EditOperation::Equal,
                    length: 2,
                },
                EditRun {
                    operation: EditOperation::Deletion,
                    length: 1,
                },
                EditRun {
                    operation: EditOperation::Equal,
                    length: 3,
                },
            ],
        )
        .unwrap();
        assert_eq!(projection.query_consumed, 6);
        assert_eq!(
            projection.deletion_intervals,
            vec![BaseInterval { start: 0, end: 1 }]
        );
        assert_eq!(
            projection.supported_intervals,
            vec![
                BaseInterval { start: 8, end: 10 },
                BaseInterval { start: 1, end: 4 }
            ]
        );
    }

    #[test]
    fn malformed_cigar_is_rejected() {
        assert!(parse_cigar("3=").is_ok());
        assert!(matches!(
            parse_cigar("=3"),
            Err(CoverageError::MissingLength { .. })
        ));
        assert!(matches!(
            parse_cigar("0="),
            Err(CoverageError::ZeroLengthRun { .. })
        ));
        assert!(matches!(
            parse_cigar("3"),
            Err(CoverageError::MissingLength { .. })
        ));
    }

    #[test]
    fn coverage_summary_unions_primary_support_and_reports_gaps() {
        let summary = summary_from_intervals(
            10,
            &[
                BaseInterval { start: 0, end: 4 },
                BaseInterval { start: 3, end: 7 },
            ],
            &[BaseInterval { start: 4, end: 9 }],
        )
        .unwrap();
        assert_eq!(summary.supported_bases, 7);
        assert_eq!(
            summary.primary_intervals,
            vec![BaseInterval { start: 0, end: 7 }]
        );
        assert_eq!(summary.gaps[0].interval, BaseInterval { start: 7, end: 10 });
        assert_eq!(
            summary.gaps[0].segments,
            vec![BaseInterval { start: 7, end: 10 }]
        );
        assert!(!summary.gaps[0].wraps_origin);
    }

    #[test]
    fn largest_circular_gap_combines_terminal_gap_pieces() {
        let gaps = vec![
            GapRecord {
                interval: BaseInterval { start: 0, end: 2 },
                segments: vec![BaseInterval { start: 0, end: 2 }],
                wraps_origin: false,
                length: 2,
            },
            GapRecord {
                interval: BaseInterval { start: 6, end: 10 },
                segments: vec![BaseInterval { start: 6, end: 10 }],
                wraps_origin: false,
                length: 4,
            },
        ];
        assert_eq!(largest_circular_gap(&gaps, 10).unwrap(), 6);
        assert_eq!(largest_circular_gap(&[], 10).unwrap(), 0);
        assert_eq!(
            largest_circular_gap(
                &[GapRecord {
                    interval: BaseInterval { start: 0, end: 10 },
                    segments: vec![BaseInterval { start: 0, end: 10 }],
                    wraps_origin: false,
                    length: 10,
                }],
                10,
            )
            .unwrap(),
            10
        );
    }
}
