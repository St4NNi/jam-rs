//! Checked interval operations used by trace coverage and mosaic reporting.
//!
//! Every interval in this module is zero-based and half-open.  Circular
//! intervals are represented as one or two ordinary [`BaseInterval`] values;
//! no operation accepts a reversed `BaseInterval` as an implicit wrap.

use crate::trace::model::BaseInterval;
use thiserror::Error;

/// Errors raised while normalising intervals against a circular sequence.
#[derive(Clone, Debug, Eq, PartialEq, Error)]
pub enum IntervalError {
    #[error("circular sequence length must be greater than zero")]
    ZeroLength,
    #[error("interval [{start}, {end}) is outside sequence length {length}")]
    OutOfBounds { start: u64, end: u64, length: u64 },
    #[error("interval [{start}, {end}) is reversed")]
    Reversed { start: u64, end: u64 },
}

/// Split a possibly origin-crossing circular interval into non-wrapping
/// intervals.  `start == end` denotes an empty interval.  Use
/// [`full_circle`] when a complete circle is intended.
pub fn split_circular(
    start: u64,
    end: u64,
    length: u64,
) -> Result<Vec<BaseInterval>, IntervalError> {
    if length == 0 {
        return Err(IntervalError::ZeroLength);
    }
    if start > length || end > length {
        return Err(IntervalError::OutOfBounds { start, end, length });
    }
    if start == end {
        return Ok(Vec::new());
    }
    if start < end {
        return Ok(vec![
            BaseInterval::new(start, end).map_err(|_| IntervalError::Reversed { start, end })?,
        ]);
    }

    let mut segments = Vec::with_capacity(2);
    if start < length {
        segments.push(
            BaseInterval::new(start, length)
                .map_err(|_| IntervalError::Reversed { start, end: length })?,
        );
    }
    if end > 0 {
        segments.push(
            BaseInterval::new(0, end).map_err(|_| IntervalError::Reversed { start: 0, end })?,
        );
    }
    Ok(segments)
}

/// Return the one non-wrapping interval covering an entire circular
/// sequence.
pub fn full_circle(length: u64) -> Result<BaseInterval, IntervalError> {
    if length == 0 {
        return Err(IntervalError::ZeroLength);
    }
    BaseInterval::new(0, length).map_err(|_| IntervalError::Reversed {
        start: 0,
        end: length,
    })
}

/// Validate ordinary intervals and merge overlaps and touching intervals.
pub fn union<I>(intervals: I) -> Vec<BaseInterval>
where
    I: IntoIterator<Item = BaseInterval>,
{
    let mut sorted: Vec<_> = intervals
        .into_iter()
        .filter(|interval| !interval.is_empty())
        .collect();
    sorted.sort_unstable_by_key(|interval| (interval.start, interval.end));

    let mut merged: Vec<BaseInterval> = Vec::with_capacity(sorted.len());
    for interval in sorted {
        if let Some(last) = merged.last_mut()
            && interval.start <= last.end
        {
            if interval.end > last.end {
                last.end = interval.end;
            }
            continue;
        }
        merged.push(interval);
    }
    merged
}

/// Alias with a descriptive name for callers handling multiple interval
/// sources.
#[must_use]
pub fn union_intervals(intervals: Vec<BaseInterval>) -> Vec<BaseInterval> {
    union(intervals)
}

/// Validate and union non-wrapping segments in a circular sequence.
pub fn circular_union(
    segments: &[BaseInterval],
    length: u64,
) -> Result<Vec<BaseInterval>, IntervalError> {
    if length == 0 {
        return if segments.iter().all(|interval| interval.is_empty()) {
            Ok(Vec::new())
        } else {
            Err(IntervalError::ZeroLength)
        };
    }
    for interval in segments {
        validate(*interval, length)?;
    }
    Ok(union(segments.iter().copied()))
}

/// Validate ordinary interval coordinates against a sequence length.
pub fn validate(interval: BaseInterval, length: u64) -> Result<(), IntervalError> {
    if interval.start > interval.end {
        return Err(IntervalError::Reversed {
            start: interval.start,
            end: interval.end,
        });
    }
    if interval.end > length {
        return Err(IntervalError::OutOfBounds {
            start: interval.start,
            end: interval.end,
            length,
        });
    }
    Ok(())
}

/// Validate and normalise a list of segments.  The returned list is sorted,
/// merged, and remains non-wrapping.
pub fn normalize_segments(
    segments: &[BaseInterval],
    length: u64,
) -> Result<Vec<BaseInterval>, IntervalError> {
    circular_union(segments, length)
}

/// Compute the length of a union of already-normalised intervals.
#[must_use]
pub fn covered_length(intervals: &[BaseInterval]) -> u64 {
    intervals.iter().map(|interval| interval.len()).sum()
}

/// Return the overlap length of two ordinary intervals.
#[must_use]
pub fn overlap_len(left: BaseInterval, right: BaseInterval) -> u64 {
    left.end
        .min(right.end)
        .saturating_sub(left.start.max(right.start))
}

/// Return the intersection of two ordinary intervals, if non-empty.
#[must_use]
pub fn intersection(left: BaseInterval, right: BaseInterval) -> Option<BaseInterval> {
    let start = left.start.max(right.start);
    let end = left.end.min(right.end);
    (start < end).then_some(BaseInterval { start, end })
}

/// Subtract a union of covered intervals from one ordinary interval.
#[must_use]
pub fn subtract(interval: BaseInterval, covered: &[BaseInterval]) -> Vec<BaseInterval> {
    if interval.is_empty() {
        return Vec::new();
    }
    let mut result = Vec::new();
    let mut cursor = interval.start;
    for covered_interval in covered {
        if covered_interval.end <= cursor {
            continue;
        }
        if covered_interval.start >= interval.end {
            break;
        }
        if covered_interval.start > cursor {
            result.push(BaseInterval {
                start: cursor,
                end: covered_interval.start.min(interval.end),
            });
        }
        cursor = cursor.max(covered_interval.end);
        if cursor >= interval.end {
            break;
        }
    }
    if cursor < interval.end {
        result.push(BaseInterval {
            start: cursor,
            end: interval.end,
        });
    }
    result
}

/// Subtract a normalised circular union from the complete sequence.  Gaps
/// are returned as ordinary intervals in coordinate order.  If an uncovered
/// gap crosses the origin it is represented by its two non-wrapping pieces,
/// matching the representation used for origin-crossing alignments.
pub fn circular_gap_complement(
    covered: &[BaseInterval],
    length: u64,
) -> Result<Vec<BaseInterval>, IntervalError> {
    if length == 0 {
        return if covered.iter().all(|interval| interval.is_empty()) {
            Ok(Vec::new())
        } else {
            Err(IntervalError::ZeroLength)
        };
    }
    let covered = circular_union(covered, length)?;
    Ok(complement(&covered, length))
}

/// Alias used by coverage callers that describe unsupported sequence as gaps.
pub fn circular_gaps(
    covered: &[BaseInterval],
    length: u64,
) -> Result<Vec<BaseInterval>, IntervalError> {
    circular_gap_complement(covered, length)
}

/// Compute the linear complement of a validated union in `[0, length)`.
#[must_use]
pub fn complement(covered: &[BaseInterval], length: u64) -> Vec<BaseInterval> {
    if length == 0 {
        return Vec::new();
    }
    let mut gaps = Vec::new();
    let mut cursor = 0;
    for interval in covered {
        if interval.start > cursor {
            gaps.push(BaseInterval {
                start: cursor,
                end: interval.start,
            });
        }
        cursor = cursor.max(interval.end);
        if cursor >= length {
            break;
        }
    }
    if cursor < length {
        gaps.push(BaseInterval {
            start: cursor,
            end: length,
        });
    }
    gaps
}

/// Compute the complement of a linear interval union in `[0, length)`.
/// Inputs are normalized internally so callers may provide overlapping or
/// out-of-order support pieces.  This is the linear counterpart to
/// [`circular_gap_complement`].
#[must_use]
pub fn linear_gap_complement(covered: &[BaseInterval], length: u64) -> Vec<BaseInterval> {
    complement(&union(covered.iter().copied()), length)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn circular_split_preserves_order_at_origin() {
        assert_eq!(
            split_circular(8, 2, 10).unwrap(),
            vec![
                BaseInterval { start: 8, end: 10 },
                BaseInterval { start: 0, end: 2 }
            ]
        );
    }

    #[test]
    fn union_deduplicates_overlapping_and_touching_intervals() {
        assert_eq!(
            union([
                BaseInterval { start: 8, end: 10 },
                BaseInterval { start: 0, end: 3 },
                BaseInterval { start: 2, end: 8 },
            ]),
            vec![BaseInterval { start: 0, end: 10 }]
        );
    }

    #[test]
    fn circular_complement_has_two_origin_pieces_when_needed() {
        assert_eq!(
            circular_gap_complement(&[BaseInterval { start: 3, end: 7 }], 10).unwrap(),
            vec![
                BaseInterval { start: 0, end: 3 },
                BaseInterval { start: 7, end: 10 }
            ]
        );
    }

    #[test]
    fn subtraction_does_not_reintroduce_covered_bases() {
        assert_eq!(
            subtract(
                BaseInterval { start: 0, end: 10 },
                &[
                    BaseInterval { start: 2, end: 5 },
                    BaseInterval { start: 7, end: 8 }
                ]
            ),
            vec![
                BaseInterval { start: 0, end: 2 },
                BaseInterval { start: 5, end: 7 },
                BaseInterval { start: 8, end: 10 }
            ]
        );
    }

    #[test]
    fn circular_operations_reject_invalid_coordinates() {
        assert!(matches!(
            circular_union(&[BaseInterval { start: 0, end: 11 }], 10),
            Err(IntervalError::OutOfBounds { .. })
        ));
        assert!(matches!(
            circular_union(&[BaseInterval { start: 8, end: 2 }], 10),
            Err(IntervalError::Reversed { .. })
        ));
        assert!(matches!(
            split_circular(11, 0, 10),
            Err(IntervalError::OutOfBounds { .. })
        ));
        assert!(matches!(full_circle(0), Err(IntervalError::ZeroLength)));
    }

    #[test]
    fn circular_union_bounds_coverage_and_handles_full_or_zero_cases() {
        let empty = circular_union(&[], 10).unwrap();
        assert_eq!(covered_length(&empty), 0);

        let full = circular_union(
            &[
                BaseInterval { start: 0, end: 6 },
                BaseInterval { start: 5, end: 10 },
            ],
            10,
        )
        .unwrap();
        assert_eq!(full, vec![BaseInterval { start: 0, end: 10 }]);
        assert_eq!(covered_length(&full), 10);
    }

    #[test]
    fn circular_gap_complement_is_origin_rotation_invariant_in_length() {
        let source = [
            BaseInterval { start: 10, end: 12 },
            BaseInterval { start: 0, end: 3 },
        ];
        let reference = circular_gap_complement(&source, 12).unwrap();
        assert_eq!(covered_length(&reference), 7);

        // Exercise every possible stored-origin rotation without a random
        // dependency.  The coordinates rotate, but covered and gap lengths
        // remain invariant.
        for shift in 0..12 {
            let mut rotated = Vec::new();
            for interval in source {
                let start = (interval.start + shift) % 12;
                let end = (start + interval.len()) % 12;
                rotated.extend(split_circular(start, end, 12).unwrap());
            }
            let gaps = circular_gap_complement(&rotated, 12).unwrap();
            assert_eq!(covered_length(&gaps), 7);
            assert_eq!(gaps.iter().map(|interval| interval.len()).sum::<u64>(), 7);
        }
    }

    #[test]
    fn linear_gap_complement_normalizes_input_order_and_overlap() {
        assert_eq!(
            linear_gap_complement(
                &[
                    BaseInterval { start: 8, end: 10 },
                    BaseInterval { start: 2, end: 6 },
                    BaseInterval { start: 4, end: 9 },
                ],
                12,
            ),
            vec![
                BaseInterval { start: 0, end: 2 },
                BaseInterval { start: 10, end: 12 },
            ]
        );
    }
}
