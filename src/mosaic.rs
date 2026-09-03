use crate::alignment::{Alignment, AlignmentError, EditOperation, Interval, Strand};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use thiserror::Error;

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct Fragment {
    pub contig_id: u32,
    pub query_segments: Vec<Interval>,
    pub alignment: Alignment,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct SelectedFragment {
    pub newly_supported_bases: u64,
    pub fragment: Fragment,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct Mosaic {
    pub query_length: u64,
    pub covered_bases: u64,
    pub covered_intervals: Vec<Interval>,
    pub gaps: Vec<Interval>,
    pub primary: Vec<SelectedFragment>,
    pub alternatives: Vec<Fragment>,
}

pub fn build_mosaic(query_length: u64, fragments: &[Fragment]) -> Result<Mosaic, MosaicError> {
    if query_length == 0 {
        return Err(MosaicError::EmptyQuery);
    }
    let mut prepared = fragments
        .iter()
        .cloned()
        .map(|fragment| {
            let support = supported_intervals(query_length, &fragment)?;
            let supported_bases = interval_bases(&support);
            Ok(PreparedFragment {
                fragment,
                support,
                supported_bases,
            })
        })
        .collect::<Result<Vec<_>, MosaicError>>()?;
    prepared.sort_by(compare_fragment);
    prepared.dedup_by(|left, right| left.fragment == right.fragment);

    let mut covered = Vec::new();
    let mut primary = Vec::new();
    let mut alternatives = Vec::new();
    for prepared in prepared {
        // ponytail: rebuild this small interval union; use an interval tree if mosaics grow large.
        let previous = interval_bases(&covered);
        let mut combined = covered.clone();
        combined.extend_from_slice(&prepared.support);
        let next = union_intervals(combined);
        let added = interval_bases(&next).saturating_sub(previous);
        if added == 0 {
            alternatives.push(prepared.fragment);
        } else {
            covered = next;
            primary.push(SelectedFragment {
                newly_supported_bases: added,
                fragment: prepared.fragment,
            });
        }
    }

    Ok(Mosaic {
        query_length,
        covered_bases: interval_bases(&covered),
        gaps: complement(&covered, query_length),
        covered_intervals: covered,
        primary,
        alternatives,
    })
}

struct PreparedFragment {
    fragment: Fragment,
    support: Vec<Interval>,
    supported_bases: u64,
}

fn compare_fragment(left: &PreparedFragment, right: &PreparedFragment) -> Ordering {
    right
        .fragment
        .alignment
        .score
        .cmp(&left.fragment.alignment.score)
        .then_with(|| right.supported_bases.cmp(&left.supported_bases))
        .then_with(|| {
            right
                .fragment
                .alignment
                .matches
                .cmp(&left.fragment.alignment.matches)
        })
        .then_with(|| {
            edit_bases(&left.fragment.alignment).cmp(&edit_bases(&right.fragment.alignment))
        })
        .then_with(|| left.fragment.contig_id.cmp(&right.fragment.contig_id))
        .then_with(|| {
            strand_key(left.fragment.alignment.strand)
                .cmp(&strand_key(right.fragment.alignment.strand))
        })
        .then_with(|| {
            left.fragment
                .alignment
                .target_interval
                .cmp(&right.fragment.alignment.target_interval)
        })
        .then_with(|| {
            left.fragment
                .alignment
                .query_interval
                .cmp(&right.fragment.alignment.query_interval)
        })
        .then_with(|| {
            left.fragment
                .query_segments
                .cmp(&right.fragment.query_segments)
        })
        .then_with(|| {
            left.fragment
                .alignment
                .cigar
                .cmp(&right.fragment.alignment.cigar)
        })
}

fn edit_bases(alignment: &Alignment) -> u64 {
    alignment
        .substitutions
        .saturating_add(alignment.insertions)
        .saturating_add(alignment.deletions)
}

fn strand_key(strand: Strand) -> u8 {
    match strand {
        Strand::Forward => 0,
        Strand::Reverse => 1,
    }
}

fn supported_intervals(
    query_length: u64,
    fragment: &Fragment,
) -> Result<Vec<Interval>, MosaicError> {
    fragment.alignment.validate_cigar()?;
    if fragment.query_segments.is_empty()
        || fragment
            .query_segments
            .iter()
            .any(|segment| segment.is_empty() || segment.end > query_length)
    {
        return Err(MosaicError::InvalidQuerySegments);
    }
    let mut ordered = fragment.query_segments.clone();
    ordered.sort_unstable();
    if ordered.windows(2).any(|pair| pair[0].end > pair[1].start) {
        return Err(MosaicError::InvalidQuerySegments);
    }
    let segment_bases = fragment
        .query_segments
        .iter()
        .try_fold(0u64, |total, segment| total.checked_add(segment.len()))
        .ok_or(MosaicError::LengthOverflow)?;
    if segment_bases != fragment.alignment.query_interval.len() {
        return Err(MosaicError::QuerySpanMismatch);
    }

    let mut query_offset = 0u64;
    let mut support = Vec::new();
    for run in &fragment.alignment.edit_script {
        let length = u64::from(run.length);
        match run.operation {
            EditOperation::Equal | EditOperation::Substitution => {
                project_segments(&fragment.query_segments, query_offset, length, &mut support)?;
                query_offset = query_offset
                    .checked_add(length)
                    .ok_or(MosaicError::LengthOverflow)?;
            }
            EditOperation::Deletion => {
                query_offset = query_offset
                    .checked_add(length)
                    .ok_or(MosaicError::LengthOverflow)?;
            }
            EditOperation::Insertion => {}
        }
    }
    if query_offset != segment_bases {
        return Err(MosaicError::QuerySpanMismatch);
    }
    Ok(union_intervals(support))
}

fn project_segments(
    segments: &[Interval],
    mut offset: u64,
    mut length: u64,
    output: &mut Vec<Interval>,
) -> Result<(), MosaicError> {
    for segment in segments {
        if offset >= segment.len() {
            offset -= segment.len();
            continue;
        }
        let count = length.min(segment.len() - offset);
        let start = segment
            .start
            .checked_add(offset)
            .ok_or(MosaicError::LengthOverflow)?;
        output.push(Interval::new(
            start,
            start
                .checked_add(count)
                .ok_or(MosaicError::LengthOverflow)?,
        )?);
        length -= count;
        offset = 0;
        if length == 0 {
            return Ok(());
        }
    }
    Err(MosaicError::QuerySpanMismatch)
}

fn union_intervals(mut intervals: Vec<Interval>) -> Vec<Interval> {
    intervals.sort_unstable();
    intervals
        .into_iter()
        .fold(Vec::new(), |mut union, interval| {
            if let Some(last) = union.last_mut()
                && interval.start <= last.end
            {
                last.end = last.end.max(interval.end);
            } else {
                union.push(interval);
            }
            union
        })
}

fn interval_bases(intervals: &[Interval]) -> u64 {
    intervals
        .iter()
        .fold(0, |total, interval| total.saturating_add(interval.len()))
}

fn complement(intervals: &[Interval], query_length: u64) -> Vec<Interval> {
    let mut gaps = Vec::new();
    let mut start = 0;
    for interval in intervals {
        if start < interval.start {
            gaps.push(Interval {
                start,
                end: interval.start,
            });
        }
        start = interval.end;
    }
    if start < query_length {
        gaps.push(Interval {
            start,
            end: query_length,
        });
    }
    gaps
}

#[derive(Debug, Error)]
pub enum MosaicError {
    #[error("mosaic query is empty")]
    EmptyQuery,
    #[error("mosaic query segments are invalid")]
    InvalidQuerySegments,
    #[error("mosaic query segments do not match the alignment CIGAR")]
    QuerySpanMismatch,
    #[error("mosaic coordinate length overflows")]
    LengthOverflow,
    #[error(transparent)]
    Alignment(#[from] AlignmentError),
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::{EditOperation, EditRun, cigar_from_runs};

    fn fragment(
        contig_id: u32,
        score: i32,
        segments: Vec<Interval>,
        runs: Vec<EditRun>,
    ) -> Fragment {
        let (matches, substitutions, insertions, deletions) = runs.iter().fold(
            (0, 0, 0, 0),
            |(matches, substitutions, insertions, deletions), run| {
                let length = u64::from(run.length);
                match run.operation {
                    EditOperation::Equal => {
                        (matches + length, substitutions, insertions, deletions)
                    }
                    EditOperation::Substitution => {
                        (matches, substitutions + length, insertions, deletions)
                    }
                    EditOperation::Insertion => {
                        (matches, substitutions, insertions + length, deletions)
                    }
                    EditOperation::Deletion => {
                        (matches, substitutions, insertions, deletions + length)
                    }
                }
            },
        );
        let query_bases = matches + substitutions + deletions;
        let target_bases = matches + substitutions + insertions;
        let cigar = cigar_from_runs(&runs).unwrap();
        Fragment {
            contig_id,
            query_segments: segments,
            alignment: Alignment {
                score,
                strand: Strand::Forward,
                query_interval: Interval::new(0, query_bases).unwrap(),
                target_interval: Interval::new(0, target_bases).unwrap(),
                matches,
                substitutions,
                insertions,
                deletions,
                cigar,
                edit_script: runs,
            },
        }
    }

    fn equal_fragment(contig_id: u32, score: i32, start: u64, end: u64) -> Fragment {
        fragment(
            contig_id,
            score,
            vec![Interval::new(start, end).unwrap()],
            vec![EditRun {
                operation: EditOperation::Equal,
                length: u32::try_from(end - start).unwrap(),
            }],
        )
    }

    #[test]
    fn deletions_are_not_reported_as_coverage() {
        let mosaic = build_mosaic(
            5,
            &[fragment(
                0,
                7,
                vec![Interval::new(0, 5).unwrap()],
                vec![
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
                        length: 2,
                    },
                ],
            )],
        )
        .unwrap();
        assert_eq!(mosaic.covered_bases, 4);
        assert_eq!(mosaic.gaps, vec![Interval::new(2, 3).unwrap()]);
    }

    #[test]
    fn overlapping_fragments_add_only_new_support() {
        let input = [
            equal_fragment(2, 7, 1, 5),
            equal_fragment(1, 8, 4, 10),
            equal_fragment(0, 10, 0, 6),
        ];
        let mosaic = build_mosaic(10, &input).unwrap();
        assert_eq!(mosaic.covered_bases, 10);
        assert_eq!(mosaic.primary.len(), 2);
        assert_eq!(mosaic.primary[0].fragment.contig_id, 0);
        assert_eq!(mosaic.primary[1].newly_supported_bases, 4);
        assert_eq!(mosaic.alternatives.len(), 1);
        let mut reversed = input.to_vec();
        reversed.reverse();
        assert_eq!(build_mosaic(10, &reversed).unwrap(), mosaic);
    }

    #[test]
    fn origin_split_segments_are_covered_without_wrapping_intervals() {
        let mosaic = build_mosaic(
            10,
            &[fragment(
                0,
                8,
                vec![Interval::new(8, 10).unwrap(), Interval::new(0, 2).unwrap()],
                vec![EditRun {
                    operation: EditOperation::Equal,
                    length: 4,
                }],
            )],
        )
        .unwrap();
        assert_eq!(
            mosaic.covered_intervals,
            vec![Interval::new(0, 2).unwrap(), Interval::new(8, 10).unwrap()]
        );
        assert_eq!(mosaic.gaps, vec![Interval::new(2, 8).unwrap()]);
    }
}
