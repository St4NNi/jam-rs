//! Deterministic nonredundant support selection for overlapping alignments.
//!
//! Alignments are ranked by evidence before selection.  An alignment that
//! contributes at least one previously unsupported plasmid base is retained
//! as primary; an alignment with no novel supported base remains visible as a
//! secondary alignment and is classified by its overlap with primary support.
//! The final coverage summary unions all intervals, so repeated contigs never
//! inflate supported bases.

use crate::trace::coverage::{CoverageError, project_alignment};
use crate::trace::intervals::{IntervalError, circular_union, intersection, subtract, union};
use crate::trace::model::{BaseAlignment, BaseInterval, CoverageSummary};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use thiserror::Error;

/// Relationship between an alignment's supported query bases and support
/// selected before it.
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum OverlapClass {
    Disjoint,
    Partial,
    Contained,
    Identical,
}

/// Whether an alignment contributes nonredundant primary support.
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum SupportClass {
    Primary,
    SecondaryOverlap,
    SecondaryRedundant,
}

/// Deterministic classification for one input alignment.  `input_index`
/// points to the alignment's position in the input slice.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct AlignmentClassification {
    pub input_index: usize,
    pub support_class: SupportClass,
    pub overlap_class: OverlapClass,
    pub supported_bases: u64,
    pub overlapping_bases: u64,
    pub novel_intervals: Vec<BaseInterval>,
}

/// Result of primary/secondary mosaic selection.
#[derive(Clone, Debug, PartialEq)]
pub struct MosaicSelection {
    pub primary_indices: Vec<usize>,
    pub secondary_indices: Vec<usize>,
    pub classifications: Vec<AlignmentClassification>,
    pub primary_intervals: Vec<BaseInterval>,
    pub secondary_intervals: Vec<BaseInterval>,
    pub coverage: CoverageSummary,
}

#[derive(Debug, Error)]
pub enum MosaicError {
    #[error("plasmid length must be greater than zero")]
    ZeroPlasmidLength,
    #[error("alignment {index} could not be projected: {source}")]
    AlignmentProjection {
        index: usize,
        #[source]
        source: CoverageError,
    },
    #[error("coverage summary failed: {0}")]
    Coverage(#[from] CoverageError),
    #[error("interval operation failed: {0}")]
    Interval(#[from] IntervalError),
}

/// Select deterministic nonredundant primary support and classify every
/// alignment.  The input order is used only as the final tie-breaker after
/// stable sequence identifiers and coordinates, so parallel callers can sort
/// their input by the same fields and obtain identical results.
pub fn select_primary(
    plasmid_length: u64,
    alignments: &[BaseAlignment],
) -> Result<MosaicSelection, MosaicError> {
    if plasmid_length == 0 {
        return Err(MosaicError::ZeroPlasmidLength);
    }

    let mut projections = Vec::with_capacity(alignments.len());
    for (index, alignment) in alignments.iter().enumerate() {
        let projection = project_alignment(alignment)
            .map_err(|source| MosaicError::AlignmentProjection { index, source })?;
        let supported = circular_union(&projection.supported_intervals, plasmid_length)
            .map_err(CoverageError::from)
            .map_err(|source| MosaicError::AlignmentProjection { index, source })?;
        // Rank by nonredundant query support.  The edit-script counter is
        // deliberately kept for per-alignment diagnostics, rather than being
        // treated as a coordinate-union measure.
        let supported_bases = covered_bases(&supported);
        projections.push((supported, supported_bases));
    }

    let mut order: Vec<usize> = (0..alignments.len()).collect();
    order.sort_unstable_by(|left, right| {
        compare_alignment(
            &alignments[*left],
            projections[*left].1,
            *left,
            &alignments[*right],
            projections[*right].1,
            *right,
        )
    });

    let mut selected_support = Vec::new();
    let mut primary_support = Vec::new();
    let mut secondary_support = Vec::new();
    let mut primary_indices = Vec::new();
    let mut secondary_indices = Vec::new();
    let mut classifications: Vec<Option<AlignmentClassification>> = vec![None; alignments.len()];

    for index in order {
        let supported = &projections[index].0;
        let supported_bases = covered_bases(supported);
        let overlapping_bases = overlap_with_union(supported, &selected_support);
        let overlap_class = overlap_class(supported, &selected_support);
        let novel_intervals = supported
            .iter()
            .flat_map(|interval| subtract(*interval, &selected_support))
            .collect::<Vec<_>>();
        let is_primary = !novel_intervals.is_empty();
        if is_primary {
            primary_indices.push(index);
            primary_support.extend(supported.iter().copied());
            selected_support = union(
                selected_support
                    .into_iter()
                    .chain(supported.iter().copied()),
            );
        } else {
            secondary_indices.push(index);
            secondary_support.extend(supported.iter().copied());
        }
        classifications[index] = Some(AlignmentClassification {
            input_index: index,
            support_class: if is_primary {
                SupportClass::Primary
            } else if overlapping_bases > 0 {
                SupportClass::SecondaryOverlap
            } else {
                SupportClass::SecondaryRedundant
            },
            overlap_class,
            supported_bases,
            overlapping_bases,
            novel_intervals,
        });
    }

    primary_indices.sort_unstable();
    secondary_indices.sort_unstable();
    let classifications = classifications
        .into_iter()
        .enumerate()
        .map(|(index, classification)| {
            classification.unwrap_or_else(|| AlignmentClassification {
                input_index: index,
                support_class: SupportClass::SecondaryRedundant,
                overlap_class: OverlapClass::Disjoint,
                supported_bases: 0,
                overlapping_bases: 0,
                novel_intervals: Vec::new(),
            })
        })
        .collect();

    let primary_alignments: Vec<_> = primary_indices
        .iter()
        .map(|index| alignments[*index].clone())
        .collect();
    let secondary_alignments: Vec<_> = secondary_indices
        .iter()
        .map(|index| alignments[*index].clone())
        .collect();
    let coverage = crate::trace::coverage::summarize_coverage(
        plasmid_length,
        &primary_alignments,
        &secondary_alignments,
    )?;
    let primary_intervals = circular_union(&primary_support, plasmid_length)?;
    let secondary_intervals = circular_union(&secondary_support, plasmid_length)?;
    Ok(MosaicSelection {
        primary_indices,
        secondary_indices,
        classifications,
        primary_intervals,
        secondary_intervals,
        coverage,
    })
}

/// Alias that makes the nonredundant behavior explicit at call sites.
pub fn nonredundant_primary_support(
    plasmid_length: u64,
    alignments: &[BaseAlignment],
) -> Result<MosaicSelection, MosaicError> {
    select_primary(plasmid_length, alignments)
}

/// Alias for callers that describe the operation as mosaic summarisation.
pub fn summarize_mosaic(
    plasmid_length: u64,
    alignments: &[BaseAlignment],
) -> Result<MosaicSelection, MosaicError> {
    select_primary(plasmid_length, alignments)
}

/// Classify one candidate alignment against already selected primary support.
pub fn classify_secondary_overlap(
    plasmid_length: u64,
    alignment: &BaseAlignment,
    selected_support: &[BaseInterval],
) -> Result<OverlapClass, MosaicError> {
    let projection = project_alignment(alignment)
        .map_err(|source| MosaicError::AlignmentProjection { index: 0, source })?;
    let supported = circular_union(&projection.supported_intervals, plasmid_length)?;
    let selected = circular_union(selected_support, plasmid_length)?;
    Ok(overlap_class(&supported, &selected))
}

/// Compute the overlap class for two alignments' supported query intervals.
pub fn classify_overlap(
    plasmid_length: u64,
    left: &BaseAlignment,
    right: &BaseAlignment,
) -> Result<OverlapClass, MosaicError> {
    let left_projection = project_alignment(left)
        .map_err(|source| MosaicError::AlignmentProjection { index: 0, source })?;
    let right_projection = project_alignment(right)
        .map_err(|source| MosaicError::AlignmentProjection { index: 1, source })?;
    let left_support = circular_union(&left_projection.supported_intervals, plasmid_length)?;
    let right_support = circular_union(&right_projection.supported_intervals, plasmid_length)?;
    Ok(overlap_class(&left_support, &right_support))
}

fn compare_alignment(
    left: &BaseAlignment,
    left_supported_bases: u64,
    left_index: usize,
    right: &BaseAlignment,
    right_supported_bases: u64,
    right_index: usize,
) -> Ordering {
    right_supported_bases
        .cmp(&left_supported_bases)
        .then_with(|| right.score.cmp(&left.score))
        .then_with(|| right.chain_score.cmp(&left.chain_score))
        .then_with(|| right.matches.cmp(&left.matches))
        .then_with(|| left.substitutions.cmp(&right.substitutions))
        .then_with(|| left.insertions.cmp(&right.insertions))
        .then_with(|| left.deletions.cmp(&right.deletions))
        .then_with(|| left.contig_id.cmp(&right.contig_id))
        .then_with(|| left.target_interval.start.cmp(&right.target_interval.start))
        .then_with(|| left.target_interval.end.cmp(&right.target_interval.end))
        .then_with(|| strand_key(left).cmp(&strand_key(right)))
        .then_with(|| compare_query_segments(&left.query_segments, &right.query_segments))
        .then_with(|| left_index.cmp(&right_index))
}

fn strand_key(alignment: &BaseAlignment) -> u8 {
    match alignment.strand {
        crate::trace::model::Strand::Forward => 0,
        crate::trace::model::Strand::Reverse => 1,
    }
}

fn compare_query_segments(left: &[BaseInterval], right: &[BaseInterval]) -> Ordering {
    for (left_segment, right_segment) in left.iter().zip(right) {
        match left_segment.start.cmp(&right_segment.start) {
            Ordering::Equal => {}
            ordering => return ordering,
        }
        match left_segment.end.cmp(&right_segment.end) {
            Ordering::Equal => {}
            ordering => return ordering,
        }
    }
    left.len().cmp(&right.len())
}

fn overlap_class(left: &[BaseInterval], right: &[BaseInterval]) -> OverlapClass {
    if left.is_empty() || right.is_empty() {
        return OverlapClass::Disjoint;
    }
    let overlap = overlap_with_union(left, right);
    if overlap == 0 {
        OverlapClass::Disjoint
    } else if intervals_equal(left, right) {
        OverlapClass::Identical
    } else if overlap == covered_bases(left) {
        OverlapClass::Contained
    } else {
        OverlapClass::Partial
    }
}

fn intervals_equal(left: &[BaseInterval], right: &[BaseInterval]) -> bool {
    left == right
}

fn overlap_with_union(left: &[BaseInterval], right: &[BaseInterval]) -> u64 {
    left.iter()
        .map(|left_interval| {
            right
                .iter()
                .filter_map(|right_interval| intersection(*left_interval, *right_interval))
                .map(|intersection| intersection.len())
                .sum::<u64>()
        })
        .sum()
}

fn covered_bases(intervals: &[BaseInterval]) -> u64 {
    intervals.iter().map(|interval| interval.len()).sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::trace::model::{EditOperation, EditRun, Strand};

    fn alignment(contig: &str, start: u64, end: u64) -> BaseAlignment {
        let length = end - start;
        BaseAlignment {
            plasmid_id: "p".to_string(),
            metagenome_id: "m".to_string(),
            contig_id: contig.to_string(),
            strand: Strand::Forward,
            query_segments: vec![BaseInterval { start, end }],
            target_interval: BaseInterval {
                start: 0,
                end: length,
            },
            query_length: length,
            target_length: length,
            origin_crossing: false,
            score: length as i64,
            matches: length,
            substitutions: 0,
            insertions: 0,
            deletions: 0,
            cigar: format!("{length}="),
            edit_script: vec![EditRun {
                operation: EditOperation::Equal,
                length: u32::try_from(length).unwrap(),
            }],
            chain_score: length as i64,
            primary: false,
        }
    }

    #[test]
    fn overlapping_support_is_union_deduplicated() {
        let alignments = vec![alignment("a", 0, 100), alignment("b", 80, 160)];
        let mosaic = select_primary(176, &alignments).unwrap();
        assert_eq!(mosaic.coverage.supported_bases, 160);
        assert_eq!(
            mosaic.primary_intervals,
            vec![BaseInterval { start: 0, end: 160 }]
        );
        assert!(
            mosaic
                .classifications
                .iter()
                .all(|classification| classification.support_class == SupportClass::Primary)
        );
    }

    #[test]
    fn fully_repeated_alignment_is_secondary_overlap() {
        let alignments = vec![alignment("a", 0, 100), alignment("b", 0, 100)];
        let mosaic = select_primary(176, &alignments).unwrap();
        assert_eq!(mosaic.primary_indices.len(), 1);
        assert_eq!(mosaic.secondary_indices.len(), 1);
        assert_eq!(
            mosaic.classifications[1].support_class,
            SupportClass::SecondaryOverlap
        );
        assert_eq!(
            mosaic.classifications[1].overlap_class,
            OverlapClass::Identical
        );
    }

    #[test]
    fn primary_order_does_not_change_union() {
        let first = vec![alignment("a", 0, 40), alignment("b", 40, 80)];
        let second = vec![first[1].clone(), first[0].clone()];
        assert_eq!(
            select_primary(100, &first).unwrap().coverage,
            select_primary(100, &second).unwrap().coverage
        );
    }
}
