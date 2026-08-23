//! Deterministic query-coordinate fragment mosaics for overlapping alignments.
//!
//! Alignments are ranked by evidence before selection.  An alignment that
//! contributes at least one previously unsupported plasmid base is retained
//! as primary; an alignment with no novel supported base remains visible as a
//! secondary alignment and is classified by its overlap with primary support.
//! The final coverage summary unions all intervals, so repeated contigs never
//! inflate supported bases.

use crate::trace::config::MOSAIC_ALGORITHM_ID;
use crate::trace::coverage::{CoverageError, project_alignment};
use crate::trace::intervals::{
    IntervalError, circular_gap_complement, circular_union, complement, covered_length,
    intersection, subtract, union,
};
use crate::trace::model::{
    AlignmentMosaicEvidence, AlignmentRole, BaseAlignment, BaseInterval, CoordinateModel,
    CoverageSummary, FragmentMosaicSummary, GapRecord, MosaicAtomicInterval,
    MosaicSelectionComponents, TopologyAssessment, TopologyEvidence, TopologyModelEvidence,
    TopologyModelSummary, TopologyRequested,
};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::fmt::Write as _;
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
    #[error("query length must be greater than zero")]
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
    #[error("alignment {index} has query coordinates outside query length {query_length}")]
    QueryCoordinateOutOfBounds { index: usize, query_length: u64 },
    #[error("alignment {index} consumes more than one query length")]
    QuerySpanExceedsLength { index: usize },
}

/// The projected evidence used to construct a query-coordinate fragment
/// mosaic.  The projection is retained separately for support, aligned span,
/// and deletion accounting because a `D` operation is an aligned query span
/// but is not base-level support.
struct FragmentAlignment<'a> {
    alignment: &'a BaseAlignment,
    alignment_id: String,
    canonical_key: String,
    supported: Vec<BaseInterval>,
    aligned: Vec<BaseInterval>,
    deletions: Vec<BaseInterval>,
}

#[derive(Clone, Copy, Debug, Default)]
struct AlignmentEvidenceCounts {
    primary: u64,
    secondary: u64,
}

/// Build a deterministic query-coordinate fragment mosaic.
///
/// `CoordinateModel::Linear` and `CoordinateModel::Undetermined` use ordinary
/// linear complements.  Origin-crossing alignments are omitted from those
/// models so a linear display can never silently wrap.  `CoordinateModel::Wrap`
/// accepts split origin segments and uses circular unions and complements.
/// Alignment IDs are generated from the complete alignment content, so the
/// summary is stable when callers provide the accepted alignments in a
/// different order.
pub fn summarize_fragment_mosaic(
    query_length: u64,
    coordinate_model: CoordinateModel,
    alignments: &[BaseAlignment],
) -> Result<FragmentMosaicSummary, MosaicError> {
    if query_length == 0 {
        return Err(MosaicError::ZeroPlasmidLength);
    }

    let mut order: Vec<usize> = (0..alignments.len()).collect();
    let canonical_keys: Vec<String> = alignments.iter().map(canonical_alignment_key).collect();
    order.sort_unstable_by(|left, right| {
        canonical_keys[*left]
            .cmp(&canonical_keys[*right])
            .then_with(|| left.cmp(right))
    });

    let mut infos = Vec::with_capacity(alignments.len());
    let mut previous_key: Option<&str> = None;
    let mut duplicate_ordinal = 0_u32;
    for index in order {
        let alignment = &alignments[index];
        let key = &canonical_keys[index];
        if previous_key == Some(key.as_str()) {
            duplicate_ordinal = duplicate_ordinal.saturating_add(1);
        } else {
            previous_key = Some(key.as_str());
            duplicate_ordinal = 0;
        }
        validate_query_coordinates(index, alignment, query_length)?;

        let is_linear_compatible = matches!(
            coordinate_model,
            CoordinateModel::Linear | CoordinateModel::Undetermined
        ) && !alignment.origin_crossing
            && query_segments_are_linear(&alignment.query_segments);
        if matches!(
            coordinate_model,
            CoordinateModel::Linear | CoordinateModel::Undetermined
        ) && !is_linear_compatible
        {
            // Auto/unknown comparisons do not treat origin-crossing evidence
            // as linear support.  It remains available in the wrap model.
            continue;
        }

        let projection = project_alignment(alignment)
            .map_err(|source| MosaicError::AlignmentProjection { index, source })?;
        let supported = coordinate_union(
            &projection.supported_intervals,
            query_length,
            coordinate_model,
        )?;
        let aligned = coordinate_union(
            &projection.aligned_intervals,
            query_length,
            coordinate_model,
        )?;
        let deletions = coordinate_union(
            &projection.deletion_intervals,
            query_length,
            coordinate_model,
        )?;
        let alignment_id = if alignment.alignment_id.is_empty() {
            format!(
                "alignment-{}-{:04}",
                crate::provenance::sha256_bytes(key.as_bytes()),
                duplicate_ordinal
            )
        } else if duplicate_ordinal == 0 {
            alignment.alignment_id.clone()
        } else {
            format!("{}-{:04}", alignment.alignment_id, duplicate_ordinal)
        };
        infos.push(FragmentAlignment {
            alignment,
            alignment_id,
            canonical_key: key.clone(),
            supported,
            aligned,
            deletions,
        });
    }

    let boundaries = atomic_boundaries(query_length, &infos);
    let mut atomic_intervals = Vec::new();
    let mut counts = vec![AlignmentEvidenceCounts::default(); infos.len()];

    for pair in boundaries.windows(2) {
        let interval = BaseInterval {
            start: pair[0],
            end: pair[1],
        };
        if interval.is_empty() {
            continue;
        }
        let mut covering: Vec<usize> = infos
            .iter()
            .enumerate()
            .filter_map(|(index, info)| {
                interval_is_covered(interval, &info.supported).then_some(index)
            })
            .collect();
        covering.sort_unstable_by(|left, right| {
            compare_fragment_alignment(&infos[*left], &infos[*right])
        });

        let primary_index = covering.first().copied();
        let primary_alignment_id = primary_index.map(|index| {
            counts[index].primary = counts[index].primary.saturating_add(interval.len());
            infos[index].alignment_id.clone()
        });
        let mut alternative_alignment_ids = Vec::new();
        let primary_skip = if primary_index.is_some() { 1 } else { 0 };
        for index in covering.into_iter().skip(primary_skip) {
            counts[index].secondary = counts[index].secondary.saturating_add(interval.len());
            alternative_alignment_ids.push(infos[index].alignment_id.clone());
        }
        atomic_intervals.push(MosaicAtomicInterval {
            interval,
            primary_alignment_id,
            alternative_alignment_ids,
        });
    }

    let all_supported: Vec<BaseInterval> = infos
        .iter()
        .flat_map(|info| info.supported.iter().copied())
        .collect();
    let all_aligned: Vec<BaseInterval> = infos
        .iter()
        .flat_map(|info| info.aligned.iter().copied())
        .collect();
    let all_deletions: Vec<BaseInterval> = infos
        .iter()
        .flat_map(|info| info.deletions.iter().copied())
        .collect();
    let nonrepetitive_support: Vec<BaseInterval> = infos
        .iter()
        .filter(|info| info.alignment.seed_evidence.nonrepetitive_anchor_count > 0)
        .flat_map(|info| info.supported.iter().copied())
        .collect();
    let common_support: Vec<BaseInterval> = infos
        .iter()
        .filter(|info| info.alignment.seed_evidence.common_anchor_count > 0)
        .flat_map(|info| info.supported.iter().copied())
        .collect();
    let repeat_only_support: Vec<BaseInterval> = infos
        .iter()
        .filter(|info| {
            info.alignment.seed_evidence.repetitive_seed_count > 0
                && info.alignment.seed_evidence.nonrepetitive_anchor_count == 0
        })
        .flat_map(|info| info.supported.iter().copied())
        .collect();
    let covered_intervals = coordinate_union(&all_supported, query_length, coordinate_model)?;
    let aligned_intervals = coordinate_union(&all_aligned, query_length, coordinate_model)?;
    let deletion_intervals = coordinate_union(&all_deletions, query_length, coordinate_model)?;
    let nonrepetitive_intervals =
        coordinate_union(&nonrepetitive_support, query_length, coordinate_model)?;
    let common_sequence_intervals =
        coordinate_union(&common_support, query_length, coordinate_model)?;
    let repeat_only_intervals =
        coordinate_union(&repeat_only_support, query_length, coordinate_model)?;
    let unsupported_gaps = coordinate_gaps(&covered_intervals, query_length, coordinate_model)?;
    let base_covered_bases = covered_length(&covered_intervals);
    let aligned_span_bases = covered_length(&aligned_intervals);
    let nonrepetitive_supported_bases = covered_length(&nonrepetitive_intervals);
    let common_sequence_supported_bases = covered_length(&common_sequence_intervals);
    let repeat_only_supported_bases = covered_length(&repeat_only_intervals);
    let total_supported_bases = infos.iter().fold(0_u64, |total, info| {
        total.saturating_add(covered_length(&info.supported))
    });
    let redundant_overlap_bases = total_supported_bases.saturating_sub(base_covered_bases);

    let mut evidence = Vec::with_capacity(infos.len());
    let mut supporting_contigs = Vec::new();
    let mut local_alignment_score_sum = 0_i64;
    let mut nonrepetitive_anchor_evidence = 0_u64;
    let mut fragment_count = 0_u64;
    let mut alternative_alignment_count = 0_u64;
    for (index, info) in infos.iter().enumerate() {
        if !info.supported.is_empty() || !info.aligned.is_empty() {
            supporting_contigs.push(info.alignment.contig_id.clone());
        }
        if counts[index].primary > 0 {
            local_alignment_score_sum =
                local_alignment_score_sum.saturating_add(info.alignment.score);
            nonrepetitive_anchor_evidence = nonrepetitive_anchor_evidence
                .saturating_add(info.alignment.seed_evidence.nonrepetitive_anchor_count);
            fragment_count = fragment_count.saturating_add(1);
        }
        if counts[index].secondary > 0 {
            alternative_alignment_count = alternative_alignment_count.saturating_add(1);
        }
        let repeat_only = info.alignment.seed_evidence.repetitive_seed_count > 0
            && info.alignment.seed_evidence.nonrepetitive_anchor_count == 0;
        let common_sequence = info.alignment.seed_evidence.common_anchor_count > 0
            && info.alignment.seed_evidence.nonrepetitive_anchor_count == 0;
        let role = if info.alignment.origin_crossing {
            AlignmentRole::OriginCrossing
        } else if repeat_only {
            AlignmentRole::RepeatOnly
        } else if common_sequence {
            AlignmentRole::CommonSequence
        } else if counts[index].primary > 0 {
            AlignmentRole::PrimaryMosaic
        } else if counts[index].secondary > 0 {
            AlignmentRole::OverlappingSupport
        } else {
            AlignmentRole::AlternativeMapping
        };
        evidence.push(AlignmentMosaicEvidence {
            alignment_id: info.alignment_id.clone(),
            primary_supported_bases: counts[index].primary,
            secondary_supported_bases: counts[index].secondary,
            newly_supported_bases: counts[index].primary,
            role,
        });
    }
    supporting_contigs.sort_unstable();
    supporting_contigs.dedup();
    evidence.sort_by(|left, right| left.alignment_id.cmp(&right.alignment_id));

    let selection_components = MosaicSelectionComponents {
        local_alignment_score_sum,
        newly_supported_query_bases: base_covered_bases,
        nonrepetitive_anchor_evidence,
        redundant_overlap_bases,
        coordinate_contradictions: 0,
        sequence_contradictions: 0,
        fragment_count,
    };

    Ok(FragmentMosaicSummary {
        mosaic_algorithm: MOSAIC_ALGORITHM_ID.to_string(),
        coordinate_model,
        base_covered_bases,
        base_coverage_fraction: base_covered_bases as f64 / query_length as f64,
        aligned_span_bases,
        aligned_span_fraction: aligned_span_bases as f64 / query_length as f64,
        covered_intervals,
        unsupported_gaps,
        alignment_deletions: deletion_intervals,
        common_sequence_intervals,
        repeat_only_intervals,
        supporting_contigs,
        accepted_alignment_count: infos.len() as u64,
        alternative_alignment_count,
        nonrepetitive_supported_bases,
        common_sequence_supported_bases,
        repeat_only_supported_bases,
        atomic_intervals,
        alignment_evidence: evidence,
        selection_components,
    })
}

/// Alias retained for callers that use the shorter model-building name.
pub fn build_fragment_mosaic(
    query_length: u64,
    coordinate_model: CoordinateModel,
    alignments: &[BaseAlignment],
) -> Result<FragmentMosaicSummary, MosaicError> {
    summarize_fragment_mosaic(query_length, coordinate_model, alignments)
}

/// Short alias for callers that treat the operation as the fragment-mosaic
/// constructor itself.
pub fn fragment_mosaic(
    query_length: u64,
    coordinate_model: CoordinateModel,
    alignments: &[BaseAlignment],
) -> Result<FragmentMosaicSummary, MosaicError> {
    summarize_fragment_mosaic(query_length, coordinate_model, alignments)
}

/// Build linear and, when requested, wrap-aware models from one accepted
/// alignment slice.  Auto selection is based on the signed difference in
/// newly supported query bases; a difference must be strictly greater than
/// `margin_bases` before a coordinate model is selected.  A selected model is
/// a coordinate representation supported by the evidence, not a biological
/// topology call.
pub fn assess_topology(
    query_length: u64,
    topology_requested: TopologyRequested,
    margin_bases: u64,
    alignments: &[BaseAlignment],
) -> Result<TopologyAssessment, MosaicError> {
    let linear_mosaic =
        summarize_fragment_mosaic(query_length, CoordinateModel::Linear, alignments)?;
    let linear_model = model_summary(&linear_mosaic);
    let needs_wrap = matches!(
        topology_requested,
        TopologyRequested::Circular | TopologyRequested::Auto | TopologyRequested::Unknown
    );
    let wrap_model = if needs_wrap {
        Some(model_summary(&summarize_fragment_mosaic(
            query_length,
            CoordinateModel::Wrap,
            alignments,
        )?))
    } else {
        None
    };

    let (coordinate_model, topology_evidence) = match topology_requested {
        TopologyRequested::Linear => (
            CoordinateModel::Linear,
            if linear_mosaic.base_covered_bases > 0 {
                TopologyEvidence::LinearSupported
            } else {
                TopologyEvidence::Insufficient
            },
        ),
        TopologyRequested::Circular => (
            CoordinateModel::Wrap,
            if wrap_model
                .as_ref()
                .is_some_and(|model| model.mosaic.base_covered_bases > 0)
            {
                TopologyEvidence::WrapSupported
            } else {
                TopologyEvidence::Insufficient
            },
        ),
        TopologyRequested::Unknown => (
            CoordinateModel::Undetermined,
            TopologyEvidence::Undetermined,
        ),
        TopologyRequested::Auto => {
            let wrap = wrap_model
                .as_ref()
                .expect("auto always builds a wrap model");
            let linear_bases = linear_mosaic.base_covered_bases;
            let wrap_bases = wrap.mosaic.base_covered_bases;
            if wrap_bases > linear_bases.saturating_add(margin_bases) {
                (CoordinateModel::Wrap, TopologyEvidence::WrapSupported)
            } else if linear_bases > wrap_bases.saturating_add(margin_bases) {
                (CoordinateModel::Linear, TopologyEvidence::LinearSupported)
            } else if wrap_bases == linear_bases && wrap_bases > 0 {
                (
                    CoordinateModel::Undetermined,
                    TopologyEvidence::BothCompatible,
                )
            } else {
                (
                    CoordinateModel::Undetermined,
                    TopologyEvidence::Undetermined,
                )
            }
        }
    };

    Ok(TopologyAssessment {
        topology_requested,
        coordinate_model,
        topology_evidence,
        selection_margin_bases: margin_bases,
        linear_model,
        wrap_model,
    })
}

/// Alias for code that names the operation as topology-model assessment.
pub fn assess_topology_models(
    query_length: u64,
    topology_requested: TopologyRequested,
    margin_bases: u64,
    alignments: &[BaseAlignment],
) -> Result<TopologyAssessment, MosaicError> {
    assess_topology(query_length, topology_requested, margin_bases, alignments)
}

fn model_summary(mosaic: &FragmentMosaicSummary) -> TopologyModelSummary {
    let terminal_gap_bases = if matches!(
        mosaic.coordinate_model,
        CoordinateModel::Linear | CoordinateModel::Undetermined
    ) {
        let query_length = mosaic.base_covered_bases.saturating_add(
            mosaic
                .unsupported_gaps
                .iter()
                .fold(0_u64, |total, gap| total.saturating_add(gap.length)),
        );
        mosaic
            .unsupported_gaps
            .iter()
            .filter(|gap| gap.interval.start == 0 || gap.interval.end == query_length)
            .fold(0_u64, |total, gap| total.saturating_add(gap.length))
    } else {
        0
    };
    let origin_crossing_alignment_count = mosaic
        .alignment_evidence
        .iter()
        .filter(|evidence| evidence.role == AlignmentRole::OriginCrossing)
        .count() as u64;
    TopologyModelSummary {
        coordinate_model: mosaic.coordinate_model,
        mosaic: mosaic.clone(),
        evidence: TopologyModelEvidence {
            newly_supported_query_bases: mosaic.base_covered_bases,
            alignment_quality_sum: mosaic.selection_components.local_alignment_score_sum,
            nonrepetitive_anchor_support: mosaic.selection_components.nonrepetitive_anchor_evidence,
            origin_crossing_alignment_count,
            unexplained_terminal_gap_bases: terminal_gap_bases,
            contradictory_alignment_count: mosaic
                .selection_components
                .coordinate_contradictions
                .saturating_add(mosaic.selection_components.sequence_contradictions),
            fragment_count: mosaic.selection_components.fragment_count,
        },
    }
}

fn validate_query_coordinates(
    index: usize,
    alignment: &BaseAlignment,
    query_length: u64,
) -> Result<(), MosaicError> {
    let query_span = alignment
        .query_segments
        .iter()
        .try_fold(0_u64, |total, segment| total.checked_add(segment.len()));
    if alignment.query_length > query_length || query_span.is_none_or(|span| span > query_length) {
        return Err(MosaicError::QuerySpanExceedsLength { index });
    }
    if alignment
        .query_segments
        .iter()
        .any(|segment| segment.end > query_length)
    {
        return Err(MosaicError::QueryCoordinateOutOfBounds {
            index,
            query_length,
        });
    }
    Ok(())
}

fn query_segments_are_linear(segments: &[BaseInterval]) -> bool {
    segments.windows(2).all(|pair| pair[0].end <= pair[1].start)
}

fn coordinate_union(
    intervals: &[BaseInterval],
    query_length: u64,
    coordinate_model: CoordinateModel,
) -> Result<Vec<BaseInterval>, MosaicError> {
    match coordinate_model {
        CoordinateModel::Wrap => Ok(circular_union(intervals, query_length)?),
        CoordinateModel::Linear | CoordinateModel::Undetermined => {
            for interval in intervals {
                if interval.end > query_length {
                    return Err(MosaicError::QueryCoordinateOutOfBounds {
                        index: 0,
                        query_length,
                    });
                }
            }
            Ok(union(intervals.iter().copied()))
        }
    }
}

fn coordinate_gaps(
    covered: &[BaseInterval],
    query_length: u64,
    coordinate_model: CoordinateModel,
) -> Result<Vec<GapRecord>, MosaicError> {
    let intervals = match coordinate_model {
        CoordinateModel::Wrap => circular_gap_complement(covered, query_length)?,
        CoordinateModel::Linear | CoordinateModel::Undetermined => {
            complement(covered, query_length)
        }
    };
    let mut gaps: Vec<GapRecord> = intervals
        .into_iter()
        .map(|interval| GapRecord {
            segments: vec![interval],
            wraps_origin: false,
            length: interval.len(),
            interval,
        })
        .collect();
    if matches!(coordinate_model, CoordinateModel::Wrap) && gaps.len() >= 2 {
        let first_is_terminal = gaps.first().is_some_and(|gap| gap.interval.start == 0);
        let last_is_terminal = gaps
            .last()
            .is_some_and(|gap| gap.interval.end == query_length);
        if first_is_terminal && last_is_terminal {
            let first = gaps.remove(0);
            let last = gaps.pop().expect("at least two terminal gaps");
            let mut segments = first.segments;
            segments.extend(last.segments);
            let length = first.length.saturating_add(last.length);
            gaps.insert(
                0,
                GapRecord {
                    interval: first.interval,
                    segments,
                    wraps_origin: true,
                    length,
                },
            );
        }
    }
    Ok(gaps)
}

fn atomic_boundaries(query_length: u64, infos: &[FragmentAlignment<'_>]) -> Vec<u64> {
    let mut boundaries = BTreeSet::from([0, query_length]);
    for info in infos {
        for interval in info
            .supported
            .iter()
            .chain(info.aligned.iter())
            .chain(info.deletions.iter())
        {
            boundaries.insert(interval.start);
            boundaries.insert(interval.end);
        }
    }
    boundaries.into_iter().collect()
}

fn interval_is_covered(interval: BaseInterval, coverage: &[BaseInterval]) -> bool {
    coverage
        .iter()
        .any(|covered| covered.start <= interval.start && covered.end >= interval.end)
}

fn compare_fragment_alignment(
    left: &FragmentAlignment<'_>,
    right: &FragmentAlignment<'_>,
) -> Ordering {
    right
        .alignment
        .score
        .cmp(&left.alignment.score)
        .then_with(|| right.alignment.identity.total_cmp(&left.alignment.identity))
        .then_with(|| {
            compare_ratio_desc(
                left.alignment.matches,
                left.alignment
                    .matches
                    .saturating_add(left.alignment.substitutions)
                    .saturating_add(left.alignment.deletions),
                right.alignment.matches,
                right
                    .alignment
                    .matches
                    .saturating_add(right.alignment.substitutions)
                    .saturating_add(right.alignment.deletions),
            )
        })
        .then_with(|| {
            compare_ratio_desc(
                left.alignment.score.unsigned_abs(),
                left.alignment.query_length.max(1),
                right.alignment.score.unsigned_abs(),
                right.alignment.query_length.max(1),
            )
        })
        .then_with(|| {
            right
                .alignment
                .seed_evidence
                .nonrepetitive_anchor_count
                .cmp(&left.alignment.seed_evidence.nonrepetitive_anchor_count)
        })
        .then_with(|| covered_length(&right.supported).cmp(&covered_length(&left.supported)))
        .then_with(|| left.alignment.contig_id.cmp(&right.alignment.contig_id))
        .then_with(|| left.alignment_id.cmp(&right.alignment_id))
        .then_with(|| left.canonical_key.cmp(&right.canonical_key))
}

fn compare_ratio_desc(
    left_numerator: u64,
    left_denominator: u64,
    right_numerator: u64,
    right_denominator: u64,
) -> Ordering {
    let left = u128::from(left_numerator) * u128::from(right_denominator.max(1));
    let right = u128::from(right_numerator) * u128::from(left_denominator.max(1));
    right.cmp(&left)
}

/// Assign deterministic IDs to accepted alignments before they enter a
/// mosaic.  IDs are derived from alignment content rather than slice order;
/// identical content receives deterministic ordinal suffixes after canonical
/// sorting.  Existing IDs are intentionally replaced so rerunning this
/// function cannot preserve an order-dependent identifier.
pub fn assign_alignment_ids(alignments: &mut [BaseAlignment]) {
    let keys: Vec<String> = alignments.iter().map(canonical_alignment_key).collect();
    let mut order: Vec<usize> = (0..alignments.len()).collect();
    order.sort_unstable_by(|left, right| {
        keys[*left].cmp(&keys[*right]).then_with(|| left.cmp(right))
    });
    let mut previous_key: Option<&str> = None;
    let mut ordinal = 0_u32;
    for index in order {
        let key = &keys[index];
        if previous_key == Some(key.as_str()) {
            ordinal = ordinal.saturating_add(1);
        } else {
            previous_key = Some(key.as_str());
            ordinal = 0;
        }
        alignments[index].alignment_id = format!(
            "alignment-{}-{:04}",
            crate::provenance::sha256_bytes(key.as_bytes()),
            ordinal
        );
    }
}

fn canonical_alignment_key(alignment: &BaseAlignment) -> String {
    let mut key = String::new();
    let _ = write!(
        key,
        "{}|{}|{}|{:?}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{:?}",
        alignment.plasmid_id,
        alignment.metagenome_id,
        alignment.contig_id,
        alignment.strand,
        alignment.query_length,
        alignment.target_length,
        alignment.origin_crossing,
        alignment.score,
        alignment.matches,
        alignment.substitutions,
        alignment.insertions,
        alignment.deletions,
        alignment.chain_score,
        alignment.target_interval.start,
        alignment.target_interval.end,
        alignment.query_segments.len(),
        alignment.query_segments,
    );
    let _ = write!(
        key,
        "|{}|{:?}|{}|{:?}",
        alignment.cigar, alignment.edit_script, alignment.identity, alignment.seed_evidence,
    );
    key
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
    use crate::trace::model::{AlignmentRole, EditOperation, EditRun, SeedEvidence, Strand};

    fn alignment(contig: &str, start: u64, end: u64) -> BaseAlignment {
        let length = end - start;
        BaseAlignment {
            alignment_id: format!("{contig}:{start}-{end}"),
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
            identity: 1.0,
            seed_evidence: SeedEvidence::default(),
            primary_supported_bases: 0,
            secondary_supported_bases: 0,
            newly_supported_bases: 0,
            role: AlignmentRole::AlternativeMapping,
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
