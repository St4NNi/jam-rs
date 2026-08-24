//! Bounded, gap-directed trace cascade planning.
//!
//! This module contains the state machine used by the positional runner.  It
//! deliberately does not know how seeds or sequence bytes are stored: the
//! runner supplies stage metrics and keeps the accepted alignment evidence.
//! A stage that makes no progress is terminal unless an explicit resource
//! budget still permits a caller-selected optional stage.

use crate::archive::SeedSchemeDescriptor;
use crate::trace::model::{BaseInterval, TraceStageMetrics};

/// Descriptor algorithm identifier reserved for variable-length exact Gear
/// fragment records. A caller must still verify returned sequence bytes.
pub const GEAR_FRAGMENT_ALGORITHM_ID: u32 = 0x4745_4152;

/// Return the advertised exact-fragment scheme, if the archive has one. The
/// key encoding flag is checked here so a conventional fixed-k seed section
/// can never accidentally activate the fragment path.
#[must_use]
pub fn advertised_exact_gear_scheme(
    descriptors: &[SeedSchemeDescriptor],
) -> Option<SeedSchemeDescriptor> {
    descriptors
        .iter()
        .copied()
        .filter(|descriptor| {
            descriptor.algorithm_id == GEAR_FRAGMENT_ALGORITHM_ID
                && descriptor.key_encoding == 2
                && descriptor.occurrence_encoding == 2
        })
        .min_by_key(|descriptor| descriptor.scheme_id)
}

/// Ordered positional stages.  The first stage is an optional exact Gear
/// fragment fast path; it is present only when the archive advertises one.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CascadeStageKind {
    ExactGearFragments,
    SparseK31,
    DenseGapK31,
    MixedGapK31K21,
    AlternativeRescue,
}

impl CascadeStageKind {
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::ExactGearFragments => "exact_gear_fragments",
            Self::SparseK31 => "sparse_k31",
            Self::DenseGapK31 => "dense_gap_k31",
            Self::MixedGapK31K21 => "mixed_gap_k31_k21",
            Self::AlternativeRescue => "alternative_rescue",
        }
    }

    #[must_use]
    pub const fn optional(self) -> bool {
        matches!(self, Self::ExactGearFragments | Self::AlternativeRescue)
    }
}

/// A planned stage with explicit bounded work.  The runner may further lower
/// each value from candidate-level configuration before executing it.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CascadeStage {
    pub ordinal: u8,
    pub kind: CascadeStageKind,
    pub max_seed_pages: u64,
    pub max_sequence_blocks: u64,
    pub max_alignment_windows: u64,
}

/// Stop reason retained for diagnostics and tests.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CascadeStopReason {
    Complete,
    NoGain,
    NoEligibleSeedPages,
    ResourceBudget,
    StageLimit,
}

/// The state shared by all cascade stages.
#[derive(Clone, Debug)]
pub struct CascadeState {
    query_length: u64,
    accepted_intervals: Vec<BaseInterval>,
    pub stages: Vec<TraceStageMetrics>,
    pub stop_reason: Option<CascadeStopReason>,
}

impl CascadeState {
    #[must_use]
    pub fn new(query_length: u64) -> Self {
        Self {
            query_length,
            accepted_intervals: Vec::new(),
            stages: Vec::new(),
            stop_reason: None,
        }
    }

    #[must_use]
    pub const fn query_length(&self) -> u64 {
        self.query_length
    }

    #[must_use]
    pub fn accepted_intervals(&self) -> &[BaseInterval] {
        &self.accepted_intervals
    }

    /// Replace the accepted interval union with a new mosaic projection.
    /// Intervals are clipped, sorted, and merged so overlapping contigs or
    /// duplicate alignments cannot inflate the supported-base count.
    pub fn retain_supported_intervals(
        &mut self,
        intervals: impl IntoIterator<Item = BaseInterval>,
    ) {
        let mut intervals = intervals
            .into_iter()
            .filter_map(|interval| {
                let start = interval.start.min(self.query_length);
                let end = interval.end.min(self.query_length);
                (start < end).then_some(BaseInterval { start, end })
            })
            .collect::<Vec<_>>();
        intervals.sort_by_key(|interval| (interval.start, interval.end));
        let mut merged: Vec<BaseInterval> = Vec::with_capacity(intervals.len());
        for interval in intervals {
            if let Some(last) = merged.last_mut()
                && interval.start <= last.end
            {
                last.end = last.end.max(interval.end);
                continue;
            }
            merged.push(interval);
        }
        self.accepted_intervals = merged;
    }

    #[must_use]
    pub fn supported_bases(&self) -> u64 {
        self.accepted_intervals
            .iter()
            .map(|interval| interval.end.saturating_sub(interval.start))
            .fold(0, u64::saturating_add)
    }

    /// Record stage metrics and return the number of newly supported bases.
    /// The caller supplies the accepted mosaic intervals after the stage has
    /// completed; previously accepted evidence is never discarded.
    pub fn record_stage(
        &mut self,
        mut metrics: TraceStageMetrics,
        supported_intervals: impl IntoIterator<Item = BaseInterval>,
    ) -> u64 {
        let before = self.supported_bases();
        let mut union = self.accepted_intervals.clone();
        union.extend(supported_intervals);
        self.retain_supported_intervals(union);
        let after = self.supported_bases();
        let gain = after.saturating_sub(before);
        metrics.new_query_bases_supported = gain;
        self.stages.push(metrics);
        gain
    }

    /// Apply the ordered stop rules after one stage.  Completion always wins;
    /// a zero-gain stage stops the default cascade before optional methods can
    /// repeat the same work unless the caller explicitly enables them.
    pub fn stop_after_stage(
        &mut self,
        gain: u64,
        query_complete: bool,
        eligible_seed_pages: bool,
        budget_available: bool,
    ) -> CascadeStopReason {
        let reason = if query_complete {
            CascadeStopReason::Complete
        } else if !eligible_seed_pages {
            CascadeStopReason::NoEligibleSeedPages
        } else if !budget_available {
            CascadeStopReason::ResourceBudget
        } else if gain == 0 {
            CascadeStopReason::NoGain
        } else {
            return CascadeStopReason::StageLimit;
        };
        self.stop_reason = Some(reason);
        reason
    }

    #[must_use]
    pub fn unresolved_gaps(&self) -> Vec<BaseInterval> {
        let mut gaps = Vec::new();
        let mut cursor = 0;
        for interval in &self.accepted_intervals {
            if cursor < interval.start {
                gaps.push(BaseInterval {
                    start: cursor,
                    end: interval.start,
                });
            }
            cursor = cursor.max(interval.end);
        }
        if cursor < self.query_length {
            gaps.push(BaseInterval {
                start: cursor,
                end: self.query_length,
            });
        }
        gaps
    }
}

/// Construct the default ordered plan.  Gear and alternative schemes are
/// opt-in based solely on archive descriptors; no target-wide index is built
/// when either capability is absent.
#[must_use]
pub fn default_stage_plan(
    has_gear_fragments: bool,
    has_dense_k31: bool,
    has_k21: bool,
    has_alternative_scheme: bool,
    max_seed_pages: u64,
    max_sequence_blocks: u64,
    max_alignment_windows: u64,
) -> Vec<CascadeStage> {
    let mut kinds = Vec::new();
    if has_gear_fragments {
        kinds.push(CascadeStageKind::ExactGearFragments);
    }
    kinds.push(CascadeStageKind::SparseK31);
    if has_dense_k31 {
        kinds.push(CascadeStageKind::DenseGapK31);
    }
    if has_k21 {
        kinds.push(CascadeStageKind::MixedGapK31K21);
    }
    if has_alternative_scheme {
        kinds.push(CascadeStageKind::AlternativeRescue);
    }
    kinds
        .into_iter()
        .enumerate()
        .map(|(index, kind)| CascadeStage {
            ordinal: u8::try_from(index + 1).unwrap_or(u8::MAX),
            kind,
            max_seed_pages,
            max_sequence_blocks,
            max_alignment_windows,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn interval(start: u64, end: u64) -> BaseInterval {
        BaseInterval { start, end }
    }

    #[test]
    fn plan_is_ordered_and_optional_methods_are_descriptor_gated() {
        let plan = default_stage_plan(true, true, true, true, 10, 20, 30);
        assert_eq!(
            plan.iter().map(|stage| stage.kind).collect::<Vec<_>>(),
            vec![
                CascadeStageKind::ExactGearFragments,
                CascadeStageKind::SparseK31,
                CascadeStageKind::DenseGapK31,
                CascadeStageKind::MixedGapK31K21,
                CascadeStageKind::AlternativeRescue,
            ]
        );
        assert!(plan[0].kind.optional());
        assert!(!plan[1].kind.optional());
    }

    #[test]
    fn accepted_evidence_is_union_and_no_gain_stops() {
        let mut state = CascadeState::new(100);
        let first = TraceStageMetrics {
            stage: 1,
            name: "primary".to_string(),
            ..TraceStageMetrics::default()
        };
        assert_eq!(state.record_stage(first, [interval(10, 40)]), 30);
        let second = TraceStageMetrics {
            stage: 2,
            name: "rescue".to_string(),
            ..TraceStageMetrics::default()
        };
        assert_eq!(state.record_stage(second, [interval(20, 40)]), 0);
        assert_eq!(
            state.stop_after_stage(0, false, true, true),
            CascadeStopReason::NoGain
        );
        assert_eq!(state.accepted_intervals(), &[interval(10, 40)]);
    }

    #[test]
    fn gaps_are_query_coordinate_intervals() {
        let mut state = CascadeState::new(100);
        state.retain_supported_intervals([interval(10, 20), interval(18, 50), interval(75, 90)]);
        assert_eq!(
            state.unresolved_gaps(),
            vec![interval(0, 10), interval(50, 75), interval(90, 100)]
        );
    }

    #[test]
    fn exact_gear_fast_path_requires_variable_fragment_encoding() {
        let descriptor = SeedSchemeDescriptor {
            scheme_id: 7,
            algorithm_id: GEAR_FRAGMENT_ALGORITHM_ID,
            span: 0,
            informative_bases: 0,
            density_parameter: 0,
            bucket_bits: 8,
            key_encoding: 2,
            occurrence_encoding: 2,
            flags: 0,
        };
        assert_eq!(
            advertised_exact_gear_scheme(&[descriptor]),
            Some(descriptor)
        );
        assert!(
            advertised_exact_gear_scheme(&[SeedSchemeDescriptor {
                key_encoding: 1,
                ..descriptor
            }])
            .is_none()
        );
    }
}
