//! Internal study of local fragment islands inside admitted parent regions.
//!
//! This is a search-scope experiment, not a lossless change and not a default. The parent
//! envelope stays the complete reference and the fallback. An island is the standard envelope of
//! a group of the parent's own anchors; islands inherit the parent's admission and are never
//! admitted again on their own.

use crate::alignment::{Alignment, AlignmentWork};
use crate::trace::{
    AlignmentTask, FragmentEnvelope, RegionAccumulator, RegionKey, SeedHit, TraceConfig,
    TraceError, scoped_fragment_envelope,
};
use serde::Serialize;

/// Seed length of the core anchors the shared index projects every context onto.
const ANCHOR_BASES: u64 = 15;
/// Anchor coordinates are written only for small regions; counts are always written.
const LEDGER_ANCHOR_LIMIT: usize = 16;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) enum RegionPolicy {
    /// Every admitted region runs one task on its full envelope.
    #[default]
    Parent,
    /// Split parents run one task per island; any inner-edge contact reruns the full parent.
    #[cfg_attr(not(any(test, feature = "bench-internals")), allow(dead_code))]
    Islands,
}

impl RegionPolicy {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Parent => "parent",
            Self::Islands => "islands",
        }
    }
}

/// Island window sides that lie strictly inside the parent window. Only these sides can hide an
/// alignment that the parent task could continue.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct InnerEdges {
    pub(crate) query_left: bool,
    pub(crate) query_right: bool,
    pub(crate) target_left: bool,
    pub(crate) target_right: bool,
}

pub(crate) struct Island {
    pub(crate) region: RegionAccumulator,
    pub(crate) envelope: FragmentEnvelope,
    pub(crate) edges: InnerEdges,
}

/// Query interval of an envelope in the unwrapped frame of `anchor`, which the envelope contains.
fn unwrapped_query(envelope: &FragmentEnvelope, anchor: u64, query_length: u64) -> (i128, i128) {
    let offset = if anchor >= envelope.query_start {
        anchor - envelope.query_start
    } else {
        anchor + query_length - envelope.query_start
    };
    let start = i128::from(anchor) - i128::from(offset);
    (start, start + i128::from(envelope.query_span))
}

fn windows_overlap(first: &Island, second: &Island, query_length: u64) -> bool {
    let target = first.envelope.target_start < second.envelope.target_end
        && second.envelope.target_start < first.envelope.target_end;
    let (first_start, first_end) =
        unwrapped_query(&first.envelope, first.region.query_start, query_length);
    let (second_start, second_end) =
        unwrapped_query(&second.envelope, second.region.query_start, query_length);
    target && first_start < second_end && second_start < first_end
}

/// Splits an admitted parent into islands. Each anchor starts its own group; adjacent groups
/// merge while their bounded envelopes share a base in both target and unwrapped query
/// coordinates. Islands never use the whole-contig window of a short contig. Returns `None` when
/// the only island window equals the parent window.
pub(crate) fn plan_islands(
    parent: &RegionAccumulator,
    parent_envelope: &FragmentEnvelope,
    anchors: &[SeedHit],
    key: RegionKey,
    query_length: u64,
    contig_length: u64,
    config: TraceConfig,
) -> Result<Option<Vec<Island>>, TraceError> {
    if anchors.is_empty() {
        return Ok(None);
    }
    let bounded = |region: &RegionAccumulator| {
        scoped_fragment_envelope(region, key, query_length, contig_length, config, false)
    };
    let mut islands = Vec::<Island>::new();
    for (ordinal, &hit) in anchors.iter().enumerate() {
        let mut region = RegionAccumulator::new(hit);
        let position = parent.support.start + ordinal;
        region.support = position..position + 1;
        let envelope = bounded(&region)?;
        islands.push(Island {
            region,
            envelope,
            edges: InnerEdges::default(),
        });
        while islands.len() >= 2
            && windows_overlap(
                &islands[islands.len() - 2],
                &islands[islands.len() - 1],
                query_length,
            )
        {
            let last = islands.pop().expect("two islands");
            let previous = islands.last_mut().expect("two islands");
            let region = &mut previous.region;
            region.query_end = last.region.query_end;
            region.target_end = last.region.target_end;
            region.diagonal_min = region.diagonal_min.min(last.region.diagonal_min);
            region.diagonal_max = region.diagonal_max.max(last.region.diagonal_max);
            region.hits = region.hits.saturating_add(last.region.hits);
            region.support.end = last.region.support.end;
            previous.envelope = bounded(region)?;
        }
    }
    if islands.len() == 1 && islands[0].envelope == *parent_envelope {
        return Ok(None);
    }
    let (parent_start, parent_end) =
        unwrapped_query(parent_envelope, parent.query_start, query_length);
    let whole_query = parent_envelope.query_span >= query_length;
    for island in &mut islands {
        let (start, end) =
            unwrapped_query(&island.envelope, island.region.query_start, query_length);
        let partial = whole_query && island.envelope.query_span < parent_envelope.query_span;
        island.edges = InnerEdges {
            query_left: partial || start > parent_start,
            query_right: partial || end < parent_end,
            target_left: island.envelope.target_start > parent_envelope.target_start,
            target_right: island.envelope.target_end < parent_envelope.target_end,
        };
    }
    Ok(Some(islands))
}

/// True when an island alignment reaches within `margin` bases of an inner window side, in query
/// window or forward target coordinates. The parent task could continue such an alignment.
pub(crate) fn touches_inner_edge(task: &AlignmentTask, alignment: &Alignment, margin: u64) -> bool {
    let Some(edges) = task.island else {
        return false;
    };
    let query = alignment.query_interval;
    let target = alignment.target_interval;
    (edges.query_left && query.start < margin)
        || (edges.query_right && query.end.saturating_add(margin) > task.query_span)
        || (edges.target_left && target.start < task.target_start.saturating_add(margin))
        || (edges.target_right && target.end.saturating_add(margin) > task.target_end)
}

#[derive(Clone, Debug, Serialize)]
pub(crate) struct WindowLedger {
    pub(crate) query_start: u64,
    pub(crate) query_span: u64,
    pub(crate) target_start: u64,
    pub(crate) target_end: u64,
    pub(crate) diagonal_offset: i64,
    pub(crate) anchors: u32,
}

impl WindowLedger {
    pub(crate) fn new(envelope: &FragmentEnvelope, anchors: u32) -> Self {
        Self {
            query_start: envelope.query_start,
            query_span: envelope.query_span,
            target_start: envelope.target_start,
            target_end: envelope.target_end,
            diagonal_offset: envelope.diagonal_offset,
            anchors,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub(crate) struct AlignmentLedger {
    pub(crate) query_start: u64,
    pub(crate) query_end: u64,
    pub(crate) target_start: u64,
    pub(crate) target_end: u64,
    pub(crate) score: i32,
    pub(crate) identity: f64,
    /// Parent anchors whose query start lies inside the aligned query interval.
    pub(crate) anchors_inside: u32,
    /// Query bases from the aligned interval to the nearest parent anchor; zero when inside.
    pub(crate) nearest_anchor_distance: u64,
}

impl AlignmentLedger {
    pub(crate) fn new(
        alignment: &Alignment,
        window_start: u64,
        anchors: &[SeedHit],
        query_length: u64,
    ) -> Self {
        let interval = alignment.query_interval;
        let mut inside = 0;
        let mut nearest = u64::MAX;
        for anchor in anchors {
            let position = if anchor.query >= window_start {
                anchor.query - window_start
            } else {
                anchor.query + query_length - window_start
            };
            if position >= interval.start && position < interval.end {
                inside += 1;
                nearest = 0;
            } else {
                let distance = if position < interval.start {
                    interval.start - position
                } else {
                    position + 1 - interval.end
                };
                nearest = nearest.min(distance);
            }
        }
        Self {
            query_start: interval.start,
            query_end: interval.end,
            target_start: alignment.target_interval.start,
            target_end: alignment.target_interval.end,
            score: alignment.score,
            identity: alignment.identity(),
            anchors_inside: inside,
            nearest_anchor_distance: nearest,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub(crate) struct TaskLedger {
    /// `parent`, `island` or `fallback`.
    pub(crate) role: &'static str,
    pub(crate) window: WindowLedger,
    pub(crate) band_cells: u64,
    pub(crate) work: AlignmentWork,
    pub(crate) core: Option<AlignmentLedger>,
    pub(crate) selected: Option<AlignmentLedger>,
    pub(crate) accepted: bool,
    pub(crate) contact: bool,
}

#[derive(Clone, Debug, Default, Serialize)]
pub(crate) struct SupportLedger {
    pub(crate) anchors: u32,
    pub(crate) anchors_with_21: u32,
    pub(crate) anchors_with_31: u32,
    pub(crate) anchor_query: Vec<u64>,
    pub(crate) anchor_target: Vec<u64>,
    pub(crate) union_query_bases: u64,
    pub(crate) union_target_bases: u64,
    pub(crate) query_span: u64,
    pub(crate) target_span: u64,
    pub(crate) max_query_gap: u64,
    pub(crate) max_target_gap: u64,
    pub(crate) max_gap_difference: u64,
    pub(crate) diagonal_spread: u64,
}

impl SupportLedger {
    /// Summarizes distinct core anchor pairs. `contexts` holds one bit per anchor for 21 (2) and
    /// 31 (4) contexts; contexts are evidence about the same pair, not further hits.
    pub(crate) fn new(anchors: &[SeedHit], contexts: &[u8]) -> Self {
        let mut ledger = Self {
            anchors: u32::try_from(anchors.len()).unwrap_or(u32::MAX),
            ..Self::default()
        };
        let (Some(first), Some(last)) = (anchors.first(), anchors.last()) else {
            return ledger;
        };
        ledger.anchors_with_21 = contexts.iter().filter(|&&mask| mask & 2 != 0).count() as u32;
        ledger.anchors_with_31 = contexts.iter().filter(|&&mask| mask & 4 != 0).count() as u32;
        if anchors.len() <= LEDGER_ANCHOR_LIMIT {
            ledger.anchor_query = anchors.iter().map(|hit| hit.query).collect();
            ledger.anchor_target = anchors.iter().map(|hit| hit.target).collect();
        }
        ledger.query_span = last.query + ANCHOR_BASES - first.query;
        ledger.target_span = last.target + ANCHOR_BASES - first.target;
        let (mut diagonal_min, mut diagonal_max) = (first.diagonal, first.diagonal);
        let mut query_end = first.query;
        let mut target_end = first.target;
        for (index, hit) in anchors.iter().enumerate() {
            ledger.union_query_bases +=
                ANCHOR_BASES - query_end.saturating_sub(hit.query).min(ANCHOR_BASES);
            ledger.union_target_bases +=
                ANCHOR_BASES - target_end.saturating_sub(hit.target).min(ANCHOR_BASES);
            if index > 0 {
                let previous = anchors[index - 1];
                ledger.max_query_gap = ledger
                    .max_query_gap
                    .max(hit.query.saturating_sub(previous.query + ANCHOR_BASES));
                ledger.max_target_gap = ledger
                    .max_target_gap
                    .max(hit.target.saturating_sub(previous.target + ANCHOR_BASES));
                ledger.max_gap_difference = ledger
                    .max_gap_difference
                    .max((hit.diagonal - previous.diagonal).unsigned_abs() as u64);
            }
            query_end = hit.query + ANCHOR_BASES;
            target_end = hit.target + ANCHOR_BASES;
            diagonal_min = diagonal_min.min(hit.diagonal);
            diagonal_max = diagonal_max.max(hit.diagonal);
        }
        ledger.diagonal_spread = (diagonal_max - diagonal_min).unsigned_abs() as u64;
        ledger
    }
}

#[derive(Clone, Debug, Serialize)]
pub(crate) struct ParentLedger {
    pub(crate) query_id: String,
    pub(crate) query_length: u64,
    pub(crate) circular: bool,
    pub(crate) policy: &'static str,
    pub(crate) parent: u32,
    pub(crate) metagenome_id: u32,
    pub(crate) contig_id: u32,
    pub(crate) contig_length: u64,
    pub(crate) strand: crate::alignment::Strand,
    pub(crate) support: SupportLedger,
    pub(crate) parent_window: WindowLedger,
    pub(crate) islands: Vec<WindowLedger>,
    pub(crate) fallback: bool,
    pub(crate) tasks: Vec<TaskLedger>,
}
