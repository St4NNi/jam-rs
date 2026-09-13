//! Internal study of local fragment islands inside admitted parent regions.
//!
//! This is a search-scope experiment, not a lossless change and not a default. The parent
//! envelope stays the complete reference and the fallback. An island is the standard envelope of
//! a group of the parent's own anchors; islands inherit the parent's admission and are never
//! admitted again on their own.

use crate::alignment::Alignment;
use crate::trace::{
    AlignmentTask, FragmentEnvelope, RegionAccumulator, RegionKey, SHORT_CONTIG_ENVELOPE_BYTES,
    SeedHit, TraceConfig, TraceError, fragment_envelope,
};

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) enum RegionPolicy {
    /// Every admitted region runs one task on its full envelope.
    #[default]
    Parent,
    /// Split parents run one task per island; any inner-edge contact reruns the full parent.
    Islands,
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
/// merge while their standard envelopes share a base in both target and unwrapped query
/// coordinates. Returns `None` when the parent stays one group or uses the short-contig window.
pub(crate) fn plan_islands(
    parent: &RegionAccumulator,
    parent_envelope: &FragmentEnvelope,
    anchors: &[SeedHit],
    key: RegionKey,
    query_length: u64,
    contig_length: u64,
    config: TraceConfig,
) -> Result<Option<Vec<Island>>, TraceError> {
    if contig_length <= SHORT_CONTIG_ENVELOPE_BYTES || anchors.len() < 2 {
        return Ok(None);
    }
    let mut islands = Vec::<Island>::new();
    for (ordinal, &hit) in anchors.iter().enumerate() {
        let mut region = RegionAccumulator::new(hit);
        let position = parent.support.start + ordinal;
        region.support = position..position + 1;
        let envelope = fragment_envelope(&region, key, query_length, contig_length, config)?;
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
            previous.envelope =
                fragment_envelope(region, key, query_length, contig_length, config)?;
        }
    }
    if islands.len() < 2 {
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
