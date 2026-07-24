//! Aligning one chain: one banded local alignment around the chain's diagonals and the query
//! bases beyond it. The window starts close around the chain and widens while the alignment runs
//! into its edge, so a short element shared with a long sequence does not cost an alignment of
//! all of it.
//!
//! Query positions of the results count from the query start without wrapping: for a circular
//! query an alignment over the origin ends beyond the query length.

use super::anchors::QueryTable;
use super::chains::Chain;
use super::{
    Alignment, AlignmentConfig, AlignmentError, AlignmentInput, AlignmentWorkspace, EditOperation,
    EditRun, Interval, complement,
};

/// The band that aligns chains and walks: scores, width and the drift of anchors it follows.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BandConfig {
    /// Scores, the narrowest band and the resident cell limit; the diagonal is set per chain.
    pub alignment: AlignmentConfig,
    /// Largest drift of a chain's anchors off its first diagonal; drifting chains are split.
    pub max_drift: u32,
}

impl Default for BandConfig {
    fn default() -> Self {
        Self {
            alignment: AlignmentConfig::default(),
            max_drift: 1_024,
        }
    }
}

/// Diagonals and bases kept beyond a drift, so that an alignment that follows the drift still
/// has room for the indels around it.
const DRIFT_MARGIN: u32 = 64;

impl BandConfig {
    /// The widest band a chain or walk alignment uses: the configured band or the drift margin.
    pub fn max_band(self) -> u32 {
        self.alignment.band_width.max(self.max_drift + DRIFT_MARGIN)
    }

    /// The largest drift and its margin: how far an alignment may leave the diagonal of the
    /// anchors or walk it follows.
    pub(super) fn drift_margin(self) -> u64 {
        u64::from(self.max_drift + DRIFT_MARGIN)
    }
}

/// The query bases beyond a chain that its first alignment may use on each side.
const FIRST_EXTENSION: u64 = 128;

/// Aligns one chain of `window` (the contig bases `offset..offset + window.len()`), allowing the
/// alignment to reach `reach` query bases beyond the chain on each side. `None` when no positive
/// alignment exists.
pub fn align_chain(
    chain: &Chain,
    window: &[u8],
    offset: u64,
    table: &QueryTable,
    config: BandConfig,
    workspace: &mut AlignmentWorkspace,
    reach: u64,
) -> Result<Option<Alignment>, AlignmentError> {
    if let Some(exact) = exact_alignment(chain, window, offset, table, config)? {
        return Ok(Some(exact));
    }
    let geometry = ChainGeometry::new(chain, window.len() as u64, table, config, reach);
    let a = u64::from(table.anchors().a);
    let length = u64::from(table.query_length());
    let (_, widest_left, widest_span, _, _) = geometry.window_at(geometry.extension);
    let widest_right = widest_span - widest_left - geometry.span.min(widest_span);
    // A window that reaches at most the margin beyond the chain always counts as run into: the
    // chain itself starts within the margin. So the first window reaches past the margin.
    let mut narrow = FIRST_EXTENSION;
    while narrow <= geometry.margin {
        narrow = narrow.saturating_mul(4);
    }
    let mut narrow = narrow.min(geometry.extension);
    loop {
        let (query_start, left, query_span, contig_low, contig_high) = geometry.window_at(narrow);
        let (forward_low, forward_high) = if chain.strand.is_reverse() {
            (
                window.len() as u64 - contig_high,
                window.len() as u64 - contig_low,
            )
        } else {
            (contig_low, contig_high)
        };
        let scores = AlignmentConfig {
            band_width: geometry.band,
            diagonal_offset: (geometry.first_contig - contig_low) as i64 - left as i64
                + geometry.centre,
            ..config.alignment
        };
        let mut bases = Vec::new();
        table.copy_bases(query_start, query_span, &mut bases)?;
        let input = AlignmentInput {
            query: &bases,
            query_offset: query_start,
            contig: &window[forward_low as usize..forward_high as usize],
            contig_offset: offset + forward_low,
        };
        let alignment = match workspace.align(input, chain.strand, scores) {
            Ok(alignment) => alignment,
            Err(AlignmentError::NoAlignment) => return Ok(None),
            Err(error) => return Err(error),
        };
        let (aligned_start, aligned_end) = (
            alignment.query_interval.start - query_start,
            alignment.query_interval.end - query_start,
        );
        let right = query_span - left - geometry.span.min(query_span);
        let runs_into_edge = (left < widest_left && aligned_start <= geometry.margin)
            || (right < widest_right && query_span - aligned_end <= geometry.margin);
        if runs_into_edge && narrow < geometry.extension {
            narrow = narrow.saturating_mul(4).min(geometry.extension);
            continue;
        }
        // An alignment of the whole query that stops within an anchor length of both ends is
        // completed through them.
        if query_span == length && aligned_start <= a && query_span - aligned_end <= a {
            return workspace
                .complete_endpoints(alignment, input, a as usize, scores)
                .map(Some);
        }
        return Ok(Some(alignment));
    }
}

/// The alignment of a linear chain whose anchors all lie on one diagonal along which the query
/// and the window are equal, without the kernel. When both are copies of one short repeat, the
/// diagonal of the chain is taken as it is.
fn exact_alignment(
    chain: &Chain,
    window: &[u8],
    offset: u64,
    table: &QueryTable,
    config: BandConfig,
) -> Result<Option<Alignment>, AlignmentError> {
    let query = table.query();
    let length = query.len() as u64;
    let (first_contig, first_query) = chain.anchors[0];
    let diagonal = i64::from(first_contig) - i64::from(first_query);
    if table.circular()
        || !chain
            .anchors
            .iter()
            .all(|&(t, q)| i64::from(t) - i64::from(q) == diagonal)
    {
        return Ok(None);
    }
    let (q, t) = ((-diagonal).max(0) as u64, diagonal.max(0) as u64);
    let paired = (length - q).min(window.len() as u64 - t);
    let equal = (0..paired as usize).all(|at| {
        let contig = if chain.strand.is_reverse() {
            complement(window[window.len() - 1 - t as usize - at])
        } else {
            window[t as usize + at]
        };
        query[q as usize + at].eq_ignore_ascii_case(&contig)
    }) || table.periodic_pair(window).is_some();
    if paired == 0 || !equal {
        return Ok(None);
    }
    let (low, high) = if chain.strand.is_reverse() {
        (window.len() as u64 - t - paired, window.len() as u64 - t)
    } else {
        (t, t + paired)
    };
    let alignment = Alignment {
        score: config
            .alignment
            .match_score
            .checked_mul(paired as i32)
            .ok_or(AlignmentError::LengthOverflow)?,
        strand: chain.strand,
        query_interval: Interval::new(q, q + paired)?,
        contig_interval: Interval::new(offset + low, offset + high)?,
        runs: vec![EditRun {
            operation: EditOperation::Match,
            length: u32::try_from(paired).map_err(|_| AlignmentError::RunTooLong)?,
        }],
    };
    Ok(Some(alignment))
}

/// The band and the window choices of one chain.
struct ChainGeometry {
    circular: bool,
    length: u64,
    a: u64,
    window_length: u64,
    first_contig: u64,
    first_query: u64,
    last_contig: u64,
    /// Query bases from the first anchor's start to the last anchor's end.
    span: u64,
    slack: u64,
    margin: u64,
    centre: i64,
    band: u32,
    extension: u64,
}

impl ChainGeometry {
    fn new(
        chain: &Chain,
        window_length: u64,
        table: &QueryTable,
        config: BandConfig,
        reach: u64,
    ) -> Self {
        let length = u64::from(table.query_length());
        let a = u64::from(table.anchors().a);
        let (first_contig, first_query) = chain.anchors[0];
        let first_contig = u64::from(first_contig);
        // The lowest and highest diagonal of the anchors, counted from the first anchor's.
        let (low, high) = chain.anchors.iter().zip(&chain.offsets).fold(
            (0i64, 0i64),
            |(low, high), (&(contig, _), &offset)| {
                let diagonal = i64::from(contig) - first_contig as i64 - i64::from(offset);
                (low.min(diagonal), high.max(diagonal))
            },
        );
        let drift = low.unsigned_abs().max(high.unsigned_abs());
        let margin = drift + u64::from(DRIFT_MARGIN);
        let band_width = config.alignment.band_width;
        // A band wider than the configured one follows the middle of the anchor diagonals: a
        // chain that drifts to one side needs half its drift range on each side, not all of it.
        let (centre, band) = if margin > u64::from(band_width) {
            let centre = low + (high - low) / 2;
            let half = (high - centre).max(centre - low).unsigned_abs() + u64::from(DRIFT_MARGIN);
            (centre, band_width.max(half as u32))
        } else {
            (0, band_width.max(margin as u32))
        };
        Self {
            circular: table.circular(),
            length,
            a,
            window_length,
            first_contig,
            first_query: u64::from(first_query),
            last_contig: u64::from(chain.anchors[chain.anchors.len() - 1].0),
            span: u64::from(*chain.offsets.last().unwrap_or(&0)) + a,
            slack: config.drift_margin(),
            margin,
            centre,
            band,
            extension: reach.min(length),
        }
    }

    /// The query start, the query bases before the chain, the query length and the contig range
    /// of an alignment that may extend by `extension` query bases beyond the chain.
    fn window_at(&self, extension: u64) -> (u64, u64, u64, u64, u64) {
        let (length, a, span) = (self.length, self.a, self.span);
        let after = self.window_length.saturating_sub(self.last_contig + a);
        let (query_start, left, query_span) = if self.circular {
            // The query bases beyond the chain are split between its sides, but a side without
            // window bases gives its share to the other, so a whole linear copy is not cut short.
            let room = length.saturating_sub(span).min(2 * extension);
            let right = (room - room / 2).min(after);
            let left = (room - right).min(self.first_contig);
            let right = room - left;
            // Query bases beyond the contig bases that the window has on a side can only be a
            // long gap, so the query window does not reach further than those and the drift
            // allowance.
            let (left, right) = (
                left.min(self.first_contig + self.slack),
                right.min(after + self.slack),
            );
            (
                (self.first_query + length - left) % length,
                left,
                span.min(length) + left + right,
            )
        } else {
            let left = extension
                .min(self.first_query)
                .min(self.first_contig + self.slack);
            let end = (self.first_query + span + extension.min(after + self.slack)).min(length);
            (
                self.first_query - left,
                left,
                end - (self.first_query - left),
            )
        };
        let contig_low = self.first_contig.saturating_sub(left + self.margin);
        let contig_high =
            (self.last_contig + a + (query_span - left - span.min(query_span)) + self.margin)
                .min(self.window_length);
        (query_start, left, query_span, contig_low, contig_high)
    }
}
