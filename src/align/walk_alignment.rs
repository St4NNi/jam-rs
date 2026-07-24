//! Alignments of dense walks: the ungapped alignment of one walk, and the probes that decide
//! whether a walk is aligned again with gaps or lies inside an alignment already made.

use super::anchors::{QueryTable, kmers};
use super::chain_alignment::BandConfig;
use super::cigar::{append_runs, counts};
use super::dense::Walk;
use super::{Alignment, AlignmentError, EditOperation, EditRun, Interval};

/// The ungapped alignment of one walk. `oriented` is the window on the walk's strand (the
/// reverse complement for the reverse strand); the window starts `offset` bases into the contig.
pub fn walk_alignment(
    walk: &Walk,
    oriented: &[u8],
    offset: u64,
    table: &QueryTable,
    config: BandConfig,
) -> Result<Alignment, AlignmentError> {
    let query = table.query();
    let length = query.len();
    let mut runs: Vec<EditRun> = Vec::new();
    let (mut matches, mut substitutions) = (0i64, 0i64);
    let (contig, mut at) = (walk.contig as usize, 0usize);
    while at < walk.length as usize {
        // The longest piece in which the query does not wrap.
        let from = (walk.query as usize + at) % length;
        let piece = (walk.length as usize - at).min(length - from);
        let (a, b) = (
            &oriented[contig + at..contig + at + piece],
            &query[from..from + piece],
        );
        let mut step = 0;
        while step < piece {
            let equal = step;
            loop {
                if step + 8 <= piece && a[step..step + 8] == b[step..step + 8] {
                    step += 8;
                } else if step < piece && a[step] == b[step] {
                    step += 1;
                } else {
                    break;
                }
            }
            if step > equal {
                append_runs(&mut runs, [run(EditOperation::Match, step - equal)])?;
                matches += (step - equal) as i64;
            }
            if step < piece {
                append_runs(&mut runs, [run(EditOperation::Substitution, 1)])?;
                (substitutions, step) = (substitutions + 1, step + 1);
            }
        }
        at += piece;
    }
    let score = i64::from(config.alignment.match_score) * matches
        + i64::from(config.alignment.mismatch_score) * substitutions;
    oriented_alignment(runs, score, walk, oriented.len(), offset)
}

fn run(operation: EditOperation, length: usize) -> EditRun {
    EditRun {
        operation,
        length: length as u32,
    }
}

/// The alignment of `runs` starting at the oriented window and query positions of `first`, in
/// contig positions. Query positions run on over the origin of a circular query.
pub(super) fn oriented_alignment(
    runs: Vec<EditRun>,
    score: i64,
    first: &Walk,
    window_length: usize,
    offset: u64,
) -> Result<Alignment, AlignmentError> {
    let counts = counts(&runs);
    let (contig_span, query_span) = (counts.contig_bases(), counts.query_bases());
    let low = u64::from(first.contig);
    let (start, end) = if first.strand.is_reverse() {
        let window_length = window_length as u64;
        (window_length - low - contig_span, window_length - low)
    } else {
        (low, low + contig_span)
    };
    let query_start = u64::from(first.query);
    Ok(Alignment {
        score: i32::try_from(score).map_err(|_| AlignmentError::LengthOverflow)?,
        strand: first.strand,
        query_interval: Interval::new(query_start, query_start + query_span)?,
        contig_interval: Interval::new(offset + start, offset + end)?,
        runs,
    })
}

/// Whether an unpaired end can reach a conserved block on a different diagonal: whether an
/// exact k-mer within `max_drift` bases beyond the end of the query pairs with one within
/// `max_drift` bases beyond the end of the contig at another distance from the end. Differences
/// next to the end, such as a substitution beside an insertion, do not hide the block. `codes`
/// is reused scratch of at most `max_drift + 1` codes.
pub fn gap_continuation(
    walk: &Walk,
    oriented: &[u8],
    table: &QueryTable,
    config: BandConfig,
    codes: &mut Vec<(u64, u32)>,
) -> Result<bool, AlignmentError> {
    let query = table.query();
    let a = usize::from(table.anchors().a);
    let reach = config.max_drift as usize + a;
    for right in [false, true] {
        let contig = walk.contig as usize + if right { walk.length as usize } else { 0 };
        let edge = walk.query as usize + if right { walk.length as usize } else { 0 };
        let room = if table.circular() {
            query.len() - walk.length as usize
        } else if right {
            query.len() - edge
        } else {
            edge
        };
        let query_bases = room.min(reach);
        let contig_bases = if right {
            oriented.len().saturating_sub(contig)
        } else {
            contig
        }
        .min(reach);
        if query_bases < a || contig_bases < a {
            continue;
        }
        // The bases beyond the end, in query order; a k-mer's distance from the end counts the
        // bases between them.
        let query_from = if right {
            edge
        } else {
            edge + query.len() - query_bases
        };
        let mut query_side = Vec::with_capacity(query_bases);
        let from = (query_from % query.len()) as u64;
        table.copy_bases(from, query_bases as u64, &mut query_side)?;
        let contig_side = if right {
            &oriented[contig..contig + contig_bases]
        } else {
            &oriented[contig - contig_bases..contig]
        };
        let distance = |at: usize, bases: usize| if right { at } else { bases - a - at };
        // The end continues on its own diagonal, as at a circular query's whole length.
        if query_side[distance(0, query_bases)..][..a]
            == contig_side[distance(0, contig_bases)..][..a]
        {
            continue;
        }
        codes.clear();
        kmers(&query_side, table.anchors().a, |at, code, _| {
            codes.push((code, distance(at, query_bases) as u32))
        });
        codes.sort_unstable();
        let mut found = false;
        kmers(contig_side, table.anchors().a, |at, code, _| {
            let at = distance(at, contig_bases) as u32;
            // Pairs of one code are ordered by distance, so another distance is at an end.
            let first = codes.partition_point(|pair| pair.0 < code);
            let last = codes.partition_point(|pair| pair.0 <= code);
            found |= first < last && (codes[first].1 != at || codes[last - 1].1 != at);
        });
        if found {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Whether `alignment` pairs the query positions of `walk` within its contig positions, on its
/// strand. The walk comes from an alignment window of `window_length` bases that starts `offset`
/// bases into the contig.
pub fn holds_walk(
    alignment: &Alignment,
    walk: &Walk,
    window_length: usize,
    offset: u64,
    query_length: u64,
) -> bool {
    if (alignment.strand == super::Strand::Reverse) != walk.strand.is_reverse() {
        return false;
    }
    let (low, high) = (u64::from(walk.contig), u64::from(walk.contig + walk.length));
    let (low, high) = if walk.strand.is_reverse() {
        (window_length as u64 - high, window_length as u64 - low)
    } else {
        (low, high)
    };
    let contig = alignment.contig_interval;
    if offset + low < contig.start || offset + high > contig.end {
        return false;
    }
    // The query positions of the alignment may run on over the origin.
    let walk_query = Interval {
        start: u64::from(walk.query),
        end: u64::from(walk.query + walk.length),
    };
    walk_query.wrap(query_length).all(|piece| {
        alignment
            .query_interval
            .wrap(query_length)
            .any(|segment| segment.start <= piece.start && piece.end <= segment.end)
    })
}
