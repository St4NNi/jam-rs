//! Alignments of dense walks: the ungapped alignment of one walk, and the probes that decide
//! whether a walk is aligned again with gaps or lies inside an alignment already made.

use super::anchors::{QueryTable, kmers};
use super::chain_alignment::{BandConfig, query_window, strand};
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
                append(&mut runs, &[run(EditOperation::Match, step - equal)]);
                matches += (step - equal) as i64;
            }
            if step < piece {
                append(&mut runs, &[run(EditOperation::Substitution, 1)]);
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

/// Appends `more` to `runs`, merging runs of one operation.
pub(super) fn append(runs: &mut Vec<EditRun>, more: &[EditRun]) {
    for run in more {
        match runs.last_mut() {
            Some(last) if last.operation == run.operation => last.length += run.length,
            _ => runs.push(*run),
        }
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
    let count = |operations: &[EditOperation]| {
        runs.iter()
            .filter(|run| operations.contains(&run.operation))
            .map(|run| u64::from(run.length))
            .sum::<u64>()
    };
    let paired = [EditOperation::Match, EditOperation::Substitution];
    let contig_span = count(&paired) + count(&[EditOperation::Deletion]);
    let query_span = count(&paired) + count(&[EditOperation::Insertion]);
    let low = u64::from(first.contig);
    let (start, end) = if first.reverse {
        let window_length = window_length as u64;
        (window_length - low - contig_span, window_length - low)
    } else {
        (low, low + contig_span)
    };
    let query_start = u64::from(first.query);
    Ok(Alignment {
        score: i32::try_from(score).map_err(|_| AlignmentError::LengthOverflow)?,
        strand: strand(first.reverse),
        query_interval: Interval::new(query_start, query_start + query_span)?,
        contig_interval: Interval::new(offset + start, offset + end)?,
        runs,
    })
}

/// Whether an unpaired end can reach a conserved block on a different diagonal: whether an
/// exact k-mer within `max_drift` bases beyond the end of the query pairs with one within
/// `max_drift` bases beyond the end of the contig at another distance from the end. Differences
/// next to the end, such as a substitution beside an insertion, do not hide the block. `codes`
/// is reused scratch of at most `max_drift + 1` entries.
pub fn gap_continuation(
    walk: &Walk,
    oriented: &[u8],
    table: &QueryTable,
    config: BandConfig,
    codes: &mut Vec<(u64, u32)>,
) -> Result<bool, AlignmentError> {
    let query = table.query();
    let k = usize::from(table.config().k);
    let reach = config.max_drift as usize + k;
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
        if query_bases < k || contig_bases < k {
            continue;
        }
        // The bases beyond the end, in query order; a k-mer's distance from the end counts the
        // bases between them.
        let query_from = if right {
            edge
        } else {
            edge + query.len() - query_bases
        };
        let query_side = query_window(
            query,
            (query_from % query.len()) as u64,
            query_bases as u64,
            table.circular(),
        )?;
        let contig_side = if right {
            &oriented[contig..contig + contig_bases]
        } else {
            &oriented[contig - contig_bases..contig]
        };
        let distance = |at: usize, bases: usize| if right { at } else { bases - k - at };
        // The end continues on its own diagonal, as at a circular query's whole length.
        if query_side[distance(0, query_bases)..][..k]
            == contig_side[distance(0, contig_bases)..][..k]
        {
            continue;
        }
        codes.clear();
        kmers(&query_side, table.config().k, |at, code, _| {
            codes.push((code, distance(at, query_bases) as u32))
        });
        codes.sort_unstable();
        let mut found = false;
        kmers(contig_side, table.config().k, |at, code, _| {
            let at = distance(at, contig_bases) as u32;
            // Entries of one code are ordered by distance, so another distance is at an end.
            let first = codes.partition_point(|entry| entry.0 < code);
            let last = codes.partition_point(|entry| entry.0 <= code);
            found |= first < last && (codes[first].1 != at || codes[last - 1].1 != at);
        });
        if found {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Whether `alignment` pairs the query positions of `walk` within its contig positions, on its
/// strand. `window_length` is the length of the window the walk was read from, which starts
/// `offset` bases into the contig; `query_length` is the length of the query.
pub fn contains(
    alignment: &Alignment,
    walk: &Walk,
    window_length: usize,
    offset: u64,
    query_length: u64,
) -> bool {
    if (alignment.strand == super::Strand::Reverse) != walk.reverse {
        return false;
    }
    let (low, high) = (u64::from(walk.contig), u64::from(walk.contig + walk.length));
    let (low, high) = if walk.reverse {
        (window_length as u64 - high, window_length as u64 - low)
    } else {
        (low, high)
    };
    let contig = alignment.contig_interval;
    if offset + low < contig.start || offset + high > contig.end {
        return false;
    }
    // The query positions of the alignment, which may run on over the origin, as pieces inside
    // the query.
    let (aligned_start, aligned_end) =
        (alignment.query_interval.start, alignment.query_interval.end);
    let segments = if aligned_end <= query_length {
        [(aligned_start, aligned_end), (0, 0)]
    } else if aligned_start >= query_length {
        [
            (aligned_start - query_length, aligned_end - query_length),
            (0, 0),
        ]
    } else {
        [
            (aligned_start, query_length),
            (0, aligned_end - query_length),
        ]
    };
    let (start, end) = (u64::from(walk.query), u64::from(walk.query + walk.length));
    [
        (start, end.min(query_length)),
        (0, end.saturating_sub(query_length)),
    ]
    .iter()
    .filter(|(from, to)| from < to)
    .all(|&(from, to)| {
        segments
            .iter()
            .any(|&(segment_start, segment_end)| segment_start <= from && to <= segment_end)
    })
}
