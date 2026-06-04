//! Reading one alignment window and finding its chains and dense walks.

use super::{AlignStageConfig, AlignStageError, AlignWork, PreparedQuery};
use crate::align::{
    Anchor, CHAINING_BYTES, Chain, Walk, dense_walks, forward_extent, window_chains,
};
use crate::index::PackedSequence;
use std::mem::size_of;
use std::ops::Range;

/// Buffers of window reads, kept across rounds and entries.
#[derive(Default)]
pub(super) struct WindowBuffers {
    anchors: Vec<Anchor>,
    bases: Vec<Vec<u8>>,
}

impl WindowBuffers {
    /// Keeps the bases buffer of a window that is done for the next read.
    pub fn recycle(&mut self, window: PreparedWindow) {
        self.bases.push(window.bases);
    }
}

/// One alignment window as read in one widening round, with its chains and dense walks.
pub(super) struct PreparedWindow {
    /// Contig position of the first base.
    pub start: u64,
    pub bases: Vec<u8>,
    /// Chains of anchors, best score first per strand.
    pub chains: Vec<Chain>,
    /// Copies in dense repeats and low-complexity runs, walked along their diagonals.
    pub walks: Vec<Walk>,
    /// The reverse complement of the window, which reverse walks are read from.
    pub reverse: Vec<u8>,
}

impl PreparedWindow {
    /// The contig positions of the window.
    pub fn range(&self) -> Range<u64> {
        self.start..self.start + self.bases.len() as u64
    }

    /// The window bases on a walk's strand.
    #[cfg(test)]
    pub fn oriented(&self, walk: &Walk) -> &[u8] {
        if walk.reverse {
            &self.reverse
        } else {
            &self.bases
        }
    }

    /// Bytes held by the bases, chains and walks.
    pub fn bytes(&self) -> usize {
        let chains: usize = self
            .chains
            .iter()
            .map(|chain| {
                size_of::<Chain>()
                    + chain.anchors.capacity() * size_of::<(u32, u32)>()
                    + chain.offsets.capacity() * size_of::<u32>()
            })
            .sum();
        self.bases.capacity()
            + self.reverse.capacity()
            + self.walks.capacity() * size_of::<Walk>()
            + chains
    }

    /// Whether a chain or walk lies within `reach` bases of the window's start and its end.
    pub fn near_edges(&self, k: u64, reach: u64) -> (bool, bool) {
        let length = self.bases.len() as u64;
        let walked = self.walks.iter().map(|walk| walk_extent(walk, length));
        let extents = self
            .chains
            .iter()
            .map(|chain| forward_extent(chain, length, k))
            .chain(walked);
        let (mut left, mut right) = (false, false);
        for (low, high) in extents {
            left |= low < reach;
            right |= high + reach > length;
        }
        (left, right)
    }
}

/// Forward window positions `start..end` of a walk.
pub(super) fn walk_extent(walk: &Walk, window_length: u64) -> (u64, u64) {
    let (low, high) = (u64::from(walk.contig), u64::from(walk.contig + walk.length));
    if walk.reverse {
        (window_length - high, window_length - low)
    } else {
        (low, high)
    }
}

/// Reads bases `range` of a contig and finds their chains and dense walks.
pub(super) fn prepare_window(
    packed: &PackedSequence,
    contig: u32,
    range: Range<u64>,
    query: &PreparedQuery,
    config: &AlignStageConfig,
    buffers: &mut WindowBuffers,
    work: &mut AlignWork,
) -> Result<PreparedWindow, AlignStageError> {
    let mut bases = buffers.bases.pop().unwrap_or_default();
    packed.read(contig as usize, range.clone(), &mut bases)?;
    work.windows += 1;
    work.bases_read += bases.len() as u64;
    let (chains, above) = window_chains(&query.table, &bases, &mut buffers.anchors)?;
    work.anchors += buffers.anchors.len() as u64;
    // K-mers too frequent for anchors, and those of low-complexity runs: their copies are
    // walked along their diagonals.
    let (mut walks, mut reverse) = match above {
        Some(above) => dense_walks(&query.table, &bases, above, config.chain, CHAINING_BYTES)?,
        None => Default::default(),
    };
    if let Some(low) = &query.low {
        let (more, copy) = dense_walks(low, &bases, 0, config.chain, CHAINING_BYTES)?;
        walks.extend(more);
        if reverse.is_empty() {
            reverse = copy;
        }
    }
    Ok(PreparedWindow {
        start: range.start,
        bases,
        chains,
        walks,
        reverse,
    })
}

/// Sorts spans and merges those that overlap or touch.
pub(super) fn merge_spans(spans: &mut Vec<Range<u64>>) {
    spans.sort_unstable_by_key(|span| (span.start, span.end));
    let mut merged: Vec<Range<u64>> = Vec::with_capacity(spans.len());
    for span in spans.drain(..) {
        match merged.last_mut() {
            Some(last) if span.start <= last.end => last.end = last.end.max(span.end),
            _ => merged.push(span),
        }
    }
    *spans = merged;
}

/// Grows `span` by `step` bases beyond the bases `read` on the sides asked, within a contig of
/// `contig_length` bases. Returns whether it grew.
pub(super) fn widen(
    span: &mut Range<u64>,
    read: &Range<u64>,
    (left, right): (bool, bool),
    step: u64,
    contig_length: u64,
) -> bool {
    let before = span.clone();
    if left {
        span.start = span.start.min(read.start.saturating_sub(step));
    }
    if right {
        span.end = span
            .end
            .max(read.end.saturating_add(step).min(contig_length));
    }
    *span != before
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::{BuildOptions, EntryInput, Parameters, Shard, build_shard};

    fn random_bases(length: usize, mut state: u64) -> Vec<u8> {
        (0..length)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state % 4) as usize]
            })
            .collect()
    }

    /// A shard of one entry with one contig.
    fn shard_of(directory: &std::path::Path, contig: &[u8]) -> Shard {
        let path = directory.join("entry.fa");
        std::fs::write(&path, [b">c\n".as_slice(), contig, b"\n"].concat()).unwrap();
        let options = BuildOptions {
            parameters: Parameters::default(),
            bucket_bits: Some(8),
            threads: 1,
            temp_dir: directory.to_path_buf(),
            entry_per_record: false,
        };
        let root = directory.join("shard");
        build_shard(&[EntryInput::from_path(path)], &root, &options).unwrap();
        Shard::open(&root).unwrap()
    }

    #[test]
    fn a_window_read_gives_its_chains_walks_and_edges() {
        let directory = tempfile::tempdir().unwrap();
        let contig = [random_bases(3000, 7), vec![b'T'; 200]].concat();
        let shard = shard_of(directory.path(), &contig);
        let packed = shard.packed(0).unwrap();
        let config = AlignStageConfig::default();
        let bases = [&contig[1000..2000], &[b'A'; 100][..]].concat();
        let query = PreparedQuery::new(&bases, false, &config).unwrap();
        let (mut buffers, mut work) = (WindowBuffers::default(), AlignWork::default());
        let mut read = |range: Range<u64>, buffers: &mut WindowBuffers| {
            prepare_window(&packed, 0, range, &query, &config, buffers, &mut work).unwrap()
        };
        let copy = read(900..2100, &mut buffers);
        assert_eq!(copy.range(), 900..2100);
        assert_eq!(copy.bases, contig[900..2100]);
        // The copy lies 100 bases inside each end of the window.
        assert!(copy.walks.is_empty());
        // Anchors start at every fourth window position: the last one at 1,084.
        assert_eq!(forward_extent(&copy.chains[0], 1200, 15), (100, 1099));
        assert_eq!(copy.near_edges(15, 200), (true, true));
        assert_eq!(copy.near_edges(15, 50), (false, false));
        buffers.recycle(copy);
        assert_eq!(buffers.bases.len(), 1);
        let run = read(3000..3200, &mut buffers);
        assert!(buffers.bases.is_empty());
        // The run of T pairs the run of A of the query on the reverse strand only, as walks.
        assert!(!run.walks.is_empty() && run.walks.iter().all(|walk| walk.reverse));
        assert_eq!(run.oriented(&run.walks[0]), vec![b'A'; 200]);
        assert!(run.bytes() >= 400);
        assert_eq!((work.windows, work.bases_read), (2, 1400));
    }

    #[test]
    fn spans_merge_when_they_overlap_or_touch() {
        let mut spans = vec![30..40, 0..10, 5..20, 20..25, 50..60];
        merge_spans(&mut spans);
        assert_eq!(spans, [0..25, 30..40, 50..60]);
    }

    #[test]
    fn spans_widen_once_per_round_within_the_contig() {
        let read = 100..200;
        let mut span = read.clone();
        assert!(widen(&mut span, &read, (true, false), 64, 1000));
        assert_eq!(span, 36..200);
        // A second request in the same round does not grow it again.
        assert!(!widen(&mut span, &read, (true, false), 64, 1000));
        assert!(widen(&mut span, &read, (true, true), 1000, 1000));
        assert_eq!(span, 0..1000);
        assert!(!widen(&mut span, &read, (true, true), 1000, 1000));
    }

    #[test]
    fn walks_on_the_reverse_strand_are_read_from_the_window_end() {
        let walk = Walk {
            reverse: true,
            contig: 10,
            query: 0,
            length: 30,
        };
        assert_eq!(walk_extent(&walk, 100), (60, 90));
        let forward = Walk {
            reverse: false,
            ..walk
        };
        assert_eq!(walk_extent(&forward, 100), (10, 40));
    }
}
