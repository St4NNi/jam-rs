//! Alignments placed as fragments, and the fragments one entry collects within its limit.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::{AlignStageError, Limit, PreparedQuery};
use crate::align::{EditRun, Interval};
use crate::decide::Fragment;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::mem::size_of;

/// Whether the query has bases before and after an alignment with this query interval. A
/// circular query has room on both sides until an alignment spans all of it.
pub(super) fn query_room(interval: Interval, query: &PreparedQuery) -> (bool, bool) {
    if query.circular() {
        let room = interval.len() < query.length();
        (room, room)
    } else {
        (interval.start > 0, interval.end < query.length())
    }
}

/// The fragments of one entry and the bytes they hold, within a limit. Past the limit no more
/// fragments are kept, but their bytes are still counted, so the error names what was needed.
/// As a last resort it keeps the fragments of the highest scores within the limit instead.
pub(super) struct Collected {
    list: Vec<Fragment>,
    /// Bytes of the fragments added since the last reset, kept or not.
    bytes: usize,
    /// The most that `bytes` reached.
    peak: usize,
    limit: usize,
    /// With `Some`, the fragments of the highest scores within the limit, in place of `list`.
    best: Option<Best>,
}

/// A count of fragments and their bytes, to go back to when a round is aligned again.
#[derive(Clone, Copy)]
pub(super) struct Mark {
    length: usize,
    bytes: usize,
    /// Fragments added before the mark.
    added: u64,
}

/// The fragments of the highest scores within a limit. Past the limit the fragments of the
/// lowest score are left out first, and of those the latest. The places they leave empty are
/// freed once they hold an eighth of the limit, and the fragments held take the rest.
#[derive(Default)]
struct Best {
    /// The fragments held in the order they came, and an empty place for each one left out.
    places: Vec<Option<Fragment>>,
    /// The number of each place: how many fragments came before it.
    numbers: Vec<u64>,
    /// Score and number of the fragments held, lowest score first, the latest first of equals.
    /// Entries of places dropped by a reset are skipped.
    lowest: BinaryHeap<Reverse<(i32, Reverse<u64>)>>,
    /// Fragments that came so far.
    added: u64,
    /// Bytes of the fragments held, with their numbers and entries by score.
    held: usize,
    /// Empty places.
    empty: usize,
}

/// Bytes of a held fragment beside the fragment: its number and its entry by score.
const HELD_BYTES: usize = size_of::<u64>() + size_of::<Reverse<(i32, Reverse<u64>)>>();
/// Bytes of an empty place and its number.
const EMPTY_BYTES: usize = size_of::<Option<Fragment>>() + size_of::<u64>();

impl Best {
    fn push(&mut self, fragment: Fragment, limit: usize) {
        let number = self.added;
        self.added += 1;
        let size = fragment_bytes(&fragment) + HELD_BYTES;
        let key = (fragment.alignment.score, Reverse(number));
        let room = limit - limit / 8;
        if size > room {
            return;
        }
        while self.held + size > room {
            match self.lowest_held() {
                Some(lowest) if lowest < key => self.leave_out(lowest.1.0, limit),
                // The new fragment has the lowest score.
                _ => return,
            }
        }
        self.places.push(Some(fragment));
        self.numbers.push(number);
        self.lowest.push(Reverse(key));
        self.held += size;
    }

    /// The key of the held fragment of the lowest score.
    fn lowest_held(&mut self) -> Option<(i32, Reverse<u64>)> {
        while let Some(&Reverse(key)) = self.lowest.peek() {
            let place = self.numbers.binary_search(&key.1.0).ok();
            if place.is_some_and(|at| self.places[at].is_some()) {
                return Some(key);
            }
            self.lowest.pop();
        }
        None
    }

    /// Leaves out the held fragment of the lowest score, numbered `number`.
    fn leave_out(&mut self, number: u64, limit: usize) {
        self.lowest.pop();
        let at = (self.numbers.binary_search(&number)).expect("the fragment is held");
        let fragment = self.places[at].take().expect("the fragment is held");
        self.held -= fragment_bytes(&fragment) + HELD_BYTES;
        self.empty += 1;
        if self.empty * EMPTY_BYTES > limit / 8 {
            self.free_places();
        }
    }

    /// Frees the empty places.
    fn free_places(&mut self) {
        let mut kept = self.places.iter().map(Option::is_some);
        self.numbers.retain(|_| kept.next() == Some(true));
        self.places.retain(Option::is_some);
        self.empty = 0;
        self.rebuild();
    }

    /// The scores of the fragments held, anew.
    fn rebuild(&mut self) {
        let places = self.places.iter().zip(&self.numbers);
        let held = places.filter_map(|(place, &number)| {
            let fragment = place.as_ref()?;
            Some(Reverse((fragment.alignment.score, Reverse(number))))
        });
        self.lowest = held.collect();
    }

    /// Drops the fragments that came after the first `added`.
    fn reset(&mut self, added: u64) {
        let from = self.numbers.partition_point(|&number| number < added);
        for place in self.places.drain(from..) {
            match place {
                Some(fragment) => self.held -= fragment_bytes(&fragment) + HELD_BYTES,
                None => self.empty -= 1,
            }
        }
        self.numbers.truncate(from);
        self.rebuild();
    }
}

/// The bytes that a fragment holds.
pub(super) fn fragment_bytes(fragment: &Fragment) -> usize {
    size_of::<Fragment>()
        + fragment.query_segments.capacity() * size_of::<Interval>()
        + fragment.alignment.runs.capacity() * size_of::<EditRun>()
}

impl Collected {
    pub(super) fn new(limit: usize) -> Self {
        Self {
            list: Vec::new(),
            bytes: 0,
            peak: 0,
            limit,
            best: None,
        }
    }

    /// A collection that keeps the fragments of the highest scores within `limit`, and never
    /// fails for the limit.
    pub(super) fn keeping_best(limit: usize) -> Self {
        Self {
            best: Some(Best::default()),
            ..Self::new(limit)
        }
    }

    /// Adds a fragment; once the fragments pass the limit, it is only counted.
    pub(super) fn push(&mut self, fragment: Fragment) {
        self.bytes += fragment_bytes(&fragment);
        self.peak = self.peak.max(self.bytes);
        if let Some(best) = &mut self.best {
            best.push(fragment, self.limit);
            return;
        }
        // Past the limit the list is incomplete, which only a reset before that point mends.
        if self.bytes <= self.limit {
            self.list.push(fragment);
        }
    }

    pub(super) fn mark(&self) -> Mark {
        Mark {
            length: self.list.len(),
            bytes: self.bytes,
            added: self.best.as_ref().map_or(0, |best| best.added),
        }
    }

    /// Drops the fragments added after `mark`.
    pub(super) fn reset(&mut self, mark: Mark) {
        self.list.truncate(mark.length);
        self.bytes = mark.bytes;
        if let Some(best) = &mut self.best {
            best.reset(mark.added);
        }
    }

    /// The bytes of the fragments added since the last reset, kept or not.
    pub(super) fn needed(&self) -> usize {
        self.bytes
    }

    /// The fragments, or the error of the limit when they passed it, with the most bytes they
    /// needed at once. Fragments kept by score are never an error.
    pub(super) fn into_list(self) -> Result<Vec<Fragment>, AlignStageError> {
        if let Some(best) = self.best {
            return Ok(best.places.into_iter().flatten().collect());
        }
        if self.peak > self.limit {
            return Err(AlignStageError::TooLarge {
                what: "the fragments of one entry",
                needed: self.peak as u64,
                limit: self.limit as u64,
                kind: Limit::FragmentBytes,
            });
        }
        Ok(self.list)
    }
}

#[cfg(test)]
mod tests {
    use super::super::AlignStageConfig;
    use super::*;
    use crate::align::{Alignment, Strand, parse_cigar};

    fn alignment(query: (u64, u64)) -> Alignment {
        let length = query.1 - query.0;
        Alignment {
            score: 1,
            strand: Strand::Forward,
            query_interval: Interval::new(query.0, query.1).unwrap(),
            contig_interval: Interval::new(0, length).unwrap(),
            runs: parse_cigar(&format!("{length}=")).unwrap(),
        }
    }

    fn interval(start: u64, end: u64) -> Interval {
        Interval::new(start, end).unwrap()
    }

    #[test]
    fn query_room_sides() {
        let bases = b"ACGTTGCAAC".repeat(10);
        let config = AlignStageConfig::default();
        let circular = PreparedQuery::new(&bases, true, &config).unwrap();
        let linear = PreparedQuery::new(&bases, false, &config).unwrap();
        assert_eq!(query_room(interval(0, 100), &linear), (false, false));
        assert_eq!(query_room(interval(5, 50), &linear), (true, true));
        assert_eq!(query_room(interval(0, 50), &circular), (true, true));
        assert_eq!(query_room(interval(50, 150), &circular), (false, false));
    }

    #[test]
    fn fragment_limit_rounds() {
        let bases = b"ACGTTGCAAC".repeat(10);
        let query = PreparedQuery::new(&bases, false, &AlignStageConfig::default()).unwrap();
        let fragment = Fragment::place(0, alignment((0, 50)), query.length(), false).unwrap();
        // A clone holds no spare capacity, so all pushed fragments have this size.
        let one = fragment_bytes(&fragment.clone());
        let mut collected = Collected::new(2 * one);
        collected.push(fragment.clone());
        let mark = collected.mark();
        collected.push(fragment.clone());
        collected.reset(mark);
        collected.push(fragment.clone());
        assert_eq!(collected.into_list().unwrap().len(), 2);
        // Past the limit the error names the most bytes needed at once, four fragments, even
        // when a reset dropped them later, as the entry held them at once.
        let mut collected = Collected::new(2 * one);
        let mark = collected.mark();
        for _ in 0..4 {
            collected.push(fragment.clone());
        }
        collected.reset(mark);
        collected.push(fragment.clone());
        let error = collected.into_list().unwrap_err();
        let needed = (4 * one) as u64;
        assert!(matches!(error, AlignStageError::TooLarge { needed: n, .. } if n == needed));
    }

    /// Kept by score, the fragments of the highest scores stay within the limit, in the order
    /// they came, and a reset drops the later ones.
    #[test]
    fn best_fragments_kept() {
        let bases = b"ACGTTGCAAC".repeat(10);
        let query = PreparedQuery::new(&bases, false, &AlignStageConfig::default()).unwrap();
        let scored = |score: i32| {
            let place = Fragment::place(0, alignment((0, 50)), query.length(), false);
            let mut fragment = place.unwrap();
            fragment.alignment.score = score;
            // A clone holds no spare capacity, so all fragments have the same size.
            fragment.clone()
        };
        let one = fragment_bytes(&scored(1)) + HELD_BYTES;
        let scores = |list: Vec<Fragment>| -> Vec<i32> {
            list.iter()
                .map(|fragment| fragment.alignment.score)
                .collect()
        };
        // The limit that holds `count` fragments, beside the eighth for empty places.
        let room_for = |count: usize| count * one * 8 / 7 + 8;
        // Room for three: of 5, 1, 7, 3, 9 and 2 the three best stay, in their order.
        let mut collected = Collected::keeping_best(room_for(3));
        for score in [5, 1, 7, 3, 9, 2] {
            collected.push(scored(score));
        }
        assert_eq!(collected.needed(), 6 * (one - HELD_BYTES));
        assert_eq!(scores(collected.into_list().unwrap()), [5, 7, 9]);
        // Of equal scores the first stay; a reset drops what came after the mark.
        let mut collected = Collected::keeping_best(room_for(2));
        collected.push(scored(4));
        let mark = collected.mark();
        for score in [4, 4, 8] {
            collected.push(scored(score));
        }
        collected.reset(mark);
        collected.push(scored(6));
        assert_eq!(scores(collected.into_list().unwrap()), [4, 6]);
        // A fragment larger than the limit is never held.
        let mut collected = Collected::keeping_best(one / 2);
        collected.push(scored(9));
        assert!(collected.into_list().unwrap().is_empty());
    }
}
