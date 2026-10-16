//! A source with one byte range held locally (pinned), and the rest read from another source.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::counts::count_pinned;
use super::{ReadAt, checked_end};
use std::io;
use std::ops::Range;
use std::sync::Arc;

/// A source whose bytes `range` come from a local copy and all others from `rest`, for example
/// S3 through the block cache.
pub struct PinnedSource {
    pinned: Arc<dyn ReadAt>,
    range: Range<u64>,
    rest: Arc<dyn ReadAt>,
}

impl PinnedSource {
    /// Reads `range` from `pinned`, which holds exactly those bytes, and the rest from `rest`.
    pub fn new(pinned: Arc<dyn ReadAt>, range: Range<u64>, rest: Arc<dyn ReadAt>) -> Self {
        Self {
            pinned,
            range,
            rest,
        }
    }

    /// The part of `offset..end` that lies in the pinned range; empty at the end of the part
    /// before it when they do not meet.
    fn inside(&self, offset: u64, end: u64) -> Range<u64> {
        let start = offset.max(self.range.start).min(end);
        start..end.min(self.range.end).max(start)
    }
}

impl ReadAt for PinnedSource {
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> io::Result<()> {
        let end = checked_end(offset, buffer.len(), self.rest.size())?;
        let inside = self.inside(offset, end);
        let (head, tail) = buffer.split_at_mut((inside.start - offset) as usize);
        let (middle, tail) = tail.split_at_mut((inside.end - inside.start) as usize);
        if !head.is_empty() {
            self.rest.read_exact_at(offset, head)?;
        }
        if !middle.is_empty() {
            count_pinned();
            self.pinned
                .read_exact_at(inside.start - self.range.start, middle)?;
        }
        if !tail.is_empty() {
            self.rest.read_exact_at(inside.end, tail)?;
        }
        Ok(())
    }

    fn size(&self) -> u64 {
        self.rest.size()
    }

    fn version(&self) -> io::Result<u64> {
        self.rest.version()
    }

    fn will_need(&self, offset: u64, length: u64) {
        let end = offset.saturating_add(length);
        let inside = self.inside(offset, end);
        if inside.start < inside.end {
            let at = inside.start - self.range.start;
            self.pinned.will_need(at, inside.end - inside.start);
        }
        if offset < inside.start {
            self.rest.will_need(offset, inside.start - offset);
        }
        if inside.end < end {
            self.rest.will_need(inside.end, end - inside.end);
        }
    }
}

/// A source whose pages of `page` bytes at `starts` come from a local copy that holds them one
/// after another, and all other bytes from `rest`. A shard pins the header page of every packed
/// sequence of a sequence file this way.
pub struct PinnedPages {
    pinned: Arc<dyn ReadAt>,
    /// Where the pages start, ascending, each at least `page` bytes after the one before.
    starts: Vec<u64>,
    page: u64,
    rest: Arc<dyn ReadAt>,
}

impl PinnedPages {
    /// Reads the pages at `starts` (ascending, not overlapping) from `pinned`, the rest from
    /// `rest`.
    pub fn new(
        pinned: Arc<dyn ReadAt>,
        starts: Vec<u64>,
        page: u64,
        rest: Arc<dyn ReadAt>,
    ) -> Self {
        Self {
            pinned,
            starts,
            page,
            rest,
        }
    }

    /// The part of `at..end` that starts at `at` and lies in one page, or before the next one:
    /// its end, and the place of its page in the pinned copy when it is pinned.
    fn piece(&self, at: u64, end: u64) -> (u64, Option<u64>) {
        let next = self.starts.partition_point(|&start| start <= at);
        if let Some(index) = next.checked_sub(1) {
            let start = self.starts[index];
            if at < start + self.page {
                let place = index as u64 * self.page + at - start;
                return ((start + self.page).min(end), Some(place));
            }
        }
        let stop = self.starts.get(next).map_or(end, |&start| start.min(end));
        (stop, None)
    }
}

impl ReadAt for PinnedPages {
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> io::Result<()> {
        let end = checked_end(offset, buffer.len(), self.rest.size())?;
        let mut at = offset;
        while at < end {
            let (stop, place) = self.piece(at, end);
            let part = &mut buffer[(at - offset) as usize..(stop - offset) as usize];
            match place {
                Some(place) => {
                    count_pinned();
                    self.pinned.read_exact_at(place, part)?;
                }
                None => self.rest.read_exact_at(at, part)?,
            }
            at = stop;
        }
        Ok(())
    }

    fn size(&self) -> u64 {
        self.rest.size()
    }

    fn version(&self) -> io::Result<u64> {
        self.rest.version()
    }

    fn will_need(&self, offset: u64, length: u64) {
        let end = offset.saturating_add(length).min(self.rest.size());
        let mut at = offset;
        while at < end {
            let (stop, place) = self.piece(at, end);
            match place {
                Some(place) => self.pinned.will_need(place, stop - at),
                None => self.rest.will_need(at, stop - at),
            }
            at = stop;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ranges_split_exactly() {
        let data: Vec<u8> = (0..=255).cycle().take(10_000).collect();
        let pinned = Arc::new(data[3_000..7_000].to_vec());
        // The rest holds zeros where the pin lies, so bytes from the wrong side show.
        let mut rest = data.clone();
        rest[3_000..7_000].fill(0);
        let source = PinnedSource::new(pinned, 3_000..7_000, Arc::new(rest));
        assert_eq!(source.size(), 10_000);
        for (offset, length) in [
            (0, 100),
            (2_990, 20),
            (3_000, 4_000),
            (6_990, 20),
            (0, 10_000),
            (8_000, 50),
        ] {
            let mut buffer = vec![0; length];
            source.read_exact_at(offset as u64, &mut buffer).unwrap();
            assert_eq!(buffer, data[offset..offset + length], "{offset} {length}");
        }
        assert!(source.read_exact_at(9_999, &mut [0; 2]).is_err());
        source.will_need(0, 10_000);
        source.will_need(9_000, u64::MAX);
    }

    #[test]
    fn pages_split_exactly() {
        let data: Vec<u8> = (0..=255).cycle().take(10_000).collect();
        let starts = vec![1_000, 4_000, 4_100];
        let page = 100;
        let pinned: Vec<u8> = (starts.iter())
            .flat_map(|&start| data[start as usize..(start + page) as usize].to_vec())
            .collect();
        let mut rest = data.clone();
        for &start in &starts {
            rest[start as usize..(start + page) as usize].fill(0);
        }
        let source = PinnedPages::new(Arc::new(pinned), starts, page, Arc::new(rest));
        assert_eq!(source.size(), 10_000);
        for (offset, length) in [
            (0, 100),
            (990, 20),
            (1_000, 100),
            (1_050, 3_100),
            (4_000, 200),
            (4_150, 10),
            (0, 10_000),
            (9_999, 1),
        ] {
            let mut buffer = vec![0; length];
            source.read_exact_at(offset as u64, &mut buffer).unwrap();
            assert_eq!(buffer, data[offset..offset + length], "{offset} {length}");
        }
        assert!(source.read_exact_at(9_999, &mut [0; 2]).is_err());
        source.will_need(0, 10_000);
        source.will_need(9_000, u64::MAX);
    }
}
