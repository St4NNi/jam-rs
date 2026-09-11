//! A source with one byte range held locally (pinned), and the rest read from another source.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

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
}
