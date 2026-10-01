//! A source that reads one byte range of another source as a whole file.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::{ReadAt, checked_end};
use std::io;
use std::ops::Range;
use std::sync::Arc;

/// The bytes `range` of `inner`, read from offset 0, for example one assembly in a sequence file.
pub struct SliceSource {
    inner: Arc<dyn ReadAt>,
    range: Range<u64>,
}

impl SliceSource {
    /// The bytes `range` of `inner`; it fails when the range does not lie within `inner`.
    pub fn new(inner: Arc<dyn ReadAt>, range: Range<u64>) -> io::Result<Self> {
        let length = (range.end.checked_sub(range.start)).ok_or_else(super::past_the_end)?;
        checked_end(range.start, length as usize, inner.size())?;
        Ok(Self { inner, range })
    }
}

impl ReadAt for SliceSource {
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> io::Result<()> {
        checked_end(offset, buffer.len(), self.size())?;
        self.inner.read_exact_at(self.range.start + offset, buffer)
    }

    fn size(&self) -> u64 {
        self.range.end - self.range.start
    }

    fn version(&self) -> io::Result<u64> {
        self.inner.version()
    }

    fn will_need(&self, offset: u64, length: u64) {
        let length = length.min(self.size().saturating_sub(offset));
        if length > 0 {
            self.inner.will_need(self.range.start + offset, length);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slice_reads_range() {
        let data: Vec<u8> = (0..=255).cycle().take(10_000).collect();
        let source = SliceSource::new(Arc::new(data.clone()), 4_096..6_000).unwrap();
        assert_eq!(source.size(), 1_904);
        let mut buffer = vec![0; 100];
        source.read_exact_at(1_804, &mut buffer).unwrap();
        assert_eq!(buffer, data[5_900..6_000]);
        assert!(source.read_exact_at(1_805, &mut buffer).is_err());
        assert!(SliceSource::new(Arc::new(data.clone()), 9_000..10_001).is_err());
        let (start, end) = (5, 4);
        assert!(SliceSource::new(Arc::new(data), start..end).is_err());
        source.will_need(1_900, 100);
    }
}
