//! The bucket table of the seeds file: where each bucket of a part lies, so a lookup reads exactly
//! one byte range. Buckets never cross a page unless they are larger than one.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::IndexError;
use super::bits::{u16_at, u32_at};
use std::ops::Range;

/// Bytes of one page row: buckets that start before the page (u32), data end in the page (u32).
pub(crate) const PAGE_ROW_BYTES: usize = 8;

/// One page of a part.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct PageRow {
    /// Buckets that start in earlier pages, so the first bucket that may start in this page.
    starts_before: u32,
    /// Bytes of bucket data in this page; the rest of the page is padding.
    data_end: u32,
}

/// Where every bucket of one part lies: per page a row, per bucket its start within its page.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct PartTable {
    page_bits: u32,
    pages: Vec<PageRow>,
    /// One start per bucket and one more for the end of the part.
    starts: Vec<u16>,
}

impl PartTable {
    /// Decodes the table of a part of `buckets` buckets with `pages` page rows from `bytes`: the
    /// page rows, then the starts. Checks that the rows ascend and fit `part_bytes`.
    pub(crate) fn decode(
        bytes: &[u8],
        page_bits: u32,
        buckets: usize,
        pages: usize,
        part_bytes: u64,
    ) -> Result<Self, IndexError> {
        let damaged = || IndexError::format("seeds bucket table");
        let rows = (0..pages)
            .map(|page| {
                let at = page * PAGE_ROW_BYTES;
                Some(PageRow {
                    starts_before: u32_at(bytes, at)?,
                    data_end: u32_at(bytes, at + 4)?,
                })
            })
            .collect::<Option<Vec<_>>>()
            .ok_or_else(damaged)?;
        let at = pages * PAGE_ROW_BYTES;
        let starts = (0..=buckets)
            .map(|bucket| u16_at(bytes, at + 2 * bucket))
            .collect::<Option<Vec<_>>>()
            .ok_or_else(damaged)?;
        let page_bytes = 1u64 << page_bits;
        let ordered = rows.first().is_some_and(|row| row.starts_before == 0)
            && rows.is_sorted_by_key(|row| row.starts_before)
            && rows.iter().all(|row| u64::from(row.data_end) <= page_bytes)
            && rows
                .last()
                .is_some_and(|row| row.starts_before as usize <= buckets)
            && starts.iter().all(|start| u64::from(*start) < page_bytes);
        let table = Self {
            page_bits,
            pages: rows,
            starts,
        };
        if !ordered || table.start(buckets) != part_bytes {
            return Err(damaged());
        }
        Ok(table)
    }

    /// Bytes of the encoded table of a part.
    pub(crate) fn encoded_bytes(buckets: usize, pages: usize) -> usize {
        pages * PAGE_ROW_BYTES + 2 * (buckets + 1)
    }

    /// Page rows of the table.
    pub(crate) fn page_count(&self) -> usize {
        self.pages.len()
    }

    /// The bytes of `bucket` within the part, or `None` for a bucket past the last.
    pub(crate) fn range(&self, bucket: usize) -> Option<Range<u64>> {
        let next = bucket
            .checked_add(1)
            .filter(|next| *next < self.starts.len())?;
        let (page, start) = (self.page_of(bucket), self.start(bucket));
        let next_page = self.page_of(next);
        let end = if next_page == page || self.starts[next] != 0 {
            self.start(next)
        } else {
            // The next bucket starts a page, so padding may lie before it.
            let last = next_page - 1;
            ((last as u64) << self.page_bits) + u64::from(self.pages[last].data_end)
        };
        (start <= end).then_some(start..end)
    }

    /// Writes the table: the page rows, then the starts.
    pub(crate) fn encode(&self, out: &mut Vec<u8>) {
        for row in &self.pages {
            out.extend_from_slice(&row.starts_before.to_le_bytes());
            out.extend_from_slice(&row.data_end.to_le_bytes());
        }
        for start in &self.starts {
            out.extend_from_slice(&start.to_le_bytes());
        }
    }

    /// The page that `bucket` starts in: the last page whose earlier pages hold no later bucket.
    fn page_of(&self, bucket: usize) -> usize {
        let bucket = bucket as u32;
        self.pages
            .partition_point(|row| row.starts_before <= bucket)
            .saturating_sub(1)
    }

    fn start(&self, bucket: usize) -> u64 {
        ((self.page_of(bucket) as u64) << self.page_bits) + u64::from(self.starts[bucket])
    }
}

/// Places the buckets of one part as they are written, in bucket order.
pub(crate) struct PartLayout {
    table: PartTable,
    /// Bytes of the part so far, padding included.
    written: u64,
}

impl PartLayout {
    /// An empty part with pages of 2^`page_bits` bytes (at most 64 KiB).
    pub(crate) fn new(page_bits: u32) -> Self {
        Self {
            table: PartTable {
                page_bits,
                pages: Vec::new(),
                starts: Vec::new(),
            },
            written: 0,
        }
    }

    /// Places the next bucket of `length` bytes and returns the padding to write before it.
    pub(crate) fn place(&mut self, length: u64) -> u64 {
        let page_bytes = 1u64 << self.table.page_bits;
        let in_page = self.written % page_bytes;
        let crosses = in_page + length > page_bytes;
        let padding = if length > 0 && in_page > 0 && crosses {
            let page = (self.written >> self.table.page_bits) as usize;
            self.add_pages_through(page);
            self.table.pages[page].data_end = in_page as u32;
            page_bytes - in_page
        } else {
            0
        };
        self.written += padding;
        self.start_bucket();
        self.written += length;
        padding
    }

    /// Bytes of the part so far.
    pub(crate) fn written(&self) -> u64 {
        self.written
    }

    /// Ends the part: adds the start that closes the last bucket.
    pub(crate) fn finish(mut self) -> PartTable {
        self.start_bucket();
        let page_bytes = 1u64 << self.table.page_bits;
        if let Some(row) = self.table.pages.last_mut() {
            row.data_end = (self.written % page_bytes) as u32;
        }
        self.table
    }

    /// Records a bucket that starts at the current end, adding rows up to its page.
    fn start_bucket(&mut self) {
        let page_bits = self.table.page_bits;
        let page = (self.written >> page_bits) as usize;
        self.add_pages_through(page);
        let start = self.written - ((page as u64) << page_bits);
        self.table.starts.push(start as u16);
    }

    /// Adds full page rows up to `page`; the buckets placed so far start before them.
    fn add_pages_through(&mut self, page: usize) {
        let bucket = self.table.starts.len() as u32;
        while self.table.pages.len() <= page {
            self.table.pages.push(PageRow {
                starts_before: bucket,
                data_end: 1 << self.table.page_bits,
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::rng;
    use rand::RngExt;

    /// Lays out buckets of these lengths and checks every range against the written bytes.
    fn check_layout(lengths: &[u64], page_bits: u32) {
        let mut layout = PartLayout::new(page_bits);
        let mut placed = Vec::new();
        for &length in lengths {
            let padding = layout.place(length);
            let start = layout.written() - length;
            assert!(padding == 0 || start.is_multiple_of(1 << page_bits));
            placed.push(start..start + length);
        }
        let bytes = layout.written();
        let table = layout.finish();
        let mut encoded = Vec::new();
        table.encode(&mut encoded);
        let pages = table.page_count();
        assert_eq!(
            encoded.len(),
            PartTable::encoded_bytes(lengths.len(), pages)
        );
        let decoded = PartTable::decode(&encoded, page_bits, lengths.len(), pages, bytes).unwrap();
        assert_eq!(decoded, table);
        let page_bytes = 1u64 << page_bits;
        for (bucket, expected) in placed.iter().enumerate() {
            let range = decoded.range(bucket).unwrap();
            assert_eq!(range, *expected, "bucket {bucket}");
            let fits = expected.end - expected.start <= page_bytes;
            if fits && !range.is_empty() {
                assert_eq!(range.start / page_bytes, (range.end - 1) / page_bytes);
            }
        }
        assert_eq!(decoded.range(lengths.len()), None);
        assert!(PartTable::decode(&encoded, page_bits, lengths.len(), pages, bytes + 1).is_err());
    }

    #[test]
    fn buckets_fit_pages() {
        let mut draws = rng(11);
        for page_bits in [6, 8, 16] {
            let page = 1u64 << page_bits;
            let lengths: Vec<u64> = (0..2_000)
                .map(|_| match draws.random_range(0..10) {
                    0 => 0,
                    1 => draws.random_range(page..4 * page),
                    2 => page,
                    _ => draws.random_range(1..page / 2),
                })
                .collect();
            check_layout(&lengths, page_bits);
        }
        check_layout(&[], 6);
        check_layout(&[0, 0, 0], 6);
        check_layout(&[64, 0, 64, 64, 0], 6);
        check_layout(&[200, 0, 10, 54, 10], 6);
    }

    #[test]
    fn damaged_table_refused() {
        let mut layout = PartLayout::new(6);
        for length in [10, 60, 100, 5] {
            layout.place(length);
        }
        let bytes = layout.written();
        let table = layout.finish();
        let mut encoded = Vec::new();
        table.encode(&mut encoded);
        let pages = table.page_count();
        assert!(PartTable::decode(&encoded[..encoded.len() - 1], 6, 4, pages, bytes).is_err());
        let mut bad = encoded.clone();
        bad[0] = 1;
        assert!(PartTable::decode(&bad, 6, 4, pages, bytes).is_err());
    }
}
