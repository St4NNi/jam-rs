//! Fragments and the query windows they come from.

use super::DecideError;
use crate::align::{Alignment, Interval};
use serde::Serialize;

/// An alignment placed on the query and on one contig of an entry.
///
/// The alignment's query interval has the length of the query segments together. Over the end of
/// a circular query a fragment has two query segments, the first up to the query end and the
/// second from the query start.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Fragment {
    /// The index of the contig within its entry.
    pub contig_index: u32,
    /// The query bases the alignment spans, in alignment order.
    pub query_segments: Vec<Interval>,
    pub alignment: Alignment,
}

/// A stretch of the query that is aligned as one piece. On a circular query it may run past the
/// query end and continue at the query start.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct QueryWindow {
    start: u64,
    length: u64,
    query_length: u64,
    circular: bool,
}

impl QueryWindow {
    /// A window of `length` bases from `start`. A linear window must lie inside the query. A
    /// circular window may wrap once, so it is at most `query_length` long.
    pub fn new(
        start: u64,
        length: u64,
        query_length: u64,
        circular: bool,
    ) -> Result<Self, DecideError> {
        let end = start.checked_add(length);
        let fits = match circular {
            false => end.is_some_and(|end| end <= query_length),
            true => start < query_length && length <= query_length,
        };
        if query_length == 0 || !fits {
            return Err(DecideError::InvalidQueryWindow);
        }
        Ok(Self {
            start,
            length,
            query_length,
            circular,
        })
    }

    /// The bases of the window, copied out of `query`, which must have the length given to `new`.
    pub fn bases(&self, query: &[u8]) -> Result<Vec<u8>, DecideError> {
        if query.len() as u64 != self.query_length {
            return Err(DecideError::InvalidQueryWindow);
        }
        let (start, length) = (self.start as usize, self.length as usize);
        let first = length.min(query.len() - start);
        let mut bases = query[start..start + first].to_vec();
        bases.extend_from_slice(&query[..length - first]);
        Ok(bases)
    }

    /// The query segments of `local`, an interval of window positions (so `0..length` is the
    /// whole window). It gives two segments when the interval crosses the query end.
    pub fn segments(&self, local: Interval) -> Result<Vec<Interval>, DecideError> {
        if local.end > self.length {
            return Err(DecideError::InvalidQueryWindow);
        }
        let start = self.start + local.start;
        let end = self.start + local.end;
        if end <= self.query_length {
            return Ok(vec![Interval::new(start, end)?]);
        }
        if start >= self.query_length {
            let (start, end) = (start - self.query_length, end - self.query_length);
            return Ok(vec![Interval::new(start, end)?]);
        }
        Ok(vec![
            Interval::new(start, self.query_length)?,
            Interval::new(0, end - self.query_length)?,
        ])
    }
}

impl Fragment {
    /// A fragment from an alignment whose query interval is given in positions of `window`.
    pub fn in_window(
        contig_index: u32,
        window: &QueryWindow,
        alignment: Alignment,
    ) -> Result<Self, DecideError> {
        Ok(Self {
            contig_index,
            query_segments: window.segments(alignment.query_interval)?,
            alignment,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn interval(start: u64, end: u64) -> Interval {
        Interval::new(start, end).unwrap()
    }

    #[test]
    fn linear_windows_copy_and_map_directly() {
        let window = QueryWindow::new(2, 5, 10, false).unwrap();
        assert_eq!(window.bases(b"ABCDEFGHIJ").unwrap(), b"CDEFG");
        assert_eq!(window.segments(interval(1, 4)).unwrap(), [interval(3, 6)]);
        assert!(QueryWindow::new(6, 5, 10, false).is_err());
        assert!(window.bases(b"ABC").is_err());
    }

    #[test]
    fn circular_windows_wrap_past_the_end() {
        let window = QueryWindow::new(7, 6, 10, true).unwrap();
        assert_eq!(window.bases(b"ABCDEFGHIJ").unwrap(), b"HIJABC");
        assert_eq!(
            window.segments(interval(0, 6)).unwrap(),
            [interval(7, 10), interval(0, 3)]
        );
        assert_eq!(window.segments(interval(1, 3)).unwrap(), [interval(8, 10)]);
        assert_eq!(window.segments(interval(3, 5)).unwrap(), [interval(0, 2)]);
        assert_eq!(
            window.segments(interval(1, 4)).unwrap(),
            [interval(8, 10), interval(0, 1)]
        );
        assert!(window.segments(interval(0, 7)).is_err());
    }

    #[test]
    fn windows_that_do_not_fit_are_rejected() {
        assert!(QueryWindow::new(0, 0, 0, true).is_err());
        assert!(QueryWindow::new(10, 1, 10, true).is_err());
        assert!(QueryWindow::new(0, 11, 10, true).is_err());
        assert!(QueryWindow::new(u64::MAX, 2, 10, false).is_err());
        assert!(QueryWindow::new(3, 0, 10, true).is_ok());
    }
}
