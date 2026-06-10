//! Hints that ask the kernel to read file ranges into the page cache before they are used. A hint
//! never changes what is read, and its errors are ignored.

use std::fs::File;

/// Sorts `(start, end)` spans and joins those that overlap, touch or lie within `gap` bytes.
pub fn merge_spans(spans: &mut Vec<(u64, u64)>, gap: u64) {
    spans.sort_unstable();
    let mut merged = Vec::<(u64, u64)>::with_capacity(spans.len());
    for &(start, end) in spans.iter() {
        if let Some((_, previous_end)) = merged.last_mut()
            && start <= previous_end.saturating_add(gap)
        {
            *previous_end = (*previous_end).max(end);
        } else {
            merged.push((start, end));
        }
    }
    *spans = merged;
}

/// Asks the kernel to read `length` bytes of the file from `offset`.
pub fn will_need(file: &File, offset: u64, length: u64) {
    #[cfg(target_os = "linux")]
    if let (Ok(offset), Ok(length)) = (libc::off_t::try_from(offset), libc::off_t::try_from(length))
    {
        use std::os::fd::AsRawFd;
        // SAFETY: the descriptor stays open for the call, and the call only hints.
        unsafe {
            libc::posix_fadvise(file.as_raw_fd(), offset, length, libc::POSIX_FADV_WILLNEED);
        }
    }
    #[cfg(not(target_os = "linux"))]
    let _ = (file, offset, length);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merges_overlapping_touching_and_near_spans() {
        let spans = [(30, 40), (0, 10), (10, 15), (5, 8), (100, 110), (42, 50)];
        let mut exact = spans.to_vec();
        merge_spans(&mut exact, 0);
        assert_eq!(exact, [(0, 15), (30, 40), (42, 50), (100, 110)]);
        let mut near = spans.to_vec();
        merge_spans(&mut near, 2);
        assert_eq!(near, [(0, 15), (30, 50), (100, 110)]);
        let mut none = Vec::new();
        merge_spans(&mut none, 0);
        assert!(none.is_empty());
    }

    #[test]
    fn hint_ignores_bad_ranges() {
        let file = tempfile::tempfile().unwrap();
        will_need(&file, 0, 4096);
        will_need(&file, u64::MAX, u64::MAX);
    }
}
