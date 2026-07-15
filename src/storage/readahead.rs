//! Hints that ask the kernel to read file ranges into the page cache before they are used. A hint
//! never changes what is read, and its errors are ignored.

use std::fs::File;

/// Asks the kernel to read `length` bytes of the file from `offset`.
pub(crate) fn will_need(file: &File, offset: u64, length: u64) {
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
    fn bad_hints_ignored() {
        let file = tempfile::tempfile().unwrap();
        will_need(&file, 0, 4096);
        will_need(&file, u64::MAX, u64::MAX);
    }
}
