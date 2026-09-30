//! The ETags that S3 objects had when this process first read them, kept in memory and, when a
//! cache directory is given, in a log file there for later runs.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::lock;
use sha2::{Digest, Sha256};
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicU64, Ordering};

/// The name of the log file in the cache directory.
const LOG_FILE: &str = "etags.log";

/// The error of a read that found its object replaced.
pub(crate) fn changed_error() -> io::Error {
    io::Error::other("the index changed in S3 while it was read; retry")
}

/// A number that names an ETag, for block keys.
pub(crate) fn version_of(etag: &str) -> u64 {
    let digest = Sha256::digest(etag.as_bytes());
    u64::from_le_bytes(digest[..8].try_into().unwrap_or_default())
}

/// The known ETag of each object by its `s3://` location.
#[derive(Default)]
pub(crate) struct EtagBook {
    inner: Mutex<Inner>,
    /// How often an object was found replaced.
    changes: AtomicU64,
}

#[derive(Default)]
struct Inner {
    etags: HashMap<String, String>,
    log: Option<(PathBuf, File)>,
}

impl EtagBook {
    /// The ETag of `location`, if it is known.
    pub(crate) fn get(&self, location: &str) -> Option<String> {
        lock(&self.inner).etags.get(location).cloned()
    }

    /// Records `seen` as the ETag of `location` if none is known. Fails with the known ETag
    /// when it differs.
    pub(crate) fn settle(&self, location: &str, seen: &str) -> Result<(), String> {
        let mut inner = lock(&self.inner);
        if let Some(known) = inner.etags.get(location) {
            return if known == seen {
                Ok(())
            } else {
                Err(known.clone())
            };
        }
        inner.etags.insert(location.to_string(), seen.to_string());
        inner.append(location, seen);
        Ok(())
    }

    /// Forgets the ETag of `location`, if it is still `old`, and counts the change.
    pub(crate) fn changed(&self, location: &str, old: Option<&str>) {
        let mut inner = lock(&self.inner);
        if inner.etags.get(location).map(String::as_str) == old {
            inner.etags.remove(location);
            inner.append(location, "");
            self.changes.fetch_add(1, Ordering::Relaxed);
        }
    }

    /// How often an object was found replaced.
    pub(crate) fn changes(&self) -> u64 {
        self.changes.load(Ordering::Relaxed)
    }

    /// Loads the ETags that earlier runs left in the log of `directory`, and appends new ones
    /// there from now on. ETags known already win.
    pub(crate) fn keep(&self, directory: &Path) -> io::Result<()> {
        let path = directory.join(LOG_FILE);
        let mut inner = lock(&self.inner);
        if inner.log.as_ref().is_some_and(|(kept, _)| *kept == path) {
            return Ok(());
        }
        std::fs::create_dir_all(directory)?;
        let file = OpenOptions::new().create(true).append(true).open(&path)?;
        let text = std::fs::read_to_string(&path).unwrap_or_default();
        let mut loaded = HashMap::new();
        for line in text.lines() {
            match line.split_once('\t') {
                Some((location, "")) => loaded.remove(location),
                Some((location, etag)) => loaded.insert(location.to_string(), etag.to_string()),
                None => None,
            };
        }
        for (location, etag) in loaded {
            inner.etags.entry(location).or_insert(etag);
        }
        inner.log = Some((path, file));
        Ok(())
    }
}

impl Inner {
    /// Appends a line to the log; an empty `etag` forgets the location. Best effort, as a lost
    /// line only means the ETag is learned again.
    fn append(&mut self, location: &str, etag: &str) {
        if let Some((_, file)) = &mut self.log {
            let _ = file.write_all(format!("{location}\t{etag}\n").as_bytes());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_survives_runs() {
        let directory = tempfile::tempdir().unwrap();
        let book = EtagBook::default();
        book.keep(directory.path()).unwrap();
        book.settle("s3://b/one", "\"1\"").unwrap();
        book.settle("s3://b/two", "\"2\"").unwrap();
        assert_eq!(book.settle("s3://b/one", "\"9\""), Err("\"1\"".into()));
        book.changed("s3://b/two", Some("\"2\""));
        assert_eq!(book.changes(), 1);
        let next = EtagBook::default();
        next.keep(directory.path()).unwrap();
        assert_eq!(next.get("s3://b/one").as_deref(), Some("\"1\""));
        assert_eq!(next.get("s3://b/two"), None);
        // A change of an ETag that is no longer the known one is not counted.
        book.changed("s3://b/one", Some("\"7\""));
        assert_eq!(book.changes(), 1);
    }
}
