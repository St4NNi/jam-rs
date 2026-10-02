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
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Condvar, Mutex, PoisonError};

/// The name of the log file in the cache directory.
const LOG_FILE: &str = "etags.log";

/// The error of a read that found its object replaced.
#[derive(Debug, thiserror::Error)]
#[error("the index changed in S3 while it was read; retry")]
pub struct ChangedError;

pub(crate) fn changed_error() -> io::Error {
    io::Error::other(ChangedError)
}

/// Whether `error` says that an object was replaced.
pub fn is_changed(error: &io::Error) -> bool {
    error
        .get_ref()
        .is_some_and(|inner| inner.is::<ChangedError>())
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

/// The check of one object, once per process: a read of its first byte that the known ETag must
/// match.
#[derive(Default)]
pub(crate) struct Latch {
    /// `None` while the check runs; then its outcome, with the reason of a failure.
    outcome: Mutex<Option<Result<(), (bool, String)>>>,
    done: Condvar,
}

impl Latch {
    /// Waits for the end of the check.
    pub(crate) fn wait(&self) -> io::Result<()> {
        let mut outcome = lock(&self.outcome);
        while outcome.is_none() {
            outcome = (self.done.wait(outcome)).unwrap_or_else(PoisonError::into_inner);
        }
        match outcome.as_ref() {
            Some(Err((true, _))) => Err(changed_error()),
            Some(Err((false, reason))) => Err(io::Error::other(reason.clone())),
            _ => Ok(()),
        }
    }
}

/// What a read finds when it asks for the check of its object.
pub(crate) enum Start {
    /// The check passed.
    Passed,
    /// Another read started the check; wait for its latch.
    Wait(Arc<Latch>),
    /// This read starts the check and must end the latch.
    Run(Arc<Latch>),
}

enum Check {
    Running(Arc<Latch>),
    Passed,
}

#[derive(Default)]
struct Inner {
    checks: HashMap<String, Check>,
    /// Each ETag with its version number.
    etags: HashMap<String, (String, u64)>,
    log: Option<(PathBuf, File)>,
}

impl EtagBook {
    /// The ETag of `location`, if it is known.
    pub(crate) fn get(&self, location: &str) -> Option<String> {
        lock(&self.inner)
            .etags
            .get(location)
            .map(|(etag, _)| etag.clone())
    }

    /// The version number of the ETag of `location`, if it is known.
    pub(crate) fn version(&self, location: &str) -> Option<u64> {
        lock(&self.inner)
            .etags
            .get(location)
            .map(|(_, version)| *version)
    }

    /// Records `seen` as the ETag of `location` if none is known. Fails with the known ETag
    /// when it differs.
    pub(crate) fn settle(&self, location: &str, seen: &str) -> Result<(), String> {
        let mut inner = lock(&self.inner);
        if let Some((known, _)) = inner.etags.get(location) {
            return if known == seen {
                Ok(())
            } else {
                Err(known.clone())
            };
        }
        inner
            .etags
            .insert(location.to_string(), (seen.to_string(), version_of(seen)));
        inner.append(location, seen);
        Ok(())
    }

    /// Forgets the ETag of `location`, if it is still `old`, and counts the change.
    pub(crate) fn changed(&self, location: &str, old: Option<&str>) {
        let mut inner = lock(&self.inner);
        if inner.etags.get(location).map(|(etag, _)| etag.as_str()) == old {
            inner.etags.remove(location);
            inner.checks.remove(location);
            inner.append(location, "");
            self.changes.fetch_add(1, Ordering::Relaxed);
        }
    }

    /// Asks for the check of `location`.
    pub(crate) fn begin_check(&self, location: &str) -> Start {
        let mut inner = lock(&self.inner);
        match inner.checks.get(location) {
            Some(Check::Passed) => Start::Passed,
            Some(Check::Running(latch)) => Start::Wait(Arc::clone(latch)),
            None => {
                let latch = Arc::<Latch>::default();
                inner
                    .checks
                    .insert(location.to_string(), Check::Running(Arc::clone(&latch)));
                Start::Run(latch)
            }
        }
    }

    /// Marks the check of `location` as passed when its ETag is known and no check runs.
    pub(crate) fn pass_check(&self, location: &str) {
        let mut inner = lock(&self.inner);
        if inner.etags.contains_key(location) && !inner.checks.contains_key(location) {
            inner.checks.insert(location.to_string(), Check::Passed);
        }
    }

    /// Ends the check of `location` that `latch` belongs to. A passed check is not run again; a
    /// failed one is run again by the next read.
    pub(crate) fn end_check(&self, location: &str, latch: &Latch, result: &io::Result<()>) {
        let mut inner = lock(&self.inner);
        let running = matches!(inner.checks.get(location), Some(Check::Running(l)) if std::ptr::eq(&**l, latch));
        if running {
            match result {
                Ok(()) => inner.checks.insert(location.to_string(), Check::Passed),
                Err(_) => inner.checks.remove(location),
            };
        }
        let failure = result
            .as_ref()
            .err()
            .map(|e| (is_changed(e), e.to_string()));
        *lock(&latch.outcome) = Some(failure.map_or(Ok(()), Err));
        latch.done.notify_all();
    }

    /// Forgets the ETag of `location` without counting a change.
    pub(crate) fn forget(&self, location: &str) {
        let mut inner = lock(&self.inner);
        if inner.etags.remove(location).is_some() {
            inner.checks.remove(location);
            inner.append(location, "");
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
            let version = version_of(&etag);
            inner.etags.entry(location).or_insert((etag, version));
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
