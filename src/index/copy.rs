//! Copies an index file and its sequence files between local directories and S3: the sequence
//! files first and the index file last, so a reader never sees an index without its files.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::manifest::sha256_range;
use super::seed_file::{HEADER_BYTES, SeedFile};
use super::sources::{IndexHead, Location};
use super::{AtPath, IndexError};
use crate::storage::{LocalFile, ReadAt};
use rayon::prelude::*;
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// Bytes of one ranged read while downloading.
const DOWNLOAD_PART: u64 = 8 << 20;
/// Ranged reads in flight while downloading.
const DOWNLOADS: usize = 16;

/// What one copy moved.
#[derive(Clone, Debug)]
pub struct CopySummary {
    /// The head of the copied index file.
    pub head: IndexHead,
    /// Bytes of all copied files.
    pub bytes: u64,
}

/// Copies the index file at `source` and its sequence files to `target`, the location of the new
/// index file; the sequence files keep their names relative to it. Downloads are checked against
/// the sizes and digests in the manifest, and an uploaded index file is read back.
pub fn copy_index(source: &Location, target: &Location) -> Result<CopySummary, IndexError> {
    let head = source.head()?;
    let files = &head.manifest.sequence_files;
    let outside = |path: &str| path.split('/').any(|part| part == ".." || part.is_empty());
    if let Some(file) = files.iter().find(|file| outside(&file.path)) {
        return Err(IndexError::Input(format!(
            "{} lies outside the directory of the index file; reshard with --copy-sequences \
             to copy it",
            file.path
        )));
    }
    let threads = rayon::ThreadPoolBuilder::new()
        .num_threads(DOWNLOADS)
        .build()?;
    for file in files {
        let check = |path: &Path| {
            let local = LocalFile::open(path)?;
            match local.size() == file.bytes && sha256_range(&local, 0..file.bytes)? == file.sha256
            {
                true => Ok(()),
                false => Err(IndexError::Checksum(format!(
                    "the copy of {} differs",
                    file.path
                ))),
            }
        };
        let from = source.open(&file.path, file.bytes)?;
        threads.install(|| copy_file(from.as_ref(), target, &file.path, check))?;
    }
    let from = source.open(&head.file, head.bytes)?;
    let check = |path: &Path| check_index(&Location::File(path.into()), &head);
    threads.install(|| copy_file(from.as_ref(), target, &target.file_name(), check))?;
    if let Location::S3 { .. } = target {
        // The reads are bound to the uploaded version, as its ETag is known.
        if target.head()?.digest != head.digest {
            return Err(IndexError::Checksum(
                "the uploaded index file differs".into(),
            ));
        }
    }
    let bytes = head.bytes + files.iter().map(|file| file.bytes).sum::<u64>();
    Ok(CopySummary { head, bytes })
}

/// Checks the index file at `location` against the head of its source: the digest of the
/// manifest, the checksums of the header and the manifest, and the digest of the seeds.
fn check_index(location: &Location, head: &IndexHead) -> Result<(), IndexError> {
    let copied = location.head()?;
    let file: Arc<dyn ReadAt> = location.open(&copied.file, copied.bytes)?;
    let seeds = HEADER_BYTES..SeedFile::open(Arc::clone(&file))?.manifest_range().start;
    let same = copied.digest == head.digest
        && sha256_range(file.as_ref(), seeds)? == head.manifest.seeds_sha256;
    match same {
        true => Ok(()),
        false => Err(IndexError::Checksum(
            "the copy of the index file differs".into(),
        )),
    }
}

/// Copies all bytes of `source` to the file `name` at `target`: an upload in parts to S3, or
/// ranged reads into a partial local file that `check` passes before it gets its name.
fn copy_file(
    source: &dyn ReadAt,
    target: &Location,
    name: &str,
    check: impl Fn(&Path) -> Result<(), IndexError>,
) -> Result<(), IndexError> {
    if target.upload(name, source)?.is_some() {
        return Ok(());
    }
    let path = target.local_path(name).unwrap_or_default();
    let parent = path.parent().unwrap_or(Path::new("."));
    std::fs::create_dir_all(parent).at(parent)?;
    let partial = PathBuf::from(format!("{}.partial-{}", path.display(), std::process::id()));
    let written = download(source, &partial).and_then(|()| check(&partial));
    match written {
        Ok(()) => std::fs::rename(&partial, &path).at(&path),
        Err(error) => {
            let _ = std::fs::remove_file(&partial);
            Err(error)
        }
    }
}

/// Writes all bytes of `source` into a new file at `path`, with ranged reads at once.
fn download(source: &dyn ReadAt, path: &Path) -> Result<(), IndexError> {
    let file = std::fs::File::create_new(path).at(path)?;
    file.set_len(source.size()).at(path)?;
    let parts = source.size().div_ceil(DOWNLOAD_PART);
    (0..parts).into_par_iter().try_for_each(|part| {
        let at = part * DOWNLOAD_PART;
        let mut buffer = vec![0; (source.size() - at).min(DOWNLOAD_PART) as usize];
        source.read_exact_at(at, &mut buffer)?;
        file.write_all_at(&buffer, at).at(path)
    })?;
    file.sync_all().at(path)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::Shard;
    use crate::storage::{S3Config, S3Limits};
    use crate::test_support::{random_bases, rng, runtime, serve_store, test_shard};

    /// A config for the store at `address`.
    fn config(address: &str) -> S3Config {
        let limits = S3Limits::default();
        let config = S3Config::new("local", Some(address), true, limits, runtime());
        config.unwrap().unsigned()
    }

    /// An index goes up to S3 and comes back the same. The ETags of the upload are kept in the
    /// cache directory, so a later search sends If-Match from its first read, and a download is
    /// checked against the digests.
    #[test]
    fn copies_both_ways() {
        let directory = tempfile::tempdir().unwrap();
        let mut rng = rng(12);
        let entries: Vec<Vec<Vec<u8>>> = (0..3)
            .map(|_| vec![random_bases(&mut rng, 200_000)])
            .collect();
        let local = directory.path().join("local");
        std::fs::create_dir(&local).unwrap();
        test_shard(&local, &entries);
        let store = directory.path().join("store");
        let (server, seen) = serve_store(&store);
        let cache = directory.path().join("cache");
        let up = config(&server.address);
        up.keep_etags(&cache).unwrap();
        let source = Location::File(local.join("shard.jam"));
        let target = Location::s3("s3://bucket/copy/shard.jam", &up).unwrap();
        copy_index(&source, &target).unwrap();
        for name in ["shard.jam", "shard.jamseq"] {
            let read = |path: PathBuf| std::fs::read(path).unwrap();
            assert!(
                read(local.join(name)) == read(store.join("copy").join(name)),
                "{name}"
            );
        }
        let puts = seen
            .lock()
            .unwrap()
            .iter()
            .filter(|r| r.method == "PUT")
            .count();
        assert_eq!(puts, 2);
        let log = std::fs::read_to_string(cache.join("etags.log")).unwrap();
        assert!(log.contains("s3://bucket/copy/shard.jam\t") && log.contains("shard.jamseq\t"));
        // A later run with the cache directory is bound to the uploaded version from the start.
        seen.lock().unwrap().clear();
        let later = config(&server.address);
        later.keep_etags(&cache).unwrap();
        let location = Location::s3("s3://bucket/copy/shard.jam", &later).unwrap();
        let shard = Shard::open_with(&location, &Default::default()).unwrap();
        shard.verify().unwrap();
        let requests = seen.lock().unwrap().clone();
        assert!(
            requests.iter().all(|r| r.if_match.is_some()),
            "{requests:?}"
        );
        // Back down, the same bytes; a changed byte in the store fails the digest check.
        let down = directory.path().join("down");
        let fresh = config(&server.address);
        let source = Location::s3("s3://bucket/copy/shard.jam", &fresh).unwrap();
        copy_index(&source, &Location::File(down.join("shard.jam"))).unwrap();
        assert!(
            std::fs::read(down.join("shard.jamseq")).unwrap()
                == std::fs::read(local.join("shard.jamseq")).unwrap()
        );
        let stored = store.join("copy/shard.jamseq");
        let mut bytes = std::fs::read(&stored).unwrap();
        bytes[5_000] ^= 1;
        std::fs::write(&stored, bytes).unwrap();
        let fresh = config(&server.address);
        let source = Location::s3("s3://bucket/copy/shard.jam", &fresh).unwrap();
        let other = directory.path().join("other");
        let error = copy_index(&source, &Location::File(other.join("shard.jam"))).unwrap_err();
        assert!(matches!(error, IndexError::Checksum(_)), "{error}");
        assert!(!other.join("shard.jamseq").exists() && !other.join("shard.jam").exists());
    }
}
