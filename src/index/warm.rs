//! Warming a cache directory: download the pinned ranges of a shard with large ranged reads,
//! many at once, check their digests, and write the plan that later opens read.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::manifest::{EntryRecord, FileRecord, Manifest, PACKED_DIRECTORY, verify_source};
use super::plan::{PLAN_FILE, Pin, Plan};
use super::sources::Location;
use super::{AtPath, IndexError};
use crate::storage::{LocalFile, ReadAt};
use rayon::prelude::*;
use std::collections::HashSet;
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};

/// Bytes of one ranged read while downloading.
const DOWNLOAD_CHUNK: usize = 8 << 20;
/// Ranged reads in flight while downloading.
const DOWNLOAD_THREADS: usize = 64;

/// Downloads every pin of `plan` from the shard at `location` into `cache`, checks the digest of
/// each pin that covers a whole file, writes the plan, and removes files of older plans of the
/// shard. Pins already in place with the right digest are kept.
pub fn warm(
    location: &Location,
    manifest: &Manifest,
    plan: &Plan,
    cache: &Path,
) -> Result<(), IndexError> {
    let threads = rayon::ThreadPoolBuilder::new()
        .num_threads(DOWNLOAD_THREADS)
        .build()?;
    threads.install(|| {
        (plan.pins.par_iter()).try_for_each(|pin| {
            let record = file_record(manifest, &pin.file)?;
            let path = plan.pin_path(cache, pin);
            if is_complete(&path, pin, record)? {
                return Ok(());
            }
            download(location, record, pin, &path)
        })
    })?;
    plan.write(cache)?;
    remove_unplanned(plan, cache)
}

/// The manifest record of a shard file named relative to the shard.
fn file_record<'a>(manifest: &'a Manifest, file: &str) -> Result<&'a FileRecord, IndexError> {
    let packed = |entry: &&EntryRecord| {
        file.strip_prefix(PACKED_DIRECTORY)
            .and_then(|rest| rest.strip_prefix('/'))
            == entry.packed.as_ref().map(|packed| packed.path.as_str())
    };
    [&manifest.seeds]
        .into_iter()
        .find(|record| record.path == file)
        .or_else(|| {
            manifest
                .entries
                .iter()
                .find(packed)
                .and_then(|entry| entry.packed.as_ref())
        })
        .ok_or_else(|| IndexError::Input(format!("the shard has no file {file}")))
}

/// Whether the pin file at `path` is there with its size and, for a whole file, its digest.
fn is_complete(path: &Path, pin: &Pin, record: &FileRecord) -> Result<bool, IndexError> {
    let Ok(file) = LocalFile::open(path) else {
        return Ok(false);
    };
    if file.size() != pin.bytes() {
        return Ok(false);
    }
    Ok(pin.bytes() != record.bytes || verify_source(&file, record).is_ok())
}

/// Downloads the range of `pin` into a partial file, checks it, and moves it to `path`.
fn download(
    location: &Location,
    record: &FileRecord,
    pin: &Pin,
    path: &Path,
) -> Result<(), IndexError> {
    let source = location.open(&pin.file, record.bytes)?;
    let parent = path.parent().unwrap_or(Path::new("."));
    std::fs::create_dir_all(parent).at(parent)?;
    let partial = PathBuf::from(format!("{}.partial", path.display()));
    let file = std::fs::File::create(&partial).at(&partial)?;
    file.set_len(pin.bytes()).at(&partial)?;
    let chunks = pin.bytes().div_ceil(DOWNLOAD_CHUNK as u64);
    (0..chunks).into_par_iter().try_for_each(|chunk| {
        let at = chunk * DOWNLOAD_CHUNK as u64;
        let mut buffer = vec![0; (pin.bytes() - at).min(DOWNLOAD_CHUNK as u64) as usize];
        source.read_exact_at(pin.start + at, &mut buffer)?;
        file.write_all_at(&buffer, at).at(&partial)
    })?;
    file.sync_all().at(&partial)?;
    if pin.bytes() == record.bytes {
        verify_source(&LocalFile::open(&partial)?, record)?;
    }
    std::fs::rename(&partial, path).at(path)
}

/// Removes the files in the shard's cache directory that the plan does not pin.
fn remove_unplanned(plan: &Plan, cache: &Path) -> Result<(), IndexError> {
    let planned: HashSet<PathBuf> = (plan.pins.iter())
        .map(|pin| plan.pin_path(cache, pin))
        .collect();
    let directory = cache.join(&plan.shard);
    for directory in [directory.clone(), directory.join(PACKED_DIRECTORY)] {
        let Ok(listing) = std::fs::read_dir(&directory) else {
            continue;
        };
        for item in listing {
            let path = item.at(&directory)?.path();
            if path.is_file()
                && !planned.contains(&path)
                && path.file_name() != Some(PLAN_FILE.as_ref())
            {
                std::fs::remove_file(&path).at(&path)?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::{BuildOptions, Part, Place, Shard, ShardOptions, build_shard};
    use crate::storage::{S3Config, S3Limits};
    use crate::test_support::{
        build_options, random_bases, rng, runtime, serve_directory, test_shard, write_entries,
    };

    #[test]
    fn warm_then_read() {
        let directory = tempfile::tempdir().unwrap();
        let mut rng = rng(3);
        let entries: Vec<Vec<Vec<u8>>> = (0..3)
            .map(|_| vec![random_bases(&mut rng, 20_000)])
            .collect();
        let local = test_shard(directory.path(), &entries);
        let (server, asked) = serve_directory(directory.path());
        let endpoint = Some(server.address.as_str());
        let config = S3Config::new("local", endpoint, true, S3Limits::default(), runtime());
        let location = Location::s3("s3://bucket/shard", &config.unwrap().unsigned()).unwrap();
        let (manifest, digest) = location.manifest().unwrap();
        assert_eq!(digest, local.digest());
        let cache = directory.path().join("cache");
        let parts = [Part::Sketch, Part::Packed(Some(vec!["entry2".into()]))];
        let plan = Plan::choose(&manifest, &digest, 1 << 30, 1 << 30, &parts).unwrap();
        assert_eq!(plan.pins[0].place, Place::Ram);
        let stale = cache.join(&digest).join("positions");
        std::fs::create_dir_all(stale.parent().unwrap()).unwrap();
        std::fs::write(&stale, b"old").unwrap();
        warm(&location, &manifest, &plan, &cache).unwrap();
        assert!(!stale.exists());
        assert_eq!(Plan::read(&cache, &digest).unwrap(), Some(plan.clone()));
        // A second warm finds every pin in place and downloads nothing.
        let before = asked.lock().unwrap().len();
        warm(&location, &manifest, &plan, &cache).unwrap();
        assert_eq!(asked.lock().unwrap().len(), before);
        // A damaged pin is downloaded again.
        let pinned = plan.pin_path(&cache, &plan.pins[1]);
        let mut bytes = std::fs::read(&pinned).unwrap();
        bytes[100] ^= 1;
        std::fs::write(&pinned, bytes).unwrap();
        warm(&location, &manifest, &plan, &cache).unwrap();
        assert!(asked.lock().unwrap().len() > before);
        let options = ShardOptions {
            pins: Some(cache),
            ..ShardOptions::default()
        };
        let shard = Shard::open_with(&location, &options).unwrap();
        shard.verify().unwrap();
        let mut bases = Vec::new();
        shard
            .packed(2)
            .unwrap()
            .read(0, 0..20_000, &mut bases)
            .unwrap();
        assert_eq!(bases, entries[2][0]);
    }

    #[test]
    fn sketch_only_warmed() {
        let directory = tempfile::tempdir().unwrap();
        let mut rng = rng(4);
        let entries: Vec<Vec<Vec<u8>>> = (0..2)
            .map(|_| vec![random_bases(&mut rng, 20_000)])
            .collect();
        let inputs = write_entries(directory.path(), &entries);
        let options = BuildOptions {
            sketch_only: true,
            ..build_options(directory.path(), 2)
        };
        build_shard(&inputs, &directory.path().join("shard"), &options).unwrap();
        let (server, asked) = serve_directory(directory.path());
        let endpoint = Some(server.address.as_str());
        let config = S3Config::new("local", endpoint, true, S3Limits::default(), runtime());
        let location = Location::s3("s3://bucket/shard", &config.unwrap().unsigned()).unwrap();
        let (manifest, digest) = location.manifest().unwrap();
        let cache = directory.path().join("cache");
        // The default plan pins the sketch part alone; there are no locators or packed files.
        let plan = Plan::choose(&manifest, &digest, 1 << 30, 1 << 30, &[]).unwrap();
        let parts: Vec<&str> = plan.pins.iter().map(|pin| pin.part.as_str()).collect();
        assert_eq!(parts, ["sketch"]);
        let packed = Plan::choose(&manifest, &digest, 1 << 30, 1 << 30, &[Part::Packed(None)]);
        assert!(packed.unwrap().pins.is_empty());
        warm(&location, &manifest, &plan, &cache).unwrap();
        let options = ShardOptions {
            pins: Some(cache),
            ..ShardOptions::default()
        };
        let before = asked.lock().unwrap().len();
        let shard = Shard::open_with(&location, &options).unwrap();
        let ids = (0..manifest.color_sets as u32).map(crate::index::ColorSetId);
        assert!(shard.seed_file().read_color_sets(ids).is_ok());
        // Opening reads the manifest; the seeds file is read from the pin alone.
        let asked = asked.lock().unwrap();
        assert!(asked[before..].iter().all(|key| !key.ends_with("seeds")));
    }
}
