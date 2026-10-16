//! Warming a cache directory: download the pinned ranges of a shard with large ranged reads,
//! many at once, check their digests, and write the plan that later opens read.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::manifest::{EntryRecord, FileRecord, sha256_range};
use super::packed::{BASES_AT, header_page_intact};
use super::plan::{PLAN_FILE, Pin, Plan, etag_file, header_starts};
use super::seed_file::{HEADER_BYTES, verify_part};
use super::sources::{IndexHead, Location};
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

/// Downloads every pin of `plan` from the index at `location` (head `head`) into `cache`, checks
/// each, and writes the plan. The old plan and its pins that `plan` lacks go first, so the disk
/// budget holds while new pins arrive; pins already in place and intact are kept.
pub fn warm(
    location: &Location,
    head: &IndexHead,
    plan: &Plan,
    cache: &Path,
) -> Result<(), IndexError> {
    location.keep_etags(cache)?;
    remove_obsolete(plan, cache)?;
    let threads = rayon::ThreadPoolBuilder::new()
        .num_threads(DOWNLOAD_THREADS)
        .build()?;
    let seeds = FileRecord {
        path: head.file.clone(),
        bytes: head.bytes,
        sha256: head.digest.clone(),
    };
    let files = &head.manifest.sequence_files;
    threads.install(|| {
        (plan.pins.par_iter()).try_for_each(|pin| {
            let file = (files.iter().position(|file| file.path == pin.file))
                .map(|place| (&head.manifest.entries[..], place));
            let record = match file {
                Some((_, place)) => &files[place],
                None if pin.file == head.file => &seeds,
                None => return Err(IndexError::Input(format!("no file {}", pin.file))),
            };
            let pages =
                (file.zip(pin.pages)).map(|((_, place), _)| header_starts(&head.manifest, place));
            let pinned = Pinned {
                pin,
                record,
                file,
                pages,
            };
            let path = plan.pin_path(cache, pin);
            if is_complete(location, &path, &pinned) {
                return Ok(());
            }
            download(location, &pinned, &path)
        })
    })?;
    plan.write(cache)
}

/// A pin with the record of its file; for a sequence file also the entries and the place of the
/// file in the manifest, and for a pin of headers where its pages start.
struct Pinned<'a> {
    pin: &'a Pin,
    record: &'a FileRecord,
    file: Option<(&'a [EntryRecord], usize)>,
    pages: Option<Vec<u64>>,
}

/// Whether the pin file at `path` is there and intact, as `verify_pin` checks it, and its
/// object on S3 still has the ETag that the pin was downloaded under.
fn is_complete(location: &Location, path: &Path, pinned: &Pinned<'_>) -> bool {
    let intact = LocalFile::open(path).is_ok_and(|f| verify_pin(location, &f, pinned).is_ok());
    let etag = std::fs::read_to_string(etag_file(path)).ok();
    let (file, bytes) = (&pinned.pin.file, pinned.record.bytes);
    intact && (location.verify_etag(file, bytes, etag.as_deref())).is_ok()
}

/// Checks a pin file: its size; for a sequence file the digest of the whole file, else of each
/// packed sequence the range holds; for the seeds file the CRC-32 of every section that the
/// range holds, with the header from the pin when it starts the file, else from `location`.
fn verify_pin(
    location: &Location,
    file: &LocalFile,
    pinned: &Pinned<'_>,
) -> Result<(), IndexError> {
    let Pinned { pin, record, .. } = *pinned;
    if file.size() != pin.bytes() {
        return Err(IndexError::format(format!(
            "the pin of {} is cut",
            pin.file
        )));
    }
    let damaged = || IndexError::Checksum(format!("the pin of {} is damaged", pin.file));
    if let Some(starts) = &pinned.pages {
        let mut page = vec![0u8; BASES_AT as usize];
        for index in 0..starts.len() as u64 {
            file.read_exact_at(index * BASES_AT, &mut page)?;
            if !header_page_intact(&page) {
                return Err(damaged());
            }
        }
        return Ok(());
    }
    if let Some((entries, place)) = pinned.file {
        if pin.bytes() == record.bytes {
            let whole = sha256_range(file, 0..pin.bytes())? == record.sha256;
            return whole.then_some(()).ok_or_else(damaged);
        }
        let inside = (entries.iter().filter_map(|entry| entry.sequence.as_ref()))
            .filter(|s| s.file == place && pin.start <= s.offset && s.range().end <= pin.end);
        for sequence in inside {
            let range = sequence.offset - pin.start..sequence.range().end - pin.start;
            if sha256_range(file, range)? != sequence.sha256 {
                return Err(damaged());
            }
        }
        return Ok(());
    }
    let mut header = [0u8; HEADER_BYTES as usize];
    match pin.start {
        0 => file.read_exact_at(0, &mut header)?,
        _ => (location.open(&pin.file, record.bytes)?).read_exact_at(0, &mut header)?,
    }
    verify_part(&header, file, pin.start..pin.end)
}

/// Downloads the range of `pin` into a partial file, checks it, and moves it to `path`.
fn download(location: &Location, pinned: &Pinned<'_>, path: &Path) -> Result<(), IndexError> {
    let Pinned { pin, record, .. } = *pinned;
    let etag_path = etag_file(path);
    let _ = std::fs::remove_file(&etag_path);
    let etag = location.refresh_etag(&pin.file, record.bytes)?;
    let source = location.open(&pin.file, record.bytes)?;
    let parent = path.parent().unwrap_or(Path::new("."));
    std::fs::create_dir_all(parent).at(parent)?;
    let partial = PathBuf::from(format!("{}.partial", path.display()));
    let file = std::fs::File::create(&partial).at(&partial)?;
    file.set_len(pin.bytes()).at(&partial)?;
    if let Some(starts) = &pinned.pages {
        (starts.par_iter().enumerate()).try_for_each(|(index, &start)| {
            let mut page = vec![0; BASES_AT as usize];
            source.read_exact_at(start, &mut page)?;
            file.write_all_at(&page, index as u64 * BASES_AT)
                .at(&partial)
        })?;
    }
    let chunks = match pinned.pages {
        Some(_) => 0,
        None => pin.bytes().div_ceil(DOWNLOAD_CHUNK as u64),
    };
    (0..chunks).into_par_iter().try_for_each(|chunk| {
        let at = chunk * DOWNLOAD_CHUNK as u64;
        let mut buffer = vec![0; (pin.bytes() - at).min(DOWNLOAD_CHUNK as u64) as usize];
        source.read_exact_at(pin.start + at, &mut buffer)?;
        file.write_all_at(&buffer, at).at(&partial)
    })?;
    file.sync_all().at(&partial)?;
    verify_pin(location, &LocalFile::open(&partial)?, pinned)?;
    std::fs::rename(&partial, path).at(path)?;
    match etag {
        Some(etag) => std::fs::write(&etag_path, etag).at(&etag_path),
        None => Ok(()),
    }
}

/// Removes the plan of the shard in `cache`, then every pin of it that `plan` does not hold, by
/// its exact path. A plan that cannot be read names no pins to remove.
fn remove_obsolete(plan: &Plan, cache: &Path) -> Result<(), IndexError> {
    let Ok(Some(previous)) = Plan::read(cache, &plan.shard) else {
        return Ok(());
    };
    let path = cache.join(&plan.shard).join(PLAN_FILE);
    std::fs::remove_file(&path).at(&path)?;
    let planned: HashSet<PathBuf> = (plan.pins.iter())
        .map(|pin| plan.pin_path(cache, pin))
        .collect();
    for pin in &previous.pins {
        let path = previous.pin_path(cache, pin);
        if planned.contains(&path) {
            continue;
        }
        let _ = std::fs::remove_file(etag_file(&path));
        match std::fs::remove_file(&path) {
            Err(error) if error.kind() != std::io::ErrorKind::NotFound => {
                return Err(error).at(&path);
            }
            _ => {}
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
        let location = Location::s3("s3://bucket/shard.jam", &config.unwrap().unsigned()).unwrap();
        let head = location.head().unwrap();
        assert_eq!(head.digest, local.digest());
        let cache = directory.path().join("cache");
        let parts = [
            Part::Sketch,
            Part::Packed(Some(vec!["entry2".into()])),
            Part::Headers,
        ];
        let plan = Plan::choose(&head, 1 << 30, 1 << 30, &parts).unwrap();
        assert_eq!(plan.pins[0].place, Place::Ram);
        assert_eq!(plan.pins[2].place, Place::Ssd);
        warm(&location, &head, &plan, &cache).unwrap();
        assert_eq!(
            Plan::read(&cache, &head.digest).unwrap(),
            Some(plan.clone())
        );
        // A second warm finds every pin in place and downloads nothing; each pin only costs its
        // check read.
        let before = asked.lock().unwrap().len();
        warm(&location, &head, &plan, &cache).unwrap();
        assert_eq!(asked.lock().unwrap().len(), before + plan.pins.len());
        // A damaged pin is downloaded again, also a pin of part of the seeds file and of headers.
        for pin in &plan.pins {
            let pinned = plan.pin_path(&cache, pin);
            let mut bytes = std::fs::read(&pinned).unwrap();
            let at = bytes.len() - 100;
            bytes[at] ^= 1;
            std::fs::write(&pinned, &bytes).unwrap();
            let before = asked.lock().unwrap().len();
            warm(&location, &head, &plan, &cache).unwrap();
            assert!(asked.lock().unwrap().len() > before, "{}", pin.file);
            assert_ne!(std::fs::read(&pinned).unwrap(), bytes);
        }
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

    /// Bytes of all files below `directory`, except the ETag files.
    fn disk_bytes(directory: &Path) -> u64 {
        let listing = std::fs::read_dir(directory).unwrap();
        (listing.map(|item| item.unwrap().path()))
            .filter(|path| {
                !path.to_string_lossy().contains(".etag") && !path.ends_with("etags.log")
            })
            .map(|path| match path.is_dir() {
                true => disk_bytes(&path),
                false => std::fs::metadata(&path).unwrap().len(),
            })
            .sum()
    }

    /// Rewarming with fewer parts removes the pins of the parts left out, so the cache holds
    /// only the pins of the new plan.
    #[test]
    fn rewarm_removes_pins() {
        let directory = tempfile::tempdir().unwrap();
        let mut rng = rng(5);
        let entries: Vec<Vec<Vec<u8>>> = (0..2)
            .map(|_| vec![random_bases(&mut rng, 20_000)])
            .collect();
        test_shard(directory.path(), &entries);
        let (server, _) = serve_directory(directory.path());
        let endpoint = Some(server.address.as_str());
        let config = S3Config::new("local", endpoint, true, S3Limits::default(), runtime());
        let location = Location::s3("s3://bucket/shard.jam", &config.unwrap().unsigned()).unwrap();
        let head = location.head().unwrap();
        let cache = directory.path().join("cache");
        let all = [Part::Sketch, Part::Locators, Part::Packed(None)];
        let wide = Plan::choose(&head, 1 << 30, 1 << 30, &all).unwrap();
        warm(&location, &head, &wide, &cache).unwrap();
        let narrow = Plan::choose(&head, 1 << 30, 1 << 30, &[Part::Sketch]).unwrap();
        warm(&location, &head, &narrow, &cache).unwrap();
        for pin in &wide.pins[1..] {
            assert!(!wide.pin_path(&cache, pin).exists(), "{}", pin.file);
        }
        assert!(narrow.pin_path(&cache, &narrow.pins[0]).exists());
        let plan_bytes = std::fs::metadata(cache.join(&head.digest).join(PLAN_FILE)).unwrap();
        let pinned: u64 = narrow.pins.iter().map(Pin::bytes).sum();
        assert_eq!(disk_bytes(&cache), pinned + plan_bytes.len());
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
        build_shard(&inputs, &directory.path().join("shard.jam"), &options).unwrap();
        let (server, asked) = serve_directory(directory.path());
        let endpoint = Some(server.address.as_str());
        let config = S3Config::new("local", endpoint, true, S3Limits::default(), runtime());
        let location = Location::s3("s3://bucket/shard.jam", &config.unwrap().unsigned()).unwrap();
        let head = location.head().unwrap();
        let cache = directory.path().join("cache");
        // The default plan pins the sketch part alone; there are no locators or packed files.
        let plan = Plan::choose(&head, 1 << 30, 1 << 30, &[]).unwrap();
        let parts: Vec<&str> = plan.pins.iter().map(|pin| pin.part.as_str()).collect();
        assert_eq!(parts, ["sketch"]);
        let packed = Plan::choose(&head, 1 << 30, 1 << 30, &[Part::Packed(None)]);
        assert!(packed.unwrap().pins.is_empty());
        warm(&location, &head, &plan, &cache).unwrap();
        let options = ShardOptions {
            pins: Some(cache),
            ..ShardOptions::default()
        };
        let before = asked.lock().unwrap().len();
        let shard = Shard::open_with(&location, &options).unwrap();
        let ids = (0..head.manifest.color_sets as u32).map(crate::index::ColorSetId);
        assert!(shard.seed_file().read_color_sets(ids).is_ok());
        // Opening reads the header and the manifest, which also check the index file; the
        // seeds are read from the pin.
        let asked = asked.lock().unwrap();
        let seeds = asked[before..].iter().filter(|key| *key == "shard.jam");
        assert_eq!(seeds.count(), 2);
    }
}
