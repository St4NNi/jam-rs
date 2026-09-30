//! Tests of shards on S3 whose objects are learned, checked and replaced.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::*;
use crate::storage::{BlockCache, CacheConfig, S3Config, S3Limits};
use crate::test_support::{Seen, random_bases, rng, runtime, serve_seen, test_shard};
use std::path::Path;
use std::sync::{Arc, Mutex};

/// Bases per entry, about a dozen blocks of a packed file.
const LENGTH: usize = 200_000;

/// Three entries of one contig each, from the generator `seed`.
fn entries_of(seed: u64) -> Vec<Vec<Vec<u8>>> {
    let mut rng = rng(seed);
    (0..3)
        .map(|_| vec![random_bases(&mut rng, LENGTH)])
        .collect()
}

/// A store over `root` and a config for it.
fn store(root: &Path) -> (crate::test_support::Server, Arc<Mutex<Vec<Seen>>>, S3Config) {
    let (server, seen) = serve_seen(root);
    let endpoint = Some(server.address.as_str());
    let config = S3Config::new("local", endpoint, true, S3Limits::default(), runtime());
    (server, seen, config.unwrap().unsigned())
}

/// Reads `range` of the first contig of entry `entry`, and checks the bases.
fn read_range(
    shard: &Shard,
    entries: &[Vec<Vec<u8>>],
    entry: usize,
    range: std::ops::Range<u64>,
) -> Result<(), IndexError> {
    let mut bases = Vec::new();
    shard
        .packed(entry as u32)?
        .read(0, range.clone(), &mut bases)?;
    assert_eq!(
        bases,
        entries[entry][0][range.start as usize..range.end as usize]
    );
    Ok(())
}

/// Reads all bases of all entries, or fails with the first error.
fn read_all(shard: &Shard, entries: &[Vec<Vec<u8>>]) -> Result<(), IndexError> {
    (0..entries.len()).try_for_each(|entry| read_range(shard, entries, entry, 0..LENGTH as u64))
}

/// Replaces the files of the shard in `root` by those of another shard built from `entries`;
/// with `packed` only the packed sequences.
fn replace(root: &Path, entries: &[Vec<Vec<u8>>], packed: bool) {
    let other = tempfile::tempdir().unwrap();
    test_shard(other.path(), entries);
    let mut stack = vec![other.path().join("shard")];
    while let Some(directory) = stack.pop() {
        for file in std::fs::read_dir(directory).unwrap() {
            let path = file.unwrap().path();
            if path.is_dir() {
                stack.push(path);
                continue;
            }
            let relative = path.strip_prefix(other.path()).unwrap();
            if packed && !relative.starts_with("shard/packed") {
                continue;
            }
            std::fs::copy(&path, root.join(relative)).unwrap();
        }
    }
}

fn cache() -> Arc<BlockCache> {
    BlockCache::open(&CacheConfig::new(8 << 20, 4096), runtime()).unwrap()
}

fn cached(cache: &Arc<BlockCache>) -> ShardOptions {
    ShardOptions {
        cache: Some(Arc::clone(cache)),
        ..ShardOptions::default()
    }
}

/// The first reads of each object learn its ETag and send no If-Match; every later read sends
/// it. No HEAD is sent, and no object is asked for beyond one learning byte per file.
#[test]
fn first_reads_learn() {
    let directory = tempfile::tempdir().unwrap();
    let entries = entries_of(1);
    test_shard(directory.path(), &entries);
    let (_server, seen, config) = store(directory.path());
    let location = Location::s3("s3://bucket/shard", &config).unwrap();
    let shard = Shard::open_with(&location, &cached(&cache())).unwrap();
    read_all(&shard, &entries).unwrap();
    read_all(&shard, &entries).unwrap();
    let seen = seen.lock().unwrap().clone();
    assert!(seen.iter().all(|request| request.method == "GET"));
    let mut keys: Vec<&str> = seen.iter().map(|request| request.key.as_str()).collect();
    keys.sort_unstable();
    keys.dedup();
    for key in keys.into_iter().filter(|key| key.contains("packed/")) {
        let of_key: Vec<_> = seen.iter().filter(|request| request.key == key).collect();
        assert!(of_key[0].if_match.is_none(), "{key}");
        assert!(of_key[1..].iter().all(|r| r.if_match.is_some()), "{key}");
    }
}

/// An object replaced between two searches fails the next read of an uncached block with the
/// retryable error. The shard then stays unusable, so no read mixes the two versions. Opening
/// it again reads the new manifest and the new bytes, and no byte of the old version comes back.
#[test]
fn replaced_between_searches() {
    let directory = tempfile::tempdir().unwrap();
    let old = entries_of(1);
    test_shard(directory.path(), &old);
    let (_server, _seen, config) = store(directory.path());
    let location = Location::s3("s3://bucket/shard", &config).unwrap();
    let cache = cache();
    let shard = Shard::open_with(&location, &cached(&cache)).unwrap();
    read_range(&shard, &old, 0, 0..1_000).unwrap();
    let new = entries_of(2);
    replace(directory.path(), &new, false);
    let error = read_range(&shard, &new, 0, 150_000..200_000).unwrap_err();
    assert!(error.to_string().contains("retry"), "{error}");
    assert_eq!(config.changes(), 1);
    assert!(read_range(&shard, &old, 0, 0..1_000).is_err());
    // The seeds file was replaced too, and its old ETag is still known, so the first open
    // learns that; the next open works.
    assert!(Shard::open_with(&location, &cached(&cache)).is_err());
    let shard = Shard::open_with(&location, &cached(&cache)).unwrap();
    read_all(&shard, &new).unwrap();
    assert_eq!(config.changes(), 2);
}

/// A later run finds the ETags in the cache directory, so its cached blocks are used, and its
/// check read is the only request of a file. The same holds when the object was replaced: the
/// check read finds it, though every block of the old version is cached.
#[test]
fn log_keeps_etags() {
    let directory = tempfile::tempdir().unwrap();
    let old = entries_of(1);
    test_shard(directory.path(), &old);
    let kept = directory.path().join("kept");
    let (_server, seen, config) = store(directory.path());
    config.keep_etags(&kept).unwrap();
    let location = Location::s3("s3://bucket/shard", &config).unwrap();
    let cache = cache();
    let (manifest, digest) = location.manifest().unwrap();
    let options = cached(&cache);
    let shard = Shard::open_with_manifest(&location, &options, manifest.clone(), digest.clone());
    read_all(&shard.unwrap(), &old).unwrap();
    assert!(seen.lock().unwrap().iter().any(|r| r.if_match.is_none()));
    // A next run: a new config over the same store, the same cached blocks.
    let next_run = || {
        let (server, seen, config) = store(directory.path());
        config.keep_etags(&kept).unwrap();
        let location = Location::s3("s3://bucket/shard", &config).unwrap();
        let (manifest, digest) = (manifest.clone(), digest.clone());
        let shard = Shard::open_with_manifest(&location, &options, manifest, digest).unwrap();
        (server, seen, shard)
    };
    let (_server, seen, shard) = next_run();
    read_all(&shard, &old).unwrap();
    let requests = seen.lock().unwrap().clone();
    assert!(
        requests.iter().all(|r| r.if_match.is_some()),
        "{requests:?}"
    );
    let packed = requests.iter().filter(|r| r.key.contains("packed/"));
    assert_eq!(packed.count(), 3, "one check read per file: {requests:?}");
    replace(directory.path(), &entries_of(2), true);
    let (_server, seen, shard) = next_run();
    let error = read_all(&shard, &old).unwrap_err();
    assert!(error.to_string().contains("retry"), "{error}");
    let requests = seen.lock().unwrap().clone();
    assert!(
        requests
            .iter()
            .filter(|r| r.key.contains("packed/"))
            .count()
            <= 3
    );
}
