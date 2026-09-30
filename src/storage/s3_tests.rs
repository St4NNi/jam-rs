//! Tests of the S3 reader against a small HTTP server on a local port, which plays the store.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::*;
use crate::storage::{BlockCache, CacheConfig, StorageCounts};
use crate::test_support::{Reply, Server, etag_of, runtime, serve};
use std::sync::Barrier;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

/// The digest the test object is cached under.
const DIGEST: &str = "1111111111111111111111111111111111111111111111111111111111111111";

/// Limits that wait only milliseconds between tries.
fn quick() -> S3Limits {
    S3Limits {
        backoff: Duration::from_millis(5),
        ..S3Limits::default()
    }
}

fn config(server: &Server, limits: S3Limits) -> S3Config {
    let address = Some(server.address.as_str());
    let config = S3Config::new("test", address, true, limits, runtime()).unwrap();
    config.unsigned()
}

/// A store with one object `bucket/entry.packed`, which refuses HEAD.
fn store(data: Vec<u8>) -> Server {
    serve(move |_, request| {
        assert_eq!(request.method, "GET", "the size must not come from a HEAD");
        assert_eq!(request.path, "/bucket/entry.packed");
        Reply::object(&data, request)
    })
}

fn pattern(length: u32) -> Vec<u8> {
    (0..length).map(|value| (value % 251) as u8).collect()
}

#[test]
fn exact_ranges_read() {
    let data = pattern(50_000);
    let server = store(data.clone());
    let location = "s3://bucket/entry.packed";
    let object = S3Object::open(&config(&server, quick()), location, 50_000).unwrap();
    assert_eq!(object.size(), 50_000);
    for (offset, length) in [(0usize, 10usize), (49_990, 10), (12_345, 1000)] {
        let mut buffer = vec![0; length];
        object.read_exact_at(offset as u64, &mut buffer).unwrap();
        assert_eq!(buffer, data[offset..offset + length]);
    }
    assert!(object.read_exact_at(49_999, &mut [0; 2]).is_err());
    object.read_exact_at(0, &mut []).unwrap();
    // Opening sends nothing; each read is one request.
    assert_eq!(server.requests.load(Ordering::SeqCst), 3);
    let whole = S3Object::read_whole(&config(&server, quick()), location).unwrap();
    assert_eq!(whole, data);
}

#[test]
fn cache_merges_runs() {
    let data = pattern(50_000);
    let server = store(data.clone());
    let location = "s3://bucket/entry.packed";
    let object = S3Object::open(&config(&server, quick()), location, 50_000).unwrap();
    let cache = BlockCache::open(&CacheConfig::new(1 << 20, 4096), runtime()).unwrap();
    let cached = cache.wrap(Arc::new(object), DIGEST, None).unwrap();
    let mut buffer = vec![0; 3000];
    for _ in 0..2 {
        cached.read_exact_at(20_000, &mut buffer).unwrap();
        assert_eq!(buffer, data[20_000..23_000]);
        cached.read_exact_at(40_000, &mut buffer[..100]).unwrap();
    }
    // One fetch per missing run.
    assert_eq!(server.requests.load(Ordering::SeqCst), 2);
}

/// The first read learns the ETag and sends no If-Match; later reads send it. No HEAD is sent.
#[test]
fn first_read_learns() {
    let data = pattern(50_000);
    let matches = Arc::new(Mutex::new(Vec::new()));
    let server = serve({
        let (data, matches) = (data.clone(), Arc::clone(&matches));
        move |_, request| {
            assert_eq!(request.method, "GET");
            lock(&matches).push(request.headers.get("if-match").cloned());
            Reply::object(&data, request)
        }
    });
    let config = config(&server, quick());
    let object = S3Object::open(&config, "s3://bucket/entry.packed", 50_000).unwrap();
    object.read_exact_at(0, &mut [0; 10]).unwrap();
    object.read_exact_at(100, &mut [0; 10]).unwrap();
    let other = S3Object::open(&config, "s3://bucket/entry.packed", 50_000).unwrap();
    other.read_exact_at(200, &mut [0; 10]).unwrap();
    let etag = Some(etag_of(&data));
    assert_eq!(*lock(&matches), vec![None, etag.clone(), etag]);
}

/// A replaced object fails one read with the retryable error; the next read learns the new
/// ETag and gives the new bytes. The hook gets the version of the old ETag.
#[test]
fn replaced_object_refused() {
    let old = pattern(50_000);
    let new: Vec<u8> = old.iter().map(|byte| byte ^ 1).collect();
    let replaced = Arc::new(AtomicBool::new(false));
    let server = serve({
        let (old, new, replaced) = (old.clone(), new.clone(), Arc::clone(&replaced));
        move |_, request| match replaced.load(Ordering::SeqCst) {
            true => Reply::object(&new, request),
            false => Reply::object(&old, request),
        }
    });
    let config = config(&server, quick());
    let seen = Arc::new(Mutex::new(Vec::new()));
    let hook: ChangeHook = {
        let seen = Arc::clone(&seen);
        Arc::new(move |version| lock(&seen).push(version))
    };
    let location = "s3://bucket/entry.packed";
    let object = S3Object::open(&config, location, 50_000)
        .unwrap()
        .on_change(hook);
    let mut buffer = vec![0; 100];
    object.read_exact_at(0, &mut buffer).unwrap();
    assert_eq!(buffer, old[..100]);
    replaced.store(true, Ordering::SeqCst);
    let error = object.read_exact_at(20_000, &mut buffer).unwrap_err();
    assert!(error.to_string().contains("changed in S3"), "{error}");
    assert!(error.to_string().contains("retry"), "{error}");
    assert_eq!(config.changes(), 1);
    assert_eq!(*lock(&seen), vec![version_of(&etag_of(&old))]);
    object.read_exact_at(20_000, &mut buffer).unwrap();
    assert_eq!(buffer, new[20_000..20_100]);
    assert_eq!(config.changes(), 1);
}

#[test]
fn bad_locations_refused() {
    let server = store(Vec::new());
    let config = config(&server, quick());
    for location in [
        "s3:/bucket",
        "s3://bucket",
        "s3://bucket/",
        "s3:///key",
        "b/key",
    ] {
        let opened = S3Object::open(&config, location, 1);
        assert!(
            matches!(opened, Err(StorageError::Location(_))),
            "{location}"
        );
        let read = S3Object::read_whole(&config, location);
        assert!(matches!(read, Err(StorageError::Location(_))), "{location}");
    }
}

#[test]
fn forbidden_not_repeated() {
    let forbidden = serve(|_, _| Reply::status(403));
    let config = config(&forbidden, quick());
    let object = S3Object::open(&config, "s3://bucket/key", 10).unwrap();
    assert!(object.read_exact_at(0, &mut [0; 4]).is_err());
    assert!(S3Object::read_whole(&config, "s3://bucket/key").is_err());
    // A final answer is not repeated.
    assert_eq!(forbidden.requests.load(Ordering::SeqCst), 2);
}

/// Every attempt counts as a request, retried, failed or not.
#[test]
fn retries_counted() {
    let data = pattern(1000);
    // Two server errors, then the object; then errors only.
    let server = serve(move |index, request| match index {
        2 => Reply::object(&data, request),
        _ => Reply::status(503),
    });
    let limits = S3Limits {
        tries: 3,
        ..quick()
    };
    let object = S3Object::open(&config(&server, limits), "s3://bucket/key", 1000).unwrap();
    let before = StorageCounts::now();
    object.read_exact_at(0, &mut [0; 10]).unwrap();
    assert!(object.read_exact_at(0, &mut [0; 10]).is_err());
    assert_eq!(server.requests.load(Ordering::SeqCst), 6);
    // Other tests of this process may count at the same time.
    let counted = StorageCounts::now().since(&before);
    assert!(
        counted.s3_requests >= 6 && counted.s3_bytes >= 10,
        "{counted:?}"
    );
}

#[test]
fn short_body_fails() {
    let data = pattern(1000);
    let server = serve(move |index, request| {
        let mut reply = Reply::object(&data, request);
        if index > 0 {
            reply.body.truncate(5);
            reply.declared = 5;
        }
        reply
    });
    let limits = S3Limits {
        tries: 1,
        ..quick()
    };
    let object = S3Object::open(&config(&server, limits), "s3://bucket/key", 1000).unwrap();
    object.read_exact_at(0, &mut [0; 10]).unwrap();
    assert!(object.read_exact_at(0, &mut [0; 10]).is_err());
}

#[test]
fn in_flight_limited() {
    let data = pattern(100_000);
    let (active, peak) = (Arc::new(AtomicUsize::new(0)), Arc::new(AtomicUsize::new(0)));
    // The twelve reads pass in groups of three, so they only finish when three run at once.
    let together = Arc::new(Barrier::new(3));
    let server = serve({
        let (data, active, peak) = (data.clone(), Arc::clone(&active), Arc::clone(&peak));
        move |_, request| {
            let now = active.fetch_add(1, Ordering::SeqCst) + 1;
            peak.fetch_max(now, Ordering::SeqCst);
            together.wait();
            active.fetch_sub(1, Ordering::SeqCst);
            Reply::object(&data, request)
        }
    });
    let limits = S3Limits {
        max_in_flight: 3,
        ..quick()
    };
    let location = "s3://bucket/entry.packed";
    let object = S3Object::open(&config(&server, limits), location, 100_000).unwrap();
    std::thread::scope(|scope| {
        for thread in 0..12usize {
            let (object, data) = (&object, &data);
            scope.spawn(move || {
                let mut buffer = vec![0; 100];
                let offset = thread as u64 * 1000;
                object.read_exact_at(offset, &mut buffer).unwrap();
                assert_eq!(buffer, data[thread * 1000..thread * 1000 + 100]);
            });
        }
    });
    assert_eq!(peak.load(Ordering::SeqCst), 3);
}

#[test]
fn bad_settings_rejected() {
    let limits = S3Limits::default();
    for endpoint in ["", "ftp://host", "host:9000", "http://"] {
        let config = S3Config::new("r", Some(endpoint), true, limits, runtime());
        assert!(config.is_err(), "{endpoint}");
    }
    assert!(S3Config::new("", None, true, limits, runtime()).is_err());
    let no_tries = S3Limits { tries: 0, ..limits };
    assert!(S3Config::new("r", None, true, no_tries, runtime()).is_err());
}
