//! Tests of the S3 reader against a small HTTP server on a local port, which plays the store.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::*;
use crate::storage::cache::{BlockCache, CacheConfig};
use crate::test_support::runtime;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::net::TcpListener;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

/// A request as the server saw it.
struct Request {
    method: String,
    path: String,
    headers: HashMap<String, String>,
}

impl Request {
    /// The first and last byte of the requested range.
    fn range(&self) -> (usize, usize) {
        let value = self.headers["range"].strip_prefix("bytes=").unwrap();
        let (start, last) = value.split_once('-').unwrap();
        (start.parse().unwrap(), last.parse().unwrap())
    }
}

/// An answer: status, headers, the length it declares and the body.
struct Reply {
    status: u16,
    headers: Vec<(String, String)>,
    declared: usize,
    body: Vec<u8>,
}

impl Reply {
    fn status(status: u16) -> Self {
        let headers = Vec::new();
        let body = Vec::new();
        Self {
            status,
            headers,
            declared: 0,
            body,
        }
    }

    /// The answer to a range request for `data` with `etag`, which honours `If-Match`.
    fn range_of(data: &[u8], etag: &str, request: &Request) -> Self {
        if request
            .headers
            .get("if-match")
            .is_some_and(|wanted| wanted != etag)
        {
            return Self::status(412);
        }
        let (start, last) = request.range();
        let body = data[start..=last].to_vec();
        let range = format!("bytes {start}-{last}/{}", data.len());
        Self {
            status: 206,
            headers: vec![
                ("Content-Range".into(), range),
                ("ETag".into(), etag.into()),
            ],
            declared: body.len(),
            body,
        }
    }
}

/// A server that answers each request with `handler(index, request)`, where `index` counts
/// requests from 0. It stops when dropped.
struct Server {
    address: String,
    requests: Arc<AtomicUsize>,
    stop: Arc<AtomicBool>,
}

impl Drop for Server {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::Relaxed);
    }
}

fn serve(handler: impl Fn(usize, &Request) -> Reply + Send + Sync + 'static) -> Server {
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    listener.set_nonblocking(true).unwrap();
    let server = Server {
        address: format!("http://{}", listener.local_addr().unwrap()),
        requests: Arc::new(AtomicUsize::new(0)),
        stop: Arc::new(AtomicBool::new(false)),
    };
    let (requests, stop) = (Arc::clone(&server.requests), Arc::clone(&server.stop));
    let handler = Arc::new(handler);
    std::thread::spawn(move || {
        while !stop.load(Ordering::Relaxed) {
            let Ok((stream, _)) = listener.accept() else {
                std::thread::sleep(Duration::from_millis(1));
                continue;
            };
            let index = requests.fetch_add(1, Ordering::SeqCst);
            let handler = Arc::clone(&handler);
            std::thread::spawn(move || answer(stream, index, &*handler));
        }
    });
    server
}

fn answer(
    mut stream: std::net::TcpStream,
    index: usize,
    handler: &dyn Fn(usize, &Request) -> Reply,
) {
    stream.set_nonblocking(false).unwrap();
    let mut reader = BufReader::new(stream.try_clone().unwrap());
    let mut first = String::new();
    reader.read_line(&mut first).unwrap();
    let mut parts = first.split_whitespace();
    let method = parts.next().unwrap().to_string();
    let path = parts.next().unwrap().to_string();
    let mut headers = HashMap::new();
    loop {
        let mut line = String::new();
        reader.read_line(&mut line).unwrap();
        let Some((name, value)) = line.trim_end().split_once(':') else {
            break;
        };
        headers.insert(name.to_ascii_lowercase(), value.trim().to_string());
    }
    let request = Request {
        method,
        path,
        headers,
    };
    let reply = handler(index, &request);
    let mut head = format!(
        "HTTP/1.1 {} reply\r\nContent-Length: {}\r\nConnection: close\r\n",
        reply.status, reply.declared
    );
    for (name, value) in &reply.headers {
        head.push_str(&format!("{name}: {value}\r\n"));
    }
    // The client may give up and close the connection, so write errors are fine.
    let _ = stream.write_all(format!("{head}\r\n").as_bytes());
    let _ = stream.write_all(&reply.body);
}

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

/// A store with one object `bucket/entry.packed`; after `changed` requests its ETag changes.
fn store(data: Vec<u8>, changed: usize) -> Server {
    serve(move |index, request| {
        assert_eq!(request.method, "GET", "the size must not come from a HEAD");
        assert_eq!(request.path, "/bucket/entry.packed");
        let etag = if index >= changed {
            "\"other\""
        } else {
            "\"fixture\""
        };
        Reply::range_of(&data, etag, request)
    })
}

fn pattern(length: u32) -> Vec<u8> {
    (0..length).map(|value| (value % 251) as u8).collect()
}

#[test]
fn exact_ranges_read() {
    let data = pattern(50_000);
    let server = store(data.clone(), usize::MAX);
    let location = "s3://bucket/entry.packed";
    let object = S3Object::open(&config(&server, quick()), location, Some(50_000)).unwrap();
    assert_eq!(object.size(), 50_000);
    for (offset, length) in [(0usize, 10usize), (49_990, 10), (12_345, 1000)] {
        let mut buffer = vec![0; length];
        object.read_exact_at(offset as u64, &mut buffer).unwrap();
        assert_eq!(buffer, data[offset..offset + length]);
    }
    assert!(object.read_exact_at(49_999, &mut [0; 2]).is_err());
    object.read_exact_at(0, &mut []).unwrap();
    // The opening request and three reads.
    assert_eq!(server.requests.load(Ordering::SeqCst), 4);
}

#[test]
fn cache_merges_runs() {
    let data = pattern(50_000);
    let server = store(data.clone(), usize::MAX);
    let location = "s3://bucket/entry.packed";
    let object = S3Object::open(&config(&server, quick()), location, None).unwrap();
    let cache = BlockCache::open(&CacheConfig::new(1 << 20, 4096), runtime()).unwrap();
    let cached = cache.wrap(object);
    let mut buffer = vec![0; 3000];
    for _ in 0..2 {
        cached.read_exact_at(20_000, &mut buffer).unwrap();
        assert_eq!(buffer, data[20_000..23_000]);
        cached.read_exact_at(40_000, &mut buffer[..100]).unwrap();
    }
    // The opening request and one fetch per missing run.
    assert_eq!(server.requests.load(Ordering::SeqCst), 3);
}

#[test]
fn changed_object_fails() {
    let server = store(vec![7u8; 1000], 1);
    let config = config(&server, quick());
    let object = S3Object::open(&config, "s3://bucket/entry.packed", Some(1000)).unwrap();
    let error = object.read_exact_at(0, &mut [0; 10]).unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
    let wrong_size = S3Object::open(&config, "s3://bucket/entry.packed", Some(999));
    assert!(matches!(wrong_size, Err(StorageError::Object { .. })));
    for location in [
        "s3:/bucket",
        "s3://bucket",
        "s3://bucket/",
        "s3:///key",
        "b/key",
    ] {
        let opened = S3Object::open(&config, location, None);
        assert!(
            matches!(opened, Err(StorageError::Location(_))),
            "{location}"
        );
    }
}

#[test]
fn weak_objects_refused() {
    let weak = serve(|_, request| Reply::range_of(&[0; 10], "W/\"weak\"", request));
    assert!(S3Object::open(&config(&weak, quick()), "s3://bucket/key", None).is_err());
    let forbidden = serve(|_, _| Reply::status(403));
    assert!(S3Object::open(&config(&forbidden, quick()), "s3://bucket/key", None).is_err());
    // A final answer is not repeated.
    assert_eq!(forbidden.requests.load(Ordering::SeqCst), 1);
}

#[test]
fn short_body_fails() {
    let data = pattern(1000);
    let server = serve(move |index, request| {
        let mut reply = Reply::range_of(&data, "\"fixture\"", request);
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
    let object = S3Object::open(&config(&server, limits), "s3://bucket/key", None).unwrap();
    assert!(object.read_exact_at(0, &mut [0; 10]).is_err());
}

#[test]
fn in_flight_limited() {
    let data = pattern(100_000);
    let (active, peak) = (Arc::new(AtomicUsize::new(0)), Arc::new(AtomicUsize::new(0)));
    let server = serve({
        let (data, active, peak) = (data.clone(), Arc::clone(&active), Arc::clone(&peak));
        move |_, request| {
            let now = active.fetch_add(1, Ordering::SeqCst) + 1;
            peak.fetch_max(now, Ordering::SeqCst);
            std::thread::sleep(Duration::from_millis(30));
            active.fetch_sub(1, Ordering::SeqCst);
            Reply::range_of(&data, "\"fixture\"", request)
        }
    });
    let limits = S3Limits {
        max_in_flight: 3,
        ..quick()
    };
    let location = "s3://bucket/entry.packed";
    let object = S3Object::open(&config(&server, limits), location, None).unwrap();
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
    let peak = peak.load(Ordering::SeqCst);
    assert!(peak > 1 && peak <= 3, "peak {peak}");
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
