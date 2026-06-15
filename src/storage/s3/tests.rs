//! Tests against a small HTTP server on a local port, which plays the store.

use super::*;
use crate::storage::cache::{BlockCache, CacheConfig};
use crate::storage::{ReadAt, StorageError};
use std::collections::HashMap;
use std::io::{self, BufRead, BufReader, Write};
use std::net::TcpListener;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::time::{Duration, Instant};

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

/// An answer: status, headers, the length it declares, the body, and a wait between the two.
struct Reply {
    status: u16,
    headers: Vec<(String, String)>,
    declared: usize,
    body: Vec<u8>,
    stall: Duration,
}

impl Reply {
    fn status(status: u16) -> Self {
        Self {
            status,
            headers: Vec::new(),
            declared: 0,
            body: Vec::new(),
            stall: Duration::ZERO,
        }
    }

    /// The right answer to a range request for `data`, with `etag`.
    fn range_of(data: &[u8], etag: &str, request: &Request) -> Self {
        let (start, last) = request.range();
        let body = data[start..=last].to_vec();
        Self {
            status: 206,
            headers: vec![
                (
                    "Content-Range".into(),
                    format!("bytes {start}-{last}/{}", data.len()),
                ),
                ("ETag".into(), etag.into()),
            ],
            declared: body.len(),
            body,
            stall: Duration::ZERO,
        }
    }
}

/// A server that answers each request with the reply of `handler(index, request)`, where `index`
/// counts requests from 0. It stops when dropped.
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
    let (requests, stop, handler) = (
        Arc::clone(&server.requests),
        Arc::clone(&server.stop),
        Arc::new(handler),
    );
    std::thread::spawn(move || {
        while !stop.load(Ordering::Relaxed) {
            let Ok((stream, _)) = listener.accept() else {
                std::thread::sleep(Duration::from_millis(1));
                continue;
            };
            let (handler, index) = (
                Arc::clone(&handler),
                requests.fetch_add(1, Ordering::SeqCst),
            );
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
    let (method, path) = (
        parts.next().unwrap().to_string(),
        parts.next().unwrap().to_string(),
    );
    let mut headers = HashMap::new();
    loop {
        let mut line = String::new();
        reader.read_line(&mut line).unwrap();
        let Some((name, value)) = line.trim_end().split_once(':') else {
            break;
        };
        headers.insert(name.to_ascii_lowercase(), value.trim().to_string());
    }
    let reply = handler(
        index,
        &Request {
            method,
            path,
            headers,
        },
    );
    let mut head = format!(
        "HTTP/1.1 {} reply\r\nContent-Length: {}\r\nConnection: close\r\n",
        reply.status, reply.declared
    );
    for (name, value) in &reply.headers {
        head.push_str(&format!("{name}: {value}\r\n"));
    }
    // The client may give up and close the connection, so write errors are fine.
    let _ = stream.write_all(format!("{head}\r\n").as_bytes());
    std::thread::sleep(reply.stall);
    let _ = stream.write_all(&reply.body);
}

fn credentials() -> Credentials {
    Credentials {
        access_key: "access".into(),
        secret_key: "secret".into(),
        session_token: None,
    }
}

fn config(server: &Server, limits: S3Limits) -> S3Config {
    S3Config::new("test", Some(&server.address), true, credentials(), limits).unwrap()
}

/// Limits that wait only milliseconds between tries.
fn quick() -> S3Limits {
    S3Limits {
        backoff: Duration::from_millis(5),
        ..S3Limits::default()
    }
}

/// A store with one object `bucket/entry.packed` whose ETag is `etag`; after `changed`
/// requests the ETag is another one.
fn store(data: Vec<u8>, changed: usize) -> Server {
    serve(move |index, request| {
        assert_eq!(
            request.method, "GET",
            "the size must not come from a HEAD request"
        );
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
fn exact_ranges_are_read_signed_and_counted() {
    let data = pattern(50_000);
    let server = serve({
        let data = data.clone();
        move |index, request| {
            let authorization = &request.headers["authorization"];
            assert!(authorization.starts_with("AWS4-HMAC-SHA256 Credential=access/"));
            assert!(authorization.contains("SignedHeaders=host;"));
            assert!(request.headers.contains_key("x-amz-date"));
            if index > 0 {
                assert_eq!(request.headers["if-match"], "\"fixture\"");
            }
            Reply::range_of(&data, "\"fixture\"", request)
        }
    });
    let object = S3Object::open(
        &config(&server, quick()),
        "s3://bucket/entry.packed",
        Some(50_000),
    )
    .unwrap();
    assert_eq!(object.size(), 50_000);
    for (offset, length) in [(0usize, 10usize), (49_990, 10), (12_345, 1000)] {
        let mut buffer = vec![0; length];
        object.read_exact_at(offset as u64, &mut buffer).unwrap();
        assert_eq!(buffer, data[offset..offset + length]);
    }
    assert!(object.read_exact_at(49_999, &mut [0; 2]).is_err());
    object.read_exact_at(0, &mut []).unwrap();
    assert_eq!(object.requests(), 4);
}

#[test]
fn a_cache_in_front_needs_one_request_per_missing_run() {
    let data = pattern(50_000);
    let server = store(data.clone(), usize::MAX);
    let object =
        S3Object::open(&config(&server, quick()), "s3://bucket/entry.packed", None).unwrap();
    let cache = BlockCache::open(&CacheConfig::new(1 << 20, 4096)).unwrap();
    let cached = cache.wrap(object);
    let mut buffer = vec![0; 3000];
    for _ in 0..2 {
        cached.read_exact_at(20_000, &mut buffer).unwrap();
        assert_eq!(buffer, data[20_000..23_000]);
        cached.read_exact_at(40_000, &mut buffer[..100]).unwrap();
    }
    assert_eq!(cache.stats().requests, 2);
    // The opening request and the two fetches.
    assert_eq!(server.requests.load(Ordering::SeqCst), 3);
}

#[test]
fn a_changed_object_and_a_wrong_size_fail() {
    let server = store(vec![7u8; 1000], 1);
    let config = config(&server, quick());
    let object = S3Object::open(&config, "s3://bucket/entry.packed", Some(1000)).unwrap();
    let error = object.read_exact_at(0, &mut [0; 10]).unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
    assert!(S3Object::open(&config, "s3://bucket/entry.packed", Some(999)).is_err());
    for location in [
        "s3:/bucket",
        "s3://bucket",
        "s3://bucket/",
        "s3:///key",
        "bucket/key",
    ] {
        assert!(
            matches!(
                S3Object::open(&config, location, None),
                Err(StorageError::Location(_))
            ),
            "{location}"
        );
    }
}

#[test]
fn objects_without_a_strong_etag_or_size_are_refused() {
    let weak = serve(|_, request| Reply::range_of(&[0; 10], "W/\"weak\"", request));
    assert!(S3Object::open(&config(&weak, quick()), "s3://bucket/key", None).is_err());
    let plain = serve(|_, _| Reply {
        status: 200,
        declared: 1,
        body: vec![0],
        ..Reply::status(200)
    });
    assert!(S3Object::open(&config(&plain, quick()), "s3://bucket/key", None).is_err());
    let forbidden = serve(|_, _| Reply::status(403));
    assert!(S3Object::open(&config(&forbidden, quick()), "s3://bucket/key", None).is_err());
    // A final answer is not repeated.
    assert_eq!(forbidden.requests.load(Ordering::SeqCst), 1);
}

#[test]
fn temporary_errors_are_retried_and_then_given_up() {
    let data = vec![9u8; 1000];
    // The opening request and the read each answer 503 first and the range at the second try.
    let server = serve(move |index, request| {
        if index % 2 == 0 {
            Reply::status(503)
        } else {
            Reply::range_of(&data, "\"fixture\"", request)
        }
    });
    let object = S3Object::open(
        &config(&server, quick()),
        "s3://bucket/entry.packed",
        Some(1000),
    )
    .unwrap();
    let mut buffer = [0; 10];
    object.read_exact_at(0, &mut buffer).unwrap();
    assert_eq!(buffer, [9; 10]);
    assert_eq!(server.requests.load(Ordering::SeqCst), 4);

    let busy = serve(|_, _| Reply::status(503));
    assert!(S3Object::open(&config(&busy, quick()), "s3://bucket/key", None).is_err());
    assert_eq!(busy.requests.load(Ordering::SeqCst), 4);
}

#[test]
fn an_answer_for_another_range_is_not_repeated() {
    let data = pattern(1000);
    // After the opening request the store answers a wrong range with status 206 first.
    let server = serve(move |index, request| match index {
        1 => Reply::range_of(
            &data,
            "\"fixture\"",
            &Request {
                method: "GET".into(),
                path: request.path.clone(),
                headers: HashMap::from([("range".to_string(), "bytes=1-10".to_string())]),
            },
        ),
        _ => Reply::range_of(&data, "\"fixture\"", request),
    });
    let object =
        S3Object::open(&config(&server, quick()), "s3://bucket/entry.packed", None).unwrap();
    let error = object.read_exact_at(0, &mut [0; 10]).unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
    // Opening, the wrong answer: not repeated.
    assert_eq!(server.requests.load(Ordering::SeqCst), 2);
}

#[test]
fn requests_in_flight_stay_within_the_limit() {
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
    let object =
        S3Object::open(&config(&server, limits), "s3://bucket/entry.packed", None).unwrap();
    std::thread::scope(|scope| {
        for thread in 0..12usize {
            let (object, data) = (&object, &data);
            scope.spawn(move || {
                let mut buffer = vec![0; 100];
                object
                    .read_exact_at(thread as u64 * 1000, &mut buffer)
                    .unwrap();
                assert_eq!(buffer, data[thread * 1000..thread * 1000 + 100]);
            });
        }
    });
    let peak = peak.load(Ordering::SeqCst);
    assert!(peak > 1 && peak <= 3, "peak {peak}");
}

#[test]
fn the_request_timeout_bounds_a_stalled_body() {
    for status in [206, 500] {
        let data = vec![1u8; 1000];
        let server = serve(move |index, request| {
            if index == 0 {
                return Reply::range_of(&data, "\"fixture\"", request);
            }
            // The headers promise a body that arrives too late.
            let mut reply = Reply::range_of(&data, "\"fixture\"", request);
            reply.status = status;
            reply.stall = Duration::from_millis(600);
            reply
        });
        let limits = S3Limits {
            request_timeout: Duration::from_millis(50),
            tries: 1,
            ..quick()
        };
        let object =
            S3Object::open(&config(&server, limits), "s3://bucket/entry.packed", None).unwrap();
        let started = Instant::now();
        assert!(object.read_exact_at(0, &mut [0; 10]).is_err());
        assert!(
            started.elapsed() < Duration::from_millis(400),
            "{:?}",
            started.elapsed()
        );
    }
}

#[test]
fn the_connect_timeout_bounds_a_stuck_connection() {
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let address = listener.local_addr().unwrap();
    // Nobody accepts: once the accept queue is full, a new connection waits for an answer to its
    // retried SYN, which the default connect timeout would allow for 30 s.
    let mut queued = Vec::new();
    while let Ok(stream) =
        std::net::TcpStream::connect_timeout(&address, Duration::from_millis(200))
    {
        queued.push(stream);
        assert!(queued.len() < 10_000);
    }
    let limits = S3Limits {
        connect_timeout: Duration::from_millis(100),
        tries: 1,
        ..quick()
    };
    let config = S3Config::new(
        "test",
        Some(&format!("http://{address}")),
        true,
        credentials(),
        limits,
    )
    .unwrap();
    let started = Instant::now();
    assert!(S3Object::open(&config, "s3://bucket/key", None).is_err());
    assert!(started.elapsed() < Duration::from_secs(2));
    drop(listener);
}
