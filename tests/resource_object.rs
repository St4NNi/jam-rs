use jam_rs::resource::object::ObjectResource;
use jam_rs::resource::{ByteRange, RangeReader, ResourceOpenOptions};
#[path = "support/remote/mod.rs"]
mod remote;
use std::io::{Read, Write};
use std::net::TcpListener;
use std::thread;

#[test]
fn http_ranges_use_a_bounded_cache() {
    // The object reader uses curl for HTTP transport.  Keep this test local and
    // deterministic so it exercises request validation without internet I/O.
    if std::process::Command::new("curl")
        .arg("--version")
        .output()
        .is_err()
    {
        return;
    }
    let body = b"0123456789abcdef".to_vec();
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let address = listener.local_addr().unwrap();
    let server_body = body.clone();
    let server = thread::spawn(move || {
        for connection in listener.incoming().take(4) {
            let mut stream = connection.unwrap();
            let mut request = Vec::new();
            let mut buffer = [0_u8; 1024];
            loop {
                let read = stream.read(&mut buffer).unwrap();
                if read == 0 {
                    break;
                }
                request.extend_from_slice(&buffer[..read]);
                if request.windows(4).any(|window| window == b"\r\n\r\n") {
                    break;
                }
            }
            let request = String::from_utf8_lossy(&request);
            if request.starts_with("HEAD ") {
                write!(
                    stream,
                    "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nAccept-Ranges: bytes\r\nETag: \"fixture\"\r\nConnection: close\r\n\r\n",
                    server_body.len()
                )
                .unwrap();
                continue;
            }
            let (start, end) = request
                .lines()
                .find_map(|line| line.strip_prefix("Range: bytes="))
                .and_then(|value| value.trim().split_once('-'))
                .map(|(start, end)| {
                    (
                        start.parse::<usize>().unwrap(),
                        end.parse::<usize>().unwrap(),
                    )
                })
                .unwrap();
            let selected = &server_body[start..=end];
            write!(
                stream,
                "HTTP/1.1 206 Partial Content\r\nContent-Length: {}\r\nContent-Range: bytes {}-{}/{}\r\nAccept-Ranges: bytes\r\nETag: \"fixture\"\r\nConnection: close\r\n\r\n",
                selected.len(),
                start,
                end,
                server_body.len()
            )
            .unwrap();
            stream.write_all(selected).unwrap();
        }
    });

    let options = ResourceOpenOptions {
        cache_block_bytes: 4,
        max_cache_bytes: 16,
        max_retries: 0,
        ..ResourceOpenOptions::default()
    };
    let resource = ObjectResource::open(format!("http://{address}/fixture"), options).unwrap();
    assert_eq!(
        resource.read_range(ByteRange::new(3, 7).unwrap()).unwrap(),
        b"3456789"
    );
    let first = resource.metrics();
    assert_eq!(
        resource.read_range(ByteRange::new(3, 7).unwrap()).unwrap(),
        b"3456789"
    );
    let second = resource.metrics();
    assert!(second.cache_hits > first.cache_hits);
    assert!(second.remote_bytes > 0);
    server.join().unwrap();
}

#[test]
fn ranged_reads_require_an_exact_content_range() {
    let fixture = remote_fixture(b"0123456789abcdef".to_vec()).with_malformed_content_range();
    let resource = ObjectResource::open(
        fixture.url("/fixture"),
        ResourceOpenOptions {
            cache_block_bytes: 16,
            max_cache_bytes: 64,
            max_retries: 0,
            ..ResourceOpenOptions::default()
        },
    )
    .unwrap();
    let error = resource
        .read_range(ByteRange::new(3, 4).unwrap())
        .expect_err("mismatched Content-Range must be rejected");
    assert!(matches!(
        error,
        jam_rs::resource::ResourceError::Transport { .. }
    ));
    let metrics = resource.metrics();
    assert!(metrics.head_requests >= 1);
    assert!(metrics.get_requests >= 1);
    assert!(metrics.requested_bytes >= 4);
    assert!(metrics.returned_bytes >= 4);
}

#[test]
fn full_object_fallback_is_counted_and_bounded() {
    let fixture = remote_fixture(b"0123456789abcdef".to_vec()).ignoring_ranges();
    let allowed = ObjectResource::open(
        fixture.url("/fixture"),
        ResourceOpenOptions {
            cache_block_bytes: 16,
            max_cache_bytes: 32,
            max_retries: 0,
            ..ResourceOpenOptions::default()
        },
    )
    .unwrap();
    assert_eq!(
        allowed.read_range(ByteRange::new(3, 4).unwrap()).unwrap(),
        b"3456"
    );
    let allowed_metrics = allowed.metrics();
    assert_eq!(allowed_metrics.full_object_fallbacks, 1);
    assert_eq!(allowed_metrics.returned_bytes, 16);
    assert_eq!(allowed_metrics.requested_bytes, 16);

    let bounded_fixture = remote_fixture(b"0123456789abcdef".to_vec()).ignoring_ranges();
    let bounded = ObjectResource::open(
        bounded_fixture.url("/fixture"),
        ResourceOpenOptions {
            cache_block_bytes: 4,
            max_cache_bytes: 8,
            max_retries: 0,
            ..ResourceOpenOptions::default()
        },
    )
    .unwrap();
    assert!(bounded.read_range(ByteRange::new(3, 4).unwrap()).is_err());
    assert_eq!(bounded.metrics().full_object_fallbacks, 0);
}

#[test]
fn http_statuses_are_structured_and_retry_only_for_server_errors() {
    let forbidden = remote_fixture(b"unavailable".to_vec()).with_status(403);
    let resource = ObjectResource::open(
        forbidden.url("/forbidden?token=secret"),
        ResourceOpenOptions {
            max_retries: 4,
            ..ResourceOpenOptions::default()
        },
    )
    .unwrap();
    let error = resource.metadata().expect_err("403 must fail");
    assert!(matches!(
        error,
        jam_rs::resource::ResourceError::HttpStatus { status: 403, .. }
    ));
    assert_eq!(forbidden.requests().len(), 1);

    let unavailable = remote_fixture(b"unavailable".to_vec()).with_status(503);
    let resource = ObjectResource::open(
        unavailable.url("/unavailable"),
        ResourceOpenOptions {
            max_retries: 2,
            ..ResourceOpenOptions::default()
        },
    )
    .unwrap();
    let error = resource.metadata().expect_err("503 must fail");
    assert!(matches!(
        error,
        jam_rs::resource::ResourceError::HttpStatus { status: 503, .. }
    ));
    assert_eq!(unavailable.requests().len(), 3);
    assert_eq!(resource.metrics().retries, 2);
}

fn remote_fixture(body: Vec<u8>) -> remote::HttpFixture {
    // Keep this test file self-contained while sharing the canonical fixture
    // behavior used by the trace remote integration tests.
    remote::HttpFixture::new(body)
}
