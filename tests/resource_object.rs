use jam_rs::resource::object::ObjectResource;
use jam_rs::resource::{ByteRange, RangeReader, ResourceOpenOptions};
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
