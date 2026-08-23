//! Small in-process HTTP fixtures for remote resource integration tests.
//!
//! The production resource layer deliberately uses `curl`, so these fixtures
//! speak just enough HTTP/1.1 to exercise range requests, full-download
//! fallback, retry behavior, and streaming without network access.

use std::io::{Read, Write};
use std::net::{SocketAddr, TcpListener, TcpStream};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};
use std::time::Duration;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RequestRecord {
    pub method: String,
    pub path: String,
    pub range: Option<(usize, usize)>,
}

#[derive(Clone, Debug)]
struct FixtureState {
    body: Vec<u8>,
    etag: String,
    accepts_ranges: bool,
    fail_all_requests: usize,
    fail_range_requests: usize,
    truncate_range_after: Option<usize>,
    requests: Vec<RequestRecord>,
}

/// A local HTTP server with mutable response behavior.
pub struct HttpFixture {
    address: SocketAddr,
    state: Arc<Mutex<FixtureState>>,
    stop: Arc<AtomicBool>,
    thread: Option<JoinHandle<()>>,
}

#[allow(dead_code)]
impl HttpFixture {
    pub fn new(body: impl Into<Vec<u8>>) -> Self {
        Self::with_state(FixtureState {
            body: body.into(),
            etag: "fixture-v1".to_string(),
            accepts_ranges: true,
            fail_all_requests: 0,
            fail_range_requests: 0,
            truncate_range_after: None,
            requests: Vec::new(),
        })
    }

    pub fn without_ranges(self) -> Self {
        self.state.lock().unwrap().accepts_ranges = false;
        self
    }

    pub fn fail_range_requests(self, count: usize) -> Self {
        self.state.lock().unwrap().fail_range_requests = count;
        self
    }

    pub fn fail_all_requests(self, count: usize) -> Self {
        self.state.lock().unwrap().fail_all_requests = count;
        self
    }

    /// Close a ranged response after this many body bytes.  A subsequent
    /// request can still succeed, which models an interrupted transfer.
    pub fn truncate_range_after(self, bytes: usize) -> Self {
        self.state.lock().unwrap().truncate_range_after = Some(bytes);
        self
    }

    pub fn set_body(&self, body: impl Into<Vec<u8>>, etag: impl Into<String>) {
        let mut state = self.state.lock().unwrap();
        state.body = body.into();
        state.etag = etag.into();
    }

    pub fn set_accepts_ranges(&self, accepts_ranges: bool) {
        self.state.lock().unwrap().accepts_ranges = accepts_ranges;
    }

    pub fn url(&self, path: &str) -> String {
        format!("http://{}{}", self.address, path)
    }

    pub fn requests(&self) -> Vec<RequestRecord> {
        self.state.lock().unwrap().requests.clone()
    }

    pub fn range_requests(&self) -> Vec<RequestRecord> {
        self.requests()
            .into_iter()
            .filter(|request| request.range.is_some())
            .collect()
    }

    fn with_state(state: FixtureState) -> Self {
        let listener = TcpListener::bind("127.0.0.1:0").expect("bind HTTP fixture");
        listener
            .set_nonblocking(true)
            .expect("configure HTTP fixture listener");
        let address = listener.local_addr().expect("HTTP fixture address");
        let state = Arc::new(Mutex::new(state));
        let stop = Arc::new(AtomicBool::new(false));
        let thread_state = Arc::clone(&state);
        let thread_stop = Arc::clone(&stop);
        let thread = thread::spawn(move || {
            while !thread_stop.load(Ordering::Acquire) {
                match listener.accept() {
                    Ok((stream, _)) => handle_connection(stream, &thread_state),
                    Err(error) if error.kind() == std::io::ErrorKind::WouldBlock => {
                        thread::sleep(Duration::from_millis(1));
                    }
                    Err(_) => break,
                }
            }
        });
        Self {
            address,
            state,
            stop,
            thread: Some(thread),
        }
    }
}

impl Drop for HttpFixture {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::Release);
        // Wake a non-blocking accept loop before joining it.
        let _ = TcpStream::connect(self.address);
        if let Some(thread) = self.thread.take() {
            let _ = thread.join();
        }
    }
}

fn handle_connection(mut stream: TcpStream, state: &Arc<Mutex<FixtureState>>) {
    let _ = stream.set_read_timeout(Some(Duration::from_secs(5)));
    let mut request = Vec::new();
    let mut buffer = [0_u8; 1024];
    loop {
        match stream.read(&mut buffer) {
            Ok(0) => return,
            Ok(count) => {
                request.extend_from_slice(&buffer[..count]);
                if request.windows(4).any(|window| window == b"\r\n\r\n") {
                    break;
                }
            }
            Err(_) => return,
        }
    }

    let request_text = String::from_utf8_lossy(&request);
    let mut request_line = request_text.lines();
    let Some(first_line) = request_line.next() else {
        return;
    };
    let mut parts = first_line.split_whitespace();
    let Some(method) = parts.next() else {
        return;
    };
    let Some(path) = parts.next() else {
        return;
    };
    let range = request_text.lines().find_map(parse_range_header);
    let record = RequestRecord {
        method: method.to_string(),
        path: path.to_string(),
        range,
    };

    let response = {
        let mut fixture = state.lock().unwrap();
        fixture.requests.push(record);
        if fixture.fail_all_requests > 0 {
            fixture.fail_all_requests -= 1;
            None
        } else if fixture.fail_range_requests > 0 && range.is_some() {
            fixture.fail_range_requests -= 1;
            None
        } else {
            Some((
                fixture.body.clone(),
                fixture.etag.clone(),
                fixture.accepts_ranges,
                fixture.truncate_range_after,
            ))
        }
    };
    let Some((body, etag, accepts_ranges, truncate_range_after)) = response else {
        // Dropping the connection makes curl report a transport failure and
        // lets ObjectResource exercise its retry loop.
        return;
    };

    let (status, selected, content_range, total_length) = if let (true, Some((start, end))) =
        (method.eq_ignore_ascii_case("GET") && accepts_ranges, range)
    {
        if start >= body.len() || end < start || end >= body.len() {
            let _ = write!(
                stream,
                "HTTP/1.1 416 Range Not Satisfiable\r\nContent-Length: 0\r\nConnection: close\r\n\r\n"
            );
            return;
        }
        (
            "206 Partial Content",
            body[start..=end].to_vec(),
            Some((start, end)),
            body.len(),
        )
    } else {
        let length = body.len();
        ("200 OK", body, None, length)
    };

    let body_length = selected.len();
    let mut headers =
        format!("HTTP/1.1 {status}\r\nContent-Length: {body_length}\r\nETag: \"{etag}\"\r\n");
    if accepts_ranges {
        headers.push_str("Accept-Ranges: bytes\r\n");
    }
    if let Some((start, end)) = content_range {
        headers.push_str(&format!(
            "Content-Range: bytes {start}-{end}/{total_length}\r\n"
        ));
    }
    headers.push_str("Connection: close\r\n\r\n");
    if stream.write_all(headers.as_bytes()).is_err() {
        return;
    }
    if method.eq_ignore_ascii_case("HEAD") {
        let _ = stream.flush();
        return;
    }
    let bytes_to_send = truncate_range_after
        .filter(|_| content_range.is_some())
        .map_or(selected.len(), |limit| limit.min(selected.len()));
    let _ = stream.write_all(&selected[..bytes_to_send]);
    let _ = stream.flush();
}

fn parse_range_header(line: &str) -> Option<(usize, usize)> {
    let value = line.strip_prefix("Range: bytes=")?.trim();
    let (start, end) = value.split_once('-')?;
    Some((start.parse().ok()?, end.parse().ok()?))
}
