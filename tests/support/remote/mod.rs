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
    /// HTTP response status. A zero status means that the fixture closed the
    /// connection before writing a response.
    pub status: u16,
    /// Response `Content-Range`, including the advertised object length.
    pub content_range: Option<(usize, usize, usize)>,
    /// Number of bytes requested by the client, or zero for an un-ranged
    /// request.
    pub requested_bytes: usize,
    /// Number of body bytes written by the fixture.
    pub returned_bytes: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct UploadRecord {
    pub method: String,
    pub path: String,
    pub body: Vec<u8>,
}

#[derive(Clone, Debug)]
struct FixtureState {
    body: Vec<u8>,
    etag: String,
    accepts_ranges: bool,
    /// Force a response status instead of normal range handling.  This is
    /// intentionally a single status so tests can model an unavailable
    /// object without introducing a second server implementation.
    response_status: Option<u16>,
    /// Return a 206 body with a malformed Content-Range header.
    malformed_content_range: bool,
    /// Ignore a requested range and return a 200 full-object response.
    ignore_range: bool,
    /// Delay the response long enough for a client with a short timeout to
    /// exercise its timeout path.
    response_delay: Option<Duration>,
    fail_all_requests: usize,
    fail_range_requests: usize,
    truncate_range_after: Option<usize>,
    upload_status: Option<u16>,
    interrupt_upload: bool,
    uploaded: Vec<UploadRecord>,
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
            response_status: None,
            malformed_content_range: false,
            ignore_range: false,
            response_delay: None,
            fail_all_requests: 0,
            fail_range_requests: 0,
            truncate_range_after: None,
            upload_status: None,
            interrupt_upload: false,
            uploaded: Vec::new(),
            requests: Vec::new(),
        })
    }

    pub fn without_ranges(self) -> Self {
        self.state.lock().unwrap().accepts_ranges = false;
        self
    }

    /// Configure a non-success response such as 403, 404, or 503.
    pub fn with_status(self, status: u16) -> Self {
        self.state.lock().unwrap().response_status = Some(status);
        self
    }

    /// Return a 206 response whose body is valid but whose Content-Range is
    /// intentionally inconsistent with the requested range.
    pub fn with_malformed_content_range(self) -> Self {
        self.state.lock().unwrap().malformed_content_range = true;
        self
    }

    /// Ignore Range and return a 200 full-object response. This models a
    /// server that advertises ranges but does not honor them.
    pub fn ignoring_ranges(self) -> Self {
        self.state.lock().unwrap().ignore_range = true;
        self
    }

    /// Sleep before writing a response. Use a duration larger than the
    /// client's request timeout to model a timeout without a busy loop.
    pub fn delay_response(self, delay: Duration) -> Self {
        self.state.lock().unwrap().response_delay = Some(delay);
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

    pub fn with_upload_status(self, status: u16) -> Self {
        self.state.lock().unwrap().upload_status = Some(status);
        self
    }

    pub fn interrupt_upload(self) -> Self {
        self.state.lock().unwrap().interrupt_upload = true;
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

    pub fn statuses(&self) -> Vec<u16> {
        self.requests()
            .into_iter()
            .map(|request| request.status)
            .collect()
    }

    pub fn range_requests(&self) -> Vec<RequestRecord> {
        self.requests()
            .into_iter()
            .filter(|request| request.range.is_some())
            .collect()
    }

    pub fn uploads(&self) -> Vec<UploadRecord> {
        self.state.lock().unwrap().uploaded.clone()
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

    let Some(header_end) = request
        .windows(4)
        .position(|window| window == b"\r\n\r\n")
        .map(|offset| offset + 4)
    else {
        return;
    };
    let header = &request[..header_end];
    let content_length = String::from_utf8_lossy(header)
        .lines()
        .find_map(|line| {
            line.strip_prefix("Content-Length:")
                .and_then(|value| value.trim().parse::<usize>().ok())
        })
        .unwrap_or(0);
    let mut request_body = request[header_end..].to_vec();
    let mut body_buffer = [0_u8; 8192];
    while request_body.len() < content_length {
        match stream.read(&mut body_buffer) {
            Ok(0) => break,
            Ok(count) => request_body.extend_from_slice(&body_buffer[..count]),
            Err(_) => return,
        }
    }
    request_body.truncate(content_length.min(request_body.len()));
    let request_text = String::from_utf8_lossy(header);
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
    // Keep fixture diagnostics safe to print: query strings may contain
    // signed URL credentials even when the production locator is redacted.
    let recorded_path = path.split_once('?').map_or(path, |(path, _)| path);
    let recorded_path = recorded_path
        .split_once('#')
        .map_or(recorded_path, |(path, _)| path);
    let range = request_text.lines().find_map(parse_range_header);
    let requested_bytes = range
        .and_then(|(start, end)| {
            end.checked_sub(start)
                .and_then(|length| length.checked_add(1))
        })
        .unwrap_or(0);
    let record = RequestRecord {
        method: method.to_string(),
        path: recorded_path.to_string(),
        range,
        status: 0,
        content_range: None,
        requested_bytes,
        returned_bytes: 0,
    };

    let (record_index, response) = {
        let mut fixture = state.lock().unwrap();
        let record_index = fixture.requests.len();
        fixture.requests.push(record);
        let response = if fixture.fail_all_requests > 0 {
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
                fixture.response_status,
                fixture.malformed_content_range,
                fixture.ignore_range,
                fixture.response_delay,
                fixture.truncate_range_after,
                fixture.upload_status,
                fixture.interrupt_upload,
            ))
        };
        (record_index, response)
    };
    let Some((
        body,
        etag,
        accepts_ranges,
        response_status,
        malformed_content_range,
        ignore_range,
        response_delay,
        truncate_range_after,
        upload_status,
        interrupt_upload,
    )) = response
    else {
        // Dropping the connection makes curl report a transport failure and
        // lets ObjectResource exercise its retry loop.
        let mut fixture = state.lock().unwrap();
        if let Some(record) = fixture.requests.get_mut(record_index) {
            record.status = 0;
        }
        return;
    };

    if let Some(delay) = response_delay {
        thread::sleep(delay);
    }

    if let Some(status) = response_status {
        let status_text = status_reason(status);
        let response_body = format!("fixture response {status}\n").into_bytes();
        let headers = format!(
            "HTTP/1.1 {status} {status_text}\r\nContent-Length: {}\r\nETag: \"{etag}\"\r\nConnection: close\r\n\r\n",
            response_body.len()
        );
        let _ = stream.write_all(headers.as_bytes());
        if !method.eq_ignore_ascii_case("HEAD") {
            let _ = stream.write_all(&response_body);
        }
        let _ = stream.flush();
        let mut fixture = state.lock().unwrap();
        if let Some(record) = fixture.requests.get_mut(record_index) {
            record.status = status;
            record.returned_bytes = if method.eq_ignore_ascii_case("HEAD") {
                0
            } else {
                response_body.len()
            };
        }
        return;
    }

    if method.eq_ignore_ascii_case("PUT") || method.eq_ignore_ascii_case("POST") {
        if interrupt_upload {
            let mut fixture = state.lock().unwrap();
            fixture.uploaded.push(UploadRecord {
                method: method.to_string(),
                path: recorded_path.to_string(),
                body: request_body,
            });
            if let Some(record) = fixture.requests.get_mut(record_index) {
                record.status = 0;
            }
            return;
        }
        let status = upload_status.unwrap_or(200);
        let status_text = status_reason(status);
        let headers = format!(
            "HTTP/1.1 {status} {status_text}\r\nContent-Length: 0\r\nETag: \"{etag}\"\r\nConnection: close\r\n\r\n"
        );
        let _ = stream.write_all(headers.as_bytes());
        let _ = stream.flush();
        let mut fixture = state.lock().unwrap();
        fixture.uploaded.push(UploadRecord {
            method: method.to_string(),
            path: recorded_path.to_string(),
            body: request_body,
        });
        if let Some(record) = fixture.requests.get_mut(record_index) {
            record.status = status;
        }
        return;
    }

    let (status, selected, content_range, total_length) = if let (true, Some((start, end))) = (
        method.eq_ignore_ascii_case("GET") && accepts_ranges && !ignore_range,
        range,
    ) {
        if start >= body.len() || end < start || end >= body.len() {
            let _ = write!(
                stream,
                "HTTP/1.1 416 Range Not Satisfiable\r\nContent-Length: 0\r\nConnection: close\r\n\r\n"
            );
            let mut fixture = state.lock().unwrap();
            if let Some(record) = fixture.requests.get_mut(record_index) {
                record.status = 416;
            }
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
    let response_content_range = content_range.map(|(start, end)| {
        if malformed_content_range {
            if start < end {
                (start.saturating_add(1), end, total_length)
            } else {
                (start, end.saturating_add(1), total_length)
            }
        } else {
            (start, end, total_length)
        }
    });
    let mut headers =
        format!("HTTP/1.1 {status}\r\nContent-Length: {body_length}\r\nETag: \"{etag}\"\r\n");
    if accepts_ranges {
        headers.push_str("Accept-Ranges: bytes\r\n");
    }
    if let Some((start, end, total_length)) = response_content_range {
        headers.push_str(&format!(
            "Content-Range: bytes {start}-{end}/{total_length}\r\n"
        ));
    }
    headers.push_str("Connection: close\r\n\r\n");
    if stream.write_all(headers.as_bytes()).is_err() {
        let mut fixture = state.lock().unwrap();
        if let Some(record) = fixture.requests.get_mut(record_index) {
            record.status = status_code(status);
            record.content_range = response_content_range;
        }
        return;
    }
    if method.eq_ignore_ascii_case("HEAD") {
        let _ = stream.flush();
        let mut fixture = state.lock().unwrap();
        if let Some(record) = fixture.requests.get_mut(record_index) {
            record.status = status_code(status);
            record.content_range = response_content_range;
            record.returned_bytes = 0;
        }
        return;
    }
    let bytes_to_send = truncate_range_after
        .filter(|_| content_range.is_some())
        .map_or(selected.len(), |limit| limit.min(selected.len()));
    let _ = stream.write_all(&selected[..bytes_to_send]);
    let _ = stream.flush();

    let mut fixture = state.lock().unwrap();
    if let Some(record) = fixture.requests.get_mut(record_index) {
        record.status = status_code(status);
        record.content_range = response_content_range;
        record.returned_bytes = if method.eq_ignore_ascii_case("HEAD") {
            0
        } else {
            bytes_to_send
        };
    }
}

fn status_code(status: &str) -> u16 {
    status
        .split_whitespace()
        .next()
        .and_then(|value| value.parse().ok())
        .unwrap_or(0)
}

fn status_reason(status: u16) -> &'static str {
    match status {
        200 => "OK",
        201 => "Created",
        204 => "No Content",
        206 => "Partial Content",
        403 => "Forbidden",
        404 => "Not Found",
        408 => "Request Timeout",
        416 => "Range Not Satisfiable",
        500 => "Internal Server Error",
        502 => "Bad Gateway",
        503 => "Service Unavailable",
        504 => "Gateway Timeout",
        _ => "Fixture Response",
    }
}

fn parse_range_header(line: &str) -> Option<(usize, usize)> {
    let value = line.strip_prefix("Range: bytes=")?.trim();
    let (start, end) = value.split_once('-')?;
    Some((start.parse().ok()?, end.parse().ok()?))
}
