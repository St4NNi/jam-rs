//! One TCP connection that reads and writes frames. Reads and writes end at a deadline, reads also
//! end soon after a stop flag is set, and a write never waits longer than `WRITE_STALL` on a
//! reader that does not read.

use super::ServiceError;
use super::frame::{Frame, read_frame, write_frame};
use std::io::{BufReader, BufWriter, ErrorKind, Read, Write};
use std::net::{SocketAddr, TcpStream, ToSocketAddrs};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::{Duration, Instant};

/// How often a waiting read looks at the stop flag.
pub const POLL: Duration = Duration::from_millis(50);

/// The longest a single write may wait for the other side to read.
pub const WRITE_STALL: Duration = Duration::from_secs(60);

/// A TCP connection that carries frames.
pub struct Connection {
    reader: BufReader<TimedStream>,
    writer: BufWriter<TimedStream>,
}

/// One direction of a connection with its deadline and stop flag.
struct TimedStream {
    stream: TcpStream,
    deadline: Option<Instant>,
    stop: Option<Arc<AtomicBool>>,
}

impl Connection {
    /// Wraps an accepted or opened stream.
    pub fn new(stream: TcpStream) -> Result<Self, ServiceError> {
        stream.set_nodelay(true)?;
        let timed = |stream| TimedStream {
            stream,
            deadline: None,
            stop: None,
        };
        Ok(Self {
            reader: BufReader::new(timed(stream.try_clone()?)),
            writer: BufWriter::new(timed(stream)),
        })
    }

    /// Opens a connection to `address` (host and port), trying each address it names until one
    /// answers within `timeout`.
    pub fn open(address: &str, timeout: Duration) -> Result<Self, ServiceError> {
        let addresses: Vec<SocketAddr> = address.to_socket_addrs()?.collect();
        let mut failure = std::io::Error::new(ErrorKind::NotFound, "the address names no host");
        for address in addresses {
            match TcpStream::connect_timeout(&address, timeout) {
                Ok(stream) => return Self::new(stream),
                Err(error) => failure = error,
            }
        }
        Err(failure.into())
    }

    /// Sets when reads stop waiting.
    pub fn set_read_deadline(&mut self, deadline: Option<Instant>) {
        self.reader.get_mut().deadline = deadline;
    }

    /// Sets when writes stop waiting.
    pub fn set_write_deadline(&mut self, deadline: Option<Instant>) {
        self.writer.get_mut().deadline = deadline;
    }

    /// Sets the flag that ends a waiting read.
    pub fn set_stop(&mut self, stop: Option<Arc<AtomicBool>>) {
        self.reader.get_mut().stop = stop;
    }

    /// The underlying stream, such as for watching it for a disconnect.
    pub fn stream(&self) -> &TcpStream {
        &self.writer.get_ref().stream
    }

    /// Reads the next frame, or `None` when the other side closed between frames.
    pub fn read(&mut self) -> Result<Option<Frame>, ServiceError> {
        read_frame(&mut self.reader)
    }

    /// Writes and flushes one frame.
    pub fn write(&mut self, frame: &Frame) -> Result<(), ServiceError> {
        write_frame(&mut self.writer, frame)
    }

    /// Waits for the first byte of the next frame until the read deadline. False when the other
    /// side closed the connection or the deadline passed first.
    pub fn wait_for_frame(&mut self) -> Result<bool, ServiceError> {
        match std::io::BufRead::fill_buf(&mut self.reader) {
            Ok(bytes) => Ok(!bytes.is_empty()),
            Err(error) if error.kind() == ErrorKind::TimedOut => Ok(false),
            Err(error) => Err(error.into()),
        }
    }
}

/// The time left until `deadline`, or an error once it passed.
fn remaining(deadline: Option<Instant>) -> std::io::Result<Option<Duration>> {
    let Some(deadline) = deadline else {
        return Ok(None);
    };
    match deadline.checked_duration_since(Instant::now()) {
        Some(left) if !left.is_zero() => Ok(Some(left)),
        _ => Err(std::io::Error::new(
            ErrorKind::TimedOut,
            "the deadline passed",
        )),
    }
}

fn is_timeout(error: &std::io::Error) -> bool {
    matches!(error.kind(), ErrorKind::WouldBlock | ErrorKind::TimedOut)
}

impl Read for TimedStream {
    fn read(&mut self, bytes: &mut [u8]) -> std::io::Result<usize> {
        loop {
            if self
                .stop
                .as_ref()
                .is_some_and(|stop| stop.load(Ordering::Relaxed))
            {
                return Err(std::io::Error::new(
                    ErrorKind::ConnectionAborted,
                    "the search was asked to stop",
                ));
            }
            let left = remaining(self.deadline)?;
            // Without a stop flag a read waits until the deadline at once.
            let wait = match (&self.stop, left) {
                (Some(_), left) => Some(left.map_or(POLL, |left| left.min(POLL))),
                (None, left) => left,
            };
            self.stream.set_read_timeout(wait)?;
            match self.stream.read(bytes) {
                Err(error) if is_timeout(&error) => {}
                read => return read,
            }
        }
    }
}

impl Write for TimedStream {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let left = remaining(self.deadline)?.map_or(WRITE_STALL, |left| left.min(WRITE_STALL));
        self.stream.set_write_timeout(Some(left))?;
        self.stream
            .write(bytes)
            .map_err(|error| match is_timeout(&error) {
                true => std::io::Error::new(ErrorKind::TimedOut, "the other side does not read"),
                false => error,
            })
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.stream.flush()
    }
}

/// Sets `stop` when the other side of `stream` closes the connection, and returns when `done` is
/// set or the other side sends a byte. A server runs it while a search runs, as the client then
/// sends nothing.
pub fn watch_disconnect(stream: &TcpStream, done: &AtomicBool, stop: &AtomicBool) {
    if stream.set_read_timeout(Some(POLL)).is_err() {
        return;
    }
    while !done.load(Ordering::Relaxed) {
        match stream.peek(&mut [0]) {
            Ok(0) => return stop.store(true, Ordering::Relaxed),
            Ok(_) => return,
            Err(error) if is_timeout(&error) || error.kind() == ErrorKind::Interrupted => {}
            Err(_) => return stop.store(true, Ordering::Relaxed),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::net::TcpListener;

    /// A connected pair: the client side and the accepted side.
    fn pair() -> (Connection, TcpStream) {
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let address = listener.local_addr().unwrap().to_string();
        let client = Connection::open(&address, Duration::from_secs(5)).unwrap();
        (client, listener.accept().unwrap().0)
    }

    #[test]
    fn frames_cross_a_connection_and_reads_end_at_the_deadline() {
        let (mut client, accepted) = pair();
        let mut server = Connection::new(accepted).unwrap();
        server.write(&Frame::Complete { complete: true }).unwrap();
        assert_eq!(
            client.read().unwrap(),
            Some(Frame::Complete { complete: true })
        );

        let started = Instant::now();
        client.set_read_deadline(Some(started + Duration::from_millis(100)));
        let error = client.read().unwrap_err();
        assert!(matches!(error, ServiceError::Io(e) if e.kind() == ErrorKind::TimedOut));
        assert!(started.elapsed() < Duration::from_secs(2));
        drop(server);
        client.set_read_deadline(None);
        assert_eq!(client.read().unwrap(), None);
    }

    #[test]
    fn a_stop_flag_ends_a_waiting_read_and_a_close_is_seen_by_the_watcher() {
        let (mut client, accepted) = pair();
        let stop = Arc::new(AtomicBool::new(false));
        client.set_stop(Some(stop.clone()));
        std::thread::scope(|scope| {
            scope.spawn(|| {
                std::thread::sleep(Duration::from_millis(100));
                stop.store(true, Ordering::Relaxed);
            });
            assert!(client.read().is_err());
        });

        let (done, closed) = (AtomicBool::new(false), AtomicBool::new(false));
        std::thread::scope(|scope| {
            scope.spawn(|| watch_disconnect(&accepted, &done, &closed));
            drop(client);
        });
        assert!(closed.load(Ordering::Relaxed));
    }
}
