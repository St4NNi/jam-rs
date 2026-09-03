use crate::jidx::sha256_reader;
use s3::Bucket;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::Path;
use thiserror::Error;

// ponytail: one 1 MiB block; add a bounded multi-block cache only if profiles show reuse.
const S3_BLOCK_BYTES: u64 = 1024 * 1024;

pub enum RangeSource {
    Local(LocalSource),
    S3(S3Source),
}

pub struct LocalSource {
    file: File,
    length: u64,
    stats: RangeStats,
}

pub struct S3Source {
    bucket: Box<Bucket>,
    key: String,
    length: u64,
    etag: Option<String>,
    position: u64,
    cache_offset: u64,
    cache: Vec<u8>,
    stats: RangeStats,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct RangeStats {
    pub metadata_requests: u64,
    pub read_requests: u64,
    pub bytes_read: u64,
}

impl RangeSource {
    pub fn local(path: impl AsRef<Path>, expected_size: u64) -> Result<Self, RangeSourceError> {
        let file = File::open(path)?;
        let length = file.metadata()?.len();
        if length != expected_size {
            return Err(RangeSourceError::SizeMismatch);
        }
        Ok(Self::Local(LocalSource {
            file,
            length,
            stats: RangeStats {
                metadata_requests: 1,
                ..RangeStats::default()
            },
        }))
    }

    pub fn s3(
        bucket: Box<Bucket>,
        key: impl Into<String>,
        expected_size: u64,
    ) -> Result<Self, RangeSourceError> {
        let key = key.into();
        if key.is_empty() {
            return Err(RangeSourceError::InvalidKey);
        }
        let (head, status) = bucket.head_object(&key).map_err(|_| RangeSourceError::S3)?;
        let length = head
            .content_length
            .and_then(|length| u64::try_from(length).ok())
            .ok_or(RangeSourceError::S3)?;
        if status != 200 || length != expected_size {
            return Err(RangeSourceError::SizeMismatch);
        }
        Ok(Self::S3(S3Source {
            bucket,
            key,
            length,
            etag: head.e_tag,
            position: 0,
            cache_offset: 0,
            cache: Vec::new(),
            stats: RangeStats {
                metadata_requests: 1,
                ..RangeStats::default()
            },
        }))
    }

    pub fn len(&self) -> u64 {
        match self {
            Self::Local(source) => source.length,
            Self::S3(source) => source.length,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn stats(&self) -> RangeStats {
        match self {
            Self::Local(source) => source.stats,
            Self::S3(source) => source.stats,
        }
    }

    pub fn verify_sha256(&mut self, expected: [u8; 32]) -> io::Result<bool> {
        let position = self.stream_position()?;
        self.seek(SeekFrom::Start(0))?;
        let actual = sha256_reader(&mut *self);
        let restored = self.seek(SeekFrom::Start(position));
        let actual = actual?;
        restored?;
        Ok(actual == expected)
    }
}

impl Read for RangeSource {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        match self {
            Self::Local(source) => {
                let position = source.file.stream_position()?;
                if position >= source.length {
                    return Ok(0);
                }
                let remaining = usize::try_from(source.length - position).unwrap_or(usize::MAX);
                let count = buffer.len().min(remaining);
                let read = source.file.read(&mut buffer[..count])?;
                if read != 0 {
                    source.stats.read_requests = source.stats.read_requests.saturating_add(1);
                    source.stats.bytes_read = source.stats.bytes_read.saturating_add(read as u64);
                }
                Ok(read)
            }
            Self::S3(source) => source.read(buffer),
        }
    }
}

impl Seek for RangeSource {
    fn seek(&mut self, position: SeekFrom) -> io::Result<u64> {
        match self {
            Self::Local(source) => {
                let current = source.file.stream_position()?;
                let target = checked_position(source.length, current, position)?;
                source.file.seek(SeekFrom::Start(target))
            }
            Self::S3(source) => {
                source.position = checked_position(source.length, source.position, position)?;
                Ok(source.position)
            }
        }
    }
}

impl S3Source {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        let mut written = 0;
        while written < buffer.len() && self.position < self.length {
            let cache_end = self
                .cache_offset
                .checked_add(self.cache.len() as u64)
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "S3 cache range"))?;
            if self.cache.is_empty()
                || self.position < self.cache_offset
                || self.position >= cache_end
            {
                self.fill_cache()?;
                continue;
            }
            let cache_start = usize::try_from(self.position - self.cache_offset)
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "S3 cache offset"))?;
            let count = (buffer.len() - written).min(self.cache.len() - cache_start);
            buffer[written..written + count]
                .copy_from_slice(&self.cache[cache_start..cache_start + count]);
            written += count;
            self.position += count as u64;
        }
        Ok(written)
    }

    fn fill_cache(&mut self) -> io::Result<()> {
        let end = self
            .position
            .saturating_add(S3_BLOCK_BYTES)
            .min(self.length)
            - 1;
        let response = self
            .bucket
            .get_object_range(&self.key, self.position, Some(end))
            .map_err(|_| io::Error::other("S3 range request failed"))?;
        let expected = usize::try_from(end - self.position + 1)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "S3 range length"))?;
        let full_object = self.position == 0 && end.checked_add(1) == Some(self.length);
        let status_valid =
            response.status_code() == 206 || (response.status_code() == 200 && full_object);
        let headers = response.headers();
        let header = |name: &str| {
            headers
                .iter()
                .find(|(key, _)| key.eq_ignore_ascii_case(name))
                .map(|(_, value)| value.as_str())
        };
        let expected_range = format!("bytes {}-{}/{}", self.position, end, self.length);
        let range_valid = response.status_code() != 206
            || header("content-range") == Some(expected_range.as_str());
        let etag_valid = self
            .etag
            .as_deref()
            .is_none_or(|etag| header("etag") == Some(etag));
        if !status_valid || !range_valid || !etag_valid || response.as_slice().len() != expected {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "S3 returned an invalid range",
            ));
        }
        self.cache_offset = self.position;
        self.cache.clear();
        self.cache.extend_from_slice(response.as_slice());
        self.stats.read_requests = self.stats.read_requests.saturating_add(1);
        self.stats.bytes_read = self.stats.bytes_read.saturating_add(expected as u64);
        Ok(())
    }
}

fn checked_position(length: u64, current: u64, position: SeekFrom) -> io::Result<u64> {
    let target = match position {
        SeekFrom::Start(position) => i128::from(position),
        SeekFrom::End(delta) => i128::from(length) + i128::from(delta),
        SeekFrom::Current(delta) => i128::from(current) + i128::from(delta),
    };
    if !(0..=i128::from(length)).contains(&target) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "range seek is outside the source",
        ));
    }
    u64::try_from(target)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "range seek overflow"))
}

#[derive(Debug, Error)]
pub enum RangeSourceError {
    #[error("range source I/O failed: {0}")]
    Io(#[from] io::Error),
    #[error("S3 request failed")]
    S3,
    #[error("range source size differs from JIDX metadata")]
    SizeMismatch,
    #[error("S3 object key is empty")]
    InvalidKey,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::jidx::sha256;
    use s3::{Region, creds::Credentials};
    use std::io::{BufRead, BufReader, Write};
    use std::net::TcpListener;
    use std::thread;

    fn serve(listener: TcpListener, data: Vec<u8>) {
        for _ in 0..2 {
            let (mut stream, _) = listener.accept().unwrap();
            let mut request = String::new();
            let mut reader = BufReader::new(stream.try_clone().unwrap());
            loop {
                let mut line = String::new();
                reader.read_line(&mut line).unwrap();
                if line == "\r\n" || line.is_empty() {
                    break;
                }
                request.push_str(&line);
            }
            if request.starts_with("HEAD ") {
                write!(
                    stream,
                    "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nETag: \"fixture\"\r\nConnection: close\r\n\r\n",
                    data.len()
                )
                .unwrap();
            } else {
                let range = request
                    .lines()
                    .find(|line| line.to_ascii_lowercase().starts_with("range:"))
                    .unwrap()
                    .split_once('=')
                    .unwrap()
                    .1;
                let (start, end) = range.split_once('-').unwrap();
                let start: usize = start.parse().unwrap();
                let end: usize = end.parse().unwrap();
                let body = &data[start..=end];
                write!(
                    stream,
                    "HTTP/1.1 206 Partial Content\r\nContent-Length: {}\r\nContent-Range: bytes {}-{}/{}\r\nETag: \"fixture\"\r\nConnection: close\r\n\r\n",
                    body.len(), start, end, data.len()
                )
                .unwrap();
                stream.write_all(body).unwrap();
            }
        }
    }

    #[test]
    fn s3_reads_checked_ranges_and_reuses_cache() {
        let data: Vec<_> = (0..100).collect();
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let endpoint = format!("http://{}", listener.local_addr().unwrap());
        let server_data = data.clone();
        let server = thread::spawn(move || serve(listener, server_data));
        let credentials =
            Credentials::new(Some("access"), Some("secret"), None, None, None).unwrap();
        let bucket = Bucket::new(
            "bucket",
            Region::Custom {
                region: "test".into(),
                endpoint,
            },
            credentials,
        )
        .unwrap()
        .with_path_style();
        let mut source = RangeSource::s3(bucket, "object", data.len() as u64).unwrap();
        let mut first = [0; 5];
        source.read_exact(&mut first).unwrap();
        assert_eq!(first, data[..5]);
        source.seek(SeekFrom::Start(90)).unwrap();
        let mut last = [0; 10];
        source.read_exact(&mut last).unwrap();
        assert_eq!(last, data[90..]);
        assert!(source.verify_sha256(sha256(&data)).unwrap());
        assert_eq!(source.stats().read_requests, 1);
        server.join().unwrap();
    }
}
