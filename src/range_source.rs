use crate::jidx::sha256_reader;
use s3::{Bucket, Region, creds::Credentials};
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::Path;
use thiserror::Error;

const S3_BLOCK_BYTES: u64 = 1024 * 1024;

/// Credential selection for the limited-credentials development build.
pub fn development_credentials() -> Result<Credentials, RangeSourceError> {
    let present = |name| std::env::var_os(name).is_some();
    if [
        "AWS_ROLE_ARN",
        "AWS_WEB_IDENTITY_TOKEN_FILE",
        "AWS_CONTAINER_CREDENTIALS_RELATIVE_URI",
        "AWS_CONTAINER_CREDENTIALS_FULL_URI",
        "AWS_EC2_METADATA_SERVICE_ENDPOINT",
    ]
    .into_iter()
    .any(present)
    {
        return Err(RangeSourceError::UnsupportedCredentialProvider);
    }
    let credentials = if [
        "AWS_ACCESS_KEY_ID",
        "AWS_SECRET_ACCESS_KEY",
        "AWS_SESSION_TOKEN",
        "AWS_SECURITY_TOKEN",
    ]
    .into_iter()
    .any(present)
    {
        // Partial environment credentials must not select another account from a profile.
        Credentials::from_env().map_err(|_| RangeSourceError::UnsupportedCredentialProvider)?
    } else {
        let profile = std::env::var_os("AWS_PROFILE")
            .or_else(|| std::env::var_os("AWS_DEFAULT_PROFILE"))
            .map(|value| value.into_string())
            .transpose()
            .map_err(|_| RangeSourceError::UnsupportedCredentialProvider)?;
        if present("AWS_SHARED_CREDENTIALS_FILE")
            && std::env::var("AWS_SHARED_CREDENTIALS_FILE").is_err()
        {
            return Err(RangeSourceError::UnsupportedCredentialProvider);
        }
        static_profile_credentials(profile.as_deref().unwrap_or("default"))?
    };
    if credentials.access_key.as_deref().is_none_or(str::is_empty)
        || credentials.secret_key.as_deref().is_none_or(str::is_empty)
    {
        return Err(RangeSourceError::UnsupportedCredentialProvider);
    }
    Ok(credentials)
}

fn static_profile_credentials(profile: &str) -> Result<Credentials, RangeSourceError> {
    let default_path = |name| {
        home::home_dir()
            .map(|path| path.join(".aws").join(name))
            .ok_or(RangeSourceError::UnsupportedCredentialProvider)
    };
    let credentials_path = match std::env::var_os("AWS_SHARED_CREDENTIALS_FILE") {
        Some(path) => path.into(),
        None => default_path("credentials")?,
    };
    let selected_config = std::env::var_os("AWS_CONFIG_FILE");
    let config_path = match selected_config.as_ref() {
        Some(path) => std::path::PathBuf::from(path),
        None => default_path("config")?,
    };
    let credentials = ini::Ini::load_from_file(credentials_path)
        .map_err(|_| RangeSourceError::UnsupportedCredentialProvider)?;
    let config = match ini::Ini::load_from_file(config_path) {
        Ok(config) => Some(config),
        Err(ini::Error::Io(error))
            if selected_config.is_none() && error.kind() == io::ErrorKind::NotFound =>
        {
            None
        }
        Err(_) => return Err(RangeSourceError::UnsupportedCredentialProvider),
    };
    let config_section = if profile == "default" {
        profile.to_owned()
    } else {
        format!("profile {profile}")
    };
    let sections = credentials.section_all(Some(profile)).chain(
        config
            .iter()
            .flat_map(|config| config.section_all(Some(&config_section))),
    );
    for section in sections {
        if section.iter().any(|(key, _)| {
            matches!(
                key,
                "role_arn"
                    | "source_profile"
                    | "credential_source"
                    | "web_identity_token_file"
                    | "credential_process"
            ) || key.starts_with("sso_")
        }) {
            return Err(RangeSourceError::UnsupportedCredentialProvider);
        }
    }
    let section = credentials
        .section(Some(profile))
        .ok_or(RangeSourceError::UnsupportedCredentialProvider)?;
    let access = section
        .get("aws_access_key_id")
        .ok_or(RangeSourceError::UnsupportedCredentialProvider)?;
    let secret = section
        .get("aws_secret_access_key")
        .ok_or(RangeSourceError::UnsupportedCredentialProvider)?;
    Credentials::new(
        Some(access),
        Some(secret),
        section.get("aws_security_token"),
        section.get("aws_session_token"),
        None,
    )
    .map_err(|_| RangeSourceError::UnsupportedCredentialProvider)
}

pub enum RangeSource {
    Local(LocalSource),
    S3(S3Source),
}

#[derive(Clone)]
pub struct S3Config {
    region: Region,
    credentials: Credentials,
    path_style: bool,
}

impl S3Config {
    pub fn new(
        region: &str,
        endpoint: Option<&str>,
        path_style: bool,
        credentials: Credentials,
    ) -> Result<Self, RangeSourceError> {
        let region = match endpoint {
            Some(endpoint) if !endpoint.is_empty() => Region::Custom {
                region: region.to_string(),
                endpoint: endpoint.to_string(),
            },
            Some(_) => return Err(RangeSourceError::S3),
            None => region.parse().map_err(|_| RangeSourceError::S3)?,
        };
        Ok(Self {
            region,
            credentials,
            path_style,
        })
    }

    fn bucket(&self, name: &str) -> Result<Box<Bucket>, RangeSourceError> {
        let bucket = Bucket::new(name, self.region.clone(), self.credentials.clone())
            .map_err(|_| RangeSourceError::S3)?;
        Ok(if self.path_style {
            bucket.with_path_style()
        } else {
            bucket
        })
    }
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
    etag: String,
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
    pub read_nanoseconds: Option<u64>,
}

impl RangeSource {
    pub fn open(
        uri: &str,
        expected_size: u64,
        s3: Option<&S3Config>,
    ) -> Result<Self, RangeSourceError> {
        if let Some(value) = uri.strip_prefix("s3://") {
            let (bucket, key) = value.split_once('/').ok_or(RangeSourceError::InvalidKey)?;
            if bucket.is_empty() || key.is_empty() {
                return Err(RangeSourceError::InvalidKey);
            }
            return Self::s3(
                s3.ok_or(RangeSourceError::MissingS3Config)?
                    .bucket(bucket)?,
                key,
                expected_size,
            );
        }
        let path = if let Some(path) = uri.strip_prefix("file://") {
            if !path.starts_with('/') {
                return Err(RangeSourceError::InvalidLocalUri);
            }
            Path::new(path)
        } else {
            if uri.contains("://") {
                return Err(RangeSourceError::InvalidLocalUri);
            }
            Path::new(uri)
        };
        Self::local(path, expected_size)
    }

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
        let etag = head.e_tag.ok_or(RangeSourceError::S3)?;
        if etag.trim().is_empty() || etag.trim_start().starts_with("W/") {
            return Err(RangeSourceError::S3);
        }
        Ok(Self::S3(S3Source {
            bucket,
            key,
            length,
            etag,
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

    pub(crate) fn enable_timing(&mut self) {
        if let Self::Local(source) = self {
            source.stats.read_nanoseconds = Some(0);
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
                let started = source
                    .stats
                    .read_nanoseconds
                    .map(|_| std::time::Instant::now());
                let position = source.file.stream_position()?;
                if position >= source.length {
                    return Ok(0);
                }
                let remaining = usize::try_from(source.length - position).unwrap_or(usize::MAX);
                let count = buffer.len().min(remaining);
                let read = source.file.read(&mut buffer[..count])?;
                if let (Some(started), Some(elapsed)) =
                    (started, &mut source.stats.read_nanoseconds)
                {
                    *elapsed = elapsed.saturating_add(started.elapsed().as_nanos() as u64);
                }
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
        let etag_valid = header("etag") == Some(self.etag.as_str());
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
    #[error("S3 configuration is required")]
    MissingS3Config,
    #[error(
        "unsupported credential provider in this limited-credentials development build: use complete explicit AWS environment or static profile credentials; STS web identity, container endpoints and instance metadata are unavailable"
    )]
    UnsupportedCredentialProvider,
    #[error("range source is not a supported local URI")]
    InvalidLocalUri,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::jidx::sha256;
    use s3::creds::Credentials;
    use std::io::{BufRead, BufReader, Write};
    use std::net::{TcpListener, TcpStream};
    use std::thread;

    fn read_request(stream: &TcpStream) -> String {
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
        request
    }

    fn serve(listener: TcpListener, data: Vec<u8>) {
        for _ in 0..2 {
            let (mut stream, _) = listener.accept().unwrap();
            let request = read_request(&stream);
            let headers = request.to_ascii_lowercase();
            assert!(headers.contains("authorization: aws4-hmac-sha256 "));
            assert!(headers.contains("x-amz-security-token: fixture-session"));
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

    fn serve_head(listener: TcpListener, length: usize, etag: Option<&str>) {
        let (mut stream, _) = listener.accept().unwrap();
        read_request(&stream);
        let etag = etag
            .map(|etag| format!("ETag: {etag}\r\n"))
            .unwrap_or_default();
        write!(
            stream,
            "HTTP/1.1 200 OK\r\nContent-Length: {length}\r\n{etag}Connection: close\r\n\r\n"
        )
        .unwrap();
    }

    #[test]
    fn s3_reads_checked_ranges_and_reuses_cache() {
        let data: Vec<_> = (0..100).collect();
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let endpoint = format!("http://{}", listener.local_addr().unwrap());
        let server_data = data.clone();
        let server = thread::spawn(move || serve(listener, server_data));
        let credentials = Credentials::new(
            Some("access"),
            Some("secret"),
            None,
            Some("fixture-session"),
            None,
        )
        .unwrap();
        let config = S3Config::new("test", Some(&endpoint), true, credentials).unwrap();
        let mut source =
            RangeSource::open("s3://bucket/object", data.len() as u64, Some(&config)).unwrap();
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

    #[test]
    fn s3_rejects_missing_empty_and_weak_validators() {
        for (etag, valid) in [
            (None, false),
            (Some(""), false),
            (Some("W/\"fixture\""), false),
            (Some("\"\""), true),
        ] {
            let listener = TcpListener::bind("127.0.0.1:0").unwrap();
            let endpoint = format!("http://{}", listener.local_addr().unwrap());
            let server = thread::spawn(move || serve_head(listener, 100, etag));
            let credentials =
                Credentials::new(Some("access"), Some("secret"), None, None, None).unwrap();
            let config = S3Config::new("test", Some(&endpoint), true, credentials).unwrap();
            assert_eq!(
                RangeSource::open("s3://bucket/object", 100, Some(&config)).is_ok(),
                valid,
            );
            server.join().unwrap();
        }
    }

    #[test]
    fn limited_credentials_child() {
        let Ok(case) = std::env::var("JAM_CREDENTIAL_TEST") else {
            return;
        };
        let result = development_credentials();
        if case == "supported" {
            let credentials = result.unwrap();
            assert!(credentials.access_key.as_deref() == Some("fixture-access"));
            assert!(credentials.session_token.as_deref() == Some("fixture-session"));
        } else {
            assert!(matches!(
                result,
                Err(RangeSourceError::UnsupportedCredentialProvider)
            ));
        }
    }

    #[test]
    fn limited_credentials_are_explicit_and_never_contact_providers() {
        let directory = tempfile::tempdir().unwrap();
        let profile = directory.path().join("credentials");
        let config = directory.path().join("config");
        std::fs::write(&config, "").unwrap();
        std::fs::write(&profile, "[default]\naws_access_key_id=fixture-access\naws_secret_access_key=fixture-secret\naws_session_token=fixture-session\n[named]\naws_access_key_id=fixture-access\naws_secret_access_key=fixture-secret\naws_session_token=fixture-session\n[role]\nrole_arn=fixture-role\nsource_profile=default\n").unwrap();
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        listener.set_nonblocking(true).unwrap();
        let endpoint = format!("http://{}", listener.local_addr().unwrap());
        for case in [
            "environment",
            "profile",
            "named",
            "partial",
            "missing",
            "role",
            "AWS_WEB_IDENTITY_TOKEN_FILE",
            "AWS_CONTAINER_CREDENTIALS_RELATIVE_URI",
            "AWS_CONTAINER_CREDENTIALS_FULL_URI",
            "AWS_EC2_METADATA_SERVICE_ENDPOINT",
        ] {
            let supported = matches!(case, "environment" | "profile" | "named");
            let mut child = std::process::Command::new(std::env::current_exe().unwrap());
            child
                .args(["--exact", "range_source::tests::limited_credentials_child"])
                .env_clear()
                .env("AWS_SHARED_CREDENTIALS_FILE", &profile)
                .env("AWS_CONFIG_FILE", &config)
                .env(
                    "JAM_CREDENTIAL_TEST",
                    if supported {
                        "supported"
                    } else {
                        "unsupported"
                    },
                );
            match case {
                "environment" => {
                    child
                        .env("AWS_ACCESS_KEY_ID", "fixture-access")
                        .env("AWS_SECRET_ACCESS_KEY", "fixture-secret")
                        .env("AWS_SESSION_TOKEN", "fixture-session");
                }
                "named" => {
                    child.env("AWS_PROFILE", "named");
                }
                "partial" => {
                    child.env("AWS_ACCESS_KEY_ID", "fixture-access");
                }
                "missing" | "role" => {
                    child.env("AWS_PROFILE", case);
                }
                "profile" => {}
                provider => {
                    child.env(provider, &endpoint);
                }
            }
            assert!(
                child.output().unwrap().status.success(),
                "credential case {case}"
            );
            assert!(
                matches!(listener.accept(), Err(error) if error.kind() == io::ErrorKind::WouldBlock)
            );
        }
    }

    #[test]
    fn limited_credentials_reject_provider_profiles_with_static_keys() {
        let directory = tempfile::tempdir().unwrap();
        let credentials = directory.path().join("credentials");
        let config = directory.path().join("config");
        let keys = "aws_access_key_id=fixture-access\naws_secret_access_key=fixture-secret\naws_session_token=fixture-session\n";
        for selected in ["default", "named"] {
            for in_config in [false, true] {
                for selected_provider in [false, true] {
                    for provider in [
                        "role_arn",
                        "source_profile",
                        "credential_source",
                        "web_identity_token_file",
                        "credential_process",
                        "sso_start_url",
                        "sso_session",
                    ] {
                        let provider_profile = if selected_provider { selected } else { "other" };
                        let mut shared = format!("[{selected}]\n{keys}");
                        let mut settings = String::new();
                        if in_config {
                            let section = if provider_profile == "default" {
                                "default".into()
                            } else {
                                format!("profile {provider_profile}")
                            };
                            settings = format!("[{section}]\n{provider}=fixture-provider\n");
                        } else {
                            if !selected_provider {
                                shared.push_str("[other]\n");
                            }
                            shared.push_str(&format!("{provider}=fixture-provider\n"));
                        }
                        std::fs::write(&credentials, shared).unwrap();
                        std::fs::write(&config, settings).unwrap();
                        let mut child =
                            std::process::Command::new(std::env::current_exe().unwrap());
                        child
                            .args(["--exact", "range_source::tests::limited_credentials_child"])
                            .env_clear()
                            .env("AWS_SHARED_CREDENTIALS_FILE", &credentials)
                            .env("AWS_CONFIG_FILE", &config)
                            .env("AWS_PROFILE", selected)
                            .env(
                                "JAM_CREDENTIAL_TEST",
                                if selected_provider {
                                    "unsupported"
                                } else {
                                    "supported"
                                },
                            );
                        assert!(
                            child.output().unwrap().status.success(),
                            "provider {provider}, config {in_config}, selected {selected_provider}"
                        );
                        // Explicit complete environment credentials keep their documented precedence.
                        child
                            .env("AWS_ACCESS_KEY_ID", "fixture-access")
                            .env("AWS_SECRET_ACCESS_KEY", "fixture-secret")
                            .env("AWS_SESSION_TOKEN", "fixture-session")
                            .env("JAM_CREDENTIAL_TEST", "supported");
                        assert!(child.output().unwrap().status.success());
                    }
                }
            }
        }
    }

    #[test]
    fn s3_xml_shaped_ranges_are_bytes_and_error_bodies_are_redacted() {
        let data = b"<broken xmlns:a='x' xmlns:a='y' attribute='one' attribute='two'".to_vec();
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let endpoint = format!("http://{}", listener.local_addr().unwrap());
        let server_data = data.clone();
        let server = thread::spawn(move || serve(listener, server_data));
        let credentials = Credentials::new(
            Some("access"),
            Some("secret"),
            None,
            Some("fixture-session"),
            None,
        )
        .unwrap();
        let config = S3Config::new("test", Some(&endpoint), true, credentials.clone()).unwrap();
        let mut source =
            RangeSource::open("s3://bucket/object", data.len() as u64, Some(&config)).unwrap();
        let mut actual = Vec::new();
        source.read_to_end(&mut actual).unwrap();
        assert_eq!(actual, data);
        server.join().unwrap();
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let endpoint = format!("http://{}", listener.local_addr().unwrap());
        let server = thread::spawn(move || {
            let (mut stream, _) = listener.accept().unwrap();
            read_request(&stream);
            write!(stream, "HTTP/1.1 200 OK\r\nContent-Length: 34\r\nETag: \"fixture\"\r\nConnection: close\r\n\r\n").unwrap();
            drop(stream);
            for _ in 0..=s3::get_retries() {
                let (mut stream, _) = listener.accept().unwrap();
                read_request(&stream);
                let body = "<Error><Message>private</Message>";
                write!(
                    stream,
                    "HTTP/1.1 403 Forbidden\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{body}",
                    body.len()
                )
                .unwrap();
            }
        });
        let config = S3Config::new("test", Some(&endpoint), true, credentials).unwrap();
        let mut source = RangeSource::open("s3://bucket/object", 34, Some(&config)).unwrap();
        assert_eq!(
            source.read(&mut [0; 34]).unwrap_err().to_string(),
            "S3 range request failed"
        );
        server.join().unwrap();
    }
}
