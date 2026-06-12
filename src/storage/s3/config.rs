//! The settings of S3 requests: where the store is, who asks, and how patient and how parallel
//! the requests are. One config is shared by all objects of a process, so the limit on requests
//! in flight covers them all.

use super::Credentials;
use super::sign::{self, Get};
use crate::storage::{StorageError, limit::InFlight};
use std::io;
use std::sync::Arc;
use std::time::{Duration, SystemTime};

/// Limits on requests to the store.
#[derive(Clone, Copy, Debug)]
pub struct S3Limits {
    /// The longest one request may take, from its start to its last byte.
    pub request_timeout: Duration,
    /// The longest a new connection may take to open. Opening one takes milliseconds, so a
    /// stuck one is cut early and tried again.
    pub connect_timeout: Duration,
    /// Tries of one request, the first included. A read fails after the last one.
    pub tries: u32,
    /// Requests in flight at once, over all objects of the config. Some stores answer "another
    /// range or object version" when 128 range reads run in parallel, so callers can set this
    /// lower; tries with backoff cover the rest.
    pub max_in_flight: usize,
}

impl Default for S3Limits {
    fn default() -> Self {
        Self {
            request_timeout: Duration::from_secs(60),
            connect_timeout: Duration::from_secs(1),
            tries: 4,
            max_in_flight: 64,
        }
    }
}

/// A custom endpoint such as `http://127.0.0.1:9000`, split into scheme and host with port.
#[derive(Clone, Debug)]
struct Endpoint {
    scheme: String,
    host: String,
}

/// Where one object lives: the URL to call, the `Host` header and the signed path.
pub struct Target {
    url: String,
    host: String,
    path: String,
}

/// The place of the store, the credentials and the limits.
#[derive(Clone)]
pub struct S3Config {
    region: String,
    endpoint: Option<Endpoint>,
    path_style: bool,
    credentials: Credentials,
    limits: S3Limits,
    /// Keeps connections to the store open between requests: a new connection costs more than a
    /// range read of a packed sequence, and a search sends many thousands of them.
    agent: ureq::Agent,
    in_flight: Arc<InFlight>,
}

impl S3Config {
    /// A config for `region`. With an `endpoint` (for example `https://s3.example.org`) the
    /// store is not AWS; `path_style` puts the bucket in the path and not in the host name.
    pub fn new(
        region: &str,
        endpoint: Option<&str>,
        path_style: bool,
        credentials: Credentials,
        limits: S3Limits,
    ) -> Result<Self, StorageError> {
        let endpoint = endpoint.map(parse_endpoint).transpose()?;
        if region.is_empty() || limits.tries == 0 {
            return Err(StorageError::Config(
                "the region and the tries must not be empty".into(),
            ));
        }
        let agent = ureq::AgentBuilder::new()
            .max_idle_connections(1024)
            .max_idle_connections_per_host(1024)
            .timeout_connect(limits.connect_timeout)
            .build();
        Ok(Self {
            region: region.to_string(),
            endpoint,
            path_style,
            credentials,
            limits,
            agent,
            in_flight: Arc::new(InFlight::new(limits.max_in_flight)),
        })
    }

    pub fn limits(&self) -> &S3Limits {
        &self.limits
    }

    pub fn in_flight(&self) -> &InFlight {
        &self.in_flight
    }

    /// The target of `key` in `bucket`.
    pub fn target(&self, bucket: &str, key: &str) -> Target {
        let encoded = encode_path(key);
        let (scheme, store) = match &self.endpoint {
            Some(endpoint) => (endpoint.scheme.as_str(), endpoint.host.clone()),
            None => ("https", format!("s3.{}.amazonaws.com", self.region)),
        };
        let (host, path) = if self.path_style {
            (store, format!("/{bucket}/{encoded}"))
        } else {
            (format!("{bucket}.{store}"), format!("/{encoded}"))
        };
        Target {
            url: format!("{scheme}://{host}{path}"),
            host,
            path,
        }
    }

    /// Sends one signed GET of the bytes `start` to `last` (both included) and returns the answer
    /// once its headers have arrived. An answer with an error status is returned too, so that
    /// the caller can tell a temporary status from a final one. A request that cannot be sent
    /// (connection, timeout) is an `Err`.
    pub fn get_range(
        &self,
        target: &Target,
        (start, last): (u64, u64),
        if_match: Option<&str>,
    ) -> io::Result<ureq::Response> {
        let range = format!("bytes={start}-{last}");
        let get = Get {
            host: &target.host,
            path: &target.path,
            range: &range,
            if_match,
            region: &self.region,
        };
        let mut request = self
            .agent
            .get(&target.url)
            .timeout(self.limits.request_timeout);
        for (name, value) in sign::headers(&get, &self.credentials, SystemTime::now()) {
            request = request.set(name, &value);
        }
        match request.call() {
            Ok(response) | Err(ureq::Error::Status(_, response)) => Ok(response),
            Err(ureq::Error::Transport(transport)) => Err(io::Error::other(format!(
                "S3 request failed: {}",
                transport.kind()
            ))),
        }
    }
}

fn parse_endpoint(endpoint: &str) -> Result<Endpoint, StorageError> {
    let invalid = || StorageError::Config(format!("{endpoint:?} is not an http or https endpoint"));
    let (scheme, host) = endpoint.split_once("://").ok_or_else(invalid)?;
    let host = host.trim_end_matches('/');
    if !matches!(scheme, "http" | "https") || host.is_empty() || host.contains('/') {
        return Err(invalid());
    }
    Ok(Endpoint {
        scheme: scheme.to_string(),
        host: host.to_string(),
    })
}

/// Percent-encodes an object key for the path of a request; `/` stays as it is.
fn encode_path(key: &str) -> String {
    let mut encoded = String::with_capacity(key.len());
    for byte in key.bytes() {
        if byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_' | b'.' | b'~' | b'/') {
            encoded.push(byte as char);
        } else {
            encoded.push_str(&format!("%{byte:02X}"));
        }
    }
    encoded
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config(region: &str, endpoint: Option<&str>, path_style: bool) -> S3Config {
        let credentials = Credentials {
            access_key: "a".into(),
            secret_key: "s".into(),
            session_token: None,
        };
        S3Config::new(
            region,
            endpoint,
            path_style,
            credentials,
            S3Limits::default(),
        )
        .unwrap()
    }

    #[test]
    fn urls_follow_the_addressing_style() {
        let aws = config("eu-west-1", None, false).target("bucket", "dir/a b.bin");
        assert_eq!(
            aws.url,
            "https://bucket.s3.eu-west-1.amazonaws.com/dir/a%20b.bin"
        );
        assert_eq!(aws.host, "bucket.s3.eu-west-1.amazonaws.com");
        let path = config("r", Some("http://127.0.0.1:9000/"), true).target("bucket", "key");
        assert_eq!(
            (path.url.as_str(), path.host.as_str(), path.path.as_str()),
            (
                "http://127.0.0.1:9000/bucket/key",
                "127.0.0.1:9000",
                "/bucket/key"
            )
        );
        let virtual_host =
            config("r", Some("https://s3.example.org"), false).target("bucket", "key");
        assert_eq!(virtual_host.url, "https://bucket.s3.example.org/key");
    }

    #[test]
    fn bad_settings_are_rejected() {
        let credentials = || Credentials {
            access_key: "a".into(),
            secret_key: "s".into(),
            session_token: None,
        };
        for endpoint in ["", "ftp://host", "host:9000", "http://", "http://host/path"] {
            assert!(
                S3Config::new(
                    "r",
                    Some(endpoint),
                    true,
                    credentials(),
                    S3Limits::default()
                )
                .is_err(),
                "{endpoint}"
            );
        }
        assert!(S3Config::new("", None, true, credentials(), S3Limits::default()).is_err());
        let no_tries = S3Limits {
            tries: 0,
            ..S3Limits::default()
        };
        assert!(S3Config::new("r", None, true, credentials(), no_tries).is_err());
    }
}
