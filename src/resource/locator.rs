//! Locator parsing and redaction helpers.

use super::{ResourceError, ResourceLocator, ResourceResult, ResourceScheme};

pub(super) fn parse(input: &str) -> ResourceResult<ResourceLocator> {
    let input = input.trim();
    if input.is_empty() {
        return Err(ResourceError::InvalidLocator("empty locator".to_string()));
    }
    let scheme = if input.starts_with("file://") {
        ResourceScheme::File
    } else if input.starts_with("http://") {
        ResourceScheme::Http
    } else if input.starts_with("https://") {
        ResourceScheme::Https
    } else if input.starts_with("s3://") {
        ResourceScheme::S3
    } else if input.contains("://") {
        return Err(ResourceError::UnsupportedScheme(
            input.split("://").next().unwrap_or_default().to_string(),
        ));
    } else {
        ResourceScheme::Local
    };
    validate(input, scheme)?;
    Ok(ResourceLocator::from_parts(scheme, input.to_string()))
}

pub(super) fn redact(locator: &ResourceLocator) -> String {
    match locator.scheme() {
        ResourceScheme::Local => locator.raw().to_string(),
        ResourceScheme::File => strip_suffix(locator.raw()).to_string(),
        ResourceScheme::Http | ResourceScheme::Https | ResourceScheme::S3 => {
            let (prefix, remainder) = locator
                .raw()
                .split_once("://")
                .unwrap_or(("", locator.raw()));
            let without_query = strip_suffix(remainder);
            let (authority, path) = without_query
                .split_once('/')
                .map_or((without_query, ""), |(authority, path)| (authority, path));
            let without_credentials = authority
                .rsplit_once('@')
                .map_or(authority, |(_, host)| host);
            let without_credentials = if path.is_empty() {
                without_credentials.to_string()
            } else {
                format!("{without_credentials}/{path}")
            };
            format!("{prefix}://{without_credentials}")
        }
    }
}

fn strip_suffix(value: &str) -> &str {
    let value = value.split_once('#').map_or(value, |(prefix, _)| prefix);
    value.split_once('?').map_or(value, |(prefix, _)| prefix)
}

fn validate(input: &str, scheme: ResourceScheme) -> ResourceResult<()> {
    match scheme {
        ResourceScheme::Local => Ok(()),
        ResourceScheme::File => {
            let remainder = strip_suffix(input.strip_prefix("file://").unwrap_or_default());
            if remainder.is_empty() {
                return Err(ResourceError::InvalidLocator(
                    "file URL has no path".to_string(),
                ));
            }
            Ok(())
        }
        ResourceScheme::Http | ResourceScheme::Https => {
            let remainder = input
                .split_once("://")
                .map(|(_, value)| value)
                .unwrap_or_default();
            let remainder = strip_suffix(remainder);
            let authority = remainder
                .split_once('/')
                .map_or(remainder, |(host, _)| host)
                .trim();
            let host = authority
                .rsplit_once('@')
                .map_or(authority, |(_, host)| host);
            if host.is_empty() || host.starts_with(':') || host.chars().any(char::is_whitespace) {
                return Err(ResourceError::InvalidLocator(
                    "HTTP URL has no host".to_string(),
                ));
            }
            Ok(())
        }
        ResourceScheme::S3 => {
            let remainder = input.strip_prefix("s3://").unwrap_or_default();
            let remainder = strip_suffix(remainder);
            let (bucket, key) = remainder.split_once('/').unwrap_or((remainder, ""));
            if bucket.is_empty() || key.is_empty() {
                return Err(ResourceError::InvalidLocator(
                    "S3 URL must contain a bucket and object key".to_string(),
                ));
            }
            let bucket = bucket.trim();
            if bucket.is_empty() || bucket.chars().any(char::is_whitespace) {
                return Err(ResourceError::InvalidLocator(
                    "S3 URL must contain a bucket and object key".to_string(),
                ));
            }
            Ok(())
        }
    }
}
