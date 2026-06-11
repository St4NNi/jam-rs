//! Static S3 credentials from the environment or the AWS files. Roles, web identity, SSO and
//! credential processes need more than this layer does, so they are refused with a clear error
//! instead of falling back to another account.

use crate::storage::StorageError;
use std::ffi::OsString;
use std::path::PathBuf;

/// An access key with its secret and an optional session token.
#[derive(Clone)]
pub struct Credentials {
    /// The access key id.
    pub access_key: String,
    /// The secret access key.
    pub secret_key: String,
    /// The session token of temporary credentials.
    pub session_token: Option<String>,
}

impl std::fmt::Debug for Credentials {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // The secrets never reach a log.
        formatter
            .debug_struct("Credentials")
            .field("access_key", &self.access_key)
            .finish_non_exhaustive()
    }
}

/// Settings that select a credential provider this layer does not implement.
const UNSUPPORTED_ENVIRONMENT: [&str; 5] = [
    "AWS_ROLE_ARN",
    "AWS_WEB_IDENTITY_TOKEN_FILE",
    "AWS_CONTAINER_CREDENTIALS_RELATIVE_URI",
    "AWS_CONTAINER_CREDENTIALS_FULL_URI",
    "AWS_EC2_METADATA_SERVICE_ENDPOINT",
];

/// Keys of a profile that select another provider.
fn unsupported_key(key: &str) -> bool {
    matches!(
        key,
        "role_arn"
            | "source_profile"
            | "credential_source"
            | "web_identity_token_file"
            | "credential_process"
    ) || key.starts_with("sso_")
}

fn refused(reason: &str) -> StorageError {
    StorageError::Credentials(reason.to_string())
}

impl Credentials {
    /// Reads credentials from `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY` (with
    /// `AWS_SESSION_TOKEN`), or else from the profile in the AWS credentials file. Partial
    /// environment credentials fail and never select a profile.
    pub fn from_environment() -> Result<Self, StorageError> {
        Self::from_variables(&|name| std::env::var_os(name))
    }

    /// `from_environment` with the variables given by `lookup`, so that tests need not change
    /// the environment of the process.
    fn from_variables(lookup: &dyn Fn(&str) -> Option<OsString>) -> Result<Self, StorageError> {
        let present = |name: &str| lookup(name).is_some();
        if UNSUPPORTED_ENVIRONMENT.into_iter().any(present) {
            return Err(refused(
                "roles, web identity and instance metadata are not supported",
            ));
        }
        let variable = |name: &str| {
            lookup(name)
                .and_then(|value| value.into_string().ok())
                .filter(|value| !value.is_empty())
        };
        let keys = [
            "AWS_ACCESS_KEY_ID",
            "AWS_SECRET_ACCESS_KEY",
            "AWS_SESSION_TOKEN",
            "AWS_SECURITY_TOKEN",
        ];
        if keys.into_iter().any(present) {
            let token = variable("AWS_SESSION_TOKEN").or_else(|| variable("AWS_SECURITY_TOKEN"));
            return match (
                variable("AWS_ACCESS_KEY_ID"),
                variable("AWS_SECRET_ACCESS_KEY"),
            ) {
                (Some(access_key), Some(secret_key)) => Ok(Self {
                    access_key,
                    secret_key,
                    session_token: token,
                }),
                _ => Err(refused(
                    "the environment holds only part of the credentials",
                )),
            };
        }
        let profile = variable("AWS_PROFILE")
            .or_else(|| variable("AWS_DEFAULT_PROFILE"))
            .unwrap_or_else(|| "default".to_string());
        Self::from_profile(&profile, lookup)
    }

    /// Reads the profile `profile` of the files named by `AWS_SHARED_CREDENTIALS_FILE` and
    /// `AWS_CONFIG_FILE`, or `~/.aws/credentials` and `~/.aws/config`.
    fn from_profile(
        profile: &str,
        lookup: &dyn Fn(&str) -> Option<OsString>,
    ) -> Result<Self, StorageError> {
        let file = |variable: &str, default: &str| {
            lookup(variable)
                .map(PathBuf::from)
                .or_else(|| Some(PathBuf::from(lookup("HOME")?).join(".aws").join(default)))
        };
        let credentials_path = file("AWS_SHARED_CREDENTIALS_FILE", "credentials")
            .ok_or_else(|| refused("no home directory for the AWS credentials file"))?;
        let credentials = std::fs::read_to_string(&credentials_path)
            .map_err(|_| refused("the AWS credentials file cannot be read"))?;
        // A missing default config file is fine; a named one must exist.
        let config = match file("AWS_CONFIG_FILE", "config").map(std::fs::read_to_string) {
            Some(Ok(text)) => text,
            Some(Err(error))
                if lookup("AWS_CONFIG_FILE").is_none()
                    && error.kind() == std::io::ErrorKind::NotFound =>
            {
                String::new()
            }
            Some(Err(_)) => return Err(refused("the AWS config file cannot be read")),
            None => String::new(),
        };
        let config_section = if profile == "default" {
            profile.to_string()
        } else {
            format!("profile {profile}")
        };
        let selected = section(&credentials, profile);
        let others = section(&config, &config_section);
        if selected
            .iter()
            .chain(&others)
            .any(|(key, _)| unsupported_key(key))
        {
            return Err(refused(
                "the profile selects roles, SSO or a credential process",
            ));
        }
        let get = |key: &str| {
            selected
                .iter()
                .find(|(name, _)| name == key)
                .map(|(_, value)| value.clone())
        };
        match (get("aws_access_key_id"), get("aws_secret_access_key")) {
            (Some(access_key), Some(secret_key))
                if !access_key.is_empty() && !secret_key.is_empty() =>
            {
                Ok(Self {
                    access_key,
                    secret_key,
                    session_token: get("aws_session_token").or_else(|| get("aws_security_token")),
                })
            }
            _ => Err(refused("the profile has no access key and secret")),
        }
    }
}

/// The `key = value` pairs of the section `name` of an INI text. Comments start with `#` or `;`.
fn section(text: &str, name: &str) -> Vec<(String, String)> {
    let mut inside = false;
    let mut pairs = Vec::new();
    for line in text
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with(['#', ';']))
    {
        if let Some(header) = line
            .strip_prefix('[')
            .and_then(|rest| rest.strip_suffix(']'))
        {
            inside = header.trim() == name;
        } else if inside && let Some((key, value)) = line.split_once('=') {
            pairs.push((key.trim().to_string(), value.trim().to_string()));
        }
    }
    pairs
}

#[cfg(test)]
mod tests {
    use super::*;

    const FILE: &str = "# comment\n[default]\naws_access_key_id = AAA\naws_secret_access_key=SSS\n\n[other]\naws_access_key_id = BBB\naws_secret_access_key = TTT\naws_session_token = ZZZ\nrole_arn = x\n";

    #[test]
    fn sections_are_read_by_name() {
        assert_eq!(section(FILE, "default").len(), 2);
        assert_eq!(
            section(FILE, "other")[2],
            ("aws_session_token".to_string(), "ZZZ".to_string())
        );
        assert!(section(FILE, "missing").is_empty());
        assert!(
            section(FILE, "other")
                .iter()
                .any(|(key, _)| unsupported_key(key))
        );
        assert!(
            !section(FILE, "default")
                .iter()
                .any(|(key, _)| unsupported_key(key))
        );
    }

    #[test]
    fn secrets_stay_out_of_debug_output() {
        let credentials = Credentials {
            access_key: "AAA".into(),
            secret_key: "SSS".into(),
            session_token: Some("TTT".into()),
        };
        let shown = format!("{credentials:?}");
        assert!(shown.contains("AAA") && !shown.contains("SSS") && !shown.contains("TTT"));
    }

    fn lookup<'a>(pairs: &'a [(&'a str, String)]) -> impl Fn(&str) -> Option<OsString> + 'a {
        move |name| {
            pairs
                .iter()
                .find(|(key, _)| *key == name)
                .map(|(_, value)| value.into())
        }
    }

    #[test]
    fn environment_keys_win_and_partial_ones_fail() {
        let full = [
            ("AWS_ACCESS_KEY_ID", "A".to_string()),
            ("AWS_SECRET_ACCESS_KEY", "S".to_string()),
            ("AWS_SESSION_TOKEN", "T".to_string()),
        ];
        let credentials = Credentials::from_variables(&lookup(&full)).unwrap();
        assert_eq!(
            (
                credentials.access_key.as_str(),
                credentials.session_token.as_deref()
            ),
            ("A", Some("T"))
        );
        let partial = [("AWS_ACCESS_KEY_ID", "A".to_string())];
        assert!(Credentials::from_variables(&lookup(&partial)).is_err());
        let role = [
            ("AWS_ROLE_ARN", "arn".to_string()),
            ("AWS_ACCESS_KEY_ID", "A".to_string()),
        ];
        assert!(Credentials::from_variables(&lookup(&role)).is_err());
    }

    #[test]
    fn profiles_come_from_the_credentials_file() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("credentials");
        std::fs::write(&path, FILE).unwrap();
        let file = path.display().to_string();
        let mut variables = vec![("AWS_SHARED_CREDENTIALS_FILE", file)];
        let credentials = Credentials::from_variables(&lookup(&variables)).unwrap();
        assert_eq!(
            (
                credentials.access_key.as_str(),
                credentials.secret_key.as_str()
            ),
            ("AAA", "SSS")
        );
        // The profile "other" selects a role, which is refused.
        variables.push(("AWS_PROFILE", "other".to_string()));
        assert!(Credentials::from_variables(&lookup(&variables)).is_err());
        variables[1].1 = "missing".to_string();
        assert!(Credentials::from_variables(&lookup(&variables)).is_err());
    }
}
