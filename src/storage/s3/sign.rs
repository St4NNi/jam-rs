//! AWS Signature Version 4 for a GET request with a range, as S3 and its compatible stores
//! check it. Every header that matters is signed, so a proxy cannot change the range.

use super::credentials::Credentials;
use hmac::{Hmac, Mac};
use sha2::{Digest, Sha256};
use std::time::{SystemTime, UNIX_EPOCH};

/// The hash of an empty body, which every GET sends.
const EMPTY_PAYLOAD_SHA256: &str =
    "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855";

/// What is signed: one GET of a byte range.
pub(super) struct Get<'a> {
    /// The host name, with the port unless it is the default one.
    pub host: &'a str,
    /// The path as sent, already percent-encoded.
    pub path: &'a str,
    /// The value of the `Range` header, for example `bytes=0-9`.
    pub range: &'a str,
    /// The value of the `If-Match` header, if the read must fail for another object version.
    pub if_match: Option<&'a str>,
    pub region: &'a str,
}

/// Returns the headers to send with the request, the authorization among them. `host` is left
/// out: the HTTP client sends it from the same address.
pub(super) fn headers(
    get: &Get<'_>,
    credentials: &Credentials,
    now: SystemTime,
) -> Vec<(&'static str, String)> {
    let (date, stamp) = date_and_stamp(now);
    // Headers in the order of their lower-case names, as the canonical request needs them.
    let mut signed: Vec<(&'static str, String)> = vec![("host", get.host.to_string())];
    if let Some(etag) = get.if_match {
        signed.push(("if-match", etag.to_string()));
    }
    signed.push(("range", get.range.to_string()));
    signed.push(("x-amz-content-sha256", EMPTY_PAYLOAD_SHA256.to_string()));
    signed.push(("x-amz-date", stamp.clone()));
    if let Some(token) = &credentials.session_token {
        signed.push(("x-amz-security-token", token.clone()));
    }
    let names = signed
        .iter()
        .map(|(name, _)| *name)
        .collect::<Vec<_>>()
        .join(";");
    let block: String = signed
        .iter()
        .map(|(name, value)| format!("{name}:{}\n", value.trim()))
        .collect();
    let request = format!(
        "GET\n{}\n\n{block}\n{names}\n{EMPTY_PAYLOAD_SHA256}",
        get.path
    );
    let scope = format!("{date}/{}/s3/aws4_request", get.region);
    let to_sign = format!(
        "AWS4-HMAC-SHA256\n{stamp}\n{scope}\n{}",
        hex(&Sha256::digest(request))
    );
    let mut key = mac(
        format!("AWS4{}", credentials.secret_key).as_bytes(),
        date.as_bytes(),
    );
    for part in [get.region, "s3", "aws4_request"] {
        key = mac(&key, part.as_bytes());
    }
    let signature = hex(&mac(&key, to_sign.as_bytes()));
    let authorization = format!(
        "AWS4-HMAC-SHA256 Credential={}/{scope}, SignedHeaders={names}, Signature={signature}",
        credentials.access_key
    );
    signed.retain(|(name, _)| *name != "host");
    signed.push(("authorization", authorization));
    signed
}

fn mac(key: &[u8], data: &[u8]) -> Vec<u8> {
    // HMAC takes keys of any length, so creating it cannot fail.
    let mut mac = Hmac::<Sha256>::new_from_slice(key).expect("HMAC accepts any key length");
    mac.update(data);
    mac.finalize().into_bytes().to_vec()
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

/// The date (`20130524`) and the time stamp (`20130524T000000Z`) of `now` in UTC.
fn date_and_stamp(now: SystemTime) -> (String, String) {
    let seconds = now
        .duration_since(UNIX_EPOCH)
        .map_or(0, |since| since.as_secs());
    let (year, month, day) = civil_from_days((seconds / 86_400) as i64);
    let of_day = seconds % 86_400;
    let date = format!("{year:04}{month:02}{day:02}");
    let stamp = format!(
        "{date}T{:02}{:02}{:02}Z",
        of_day / 3600,
        of_day % 3600 / 60,
        of_day % 60
    );
    (date, stamp)
}

/// The year, month and day of a day count since 1970-01-01 (proleptic Gregorian calendar).
fn civil_from_days(days: i64) -> (i64, i64, i64) {
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let day_of_era = z.rem_euclid(146_097);
    let year_of_era =
        (day_of_era - day_of_era / 1460 + day_of_era / 36_524 - day_of_era / 146_096) / 365;
    let day_of_year = day_of_era - (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
    let shifted_month = (5 * day_of_year + 2) / 153;
    let day = day_of_year - (153 * shifted_month + 2) / 5 + 1;
    let month = if shifted_month < 10 {
        shifted_month + 3
    } else {
        shifted_month - 9
    };
    let year = year_of_era + era * 400 + i64::from(month <= 2);
    (year, month, day)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    /// The ranged GET example of the AWS documentation for Signature Version 4.
    #[test]
    fn matches_the_example_of_the_aws_documentation() {
        let credentials = Credentials {
            access_key: "AKIAIOSFODNN7EXAMPLE".into(),
            secret_key: "wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY".into(),
            session_token: None,
        };
        let get = Get {
            host: "examplebucket.s3.amazonaws.com",
            path: "/test.txt",
            range: "bytes=0-9",
            if_match: None,
            region: "us-east-1",
        };
        let now = UNIX_EPOCH + Duration::from_secs(1_369_353_600);
        let headers = headers(&get, &credentials, now);
        let header = |name: &str| {
            headers
                .iter()
                .find(|(key, _)| *key == name)
                .unwrap()
                .1
                .clone()
        };
        assert_eq!(header("x-amz-date"), "20130524T000000Z");
        assert_eq!(
            header("authorization"),
            "AWS4-HMAC-SHA256 Credential=AKIAIOSFODNN7EXAMPLE/20130524/us-east-1/s3/aws4_request, \
             SignedHeaders=host;range;x-amz-content-sha256;x-amz-date, \
             Signature=f0e8bdb87c964420e857bd35b5d6ed310bd44f0170aba48dd91039c6036bdb41"
        );
    }

    #[test]
    fn session_token_and_etag_are_signed() {
        let credentials = Credentials {
            access_key: "key".into(),
            secret_key: "secret".into(),
            session_token: Some("token".into()),
        };
        let get = Get {
            host: "h",
            path: "/b/k",
            range: "bytes=1-2",
            if_match: Some("\"tag\""),
            region: "r",
        };
        let headers = headers(&get, &credentials, SystemTime::now());
        let authorization = &headers
            .iter()
            .find(|(key, _)| *key == "authorization")
            .unwrap()
            .1;
        assert!(authorization.contains(
            "SignedHeaders=host;if-match;range;x-amz-content-sha256;x-amz-date;x-amz-security-token"
        ));
        assert!(
            headers
                .iter()
                .any(|(key, value)| *key == "if-match" && value == "\"tag\"")
        );
        assert!(headers.iter().all(|(key, _)| *key != "host"));
    }

    #[test]
    fn dates_come_out_right() {
        assert_eq!(civil_from_days(0), (1970, 1, 1));
        assert_eq!(civil_from_days(19_782), (2024, 2, 29));
        assert_eq!(
            date_and_stamp(UNIX_EPOCH + Duration::from_secs(86_399)).1,
            "19700101T235959Z"
        );
    }
}
