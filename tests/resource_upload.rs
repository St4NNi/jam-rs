#[path = "support/remote/mod.rs"]
mod remote;

use jam_rs::resource::ResourceError;
use jam_rs::resource::upload::{UploadOptions, upload_file};
use remote::HttpFixture;
use std::sync::{Mutex, OnceLock};
use std::time::Duration;

fn upload_file_fixture(bytes: &[u8]) -> tempfile::NamedTempFile {
    let mut file = tempfile::NamedTempFile::new().expect("create upload fixture");
    std::io::Write::write_all(&mut file, bytes).expect("write upload fixture");
    file
}

#[test]
fn uploads_a_complete_http_file_once() {
    let payload = b"{\"status\":\"complete\"}\n";
    let file = upload_file_fixture(payload);
    let fixture = HttpFixture::new(Vec::<u8>::new());
    let result = upload_file(
        fixture.url("/results/run.jsonl"),
        file.path(),
        UploadOptions {
            max_retries: 0,
            ..UploadOptions::default()
        },
    )
    .expect("upload must succeed");

    assert_eq!(result.bytes_uploaded, payload.len() as u64);
    assert_eq!(result.attempts, 1);
    assert_eq!(result.status, 200);
    assert_eq!(fixture.uploads().len(), 1);
    assert_eq!(fixture.uploads()[0].method, "PUT");
    assert_eq!(fixture.uploads()[0].path, "/results/run.jsonl");
    assert_eq!(fixture.uploads()[0].body, payload);
}

#[test]
fn upload_statuses_are_structured_and_retry_only_for_server_errors() {
    let file = upload_file_fixture(b"result\n");
    for status in [403, 404] {
        let fixture = HttpFixture::new(Vec::<u8>::new()).with_upload_status(status);
        let result = upload_file(
            fixture.url("/results/run.jsonl?token=do-not-leak"),
            file.path(),
            UploadOptions {
                max_retries: 3,
                ..UploadOptions::default()
            },
        );
        assert!(matches!(
            result,
            Err(ResourceError::HttpStatus {
                status: error_status,
                ..
            }) if error_status == status
        ));
        assert_eq!(fixture.uploads().len(), 1);
    }

    let fixture = HttpFixture::new(Vec::<u8>::new()).with_upload_status(503);
    let result = upload_file(
        fixture.url("/results/retry.jsonl"),
        file.path(),
        UploadOptions {
            max_retries: 2,
            ..UploadOptions::default()
        },
    );
    assert!(matches!(
        result,
        Err(ResourceError::HttpStatus { status: 503, .. })
    ));
    assert_eq!(fixture.uploads().len(), 3);
}

#[test]
fn upload_errors_redact_credentials_and_timeout_is_structured() {
    let file = upload_file_fixture(b"result\n");
    let fixture = HttpFixture::new(Vec::<u8>::new()).with_upload_status(403);
    let host_path = fixture
        .url("/results/run.jsonl")
        .strip_prefix("http://")
        .unwrap()
        .to_string();
    let result = upload_file(
        format!("http://user:password@{host_path}?X-Amz-Signature=secret"),
        file.path(),
        UploadOptions {
            max_retries: 0,
            ..UploadOptions::default()
        },
    );
    let error = result.expect_err("forbidden upload must fail");
    let message = error.to_string();
    assert!(!message.contains("password"));
    assert!(!message.contains("secret"));
    assert!(message.contains("/results/run.jsonl"));

    let slow = HttpFixture::new(Vec::<u8>::new()).delay_response(Duration::from_secs(2));
    let result = upload_file(
        slow.url("/results/slow.jsonl"),
        file.path(),
        UploadOptions {
            request_timeout_seconds: 1,
            max_retries: 0,
        },
    );
    assert!(matches!(result, Err(ResourceError::Timeout { .. })));
}

#[test]
fn interrupted_upload_is_not_reported_as_success() {
    let file = upload_file_fixture(b"interrupted-result\n");
    let fixture = HttpFixture::new(Vec::<u8>::new()).interrupt_upload();
    let result = upload_file(
        fixture.url("/results/interrupted.jsonl"),
        file.path(),
        UploadOptions {
            max_retries: 0,
            ..UploadOptions::default()
        },
    );
    assert!(matches!(result, Err(ResourceError::Transport { .. })));
    assert_eq!(fixture.uploads().len(), 1);
    assert_eq!(fixture.requests()[0].status, 0);
}

#[test]
fn uploads_map_s3_locators_to_the_configured_endpoint() {
    static ENV_LOCK: OnceLock<Mutex<()>> = OnceLock::new();
    let _guard = ENV_LOCK.get_or_init(|| Mutex::new(())).lock().unwrap();
    let fixture = HttpFixture::new(Vec::<u8>::new());
    let previous = std::env::var_os("JAM_S3_ENDPOINT");
    unsafe { std::env::set_var("JAM_S3_ENDPOINT", fixture.url("")) };
    let file = upload_file_fixture(b"s3-result\n");
    let result = upload_file(
        "s3://bucket/results/run.jsonl",
        file.path(),
        UploadOptions {
            max_retries: 0,
            ..UploadOptions::default()
        },
    );
    match previous {
        Some(value) => unsafe { std::env::set_var("JAM_S3_ENDPOINT", value) },
        None => unsafe { std::env::remove_var("JAM_S3_ENDPOINT") },
    }
    result.expect("S3-compatible upload must succeed");
    assert_eq!(fixture.uploads()[0].path, "/bucket/results/run.jsonl");
}
