#[path = "support/remote/mod.rs"]
mod remote;

use jam_rs::resource::object::ObjectResource;
use jam_rs::resource::{ByteRange, RangeReader, ResourceOpenOptions, ResourceScheme};
use remote::HttpFixture;
use std::sync::{Mutex, OnceLock};

static S3_ENV_LOCK: OnceLock<Mutex<()>> = OnceLock::new();

#[test]
fn optional_s3_endpoint_maps_bucket_and_key_to_checked_remote_reads() {
    let _guard = S3_ENV_LOCK.get_or_init(|| Mutex::new(())).lock().unwrap();
    let fixture = HttpFixture::new(b"s3-compatible-object".to_vec());
    let endpoint = fixture.url("");
    let previous = std::env::var_os("JAM_S3_ENDPOINT");
    // Rust 2024 makes process-environment mutation explicitly unsafe because
    // other threads may observe it.  The test lock keeps this fixture's
    // mutation isolated within this integration-test process.
    unsafe { std::env::set_var("JAM_S3_ENDPOINT", &endpoint) };

    let opened = ObjectResource::open(
        "s3://test-bucket/assemblies/assembly.jma",
        ResourceOpenOptions {
            cache_block_bytes: 4,
            max_cache_bytes: 32,
            max_retries: 0,
            ..ResourceOpenOptions::default()
        },
    );
    let resource = opened.expect("open S3-compatible fixture");
    let bytes = resource.read_range(ByteRange::new(4, 7).unwrap());

    match previous {
        Some(value) => unsafe { std::env::set_var("JAM_S3_ENDPOINT", value) },
        None => unsafe { std::env::remove_var("JAM_S3_ENDPOINT") },
    }

    assert_eq!(resource.locator().scheme(), ResourceScheme::S3);
    assert_eq!(
        resource.locator().redacted(),
        "s3://test-bucket/assemblies/assembly.jma"
    );
    assert_eq!(bytes.unwrap(), b"ompatib");
    assert!(fixture.requests().iter().any(|request| {
        request.method == "GET" && request.path == "/test-bucket/assemblies/assembly.jma"
    }));
}
