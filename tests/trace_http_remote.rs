#[path = "support/remote/mod.rs"]
mod remote;

use jam_rs::jma::ArchiveReader;
use jam_rs::jma::SequenceRange;
use jam_rs::jma::builder::{ArchiveBuildConfig, write_archive_from_fasta};
use jam_rs::jma::reader::JmaReader;
use jam_rs::resource::cache::BlockCache;
use jam_rs::resource::object::ObjectResource;
use jam_rs::resource::{ByteRange, CacheIdentity, RangeReader, ResourceError, ResourceOpenOptions};
use jam_rs::trace::raw::RawAssembly;
use remote::HttpFixture;
use std::fs;
use std::time::Duration;

fn temporary_directory(prefix: &str) -> tempfile::TempDir {
    tempfile::Builder::new()
        .prefix(prefix)
        .tempdir_in("target")
        .expect("temporary test directory")
}

fn archive_bytes() -> Vec<u8> {
    let directory = temporary_directory("trace-http-jma-");
    let sequence: String = (0..240)
        .map(|index| b"ACGT"[(index * 17 + index / 5) % 4] as char)
        .collect();
    let input = directory.path().join("assembly.fa");
    let output = directory.path().join("assembly.jma");
    fs::write(&input, format!(">remote-contig\n{sequence}\n")).expect("write JMA input");
    write_archive_from_fasta(
        &input,
        &output,
        ArchiveBuildConfig {
            block_bases: 32,
            k31_scale: 1,
            k21_scale: Some(1),
            ..ArchiveBuildConfig::default()
        },
    )
    .expect("build JMA fixture");
    fs::read(output).expect("read JMA fixture")
}

fn small_options() -> ResourceOpenOptions {
    ResourceOpenOptions {
        cache_block_bytes: 4,
        max_cache_bytes: 64,
        max_retries: 0,
        ..ResourceOpenOptions::default()
    }
}

#[test]
fn http_jma_reader_uses_checked_ranges_and_reusable_cache() {
    let fixture = HttpFixture::new(archive_bytes());
    let resource = ObjectResource::open(
        fixture.url("/assembly.jma"),
        ResourceOpenOptions {
            cache_block_bytes: 64,
            max_cache_bytes: 1024 * 1024,
            max_retries: 0,
            ..ResourceOpenOptions::default()
        },
    )
    .expect("open remote JMA");
    let reader = JmaReader::open(resource).expect("read remote JMA");

    assert_eq!(reader.header().format_version, 1);
    assert_eq!(reader.contigs()[0].name, "remote-contig");
    assert_eq!(
        reader
            .read_sequence(0, SequenceRange::new(11, 43).unwrap())
            .unwrap()
            .len(),
        32
    );

    let metrics = reader.resource().metrics();
    assert!(metrics.range_requests > 0);
    assert!(metrics.remote_bytes > 0);
    let requests = fixture.requests();
    assert!(requests.iter().any(|request| request.method == "HEAD"));
    assert!(
        requests
            .iter()
            .any(|request| request.method == "GET" && request.range.is_some())
    );

    // JMA construction reads complete required sections today.  Repeating a
    // direct range read still demonstrates that remote bytes are served from
    // the bounded block cache rather than fetched a second time.
    let cached_resource = ObjectResource::open(
        fixture.url("/assembly.jma"),
        ResourceOpenOptions {
            cache_block_bytes: 16,
            max_cache_bytes: 64,
            max_retries: 0,
            ..ResourceOpenOptions::default()
        },
    )
    .expect("open cache fixture");
    let range = ByteRange::new(0, 12).unwrap();
    let first = cached_resource.read_range(range).unwrap();
    let request_count = fixture.requests().len();
    let second = cached_resource.read_range(range).unwrap();
    assert_eq!(first, second);
    assert_eq!(fixture.requests().len(), request_count);
    assert!(cached_resource.metrics().cache_hits > 0);
}

#[test]
fn http_no_range_servers_use_explicit_full_download_fallback() {
    let body = b"0123456789abcdef".to_vec();
    let fixture = HttpFixture::new(body.clone()).without_ranges();
    let resource = ObjectResource::open(
        fixture.url("/assembly.fa"),
        ResourceOpenOptions {
            cache_block_bytes: 8,
            max_cache_bytes: 32,
            allow_full_download_fallback: true,
            ..small_options()
        },
    )
    .unwrap();
    assert!(!resource.metadata().unwrap().accepts_ranges);
    assert_eq!(
        resource.read_range(ByteRange::new(3, 5).unwrap()).unwrap(),
        b"34567"
    );

    let strict = ObjectResource::open(
        fixture.url("/assembly.fa"),
        ResourceOpenOptions {
            allow_full_download_fallback: false,
            ..small_options()
        },
    )
    .unwrap();
    assert!(matches!(
        strict.read_range(ByteRange::new(3, 5).unwrap()),
        Err(ResourceError::RangeUnsupported(_))
    ));
    assert!(
        fixture
            .requests()
            .iter()
            .any(|request| request.method == "GET" && request.range.is_none())
    );
}

#[test]
fn http_range_failures_are_retried_and_interrupted_reads_are_reported() {
    let retry_fixture = HttpFixture::new(b"retryable-body".to_vec()).fail_range_requests(1);
    let retry_resource = ObjectResource::open(
        retry_fixture.url("/object"),
        ResourceOpenOptions {
            max_retries: 2,
            ..small_options()
        },
    )
    .unwrap();
    assert_eq!(
        retry_resource
            .read_range(ByteRange::new(0, 4).unwrap())
            .unwrap(),
        b"retr"
    );
    assert!(retry_resource.metrics().retries > 0);
    assert!(retry_fixture.range_requests().len() >= 2);

    let interrupted_fixture =
        HttpFixture::new(b"interrupted-body".to_vec()).truncate_range_after(1);
    let interrupted_resource = ObjectResource::open(
        interrupted_fixture.url("/object"),
        ResourceOpenOptions {
            max_retries: 0,
            ..small_options()
        },
    )
    .unwrap();
    let error = interrupted_resource
        .read_range(ByteRange::new(0, 4).unwrap())
        .expect_err("short HTTP body must not be accepted as a complete range");
    assert!(matches!(error, ResourceError::Transport { .. }));
}

#[test]
fn changed_remote_identity_does_not_reuse_old_blocks() {
    let fixture = HttpFixture::new(b"old-version".to_vec());
    let first = ObjectResource::open(fixture.url("/object"), small_options()).unwrap();
    assert_eq!(
        first.read_range(ByteRange::new(0, 3).unwrap()).unwrap(),
        b"old"
    );

    fixture.set_body(b"new-version".to_vec(), "fixture-v2");
    let second = ObjectResource::open(fixture.url("/object"), small_options()).unwrap();
    assert_eq!(
        second.read_range(ByteRange::new(0, 3).unwrap()).unwrap(),
        b"new"
    );

    let cache = BlockCache::new(4, 16).unwrap();
    let version_one = CacheIdentity {
        redacted_locator: "http://fixture/object".to_string(),
        version: "fixture-v1".to_string(),
        size: 11,
    };
    let version_two = CacheIdentity {
        version: "fixture-v2".to_string(),
        ..version_one.clone()
    };
    cache.prepare(&version_one).unwrap();
    cache.insert(&version_one, 0, b"old-".to_vec()).unwrap();
    assert!(cache.get(&version_one, 0).unwrap().is_some());
    cache.prepare(&version_two).unwrap();
    assert!(cache.get(&version_two, 0).unwrap().is_none());
}

#[test]
fn remote_raw_streaming_preserves_contigs_and_redacts_query_tokens() {
    let body = b">one\nACGTACGTACGT\n>two\nTTGCAATTGCAA\n".to_vec();
    let fixture = HttpFixture::new(body);
    let resource = RawAssembly::open(
        fixture.url("/assembly.fa?token=do-not-log"),
        ResourceOpenOptions::default(),
    )
    .unwrap();
    assert_eq!(resource.contigs.len(), 2);
    assert_eq!(resource.contig("one").unwrap().sequence, b"ACGTACGTACGT");
    assert_eq!(resource.contig("two").unwrap().sequence, b"TTGCAATTGCAA");
    assert_eq!(resource.redacted_locator, fixture.url("/assembly.fa"));
    assert_eq!(resource.metrics.stream_requests, 1);
    assert!(
        fixture
            .requests()
            .iter()
            .any(|request| request.method == "GET" && request.range.is_none())
    );
}

#[test]
fn signed_http_failures_do_not_expose_credentials_or_query_tokens() {
    // Metadata may try HEAD and then a ranged GET. Keep every bounded attempt
    // failing so this test exercises redaction on an actual transport error.
    let fixture = HttpFixture::new(b"unused".to_vec()).fail_all_requests(16);
    let host_path = fixture
        .url("/missing")
        .strip_prefix("http://")
        .unwrap()
        .to_string();
    let locator = format!("http://user:password@{host_path}?X-Amz-Signature=top-secret-token");
    let resource = ObjectResource::open(
        locator,
        ResourceOpenOptions {
            max_retries: 0,
            ..ResourceOpenOptions::default()
        },
    )
    .unwrap();
    let error = resource
        .metadata()
        .expect_err("fixture must fail transport");
    let message = error.to_string();
    assert!(!message.contains("password"));
    assert!(!message.contains("top-secret-token"));
    assert!(message.contains("/missing"));
}

#[test]
fn http_fixture_records_status_ranges_and_byte_counts() {
    let fixture = HttpFixture::new(b"0123456789abcdef".to_vec());
    let resource = ObjectResource::open(
        fixture.url("/object"),
        ResourceOpenOptions {
            cache_block_bytes: 8,
            max_cache_bytes: 64,
            max_retries: 0,
            ..ResourceOpenOptions::default()
        },
    )
    .unwrap();
    assert_eq!(
        resource.read_range(ByteRange::new(3, 5).unwrap()).unwrap(),
        b"34567"
    );
    let range = fixture
        .range_requests()
        .into_iter()
        .find(|request| request.method == "GET")
        .expect("range GET");
    assert_eq!(range.status, 206);
    assert_eq!(range.range, Some((0, 7)));
    assert_eq!(range.content_range, Some((0, 7, 16)));
    assert_eq!(range.requested_bytes, 8);
    assert_eq!(range.returned_bytes, 8);
}

#[test]
fn malformed_partial_responses_are_visible_to_the_client() {
    let fixture = HttpFixture::new(b"0123456789abcdef".to_vec()).with_malformed_content_range();
    let resource = ObjectResource::open(
        fixture.url("/object"),
        ResourceOpenOptions {
            cache_block_bytes: 8,
            max_cache_bytes: 64,
            max_retries: 0,
            ..ResourceOpenOptions::default()
        },
    )
    .unwrap();
    let result = resource.read_range(ByteRange::new(3, 5).unwrap());
    assert!(result.is_err(), "incorrect Content-Range must be rejected");
    assert!(
        fixture
            .range_requests()
            .iter()
            .any(|request| { request.status == 206 && request.content_range != Some((0, 7, 16)) })
    );
}

#[test]
fn advertised_ranges_that_return_200_are_explicitly_fallback_or_error() {
    let fixture = HttpFixture::new(b"0123456789abcdef".to_vec()).ignoring_ranges();
    let resource = ObjectResource::open(
        fixture.url("/object"),
        ResourceOpenOptions {
            cache_block_bytes: 8,
            max_cache_bytes: 64,
            max_retries: 0,
            allow_full_download_fallback: true,
            ..ResourceOpenOptions::default()
        },
    )
    .unwrap();
    assert_eq!(
        resource.read_range(ByteRange::new(3, 5).unwrap()).unwrap(),
        b"34567"
    );
    assert!(fixture.range_requests().iter().any(|request| {
        request.status == 200 && request.range == Some((0, 7)) && request.returned_bytes == 16
    }));

    let strict = ObjectResource::open(
        fixture.url("/strict"),
        ResourceOpenOptions {
            cache_block_bytes: 8,
            max_cache_bytes: 64,
            max_retries: 0,
            allow_full_download_fallback: false,
            ..ResourceOpenOptions::default()
        },
    )
    .unwrap();
    assert!(strict.read_range(ByteRange::new(3, 5).unwrap()).is_err());
}

#[test]
fn http_status_failures_are_recorded_without_leaking_locator_credentials() {
    for status in [403, 404, 503] {
        let fixture = HttpFixture::new(b"unavailable".to_vec()).with_status(status);
        let host_path = fixture
            .url("/object")
            .strip_prefix("http://")
            .unwrap()
            .to_string();
        let resource = ObjectResource::open(
            format!("http://user:secret@{host_path}?token=do-not-leak"),
            ResourceOpenOptions {
                max_retries: 0,
                ..ResourceOpenOptions::default()
            },
        )
        .unwrap();
        let error = resource.metadata().expect_err("status must fail");
        let message = error.to_string();
        assert!(!message.contains("secret"));
        assert!(!message.contains("do-not-leak"));
        assert_eq!(fixture.statuses().last().copied(), Some(status));
    }
}

#[test]
fn delayed_http_response_exercises_timeout_path() {
    let fixture = HttpFixture::new(b"slow".to_vec()).delay_response(Duration::from_secs(2));
    let resource = ObjectResource::open(
        fixture.url("/slow"),
        ResourceOpenOptions {
            request_timeout_seconds: 1,
            max_retries: 0,
            ..ResourceOpenOptions::default()
        },
    )
    .unwrap();
    assert!(resource.metadata().is_err());
    assert!(!fixture.requests().is_empty());
}
