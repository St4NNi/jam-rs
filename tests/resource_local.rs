use jam_rs::resource::local::LocalResource;
use jam_rs::resource::{ByteRange, RangeReader, ResourceOpenOptions};
use std::io::Read;
use std::io::Write;
use tempfile::NamedTempFile;

#[test]
fn local_ranges_are_checked_and_cached() {
    let mut file = NamedTempFile::new().unwrap();
    file.write_all(b"0123456789abcdef").unwrap();
    let options = ResourceOpenOptions {
        cache_block_bytes: 4,
        max_cache_bytes: 12,
        ..ResourceOpenOptions::default()
    };
    let resource = LocalResource::open(file.path().to_string_lossy(), options).unwrap();
    assert_eq!(resource.metadata().unwrap().size, 16);
    assert_eq!(
        resource.read_range(ByteRange::new(3, 7).unwrap()).unwrap(),
        b"3456789"
    );
    let before = resource.metrics();
    assert_eq!(
        resource.read_range(ByteRange::new(3, 7).unwrap()).unwrap(),
        b"3456789"
    );
    let after = resource.metrics();
    assert!(after.cache_hits > before.cache_hits);
    assert!(after.cache_bytes > before.cache_bytes);
    assert!(resource.read_range(ByteRange::new(16, 1).unwrap()).is_err());
    assert!(
        resource
            .read_range(ByteRange {
                offset: u64::MAX,
                length: 1,
            })
            .is_err()
    );

    let mut stream = resource.stream().unwrap();
    let mut bytes = Vec::new();
    stream.read_to_end(&mut bytes).unwrap();
    assert_eq!(bytes, b"0123456789abcdef");
}

#[test]
fn file_url_reads_the_same_bytes() {
    let mut file = NamedTempFile::new().unwrap();
    file.write_all(b"file-url").unwrap();
    let input = format!("file://{}", file.path().display());
    let resource = LocalResource::open(input, ResourceOpenOptions::default()).unwrap();
    assert_eq!(
        resource.read_range(ByteRange::new(0, 8).unwrap()).unwrap(),
        b"file-url"
    );
}
