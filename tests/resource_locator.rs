use jam_rs::resource::{ResourceLocator, ResourceScheme};

#[test]
fn parses_supported_locator_schemes() {
    assert_eq!(
        ResourceLocator::parse("relative.jma").unwrap().scheme(),
        ResourceScheme::Local
    );
    assert_eq!(
        ResourceLocator::parse("file:///tmp/archive.jma")
            .unwrap()
            .scheme(),
        ResourceScheme::File
    );
    assert_eq!(
        ResourceLocator::parse("https://example.org/archive.jma")
            .unwrap()
            .scheme(),
        ResourceScheme::Https
    );
    assert_eq!(
        ResourceLocator::parse("s3://bucket/archive.jma")
            .unwrap()
            .scheme(),
        ResourceScheme::S3
    );
}

#[test]
fn rejects_malformed_remote_locators() {
    assert!(ResourceLocator::parse("").is_err());
    assert!(ResourceLocator::parse("https:///archive.jma").is_err());
    assert!(ResourceLocator::parse("s3://bucket").is_err());
    assert!(ResourceLocator::parse("gopher://example.org/archive").is_err());
}

#[test]
fn redacts_credentials_queries_and_fragments() {
    let locator = ResourceLocator::parse(
        "https://user:secret@example.org/archive.jma?X-Amz-Signature=secret#fragment",
    )
    .unwrap();
    assert_eq!(locator.redacted(), "https://example.org/archive.jma");

    let s3 = ResourceLocator::parse("s3://bucket/archive.jma?signature=secret").unwrap();
    assert_eq!(s3.redacted(), "s3://bucket/archive.jma");
}
