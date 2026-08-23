use jam_rs::provenance::{command_line, redact_command_arguments, redacted_command_line};

fn strings(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| (*value).to_string()).collect()
}

#[test]
fn preserves_ordinary_arguments_and_option_structure() {
    let input = strings(&["jam", "trace", "--threads", "8", "assembly.fa"]);
    assert_eq!(redact_command_arguments(&input), input);
}

#[test]
fn redacts_signed_http_urls_and_userinfo() {
    let input = strings(&[
        "jam",
        "https://user:password@example.org/archive.jma?X-Amz-Signature=secret#fragment",
        "--database=https://user:password@example.org/archive.jma?token=secret",
    ]);
    assert_eq!(
        redact_command_arguments(&input),
        strings(&[
            "jam",
            "https://example.org/archive.jma",
            "--database=https://example.org/archive.jma",
        ])
    );
}

#[test]
fn redacts_signed_s3_urls_in_both_forms() {
    let input = strings(&[
        "s3://access:secret@bucket/archive.jma?X-Amz-Credential=secret",
        "--database=s3://bucket/archive.jma?signature=secret#fragment",
    ]);
    assert_eq!(
        redact_command_arguments(&input),
        strings(&[
            "s3://bucket/archive.jma",
            "--database=s3://bucket/archive.jma"
        ])
    );
}

#[test]
fn redacts_file_url_suffix_without_changing_local_paths() {
    let input = strings(&["file:///tmp/archive.jma?secret=secret", "./archive.jma"]);
    assert_eq!(
        redact_command_arguments(&input),
        strings(&["file:///tmp/archive.jma", "./archive.jma"])
    );
}

#[test]
fn malformed_urls_do_not_echo_sensitive_suffixes() {
    let input = strings(&[
        "https://user:password@[invalid.example/archive.jma?token=secret",
        "--database=https://user:password@/archive.jma?signature=secret",
    ]);
    let output = redact_command_arguments(&input);
    assert!(output.iter().all(|value| !value.contains("password")));
    assert!(output.iter().all(|value| !value.contains("secret")));
}

#[test]
fn legacy_command_line_remains_available() {
    assert!(!command_line().is_empty());
    assert!(!redacted_command_line().is_empty());
}
