use jam_rs::format::{VERSION, bucket_id};
use jam_rs::reader::JamReader;
use jam_rs::writer::{HashSampleInput, build_from_hash_samples};
use tempfile::Builder;

#[test]
fn caller_selected_screen_is_a_pure_format_three_jam_database() {
    let directory = Builder::new()
        .prefix("jam-index-screen-")
        .tempdir_in("target")
        .unwrap();
    let output = directory.path().join("screen.jam");
    let samples = vec![
        HashSampleInput {
            sample_name: "mg-a".to_string(),
            hashes: vec![11, 7, 11, 19],
        },
        HashSampleInput {
            sample_name: "mg-b".to_string(),
            hashes: vec![7, 23],
        },
    ];
    let stats = build_from_hash_samples(&output, &samples, 21, 1).unwrap();
    let reader = JamReader::open(&output).unwrap();
    assert_eq!(VERSION, 3);
    assert_eq!(reader.kmer_size(), 21);
    assert_eq!(reader.threshold(), u64::MAX);
    assert_eq!(reader.sample_names(), ["mg-a", "mg-b"]);
    assert_eq!(reader.sample_sizes(), [3, 2]);
    assert_eq!(stats.total_entries, 5);
    let entries = reader.bucket_entries(bucket_id(7));
    assert!(
        entries
            .iter()
            .any(|entry| entry.hash == 7 && entry.sample_id == 0)
    );
    assert!(
        entries
            .iter()
            .any(|entry| entry.hash == 7 && entry.sample_id == 1)
    );
}

#[test]
fn selected_screen_rejects_zero_hash_and_duplicate_names() {
    let directory = Builder::new()
        .prefix("jam-index-screen-invalid-")
        .tempdir_in("target")
        .unwrap();
    let output = directory.path().join("screen.jam");
    assert!(
        build_from_hash_samples(
            &output,
            &[HashSampleInput {
                sample_name: "mg".to_string(),
                hashes: vec![0],
            }],
            21,
            1,
        )
        .is_err()
    );
}
