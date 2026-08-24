use jam_rs::jma::builder::{ArchiveBuildConfig, write_archive_from_fasta};
use jam_rs::jma::reader::JmaReader;
use jam_rs::resource::ResourceOpenOptions;
use jam_rs::resource::local::LocalResource;
use jam_rs::trace::config::SeedSensitivity;
use jam_rs::trace::seeds::extract_seed_level;
use std::fs;

#[test]
fn unequal_packed_kmers_with_the_same_hash_cannot_become_anchors() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("assembly.fa");
    let archive = directory.path().join("assembly.jma");
    let sequence = b"ACGTTGCATGTCAGTAGGCATCAGTACCGATGCTAGCTAGGCTAACGTTACGATCGATGCA";
    fs::write(
        &input,
        format!(">contig\n{}\n", String::from_utf8_lossy(sequence)),
    )
    .unwrap();
    write_archive_from_fasta(
        &input,
        &archive,
        ArchiveBuildConfig {
            k31_scale: 1,
            k21_scale: None,
            ..ArchiveBuildConfig::default()
        },
    )
    .unwrap();

    let resource = LocalResource::from_path(&archive, ResourceOpenOptions::default()).unwrap();
    let reader = JmaReader::from_resource(resource).unwrap();
    let seed = extract_seed_level(
        sequence,
        SeedSensitivity {
            k: 31,
            scale: 1,
            max_occurrences: 4,
        },
    )
    .unwrap()
    .seeds
    .into_iter()
    .next()
    .expect("fixture has a retained k=31 seed");

    let exact = reader.seed_occurrences_at_scale(seed.query(31), 1).unwrap();
    assert!(!exact.is_empty());

    // Reuse the stored hash while changing the packed k-mer. This models a
    // hash collision at the lookup boundary. The JMA reader must require both
    // fields before an occurrence can reach anchor generation.
    let collision_query = jam_rs::jma::SeedQuery {
        k: 31,
        hash: seed.hash,
        canonical_kmer: seed.canonical_kmer ^ 1,
    };
    let error = reader
        .seed_occurrences_at_scale(collision_query, 1)
        .unwrap_err();
    assert!(error.to_string().contains("hash does not match"));
}
