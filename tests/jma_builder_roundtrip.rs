use jam_rs::jma::builder::{
    ArchiveBuildConfig, build_archive_from_fasta, write_archive_from_fasta,
};
use jam_rs::jma::reader::JmaReader;
use jam_rs::jma::{ArchiveReader, JMA_MAGIC};
use jam_rs::resource::ResourceOpenOptions;
use jam_rs::resource::local::LocalResource;
use std::io::Write;
use tempfile::NamedTempFile;

fn fasta(contents: &[u8]) -> NamedTempFile {
    let mut file = NamedTempFile::new().unwrap();
    file.write_all(contents).unwrap();
    file
}

#[test]
fn builder_reports_provenance_and_both_seed_levels() {
    let input =
        fasta(b">contig-a\nACGTTGCATGTCAGTAGGCATCAGTACCGATGCTAGCTAGGCTAACGTTACGATCGATGCA\n");
    let output = build_archive_from_fasta(
        input.path(),
        ArchiveBuildConfig {
            k31_scale: 1,
            k21_scale: Some(1),
            ..ArchiveBuildConfig::default()
        },
    )
    .unwrap();

    assert_eq!(output.stats.contig_count, 1);
    assert!(output.stats.total_bases > 0);
    assert_ne!(output.stats.source_sha256, [0u8; 32]);
    assert_eq!(output.parts.seed_sections.len(), 2);
    assert_eq!(output.parts.seed_sections[0].k, 31);
    assert_eq!(output.parts.seed_sections[1].k, 21);
    assert_eq!(output.parts.contigs[0].name, "contig-a");
}

#[test]
fn writer_emits_jma_magic_and_hash_zero_is_excluded() {
    let input = fasta(b">all-a\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n");
    let output_path = tempfile::NamedTempFile::new().unwrap();
    let stats = write_archive_from_fasta(
        input.path(),
        output_path.path(),
        ArchiveBuildConfig {
            k31_scale: 1,
            k21_scale: Some(1),
            ..ArchiveBuildConfig::default()
        },
    )
    .unwrap();

    assert_eq!(stats.k31_seed_count, 0);
    assert_eq!(stats.k21_seed_count, 0);
    let bytes = std::fs::read(output_path.path()).unwrap();
    assert_eq!(&bytes[..JMA_MAGIC.len()], &JMA_MAGIC);
}

#[test]
fn writer_records_entropy_metadata_used_by_seed_builder() {
    let input =
        fasta(b">contig-a\nACGTTGCATGTCAGTAGGCATCAGTACCGATGCTAGCTAGGCTAACGTTACGATCGATGCA\n");
    let output = tempfile::NamedTempFile::new().unwrap();
    write_archive_from_fasta(
        input.path(),
        output.path(),
        ArchiveBuildConfig {
            k31_scale: 1,
            k21_scale: Some(1),
            min_entropy: Some(1.25),
            ..ArchiveBuildConfig::default()
        },
    )
    .unwrap();
    let resource = LocalResource::from_path(output.path(), ResourceOpenOptions::default()).unwrap();
    let reader = JmaReader::open(resource).unwrap();
    assert_eq!(
        reader.header().algorithm_id.as_deref(),
        Some("jam-seed-chain-align-v1")
    );
    assert_eq!(reader.header().algorithm_version, Some(1));
    assert_eq!(reader.header().min_entropy, Some(1.25));
}
