use jam_rs::trace::catalog::{CatalogEntry, CatalogError, TraceCatalog};
use jam_rs::trace::raw::parse_stream;
use jam_rs::trace::runner::{TraceRunner, TraceRunnerConfig};
use std::io::Cursor;
use tempfile::tempdir;

fn entry(metagenome_id: &str, resource_uri: &str) -> CatalogEntry {
    CatalogEntry {
        metagenome_id: metagenome_id.to_string(),
        resource_uri: resource_uri.to_string(),
        sha256: "00".repeat(32),
        etag: None,
        object_version: None,
        label: None,
        source_assembly_sha256: None,
    }
}

#[test]
fn catalog_tsv_resolves_relative_resources_and_preserves_order() {
    let directory = tempdir().unwrap();
    let catalog_path = directory.path().join("catalog.tsv");
    std::fs::write(
        &catalog_path,
        concat!(
            "metagenome_id\tresource_uri\tsha256\n",
            "mg-b\tarchives/b.jma\t0000000000000000000000000000000000000000000000000000000000000000\n",
            "mg-a\tarchives/a.jma\t1111111111111111111111111111111111111111111111111111111111111111\n",
        ),
    )
    .unwrap();

    let catalog = TraceCatalog::from_path(&catalog_path).unwrap();
    assert_eq!(catalog.len(), 2);
    assert_eq!(catalog.entries()[0].metagenome_id, "mg-b");
    assert!(
        catalog.entries()[0]
            .resource_uri
            .ends_with("archives/b.jma")
    );
    assert_eq!(catalog.get("mg-a").unwrap().sha256, "11".repeat(32));
}

#[test]
fn catalog_accepts_optional_object_identity_fields() {
    let directory = tempdir().unwrap();
    let catalog_path = directory.path().join("catalog.tsv");
    std::fs::write(
        &catalog_path,
        concat!(
            "metagenome_id\tresource_uri\tsha256\tetag\tobject_version\tlabel\tsource_assembly_sha256\n",
            "mg\tarchives/mg.jma\taaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\tobject-etag\tv3\tlabel\tbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb\n",
        ),
    )
    .unwrap();

    let catalog = TraceCatalog::from_path(&catalog_path).unwrap();
    let entry = catalog.get("mg").unwrap();
    assert!(entry.resource_uri.ends_with("archives/mg.jma"));
    assert_eq!(entry.etag.as_deref(), Some("object-etag"));
    assert_eq!(entry.object_version.as_deref(), Some("v3"));
}

#[test]
fn catalog_rejects_removed_resource_columns() {
    let directory = tempdir().unwrap();
    let catalog_path = directory.path().join("catalog.tsv");
    std::fs::write(
        &catalog_path,
        "metagenome_id\tresource_uri\tsha256\tjma_index\nmg\tmg.jma\t0000000000000000000000000000000000000000000000000000000000000000\tmg.idx\n",
    )
    .unwrap();
    let result = TraceCatalog::from_path(catalog_path);
    assert!(matches!(result, Err(CatalogError::Invalid { .. })));
}

#[test]
fn catalog_rejects_duplicate_ids_and_missing_resources() {
    let duplicate = TraceCatalog::from_entries(vec![entry("mg", "a.jma"), entry("mg", "b.jma")]);
    assert!(matches!(duplicate, Err(CatalogError::Invalid { .. })));

    let mut invalid = entry("mg", "");
    invalid.sha256 = "not-a-checksum".to_string();
    let missing = TraceCatalog::from_entries(vec![invalid]);
    assert!(matches!(missing, Err(CatalogError::Invalid { .. })));
}

#[test]
fn raw_parser_accepts_empty_and_fasta_streams() {
    assert!(
        parse_stream(Cursor::new(Vec::<u8>::new()), "file:///empty.fa")
            .unwrap()
            .is_empty()
    );
    let contigs = parse_stream(
        Cursor::new(b">contig-a description\nACGTN\n>contig-b\nTTAA\n".to_vec()),
        "file:///assembly.fa",
    )
    .unwrap();
    assert_eq!(contigs.len(), 2);
    assert_eq!(contigs[0].id, "contig-a");
    assert_eq!(contigs[0].sequence, b"ACGTN");
}

#[test]
fn runner_configuration_is_explicit_and_validated() {
    let runner = TraceRunner::new(TraceRunnerConfig::default()).unwrap();
    assert_eq!(runner.config().threads, 1);
    assert!(runner.config().sensitivity.primary.k == 31);
}
