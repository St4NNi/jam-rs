use jam_rs::trace::catalog::{CatalogEntry, CatalogError, TraceCatalog};
use jam_rs::trace::raw::parse_stream;
use jam_rs::trace::runner::{TraceRunner, TraceRunnerConfig};
use std::io::Cursor;
use tempfile::tempdir;

#[test]
fn catalog_tsv_resolves_relative_resources_and_preserves_order() {
    let directory = tempdir().unwrap();
    let catalog_path = directory.path().join("catalog.tsv");
    std::fs::write(
        &catalog_path,
        "metagenome_id\tjma\traw\nmg-b\tarchives/b.jma\tassemblies/b.fa\nmg-a\t\tassemblies/a.fa\n",
    )
    .unwrap();

    let catalog = TraceCatalog::from_path(&catalog_path).unwrap();
    assert_eq!(catalog.len(), 2);
    assert_eq!(catalog.entries()[0].metagenome_id, "mg-b");
    assert!(
        catalog.entries()[0]
            .jma
            .as_deref()
            .unwrap()
            .ends_with("archives/b.jma")
    );
    assert!(catalog.entries()[1].jma.is_none());
    assert_eq!(
        catalog
            .get("mg-a")
            .unwrap()
            .raw
            .as_deref()
            .unwrap()
            .rsplit('/')
            .next(),
        Some("a.fa")
    );
}

#[test]
fn catalog_rejects_duplicate_ids_and_missing_resources() {
    let duplicate = TraceCatalog::from_entries(vec![
        CatalogEntry {
            metagenome_id: "mg".to_string(),
            jma: Some("a.jma".to_string()),
            raw: None,
        },
        CatalogEntry {
            metagenome_id: "mg".to_string(),
            jma: None,
            raw: Some("a.fa".to_string()),
        },
    ]);
    assert!(matches!(duplicate, Err(CatalogError::Invalid { .. })));

    let missing = TraceCatalog::from_entries(vec![CatalogEntry {
        metagenome_id: "mg".to_string(),
        jma: None,
        raw: None,
    }]);
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
