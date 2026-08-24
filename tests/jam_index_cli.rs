use jam_rs::jam_index::load_manifest;
use jam_rs::provenance::sha256_file;
use std::fs;
use std::process::Command;

fn catalog(
    directory: &std::path::Path,
    id: &str,
    base: u8,
) -> (std::path::PathBuf, std::path::PathBuf) {
    let source = directory.join(format!("{id}.fasta"));
    fs::write(
        &source,
        format!(">contig\n{}\n", char::from(base).to_string().repeat(200)),
    )
    .unwrap();
    let catalog = directory.join(format!("{id}.tsv"));
    fs::write(
        &catalog,
        format!(
            "metagenome_id\tresource_uri\tsha256\n{id}\t{}\t{}\n",
            source.display(),
            sha256_file(&source).unwrap()
        ),
    )
    .unwrap();
    (source, catalog)
}

#[test]
fn index_builds_parts() {
    let directory = tempfile::Builder::new()
        .prefix("jam-index-cli-")
        .tempdir_in("target")
        .unwrap();
    let (_first, first_catalog) = catalog(directory.path(), "mg-a", b'A');
    let index = directory.path().join("index");
    let build = Command::new(env!("CARGO_BIN_EXE_jam"))
        .args(["--silent", "index", "build", "--metagenomes"])
        .arg(&first_catalog)
        .args(["--output"])
        .arg(&index)
        .args(["--max-part-bases", "1000", "--parallel-parts", "2"])
        .output()
        .unwrap();
    assert!(
        build.status.success(),
        "{}",
        String::from_utf8_lossy(&build.stderr)
    );
    assert_eq!(load_manifest(&index).unwrap().parts.len(), 1);

    let (_second, second_catalog) = catalog(directory.path(), "mg-b", b'C');
    let append = Command::new(env!("CARGO_BIN_EXE_jam"))
        .args(["--silent", "index", "append", "--metagenomes"])
        .arg(&second_catalog)
        .args(["--output"])
        .arg(&index)
        .output()
        .unwrap();
    assert!(
        append.status.success(),
        "{}",
        String::from_utf8_lossy(&append.stderr)
    );
    let manifest = load_manifest(&index).unwrap();
    assert_eq!(manifest.parts.len(), 2);
    assert_eq!(manifest.total_metagenomes, 2);
}
