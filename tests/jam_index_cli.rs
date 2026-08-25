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
    let offset = usize::from(base & 3);
    let mut sequence = (0..200)
        .map(|index| b"ACGT"[(index * 13 + index / 7 + offset) & 3] as char)
        .collect::<String>();
    sequence.replace_range(73..74, "R");
    fs::write(&source, format!(">contig\n{sequence}\n")).unwrap();
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
    let (first, first_catalog) = catalog(directory.path(), "mg-a", b'A');
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

    let (second, second_catalog) = catalog(directory.path(), "mg-b", b'C');
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

    let trace = directory.path().join("trace.jsonl");
    let run = Command::new(env!("CARGO_BIN_EXE_jam"))
        .args(["--silent", "trace", "--query"])
        .arg(&first)
        .args(["--index"])
        .arg(&index)
        .args(["--output"])
        .arg(&trace)
        .args([
            "--sensitivity",
            "sensitive",
            "--min-shared",
            "1",
            "--initial-contigs",
            "1",
        ])
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let records = fs::read_to_string(trace).unwrap();
    assert!(records.contains("\"record_type\":\"metagenome_result\""));
    assert!(records.contains("\"metagenome_id\":\"mg-a\""));
    assert!(records.contains("\"cigar\":\"200=\""));

    let query_a = directory.path().join("query-a.fasta");
    let query_b = directory.path().join("query-b.fasta");
    fs::write(
        &query_a,
        fs::read_to_string(&first)
            .unwrap()
            .replacen(">contig", ">query-a description retained", 1),
    )
    .unwrap();
    fs::write(
        &query_b,
        fs::read_to_string(&second)
            .unwrap()
            .replacen(">contig", ">query-b", 1),
    )
    .unwrap();
    let staged_first = directory.path().join("staged-mg-a.fasta");
    let staged_second = directory.path().join("staged-mg-b.fasta");
    fs::copy(&first, &staged_first).unwrap();
    fs::copy(&second, &staged_second).unwrap();
    let staged_catalog = directory.path().join("staged.tsv");
    fs::write(
        &staged_catalog,
        format!(
            "metagenome_id\tresource_uri\tsha256\nmg-a\t{}\t{}\nmg-b\t{}\t{}\n",
            staged_first.display(),
            sha256_file(&staged_first).unwrap(),
            staged_second.display(),
            sha256_file(&staged_second).unwrap(),
        ),
    )
    .unwrap();
    fs::rename(&first, directory.path().join("unavailable-mg-a.fasta")).unwrap();
    fs::rename(&second, directory.path().join("unavailable-mg-b.fasta")).unwrap();
    let batch = directory.path().join("batch.jsonl.zst");
    let run = Command::new(env!("CARGO_BIN_EXE_jam"))
        .args(["--silent", "trace", "--query"])
        .arg(&query_a)
        .arg(&query_b)
        .args(["--index"])
        .arg(&index)
        .args(["--metagenomes"])
        .arg(&staged_catalog)
        .args(["--output"])
        .arg(&batch)
        .args([
            "--sensitivity",
            "sensitive",
            "--min-shared",
            "1",
            "--initial-contigs",
            "1",
        ])
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let records =
        String::from_utf8(zstd::stream::decode_all(std::fs::File::open(&batch).unwrap()).unwrap())
            .unwrap();
    assert!(records.contains("\"query_id\":\"query-a\""));
    assert!(records.contains("\"query_id\":\"query-b\""));
    for suffix in ["candidates", "work", "status"] {
        let path = directory.path().join(format!("batch.{suffix}.tsv.zst"));
        let decoded = zstd::stream::decode_all(std::fs::File::open(path).unwrap()).unwrap();
        assert!(!decoded.is_empty());
        if suffix == "status" {
            assert!(String::from_utf8_lossy(&decoded).contains("query-a description retained"));
        } else if suffix == "work" {
            let decoded = String::from_utf8_lossy(&decoded);
            assert!(decoded.contains(staged_first.to_string_lossy().as_ref()));
            assert!(decoded.contains(staged_second.to_string_lossy().as_ref()));
        }
    }
    let metrics: serde_json::Value = serde_json::from_reader(std::io::BufReader::new(
        std::fs::File::open(directory.path().join("batch.metrics.json")).unwrap(),
    ))
    .unwrap();
    assert_eq!(metrics["execution"]["metrics"]["queries"], 2);
    assert_eq!(metrics["execution"]["metrics"]["source_open_count"], 2);

    let screen_only = directory.path().join("screen-only.jsonl.zst");
    let run = Command::new(env!("CARGO_BIN_EXE_jam"))
        .args(["--silent", "trace", "--query"])
        .arg(&query_a)
        .arg(&query_b)
        .args(["--index"])
        .arg(&index)
        .args(["--output"])
        .arg(&screen_only)
        .args(["--screen-only", "--min-shared", "1"])
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let records = String::from_utf8(
        zstd::stream::decode_all(std::fs::File::open(&screen_only).unwrap()).unwrap(),
    )
    .unwrap();
    assert!(!records.contains("\"record_type\":\"metagenome_result\""));
    let metrics: serde_json::Value = serde_json::from_reader(std::io::BufReader::new(
        std::fs::File::open(directory.path().join("screen-only.metrics.json")).unwrap(),
    ))
    .unwrap();
    assert_eq!(metrics["screen_only"], true);
    assert_eq!(metrics["execution"]["metrics"]["source_open_count"], 0);
    assert!(metrics["execution"]["metrics"]["candidate_pairs"] != 0);
    let help = Command::new(env!("CARGO_BIN_EXE_jam"))
        .arg("--help")
        .output()
        .unwrap();
    assert!(!String::from_utf8_lossy(&help.stdout).contains("trace-batch"));
}
