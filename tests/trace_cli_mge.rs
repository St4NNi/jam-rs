use serde_json::Value;
use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn jam() -> Command {
    Command::new(env!("CARGO_BIN_EXE_jam"))
}

fn fixture(file: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/trace_mge")
        .join(file)
}

fn run_ok(command: &mut Command) -> Output {
    let output = command.output().expect("run jam command");
    assert!(
        output.status.success(),
        "command failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn directory_entries(directory: &Path) -> BTreeSet<String> {
    fs::read_dir(directory)
        .expect("read target-local trace directory")
        .map(|entry| {
            entry
                .expect("read target-local trace entry")
                .file_name()
                .to_string_lossy()
                .into_owned()
        })
        .collect()
}

fn write_catalog(directory: &Path, raw: &Path) -> PathBuf {
    let catalog = directory.join("metagenomes.tsv");
    fs::write(
        &catalog,
        format!("metagenome_id\traw\nassembly.fa\t{}\n", raw.display()),
    )
    .expect("write metagenome catalog");
    catalog
}

fn trace_with_missing_database(
    query: &Path,
    catalog: &Path,
    output: &Path,
    query_kind: &str,
    topology: &str,
) -> Output {
    jam()
        .args(["--silent", "trace", "--query"])
        .arg(query)
        .args(["--query-kind", query_kind, "--topology", topology])
        .args(["--database"])
        .arg(query.with_file_name("missing.jam"))
        .args(["--metagenomes"])
        .arg(catalog)
        .args(["--output"])
        .arg(output)
        .output()
        .expect("run trace parser contract command")
}

#[test]
fn query_primary_accepts_all_query_kinds_and_topologies() {
    let directory = tempfile::Builder::new()
        .prefix("trace-cli-contract-")
        .tempdir_in("target")
        .expect("create target-local tempdir");
    let query = fixture("query_unknown_element.fa");
    let catalog = write_catalog(directory.path(), &fixture("assembly_circular.fa"));

    for query_kind in ["plasmid", "phage", "other", "unknown"] {
        for topology in ["linear", "circular", "auto", "unknown"] {
            let output = directory
                .path()
                .join(format!("{query_kind}-{topology}.jsonl"));
            let result =
                trace_with_missing_database(&query, &catalog, &output, query_kind, topology);
            assert!(!result.status.success());
            let stderr = String::from_utf8_lossy(&result.stderr);
            assert!(
                !stderr.contains("invalid value") && !stderr.contains("possible values"),
                "CLI rejected supported query contract value {query_kind}/{topology}: {stderr}"
            );
        }
    }
}

#[test]
fn plasmid_alias_remains_compatible_and_emits_deprecation_notice() {
    let directory = tempfile::Builder::new()
        .prefix("trace-cli-alias-")
        .tempdir_in("target")
        .expect("create target-local tempdir");
    let catalog = write_catalog(directory.path(), &fixture("assembly_circular.fa"));
    let output = directory.path().join("alias.jsonl");
    let command_output = jam()
        .args(["--silent", "trace", "--plasmid"])
        .arg(fixture("query_circular_plasmid.fa"))
        .args(["--database"])
        .arg(directory.path().join("missing.jam"))
        .args(["--metagenomes"])
        .arg(&catalog)
        .args(["--output"])
        .arg(&output)
        .output()
        .expect("run plasmid compatibility command");
    assert!(!command_output.status.success());
    let stderr = String::from_utf8_lossy(&command_output.stderr);
    assert!(
        stderr.contains("--plasmid is deprecated"),
        "stderr: {stderr}"
    );
    assert!(stderr.contains("--query-kind plasmid"), "stderr: {stderr}");
    assert!(!stderr.contains("invalid value"), "stderr: {stderr}");
}

#[test]
fn zero_and_multiple_query_records_are_rejected_before_search() {
    let directory = tempfile::Builder::new()
        .prefix("trace-cli-query-records-")
        .tempdir_in("target")
        .expect("create target-local tempdir");
    let catalog = write_catalog(directory.path(), &fixture("assembly_circular.fa"));

    let empty_query = directory.path().join("empty.fa");
    fs::write(&empty_query, "").expect("write empty query");
    let empty_output = directory.path().join("empty.jsonl");
    let empty_result =
        trace_with_missing_database(&empty_query, &catalog, &empty_output, "unknown", "auto");
    assert!(!empty_result.status.success());
    assert!(
        String::from_utf8_lossy(&empty_result.stderr)
            .contains("exactly one FASTA/FASTQ record, got 0")
    );
    assert!(!empty_output.exists());

    let multiple_query = directory.path().join("multiple.fa");
    fs::write(
        &multiple_query,
        ">first\nACGTACGTACGTACGTACGTACGTACGTACGT\n>second\nTGCATGCATGCATGCATGCATGCATGCATGCA\n",
    )
    .expect("write multiple-record query");
    let multiple_output = directory.path().join("multiple.jsonl");
    let multiple_result = trace_with_missing_database(
        &multiple_query,
        &catalog,
        &multiple_output,
        "unknown",
        "auto",
    );
    assert!(!multiple_result.status.success());
    assert!(
        String::from_utf8_lossy(&multiple_result.stderr)
            .contains("exactly one FASTA/FASTQ record, got 2")
    );
    assert!(!multiple_output.exists());
}

#[test]
fn query_trace_writes_one_v2_jsonl_stream_without_default_side_files() {
    let directory = tempfile::Builder::new()
        .prefix("trace-cli-output-")
        .tempdir_in("target")
        .expect("create target-local tempdir");
    let query = directory.path().join("query.fa");
    let assembly = directory.path().join("assembly.fa");
    fs::copy(fixture("query_circular_plasmid.fa"), &query).expect("copy query fixture");
    fs::copy(fixture("assembly_circular.fa"), &assembly).expect("copy assembly fixture");
    let database = directory.path().join("metagenomes.jam");
    run_ok(
        jam()
            .args(["--silent", "sketch"])
            .arg(&assembly)
            .args(["--output"])
            .arg(&database)
            .args(["--kmer-size", "21", "--fscale", "1"]),
    );
    let catalog = write_catalog(directory.path(), &assembly);
    let output = directory.path().join("result.jsonl");
    let before_trace = directory_entries(directory.path());

    run_ok(
        jam()
            .args(["--silent", "--threads", "1", "trace", "--query"])
            .arg(&query)
            .args([
                "--query-kind",
                "plasmid",
                "--topology",
                "circular",
                "--database",
            ])
            .arg(&database)
            .args(["--metagenomes"])
            .arg(&catalog)
            .args(["--output"])
            .arg(&output)
            .args([
                "--sensitivity",
                "sensitive",
                "--min-shared",
                "1",
                "--top-candidates",
                "1",
                "--io-concurrency",
                "1",
            ]),
    );

    let records: Vec<Value> = fs::read_to_string(&output)
        .expect("read trace JSONL output")
        .lines()
        .map(|line| serde_json::from_str(line).expect("parse trace JSONL record"))
        .collect();
    assert!(records.len() >= 3, "expected header, result, and footer");
    assert_eq!(records.first().unwrap()["record_type"], "run_header");
    assert_eq!(records.last().unwrap()["record_type"], "run_footer");
    assert_eq!(records.first().unwrap()["schema_version"], "2.0.0");
    assert_eq!(records.last().unwrap()["schema_version"], "2.0.0");
    let result_records: Vec<&Value> = records
        .iter()
        .filter(|record| record["record_type"] == "metagenome_result")
        .collect();
    assert_eq!(result_records.len(), 1);
    assert_eq!(result_records[0]["schema_version"], "2.0.0");
    assert_eq!(result_records[0]["query_kind"], "plasmid");
    assert_eq!(result_records[0]["topology_requested"], "circular");
    assert!(result_records[0].get("query_id").is_some());
    assert!(result_records[0].get("plasmid_id").is_none());

    let mut expected_entries = before_trace;
    expected_entries.insert("result.jsonl".to_string());
    assert_eq!(directory_entries(directory.path()), expected_entries);
    assert!(!directory.path().join("result.jsonl.zst").exists());
    assert!(!directory.path().join("result.jsonl.json").exists());
}
