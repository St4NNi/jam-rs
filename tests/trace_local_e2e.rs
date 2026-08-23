use serde_json::Value;
use std::fs;
use std::path::Path;
use std::process::{Command, Output};

fn jam() -> Command {
    Command::new(env!("CARGO_BIN_EXE_jam"))
}

fn run_ok(command: &mut Command) -> Output {
    let output = command.output().unwrap();
    assert!(
        output.status.success(),
        "command failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn deterministic_dna(seed: u64, length: usize) -> String {
    let mut state = seed;
    let bases = *b"ACGT";
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            bases[(state >> 62) as usize] as char
        })
        .collect()
}

fn records(path: &Path) -> Vec<Value> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect()
}

fn assert_supported_trace(path: &Path, expect_fallback: bool) {
    let records = records(path);
    assert_eq!(records.first().unwrap()["record_type"], "run_header");
    assert_eq!(records.last().unwrap()["record_type"], "run_footer");
    let result = records
        .iter()
        .find(|record| record["record_type"] == "metagenome_result")
        .unwrap();
    assert_eq!(result["candidate"]["score_mode"], "uniform");
    assert!(result["candidate"]["uniform_hash_e_value"].is_number());
    assert!(!result["alignments"].as_array().unwrap().is_empty());
    assert!(result["coverage"]["supported_fraction"].as_f64().unwrap() > 0.95);
    if expect_fallback {
        assert_eq!(result["status"], "partial");
        assert_eq!(result["failures"][0]["code"], "jma_full_download_fallback");
        assert_eq!(result["resource_metrics"]["full_object_fallbacks"], 1);
    } else {
        assert_eq!(result["status"], "complete");
        assert!(result["failures"].as_array().unwrap().is_empty());
    }
}

#[test]
fn local_jma_and_raw_trace_paths_emit_coverage_jsonl() {
    let directory = tempfile::Builder::new()
        .prefix("trace-local-")
        .tempdir_in("target")
        .unwrap();
    let plasmid_sequence = deterministic_dna(77, 2_000);
    let plasmid = directory.path().join("plasmid.fa");
    fs::write(&plasmid, format!(">plasmid-one\n{plasmid_sequence}\n")).unwrap();
    let assembly = directory.path().join("assembly.fa");
    fs::write(
        &assembly,
        format!(
            ">exact-plasmid\n{plasmid_sequence}\n>chromosome-like\n{}\n",
            deterministic_dna(9_001, 800)
        ),
    )
    .unwrap();

    let database = directory.path().join("metagenomes.jam");
    run_ok(
        jam()
            .arg("--silent")
            .arg("sketch")
            .arg(&assembly)
            .arg("--output")
            .arg(&database)
            .arg("--kmer-size")
            .arg("21")
            .arg("--fscale")
            .arg("1"),
    );
    let archive = directory.path().join("assembly.jma");
    run_ok(
        jam()
            .arg("--silent")
            .arg("archive")
            .arg("--input")
            .arg(&assembly)
            .arg("--output")
            .arg(&archive)
            .arg("--primary-scale")
            .arg("1")
            .arg("--rescue-scale")
            .arg("1")
            .arg("--block-bases")
            .arg("512"),
    );

    let indexed_catalog = directory.path().join("indexed.tsv");
    fs::write(
        &indexed_catalog,
        "metagenome_id\tjma\traw\nassembly.fa\tassembly.jma\tassembly.fa\n",
    )
    .unwrap();
    let indexed_output = directory.path().join("indexed.jsonl");
    run_ok(
        jam()
            .arg("--silent")
            .arg("trace")
            .arg("--plasmid")
            .arg(&plasmid)
            .arg("--database")
            .arg(&database)
            .arg("--catalog")
            .arg(&indexed_catalog)
            .arg("--output")
            .arg(&indexed_output)
            .arg("--sensitivity")
            .arg("sensitive")
            .arg("--min-shared")
            .arg("1")
            .arg("--top-candidates")
            .arg("1"),
    );
    assert_supported_trace(&indexed_output, true);

    let strict_output = directory.path().join("strict.jsonl");
    run_ok(
        jam()
            .arg("--silent")
            .arg("trace")
            .arg("--plasmid")
            .arg(&plasmid)
            .arg("--database")
            .arg(&database)
            .arg("--catalog")
            .arg(&indexed_catalog)
            .arg("--output")
            .arg(&strict_output)
            .arg("--sensitivity")
            .arg("sensitive")
            .arg("--min-shared")
            .arg("1")
            .arg("--top-candidates")
            .arg("1")
            .arg("--no-full-download-fallback"),
    );
    let strict_records = records(&strict_output);
    let strict_result = strict_records
        .iter()
        .find(|record| record["record_type"] == "metagenome_result")
        .unwrap();
    assert_eq!(strict_result["status"], "failed");
    assert_eq!(strict_result["failures"][0]["code"], "jma_index_required");

    let raw_catalog = directory.path().join("raw.tsv");
    fs::write(
        &raw_catalog,
        "metagenome_id\traw\nassembly.fa\tassembly.fa\n",
    )
    .unwrap();
    let raw_output = directory.path().join("raw.jsonl");
    run_ok(
        jam()
            .arg("--silent")
            .arg("trace")
            .arg("--plasmid")
            .arg(&plasmid)
            .arg("--database")
            .arg(&database)
            .arg("--catalog")
            .arg(&raw_catalog)
            .arg("--output")
            .arg(&raw_output)
            .arg("--sensitivity")
            .arg("sensitive")
            .arg("--min-shared")
            .arg("1")
            .arg("--top-candidates")
            .arg("1"),
    );
    assert_supported_trace(&raw_output, false);

    let raw_parallel_output = directory.path().join("raw-parallel.jsonl");
    run_ok(
        jam()
            .arg("--silent")
            .arg("--threads")
            .arg("4")
            .arg("trace")
            .arg("--plasmid")
            .arg(&plasmid)
            .arg("--database")
            .arg(&database)
            .arg("--catalog")
            .arg(&raw_catalog)
            .arg("--output")
            .arg(&raw_parallel_output)
            .arg("--sensitivity")
            .arg("sensitive")
            .arg("--min-shared")
            .arg("1")
            .arg("--top-candidates")
            .arg("1"),
    );
    let raw_result = records(&raw_output)
        .into_iter()
        .find(|record| record["record_type"] == "metagenome_result")
        .unwrap();
    let raw_parallel_result = records(&raw_parallel_output)
        .into_iter()
        .find(|record| record["record_type"] == "metagenome_result")
        .unwrap();
    for field in [
        "schema_version",
        "plasmid_id",
        "metagenome_id",
        "algorithm",
        "status",
        "candidate",
        "alignments",
        "coverage",
        "failures",
    ] {
        assert_eq!(
            raw_result[field], raw_parallel_result[field],
            "parallel trace changed biological field {field}"
        );
    }
}
