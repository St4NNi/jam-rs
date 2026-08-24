use serde_json::Value;
use sha2::{Digest, Sha256};
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

fn assert_supported_trace(path: &Path) {
    let records = records(path);
    assert_eq!(records.first().unwrap()["record_type"], "run_header");
    assert_eq!(records.last().unwrap()["record_type"], "run_footer");
    let result = records
        .iter()
        .find(|record| record["record_type"] == "metagenome_result")
        .unwrap();
    match result["candidate"]["score_mode"].as_str() {
        Some("uniform") => assert!(result["candidate"]["uniform_hash_e_value"].is_number()),
        Some("witness") => assert!(result["router_candidate"].is_object()),
        mode => panic!("unexpected candidate score mode: {mode:?}"),
    }
    assert!(
        !result["alignments"].as_array().unwrap().is_empty(),
        "trace result had no alignments: {result}"
    );
    assert!(result["coverage"]["supported_fraction"].as_f64().unwrap() > 0.95);
    assert_eq!(result["performance_counters"]["alignments_attempted"], 0);
    assert!(
        result["alignments"]
            .as_array()
            .unwrap()
            .iter()
            .all(|alignment| {
                alignment["cigar"]
                    .as_str()
                    .is_some_and(|cigar| cigar.ends_with('='))
            })
    );
    assert_eq!(result["status"], "complete");
    assert!(result["failures"].as_array().unwrap().is_empty());
    assert_eq!(result["resource_metrics"]["full_object_fallbacks"], 0);
}

#[test]
fn local_self_contained_jma_trace_is_thread_deterministic() {
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
    let archive_sha256 = format!("{:x}", Sha256::digest(fs::read(&archive).unwrap()));
    fs::write(
        &indexed_catalog,
        format!(
            "metagenome_id\tresource_uri\tsha256\nassembly.fa\tassembly.jma\t{archive_sha256}\n"
        ),
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
    assert_supported_trace(&indexed_output);
    assert!(!directory.path().join("assembly.jma.idx.json").exists());

    let router = directory.path().join("metagenomes.jwr");
    run_ok(
        jam()
            .arg("--silent")
            .arg("router")
            .arg("build")
            .arg("--metagenomes")
            .arg(&indexed_catalog)
            .arg("--output")
            .arg(&router)
            .arg("--k21-base-scale")
            .arg("1")
            .arg("--tiers")
            .arg("1"),
    );
    let router_output = directory.path().join("router.jsonl");
    run_ok(
        jam()
            .arg("--silent")
            .arg("trace")
            .arg("--plasmid")
            .arg(&plasmid)
            .arg("--router")
            .arg(&router)
            .arg("--catalog")
            .arg(&indexed_catalog)
            .arg("--output")
            .arg(&router_output)
            .arg("--sensitivity")
            .arg("sensitive")
            .arg("--top-candidates")
            .arg("1"),
    );
    assert_supported_trace(&router_output);

    let parallel_output = directory.path().join("parallel.jsonl");
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
            .arg(&indexed_catalog)
            .arg("--output")
            .arg(&parallel_output)
            .arg("--sensitivity")
            .arg("sensitive")
            .arg("--min-shared")
            .arg("1")
            .arg("--top-candidates")
            .arg("1"),
    );
    assert_supported_trace(&parallel_output);
    let indexed_result = records(&indexed_output)
        .into_iter()
        .find(|record| record["record_type"] == "metagenome_result")
        .unwrap();
    let parallel_result = records(&parallel_output)
        .into_iter()
        .find(|record| record["record_type"] == "metagenome_result")
        .unwrap();
    let router_result = records(&router_output)
        .into_iter()
        .find(|record| record["record_type"] == "metagenome_result")
        .unwrap();
    for field in ["alignments", "coverage", "failures"] {
        assert_eq!(
            indexed_result[field], router_result[field],
            "router trace changed biological field {field}"
        );
    }
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
            indexed_result[field], parallel_result[field],
            "parallel trace changed biological field {field}"
        );
    }
}
