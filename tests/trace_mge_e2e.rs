use serde_json::{Value, json};
use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use tempfile::TempDir;

struct LocalResources {
    directory: TempDir,
    database: PathBuf,
    metagenomes: PathBuf,
}

fn fixture(file: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/trace_mge")
        .join(file)
}

fn command_ok(command: &mut Command) -> Output {
    let output = command.output().expect("run jam command");
    assert!(
        output.status.success(),
        "command failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn local_resources() -> LocalResources {
    let directory = tempfile::Builder::new()
        .prefix("trace-mge-e2e-")
        .tempdir_in("target")
        .expect("create target-local fixture directory");
    let database = directory.path().join("metagenomes.jam");
    let assemblies = [
        ("assembly_circular.fa", "circular"),
        ("assembly_phages.fa", "phages"),
        ("assembly_common_elements.fa", "common"),
    ];
    let assembly_paths = assemblies
        .iter()
        .map(|(file, _)| fixture(file))
        .collect::<Vec<_>>();

    let mut sketch = Command::new(env!("CARGO_BIN_EXE_jam"));
    sketch
        .arg("--silent")
        .arg("sketch")
        .args(&assembly_paths)
        .arg("--output")
        .arg(&database)
        .arg("--kmer-size")
        .arg("21")
        .arg("--fscale")
        .arg("1");
    command_ok(&mut sketch);

    let mut catalog = String::from("metagenome_id\tjma\tjma_index\traw\n");
    for (file, id) in assemblies {
        let input = fixture(file);
        let archive = directory.path().join(format!("{id}.jma"));
        let mut archive_command = Command::new(env!("CARGO_BIN_EXE_jam"));
        archive_command
            .arg("--silent")
            .arg("archive")
            .arg("--input")
            .arg(&input)
            .arg("--output")
            .arg(&archive)
            .arg("--primary-scale")
            .arg("1")
            .arg("--rescue-scale")
            .arg("1")
            .arg("--block-bases")
            .arg("64");
        command_ok(&mut archive_command);
        let index = format!("{}.idx.json", archive.display());
        catalog.push_str(&format!(
            "{file}\t{}\t{index}\t{}\n",
            archive.display(),
            input.display()
        ));
    }
    let metagenomes = directory.path().join("metagenomes.tsv");
    fs::write(&metagenomes, catalog).expect("write local metagenome catalog");

    LocalResources {
        directory,
        database,
        metagenomes,
    }
}

fn run_trace(
    resources: &LocalResources,
    query: &str,
    query_kind: &str,
    topology: &str,
    threads: usize,
    output_name: &str,
) -> Vec<Value> {
    let output = resources.directory.path().join(output_name);
    let mut command = Command::new(env!("CARGO_BIN_EXE_jam"));
    command
        .arg("--silent")
        .arg("--threads")
        .arg(threads.to_string())
        .arg("trace")
        .arg("--query")
        .arg(fixture(query))
        .arg("--query-kind")
        .arg(query_kind)
        .arg("--topology")
        .arg(topology)
        .arg("--database")
        .arg(&resources.database)
        .arg("--metagenomes")
        .arg(&resources.metagenomes)
        .arg("--output")
        .arg(&output)
        .arg("--sensitivity")
        .arg("sensitive")
        .arg("--min-shared")
        .arg("1")
        .arg("--top-candidates")
        .arg("3")
        .arg("--max-alignments")
        .arg("256")
        .arg("--io-concurrency")
        .arg(threads.max(1).to_string());
    command_ok(&mut command);
    fs::read_to_string(output)
        .expect("read trace JSONL")
        .lines()
        .map(|line| serde_json::from_str(line).expect("valid trace JSONL record"))
        .collect()
}

fn result_for<'a>(records: &'a [Value], metagenome_id: &str) -> &'a Value {
    records
        .iter()
        .find(|record| {
            record["record_type"] == "metagenome_result" && record["metagenome_id"] == metagenome_id
        })
        .unwrap_or_else(|| panic!("missing metagenome result for {metagenome_id}"))
}

fn assert_bounded_coverage(result: &Value, query_length: u64) {
    let coverage = &result["coverage"];
    let supported = coverage["supported_bases"]
        .as_u64()
        .expect("supported_bases is an integer");
    assert!(supported <= query_length);
    assert!(
        coverage["supported_fraction"]
            .as_f64()
            .expect("supported_fraction is numeric")
            <= 1.0
    );
    for interval in coverage["primary_intervals"]
        .as_array()
        .expect("primary intervals is an array")
    {
        assert!(interval["start"].as_u64().unwrap() <= query_length);
        assert!(interval["end"].as_u64().unwrap() <= query_length);
    }
    for alignment in result["alignments"]
        .as_array()
        .expect("alignments is an array")
    {
        for segment in alignment["query_segments"]
            .as_array()
            .expect("query segments is an array")
        {
            assert!(segment["start"].as_u64().unwrap() <= query_length);
            assert!(segment["end"].as_u64().unwrap() <= query_length);
        }
    }
}

fn assert_accepted_alignments_are_retained(result: &Value) {
    let alignments = result["alignments"]
        .as_array()
        .expect("alignments is an array");
    assert!(
        !alignments.is_empty(),
        "expected at least one accepted alignment"
    );
    let ids = alignments
        .iter()
        .map(|alignment| {
            alignment["alignment_id"]
                .as_str()
                .expect("alignment has an ID")
                .to_string()
        })
        .collect::<BTreeSet<_>>();
    assert_eq!(
        ids.len(),
        alignments.len(),
        "alignment IDs are deterministic"
    );

    let mosaic = &result["primary_fragment_mosaic"];
    let accepted = mosaic["accepted_alignment_count"]
        .as_u64()
        .expect("accepted alignment count is an integer");
    assert!(accepted <= alignments.len() as u64);
    for evidence in mosaic["alignment_evidence"]
        .as_array()
        .expect("mosaic alignment evidence is an array")
    {
        assert!(ids.contains(evidence["alignment_id"].as_str().unwrap()));
    }
}

fn biological_result(result: &Value) -> Value {
    json!({
        "query_kind": result["query_kind"],
        "topology_requested": result["topology_requested"],
        "coordinate_model": result["coordinate_model"],
        "topology_evidence": result["topology_evidence"],
        "alignments": result["alignments"],
        "primary_fragment_mosaic": result["primary_fragment_mosaic"],
        "topology": result["topology"],
        "coverage": result["coverage"],
    })
}

#[test]
fn generic_queries_trace_local_jma_with_topology_and_thread_invariants() {
    let resources = local_resources();

    let phage_single_records = run_trace(
        &resources,
        "query_linear_phage.fa",
        "phage",
        "linear",
        1,
        "phage-single.jsonl",
    );
    let phage_single = result_for(&phage_single_records, "assembly_phages.fa");
    assert_eq!(phage_single["query_kind"], "phage");
    assert_eq!(phage_single["topology_requested"], "linear");
    assert_eq!(phage_single["coordinate_model"], "linear");
    assert_eq!(phage_single["topology_evidence"], "linear_supported");
    assert_accepted_alignments_are_retained(phage_single);
    assert_bounded_coverage(phage_single, 256);
    for alignment in phage_single["alignments"].as_array().unwrap() {
        assert!(!alignment["origin_crossing"].as_bool().unwrap());
        assert_eq!(alignment["query_segments"].as_array().unwrap().len(), 1);
    }

    let phage_parallel_records = run_trace(
        &resources,
        "query_linear_phage.fa",
        "phage",
        "linear",
        4,
        "phage-parallel.jsonl",
    );
    let phage_parallel = result_for(&phage_parallel_records, "assembly_phages.fa");
    assert_eq!(
        biological_result(phage_single),
        biological_result(phage_parallel),
        "thread count must not change biological records"
    );

    let circular_records = run_trace(
        &resources,
        "query_circular_plasmid_rotated.fa",
        "plasmid",
        "circular",
        1,
        "circular.jsonl",
    );
    let circular = result_for(&circular_records, "assembly_circular.fa");
    assert_eq!(circular["topology_requested"], "circular");
    assert_eq!(circular["coordinate_model"], "wrap");
    assert_eq!(circular["topology_evidence"], "wrap_supported");
    assert_accepted_alignments_are_retained(circular);
    assert_bounded_coverage(circular, 240);
    assert_eq!(
        circular["topology"]["wrap_model"]["coordinate_model"],
        "wrap"
    );

    let auto_records = run_trace(
        &resources,
        "query_circular_plasmid_auto.fa",
        "plasmid",
        "auto",
        1,
        "auto.jsonl",
    );
    let auto = result_for(&auto_records, "assembly_circular.fa");
    assert_eq!(auto["topology_requested"], "auto");
    assert!(auto["topology"]["linear_model"].is_object());
    assert!(auto["topology"]["wrap_model"].is_object());
    assert_accepted_alignments_are_retained(auto);
    assert_bounded_coverage(auto, 240);

    let unknown_records = run_trace(
        &resources,
        "query_unknown_element.fa",
        "unknown",
        "unknown",
        1,
        "unknown.jsonl",
    );
    let unknown = result_for(&unknown_records, "assembly_common_elements.fa");
    assert_eq!(unknown["query_kind"], "unknown");
    assert_eq!(unknown["topology_requested"], "unknown");
    assert_eq!(unknown["coordinate_model"], "undetermined");
    assert_eq!(unknown["topology_evidence"], "undetermined");
    assert!(
        unknown["warnings"]
            .as_array()
            .unwrap()
            .iter()
            .any(|warning| warning
                .as_str()
                .unwrap()
                .contains("wrapped coordinates do not establish biological topology"))
    );
    assert_accepted_alignments_are_retained(unknown);
    assert_bounded_coverage(unknown, 180);
}
