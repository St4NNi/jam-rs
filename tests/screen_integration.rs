use jam_rs::provenance::{self, DatabaseManifest, HASH_ID, HASH_ZERO_POLICY};
use jam_rs::screen::{CONTIG_HEADER, SUMMARY_HEADER};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

struct SyntheticData {
    catalog: PathBuf,
    assembly: PathBuf,
    no_hit: PathBuf,
    empty: PathBuf,
    malformed: PathBuf,
}

struct Tsv {
    header: Vec<String>,
    rows: Vec<Vec<String>>,
}

impl Tsv {
    fn read(path: &Path) -> Self {
        let text = fs::read_to_string(path).unwrap();
        let mut lines = text.lines();
        let header = lines
            .next()
            .unwrap()
            .split('\t')
            .map(str::to_string)
            .collect();
        let rows = lines
            .map(|line| line.split('\t').map(str::to_string).collect())
            .collect();
        Self { header, rows }
    }

    fn index(&self, name: &str) -> usize {
        self.header.iter().position(|field| field == name).unwrap()
    }

    fn value<'a>(&self, row: &'a [String], name: &str) -> &'a str {
        &row[self.index(name)]
    }

    fn row(&self, query: &str, reference: &str) -> &[String] {
        let query_index = self.index("query_contig");
        let reference_index = self.index("reference");
        self.rows
            .iter()
            .find(|row| row[query_index] == query && row[reference_index] == reference)
            .unwrap_or_else(|| panic!("missing row {query} x {reference}"))
    }

    fn summary_row(&self, reference: &str) -> &[String] {
        let reference_index = self.index("reference");
        self.rows
            .iter()
            .find(|row| row[reference_index] == reference)
            .unwrap_or_else(|| panic!("missing summary row for {reference}"))
    }
}

fn jam() -> Command {
    Command::new(env!("CARGO_BIN_EXE_jam"))
}

fn run(command: &mut Command) -> Output {
    command.output().unwrap()
}

fn run_ok(command: &mut Command) -> Output {
    let output = run(command);
    assert!(
        output.status.success(),
        "command failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn deterministic_dna(seed: u64, len: usize) -> String {
    let mut state = seed;
    let bases = *b"ACGT";
    (0..len)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            bases[(state >> 62) as usize] as char
        })
        .collect()
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => unreachable!(),
        })
        .collect()
}

fn fasta_record(name: &str, sequence: &str) -> String {
    format!(">{name}\n{sequence}\n")
}

fn write_synthetic_data(directory: &Path) -> SyntheticData {
    let mobile = deterministic_dna(91, 80);
    let plasmid_a = deterministic_dna(1, 280);
    let plasmid_b = format!(
        "{}{}{}",
        deterministic_dna(2, 100),
        mobile,
        deterministic_dna(3, 100)
    );
    let plasmid_c = format!(
        "{}{}{}",
        deterministic_dna(4, 80),
        mobile,
        deterministic_dna(5, 120)
    );

    let catalog = directory.join("catalog.fa");
    fs::write(
        &catalog,
        [
            fasta_record("plasmid_a", &plasmid_a),
            fasta_record("plasmid_b", &plasmid_b),
            fasta_record("plasmid_c", &plasmid_c),
        ]
        .concat(),
    )
    .unwrap();

    let assembly = directory.join("assembly.fa");
    fs::write(
        &assembly,
        [
            fasta_record("exact_a", &plasmid_a),
            fasta_record("reverse_b", &reverse_complement(&plasmid_b)),
            fasta_record("partial_a", &plasmid_a[70..180]),
            fasta_record("split_c_1", &plasmid_c[0..130]),
            fasta_record("split_c_2", &plasmid_c[100..220]),
            fasta_record("split_c_3", &plasmid_c[190..280]),
            fasta_record("ambiguous_mobile", &mobile),
            fasta_record("chromosome_like", &deterministic_dna(8_888, 360)),
        ]
        .concat(),
    )
    .unwrap();

    let no_hit = directory.join("no_hit.fa");
    fs::write(
        &no_hit,
        fasta_record("unrelated", &deterministic_dna(99_999, 360)),
    )
    .unwrap();
    let empty = directory.join("empty.fa");
    fs::write(&empty, "").unwrap();
    let malformed = directory.join("malformed.fa");
    fs::write(&malformed, "this is not FASTA\nACGT\n").unwrap();
    fs::write(
        directory.join("truth.tsv"),
        "query_contig\treference\texpectation\nexact_a\tplasmid_a\texact\nreverse_b\tplasmid_b\treverse_complement\npartial_a\tplasmid_a\tpartial\nsplit_c_1,split_c_2,split_c_3\tplasmid_c\tsplit_union\nambiguous_mobile\tplasmid_b,plasmid_c\tambiguous_repeat\nchromosome_like\tNA\tnegative\n",
    )
    .unwrap();

    SyntheticData {
        catalog,
        assembly,
        no_hit,
        empty,
        malformed,
    }
}

fn build_database(catalog: &Path, database: &Path, kmer_size: u8, bias: Option<&Path>) {
    let mut command = jam();
    command
        .args(["--threads", "2", "--memory", "1", "--silent", "--force"])
        .arg("sketch")
        .arg(catalog)
        .arg("--output")
        .arg(database)
        .arg("--kmer-size")
        .arg(kmer_size.to_string())
        .arg("--singleton");
    if let Some(bias) = bias {
        command.arg("--bias-table").arg(bias);
    } else {
        command.args(["--fscale", "1"]);
    }
    run_ok(&mut command);
}

#[allow(clippy::too_many_arguments)]
fn screen(
    input: &Path,
    database: &Path,
    output: &Path,
    summary: &Path,
    metadata: Option<&Path>,
    threads: usize,
    top_per_contig: usize,
    min_shared: usize,
) -> Output {
    let mut command = jam();
    command
        .arg("--threads")
        .arg(threads.to_string())
        .args(["--memory", "1", "--silent", "--force"])
        .arg("screen")
        .arg("--input")
        .arg(input)
        .arg("--database")
        .arg(database)
        .arg("--output")
        .arg(output)
        .arg("--summary")
        .arg(summary)
        .arg("--assembly-name")
        .arg("synthetic")
        .arg("--min-shared")
        .arg(min_shared.to_string())
        .arg("--top-per-contig")
        .arg(top_per_contig.to_string())
        .args(["--top-references", "100"]);
    if let Some(metadata) = metadata {
        command.arg("--metadata").arg(metadata);
    }
    run_ok(&mut command)
}

#[test]
fn synthetic_truth_and_union_containment_are_reported() {
    let directory = tempfile::tempdir().unwrap();
    let data = write_synthetic_data(directory.path());
    let database = directory.path().join("catalog.jam");
    build_database(&data.catalog, &database, 21, None);

    let contigs_1 = directory.path().join("contigs-1.tsv");
    let summary_1 = directory.path().join("summary-1.tsv");
    let metadata = directory.path().join("run.json");
    screen(
        &data.assembly,
        &database,
        &contigs_1,
        &summary_1,
        Some(&metadata),
        1,
        10,
        1,
    );
    let contigs_4 = directory.path().join("contigs-4.tsv");
    let summary_4 = directory.path().join("summary-4.tsv");
    screen(
        &data.assembly,
        &database,
        &contigs_4,
        &summary_4,
        None,
        4,
        10,
        1,
    );
    assert_eq!(fs::read(&contigs_1).unwrap(), fs::read(&contigs_4).unwrap());
    assert_eq!(fs::read(&summary_1).unwrap(), fs::read(&summary_4).unwrap());

    let contigs = Tsv::read(&contigs_1);
    let summary = Tsv::read(&summary_1);
    assert_eq!(contigs.header.join("\t"), CONTIG_HEADER);
    assert_eq!(summary.header.join("\t"), SUMMARY_HEADER);

    for (query, reference) in [("exact_a", "plasmid_a"), ("reverse_b", "plasmid_b")] {
        let row = contigs.row(query, reference);
        assert_eq!(contigs.value(row, "query_containment"), "1.000000");
        assert_eq!(contigs.value(row, "reference_containment"), "1.000000");
    }

    let partial = contigs.row("partial_a", "plasmid_a");
    assert_eq!(contigs.value(partial, "query_containment"), "1.000000");
    assert!(
        contigs
            .value(partial, "reference_containment")
            .parse::<f64>()
            .unwrap()
            < 1.0
    );

    contigs.row("ambiguous_mobile", "plasmid_b");
    contigs.row("ambiguous_mobile", "plasmid_c");
    let query_index = contigs.index("query_contig");
    assert!(
        contigs
            .rows
            .iter()
            .all(|row| row[query_index] != "chromosome_like")
    );

    let split = summary.summary_row("plasmid_c");
    let union = summary
        .value(split, "shared_hashes_union")
        .parse::<u64>()
        .unwrap();
    let reference_hashes = summary
        .value(split, "reference_hashes")
        .parse::<u64>()
        .unwrap();
    assert_eq!(
        union, reference_hashes,
        "overlapping contigs inflated union"
    );
    assert_eq!(
        summary.value(split, "aggregate_reference_containment"),
        "1.000000"
    );
    assert!(
        summary
            .value(split, "supporting_contigs")
            .parse::<usize>()
            .unwrap()
            >= 3
    );

    let metadata: serde_json::Value = serde_json::from_slice(&fs::read(metadata).unwrap()).unwrap();
    assert_eq!(metadata["schema_version"], "1.0.0");
    assert_eq!(metadata["score_mode"], "uniform");
    assert_eq!(metadata["results"]["contig_count"], 8);

    let filtered_output = directory.path().join("filtered.tsv");
    let filtered_summary = directory.path().join("filtered-summary.tsv");
    let mut filtered_command = jam();
    filtered_command
        .args(["--silent", "--force", "screen", "--input"])
        .arg(&data.assembly)
        .arg("--database")
        .arg(&database)
        .arg("--output")
        .arg(&filtered_output)
        .arg("--summary")
        .arg(&filtered_summary)
        .args([
            "--min-shared",
            "1",
            "--min-query-containment",
            "0.9",
            "--min-reference-containment",
            "0.9",
        ]);
    run_ok(&mut filtered_command);
    let filtered = Tsv::read(&filtered_output);
    filtered.row("exact_a", "plasmid_a");
    assert!(
        filtered
            .rows
            .iter()
            .all(|row| filtered.value(row, "query_contig") != "partial_a")
    );
}

#[test]
fn top_k_empty_no_hit_and_malformed_inputs_are_stable() {
    let directory = tempfile::tempdir().unwrap();
    let data = write_synthetic_data(directory.path());
    let database = directory.path().join("catalog.jam");
    build_database(&data.catalog, &database, 21, None);

    let top_one = directory.path().join("top-one.tsv");
    let top_one_summary = directory.path().join("top-one-summary.tsv");
    screen(
        &data.assembly,
        &database,
        &top_one,
        &top_one_summary,
        None,
        1,
        1,
        1,
    );
    let first = fs::read(&top_one).unwrap();
    screen(
        &data.assembly,
        &database,
        &top_one,
        &top_one_summary,
        None,
        4,
        1,
        1,
    );
    assert_eq!(first, fs::read(&top_one).unwrap());
    let top = Tsv::read(&top_one);
    let query_index = top.index("query_contig");
    assert_eq!(
        top.rows
            .iter()
            .filter(|row| row[query_index] == "ambiguous_mobile")
            .count(),
        1
    );

    for (name, input) in [("empty", &data.empty), ("no-hit", &data.no_hit)] {
        let output = directory.path().join(format!("{name}.tsv"));
        let summary = directory.path().join(format!("{name}-summary.tsv"));
        screen(input, &database, &output, &summary, None, 2, 10, 1);
        assert_eq!(fs::read_to_string(output).unwrap().trim(), CONTIG_HEADER);
        assert_eq!(fs::read_to_string(summary).unwrap().trim(), SUMMARY_HEADER);
    }

    let malformed_output = directory.path().join("malformed.tsv");
    let malformed_summary = directory.path().join("malformed-summary.tsv");
    let mut command = jam();
    command
        .args(["--silent", "screen", "--input"])
        .arg(&data.malformed)
        .arg("--database")
        .arg(&database)
        .arg("--output")
        .arg(&malformed_output)
        .arg("--summary")
        .arg(&malformed_summary);
    let failure = run(&mut command);
    assert!(!failure.status.success());
    assert!(String::from_utf8_lossy(&failure.stderr).contains("Parse error"));
    assert!(!malformed_output.exists());
    assert!(!malformed_summary.exists());
}

#[test]
fn manifest_stats_legacy_fixture_and_hash_zero_are_compatible() {
    let directory = tempfile::tempdir().unwrap();
    let data = write_synthetic_data(directory.path());
    let database = directory.path().join("catalog.jam");
    build_database(&data.catalog, &database, 21, None);

    let manifest_path = provenance::sidecar_path(&database);
    let manifest: DatabaseManifest =
        serde_json::from_slice(&fs::read(&manifest_path).unwrap()).unwrap();
    assert_eq!(manifest.database_format_version, 3);
    assert_eq!(manifest.hash_id, HASH_ID);
    assert_eq!(manifest.hash_zero_policy, HASH_ZERO_POLICY);
    assert_eq!(manifest.input_catalog_files.len(), 1);
    assert_eq!(manifest.sample_count, 3);

    let mut stats_command = jam();
    stats_command
        .args(["--silent", "stats", "--input"])
        .arg(&database)
        .arg("--json");
    let stats_output = run_ok(&mut stats_command);
    let stats: serde_json::Value = serde_json::from_slice(&stats_output.stdout).unwrap();
    assert_eq!(stats["database_format_version"], 3);
    assert_eq!(stats["hash_id"], HASH_ID);
    assert_eq!(stats["sample_count"], 3);

    let legacy = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/legacy_v3_0_9_11.jam");
    let legacy_input = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/testfiles/short.fa");
    let legacy_output = directory.path().join("legacy.tsv");
    let legacy_summary = directory.path().join("legacy-summary.tsv");
    screen(
        &legacy_input,
        &legacy,
        &legacy_output,
        &legacy_summary,
        None,
        2,
        10,
        1,
    );
    assert_eq!(Tsv::read(&legacy_output).rows.len(), 1);

    let zero_fasta = directory.path().join("zero.fa");
    fs::write(&zero_fasta, fasta_record("all_a", &"A".repeat(128))).unwrap();
    let zero_database = directory.path().join("zero.jam");
    build_database(&zero_fasta, &zero_database, 21, None);
    let mut zero_stats_command = jam();
    zero_stats_command
        .args(["--silent", "stats", "--input"])
        .arg(&zero_database)
        .arg("--json");
    let zero_stats_output = run_ok(&mut zero_stats_command);
    let zero_stats: serde_json::Value = serde_json::from_slice(&zero_stats_output.stdout).unwrap();
    assert_eq!(zero_stats["entry_count"], 0);
    assert_eq!(zero_stats["unique_hash_count"], 0);

    let zero_output = directory.path().join("zero.tsv");
    let zero_summary = directory.path().join("zero-summary.tsv");
    screen(
        &zero_fasta,
        &zero_database,
        &zero_output,
        &zero_summary,
        None,
        1,
        10,
        1,
    );
    assert!(Tsv::read(&zero_output).rows.is_empty());

    let vectors = [
        (0, 0),
        (1, 0x2d67_3a90_c3c1_4c02),
        (2, 0x2990_8a8b_0e98_8342),
        (3, 0xe40d_1ef9_0d53_5f7b),
        (0x1234_5678_90ab_cdef, 0x8da8_0ccc_204d_8d51),
        (u64::MAX, u64::MAX),
    ];
    for (input, expected) in vectors {
        assert_eq!(jam_rs::jamhash_u64_v1(input), expected);
    }
}

#[test]
fn bias_labels_e_values_and_parameter_mismatches_are_explicit() {
    let directory = tempfile::tempdir().unwrap();
    let data = write_synthetic_data(directory.path());
    let bias = directory.path().join("catalog.bias");
    let mut bias_command = jam();
    bias_command
        .args(["--threads", "2", "--silent", "--force", "bias", "create"])
        .arg("--positive")
        .arg(&data.catalog)
        .arg("--negative")
        .arg(&data.no_hit)
        .arg("--output")
        .arg(&bias)
        .args([
            "--kmer-size",
            "21",
            "--fscale",
            "1",
            "--cms-width",
            "1024",
            "--cms-depth",
            "3",
            "--min-positive-retention",
            "0",
        ]);
    run_ok(&mut bias_command);
    assert!(provenance::sidecar_path(&bias).is_file());

    let bias_database = directory.path().join("bias.jam");
    build_database(&data.catalog, &bias_database, 21, Some(&bias));
    let bias_output = directory.path().join("bias.tsv");
    let bias_summary = directory.path().join("bias-summary.tsv");
    screen(
        &data.assembly,
        &bias_database,
        &bias_output,
        &bias_summary,
        None,
        2,
        10,
        1,
    );
    let bias_rows = Tsv::read(&bias_output);
    assert!(!bias_rows.rows.is_empty());
    for row in &bias_rows.rows {
        assert_eq!(bias_rows.value(row, "query_containment"), "NA");
        assert_eq!(bias_rows.value(row, "reference_containment"), "NA");
        assert_ne!(bias_rows.value(row, "retained_query_containment"), "NA");
        assert_eq!(bias_rows.value(row, "uniform_hash_e_value"), "NA");
        assert_eq!(bias_rows.value(row, "score_mode"), "bias");
        assert!(bias_rows.value(row, "bias_table_id").starts_with("sha256:"));
    }

    let bias_dist = directory.path().join("bias-dist.tsv");
    let mut bias_dist_command = jam();
    bias_dist_command
        .args(["--silent", "--force", "dist", "--input"])
        .arg(&bias_database)
        .arg("--database")
        .arg(&bias_database)
        .arg("--output")
        .arg(&bias_dist);
    run_ok(&mut bias_dist_command);
    let dist = Tsv::read(&bias_dist);
    assert!(
        dist.header
            .contains(&"retained_query_containment".to_string())
    );
    assert!(
        dist.rows
            .iter()
            .all(|row| dist.value(row, "uniform_hash_e_value") == "NA")
    );

    let uniform_21 = directory.path().join("uniform-21.jam");
    let uniform_19 = directory.path().join("uniform-19.jam");
    build_database(&data.catalog, &uniform_21, 21, None);
    build_database(&data.catalog, &uniform_19, 19, None);

    let mut kmer_mismatch = jam();
    kmer_mismatch
        .args(["--silent", "dist", "--input"])
        .arg(&uniform_19)
        .arg("--database")
        .arg(&uniform_21);
    let output = run(&mut kmer_mismatch);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("k-mer size"));

    let mut bias_mismatch = jam();
    bias_mismatch
        .args(["--silent", "dist", "--input"])
        .arg(&uniform_21)
        .arg("--database")
        .arg(&bias_database);
    let output = run(&mut bias_mismatch);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("bias table"));
}

#[test]
fn bounded_memory_screen_smoke_is_incremental() {
    let directory = tempfile::tempdir().unwrap();
    let data = write_synthetic_data(directory.path());
    let database = directory.path().join("catalog.jam");
    build_database(&data.catalog, &database, 21, None);

    let sequence = deterministic_dna(1, 280);
    let many_contigs: String = (0..128)
        .map(|index| fasta_record(&format!("copy_{index:03}"), &sequence))
        .collect();
    let assembly = directory.path().join("many-contigs.fa");
    fs::write(&assembly, many_contigs).unwrap();
    let output = directory.path().join("many.tsv");
    let summary = directory.path().join("many-summary.tsv");
    let metadata = directory.path().join("many.json");
    screen(
        &assembly,
        &database,
        &output,
        &summary,
        Some(&metadata),
        4,
        1,
        1,
    );
    assert_eq!(Tsv::read(&output).rows.len(), 128);
    let metadata: serde_json::Value = serde_json::from_slice(&fs::read(metadata).unwrap()).unwrap();
    assert_eq!(metadata["parameters"]["memory_gb"], 1);
    assert_eq!(metadata["results"]["contig_count"], 128);
}
