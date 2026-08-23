use std::process::Command;

#[test]
fn trace_help_exposes_generic_query_contract() {
    let output = Command::new(env!("CARGO_BIN_EXE_jam"))
        .args(["trace", "--help"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    for flag in [
        "--query",
        "--query-kind",
        "--topology",
        "--metagenomes",
        "--io-concurrency",
    ] {
        assert!(stdout.contains(flag), "trace help is missing {flag}");
    }
    assert!(stdout.contains("--plasmid"));
}

#[test]
fn plasmid_alias_warns_and_rejects_a_conflicting_kind() {
    let output = Command::new(env!("CARGO_BIN_EXE_jam"))
        .args([
            "trace",
            "--plasmid",
            "missing.fa",
            "--query-kind",
            "phage",
            "--database",
            "missing.jam",
            "--metagenomes",
            "missing.tsv",
            "--output",
            "missing.jsonl",
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("--plasmid implies --query-kind plasmid"));
}

#[test]
fn plasmid_alias_emits_the_one_cycle_deprecation_notice() {
    let output = Command::new(env!("CARGO_BIN_EXE_jam"))
        .args([
            "trace",
            "--plasmid",
            "missing.fa",
            "--database",
            "missing.jam",
            "--metagenomes",
            "missing.tsv",
            "--output",
            "missing.jsonl",
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("--plasmid is deprecated"));
    assert!(stderr.contains("--query-kind plasmid"));
}
