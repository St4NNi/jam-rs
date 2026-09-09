use crate::cli::{
    Cli, Commands,
    handlers::{TraceArgs, TraceInput, handle_trace_command},
};
use crate::owner_tests::{build_fixture_at, random_sequence};
use crate::trace::{TraceConfig, TraceEngine};
use clap::Parser;
use std::io::Write;
use std::path::Path;

#[test]
fn owner_cli_routes_without_legacy_arguments() {
    let arguments =
        "jam trace --query query.fa --owner-index owners.json --output trace.jsonl --no-sketch";
    let parsed = Cli::try_parse_from(arguments.split_whitespace()).unwrap();
    assert!(matches!(
        parsed.command,
        Commands::Trace {
            owner_index: Some(_),
            ..
        }
    ));
    let fixture = build_fixture_at(None);
    let query = fixture.root.join("query-cli.fa");
    let output = fixture.root.join("trace-cli.jsonl");
    write_records(&query, &[("synthetic-plasmid", &fixture.query)]);
    handle_trace_command(TraceArgs {
        query,
        input: TraceInput::Owner(fixture.owner_manifest),
        audit_index: true,
        output: output.clone(),
        query_id: None,
        config: TraceConfig {
            use_sketch: false,
            circular: false,
            ..TraceConfig::default()
        },
        s3: None,
        force: false,
    })
    .unwrap();
    assert!(
        std::fs::read_to_string(output)
            .unwrap()
            .contains("synthetic-plasmid")
    );
}
#[test]
#[ignore = "writes a retained synthetic benchmark fixture"]
fn export_owner_latency_fixture() {
    let output = std::env::var_os("JAM_OWNER_FIXTURE_OUTPUT")
        .map(std::path::PathBuf::from)
        .expect("JAM_OWNER_FIXTURE_OUTPUT must be explicit");
    assert!(!output.try_exists().unwrap(), "output already exists");
    let fixture = build_fixture_at(Some(output.clone()));
    assert_eq!(fixture.root, output);
    write_records(
        &output.join("query-hit.fa"),
        &[("synthetic-plasmid", fixture.query.as_slice())],
    );
    let unrelated = random_sequence(90_001, fixture.query.len());
    write_records(
        &output.join("query-unrelated.fa"),
        &[("unrelated", unrelated.as_slice())],
    );
    let batch = (0..8)
        .map(|ordinal| {
            (
                format!("unrelated-{ordinal}"),
                random_sequence(100_003 + ordinal, fixture.query.len()),
            )
        })
        .collect::<Vec<_>>();
    let borrowed = batch
        .iter()
        .map(|(name, sequence)| (name.as_str(), sequence.as_slice()))
        .collect::<Vec<_>>();
    write_records(&output.join("query-low-reuse-batch.fa"), &borrowed);
    let owner = TraceEngine::open_owner(&fixture.owner_manifest, None).unwrap();
    owner
        .search(
            "synthetic-plasmid",
            &fixture.query,
            TraceConfig {
                use_sketch: false,
                circular: false,
                ..TraceConfig::default()
            },
        )
        .unwrap();
}

fn write_records(path: &Path, records: &[(&str, &[u8])]) {
    let mut output = std::fs::File::create(path).unwrap();
    for (name, sequence) in records {
        writeln!(output, ">{name}").unwrap();
        output.write_all(sequence).unwrap();
        output.write_all(b"\n").unwrap();
    }
}
