use crate::cli::{
    Cli, Commands,
    handlers::{TraceArgs, TraceInput, handle_trace_command},
};
use crate::jidx_reader::JidxReader;
use crate::owner_tests::{build_fixture_at, build_scale_fixture_at, random_sequence};
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

#[test]
#[ignore = "writes a retained roughly 100k-key query-ready fixture"]
fn export_owner_scale_fixture() {
    let output = std::env::var_os("JAM_OWNER_SCALE_FIXTURE_OUTPUT")
        .map(std::path::PathBuf::from)
        .expect("JAM_OWNER_SCALE_FIXTURE_OUTPUT must be explicit");
    assert!(!output.try_exists().unwrap(), "output already exists");
    let scale = build_scale_fixture_at(output.clone(), 250_000);
    let source = JidxReader::open(&scale.fixture.jidx).unwrap();
    let keys = source.header().seed_count;
    assert!(
        (50_000..=200_000).contains(&keys),
        "unexpected key count {keys}"
    );
    write_records(
        &output.join("query-primary-hit.fa"),
        &[("synthetic-plasmid", scale.fixture.query.as_slice())],
    );
    let independent = scale
        .independent_queries
        .iter()
        .map(|(name, sequence)| (name.as_str(), sequence.as_slice()))
        .collect::<Vec<_>>();
    write_records(&output.join("query-independent-hits.fa"), &independent);
    let unrelated = random_sequence(900_001, scale.fixture.query.len());
    write_records(
        &output.join("query-unrelated.fa"),
        &[("unrelated", unrelated.as_slice())],
    );
    let low_reuse = (0..8)
        .map(|ordinal| {
            (
                format!("low-reuse-{ordinal}"),
                random_sequence(910_001 + ordinal, scale.fixture.query.len()),
            )
        })
        .collect::<Vec<_>>();
    let low_reuse_records = low_reuse
        .iter()
        .map(|(name, sequence)| (name.as_str(), sequence.as_slice()))
        .collect::<Vec<_>>();
    write_records(&output.join("query-low-reuse.fa"), &low_reuse_records);
    let owner = TraceEngine::open_owner(&scale.fixture.owner_manifest, None).unwrap();
    let config = TraceConfig {
        use_sketch: false,
        circular: false,
        ..TraceConfig::default()
    };
    let queries = scale.independent_queries.clone();
    let results = owner.search_batch(&queries, config).unwrap();
    assert!(results.iter().all(|result| {
        result
            .metagenomes
            .iter()
            .any(|metagenome| metagenome.name == "doc-c")
    }));
    let unavailable = output.join("legacy-unavailable");
    std::fs::create_dir(&unavailable).unwrap();
    let legacy = [
        "synthetic.jam",
        "synthetic.jidx",
        "source.json",
        "doc-a.fa",
        "doc-b.fa",
        "doc-c.fa",
        "doc-a.bgz.fai",
        "doc-b.bgz.fai",
        "doc-c.bgz.fai",
        "doc-a.bgz.gzi",
        "doc-b.bgz.gzi",
        "doc-c.bgz.gzi",
    ];
    for name in legacy {
        std::fs::rename(output.join(name), unavailable.join(name)).unwrap();
    }
    let isolated = TraceEngine::open_owner(&scale.fixture.owner_manifest, None).unwrap();
    assert_eq!(isolated.search_batch(&queries, config).unwrap(), results);
    for name in legacy {
        std::fs::rename(unavailable.join(name), output.join(name)).unwrap();
    }
    std::fs::remove_dir(unavailable).unwrap();
    std::fs::write(
        output.join("scale-summary.json"),
        serde_json::to_vec_pretty(&serde_json::json!({
            "status": "complete",
            "keys": keys,
            "documents": source.header().document_count,
            "contigs": source.header().contig_count,
            "background_bases": 250000,
            "independent_positive_queries": results.len(),
            "isolated_owner_search": true
        }))
        .unwrap(),
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
