//! Shared loader and deterministic sequence helpers for trace fixture tests.
//!
//! Integration tests can include this module with:
//!
//! ```text
//! #[path = "trace_fixture_support/mod.rs"]
//! mod trace_fixture_support;
//! ```
//!
//! Expected relationships come from `truth.json`; this module never calls Jam
//! and therefore cannot accidentally turn a current result into test truth.

use serde::Deserialize;
use std::fs;
use std::path::{Path, PathBuf};

pub const FIXTURE_SCHEMA_VERSION: &str = "jam-rs-trace-fixtures-v1";

#[derive(Clone, Debug, Deserialize)]
pub struct FixtureManifest {
    pub schema_version: String,
    pub generator: String,
    pub catalog_file: String,
    pub datasets: Vec<DatasetTruth>,
    pub catalog: Vec<CatalogTruth>,
    pub cases: Vec<CaseTruth>,
}

#[derive(Clone, Debug, Deserialize)]
pub struct DatasetTruth {
    pub id: String,
    pub file: String,
    pub kind: String,
    pub contig_count: Option<usize>,
}

#[derive(Clone, Debug, Deserialize)]
pub struct CatalogTruth {
    pub reference: String,
    pub length: usize,
}

#[derive(Clone, Debug, Deserialize)]
pub struct CaseTruth {
    pub id: String,
    pub dataset: String,
    pub contigs: Vec<String>,
    pub references: Vec<String>,
    pub relation: String,
    pub query_intervals: Vec<IntervalTruth>,
    pub expected: ExpectedTruth,
    pub edits: Vec<EditTruth>,
    pub notes: String,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq)]
pub struct IntervalTruth {
    pub start: usize,
    pub end: usize,
}

impl IntervalTruth {
    #[must_use]
    pub fn len(self) -> usize {
        self.end - self.start
    }
}

#[derive(Clone, Debug, Deserialize)]
pub struct ExpectedTruth {
    pub status: String,
    pub reverse_complement: bool,
    pub origin_crossing: bool,
    pub minimum_query_coverage: f64,
    pub minimum_reference_coverage: f64,
    pub minimum_identity: f64,
    pub unique_reference: bool,
    pub union_deduplicates: bool,
}

#[derive(Clone, Debug, Deserialize)]
pub struct EditTruth {
    pub operation: String,
    pub query_offset: Option<usize>,
    pub query_offsets: Option<Vec<usize>>,
    pub from: Option<String>,
    pub to: Option<String>,
    pub sequence: Option<String>,
    pub length: Option<usize>,
    pub base: Option<String>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FastaRecord {
    pub id: String,
    pub sequence: String,
}

#[must_use]
pub fn fixture_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/trace")
}

#[must_use]
pub fn fixture_path(file: &str) -> PathBuf {
    fixture_dir().join(file)
}

pub fn load_manifest() -> FixtureManifest {
    let bytes = fs::read(fixture_path("truth.json")).expect("trace fixture truth.json");
    let manifest: FixtureManifest =
        serde_json::from_slice(&bytes).expect("valid trace fixture truth.json");
    assert_eq!(manifest.schema_version, FIXTURE_SCHEMA_VERSION);
    manifest
}

pub fn read_fixture(file: &str) -> String {
    fs::read_to_string(fixture_path(file)).expect("trace fixture file")
}

pub fn dataset<'a>(manifest: &'a FixtureManifest, id: &str) -> &'a DatasetTruth {
    manifest
        .datasets
        .iter()
        .find(|dataset| dataset.id == id)
        .unwrap_or_else(|| panic!("unknown trace fixture dataset {id}"))
}

pub fn case<'a>(manifest: &'a FixtureManifest, id: &str) -> &'a CaseTruth {
    manifest
        .cases
        .iter()
        .find(|case| case.id == id)
        .unwrap_or_else(|| panic!("unknown trace fixture case {id}"))
}

pub fn parse_fasta(input: &str) -> Result<Vec<FastaRecord>, String> {
    let mut records = Vec::new();
    let mut current_id: Option<String> = None;
    let mut current_sequence = String::new();

    for (line_number, raw_line) in input.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                if current_sequence.is_empty() {
                    return Err(format!("record {id} has an empty sequence"));
                }
                records.push(FastaRecord {
                    id,
                    sequence: std::mem::take(&mut current_sequence),
                });
            }
            let id = header.split_whitespace().next().unwrap_or_default();
            if id.is_empty() {
                return Err(format!("line {} has an empty FASTA identifier", line_number + 1));
            }
            current_id = Some(id.to_string());
        } else if current_id.is_none() {
            return Err(format!("line {} is sequence data before a FASTA header", line_number + 1));
        } else {
            if line.bytes().any(|base| !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n')) {
                return Err(format!("line {} contains a non-DNA symbol", line_number + 1));
            }
            current_sequence.push_str(line);
        }
    }

    if let Some(id) = current_id {
        if current_sequence.is_empty() {
            return Err(format!("record {id} has an empty sequence"));
        }
        records.push(FastaRecord {
            id,
            sequence: current_sequence,
        });
    }

    Ok(records)
}

pub fn read_fasta(file: &str) -> Result<Vec<FastaRecord>, String> {
    parse_fasta(&read_fixture(file))
}

#[must_use]
pub fn reverse_complement(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            b'N' => 'N',
            other => panic!("cannot reverse-complement non-DNA base {other:?}"),
        })
        .collect()
}

pub fn sequence_for_intervals(sequence: &str, intervals: &[IntervalTruth]) -> String {
    intervals
        .iter()
        .map(|interval| {
            sequence
                .get(interval.start..interval.end)
                .unwrap_or_else(|| panic!("fixture interval out of bounds: {interval:?}"))
        })
        .collect()
}
