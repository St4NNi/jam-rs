//! Loader for the deterministic mobile-genetic-element trace fixtures.
//!
//! The expected relationships in `truth.json` are construction metadata.  This
//! helper only reads those records and never consults a Jam result.

#![allow(dead_code)]

use serde::Deserialize;
use std::fs;
use std::path::{Path, PathBuf};

pub const FIXTURE_SCHEMA_VERSION: &str = "jam-rs-trace-mge-v1";

#[derive(Clone, Debug, Deserialize)]
pub struct FixtureManifest {
    pub schema_version: String,
    pub generator: String,
    pub truth_source: String,
    pub queries: Vec<QueryTruth>,
    pub datasets: Vec<DatasetTruth>,
    pub cases: Vec<CaseTruth>,
}

#[derive(Clone, Debug, Deserialize)]
pub struct QueryTruth {
    pub id: String,
    pub file: String,
    pub kind: String,
    pub topology_requested: String,
    pub length: usize,
    pub sha256: String,
    pub notes: String,
}

#[derive(Clone, Debug, Deserialize)]
pub struct DatasetTruth {
    pub id: String,
    pub file: String,
    pub kind: String,
    pub contig_count: usize,
    pub sha256: String,
    pub notes: String,
}

#[derive(Clone, Debug, Deserialize)]
pub struct CaseTruth {
    pub id: String,
    pub query: String,
    pub query_kind: String,
    pub topology_requested: String,
    pub dataset: String,
    pub relation: String,
    pub contigs: Vec<String>,
    pub query_intervals: Vec<IntervalTruth>,
    pub expected: ExpectedTruth,
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
    pub coordinate_model: String,
    pub topology_evidence: String,
    pub wrap_compatible: bool,
    pub unique_query_coverage_bases: usize,
    pub common_sequence: bool,
    pub repeat_label: String,
    pub independent_contig_identity: bool,
    pub contig_identity: Vec<ContigIdentityTruth>,
}

#[derive(Clone, Debug, Deserialize)]
pub struct ContigIdentityTruth {
    pub contig: String,
    pub independent: bool,
    pub length: usize,
    pub sequence_sha256: String,
    pub common_sequence: bool,
    pub repeat_label: String,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FastaRecord {
    pub id: String,
    pub sequence: String,
}

#[must_use]
pub fn fixture_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/trace_mge")
}

#[must_use]
pub fn fixture_path(file: &str) -> PathBuf {
    fixture_dir().join(file)
}

pub fn load_manifest() -> FixtureManifest {
    let bytes = fs::read(fixture_path("truth.json")).expect("MGE fixture truth.json");
    let manifest: FixtureManifest =
        serde_json::from_slice(&bytes).expect("valid MGE fixture truth.json");
    assert_eq!(manifest.schema_version, FIXTURE_SCHEMA_VERSION);
    assert_eq!(manifest.truth_source, "construction_coordinates");
    manifest
}

pub fn read_fixture(file: &str) -> String {
    fs::read_to_string(fixture_path(file)).expect("MGE fixture file")
}

pub fn read_bytes(file: &str) -> Vec<u8> {
    fs::read(fixture_path(file)).expect("MGE fixture bytes")
}

pub fn query<'a>(manifest: &'a FixtureManifest, id: &str) -> &'a QueryTruth {
    manifest
        .queries
        .iter()
        .find(|query| query.id == id)
        .unwrap_or_else(|| panic!("unknown MGE query {id}"))
}

pub fn dataset<'a>(manifest: &'a FixtureManifest, id: &str) -> &'a DatasetTruth {
    manifest
        .datasets
        .iter()
        .find(|dataset| dataset.id == id)
        .unwrap_or_else(|| panic!("unknown MGE dataset {id}"))
}

pub fn case<'a>(manifest: &'a FixtureManifest, id: &str) -> &'a CaseTruth {
    manifest
        .cases
        .iter()
        .find(|case| case.id == id)
        .unwrap_or_else(|| panic!("unknown MGE fixture case {id}"))
}

pub fn parse_fasta(input: &str) -> Result<Vec<FastaRecord>, String> {
    let mut records = Vec::new();
    let mut current_id = None;
    let mut sequence = String::new();
    for (line_number, raw_line) in input.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                if sequence.is_empty() {
                    return Err(format!("record {id} has an empty sequence"));
                }
                records.push(FastaRecord {
                    id,
                    sequence: std::mem::take(&mut sequence),
                });
            }
            let id = header.split_whitespace().next().unwrap_or_default();
            if id.is_empty() {
                return Err(format!(
                    "line {} has an empty FASTA identifier",
                    line_number + 1
                ));
            }
            current_id = Some(id.to_owned());
        } else {
            if current_id.is_none() {
                return Err(format!(
                    "line {} has sequence before a FASTA header",
                    line_number + 1
                ));
            }
            if line.bytes().any(|base| {
                !matches!(
                    base,
                    b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n'
                )
            }) {
                return Err(format!(
                    "line {} contains a non-DNA symbol",
                    line_number + 1
                ));
            }
            sequence.push_str(line);
        }
    }
    if let Some(id) = current_id {
        if sequence.is_empty() {
            return Err(format!("record {id} has an empty sequence"));
        }
        records.push(FastaRecord { id, sequence });
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
            other => panic!("cannot reverse-complement base {other:?}"),
        })
        .collect()
}

pub fn sequence_for_intervals(sequence: &str, intervals: &[IntervalTruth]) -> String {
    intervals
        .iter()
        .map(|interval| {
            sequence
                .get(interval.start..interval.end)
                .unwrap_or_else(|| panic!("MGE fixture interval out of bounds: {interval:?}"))
        })
        .collect()
}
