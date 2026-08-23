//! Generate the deterministic artificial inputs used by trace integration tests.
//!
//! This helper intentionally has no crate dependencies so it can be run from a
//! checkout without building the whole project:
//!
//!     rustc tools/generate_trace_fixtures.rs -o /tmp/generate-trace-fixtures
//!     /tmp/generate-trace-fixtures tests/data/trace

use std::env;
use std::fs;
use std::path::{Path, PathBuf};

const TRUTH_JSON: &str = include_str!("../tests/data/trace/truth.json");

fn deterministic_dna(seed: u64, length: usize) -> String {
    let mut state = seed;
    let bases = *b"ACGT";
    let mut sequence = String::with_capacity(length);
    for _ in 0..length {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        sequence.push(bases[(state >> 62) as usize] as char);
    }
    sequence
}

fn fasta_record(name: &str, sequence: &str) -> String {
    format!(">{name}\n{sequence}\n")
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
            other => panic!("cannot reverse-complement non-DNA base {other:?}"),
        })
        .collect()
}

fn write_file(directory: &Path, name: &str, contents: &str) {
    let path = directory.join(name);
    fs::write(&path, contents).unwrap_or_else(|error| {
        panic!("failed to write {}: {error}", path.display());
    });
}

fn build_catalog() -> (String, String, String) {
    let alpha = deterministic_dna(11, 192);
    let mobile = deterministic_dna(77, 48);
    let beta = format!(
        "{}{}{}",
        deterministic_dna(22, 64),
        mobile,
        deterministic_dna(23, 64)
    );
    let gamma = format!(
        "{}{}{}",
        deterministic_dna(33, 64),
        mobile,
        deterministic_dna(34, 64)
    );
    (alpha, beta, gamma)
}

fn build_trace_assembly(alpha: &str, beta: &str) -> String {
    let reverse_fragment = reverse_complement(&alpha[32..128]);
    let origin_crossing = format!("{}{}", &alpha[160..192], &alpha[..48]);

    let mut substitutions = alpha[48..144].as_bytes().to_vec();
    for (offset, replacement) in [(10, b'C'), (50, b'G'), (90, b'A')] {
        substitutions[offset] = replacement;
    }
    let substitutions = String::from_utf8(substitutions).unwrap();

    let insertion_source = &alpha[64..128];
    let insertion = format!(
        "{}GGGCCC{}",
        &insertion_source[..30],
        &insertion_source[30..]
    );

    let deletion_source = &alpha[80..144];
    let deletion = format!("{}{}", &deletion_source[..20], &deletion_source[25..]);

    let mut ambiguous = alpha[..80].as_bytes().to_vec();
    for offset in [20, 21, 22, 50] {
        ambiguous[offset] = b'N';
    }
    let ambiguous = String::from_utf8(ambiguous).unwrap();

    [
        fasta_record("alpha_exact", alpha),
        fasta_record("alpha_rc_fragment", &reverse_fragment),
        fasta_record("alpha_origin_crossing", &origin_crossing),
        fasta_record("alpha_substitutions", &substitutions),
        fasta_record("alpha_insertion", &insertion),
        fasta_record("alpha_deletion", &deletion),
        fasta_record("beta_split_1", &beta[..80]),
        fasta_record("beta_split_2", &beta[80..]),
        fasta_record("beta_overlap_1", &beta[..100]),
        fasta_record("beta_overlap_2", &beta[80..160]),
        fasta_record("shared_mobile", &beta[64..112]),
        fasta_record("unrelated_chromosome", &deterministic_dna(900, 192)),
        fasta_record("alpha_ambiguous", &ambiguous),
    ]
    .concat()
}

fn generate(directory: &Path) {
    fs::create_dir_all(directory).unwrap_or_else(|error| {
        panic!("failed to create {}: {error}", directory.display());
    });

    let (alpha, beta, gamma) = build_catalog();
    write_file(
        directory,
        "catalog.fa",
        &[
            fasta_record("plasmid_alpha", &alpha),
            fasta_record("plasmid_beta", &beta),
            fasta_record("plasmid_gamma", &gamma),
        ]
        .concat(),
    );
    write_file(directory, "assembly_trace.fa", &build_trace_assembly(&alpha, &beta));
    write_file(directory, "assembly_empty.fa", "");
    write_file(
        directory,
        "assembly_negative.fa",
        &fasta_record("unrelated_only", &deterministic_dna(900, 192)),
    );
    write_file(directory, "assembly_malformed.fa", "not-a-fasta-header\nACGTACGTACGT\n");
    write_file(directory, "truth.json", TRUTH_JSON);
}

fn main() {
    let mut arguments = env::args_os();
    let program = arguments.next().unwrap_or_else(|| "generate_trace_fixtures".into());
    let directory = arguments.next().map(PathBuf::from).unwrap_or_else(|| {
        eprintln!("usage: {} OUTPUT_DIRECTORY", Path::new(&program).display());
        std::process::exit(2);
    });
    if arguments.next().is_some() {
        eprintln!("usage: {} OUTPUT_DIRECTORY", Path::new(&program).display());
        std::process::exit(2);
    }
    generate(&directory);
}
