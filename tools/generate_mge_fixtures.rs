//! Generate the deterministic mobile-genetic-element trace fixtures.
//!
//! The generator has no crate dependencies.  It records the construction
//! coordinates directly in `truth.json`; no Jam output is read while making
//! the expected relationships.  Regenerate with:
//!
//! ```text
//! rustc tools/generate_mge_fixtures.rs -o tools/out/generate-mge-fixtures
//! tools/out/generate-mge-fixtures tests/data/trace_mge
//! ```

use std::env;
use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};

const SCHEMA_VERSION: &str = "jam-rs-trace-mge-v1";

#[derive(Clone)]
struct Query {
    id: &'static str,
    file: &'static str,
    kind: &'static str,
    topology: &'static str,
    sequence: String,
    notes: &'static str,
}

#[derive(Clone)]
struct Contig {
    id: &'static str,
    sequence: String,
}

#[derive(Clone)]
struct Dataset {
    id: &'static str,
    file: &'static str,
    kind: &'static str,
    contigs: Vec<Contig>,
    notes: &'static str,
}

#[derive(Clone, Copy)]
struct Interval {
    start: usize,
    end: usize,
}

#[derive(Clone)]
struct Case {
    id: &'static str,
    query: &'static str,
    dataset: &'static str,
    relation: &'static str,
    contigs: &'static [&'static str],
    intervals: &'static [Interval],
    status: &'static str,
    coordinate_model: &'static str,
    topology_evidence: &'static str,
    wrap_compatible: bool,
    unique_coverage: usize,
    common_sequence: bool,
    repeat_label: &'static str,
    notes: &'static str,
}

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

fn reverse_complement(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            other => panic!("cannot reverse-complement base {:?}", other),
        })
        .collect()
}

fn rotate(sequence: &str, offset: usize) -> String {
    let offset = offset % sequence.len();
    format!("{}{}", &sequence[offset..], &sequence[..offset])
}

fn fasta_record(id: &str, sequence: &str) -> String {
    let mut output = format!(">{id}\n");
    for chunk in sequence.as_bytes().chunks(80) {
        output.push_str(std::str::from_utf8(chunk).expect("fixture DNA is UTF-8"));
        output.push('\n');
    }
    output
}

fn fasta_records(records: &[(&str, &str)]) -> String {
    records
        .iter()
        .map(|(id, sequence)| fasta_record(id, sequence))
        .collect()
}

fn query_record_id(query: &Query) -> &str {
    if query.id == "circular_plasmid_auto" {
        "circular_plasmid_rotated"
    } else {
        query.id
    }
}

fn write_file(directory: &Path, name: &str, contents: &str) {
    fs::write(directory.join(name), contents)
        .unwrap_or_else(|error| panic!("failed to write {}: {}", name, error));
}

// Small, dependency-free SHA-256 implementation used only to bind the checked
// fixture files to their independently specified truth manifest.
fn sha256(bytes: &[u8]) -> [u8; 32] {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1,
        0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
        0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786,
        0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147,
        0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
        0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a,
        0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
        0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
    ];
    let mut message = bytes.to_vec();
    let bit_length = (message.len() as u64) * 8;
    message.push(0x80);
    while message.len() % 64 != 56 {
        message.push(0);
    }
    message.extend_from_slice(&bit_length.to_be_bytes());

    let mut state = [
        0x6a09e667_u32,
        0xbb67ae85,
        0x3c6ef372,
        0xa54ff53a,
        0x510e527f,
        0x9b05688c,
        0x1f83d9ab,
        0x5be0cd19,
    ];
    for block in message.chunks_exact(64) {
        let mut schedule = [0_u32; 64];
        for (index, word) in block.chunks_exact(4).take(16).enumerate() {
            schedule[index] = u32::from_be_bytes([word[0], word[1], word[2], word[3]]);
        }
        for index in 16..64 {
            let s0 = schedule[index - 15].rotate_right(7)
                ^ schedule[index - 15].rotate_right(18)
                ^ (schedule[index - 15] >> 3);
            let s1 = schedule[index - 2].rotate_right(17)
                ^ schedule[index - 2].rotate_right(19)
                ^ (schedule[index - 2] >> 10);
            schedule[index] = schedule[index - 16]
                .wrapping_add(s0)
                .wrapping_add(schedule[index - 7])
                .wrapping_add(s1);
        }
        let mut working = state;
        for index in 0..64 {
            let sigma1 = working[4].rotate_right(6)
                ^ working[4].rotate_right(11)
                ^ working[4].rotate_right(25);
            let choice = (working[4] & working[5]) ^ ((!working[4]) & working[6]);
            let temp1 = working[7]
                .wrapping_add(sigma1)
                .wrapping_add(choice)
                .wrapping_add(K[index])
                .wrapping_add(schedule[index]);
            let sigma0 = working[0].rotate_right(2)
                ^ working[0].rotate_right(13)
                ^ working[0].rotate_right(22);
            let majority = (working[0] & working[1])
                ^ (working[0] & working[2])
                ^ (working[1] & working[2]);
            let temp2 = sigma0.wrapping_add(majority);
            working = [
                temp1.wrapping_add(temp2),
                working[0],
                working[1],
                working[2],
                working[3].wrapping_add(temp1),
                working[4],
                working[5],
                working[6],
            ];
        }
        for index in 0..8 {
            state[index] = state[index].wrapping_add(working[index]);
        }
    }
    let mut digest = [0_u8; 32];
    for (index, word) in state.iter().enumerate() {
        digest[index * 4..index * 4 + 4].copy_from_slice(&word.to_be_bytes());
    }
    digest
}

fn hex_digest(bytes: &[u8]) -> String {
    let digest = sha256(bytes);
    digest.iter().fold(String::with_capacity(64), |mut result, byte| {
        write!(result, "{byte:02x}").expect("writing to String cannot fail");
        result
    })
}

fn json_string(value: &str) -> String {
    let mut result = String::with_capacity(value.len() + 2);
    result.push('"');
    for character in value.chars() {
        match character {
            '"' => result.push_str("\\\""),
            '\\' => result.push_str("\\\\"),
            '\n' => result.push_str("\\n"),
            '\r' => result.push_str("\\r"),
            '\t' => result.push_str("\\t"),
            other if other.is_control() => {
                write!(result, "\\u{:04x}", other as u32).expect("writing to String cannot fail");
            }
            other => result.push(other),
        }
    }
    result.push('"');
    result
}

fn json_strings(values: &[&str]) -> String {
    values
        .iter()
        .map(|value| json_string(value))
        .collect::<Vec<_>>()
        .join(", ")
}

fn json_intervals(values: &[Interval]) -> String {
    values
        .iter()
        .map(|interval| format!("{{\"start\": {}, \"end\": {}}}", interval.start, interval.end))
        .collect::<Vec<_>>()
        .join(", ")
}

fn json_contig_identity(case: &Case, datasets: &[Dataset]) -> String {
    let dataset = datasets
        .iter()
        .find(|dataset| dataset.id == case.dataset)
        .expect("case dataset exists");
    case.contigs
        .iter()
        .map(|contig| {
            let source = dataset
                .contigs
                .iter()
                .find(|candidate| candidate.id == *contig)
                .expect("case contig exists");
            format!(
                "{{\"contig\": {}, \"independent\": true, \"length\": {}, \"sequence_sha256\": {}, \"common_sequence\": {}, \"repeat_label\": {}}}",
                json_string(contig),
                source.sequence.len(),
                json_string(&hex_digest(source.sequence.as_bytes())),
                case.common_sequence,
                json_string(case.repeat_label)
            )
        })
        .collect::<Vec<_>>()
        .join(", ")
}

fn build_inputs() -> (Vec<Query>, Vec<Dataset>) {
    let transposon = deterministic_dna(0x6006, 40);
    let resistance = deterministic_dna(0x7007, 48);

    let circular_plasmid = format!(
        "{}{}{}",
        deterministic_dna(0x1010, 92),
        transposon,
        deterministic_dna(0x1011, 108)
    );
    let linear_plasmid = format!(
        "{}{}{}",
        deterministic_dna(0x2020, 96),
        resistance,
        deterministic_dna(0x2021, 80)
    );
    let phage_linear = deterministic_dna(0x3030, 256);
    let terminal_repeat = deterministic_dna(0x3031, 32);
    let phage_direct_repeat = format!(
        "{}{}{}",
        terminal_repeat,
        deterministic_dna(0x3032, 224),
        terminal_repeat
    );
    let phage_inverted_repeat = format!(
        "{}{}{}",
        terminal_repeat,
        deterministic_dna(0x3033, 224),
        reverse_complement(&terminal_repeat)
    );
    let synthetic = format!(
        "{}{}{}{}",
        deterministic_dna(0x4040, 88),
        transposon,
        deterministic_dna(0x4041, 100),
        resistance
    );
    let unknown = deterministic_dna(0x5050, 180);

    let queries = vec![
        Query {
            id: "circular_plasmid",
            file: "query_circular_plasmid.fa",
            kind: "plasmid",
            topology: "circular",
            sequence: circular_plasmid.clone(),
            notes: "Synthetic circular plasmid with a shared transposon segment.",
        },
        Query {
            id: "circular_plasmid_rotated",
            file: "query_circular_plasmid_rotated.fa",
            kind: "plasmid",
            topology: "circular",
            sequence: rotate(&circular_plasmid, 73),
            notes: "The same circular molecule represented with a different stored origin.",
        },
        Query {
            id: "circular_plasmid_auto",
            file: "query_circular_plasmid_auto.fa",
            kind: "plasmid",
            topology: "auto",
            sequence: rotate(&circular_plasmid, 73),
            notes: "The same rotated representation evaluated without a biological topology claim.",
        },
        Query {
            id: "linear_plasmid",
            file: "query_linear_plasmid.fa",
            kind: "plasmid",
            topology: "linear",
            sequence: linear_plasmid.clone(),
            notes: "Synthetic linear plasmid; terminal query gaps must remain distinct.",
        },
        Query {
            id: "linear_phage",
            file: "query_linear_phage.fa",
            kind: "phage",
            topology: "linear",
            sequence: phage_linear.clone(),
            notes: "Synthetic linear phage used for partial and integrated-fragment cases.",
        },
        Query {
            id: "phage_terminal_direct_repeat",
            file: "query_phage_direct_repeat.fa",
            kind: "phage",
            topology: "linear",
            sequence: phage_direct_repeat.clone(),
            notes: "Linear phage with identical terminal direct repeats.",
        },
        Query {
            id: "phage_terminal_inverted_repeat",
            file: "query_phage_inverted_repeat.fa",
            kind: "phage",
            topology: "linear",
            sequence: phage_inverted_repeat.clone(),
            notes: "Linear phage with reverse-complement terminal repeats.",
        },
        Query {
            id: "synthetic_construct",
            file: "query_synthetic_construct.fa",
            kind: "other",
            topology: "linear",
            sequence: synthetic.clone(),
            notes: "Synthetic construct containing shared transposon and resistance-cassette segments.",
        },
        Query {
            id: "unknown_element",
            file: "query_unknown_element.fa",
            kind: "unknown",
            topology: "unknown",
            sequence: unknown.clone(),
            notes: "Unknown sequence class and topology; no biological label is implied.",
        },
    ];

    let circular_contigs = vec![
        Contig { id: "circ_full", sequence: circular_plasmid.clone() },
        Contig { id: "circ_reverse_fragment", sequence: reverse_complement(&circular_plasmid[32..128]) },
        Contig { id: "circ_split_left", sequence: circular_plasmid[..140].to_string() },
        Contig { id: "circ_split_right", sequence: circular_plasmid[120..].to_string() },
        Contig { id: "circ_duplicate", sequence: circular_plasmid[..140].to_string() },
        Contig { id: "circ_origin_crossing", sequence: format!("{}{}", &circular_plasmid[200..], &circular_plasmid[..64]) },
        Contig { id: "circ_shared_transposon", sequence: transposon.clone() },
        Contig { id: "circ_unrelated_chromosome", sequence: deterministic_dna(0x9090, 240) },
    ];

    let linear_contigs = vec![
        Contig { id: "linear_full", sequence: linear_plasmid.clone() },
        Contig { id: "linear_overlap_left", sequence: linear_plasmid[..140].to_string() },
        Contig { id: "linear_overlap_right", sequence: linear_plasmid[120..].to_string() },
        Contig { id: "linear_duplicate", sequence: linear_plasmid[..140].to_string() },
        Contig { id: "linear_middle_fragment", sequence: linear_plasmid[40..180].to_string() },
        Contig { id: "linear_shared_resistance", sequence: resistance.clone() },
        Contig { id: "linear_unrelated_chromosome", sequence: deterministic_dna(0x9091, 224) },
    ];

    let phage_contigs = vec![
        Contig { id: "phage_full", sequence: phage_linear.clone() },
        Contig { id: "phage_reverse_fragment", sequence: reverse_complement(&phage_linear[80..208]) },
        Contig { id: "phage_integrated_fragment", sequence: format!("{}{}{}", deterministic_dna(0x8080, 40), &phage_linear[48..176], deterministic_dna(0x8081, 40)) },
        Contig { id: "phage_direct_full", sequence: phage_direct_repeat.clone() },
        Contig { id: "phage_direct_one_terminal", sequence: format!("{}{}", &phage_direct_repeat[..256], "") },
        Contig { id: "phage_inverted_full", sequence: phage_inverted_repeat.clone() },
        Contig { id: "phage_inverted_core", sequence: phage_inverted_repeat[32..256].to_string() },
        Contig { id: "phage_unrelated_chromosome", sequence: deterministic_dna(0x9092, 256) },
    ];

    let common_contigs = vec![
        Contig { id: "common_transposon_a", sequence: transposon.clone() },
        Contig { id: "common_transposon_b", sequence: transposon },
        Contig { id: "common_resistance_a", sequence: resistance.clone() },
        Contig { id: "common_resistance_b", sequence: resistance },
        Contig { id: "synthetic_full", sequence: synthetic.clone() },
        Contig { id: "synthetic_integrated", sequence: format!("{}{}{}", deterministic_dna(0x8180, 36), &synthetic[50..230], deterministic_dna(0x8181, 36)) },
        Contig { id: "unknown_full", sequence: unknown.clone() },
        Contig { id: "common_unrelated_chromosome", sequence: deterministic_dna(0x9093, 300) },
    ];

    let datasets = vec![
        Dataset { id: "circular", file: "assembly_circular.fa", kind: "circular_plasmid_mixed", contigs: circular_contigs, notes: "Exact, reverse, split, duplicate, overlap, wrap, shared-element, and negative contigs." },
        Dataset { id: "linear", file: "assembly_linear.fa", kind: "linear_plasmid_mixed", contigs: linear_contigs, notes: "Linear exact, overlapping, duplicate, partial, shared-cassette, and negative contigs." },
        Dataset { id: "phages", file: "assembly_phages.fa", kind: "phage_mixed", contigs: phage_contigs, notes: "Linear phage, terminal repeats, reverse fragment, and integrated element cases." },
        Dataset { id: "common", file: "assembly_common_elements.fa", kind: "common_mobile_elements", contigs: common_contigs, notes: "Shared transposon/cassette and integrated synthetic construct evidence." },
        Dataset { id: "empty", file: "assembly_empty.fa", kind: "empty", contigs: Vec::new(), notes: "Empty FASTA input must still produce output headers." },
        Dataset { id: "malformed", file: "assembly_malformed.fa", kind: "malformed", contigs: Vec::new(), notes: "Sequence data appears before a FASTA header." },
    ];

    (queries, datasets)
}

fn make_case(
    id: &'static str,
    query: &'static str,
    dataset: &'static str,
    relation: &'static str,
    contigs: &'static [&'static str],
    intervals: &'static [Interval],
    status: &'static str,
    coordinate_model: &'static str,
    topology_evidence: &'static str,
    wrap_compatible: bool,
    unique_coverage: usize,
    common_sequence: bool,
    repeat_label: &'static str,
    notes: &'static str,
) -> Case {
    Case { id, query, dataset, relation, contigs, intervals, status, coordinate_model, topology_evidence, wrap_compatible, unique_coverage, common_sequence, repeat_label, notes }
}

fn build_cases() -> Vec<Case> {
    vec![
        make_case("circular_exact", "circular_plasmid", "circular", "exact", &["circ_full"], &[Interval { start: 0, end: 240 }], "supported", "linear", "linear_supported", true, 240, false, "none", "The complete circular query is present as an exact forward contig."),
        make_case("circular_reverse_fragment", "circular_plasmid", "circular", "reverse_complement_fragment", &["circ_reverse_fragment"], &[Interval { start: 32, end: 128 }], "supported", "linear", "linear_supported", true, 96, false, "none", "Reverse-complement evidence must map to the same query interval."),
        make_case("circular_split_overlap", "circular_plasmid", "circular", "overlapping_contigs", &["circ_split_left", "circ_split_right"], &[Interval { start: 0, end: 240 }], "supported", "linear", "linear_supported", true, 240, false, "none", "The 20-base overlap is counted once in the query-coordinate union."),
        make_case("circular_duplicate", "circular_plasmid", "circular", "duplicate_contigs", &["circ_split_left", "circ_duplicate"], &[Interval { start: 0, end: 140 }], "supported", "linear", "linear_supported", true, 140, false, "none", "Duplicate contigs are independent evidence but cannot inflate coverage."),
        make_case("circular_origin_crossing", "circular_plasmid", "circular", "origin_crossing", &["circ_origin_crossing"], &[Interval { start: 200, end: 240 }, Interval { start: 0, end: 64 }], "supported", "wrap", "wrap_supported", true, 104, false, "none", "One contig crosses the stored query origin; its two query intervals are disjoint after normalization."),
        make_case("circular_rotated_reference", "circular_plasmid_rotated", "circular", "arbitrary_origin_rotation", &["circ_full"], &[Interval { start: 167, end: 240 }, Interval { start: 0, end: 167 }], "supported", "wrap", "wrap_supported", true, 240, false, "none", "The query is the same circular sequence rotated by 73 bases; total evidence is origin invariant."),
        make_case("circular_auto_origin_rotation", "circular_plasmid_auto", "circular", "auto_origin_rotation", &["circ_full"], &[Interval { start: 167, end: 240 }, Interval { start: 0, end: 167 }], "supported", "undetermined", "both_compatible", true, 240, false, "none", "Auto mode reuses the main alignment set and retains both coordinate models without inferring biological topology."),
        make_case("circular_shared_transposon", "circular_plasmid", "circular", "shared_transposon", &["circ_shared_transposon"], &[Interval { start: 92, end: 132 }], "ambiguous_candidate", "linear", "insufficient", true, 40, true, "common_sequence", "A shared transposon is candidate evidence only and must be labeled common sequence."),
        make_case("circular_unrelated", "circular_plasmid", "circular", "unrelated_chromosome", &["circ_unrelated_chromosome"], &[], "no_candidate", "linear", "insufficient", true, 0, false, "none", "The unrelated chromosome-like contig is a negative control."),
        make_case("linear_exact", "linear_plasmid", "linear", "exact", &["linear_full"], &[Interval { start: 0, end: 224 }], "supported", "linear", "linear_supported", false, 224, false, "none", "A linear query must not wrap or merge terminal gaps."),
        make_case("linear_split_overlap", "linear_plasmid", "linear", "overlapping_contigs", &["linear_overlap_left", "linear_overlap_right"], &[Interval { start: 0, end: 224 }], "supported", "linear", "linear_supported", false, 224, false, "none", "Overlapping linear contigs provide a query-coordinate union."),
        make_case("linear_partial", "linear_plasmid", "linear", "partial_fragment", &["linear_middle_fragment"], &[Interval { start: 40, end: 180 }], "supported", "linear", "linear_supported", false, 140, false, "none", "Leading and trailing unsupported query intervals remain separate."),
        make_case("linear_shared_resistance", "linear_plasmid", "linear", "shared_resistance_cassette", &["linear_shared_resistance"], &[Interval { start: 96, end: 144 }], "ambiguous_candidate", "linear", "insufficient", false, 48, true, "common_sequence", "A resistance cassette can map to multiple elements and is not a unique call."),
        make_case("linear_unrelated", "linear_plasmid", "linear", "unrelated_chromosome", &["linear_unrelated_chromosome"], &[], "no_candidate", "linear", "insufficient", false, 0, false, "none", "Chromosome-like sequence without designed query support."),
        make_case("phage_exact", "linear_phage", "phages", "exact", &["phage_full"], &[Interval { start: 0, end: 256 }], "supported", "linear", "linear_supported", false, 256, false, "none", "Exact linear phage evidence."),
        make_case("phage_reverse_fragment", "linear_phage", "phages", "reverse_complement_fragment", &["phage_reverse_fragment"], &[Interval { start: 80, end: 208 }], "supported", "linear", "linear_supported", false, 128, false, "none", "Reverse strand fragment with independent contig identity."),
        make_case("phage_integrated_fragment", "linear_phage", "phages", "integrated_element_fragment", &["phage_integrated_fragment"], &[Interval { start: 48, end: 176 }], "supported", "linear", "linear_supported", false, 128, false, "integrated", "Flanking chromosome sequence indicates integration context, not an autonomous phage."),
        make_case("phage_direct_repeat_exact", "phage_terminal_direct_repeat", "phages", "terminal_direct_repeat", &["phage_direct_full"], &[Interval { start: 0, end: 288 }], "supported", "linear", "linear_supported", false, 288, false, "terminal_direct_repeat", "The exact contig covers both direct-repeat copies, but unique coverage still uses a union."),
        make_case("phage_direct_repeat_one_copy", "phage_terminal_direct_repeat", "phages", "terminal_direct_repeat_one_copy", &["phage_direct_one_terminal"], &[Interval { start: 0, end: 256 }], "supported", "linear", "linear_supported", false, 256, false, "terminal_direct_repeat", "One terminal repeat plus the core must not be duplicated onto both query ends."),
        make_case("phage_inverted_repeat_exact", "phage_terminal_inverted_repeat", "phages", "terminal_inverted_repeat", &["phage_inverted_full"], &[Interval { start: 0, end: 288 }], "supported", "linear", "linear_supported", false, 288, false, "terminal_inverted_repeat", "Terminal inverted repeats retain strand-aware alternatives without inflating coverage."),
        make_case("phage_inverted_repeat_core", "phage_terminal_inverted_repeat", "phages", "terminal_repeat_core", &["phage_inverted_core"], &[Interval { start: 32, end: 256 }], "supported", "linear", "linear_supported", false, 224, false, "terminal_inverted_repeat", "The core is supported while terminal repeat evidence is absent."),
        make_case("synthetic_exact", "synthetic_construct", "common", "synthetic_construct", &["synthetic_full"], &[Interval { start: 0, end: 276 }], "supported", "linear", "linear_supported", false, 276, false, "none", "Synthetic construct evidence is supported without assigning a biological class."),
        make_case("synthetic_transposon_only", "synthetic_construct", "common", "shared_transposon", &["common_transposon_a", "common_transposon_b"], &[Interval { start: 88, end: 128 }], "ambiguous_candidate", "linear", "insufficient", false, 40, true, "common_sequence", "Two independent contigs contain only the shared transposon and remain ambiguous."),
        make_case("synthetic_resistance_only", "synthetic_construct", "common", "shared_resistance_cassette", &["common_resistance_a", "common_resistance_b"], &[Interval { start: 228, end: 276 }], "ambiguous_candidate", "linear", "insufficient", false, 48, true, "common_sequence", "A shared resistance cassette is not unique support for the construct."),
        make_case("synthetic_integrated", "synthetic_construct", "common", "integrated_element_fragment", &["synthetic_integrated"], &[Interval { start: 50, end: 230 }], "supported", "linear", "linear_supported", false, 180, false, "integrated", "An integrated construct fragment has chromosome-like flanks and is not a linkage proof."),
        make_case("unknown_exact", "unknown_element", "common", "unknown_query_exact", &["unknown_full"], &[Interval { start: 0, end: 180 }], "supported", "undetermined", "undetermined", false, 180, false, "none", "Exact sequence evidence is retained while unknown query class and topology remain undetermined."),
        make_case("empty_input", "unknown_element", "empty", "empty_assembly", &[], &[], "empty", "undetermined", "undetermined", false, 0, false, "none", "Empty input must return a valid run/result record with no topology claim."),
        make_case("malformed_input", "unknown_element", "malformed", "malformed_fasta", &[], &[], "error", "undetermined", "undetermined", false, 0, false, "none", "Malformed FASTA is an explicit structured error case."),
    ]
}

fn render_truth(
    output: &Path,
    queries: &[Query],
    datasets: &[Dataset],
    cases: &[Case],
) -> String {
    let mut truth = String::new();
    truth.push_str("{\n");
    writeln!(truth, "  \"schema_version\": {},", json_string(SCHEMA_VERSION)).unwrap();
    truth.push_str("  \"generator\": \"tools/generate_mge_fixtures.rs\",\n");
    truth.push_str("  \"truth_source\": \"construction_coordinates\",\n");
    truth.push_str("  \"queries\": [\n");
    for (index, query) in queries.iter().enumerate() {
        let bytes = fs::read(output.join(query.file)).expect("generated query");
        let comma = if index + 1 == queries.len() { "" } else { "," };
        writeln!(
            truth,
            "    {{\"id\": {}, \"file\": {}, \"kind\": {}, \"topology_requested\": {}, \"length\": {}, \"sha256\": {}, \"notes\": {}}}{}",
            json_string(query.id),
            json_string(query.file),
            json_string(query.kind),
            json_string(query.topology),
            query.sequence.len(),
            json_string(&hex_digest(&bytes)),
            json_string(query.notes),
            comma
        )
        .unwrap();
    }
    truth.push_str("  ],\n  \"datasets\": [\n");
    for (index, dataset) in datasets.iter().enumerate() {
        let bytes = fs::read(output.join(dataset.file)).expect("generated dataset");
        let comma = if index + 1 == datasets.len() { "" } else { "," };
        writeln!(
            truth,
            "    {{\"id\": {}, \"file\": {}, \"kind\": {}, \"contig_count\": {}, \"sha256\": {}, \"notes\": {}}}{}",
            json_string(dataset.id),
            json_string(dataset.file),
            json_string(dataset.kind),
            dataset.contigs.len(),
            json_string(&hex_digest(&bytes)),
            json_string(dataset.notes),
            comma
        )
        .unwrap();
    }
    truth.push_str("  ],\n  \"cases\": [\n");
    for (index, case) in cases.iter().enumerate() {
        let comma = if index + 1 == cases.len() { "" } else { "," };
        let contigs = json_strings(case.contigs);
        let intervals = json_intervals(case.intervals);
        let contig_identity = json_contig_identity(case, datasets);
        let query = queries
            .iter()
            .find(|query| query.id == case.query)
            .expect("case query exists");
        writeln!(
            truth,
            "    {{\"id\": {}, \"query\": {}, \"query_kind\": {}, \"topology_requested\": {}, \"dataset\": {}, \"relation\": {}, \"contigs\": [{}], \"query_intervals\": [{}], \"expected\": {{\"status\": {}, \"coordinate_model\": {}, \"topology_evidence\": {}, \"wrap_compatible\": {}, \"unique_query_coverage_bases\": {}, \"common_sequence\": {}, \"repeat_label\": {}, \"independent_contig_identity\": true, \"contig_identity\": [{}]}}, \"notes\": {}}}{}",
            json_string(case.id),
            json_string(case.query),
            json_string(query.kind),
            json_string(query.topology),
            json_string(case.dataset),
            json_string(case.relation),
            contigs,
            intervals,
            json_string(case.status),
            json_string(case.coordinate_model),
            json_string(case.topology_evidence),
            case.wrap_compatible,
            case.unique_coverage,
            case.common_sequence,
            json_string(case.repeat_label),
            contig_identity,
            json_string(case.notes),
            comma
        )
        .unwrap();
    }
    truth.push_str("  ]\n}\n");
    truth
}

fn generate(directory: &Path) {
    fs::create_dir_all(directory)
        .unwrap_or_else(|error| panic!("failed to create {}: {error}", directory.display()));
    let (queries, datasets) = build_inputs();
    for query in &queries {
        write_file(directory, query.file, &fasta_record(query_record_id(query), &query.sequence));
    }
    for dataset in &datasets {
        let records: Vec<(&str, &str)> = dataset
            .contigs
            .iter()
            .map(|contig| (contig.id, contig.sequence.as_str()))
            .collect();
        let contents = if dataset.id == "malformed" {
            "sequence-before-header\nACGTACGT\n".to_string()
        } else {
            fasta_records(&records)
        };
        write_file(directory, dataset.file, &contents);
    }
    let cases = build_cases();
    write_file(directory, "truth.json", &render_truth(directory, &queries, &datasets, &cases));

    let mut checksums = String::new();
    let mut files: Vec<&str> = queries.iter().map(|query| query.file).collect();
    files.extend(datasets.iter().map(|dataset| dataset.file));
    files.sort_unstable();
    for file in files {
        let bytes = fs::read(directory.join(file)).expect("generated fixture file");
        writeln!(checksums, "{}  {}", hex_digest(&bytes), file).unwrap();
    }
    checksums.push_str("# The checksum file intentionally does not include itself.\n");
    write_file(directory, "CHECKSUMS.sha256", &checksums);
    write_file(
        directory,
        "README.md",
        "# Mobile genetic element trace fixtures\n\nThese small deterministic FASTA files cover query-centered trace behavior for plasmids, phages, synthetic constructs, and unknown sequence classes. `truth.json` records construction coordinates and exact input checksums independently of Jam output. `CHECKSUMS.sha256` binds the generated FASTA inputs.\n\nCases include linear and circular coordinate handling, auto topology with origin rotations, reverse-complement fragments, terminal direct and inverted repeats, overlapping and duplicate contigs, shared transposons and resistance cassettes, integrated fragments, unrelated chromosome negatives, empty input, and malformed FASTA. A wrapped coordinate model is a coordinate convention; it is not biological topology proof.\n\nRegenerate with `rustc tools/generate_mge_fixtures.rs -o tools/out/generate-mge-fixtures` followed by `tools/out/generate-mge-fixtures tests/data/trace_mge`.\n",
    );
}

fn main() {
    let mut arguments = env::args_os();
    let program = arguments.next().unwrap_or_else(|| "generate_mge_fixtures".into());
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
