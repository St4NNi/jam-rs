//! Deterministic candidate-routing matrix for the witness-selection paths.
//!
//! This executable is an engineering measurement helper.  It reports exact
//! packed-key candidate recovery independently from witness counts and density;
//! it does not infer biological sensitivity from the routing result.

use std::collections::{BTreeMap, BTreeSet};
use std::error::Error;

use jam_rs::router::postings::{InMemoryPostingSource, PostingInput};
use jam_rs::router::search::{CandidateRouter, RouterSearchConfig, TieredQueryWitness};
use jam_rs::router::witness::{
    DEFAULT_QUERY_WINDOW_SIZES, FixedBottomWitness, extract_nested_witnesses,
    fixed_bottom_k_witnesses, local_bottom_r_witnesses,
};
use jam_rs::router::{HashAlgorithmId, WITNESS_K, WitnessKey, WitnessScheme};
use jam_rs::trace::contracts::WitnessPlanningRequest;
use jam_rs::trace::query_plan::plan_witness_tier;
use serde_json::{Value, json};

const BASE_SCALES: [u32; 5] = [10, 20, 25, 50, 100];
const QUERY_SCALES: [u32; 5] = [20, 50, 100, 200, 500];
const TRACE_LENGTHS: [usize; 6] = [80, 100, 160, 250, 500, 1_000];
const IDENTITIES: [u32; 5] = [100, 99, 97, 95, 90];
const FIXED_BOTTOM_K: [usize; 3] = [16, 32, 64];
const LOCAL_R: [usize; 3] = [1, 2, 4];
const POSTING_WINDOW: u32 = 256;
const BACKGROUND_LENGTH: usize = 64_000;
const PLASMID_LENGTH: usize = 12_000;

fn main() -> Result<(), Box<dyn Error>> {
    let plasmid = deterministic_dna(PLASMID_LENGTH, 0x0a11_ce55_1234_5678);
    let (names, metagenomes) = collection(&plasmid);

    let nested_sources = BASE_SCALES
        .iter()
        .copied()
        .map(|base_scale| {
            let scheme = scheme_for_base(base_scale);
            let source = build_nested_source(&names, &metagenomes, &scheme)?;
            Ok::<_, Box<dyn Error>>((base_scale, scheme, source))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let fixed_sources = FIXED_BOTTOM_K
        .iter()
        .copied()
        .map(|bottom_k| {
            let tier = fixed_tier(bottom_k);
            let source = build_fixed_source(&names, &metagenomes, bottom_k)?;
            Ok::<_, Box<dyn Error>>((bottom_k, tier, source))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let local_sources = DEFAULT_QUERY_WINDOW_SIZES
        .iter()
        .copied()
        .flat_map(|window_size| LOCAL_R.iter().copied().map(move |r| (window_size, r)))
        .map(|(window_size, r)| {
            let tier = local_tier(window_size, r);
            let source = build_local_source(&names, &metagenomes, window_size, r)?;
            Ok::<_, Box<dyn Error>>((window_size, r, tier, source))
        })
        .collect::<Result<Vec<_>, _>>()?;

    let mut records = Vec::new();
    for &trace_length in &TRACE_LENGTHS {
        for &identity_percent in &IDENTITIES {
            for layout in [
                Layout::Single,
                Layout::CircularOriginCrossing,
                Layout::MultiFragment,
            ] {
                let base_trace = layout.make_trace(&plasmid, trace_length);
                for orientation in [Orientation::Forward, Orientation::ReverseComplement] {
                    let oriented_trace = orientation.apply(&base_trace);
                    let query = mutate_to_identity(&oriented_trace, identity_percent);
                    let case_id = format!(
                        "{}-{}bp-{}pct-{}",
                        layout.as_str(),
                        trace_length,
                        identity_percent,
                        orientation.as_str()
                    );

                    for (base_scale, scheme, source) in &nested_sources {
                        let nested = extract_nested_witnesses(&query, scheme, POSTING_WINDOW)?;
                        let plan = plan_witness_tier(
                            WitnessPlanningRequest {
                                min_trace_length: trace_length as u64,
                                target_identity: f64::from(identity_percent) / 100.0,
                                max_zero_witness_risk: 0.01,
                                strict: false,
                            },
                            scheme,
                        )?;
                        for scale in &scheme.available_scales {
                            let witnesses = nested.at_scale(*scale)?;
                            records.push(run_record(
                                &case_id,
                                "nested_k21",
                                json!({
                                    "base_scale": base_scale,
                                    "tier": scale,
                                    "available_tiers": scheme.available_scales,
                                    "planner_selected_tier": plan.selected_witness_tier,
                                    "planner_estimated_zero_witness_risk": plan.estimated_zero_witness_risk,
                                    "planner_warning": plan.warning,
                                }),
                                &witnesses,
                                nested.eligible_window_count(),
                                *scale,
                                source,
                            )?);
                        }
                    }

                    for (bottom_k, tier, source) in &fixed_sources {
                        let selected = fixed_bottom_k_witnesses(&query, *bottom_k, POSTING_WINDOW)?;
                        let witnesses = selected
                            .iter()
                            .map(FixedBottomWitness::query_witness)
                            .collect::<Vec<_>>();
                        records.push(run_record(
                            &case_id,
                            "fixed_bottom_k",
                            json!({ "bottom_k": bottom_k, "tier": tier }),
                            &witnesses,
                            query_window_count(query.len(), POSTING_WINDOW),
                            *tier,
                            source,
                        )?);
                    }

                    for (window_size, r, tier, source) in &local_sources {
                        let selected = local_bottom_r_witnesses(&query, *window_size, *r)?;
                        let witnesses = selected.query_witnesses();
                        records.push(run_record(
                            &case_id,
                            "local_bottom_r",
                            json!({ "window_size": window_size, "r": r, "tier": tier }),
                            &witnesses,
                            selected.eligible_window_count,
                            *tier,
                            source,
                        )?);
                    }
                }
            }
        }
    }

    let output = json!({
        "schema": "jam_witness_matrix_v1",
        "claims_limited": true,
        "hash_identity": "jamhash_u64_v1",
        "exact_match_rule": "packed canonical k=21 equality is required alongside jamhash",
        "k": WITNESS_K,
        "base_scales": BASE_SCALES,
        "query_scales": QUERY_SCALES,
        "trace_lengths": TRACE_LENGTHS,
        "identities_percent": IDENTITIES,
        "fixed_bottom_k": FIXED_BOTTOM_K,
        "local_windows": DEFAULT_QUERY_WINDOW_SIZES,
        "local_r": LOCAL_R,
        "posting_window": POSTING_WINDOW,
        "collection": {
            "metagenomes": names,
            "background_bases_per_negative": BACKGROUND_LENGTH,
            "positive_bases": metagenomes[0].len(),
            "embedded_plasmid_bases": PLASMID_LENGTH,
            "construction": "deterministic xorshift DNA backgrounds with one embedded plasmid",
        },
        "records": records,
    });
    println!("{}", serde_json::to_string_pretty(&output)?);
    Ok(())
}

#[derive(Clone, Copy)]
enum Layout {
    Single,
    CircularOriginCrossing,
    MultiFragment,
}

impl Layout {
    const fn as_str(self) -> &'static str {
        match self {
            Self::Single => "single_fragment",
            Self::CircularOriginCrossing => "circular_origin_crossing",
            Self::MultiFragment => "multi_fragment",
        }
    }

    fn make_trace(self, plasmid: &[u8], length: usize) -> Vec<u8> {
        match self {
            Self::Single => plasmid[1_500..1_500 + length].to_vec(),
            Self::CircularOriginCrossing => {
                let tail = length / 2;
                let mut result = plasmid[plasmid.len() - tail..].to_vec();
                result.extend_from_slice(&plasmid[..length - tail]);
                result
            }
            Self::MultiFragment => {
                let gap = 16.min(length / 4).max(1);
                let fragment_length = length - gap;
                let first = fragment_length.div_ceil(2);
                let second = fragment_length - first;
                let mut result = plasmid[2_500..2_500 + first].to_vec();
                result.extend(std::iter::repeat_n(b'N', gap));
                result.extend_from_slice(&plasmid[8_000..8_000 + second]);
                result
            }
        }
    }
}

#[derive(Clone, Copy)]
enum Orientation {
    Forward,
    ReverseComplement,
}

impl Orientation {
    const fn as_str(self) -> &'static str {
        match self {
            Self::Forward => "forward",
            Self::ReverseComplement => "reverse_complement",
        }
    }

    fn apply(self, sequence: &[u8]) -> Vec<u8> {
        match self {
            Self::Forward => sequence.to_vec(),
            Self::ReverseComplement => sequence
                .iter()
                .rev()
                .map(|base| match base {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'T' => b'A',
                    _ => b'N',
                })
                .collect(),
        }
    }
}

fn collection(plasmid: &[u8]) -> (Vec<String>, Vec<Vec<u8>>) {
    let names = vec![
        "positive".to_string(),
        "negative_0".to_string(),
        "negative_1".to_string(),
        "negative_2".to_string(),
        "negative_3".to_string(),
    ];
    let mut positive = deterministic_dna(BACKGROUND_LENGTH, 0x0bad_f00d_0000_0001);
    let insert_at = 20_000;
    positive.splice(insert_at..insert_at, plasmid.iter().copied());
    let mut metagenomes = vec![positive];
    for index in 0..4_u64 {
        metagenomes.push(deterministic_dna(
            BACKGROUND_LENGTH,
            0x0bad_f00d_0000_0100 + index,
        ));
    }
    (names, metagenomes)
}

fn deterministic_dna(length: usize, mut state: u64) -> Vec<u8> {
    let alphabet = *b"ACGT";
    (0..length)
        .map(|_| {
            state ^= state << 7;
            state ^= state >> 9;
            state ^= state << 8;
            alphabet[(state as usize) & 3]
        })
        .collect()
}

fn mutate_to_identity(sequence: &[u8], identity_percent: u32) -> Vec<u8> {
    if identity_percent >= 100 {
        return sequence.to_vec();
    }
    let mut result = sequence.to_vec();
    let editable = result
        .iter()
        .enumerate()
        .filter_map(|(index, base)| matches!(base, b'A' | b'C' | b'G' | b'T').then_some(index))
        .collect::<Vec<_>>();
    let mismatch_count =
        ((editable.len() as f64) * (100 - identity_percent) as f64 / 100.0).round() as usize;
    for index in 0..mismatch_count.min(editable.len()) {
        let position = editable[(index * 7919 + 17) % editable.len()];
        result[position] = match result[position] {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            b'T' => b'A',
            base => base,
        };
    }
    result
}

fn scheme_for_base(base_scale: u32) -> WitnessScheme {
    let mut available_scales = vec![base_scale];
    available_scales.extend(QUERY_SCALES.into_iter().filter(|scale| *scale > base_scale));
    available_scales.sort_unstable();
    available_scales.dedup();
    WitnessScheme {
        scheme_id: base_scale,
        k: WITNESS_K,
        base_scale,
        available_scales,
        hash_id: HashAlgorithmId::JamhashU64V1,
        zero_excluded: true,
    }
}

fn build_nested_source(
    names: &[String],
    metagenomes: &[Vec<u8>],
    scheme: &WitnessScheme,
) -> Result<InMemoryPostingSource, Box<dyn Error>> {
    let mut key_documents = BTreeMap::<u32, BTreeMap<WitnessKey, BTreeSet<u32>>>::new();
    for (document_id, sequence) in metagenomes.iter().enumerate() {
        let set = extract_nested_witnesses(sequence, scheme, POSTING_WINDOW)?;
        for witness in set.witnesses {
            for scale in witness.retained_scales {
                key_documents
                    .entry(scale)
                    .or_default()
                    .entry(witness.key)
                    .or_default()
                    .insert(document_id as u32);
            }
        }
    }
    source_from_key_documents(names, key_documents)
}

fn build_fixed_source(
    names: &[String],
    metagenomes: &[Vec<u8>],
    bottom_k: usize,
) -> Result<InMemoryPostingSource, Box<dyn Error>> {
    let tier = fixed_tier(bottom_k);
    let mut key_documents = BTreeMap::<u32, BTreeMap<WitnessKey, BTreeSet<u32>>>::new();
    for (document_id, sequence) in metagenomes.iter().enumerate() {
        for witness in fixed_bottom_k_witnesses(sequence, bottom_k, POSTING_WINDOW)? {
            key_documents
                .entry(tier)
                .or_default()
                .entry(witness.key)
                .or_default()
                .insert(document_id as u32);
        }
    }
    source_from_key_documents(names, key_documents)
}

fn build_local_source(
    names: &[String],
    metagenomes: &[Vec<u8>],
    window_size: u32,
    r: usize,
) -> Result<InMemoryPostingSource, Box<dyn Error>> {
    let tier = local_tier(window_size, r);
    let mut key_documents = BTreeMap::<u32, BTreeMap<WitnessKey, BTreeSet<u32>>>::new();
    for (document_id, sequence) in metagenomes.iter().enumerate() {
        for witness in local_bottom_r_witnesses(sequence, window_size, r)?.witnesses {
            key_documents
                .entry(tier)
                .or_default()
                .entry(witness.key)
                .or_default()
                .insert(document_id as u32);
        }
    }
    source_from_key_documents(names, key_documents)
}

fn source_from_key_documents(
    names: &[String],
    key_documents: BTreeMap<u32, BTreeMap<WitnessKey, BTreeSet<u32>>>,
) -> Result<InMemoryPostingSource, Box<dyn Error>> {
    let mut source = InMemoryPostingSource::new(names.to_vec());
    for (tier, keys) in key_documents {
        for (key, documents) in keys {
            source.insert(PostingInput::new(
                tier,
                key,
                documents.into_iter().collect(),
            ))?;
        }
    }
    Ok(source)
}

fn run_record(
    case_id: &str,
    method: &str,
    method_parameters: Value,
    witnesses: &[jam_rs::router::QueryWitness],
    total_windows: u32,
    tier: u32,
    source: &InMemoryPostingSource,
) -> Result<Value, Box<dyn Error>> {
    let mut config = RouterSearchConfig::default();
    config.top_k = source.document_names().len();
    config.accumulator_capacity = config.top_k.saturating_mul(2).max(config.top_k);
    config.max_shared_witnesses_per_candidate = 256;
    config.max_positional_witnesses_per_candidate = 1;
    config.common_query_window_requirement = 1;
    config.total_query_windows = Some(total_windows.max(1));
    let tiered = witnesses
        .iter()
        .cloned()
        .map(|witness| TieredQueryWitness::new(tier, witness))
        .collect::<Vec<_>>();
    let candidates = CandidateRouter::new(source, config).search(&tiered)?;
    let candidate_ids = candidates
        .iter()
        .map(|candidate| candidate.metagenome_id.clone())
        .collect::<Vec<_>>();
    let mut unique_keys = witnesses
        .iter()
        .map(|witness| witness.key)
        .collect::<Vec<_>>();
    unique_keys.sort_unstable();
    unique_keys.dedup();
    Ok(json!({
        "case_id": case_id,
        "method": method,
        "parameters": method_parameters,
        "tier": tier,
        "exact_witness_occurrences": witnesses.len(),
        "exact_unique_packed_keys": unique_keys.len(),
        "eligible_query_windows": total_windows,
        "witness_density_per_query_window": witnesses.len() as f64 / total_windows.max(1) as f64,
        "candidate_count": candidate_ids.len(),
        "candidate_ids": candidate_ids,
        "candidate_recovered": candidates.iter().any(|candidate| candidate.metagenome_id == "positive"),
        "false_candidate_count": candidates.iter().filter(|candidate| candidate.metagenome_id != "positive").count(),
        "shared_witness_count": candidates.iter().map(|candidate| candidate.window_evidence.total_shared_witness_count).sum::<u32>(),
    }))
}

const fn fixed_tier(bottom_k: usize) -> u32 {
    10_000 + bottom_k as u32
}

const fn local_tier(window_size: u32, r: usize) -> u32 {
    20_000 + window_size * 10 + r as u32
}

const fn query_window_count(sequence_length: usize, window_size: u32) -> u32 {
    if sequence_length == 0 {
        0
    } else {
        ((sequence_length - 1) / window_size as usize + 1) as u32
    }
}
