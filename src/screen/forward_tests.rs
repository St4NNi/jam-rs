//! The screen agrees with a forward search: the assembly as the index entry and the plasmid as
//! the query give the same query coverage, identity and contigs.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::tests::{index_of, plasmids, screen};
use super::*;
use crate::index::Shard;
use crate::loss::LossBudget;
use crate::search::CoverageBound;
use crate::search::QueryKmers;
use crate::search::{AlignStageWorkspace, PreparedQuery, align_entry};
use crate::search::{LocateOutcome, LocateRequest, LocateWorkspace, locate_entry};
use crate::search::{SketchLimits, sketch_stage};
use crate::seeds::{Seeder, reverse_complement};
use crate::test_support::{miss_model, no_stop, random_bases, rng, with_substitutions};

/// The decision of a forward search of `query` (linear) in the only entry of `shard`: covered
/// bases, identity and the contigs of the mosaic, or `None` when it does not qualify.
fn forward(shard: &Shard, query: &[u8], config: &ScreenConfig) -> Option<(u64, f64, Vec<u32>)> {
    let control = no_stop();
    let parameters = shard.manifest().parameters;
    let seeder = Seeder::new(parameters.scheme).unwrap();
    let model = miss_model();
    let rule = &config.rule;
    let length = query.len() as u64;
    let kmers = QueryKmers::new(&seeder, query, false);
    let required = rule.min_query_coverage.required_bases(length);
    let bound = CoverageBound::new(
        model,
        LossBudget::default(),
        rule.min_identity,
        length,
        false,
        required,
    );
    let limits = SketchLimits {
        frequent_seed_entries: 1_000,
        max_seed_hits: 1 << 20,
    };
    let sketched = sketch_stage(shard.sketch(), &kmers, &bound, limits, control).unwrap();
    let possible = sketched.candidates.first()?;
    let request = LocateRequest {
        query: &kmers,
        bound: &bound,
        frequent_seeds: &sketched.frequent_seeds,
        flank: config.flank,
    };
    let mut work = Default::default();
    let outcome = locate_entry(
        shard,
        &request,
        possible,
        &mut LocateWorkspace::default(),
        &mut work,
        control,
    )
    .unwrap();
    let LocateOutcome::Located(located) = outcome else {
        return None;
    };
    let prepared = PreparedQuery::new(query, false, &config.align).unwrap();
    let mut workspace = AlignStageWorkspace::default();
    let aligned = align_entry(
        shard,
        &prepared,
        &located,
        &config.align,
        rule,
        &mut workspace,
        control,
    )
    .unwrap();
    let decision = aligned.decision;
    if !decision.qualifies {
        return None;
    }
    let mut contigs: Vec<u32> = decision
        .mosaic
        .primary
        .iter()
        .map(|selected| selected.fragment.contig_index)
        .collect();
    contigs.sort_unstable();
    contigs.dedup();
    Some((decision.mosaic.covered_bases, decision.identity(), contigs))
}

#[test]
fn screen_matches_search() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_bases(&mut rng(101), 9_000);
    let gene = random_bases(&mut rng(102), 1_500);
    let shard = plasmids(directory.path(), &[("plasmid", &plasmid), ("gene", &gene)]);
    // A relative at about 97.5% identity with small gaps, in three contigs, one reversed, and a
    // chromosome contig with the gene of the other entry.
    let mut relative = with_substitutions(&mut rng(103), &plasmid, 0.975);
    relative.drain(4_000..4_004);
    relative.splice(6_000..6_000, *b"ACG");
    let chromosome = [
        &random_bases(&mut rng(104), 6_000)[..],
        &gene,
        &random_bases(&mut rng(105), 6_000),
    ]
    .concat();
    let contigs = [
        ("one", relative[..3_000].to_vec()),
        ("chromosome", chromosome),
        ("two", reverse_complement(&relative[3_000..6_500])),
        ("three", relative[6_600..].to_vec()),
    ];
    for config in [
        ScreenConfig::default(),
        ScreenConfig::new(DecisionRule {
            min_identity: 0.95,
            min_contig_coverage: 0.5,
            min_contig_identity: 0.95,
            ..DecisionRule::default()
        }),
    ] {
        let (reported, _) = screen(&shard, directory.path(), &contigs, &config);
        let assembly = tempfile::tempdir().unwrap();
        let sequences: Vec<&[u8]> = contigs
            .iter()
            .map(|(_, sequence)| sequence.as_slice())
            .collect();
        let forward_shard = index_of(assembly.path(), &[("assembly", &sequences)]);
        for (query, name) in [(&plasmid, "plasmid"), (&gene, "gene")] {
            let screened = reported.iter().find(|report| report.entry == name);
            let searched = forward(&forward_shard, query, &config);
            assert_eq!(screened.is_some(), searched.is_some(), "{name} {config:?}");
            // The gene alone is reported only without the contig rule.
            let expected = name == "plasmid" || config.rule.min_contig_coverage == 0.0;
            assert_eq!(screened.is_some(), expected, "{name}");
            let (Some(report), Some((covered, identity, contigs))) = (screened, searched) else {
                continue;
            };
            assert!(
                report.covered_bases.abs_diff(covered) <= 10,
                "{name}: {} {covered}",
                report.covered_bases
            );
            assert!(
                (report.identity - identity).abs() < 0.002,
                "{name}: {} {identity}",
                report.identity
            );
            let ids: Vec<u32> = report.contigs.iter().map(|contig| contig.id).collect();
            assert_eq!(ids, contigs, "{name}");
        }
    }
}
