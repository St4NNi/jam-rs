//! The screen agrees with a forward search: the assembly as the index entry and the plasmid as
//! the query give the same query coverage, identity and contigs.

use super::tests::{diverged, index_of, plasmids, random_sequence, reverse_complement, screen};
use super::*;
use crate::index::Shard;
use crate::loss::{LossBudget, MissModel};
use crate::search::SearchControl;
use crate::search::align_stage::{AlignStageWorkspace, PreparedQuery, align_entry};
use crate::search::bound::CoverageBound;
use crate::search::locate_stage::{LocateOutcome, LocateRequest, LocateWorkspace, locate_entry};
use crate::search::query::QueryKmers;
use crate::search::sketch_stage::{SketchLimits, sketch_stage};
use crate::seeds::Seeder;
use std::sync::atomic::AtomicBool;

/// The decision of a forward search of `query` (linear) in the only entry of `shard`: covered
/// bases, identity and the contigs of the mosaic, or `None` when it does not qualify.
fn forward(shard: &Shard, query: &[u8], config: &ScreenConfig) -> Option<(u64, f64, Vec<u32>)> {
    let stop = AtomicBool::new(false);
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    let parameters = shard.manifest().parameters;
    let seeder = Seeder::new(parameters.scheme).unwrap();
    let model = MissModel::with_layouts(parameters.scheme, 20).unwrap();
    let rule = &config.rule;
    let length = query.len() as u64;
    let kmers = QueryKmers::new(&seeder, query, false);
    let required = rule.min_query_coverage.required_bases(length);
    let bound = CoverageBound::new(
        &model,
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
fn the_screen_and_a_forward_search_report_the_same_coverage_and_identity() {
    let directory = tempfile::tempdir().unwrap();
    let plasmid = random_sequence(9_000, 101);
    let gene = random_sequence(1_500, 102);
    let shard = plasmids(directory.path(), &[("plasmid", &plasmid), ("gene", &gene)]);
    // A relative at about 97.5% identity with small gaps, in three contigs, one reversed, and a
    // chromosome contig with the gene of the other entry.
    let mut relative = diverged(&plasmid, 40, 103);
    relative.drain(4_000..4_004);
    relative.splice(6_000..6_000, *b"ACG");
    let chromosome = [
        &random_sequence(6_000, 104)[..],
        &gene,
        &random_sequence(6_000, 105),
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
