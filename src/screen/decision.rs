//! The decision of each entry contig over the contigs of an assembly: the contig rule, the
//! closest margin, and the entry records.

use super::ScreenError;
use super::batch::AssemblyContig;
use super::fragments::EntryContig;
use super::records::{ScreenContig, ScreenEntry};
use crate::decide::{DecisionRule, Fragment};
use std::collections::{HashMap, HashSet};

/// How an entry contig aligns to one assembly contig.
#[derive(Clone, Copy, Debug)]
struct ContigScore {
    identity: f64,
    /// The aligned share of the assembly contig times the identity.
    score: f64,
}

/// The contig rule and the closest margin applied to the fragments of every entry contig.
struct ContigScores {
    /// The assembly contigs that count for each entry contig; `None` when every contig counts.
    allowed: Option<HashMap<EntryContig, HashSet<u32>>>,
    scores: HashMap<(EntryContig, u32), ContigScore>,
    /// The best score any entry contig reaches for an assembly contig.
    best: HashMap<u32, f64>,
    /// For how many entry contigs each assembly contig counts.
    entries: HashMap<u32, u32>,
}

/// The bases in the union of half-open spans.
fn union_length(mut spans: Vec<(u64, u64)>) -> u64 {
    spans.sort_unstable();
    let (mut total, mut reached) = (0, 0);
    for (start, end) in spans {
        let start = start.max(reached);
        if end > start {
            total += end - start;
            reached = end;
        }
    }
    total
}

/// The bases of each assembly contig that the fragments align, counted once.
fn aligned_bases<'a>(fragments: impl Iterator<Item = &'a Fragment>) -> HashMap<u32, u64> {
    let mut spans: HashMap<u32, Vec<(u64, u64)>> = HashMap::new();
    for fragment in fragments {
        let interval = fragment.alignment.contig_interval;
        spans
            .entry(fragment.contig_index)
            .or_default()
            .push((interval.start, interval.end));
    }
    spans
        .into_iter()
        .map(|(contig, spans)| (contig, union_length(spans)))
        .collect()
}

/// An assembly contig counts for an entry contig when enough of its bases align to it with
/// enough identity (the contig rule) and, with a margin, when no other entry contig is closer to
/// it by more than the margin.
fn contig_scores(
    found: &HashMap<EntryContig, Vec<Fragment>>,
    contigs: &[AssemblyContig],
    rule: &DecisionRule,
    closest_margin: Option<f64>,
) -> ContigScores {
    let mut scores = HashMap::new();
    let mut best: HashMap<u32, f64> = HashMap::new();
    for (&key, fragments) in found {
        let aligned = aligned_bases(fragments.iter());
        let mut pairing: HashMap<u32, (u64, u64)> = HashMap::new();
        for fragment in fragments {
            let counts = fragment.alignment.counts();
            let total = pairing.entry(fragment.contig_index).or_default();
            total.0 += counts.matches;
            total.1 += counts.columns();
        }
        for (contig, (matches, columns)) in pairing {
            let length = contigs[contig as usize].length as f64;
            let aligned = aligned[&contig] as f64;
            let identity = matches as f64 / columns.max(1) as f64;
            if aligned < rule.min_contig_coverage * length || identity < rule.min_contig_identity {
                continue;
            }
            let score = aligned / length * identity;
            scores.insert((key, contig), ContigScore { identity, score });
            let top = best.entry(contig).or_insert(score);
            *top = top.max(score);
        }
    }
    let restricted = rule.min_contig_coverage > 0.0
        || rule.min_contig_identity > 0.0
        || closest_margin.is_some();
    let allowed = restricted.then(|| {
        let margin = closest_margin.unwrap_or(f64::INFINITY);
        let mut allowed: HashMap<EntryContig, HashSet<u32>> = HashMap::new();
        for (&(key, contig), found) in &scores {
            if found.score >= best[&contig] - margin {
                allowed.entry(key).or_default().insert(contig);
            }
        }
        allowed
    });
    let mut entries: HashMap<u32, u32> = HashMap::new();
    for &(key, contig) in scores.keys() {
        let counts = allowed
            .as_ref()
            .is_none_or(|allowed| allowed.get(&key).is_some_and(|set| set.contains(&contig)));
        if counts {
            *entries.entry(contig).or_default() += 1;
        }
    }
    ContigScores {
        allowed,
        scores,
        best,
        entries,
    }
}

/// The name and length of an entry contig.
pub(super) type EntryContigName<'a> =
    dyn Fn(EntryContig) -> Result<(String, u64), ScreenError> + 'a;

/// Decides every entry contig from its fragments over the contigs of the assembly and returns
/// those that are reported, by entry and contig.
pub(super) fn decide_entries(
    mut found: HashMap<EntryContig, Vec<Fragment>>,
    contigs: &[AssemblyContig],
    rule: &DecisionRule,
    closest_margin: Option<f64>,
    assembly: &str,
    entry_contig: &EntryContigName<'_>,
) -> Result<Vec<ScreenEntry>, ScreenError> {
    let scores = contig_scores(&found, contigs, rule, closest_margin);
    let mut keys: Vec<EntryContig> = found.keys().copied().collect();
    keys.sort_unstable();
    let mut reported = Vec::new();
    for key in keys {
        let mut fragments = found.remove(&key).unwrap_or_default();
        if let Some(allowed) = &scores.allowed {
            let allowed = allowed.get(&key);
            fragments
                .retain(|fragment| allowed.is_some_and(|set| set.contains(&fragment.contig_index)));
            if fragments.is_empty() {
                continue;
            }
        }
        let (name, length) = entry_contig(key)?;
        let decision = rule.decide(length, fragments, |contig| contigs[contig as usize].length)?;
        if !decision.qualifies {
            continue;
        }
        let mosaic = decision.mosaic;
        let aligned = aligned_bases(
            mosaic
                .primary
                .iter()
                .map(|selected| &selected.fragment)
                .chain(&mosaic.alternatives),
        );
        let mut ids: Vec<u32> = aligned.keys().copied().collect();
        ids.sort_unstable();
        let reported_contigs = ids
            .into_iter()
            .map(|id| {
                let score = scores.scores.get(&(key, id));
                ScreenContig {
                    id,
                    name: contigs[id as usize].name.clone(),
                    length: contigs[id as usize].length,
                    aligned: aligned[&id],
                    identity: score.map_or(0.0, |score| score.identity),
                    score: score.map_or(0.0, |score| score.score),
                    best_score: scores.best.get(&id).copied().unwrap_or(0.0),
                    entries: scores.entries.get(&id).copied().unwrap_or(0),
                }
            })
            .collect();
        let counts = decision.counts;
        reported.push(ScreenEntry {
            kind: "entry",
            entry: name,
            entry_length: length,
            assembly: assembly.to_string(),
            covered_bases: mosaic.covered_bases,
            query_coverage: mosaic.covered_bases as f64 / length.max(1) as f64,
            identity: counts.identity(),
            matches: counts.matches,
            substitutions: counts.substitutions,
            insertions: counts.insertions,
            deletions: counts.deletions,
            contigs: reported_contigs,
            primary: mosaic.primary,
            alternatives: mosaic.alternatives,
            gaps: mosaic.gaps,
        });
    }
    Ok(reported)
}
