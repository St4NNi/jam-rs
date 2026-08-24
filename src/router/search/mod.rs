//! Bounded, deterministic candidate routing over witness postings.
//!
//! Search consumes [`PostingSource`] rather than a concrete file reader.  A
//! local mmap reader, a range reader, and the in-memory fixture can therefore
//! share exactly the same candidate evidence and handoff behavior.

use crate::router::postings::{PositionPosting, PostingError, PostingHeader, PostingSource};
use crate::router::{
    CandidateWindowEvidence, PositionalWitnessOccurrence, QueryWitness, RoutedCandidate,
    SharedWitness, WitnessClass, WitnessHandoffMode, WitnessKey,
};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use thiserror::Error;

/// One query witness together with the nested tier at which it was searched.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TieredQueryWitness {
    pub tier: u32,
    pub witness: QueryWitness,
}

impl TieredQueryWitness {
    #[must_use]
    pub const fn new(tier: u32, witness: QueryWitness) -> Self {
        Self { tier, witness }
    }
}

/// Controls evidence thresholds, bounded state, and handoff detail.
#[derive(Clone, Debug, PartialEq)]
pub struct RouterSearchConfig {
    pub top_k: usize,
    /// Capacity hint retained for compatibility with the memory planner.  A
    /// two-pass search does not use it to evict summaries, because doing so
    /// would make top-k selection order-dependent.
    pub accumulator_capacity: usize,
    pub max_shared_witnesses_per_candidate: usize,
    pub max_positional_witnesses_per_candidate: usize,
    /// A header above this document-frequency limit is suppressed before its
    /// document payload is decoded.
    pub max_document_frequency: Option<u32>,
    /// Independent occurrence cap for readers whose document frequency and
    /// occurrence count are not identical.
    pub max_posting_count: Option<u32>,
    pub rare_max_document_frequency: u32,
    pub moderately_common_max_document_frequency: u32,
    /// A candidate with no rare evidence must cover this many independent
    /// query windows to remain in the result set.
    pub common_query_window_requirement: u32,
    pub handoff_mode: WitnessHandoffMode,
    pub max_positions_per_witness: usize,
    pub total_query_windows: Option<u32>,
}

impl Default for RouterSearchConfig {
    fn default() -> Self {
        Self {
            top_k: 100,
            accumulator_capacity: 512,
            max_shared_witnesses_per_candidate: 64,
            max_positional_witnesses_per_candidate: 128,
            max_document_frequency: None,
            max_posting_count: None,
            rare_max_document_frequency: 4,
            moderately_common_max_document_frequency: 64,
            common_query_window_requirement: 2,
            handoff_mode: WitnessHandoffMode::SampleOnly,
            max_positions_per_witness: 16,
            total_query_windows: None,
        }
    }
}

impl RouterSearchConfig {
    pub fn validate(&self) -> Result<(), SearchError> {
        if self.top_k == 0 {
            return Err(SearchError::InvalidConfig(
                "top_k must be greater than zero".to_string(),
            ));
        }
        if self.accumulator_capacity < self.top_k {
            return Err(SearchError::InvalidConfig(
                "accumulator_capacity must be at least top_k".to_string(),
            ));
        }
        if self.max_shared_witnesses_per_candidate == 0 {
            return Err(SearchError::InvalidConfig(
                "max_shared_witnesses_per_candidate must be greater than zero".to_string(),
            ));
        }
        if self.max_positional_witnesses_per_candidate == 0 {
            return Err(SearchError::InvalidConfig(
                "max_positional_witnesses_per_candidate must be greater than zero".to_string(),
            ));
        }
        if self.max_positions_per_witness == 0 {
            return Err(SearchError::InvalidConfig(
                "max_positions_per_witness must be greater than zero".to_string(),
            ));
        }
        if self.rare_max_document_frequency > self.moderately_common_max_document_frequency {
            return Err(SearchError::InvalidConfig(
                "rare_max_document_frequency must not exceed moderately_common_max_document_frequency"
                    .to_string(),
            ));
        }
        Ok(())
    }
}

/// A reusable candidate router over any exact witness posting source.
pub struct CandidateRouter<'a, S: PostingSource> {
    source: &'a S,
    config: RouterSearchConfig,
}

impl<'a, S: PostingSource> CandidateRouter<'a, S> {
    #[must_use]
    pub fn new(source: &'a S, config: RouterSearchConfig) -> Self {
        Self { source, config }
    }

    #[must_use]
    pub fn source(&self) -> &'a S {
        self.source
    }

    #[must_use]
    pub fn config(&self) -> &RouterSearchConfig {
        &self.config
    }

    /// Search every supplied tier and return their union.  Scale is inverse
    /// density, so routing starts with the largest (sparsest) tier and then
    /// escalates toward denser tiers.  Pass one keeps only fixed-size
    /// per-metagenome summaries; pass two revisits postings for exact selected
    /// IDs and is the only pass which materializes witness/position details.
    pub fn search(
        &self,
        query_witnesses: &[TieredQueryWitness],
    ) -> Result<Vec<RoutedCandidate>, SearchError> {
        self.config.validate()?;
        let mut witnesses = query_witnesses.to_vec();
        normalize_query_witnesses(&mut witnesses);
        if witnesses.is_empty() {
            return Ok(Vec::new());
        }

        let eligible_windows = query_window_set(&witnesses);
        let total_windows = self
            .config
            .total_query_windows
            .unwrap_or_else(|| u32::try_from(eligible_windows.len()).unwrap_or(u32::MAX));
        let query_units = query_witness_units(&witnesses).max(1) as f64;
        let collection_size = self.source.collection_size();
        let mut summaries = BTreeMap::<u32, CandidateSummary>::new();

        // Pass one: exact counters and only enough distinct windows to decide
        // common-witness admission.  No shared-witness or position vectors
        // are allocated here.
        for tiered in &witnesses {
            let Some(descriptor) =
                posting_descriptor(self.source, tiered, &self.config, collection_size)?
            else {
                continue;
            };
            let mut missing_metagenome = None;
            let source = self.source;
            source.for_each_document_id(tiered.tier, &tiered.witness.key, &mut |document_id| {
                if missing_metagenome.is_some() {
                    return Ok(());
                }
                let Some(metagenome_id) = source.metagenome_id(document_id) else {
                    missing_metagenome = Some(document_id);
                    return Ok(());
                };
                let detail = descriptor.detail(tiered);
                summaries
                    .entry(document_id)
                    .or_insert_with(|| CandidateSummary::new(metagenome_id, tiered.tier))
                    .add(detail, self.config.common_query_window_requirement);
                Ok(())
            })?;
            if let Some(document_id) = missing_metagenome {
                return Err(SearchError::MissingMetagenome(document_id));
            }
        }

        let mut selected_ids = summaries
            .iter()
            .filter(|(_, summary)| summary.qualifies(&self.config))
            .map(|(document_id, _)| *document_id)
            .collect::<Vec<_>>();
        selected_ids.sort_by(|left, right| {
            summary_cmp(&summaries[left], &summaries[right]).then_with(|| left.cmp(right))
        });
        selected_ids.truncate(self.config.top_k);
        let selected_ids = selected_ids.into_iter().collect::<BTreeSet<_>>();
        if selected_ids.is_empty() {
            return Ok(Vec::new());
        }

        // Pass two: collect full details only for the exact selected IDs.
        let mut states = BTreeMap::<u32, CandidateState>::new();
        let mut position_cache = BTreeMap::<(u32, WitnessKey, u32), Vec<PositionPosting>>::new();
        let position_cache_capacity = self
            .config
            .top_k
            .saturating_mul(self.config.max_shared_witnesses_per_candidate)
            .max(1);
        for tiered in &witnesses {
            let Some(descriptor) =
                posting_descriptor(self.source, tiered, &self.config, collection_size)?
            else {
                continue;
            };
            let mut missing_metagenome = None;
            let source = self.source;
            source.for_each_document_id(tiered.tier, &tiered.witness.key, &mut |document_id| {
                if missing_metagenome.is_some() || !selected_ids.contains(&document_id) {
                    return Ok(());
                }
                let Some(metagenome_id) = source.metagenome_id(document_id) else {
                    missing_metagenome = Some(document_id);
                    return Ok(());
                };
                let detail = descriptor.detail(tiered);
                let state = states
                    .entry(document_id)
                    .or_insert_with(|| CandidateState::new(metagenome_id, tiered.tier));
                let Some(shared) = state.add(detail, &self.config) else {
                    return Ok(());
                };
                if !should_handoff_positions(
                    self.config.handoff_mode,
                    descriptor.class,
                    descriptor.header,
                ) {
                    return Ok(());
                }
                let cache_key = (tiered.tier, tiered.witness.key, document_id);
                let positions = if let Some(cached) = position_cache.get(&cache_key) {
                    cached.clone()
                } else {
                    let mut decoded = source.decode_positions(
                        tiered.tier,
                        &tiered.witness.key,
                        document_id,
                        self.config.max_positions_per_witness,
                    )?;
                    decoded.sort_unstable();
                    decoded.dedup();
                    decoded.truncate(self.config.max_positions_per_witness);
                    if position_cache.len() < position_cache_capacity {
                        position_cache.insert(cache_key, decoded.clone());
                    }
                    decoded
                };
                state.add_positions(&shared, detail.scheme_id, &positions, &self.config);
                Ok(())
            })?;
            if let Some(document_id) = missing_metagenome {
                return Err(SearchError::MissingMetagenome(document_id));
            }
        }

        let mut results = states
            .into_values()
            .filter(|state| state.qualifies(&self.config))
            .map(|state| state.into_result(total_windows, query_units, &eligible_windows))
            .collect::<Vec<_>>();
        results.sort_by(routed_candidate_cmp);
        results.truncate(self.config.top_k);
        Ok(results)
    }
}

/// Convenience function for callers which do not need to retain a router.
pub fn search<S: PostingSource>(
    source: &S,
    query_witnesses: &[TieredQueryWitness],
    config: &RouterSearchConfig,
) -> Result<Vec<RoutedCandidate>, SearchError> {
    CandidateRouter::new(source, config.clone()).search(query_witnesses)
}

#[derive(Clone, Copy)]
struct SharedDetail<'a> {
    key: WitnessKey,
    query_position: u64,
    query_window_ids: &'a [u32],
    query_reverse: bool,
    document_frequency: u32,
    tier: u32,
    scheme_id: u32,
    class: WitnessClass,
    weight: f64,
}

#[derive(Clone, Copy)]
struct PostingDescriptor {
    header: PostingHeader,
    class: WitnessClass,
    weight: f64,
}

impl PostingDescriptor {
    fn detail<'a>(&self, tiered: &'a TieredQueryWitness) -> SharedDetail<'a> {
        SharedDetail {
            key: tiered.witness.key,
            query_position: tiered.witness.query_position,
            query_window_ids: &tiered.witness.query_window_ids,
            query_reverse: tiered.witness.query_reverse,
            document_frequency: self.header.document_frequency,
            tier: tiered.tier,
            scheme_id: self.header.scheme_id,
            class: self.class,
            weight: self.weight,
        }
    }
}

struct CandidateState {
    metagenome_id: String,
    rare_shared_witnesses: u32,
    common_shared_witnesses: u32,
    weighted_witness_sum: f64,
    candidate_tier: u32,
    windows: BTreeSet<u32>,
    shared_witnesses: Vec<SharedWitness>,
    positions: Vec<PositionalWitnessOccurrence>,
    position_identities: Vec<PositionIdentity>,
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct PositionIdentity {
    key: WitnessKey,
    query_position: u64,
    query_reverse: bool,
    query_window_id: u32,
    contig_id: u32,
    position: u64,
    reverse: bool,
    scheme_id: u32,
}

impl CandidateState {
    fn new(metagenome_id: &str, tier: u32) -> Self {
        Self {
            metagenome_id: metagenome_id.to_string(),
            rare_shared_witnesses: 0,
            common_shared_witnesses: 0,
            weighted_witness_sum: 0.0,
            candidate_tier: tier,
            windows: BTreeSet::new(),
            shared_witnesses: Vec::new(),
            positions: Vec::new(),
            position_identities: Vec::new(),
        }
    }

    fn add(
        &mut self,
        detail: SharedDetail<'_>,
        config: &RouterSearchConfig,
    ) -> Option<SharedWitness> {
        // This is one exact query occurrence.  Overlapping query windows are
        // support labels only and must not multiply its counters or weight.
        for &window_id in detail.query_window_ids {
            self.windows.insert(window_id);
        }
        self.weighted_witness_sum += detail.weight;
        match detail.class {
            WitnessClass::Rare => self.rare_shared_witnesses += 1,
            WitnessClass::ModeratelyCommon | WitnessClass::CollectionCommon => {
                self.common_shared_witnesses += 1
            }
            WitnessClass::Suppressed => {}
        }
        if self.shared_witnesses.len() >= config.max_shared_witnesses_per_candidate {
            return None;
        }
        let shared = make_shared_witness(&detail);
        self.shared_witnesses.push(shared);
        self.shared_witnesses.sort_by(shared_witness_cmp);
        Some(shared)
    }

    fn add_positions(
        &mut self,
        shared: &SharedWitness,
        scheme_id: u32,
        positions: &[PositionPosting],
        config: &RouterSearchConfig,
    ) {
        for occurrence in positions {
            let identity = PositionIdentity {
                key: shared.key,
                query_position: shared.query_position,
                query_reverse: shared.query_reverse,
                query_window_id: shared.query_window_id,
                contig_id: occurrence.contig_id,
                position: occurrence.position,
                reverse: occurrence.reverse,
                scheme_id,
            };
            if self.position_identities.len() >= config.max_positional_witnesses_per_candidate
                || self.position_identities.contains(&identity)
            {
                continue;
            }
            self.position_identities.push(identity);
            self.positions
                .push(make_positional_occurrence(occurrence, shared, scheme_id));
        }
    }

    fn qualifies(&self, config: &RouterSearchConfig) -> bool {
        self.rare_shared_witnesses > 0
            || u32::try_from(self.windows.len()).unwrap_or(u32::MAX)
                >= config.common_query_window_requirement
    }

    fn into_result(
        self,
        total_windows: u32,
        query_units: f64,
        eligible_windows: &BTreeSet<u32>,
    ) -> RoutedCandidate {
        let matched_query_windows = u32::try_from(self.windows.len()).unwrap_or(u32::MAX);
        let longest_supported_window_run = longest_run(&self.windows);
        let total_shared = self.rare_shared_witnesses + self.common_shared_witnesses;
        RoutedCandidate {
            metagenome_id: self.metagenome_id,
            rare_shared_witnesses: self.rare_shared_witnesses,
            common_shared_witnesses: self.common_shared_witnesses,
            supported_query_windows: matched_query_windows,
            weighted_witness_sum: self.weighted_witness_sum,
            estimated_query_containment: (f64::from(total_shared) / query_units.max(1.0)).min(1.0),
            candidate_tier: self.candidate_tier,
            window_evidence: CandidateWindowEvidence {
                matched_query_windows,
                total_eligible_query_windows: total_windows
                    .max(u32::try_from(eligible_windows.len()).unwrap_or(u32::MAX)),
                longest_supported_window_run,
                rare_witness_count: self.rare_shared_witnesses,
                common_witness_count: self.common_shared_witnesses,
                total_shared_witness_count: total_shared,
            },
            shared_witnesses: self.shared_witnesses,
            positional_witnesses: self.positions,
        }
    }
}

struct CandidateSummary {
    metagenome_id: String,
    rare_shared_witnesses: u32,
    common_shared_witnesses: u32,
    weighted_witness_sum: f64,
    candidate_tier: u32,
    windows: BTreeSet<u32>,
}

impl CandidateSummary {
    fn new(metagenome_id: &str, tier: u32) -> Self {
        Self {
            metagenome_id: metagenome_id.to_string(),
            rare_shared_witnesses: 0,
            common_shared_witnesses: 0,
            weighted_witness_sum: 0.0,
            candidate_tier: tier,
            windows: BTreeSet::new(),
        }
    }

    fn add(&mut self, detail: SharedDetail<'_>, common_window_requirement: u32) {
        self.weighted_witness_sum += detail.weight;
        match detail.class {
            WitnessClass::Rare => self.rare_shared_witnesses += 1,
            WitnessClass::ModeratelyCommon | WitnessClass::CollectionCommon => {
                self.common_shared_witnesses += 1
            }
            WitnessClass::Suppressed => {}
        }
        let window_cap = usize::try_from(common_window_requirement).unwrap_or(usize::MAX);
        if self.windows.len() < window_cap {
            for &window_id in detail.query_window_ids {
                self.windows.insert(window_id);
                if self.windows.len() >= window_cap {
                    break;
                }
            }
        }
    }

    fn qualifies(&self, config: &RouterSearchConfig) -> bool {
        self.rare_shared_witnesses > 0
            || u32::try_from(self.windows.len()).unwrap_or(u32::MAX)
                >= config.common_query_window_requirement
    }
}

fn summary_cmp(left: &CandidateSummary, right: &CandidateSummary) -> Ordering {
    right
        .weighted_witness_sum
        .total_cmp(&left.weighted_witness_sum)
        .then_with(|| right.rare_shared_witnesses.cmp(&left.rare_shared_witnesses))
        .then_with(|| {
            let left_total = left.rare_shared_witnesses + left.common_shared_witnesses;
            let right_total = right.rare_shared_witnesses + right.common_shared_witnesses;
            right_total.cmp(&left_total)
        })
        // Larger scale is the earlier sparse routing tier.
        .then_with(|| right.candidate_tier.cmp(&left.candidate_tier))
        .then_with(|| left.metagenome_id.cmp(&right.metagenome_id))
}

fn routed_candidate_cmp(left: &RoutedCandidate, right: &RoutedCandidate) -> Ordering {
    right
        .weighted_witness_sum
        .total_cmp(&left.weighted_witness_sum)
        .then_with(|| right.rare_shared_witnesses.cmp(&left.rare_shared_witnesses))
        .then_with(|| {
            right
                .window_evidence
                .total_shared_witness_count
                .cmp(&left.window_evidence.total_shared_witness_count)
        })
        .then_with(|| right.candidate_tier.cmp(&left.candidate_tier))
        .then_with(|| left.metagenome_id.cmp(&right.metagenome_id))
}

fn shared_witness_cmp(left: &SharedWitness, right: &SharedWitness) -> Ordering {
    left.witness_tier
        .cmp(&right.witness_tier)
        .then_with(|| left.key.cmp(&right.key))
        .then_with(|| left.query_position.cmp(&right.query_position))
        .then_with(|| left.query_reverse.cmp(&right.query_reverse))
        .then_with(|| left.query_window_id.cmp(&right.query_window_id))
}

fn make_shared_witness(detail: &SharedDetail<'_>) -> SharedWitness {
    SharedWitness {
        key: detail.key,
        query_position: detail.query_position,
        query_reverse: detail.query_reverse,
        query_window_id: detail.query_window_ids.first().copied().unwrap_or(0),
        document_frequency: detail.document_frequency,
        witness_tier: detail.tier,
        witness_class: detail.class,
    }
}

fn make_positional_occurrence(
    occurrence: &PositionPosting,
    shared: &SharedWitness,
    scheme_id: u32,
) -> PositionalWitnessOccurrence {
    PositionalWitnessOccurrence {
        witness: *shared,
        contig_id: occurrence.contig_id,
        position: occurrence.position,
        reverse: occurrence.reverse,
        scheme_id,
    }
}

fn posting_descriptor<S: PostingSource>(
    source: &S,
    tiered: &TieredQueryWitness,
    config: &RouterSearchConfig,
    collection_size: u32,
) -> Result<Option<PostingDescriptor>, SearchError> {
    if tiered.witness.query_window_ids.is_empty() {
        return Ok(None);
    }
    let Some(header) = source.try_header(tiered.tier, &tiered.witness.key)? else {
        return Ok(None);
    };
    // Header checks happen before either pass asks the source for payload
    // bytes.  A repetitive/suppressed key is never evidence by itself.
    if header.suppressed
        || config
            .max_document_frequency
            .is_some_and(|maximum| header.document_frequency > maximum)
        || config
            .max_posting_count
            .is_some_and(|maximum| header.posting_count > maximum)
        || header.is_empty()
    {
        return Ok(None);
    }
    Ok(Some(PostingDescriptor {
        class: classify(
            header,
            config.rare_max_document_frequency,
            config.moderately_common_max_document_frequency,
        ),
        weight: witness_weight(collection_size, header.document_frequency),
        header,
    }))
}

fn classify(header: PostingHeader, rare_max: u32, moderate_max: u32) -> WitnessClass {
    if header.suppressed {
        WitnessClass::Suppressed
    } else if header.common_or_repetitive {
        WitnessClass::CollectionCommon
    } else if header.document_frequency <= rare_max {
        WitnessClass::Rare
    } else if header.document_frequency <= moderate_max {
        WitnessClass::ModeratelyCommon
    } else {
        WitnessClass::CollectionCommon
    }
}

fn should_handoff_positions(
    mode: WitnessHandoffMode,
    class: WitnessClass,
    header: PostingHeader,
) -> bool {
    if !header.position_bearing {
        return false;
    }
    match mode {
        WitnessHandoffMode::SampleOnly => false,
        WitnessHandoffMode::PositionBearing => class != WitnessClass::Suppressed,
        WitnessHandoffMode::Hybrid => class == WitnessClass::Rare,
    }
}

fn witness_weight(collection_size: u32, document_frequency: u32) -> f64 {
    (f64::from(collection_size.saturating_add(1)) / f64::from(document_frequency.saturating_add(1)))
        .ln()
        .max(0.0)
}

fn normalize_query_witnesses(witnesses: &mut Vec<TieredQueryWitness>) {
    let mut merged = BTreeMap::<(WitnessKey, u64, bool), (u32, BTreeSet<u32>)>::new();
    for tiered in witnesses.drain(..) {
        if tiered.tier == 0 || tiered.witness.query_window_ids.is_empty() {
            continue;
        }
        let entry = merged
            .entry((
                tiered.witness.key,
                tiered.witness.query_position,
                tiered.witness.query_reverse,
            ))
            .or_insert_with(|| (tiered.tier, BTreeSet::new()));
        // Scale is inverse density: retain the sparse tier as the first
        // exact lookup, then let denser tiers add new query occurrences.
        entry.0 = entry.0.max(tiered.tier);
        entry.1.extend(tiered.witness.query_window_ids);
    }
    *witnesses = merged
        .into_iter()
        .map(
            |((key, query_position, query_reverse), (tier, windows))| TieredQueryWitness {
                tier,
                witness: QueryWitness {
                    key,
                    query_position,
                    query_reverse,
                    query_window_ids: windows.into_iter().collect(),
                },
            },
        )
        .collect();
    witnesses.sort_by(|left, right| {
        right
            .tier
            .cmp(&left.tier)
            .then_with(|| left.witness.key.cmp(&right.witness.key))
            .then_with(|| {
                left.witness
                    .query_position
                    .cmp(&right.witness.query_position)
            })
            .then_with(|| left.witness.query_reverse.cmp(&right.witness.query_reverse))
            .then_with(|| {
                left.witness
                    .query_window_ids
                    .cmp(&right.witness.query_window_ids)
            })
    });
}

fn query_window_set(witnesses: &[TieredQueryWitness]) -> BTreeSet<u32> {
    witnesses
        .iter()
        .flat_map(|tiered| tiered.witness.query_window_ids.iter().copied())
        .collect()
}

fn query_witness_units(witnesses: &[TieredQueryWitness]) -> usize {
    witnesses
        .iter()
        .map(|tiered| {
            (
                tiered.witness.key,
                tiered.witness.query_position,
                tiered.witness.query_reverse,
            )
        })
        .collect::<BTreeSet<_>>()
        .len()
}

fn longest_run(windows: &BTreeSet<u32>) -> u32 {
    let mut longest = 0_u32;
    let mut current = 0_u32;
    let mut previous: Option<u32> = None;
    for window in windows {
        if previous.is_some_and(|previous| *window == previous.saturating_add(1)) {
            current = current.saturating_add(1);
        } else {
            current = 1;
        }
        longest = longest.max(current);
        previous = Some(*window);
    }
    longest
}

/// Candidate search failures, including source payload errors.
#[derive(Debug, Error, Eq, PartialEq)]
pub enum SearchError {
    #[error("invalid router search configuration: {0}")]
    InvalidConfig(String),
    #[error("metagenome document {0} has no catalog identifier")]
    MissingMetagenome(u32),
    #[error(transparent)]
    Posting(#[from] PostingError),
}
