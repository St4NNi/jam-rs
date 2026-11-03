//! The settings of stage 3, the memory limits of one entry, and how they are raised for an entry
//! that passed one.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::{AlignStageError, KEPT_PREPARED_BYTES};
use crate::align::{
    AlignmentConfig, AnchorConfig, BandConfig, CHAINING_BYTES, ChainConfig, MAX_EXAMINED,
};

/// How stage 3 finds anchors, chains them, aligns them and widens alignment windows, and how
/// much memory one entry may use (window bases of a contig, prepared windows, fragments).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AlignStageConfig {
    /// Anchor length a (15), lookup stride and entropy limit.
    pub anchors: AnchorConfig,
    /// Chaining limits: the largest gap, the fewest anchors and the lowest score.
    pub chains: ChainConfig,
    /// Alignment scores (2, -3, -5, -1), the resident cell limit and the drift limit (1,024).
    pub band: BandConfig,
    /// Query bases an alignment may reach beyond the outermost anchors of its chain. A window is
    /// widened while a chain lies within this many bases (plus the chaining gap) of a cut edge.
    pub extension: u32,
    /// Bases a window first grows by on each side; the step doubles every widening round, so
    /// the bases read again stay within twice the final windows.
    pub widening: u64,
    /// A chain is not aligned when its anchors are too sparse for its stretch to have this
    /// identity (0 aligns every chain).
    pub min_chain_identity: f64,
    /// The most window bases of one contig held at once. Windows are read in full, so this
    /// bounds the bases, anchors and chains of one widening round.
    pub max_window_bases: u64,
    /// The most bytes of the prepared windows of one contig: their bases, chains and walks.
    pub max_prepared_bytes: usize,
    /// The most bytes of fragments that one entry collects before its decision.
    pub max_fragment_bytes: usize,
    /// The most bytes of the dense walk starts of one window, and of its chains when the window
    /// is a copy of a periodic query.
    pub max_dense_bytes: usize,
    /// The most predecessors that chaining one window examines: work, not memory.
    pub max_examined: u64,
    /// The memory that an entry may hold when it passed a limit and is aligned again alone: its
    /// limits are raised only as far as `footprint` stays within it (0: never raised).
    pub memory_budget: u64,
}

/// A limit of stage 3 that an entry can pass, and that is raised when the entry is aligned again.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Limit {
    /// `max_window_bases`.
    WindowBases,
    /// `max_prepared_bytes`.
    PreparedBytes,
    /// `max_fragment_bytes`.
    FragmentBytes,
    /// The resident cells of one alignment, `band.alignment.max_cells`.
    Cells,
    /// `max_dense_bytes`.
    DenseBytes,
    /// `max_examined`.
    Examined,
}

/// How much higher a limit gets each time an entry passes it.
const RAISE: u64 = 4;
/// The most predecessors examined in one window after raising, as a multiple of the default.
const MAX_EXAMINED_RAISE: u64 = 64;
/// Bytes held per resident cell of an alignment, counted high.
const CELL_BYTES: u64 = 4;

impl Default for AlignStageConfig {
    fn default() -> Self {
        Self {
            anchors: AnchorConfig::default(),
            chains: ChainConfig::default(),
            band: BandConfig {
                alignment: AlignmentConfig {
                    max_cells: 1 << 22,
                    ..AlignmentConfig::default()
                },
                ..BandConfig::default()
            },
            extension: 2_000,
            widening: 8_192,
            min_chain_identity: 0.0,
            max_window_bases: 32 << 20,
            max_prepared_bytes: 64 << 20,
            max_fragment_bytes: 512 << 20,
            max_dense_bytes: CHAINING_BYTES,
            max_examined: MAX_EXAMINED,
            memory_budget: 4 << 30,
        }
    }
}

impl AlignStageConfig {
    /// Checks the anchors, the scores and the limits.
    pub fn validate(&self) -> Result<(), AlignStageError> {
        self.anchors.validate()?;
        self.chains.validate()?;
        self.band.alignment.validate()?;
        let identity = self.min_chain_identity;
        if !identity.is_finite() || !(0.0..=1.0).contains(&identity) {
            return Err(AlignStageError::InvalidConfig("chain identity"));
        }
        if self.extension == 0 || self.widening == 0 || self.band.max_drift >= 1 << 20 {
            return Err(AlignStageError::InvalidConfig(
                "extension, widening or drift",
            ));
        }
        if self.max_window_bases == 0
            || self.max_prepared_bytes == 0
            || self.max_fragment_bytes == 0
            || self.max_dense_bytes == 0
            || self.max_examined == 0
        {
            return Err(AlignStageError::InvalidConfig("memory limits"));
        }
        Ok(())
    }

    /// The bytes that one entry may hold under these limits when it is aligned on `threads`
    /// threads: its fragments, the prepared windows it keeps, and per thread the windows of a
    /// contig (two bytes per base, for both strands), its prepared windows, the dense walks and
    /// the anchors of a window, and the cells of an alignment.
    pub fn footprint(&self, threads: usize) -> u64 {
        let threads = threads.max(1) as u64;
        let per_thread = (2 * self.max_window_bases)
            .saturating_add(self.max_prepared_bytes as u64)
            .saturating_add(self.max_dense_bytes as u64)
            .saturating_add(CHAINING_BYTES as u64)
            .saturating_add(CELL_BYTES.saturating_mul(self.band.alignment.max_cells as u64));
        (self.max_fragment_bytes as u64)
            .saturating_add(KEPT_PREPARED_BYTES as u64)
            .saturating_add(threads.saturating_mul(per_thread))
    }

    /// The settings to align an entry again that passed `limit` and needed `needed` (when known):
    /// that limit four times higher, or `needed` when that is more, but only as high as the
    /// memory budget holds the footprint on `threads` threads. `None` when it cannot be raised.
    pub fn raised(&self, limit: Limit, needed: Option<u64>, threads: usize) -> Option<Self> {
        let per_thread = threads.max(1) as u64;
        // The bytes that one more unit of the limit adds to the footprint.
        let weight = match limit {
            Limit::WindowBases => 2 * per_thread,
            Limit::PreparedBytes | Limit::DenseBytes => per_thread,
            Limit::FragmentBytes => 1,
            Limit::Cells => CELL_BYTES * per_thread,
            Limit::Examined => 0,
        };
        let value = self.limit_value(limit);
        let highest = match weight {
            0 => MAX_EXAMINED.saturating_mul(MAX_EXAMINED_RAISE),
            weight => {
                let others = (self.footprint(threads)).saturating_sub(value.saturating_mul(weight));
                self.memory_budget.saturating_sub(others) / weight
            }
        };
        let raised = (value.saturating_mul(RAISE).max(needed.unwrap_or(0))).min(highest);
        (raised > value).then(|| self.with_limit(limit, raised))
    }

    fn limit_value(&self, limit: Limit) -> u64 {
        match limit {
            Limit::WindowBases => self.max_window_bases,
            Limit::PreparedBytes => self.max_prepared_bytes as u64,
            Limit::FragmentBytes => self.max_fragment_bytes as u64,
            Limit::Cells => self.band.alignment.max_cells as u64,
            Limit::DenseBytes => self.max_dense_bytes as u64,
            Limit::Examined => self.max_examined,
        }
    }

    fn with_limit(&self, limit: Limit, value: u64) -> Self {
        let mut config = *self;
        let bytes = usize::try_from(value).unwrap_or(usize::MAX);
        match limit {
            Limit::WindowBases => config.max_window_bases = value,
            Limit::PreparedBytes => config.max_prepared_bytes = bytes,
            Limit::FragmentBytes => config.max_fragment_bytes = bytes,
            Limit::Cells => config.band.alignment.max_cells = bytes,
            Limit::DenseBytes => config.max_dense_bytes = bytes,
            Limit::Examined => config.max_examined = value,
        }
        config
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn documented_defaults() {
        let config = AlignStageConfig::default();
        config.validate().unwrap();
        assert_eq!(config.anchors.a, 15);
        assert_eq!((config.extension, config.band.max_drift), (2_000, 1_024));
        let scores = config.band.alignment;
        let terms = (
            scores.match_score,
            scores.mismatch_score,
            scores.gap_open_score,
            scores.gap_extend_score,
        );
        assert_eq!(terms, (2, -3, -5, -1));
    }

    #[test]
    fn footprint_of_defaults() {
        // 512 MiB of fragments and 256 MiB of kept windows, then per thread 64 MiB of window
        // bases on both strands, 64 MiB of prepared windows, 64 MiB of dense walks, 64 MiB of
        // anchors and 16 MiB of alignment cells.
        let config = AlignStageConfig::default();
        assert_eq!(config.footprint(4), (768 + 4 * 272) << 20);
        assert_eq!(config.footprint(64), (768 + 64 * 272) << 20);
    }

    #[test]
    fn limits_raised_within_budget() {
        let config = AlignStageConfig {
            memory_budget: 4 << 30,
            ..AlignStageConfig::default()
        };
        // Four times higher, or what the entry needed when that is more, within the budget.
        let raised = config
            .raised(Limit::FragmentBytes, Some(1 << 20), 4)
            .unwrap();
        assert_eq!(raised.max_fragment_bytes, 2 << 30);
        let raised = config
            .raised(Limit::FragmentBytes, Some(3 << 30), 4)
            .unwrap();
        assert_eq!(raised.max_fragment_bytes, (4096 - 1344) << 20);
        assert_eq!(raised.footprint(4), config.memory_budget);
        assert!(raised.raised(Limit::FragmentBytes, None, 4).is_none());
        // A limit that each thread may hold weighs once per thread.
        let prepared = config.raised(Limit::PreparedBytes, None, 4).unwrap();
        assert_eq!(prepared.max_prepared_bytes, 256 << 20);
        let prepared = config.raised(Limit::PreparedBytes, None, 64);
        assert!(prepared.is_none());
        // Chaining work is not memory: it rises up to 64 times its default.
        let mut examined = config;
        for _ in 0..3 {
            examined = examined.raised(Limit::Examined, None, 4).unwrap();
        }
        assert_eq!(examined.max_examined, 64 * MAX_EXAMINED);
        assert!(examined.raised(Limit::Examined, None, 4).is_none());
        // Only the limit passed changes.
        let dense = config.raised(Limit::DenseBytes, None, 4).unwrap();
        let expected = AlignStageConfig {
            max_dense_bytes: 4 * CHAINING_BYTES,
            ..config
        };
        assert_eq!(dense, expected);
        let none = AlignStageConfig {
            memory_budget: 0,
            ..config
        };
        assert!(none.raised(Limit::FragmentBytes, None, 4).is_none());
    }

    #[test]
    fn invalid_settings_rejected() {
        let bad = [
            AlignStageConfig {
                min_chain_identity: 1.5,
                ..AlignStageConfig::default()
            },
            AlignStageConfig {
                widening: 0,
                ..AlignStageConfig::default()
            },
            AlignStageConfig {
                max_window_bases: 0,
                ..AlignStageConfig::default()
            },
        ];
        assert!(bad.iter().all(|config| config.validate().is_err()));
    }
}
