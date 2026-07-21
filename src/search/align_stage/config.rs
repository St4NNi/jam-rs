//! The settings of stage 3 and the memory limits of one entry.

use super::AlignStageError;
use crate::align::{AlignmentConfig, AnchorConfig, BandConfig};

/// How stage 3 finds anchors, chains them, aligns them and widens alignment windows, and how
/// much memory one entry may use.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AlignStageConfig {
    /// Anchor length a (15), lookup stride, entropy limit and chaining limits.
    pub anchors: AnchorConfig,
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
}

impl Default for AlignStageConfig {
    fn default() -> Self {
        Self {
            anchors: AnchorConfig::default(),
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
            max_fragment_bytes: 112 << 20,
        }
    }
}

impl AlignStageConfig {
    /// Checks the anchors, the scores and the limits.
    pub fn validate(&self) -> Result<(), AlignStageError> {
        self.anchors.validate()?;
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
        {
            return Err(AlignStageError::InvalidConfig("memory limits"));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_defaults_are_those_of_jam_0_10() {
        let config = AlignStageConfig::default();
        config.validate().unwrap();
        assert_eq!(config.anchors.k, 15);
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
    fn invalid_settings_are_rejected() {
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
