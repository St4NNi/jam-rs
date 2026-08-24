//! Experimental Gear-mer positional seed paths.
//!
//! The historical Gear proposal is recovered from commit `c671736` and its
//! `tests/cnn_plan.md`: use independent fine (`64/192/384`) and coarse
//! (`256/768/1536`) streams, evaluate both forward and reverse-complement
//! sequence views, and keep the Gear integer out of biological evidence.  The
//! positional implementations in this module are new search variants built
//! around that recovered chunking rule.  A Gear digest is always a lookup
//! candidate; callers must perform exact fragment or packed-seed verification
//! before retaining an anchor.

mod anchors;
mod fragments;
mod mini_sketch;
mod tables;

pub use anchors::{
    GearAnchor, GearAnchorConfig, GearAnchorIndex, GearAnchorScheme, GearAnchorVerification,
    GearSelectedAnchorError, SPACED31_WEIGHT21_MASK, canonical_anchor_at, gear_selected_anchors,
    verify_anchor_candidate, verify_anchor_sequences,
};
pub use fragments::{
    ExactFragment, ExactFragmentIndex, ExactRun, FragmentLengthStats, FragmentMatch,
    FragmentOrientation, FragmentationMode, GearConfig, GearError, GearStream,
    StrandSymmetricDistribution, VerifiedFragment, boundary_distribution, digest_fragment,
    fragment_boundaries, fragment_bytes, fragment_sequence, merge_exact_runs,
    strand_symmetric_distribution, verify_exact_fragment,
};
pub use mini_sketch::{
    FragmentMiniSketch, FragmentMiniSketchIndex, MiniSketchConfig, MiniSketchError,
    MiniSketchLocalPair, MiniSketchSeed, RelatedFragmentCandidate, build_fragment_mini_sketch,
    build_fragment_mini_sketches, lookup_related_fragments,
};
pub use tables::{GearTable, GearTableKind, INDEPENDENT_TABLE_SEEDS};
