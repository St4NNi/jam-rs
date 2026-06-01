//! The three stages of one search against one shard: sketch, locate and align.
//!
//! Stage 1 (`sketch_stage`) and stage 2 (`locate_stage`) choose candidate entries and their
//! alignment windows without reading any sequence; stage 3 (`align_stage`) reads the windows,
//! aligns the query there and decides each entry. `located` is what stage 2 hands to stage 3.

pub mod align_stage;
pub mod control;
pub mod locate_stage;
pub mod located;
pub mod sketch_stage;

pub use control::{SearchControl, Stopped};
pub use located::{ContigWindows, LocatedEntry};
