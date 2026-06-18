//! jam (just another minimizer) finds small sequences in very large sets of assemblies and aligns
//! them with CIGAR strings. `ARCHITECTURE.md` describes the stages, the terms and the formats.
//!
//! A module uses only modules listed above it in section 4 of `ARCHITECTURE.md`.

pub mod align;
pub mod cli;
pub mod decide;
pub mod index;
pub mod loss;
pub mod report;
pub mod screen;
pub mod search;
pub mod seeds;
pub mod service;
pub mod storage;

pub use cli::run;

#[cfg(test)]
mod test_support;
