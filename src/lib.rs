//! jam (just another minimizer) finds small sequences in very large sets of assemblies and aligns
//! them with CIGAR strings; `ARCHITECTURE.md` describes the stages, terms and module order.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

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
