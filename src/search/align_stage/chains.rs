//! Aligning the chains of one window. A chain whose anchors are already on an alignment of the
//! window is skipped. Anchors that an alignment leaves out, such as a second conserved block
//! beyond a costly insertion, are aligned again on their own.

use super::contig::{ContigJob, WindowOutput};
use super::paths::AlignedPath;
use super::window::PreparedWindow;
use super::{AlignStageConfig, AlignStageError};
use crate::align::{AlignmentWorkspace, Chain, align_chain, split_by_drift};

/// Query bases that the alignment of leftover anchors may reach beyond them on each side. Two
/// blocks that one alignment would have joined are less than this far apart (bridging them costs
/// less than the smaller block earns), so the window of the second block holds no part of the
/// first.
const LEFTOVER_REACH: u64 = 64;

/// Aligns every chain of `window`, with the leftover anchors of each alignment.
pub(super) fn align_chains(
    job: &ContigJob<'_>,
    window: &PreparedWindow,
    output: &mut WindowOutput<'_>,
    workspace: &mut AlignmentWorkspace,
) -> Result<(), AlignStageError> {
    let config = job.config;
    let length = job.query.length();
    for chain in &window.chains {
        let parts = split_by_drift(chain, config.band.max_drift, config.chains.min_anchors);
        // The parts in chain order, each allowed to reach the whole query.
        let mut pending: Vec<(Chain, u64)> =
            parts.into_iter().rev().map(|part| (part, length)).collect();
        while let Some((part, reach)) = pending.pop() {
            job.control.check()?;
            if output.paths.redundant(&part, length) {
                output.work.redundant_chains += 1;
                continue;
            }
            if chain_too_sparse(config, &part) {
                continue;
            }
            if let Some(rest) = align_part(job, window, output, workspace, &part, reach)? {
                pending.push((rest, LEFTOVER_REACH));
            }
        }
    }
    Ok(())
}

/// Whether the anchors of a chain are too few for its stretch to have the chain identity of
/// `config`. Anchors are the k-mers at every `stride`-th window position that the query has too,
/// and at identity `p` a k-mer is unchanged with chance `p^k`. Short chains are never judged:
/// too few anchors tell nothing.
fn chain_too_sparse(config: &AlignStageConfig, chain: &Chain) -> bool {
    const MIN_ANCHORS: usize = 24;
    if config.min_chain_identity <= 0.0 || chain.anchors.len() < MIN_ANCHORS {
        return false;
    }
    let k = u64::from(config.anchors.a);
    let stride = u64::from(config.anchors.stride);
    let first = u64::from(chain.anchors[0].0);
    let last = u64::from(chain.anchors[chain.anchors.len() - 1].0);
    let sampled = ((last + k).saturating_sub(first) / stride).max(1);
    let density = (chain.anchors.len() as f64 / sampled as f64).min(1.0);
    density.powf(1.0 / k as f64) < config.min_chain_identity
}

/// Aligns one chain part. When its alignment is kept, returns the anchors it leaves out.
fn align_part(
    job: &ContigJob<'_>,
    window: &PreparedWindow,
    output: &mut WindowOutput<'_>,
    workspace: &mut AlignmentWorkspace,
    part: &Chain,
    reach: u64,
) -> Result<Option<Chain>, AlignStageError> {
    let config = job.config;
    output.work.chains += 1;
    let aligned = align_chain(
        part,
        &window.bases,
        window.start,
        &job.query.table,
        config.band,
        workspace,
        reach,
    )?;
    let Some(alignment) = aligned else {
        return Ok(None);
    };
    output.work.alignments += 1;
    let path = AlignedPath::of(&alignment, window.bases.len() as u64, window.start);
    if !output.keep_alignment(job, window, alignment)? {
        return Ok(None);
    }
    output.paths.add(&path);
    Ok(output
        .paths
        .uncovered(part, job.query.length(), config.chains.min_anchors))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::Strand;

    fn chain(anchors: usize, spacing: u32) -> Chain {
        Chain {
            strand: Strand::Forward,
            anchors: (0..anchors as u32)
                .map(|i| (i * spacing, i * spacing))
                .collect(),
            offsets: (0..anchors as u32).map(|i| i * spacing).collect(),
            score: 100,
        }
    }

    #[test]
    fn sparse_chains_are_skipped_only_when_long_and_asked_for() {
        let asked = AlignStageConfig {
            min_chain_identity: 0.9,
            ..AlignStageConfig::default()
        };
        // An anchor at every looked-up position: identity 1.
        assert!(!chain_too_sparse(&asked, &chain(30, 4)));
        // One anchor in 100 looked-up positions: about 0.74.
        assert!(chain_too_sparse(&asked, &chain(30, 400)));
        assert!(!chain_too_sparse(&asked, &chain(20, 400)));
        assert!(!chain_too_sparse(
            &AlignStageConfig::default(),
            &chain(30, 400)
        ));
    }
}
