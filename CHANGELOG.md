# Changelog

All notable changes to this project are documented here.

## [0.10.0] - Unreleased

### Added

- `jam screen` for per-contig catalog search and assembly-reference union summaries.
- Stable `1.0.0` contig TSV, assembly TSV, and run-metadata JSON schemas.
- Automatic database and bias-table provenance sidecars with SHA-256 identities.
- `jam stats --json` for machine-readable database inspection.
- Synthetic truth tests for exact, reverse-complement, partial, split, overlapping, ambiguous-repeat, negative, empty, malformed, mismatch, top-k, thread-determinism, and hash-zero cases.
- `jam archive` for deterministic JMA v1 positional archives and `jam trace` for one-query-to-many-metagenome positional follow-up.
- Incremental trace JSONL/Zstandard output with alignment, circular coverage, gap, failure, and resource-metric records.
- Local, `file://`, HTTP(S), and S3-compatible JMA/raw resources, plus identity-checked remote `.jam` materialization into an explicit cache directory.
- Trace component and end-to-end local/mock-remote benchmark harnesses with optional BLASTn and minimap2 comparison stages.
- Stable `jam-seed-chain-align-v1` algorithm metadata and resolved parameters in trace JSON records and JMA provenance.
- Checksum-bound JMA v1 range indexes for partial seed-bucket and sequence-block reads, with explicit eager fallback.
- `--memory-target` as the canonical name for the internal memory budget; `--memory` remains a compatibility alias.
- Generic `jam trace --query` input with plasmid, phage, other, and unknown
  query kinds; linear, circular, auto, and unknown coordinate handling;
  `--plasmid` remains a deprecated one-release compatibility alias.
- Deterministic query-coordinate fragment mosaics, primary/secondary alignment
  roles, gap-directed k=31/k=21 rescue, common/repeat evidence, and schema 2.0
  JSONL records retaining all accepted alignments.
- Plasmid-plus-phage controlled evaluation inputs and pinned LexicMap v0.9.0
  comparison support alongside BLASTn and minimap2.

### Changed

- Pinned `jamhash = 0.1.2` and named its released algorithm identity `jamhash_u64_v1`.
- Excluded hash zero consistently during database and query sketch construction.
- Renamed bias evidence to retained/weighted containment, added bias-table IDs and score mode, and made the bias-mode uniform E-value `NA`.
- Added a default 25% positive-retention guard to bias-table calibration.
- Bounded sketch writer file descriptors by opening bucket files only while flushing.
- Corrected indexed trace alignment to preserve chain-derived circular query coordinates for split fragments.
- Seed-level JMA incompatibility now fails closed instead of producing
  unanchored raw-fallback alignments.
- Structural/index JMA failures no longer trigger unanchored raw fallback;
  rescue rounds use the resolved seed scale and preserve earlier evidence.

### Compatibility

- `.jam` binary format remains version 3.
- A checked fixture created by jam-rs 0.9.11 remains readable.
- Uniform `jam dist` columns remain compatible. Bias `dist` column names changed to correct previously ambiguous semantics.
- JSON manifests are additive sidecars; legacy version-3 databases without sidecars remain readable.
- JMA is a separate format at version 1; trace JSON records use schema version `2.0.0` and do not change `.jam` v3.

### Validation status

- Deterministic synthetic truth and local/mock-remote integration tests pass.
- Accession-backed controlled spike-ins, local range/memory measurements, and
  pinned BLASTn/minimap2/LexicMap comparisons are recorded with bounded claims.
- Actual S3, read-derived assemblies, independently supported natural
  positives, and the 1,000-assembly/100-query release scale remain pending; no
  general production performance or accuracy claim is made.
