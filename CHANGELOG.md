# Changelog

All notable changes to this project are documented here.

## [0.10.0] - Unreleased

### Added

- `jam screen` for per-contig catalog search and assembly-reference union summaries.
- Stable `1.0.0` contig TSV, assembly TSV, and run-metadata JSON schemas.
- Automatic database and bias-table provenance sidecars with SHA-256 identities.
- `jam stats --json` for machine-readable database inspection.
- Synthetic truth tests for exact, reverse-complement, partial, split, overlapping, ambiguous-repeat, negative, empty, malformed, mismatch, top-k, thread-determinism, and hash-zero cases.
- `jam archive` for deterministic JMA v1 positional archives and `jam trace` for one-plasmid-to-many-metagenome positional follow-up.
- Incremental trace JSONL/Zstandard output with alignment, circular coverage, gap, failure, and resource-metric records.
- Local, `file://`, HTTP(S), and S3-compatible JMA/raw resources, plus identity-checked remote `.jam` materialization into an explicit cache directory.
- Trace component and end-to-end local/mock-remote benchmark harnesses with optional BLASTn and minimap2 comparison stages.

### Changed

- Pinned `jamhash = 0.1.2` and named its released algorithm identity `jamhash_u64_v1`.
- Excluded hash zero consistently during database and query sketch construction.
- Renamed bias evidence to retained/weighted containment, added bias-table IDs and score mode, and made the bias-mode uniform E-value `NA`.
- Added a default 25% positive-retention guard to bias-table calibration.
- Bounded sketch writer file descriptors by opening bucket files only while flushing.
- Corrected indexed trace alignment to preserve chain-derived circular query coordinates for split fragments.

### Compatibility

- `.jam` binary format remains version 3.
- A checked fixture created by jam-rs 0.9.11 remains readable.
- Uniform `jam dist` columns remain compatible. Bias `dist` column names changed to correct previously ambiguous semantics.
- JSON manifests are additive sidecars; legacy version-3 databases without sidecars remain readable.
- JMA is a separate format at version 1; trace JSON records use schema version `1.0.0` and do not change `.jam` v3.

### Validation status

- Deterministic synthetic truth and local/mock-remote integration tests pass.
- Production-scale real-background, actual-S3, BLASTn, and minimap2 validation is pending and no production performance or accuracy claim is made.
