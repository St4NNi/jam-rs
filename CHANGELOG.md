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
- Local Jam Index format 1 datasets with append-only independent parts. Each
  part combines a pure version-3 `.jam` screening shard, position-free
  contig signatures, a compact contig directory, and complete two-bit contig
  sequences.
- `jam index build`, `jam index append`, and `jam trace --index` for parallel
  part screening, bounded contig selection, on-demand dense k=31/k=21 seed
  generation, and one JSONL/JSONL.zst result.
- Deterministic anonymous short-trace generation and scoring across 80 to
  1,000 bases, five identity levels, both strands, origin crossing, separate
  fragments, overlap, repeat-shared, integrated, rare, and unrelated cases.

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
- Alignment refinement reuses one packed traceback matrix, and Jam Index
  candidate concurrency is admitted before work according to query length.
- Jam Index query parsing preserves IUPAC symbols so verified complete-contig
  matches can bypass dense chaining and emit direct `=` alignments.

### Compatibility

- `.jam` binary format remains version 3.
- A checked fixture created by jam-rs 0.9.11 remains readable.
- Uniform `jam dist` columns remain compatible. Bias `dist` column names changed to correct previously ambiguous semantics.
- JSON manifests are additive sidecars; legacy version-3 databases without sidecars remain readable.
- JMA is a separate format at version 1; trace JSON records use schema version `2.0.0` and do not change `.jam` v3.
- Jam Index format 1 is local-only and has no compatibility bridge to another
  Jam Index layout.

### Validation status

- Deterministic synthetic truth and local/mock-remote integration tests pass.
- The selected anonymous 40-case Jam Index run recovered 40/40 candidates at
  1.0 base precision, 0.981834375 base recall, and 1.0 interval precision and
  recall. Three paired process-cold four-thread runs had a 4.364 s wall median
  and 313,300 KiB peak-RSS median.
- Sequence-backed queries are reported under anonymous aliases. Jam Index
  operation is local-only; read-derived assemblies, independently supported
  natural positives, and the 1,000-assembly/100-query release scale remain
  unmeasured; no general production performance or accuracy claim is made.
