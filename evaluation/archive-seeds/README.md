# Archive and seed measurements

This collection describes the reproducible diagnostics and bounded experiment
entry points for the archive and trace work.  Diagnostic input is the JSON form
of `TraceDiagnosticReport` and uses zero-based, half-open query coordinates.
`tools/trace_failure_analysis/intervals.py` is the shared projection and union
implementation used by these checks.

## Stage accounting

Each stage row records seed pages, keys, decoded occurrences, retained anchors,
retained chains, sequence blocks, alignment attempts, newly supported query
bases, wall and CPU microseconds, and bytes read.  Each truth interval records
candidate selection, per-scheme seed and occurrence counts, chain spans,
sequence-window work, band retries and edge contact, acceptance, mosaic support,
and the output cap state.

Missing truth bases are assigned exactly one of:

`candidate_miss`, `no_retained_seed`, `no_matching_seed`,
`repetitive_seed_suppression`, `anchor_cap`, `no_valid_chain`, `chain_limit`,
`sequence_budget`, `alignment_band_failure`, `alignment_rejection`,
`alignment_cap`, `mosaic_selection`, or `other`.

Reports with supported coordinates are unioned before accounting.  Count-only
reports are accepted only when truth intervals are disjoint; this prevents an
overlapping observation from inflating coverage.  The analysis fails closed if
an attribution total does not equal `truth_bases - supported_bases`.

## Commands

Use a new ignored workspace output for each run:

```text
python3 tools/trace_failure_analysis/analyze.py \
  --input path/to/diagnostics.jsonl \
  --output tools/out/archive-seeds-run-<label>.json \
  --pretty

python3 tools/trace_failure_analysis/assert_baseline.py \
  --normalized path/to/normalized/jam_fragments_sensitive_s1_t8.json
```

The five manifest runners default to planning and never execute a command
unless `--execute` is supplied.  Their manifests are validated against
`schemas/experiment-manifest-v1.schema.json` and use argument arrays rather
than shell strings:

```text
run_seed_screen.py       small correctness and 40-case seed ablations
run_archive_matrix.py    native archive and optional Pithos policy matrix
run_local_matrix.py      mmap process-cold and process-warm measurements
run_remote_matrix.py     bounded mock-HTTP or object-store range measurements
run_comparators.py       matched-recall BLASTn and minimap2 measurements
normalize_comparator.py  shared CIGAR/edit-script query projection
```

For example:

```text
python3 evaluation/archive-seeds/scripts/run_seed_screen.py \
  --manifest evaluation/archive-seeds/manifests/seed-screen.example.json \
  --output tools/out/archive-seeds-plan-<label>.json
```

Execution rows retain command metadata, wall and child CPU time, return code,
stdout size and checksum, parsed JSON when available, and a bounded stderr
tail.  Long collections are intentionally explicit rather than run by these
focused checks.

Comparator output is normalized before comparison.  `normalize_comparator.py`
projects query segments, CIGARs, or Jam edit scripts through the same interval
union code and reports supported bases, aligned deletions, and interval overlap
metrics in `normalized-comparator-v1`.

## Frozen baseline

The checked-in fixture at `tests/data/archive_seeds/baseline_40_case.json`
asserts the normalized sensitive 40-case fragment and indel matrix:

```text
selected candidates  40 / 40
base precision       1.00000000
base recall          0.79113125
interval precision   1.00000000
interval recall      0.95000000
missing truth bases  66,838
truth bases          320,000
observed bases       253,162
```

This is an accuracy accounting fixture, not a claim about a new seed or archive
default.  The normalized artifact and its source checksums remain the evidence
for a particular run.

## Gear-mer history and study plan

No complete Gear-mer implementation or design note was found in the jam-rs
source or reachable repository history.  The recovered research direction is
to compare independent Gear streams (one, two, four, and eight tables or
threshold streams) under equal token-count and serialized-byte budgets.  The
screen should use family/source-held-out 500--3000 bp fragments, substitutions
and indels, boundary stability, exact sharing, rare-tail recall, and equal-size
k=21 retention.  Exact fragments, Gear-selected anchors, and fragment
mini-sketches are separate variants and must not be conflated with that
recovered direction.

## Measurement boundaries

Archive size, seed-index size, decoded bytes, range count, coalesced ranges,
request counts, cache state, and peak memory belong in the command-produced
measurement JSON.  Comparator claims are valid only after normalizing all output
onto the shared query-coordinate intervals and matching recall.  No result is
reported here for a collection that has not been executed.
