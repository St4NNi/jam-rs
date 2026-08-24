# `jam trace` JSONL schema

`jam trace` emits schema `2.0.0` as one newline-delimited JSON record stream:

```text
run_header
metagenome_result*
run_footer
```

The default output is one UTF-8 `result.jsonl` file. A `.zst` or `.zstd`
suffix selects a Zstandard stream containing the same records. Output is
written incrementally and published atomically only after the footer has been
written and the JSONL or Zstandard stream has been finalized. A flushed prefix
without a footer is incomplete and must not be interpreted as a negative
result.

Validate records against the repository schema at
[`schemas/jam-trace-v2.schema.json`](../schemas/jam-trace-v2.schema.json). The
v1 schema remains available for older `1.x` streams; it is not the current
`2.0.0` record contract.

## Query and topology vocabulary

The stream is centered on one query sequence element searched against many
metagenome resources. The query may be a plasmid, phage, other mobile genetic
element, synthetic construct, or an unknown nucleotide reference.

`query_kind` is one of `plasmid`, `phage`, `other`, or `unknown`.
`topology_requested` is one of `linear`, `circular`, `auto`, or `unknown`.
`coordinate_model` describes coordinate handling, not molecule biology:

| Value | Meaning |
| --- | --- |
| `linear` | Ordinary `[0, query_length)` coordinates and separate terminal gaps. |
| `wrap` | Coordinates may cross the stored origin and are normalized on a circle. |
| `undetermined` | Linear coordinates are the primary display while alternative wrap evidence is retained when evaluated. |

`topology_evidence` is one of `linear_supported`, `wrap_supported`,
`both_compatible`, `insufficient`, or `undetermined`. A `wrap` coordinate
model is not proof that the biological element is circular. It can also
represent a rotated or circularly permuted reference. Separate contigs remain
physically independent, and sequence support does not prove an autonomous
element.

## Algorithm identity

The header and each metagenome result include `algorithms` with these stable
workflow-layer identifiers:

```json
{
  "screen_algorithm": "jam-fracminhash-screen-v1",
  "local_alignment_algorithm": "jam-exact-seed-chain-banded-v1",
  "mosaic_algorithm": "jam-fragment-mosaic-v1",
  "trace_workflow": "jam-trace-v1"
}
```

The `algorithm` object records the resolved implementation parameters,
including `hash_id`, seed selection, sensitivity profile, k-mer layers,
chain limits, alignment mode/scoring, and coverage mode. The screening path
uses `jamhash_u64_v1`; k=31 is the primary layer and k=21 is conditional rescue.
No seed size below 21 is represented by this contract. The identifiers describe
the implemented code and are not a claim of biological classification.

## Run header

The first record has `record_type="run_header"` and includes:

| Field | Meaning |
| --- | --- |
| `schema_version` | `2.0.0`. |
| `run_id` | Identifier shared by all records in the stream. |
| `jam_rs_version`, `source_commit` | Software version and source identity, when available. |
| `started_at_utc` | Run start timestamp. |
| `command` | Redacted argument vector used to start the run. Credentials and signed URL query strings are not emitted. |
| `query_id`, `query_length` | Stable query identifier and length in bases. Legacy `plasmid_id`/`plasmid_length` aliases may be accepted by compatibility deserializers. |
| `query_kind`, `topology_requested` | Query class and requested coordinate/topology handling. |
| `threads`, `io_concurrency` | Resolved CPU and resource-I/O concurrency. |
| `sensitivity` | Fully resolved fast, balanced, or sensitive seed, rescue, candidate, chain, window, and memory-related limits. |
| `algorithms`, `algorithm` | Layered identifiers plus resolved parameters. |
| `inputs` | Redacted query, database, catalog, and metagenome locators with optional checksums. |

## Metagenome result

Each `metagenome_result` is independently consumable and contains:

| Field | Meaning |
| --- | --- |
| `query_id`, `metagenome_id` | Query and candidate resource identifiers. |
| `query_kind`, `topology_requested` | Generic query/topology context repeated for the result. |
| `coordinate_model`, `topology_evidence` | Selected coordinate representation and evidence interpretation. |
| `algorithms`, `algorithm` | Resolved algorithm identity and parameters. |
| `status` | `complete`, `partial`, `no_candidate`, or `failed`. |
| `candidate` | JAM sketch-screening evidence, or `null`. |
| `alignments` | Every accepted local alignment, including overlapping and alternative mappings. |
| `primary_fragment_mosaic` | Query-coordinate nonredundant mosaic, or `null` when no alignment was accepted. |
| `topology` | Linear and optional wrap model summaries used for topology-coordinate assessment. |
| `rescue_rounds` | Per-round gap-directed seed and alignment counters. |
| `performance_counters` | Candidate, contig, window, alignment, failure, and elapsed-time counters. |
| `coverage` | Legacy compact coverage summary when produced; the fragment mosaic is the detailed schema-2 evidence. |
| `warnings` | Nonfatal interpretation or fallback warnings. |
| `failures` | Structured stage/code/message/resource/retryability records. |
| `resource_metrics` | Resource request, byte, cache, retry, fallback, and indexed-read counters. |

Candidate containment fields retain explicit denominators:

```text
query_containment       = shared_hashes / query_hashes
metagenome_containment  = shared_hashes / metagenome_hashes
```

Legacy serialized names `plasmid_containment` and `plasmid_hashes` may be
accepted as aliases by compatibility deserializers. In bias mode,
`uniform_hash_e_value` is `null` unless an explicit approximation is requested;
the approximation is not a calibrated significance value and must never be
used to make a final presence call. Bias-weighted values remain separate from
raw genomic containment.

## Accepted alignments

An alignment record includes an `alignment_id`, `contig_id`, `strand`, query
segments, forward contig target coordinates, score and identity, edit counts,
CIGAR, lossless edit script, chain score, seed evidence, and mosaic role.
Coordinates are zero-based half-open. A reverse strand reports the interval in
forward contig coordinates; it does not reverse the coordinate convention.

`query_segments` has one ordinary interval for a linear alignment and may have
two ordered non-wrapping intervals for an origin-crossing alignment. A deletion
consumes query coordinates but is not base support. An insertion consumes only
target coordinates. CIGAR and edit-script counts must reconstruct the same
query/target consumption.

`seed_evidence` records primary and rescue anchor counts, nonrepetitive and
common anchor counts, repetitive seed counts, and occurrence/document-frequency
limits. A hash lookup is not itself an alignment anchor; exact packed-k-mer
verification occurs before an anchor is accepted.

Allowed `role` values are:

```text
primary_mosaic
overlapping_support
alternative_mapping
common_sequence
repeat_only
origin_crossing
```

Common or repeat-only labels annotate evidence; they do not discard a valid
alignment and do not establish element origin. `origin_crossing` describes
coordinate evidence only.

## Fragment mosaic and atomic intervals

`primary_fragment_mosaic` is the query-centered evidence summary. It splits
the query at every accepted alignment boundary into `atomic_intervals`. Each
atomic interval has one deterministic `primary_alignment_id` when base support
exists and a list of retained `alternative_alignment_ids`.

The primary selection order uses alignment quality, identity, score per query
base, nonrepetitive anchor evidence, supported length, contig identifier, and
stable alignment identifier. The order is deterministic across thread counts
and input ordering after canonical alignment IDs are assigned.

The mosaic reports:

| Field | Meaning |
| --- | --- |
| `base_covered_bases`, `base_coverage_fraction` | Query bases paired to a metagenome base by `=` or `X`, counted at most once. |
| `aligned_span_bases`, `aligned_span_fraction` | Query bases inside accepted alignment spans, including `D` operations. |
| `covered_intervals` | Nonredundant supported query union. |
| `unsupported_gaps` | Query intervals without accepted base support. |
| `alignment_deletions` | Query intervals crossed by an accepted alignment `D` operation. These are distinct from unsupported gaps. |
| `supporting_contigs` | Unique contigs contributing accepted evidence. |
| `accepted_alignment_count` | Coordinate-compatible accepted alignments represented by the mosaic. |
| `alternative_alignment_count` | Alignments retained as overlapping or alternative support. |
| `common_sequence_intervals`, `common_sequence_supported_bases` | Support involving common-sequence anchors. |
| `repeat_only_intervals`, `repeat_only_supported_bases` | Support with repetitive evidence and no nonrepetitive anchors. |
| `nonrepetitive_supported_bases` | Union of support carrying nonrepetitive anchor evidence. |
| `alignment_evidence` | Per-alignment primary/secondary/newly-supported bases and role. |
| `selection_components` | Local score sum, newly supported bases, nonrepetitive evidence, redundant overlap, contradiction counts, and fragment count. |

Overlap, duplicate contigs, and terminal repeats cannot increase a query base
twice. Separate supporting contigs are not joined into one physical molecule.

## Gaps and coordinate models

Every gap has:

```json
{
  "interval": {"start": 0, "end": 10},
  "segments": [{"start": 0, "end": 10}],
  "wraps_origin": false,
  "length": 10
}
```

For a linear model, terminal gaps remain separate one-segment records. For a
wrap model, the first segment beginning at zero and the final segment ending at
the query length are merged into one record with two `segments`,
`wraps_origin=true`, and the combined length. `interval` remains the first
ordinary segment for compatibility. A wrapped gap record is a coordinate
representation, not evidence of circular biological topology.

`alignment_deletions` are never folded into `unsupported_gaps`: a deletion is
inside an accepted alignment, while an unsupported gap has no accepted base
support. Circular complements are origin-invariant; linear complements retain
leading and trailing gaps separately.

## Topology model assessment

`topology` contains a `linear_model` and, when evaluated, a `wrap_model`. Auto
assessment reuses the accepted alignment slice and compares newly supported
query bases using the configured margin, with alignment
quality, nonrepetitive anchors, origin-crossing count, terminal gaps, and
fragment count retained in the model evidence. Those additional components do
not currently select the model. A model is selected only when the supported-
base difference exceeds the margin. Ties or insufficient evidence
remain `coordinate_model="undetermined"`.

`topology_requested="unknown"` keeps linear coordinates as the primary
display, exposes a wrap alternative when evaluated, and emits
`topology_evidence="undetermined"`. No topology setting turns sequence
evidence into proof of a biological circle, phage genome structure, or
autonomous mobile element.

## Gap-directed rescue and counters

`rescue_rounds` records the monotonic gap-directed extension of the initial
search. Each round contains:

```text
round, seed_k, seed_scale, target_gaps,
seed_buckets_requested, seed_keys_tested, anchors_created,
chains_accepted, sequence_blocks_fetched, alignment_windows_attempted,
new_query_bases_supported, elapsed_millis
```

Later rounds may add support but do not remove accepted evidence. Balanced and
sensitive profiles may use dense k=31 seeds in unresolved gaps and conditional
k=21 rescue; fast mode may stop after the initial round.

`performance_counters` records candidate/contig/window/alignment counts,
failures, and elapsed milliseconds. `resource_metrics` records:

```text
metadata_requests, head_requests, get_requests, range_requests,
stream_requests, requested_bytes, returned_bytes, decoded_bytes,
remote_bytes, cache_bytes, cache_hits, cache_misses, cache_evictions,
stale_cache_rejections, retries, full_object_fallbacks,
seed_buckets_read, sequence_blocks_read
```

Indexed local or remote JMA reads should be selective. A full-object fallback
is explicit in the counters. Noncandidate resources must not be opened.

## Structured failures and completion

Each failure has `stage`, stable `code`, human-readable `message`, optional
redacted `resource`, and `retryable`. A failed candidate does not hide
successful results for other metagenomes. Resource errors must not expose
credentials, signed query strings, or secret headers.

The footer contains stream totals, completion time, and aggregate resource
metrics. Consumers should require exactly one header first and one footer last
before treating the run as complete.

## Jam Index records

Jam Index uses the same schema-2.0 JSONL record sequence. Its run header has an
input with role `jam_index` bound to the root `manifest.json`; it does not add
JMA or remote-resource inputs. Candidate evidence uses
`candidate.score_mode="jam_index"`. The candidate's shared-hash counts and
containment values remain routing evidence, not alignment evidence.

Ordinary dense chaining stages retain their existing names. A verified
complete selected contig can instead emit stage 0 with
`name="exact_contig"`, zero DP alignment attempts, and direct `=` alignment
evidence. `archive_metrics.seed_bytes_read` is zero for Jam Index because no
persistent positional seed section is read. `sequence_bytes_read` and
`decoded_sequence_bases` account selected complete contigs.

The output remains one JSONL or JSONL.zst stream. Separate contig IDs remain
attached to their alignments and mosaic evidence; their appearance in one
query-coordinate mosaic does not assert physical linkage.

Evidence labels such as `alert` and `supported` are suitable downstream
reports. `confirmed` requires independent mapping, marker, graph, or
contextual evidence; `jam trace` does not emit a confirmed biological presence
call.

## Compatibility boundaries

Schema `2.0.0` is the generic query-centered contract. A `1.x` JSONL stream is
not silently upgraded: use a v1 schema/adapter or regenerate the result. The
query/plasmid field aliases retained in Rust deserializers are compatibility
bridges, not permission to mix v1 and v2 records in one stream. `.jam` version 3
database compatibility and `jamhash_u64_v1` remain separate from JSON schema
compatibility. The CLI `--plasmid` alias is likewise a migration convenience;
new commands should use `--query`.
