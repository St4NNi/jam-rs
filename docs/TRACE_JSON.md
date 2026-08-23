# Trace JSONL schema and consumer guide

`jam trace` writes schema version `1.2.0` as newline-delimited JSON. The
logical stream is:

```text
run_header
metagenome_result*
run_footer
```

Plain output is UTF-8 JSONL. A path ending in `.zst` or `.zstd` contains the
same records in a Zstandard frame. The writer flushes each record and does not
retain a whole run's formatted output in memory. Validate field names and
types with [jam-trace-v1.schema.json](../schemas/jam-trace-v1.schema.json).

## Run header

The first record has `record_type="run_header"` and contains:

| Field | Meaning |
| --- | --- |
| `schema_version` | Trace record schema, currently `1.2.0`. |
| `run_id` | Identifier shared by all records in the stream. |
| `jam_rs_version` | Binary version. |
| `source_commit` | Source commit when available, otherwise `null`. |
| `started_at_utc` | Run start value recorded by the command. |
| `command` | Argument vector used to start the run. |
| `plasmid_id` | Query plasmid identifier. |
| `plasmid_length` | Query length in bases. |
| `sensitivity` | Profile, fixed k=31 primary/k=21 rescue settings, scoring, and resource bounds. |
| `algorithm` | Stable algorithm identifier, version, and all resolved seed, chain, alignment, and coverage parameters. |
| `inputs` | Redacted database, catalog, and plasmid locators with optional SHA-256 values. |

Locators in `inputs` are redacted. Do not infer a credential, signed URL, or
resource version from a redacted locator. Preserve the original access
configuration separately under pipeline controls.

## Metagenome result

Each `metagenome_result` is independently consumable:

| Field | Meaning |
| --- | --- |
| `status` | `complete`, `partial`, `no_candidate`, or `failed`. |
| `plasmid_id`, `metagenome_id` | Stable query and candidate identifiers. |
| `algorithm` | The same resolved algorithm identity and parameters as the run header. |
| `candidate` | Existing `.jam` sketch evidence, or `null`. |
| `alignments` | Positional alignment records, possibly empty. |
| `coverage` | Nonredundant query coverage and circular gaps, or `null`. |
| `failures` | Structured stage/code/message records; resources are redacted. |
| `resource_metrics` | HEAD/GET/range/stream, requested/returned/decoded byte, cache, retry, fallback, seed-bucket, and sequence-block counters. |

Candidate fields use explicit denominators:

```text
plasmid_containment    = shared_hashes / plasmid_hashes
metagenome_containment = shared_hashes / metagenome_hashes
```

`shared_hashes`, `plasmid_hashes`, and `metagenome_hashes` are counts of
retained distinct hashes. `rank` is one-based after configured filters and
deterministic sorting. `score_mode` is `uniform` or `bias`.

In `uniform` mode, `uniform_hash_e_value` is the value produced by the
uniform-hash candidate stage. In `bias` mode, it is `null`; the retained and
catalog-specific `bias_weighted_plasmid_containment` value is separate. A
bias-mode value is not raw genomic containment, a universal plasmid
probability, or calibrated significance. Consumers must never filter a final
result by an invented or approximated E-value in bias mode.

## Alignment fields and coordinates

An alignment contains the following evidence components:

| Field | Meaning |
| --- | --- |
| `contig_id` | Supporting assembly contig. |
| `strand` | Orientation relative to the forward stored contig. |
| `query_segments` | One query interval, or two ordered intervals for an origin crossing. |
| `target_interval` | Forward contig interval, always zero-based half-open. |
| `query_length`, `target_length` | Consumed query and target lengths. |
| `origin_crossing` | Whether the circular query crosses coordinate zero. |
| `score`, `chain_score` | Alignment and retained-chain scores. |
| `matches`, `substitutions`, `insertions`, `deletions` | Aggregate edit evidence. |
| `cigar` | Compact edit operations. |
| `edit_script` | Structured edit runs; authoritative when checking CIGAR and counters. |
| `primary` | Whether this alignment contributes to nonredundant primary coverage. |

All base coordinates are unsigned, zero-based, half-open `[start,end)`. A
reverse alignment still reports the target interval in forward stored-contig
coordinates; `strand="reverse"` says that the reverse complement of that
interval aligns to the plasmid query. A circular origin crossing has exactly
two non-wrapping query segments. Insertions consume target sequence but do not
add query coverage; deletions leave a query gap.

## Coverage fields

`coverage.plasmid_length` is the circular query length. `primary_intervals` are
the deterministic nonredundant union contributors. `secondary_intervals` show
overlapping evidence that was retained for review but does not increase
`supported_bases`. `supported_fraction` is:

```text
supported_bases / plasmid_length
```

`gaps` are unsupported circular query intervals, and `largest_gap` is the
largest gap length. Union is performed after origin segments are normalized to
`[0, plasmid_length)`, so overlapping contigs and repeated alignments cannot
double-count a base. Coverage does not establish that fragments are linked on
one molecule.

## Stream completion and examples

The footer contains bounded totals:

```json
{
  "record_type": "run_footer",
  "schema_version": "1.2.0",
  "run_id": "trace-...",
  "completed_at_utc": "...",
  "metagenomes_total": 2,
  "metagenomes_with_candidates": 1,
  "metagenomes_aligned": 1,
  "metagenomes_failed": 0,
  "alignments_total": 2,
  "resource_metrics": {}
}
```

An interrupted process can leave a flushed prefix without a footer. Consumers
must mark that stream incomplete instead of treating the prefix as a finished
negative result. A complete stream can be inspected without loading it all:

```bash
zstdcat p_001.trace.jsonl.zst | \
  jq -r 'select(.record_type == "metagenome_result") |
         [.metagenome_id,
          (if .candidate == null then "no_candidate"
           elif (.alignments | length) == 0 then "alert"
           else "supported" end),
          ([.alignments[].contig_id] | unique | join(",")),
          (.coverage.supported_fraction // 0)] | @tsv'
```

The labels above remain evidence labels. `confirmed` must be added by an
independent downstream mapping, marker, or contextual review step. jam-rs
does not emit a confirmed plasmid-presence call.
