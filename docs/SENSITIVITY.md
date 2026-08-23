# Trace sensitivity profiles

Trace sensitivity is an execution profile for searching one query element
across metagenomes with explicit resource bounds. It is not a calibrated
probability of plasmid, phage, or mobile-element recovery and it should not be
reported as a biological sensitivity estimate. The resolved profile is
recorded in the JSON `run_header` so a result can be reproduced.

All profiles use `jamhash_u64_v1`, FracMinHash retention, and fixed seed
identities:

- primary seeds: `k=31`;
- optional rescue seeds: `k=21`;
- hash zero is skipped at both archive-build and query time;
- exact packed canonical k-mers are checked after hash lookup.

| Profile | Primary | Rescue | Max candidates | Max anchors | Max chains | Minimum anchors | Alignment band | Window bases | Concurrent candidates |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `fast` | k=31, scale 500 | none | 25 | 25,000 | 4 | 3 | 64 | 250,000 | 4 |
| `balanced` | k=31, scale 200 | k=21, scale 500 | 100 | 100,000 | 8 | 3 | 128 | 1,000,000 | 4 |
| `sensitive` | k=31, scale 100 | k=21, scale 200 | 250 | 250,000 | 16 | 2 | 256 | 4,000,000 | 2 |

The resolved profiles also set these gap-rescue values:

| Profile | Rounds | Dense k=31 | Minimum gap | Flank | Max rescue buckets | Max sequence blocks | Max alignment windows |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| `fast` | 1 | none | 250 | 64 | 0 | 0 | 0 |
| `balanced` | 3 | scale 100, max occurrences 96 | 200 | 96 | 50,000 | 1,024 | 256 |
| `sensitive` | 3 | scale 100, max occurrences 192 | 100 | 128 | 100,000 | 2,048 | 512 |

The profile limits are resolved into the run header, including seed scales,
candidate and anchor caps, chain limits, alignment scoring/band width, gap
rescue rounds, and concurrent candidate limits. They are operational bounds,
not hidden classifier thresholds.

The primary scale is a FracMinHash scale: smaller values retain a denser
subset. The rescue level is attempted only when primary positional evidence is
insufficient or an unresolved query gap qualifies for gap-directed rescue. A
denser archive can serve a sparser query level, subject to the archive's
available k and scale sections. A missing or incompatible seed level is an
explicit error, not a zero-evidence result.

A seed-level mismatch is not eligible for raw fallback. This prevents an
incompatible indexed request from becoming unanchored local-alignment evidence.

## Choosing a profile

Use `fast` for exploratory work where lower I/O and fewer candidates matter.
Use `balanced` as the documented starting point for a thesis or surveillance
run. Use `sensitive` when short, fragmented, or more divergent represented
traces justify extra work and memory. Benchmark the chosen profile on the
versioned catalog and assembly set; profile names alone do not establish
recall, specificity, or a universal speed ranking.

The profile's `band_width` controls the alignment kernel's local search band,
not a biological divergence threshold. The chain and window limits bound work;
they can truncate evidence for unusually repetitive or long candidates. A
truncated or failed candidate is represented in the JSON `failures` field and
must not be silently treated as a clean negative.

## Gap-directed rescue

Rescue is query-centered and monotonic. The runner first searches the complete
query with the primary k=31 level, builds local alignments, and constructs a
nonredundant query-coordinate mosaic. It then identifies unsupported query
gaps that meet the profile's minimum length and requests denser seeds only in
those gaps plus the configured flanks. The next eligible round can use dense
k=31; when the profile supplies it, the final rescue round can use k=21.

Each round is bounded by its seed-bucket, sequence-block, and alignment-window
limits. It stops when no new query bases are supported, no eligible buckets
remain, the round limit is reached, or the candidate resource budget is
exhausted. Later rounds add evidence and do not remove earlier accepted
alignments. `rescue_rounds` records the requested gaps, seed keys, anchors,
chains, fetched blocks, attempted windows, newly supported bases, and elapsed
time.

`fast` has only the initial sparse k=31 pass. `balanced` adds dense k=31 and
conditional k=21 work for unresolved gaps. `sensitive` uses the same stages
with denser seed retention, lower minimum chain size, wider alignment bands,
and larger rescue budgets. k=21 remains conditional rescue evidence, not a
replacement for the k=31 baseline.

## Archive settings

Build JMA v1 with seed levels that cover the profiles you intend to run. For
`balanced`:

```bash
jam archive \
  --input assembly.fasta \
  --output assembly.jma \
  --primary-scale 200 \
  --rescue-scale 500
```

Use `--primary-scale 100 --rescue-scale 200` when one archive must also support
`sensitive`; those denser nested levels can serve the sparser profiles.

The archive command always creates the fixed k=31 primary section and creates
k=21 rescue data unless `--no-rescue` is supplied. An archive built without
rescue data cannot satisfy a profile that requests k=21 rescue. An archive
with only a sparse k=31 level cannot satisfy a profile requiring a denser
level. The command reports archive counts; retain those counts and the input
checksum with the archive for provenance.

`--complexity` applies the optional entropy filter during archive seed
construction. Candidate `.jam` k-mer, scale, entropy, and optional bias
settings belong to the separate sketch retrieval stage and must match between
that database and its query; they do not have to equal the JMA positional seed
settings. A bias table is part of the `.jam` identity, not a sensitivity
profile, and is tied to its training inputs and sampling configuration.

## k-mer rationale and evidence boundary

k=31 primary seeds provide longer exact anchors and therefore reduce common
repeat collisions in positional chaining. k=21 rescue seeds retain more
divergent or short represented query sequence at the cost of more repetitive
occurrences. The two levels serve different retrieval roles; neither makes a
sketch or alignment a final element-presence call or a biological topology
call. This release preserves the fixed k=31/k=21 contract and does not restart
classifier training or add new collision analysis below k=21.

## Query class and coordinate handling

Sensitivity controls how much evidence is searched; it does not infer the
query class. `--query-kind` may be `plasmid`, `phage`, `other`, or `unknown`.
The same seed and rescue profiles apply to each class.

`--topology linear` uses ordinary query coordinates and never permits a chain
to cross the stored origin. `--topology circular` permits a wrap coordinate
model. `--topology auto` evaluates linear and wrap mosaics from the same local
alignment set and applies `--topology-margin-bases` (the profile default is
250 for `fast`, 200 for `balanced`, and 100 for `sensitive`) when selecting a
coordinate model. `--topology unknown` reports an undetermined coordinate and
topology evidence; it does not turn a wrapped coordinate representation into a
biological topology call.

Changing a profile or topology setting can change which candidate resources,
seed buckets, sequence blocks, and alignment windows are examined. Preserve
the resolved run header and the catalog/database identities when comparing
runs. A candidate or supported fragment remains evidence requiring independent
confirmation, whether the query is a plasmid, phage, or another sequence
element.

## Resource controls

Global `--threads` bounds parallel candidate processing and trace
`--io-concurrency` bounds candidate resource tasks. Global `--memory-target`
sets the resource cache budget used by `jam trace`. Retained anchors, chains,
alignment windows, alignments, and output records have profile/runner limits.
The value is an internal target, not a hard process RSS limit: raw candidates
can retain a complete parsed assembly, and temporary candidate-match, seed-
bucket, or anchor collections can peak before their retained-output limits.
Mmap pages, allocator arenas, thread stacks, parser buffers, and the runtime
also add overhead.

For a reproducible run, preserve the complete `run_header`, including the
profile fields, command, database checksum, catalog checksum, source commit,
and input checksums. Candidate ordering, chain tie-breaking, mosaic assignment,
alignment IDs, and output record ordering are deterministic across the
configured thread counts. See [TRACE_JSON.md](TRACE_JSON.md).
