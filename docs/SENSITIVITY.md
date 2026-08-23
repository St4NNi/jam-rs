# Trace sensitivity profiles

Trace sensitivity is an execution profile with explicit resource bounds. It
is not a calibrated probability of plasmid recovery and it should not be
reported as a biological sensitivity estimate. The profile is recorded in the
`run_header` so a result can be reproduced.

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

The primary scale is a FracMinHash scale: smaller values retain a denser
subset. The rescue level is attempted only when primary positional evidence is
insufficient. A denser archive can serve a sparser query level, subject to the
archive's available k and scale sections. A missing or incompatible seed
level is an explicit error, not a zero-evidence result.

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

## Archive settings

Build JMA v1 with seed levels that cover the profiles you intend to run. For
the shipped defaults:

```bash
jam archive \
  --input assembly.fasta \
  --output assembly.jma \
  --primary-scale 200 \
  --rescue-scale 500
```

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
divergent or short represented sequence at the cost of more repetitive
occurrences. The two levels serve different retrieval roles; neither makes a
sketch or alignment a final plasmid-origin call. This release preserves the
fixed k=31/k=21 contract and does not restart classifier training or add new
collision analysis below k=21.

## Resource controls

Global `--threads` bounds parallel candidate processing. Global `--memory`
sets the resource cache budget used by `jam trace`; seed occurrence, anchor,
chain, alignment-window, retained-alignment, and JSON output buffers are
bounded by the selected profile and runner configuration. The value is an
operational budget, not a hard process RSS limit because mmap pages, allocator
arenas, thread stacks, parser buffers, and the runtime can add overhead.

For a reproducible run, preserve the complete `run_header`, including the
profile fields, command, database checksum, catalog checksum, source commit,
and input checksums. See [TRACE_JSON.md](TRACE_JSON.md).
