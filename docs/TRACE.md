# Query-centered fragment trace search

`jam trace` searches one query sequence element across many metagenomes. The
query can be a plasmid, phage, synthetic construct, other mobile genetic
element, or an element whose class is unknown. The command is a reference-
guided candidate and fragment-evidence workflow, not a plasmid or phage
classifier, assembler, binning method, or proof that an element is autonomous.

The workflow is:

```text
one query FASTA record
    -> existing .jam candidate index
    -> metagenome catalog lookup
    -> indexed JMA resource or raw assembly fallback
    -> exact positional seeds and contig-specific chains
    -> local base-level alignments
    -> nonredundant query-coordinate fragment mosaic
    -> JSONL/JSONL.zst evidence for independent confirmation
```

The query is exactly one FASTA/FASTQ record; zero or multiple records are
rejected with an input error. A metagenome catalog row maps the sample
identifier stored in the `.jam` index to a JMA archive, its checksum-bound
index, a raw assembly, or both. Separate supporting contigs remain separate
observations; a mosaic does not assert that they are physically linked.

## Quick start

Build a candidate index from the assemblies. With the default non-singleton
mode, each input file is one `.jam` sample, and the sample identifier is its
file name. This is appropriate when each file contains all contigs for one
metagenome:

```bash
jam sketch assemblies/sample_001.fasta assemblies/sample_002.fasta \
  --output metagenomes.jam \
  --kmer-size 21 \
  --fscale 100 \
  --threads 8 \
  --memory-target 4
```

Keep the generated `metagenomes.jam.json` sidecar with the `.jam` file. It
records the hash identity, k-mer and scale parameters, input identities, and
database checksum. The trace candidate query must use the same database
parameters; a mismatch is an error.

Build a positional archive for each assembly when indexed sequence access is
useful:

```bash
jam archive \
  --input assemblies/sample_001.fasta \
  --output archives/sample_001.jma \
  --primary-scale 200 \
  --rescue-scale 500 \
  --complexity 0.0
```

This writes `archives/sample_001.jma` and the checksum-bound range index
`archives/sample_001.jma.idx.json`. Keep and version them together.

Create a catalog beside the archive directory. Relative resource paths are
resolved relative to the catalog file:

```text
metagenome_id	jma	jma_index	raw
sample_001.fasta	archives/sample_001.jma	archives/sample_001.jma.idx.json	../assemblies/sample_001.fasta
sample_002.fasta	archives/sample_002.jma	archives/sample_002.jma.idx.json	../assemblies/sample_002.fasta
```

Run one query element across all catalog candidates:

```bash
jam trace \
  --query elements/p_001.fasta \
  --query-kind plasmid \
  --topology auto \
  --database metagenomes.jam \
  --metagenomes metagenomes.tsv \
  --output p_001.trace.jsonl.zst \
  --sensitivity balanced \
  --min-shared 3 \
  --min-query-containment 0.0 \
  --min-metagenome-containment 0.0 \
  --threads 8 \
  --io-concurrency 8 \
  --memory-target 4
```

`--metagenomes` is the primary catalog option; `--catalog` remains a visible
compatibility alias. `--query-kind` accepts `plasmid`, `phage`, `other`, and
`unknown`. It records the declared class and is not a classifier result.
`--topology` accepts `linear`, `circular`, `auto`, and `unknown`.

`--plasmid` remains a compatibility alias for one release cycle. It maps to
`--query <value> --query-kind plasmid` and emits a deprecation warning. New
scripts should use `--query`.

The `.zst` suffix selects the same JSONL record stream in a Zstandard frame.
Use a non-`.zst` path for plain UTF-8 JSONL. The run header records the command,
software version, source commit when available, sensitivity configuration,
and checksums for the local input files. See [TRACE_JSON.md](TRACE_JSON.md).

The repository example packages this sequence, archive construction, catalog
generation, candidate extraction, and a downstream confirmation command:

```bash
bash examples/trace/run_trace.sh plasmids/p_001.fasta results \
  assemblies/sample_001.fasta assemblies/sample_002.fasta
```

The script emits candidate metagenome IDs and supporting contig IDs. Its
`alert` and `supported` labels are jam-rs evidence labels; `confirmed` is
reserved for the downstream mapping step shown by the script.

Plasmid surveillance is the primary use case. The same query-centered
contract can be used for a linear phage reference by changing the declared
kind and topology:

```bash
jam trace \
  --query phages/phage_001.fasta \
  --query-kind phage \
  --topology linear \
  --database metagenomes.jam \
  --metagenomes metagenomes.tsv \
  --output phage_001.trace.jsonl \
  --sensitivity balanced \
  --threads 8 \
  --io-concurrency 8
```

The `phage` value records the query class; it does not establish that an
aligned fragment is an autonomous phage.

## Candidate and positional stages

The existing `.jam` index is the first, fast retrieval stage. It reports the
number of shared retained hashes and both containment directions:

```text
query_containment       = shared_hashes / query_hashes
metagenome_containment  = shared_hashes / metagenome_hashes
```

These denominators answer different questions and are not interchangeable.
`--min-shared`, `--min-query-containment`, and
`--min-metagenome-containment` are independent filters applied before the
candidate limit. Candidates sort deterministically by the configured evidence,
shared hashes, containment, and metagenome identifier; `rank` is one-based
within the retained candidate list. The `.jam` stage remains sketch evidence.

For each retained candidate, the runner prefers JMA when a catalog row has
one. It validates the JMA header, section directory, checksums, seed identity,
and sequence ranges, then performs bounded positional seeding, chaining, and
banded alignment. An incompatible JMA seed level is a structured
`seed_level_mismatch` failure and never falls back to unanchored raw alignment.
Permitted resource/transport failures may retain the JMA failure and retry
through a raw resource when one is present. Structural, checksum, format,
sidecar-identity, and seed-level failures remain structured failures without
unanchored raw fallback. A catalog row with only raw data uses that path
directly.

JMA v1 stores k=31 primary and optional k=21 rescue positional seeds. The
primary k=31 value is a longer exact-anchor workflow contract, not a claim that
k=31 is universally optimal. The k=21 rescue level helps short or more
divergent represented traces. Profiles and their bounds are described in
[SENSITIVITY.md](SENSITIVITY.md), and the exact seed, chain, alignment, and
coverage rules are specified in [ALGORITHM.md](ALGORITHM.md). No seed size below
21 is used by the trace workflow.

An archive must be at least as dense as the requested profile. Build k=31 at
scale 100 and k=21 at scale 200 when the same archive must support all three
shipped profiles; a sparser archive cannot satisfy `sensitive` and fails
closed.

## Query coordinates and topology

`topology_requested` records the command option. `coordinate_model` records
how query coordinates were handled:

| Coordinate model | Meaning |
| --- | --- |
| `linear` | Ordinary `[0, query_length)` coordinates; terminal gaps remain separate and no origin crossing is accepted. |
| `wrap` | A coordinate convention that may use unwrapped `[0, 2L)` chain coordinates and normalize the final union on a circle. |
| `undetermined` | Linear coordinates are retained for display, while the result does not select a topology model. |

`topology_evidence` is a separate evidence label: `linear_supported`,
`wrap_supported`, `both_compatible`, `insufficient`, or `undetermined`.
`coordinate_model=wrap` is never a biological claim that the molecule is
circular. It can also represent a circularly permuted or arbitrarily rotated
query representation.

`linear` prohibits origin-crossing chains and keeps leading/trailing gaps
distinct. `circular` selects the wrap coordinate model. `auto` builds linear
and wrap summaries from the same accepted local alignments and selects a model
only when the configured evidence margin is exceeded; otherwise it retains an
undetermined result. `unknown` keeps the primary result undetermined and does
not claim a topology, while possible wrap evidence can remain in the model
summary for inspection.

## Fragment mosaic and evidence fields

The final unit of evidence is one metagenome-level query mosaic. Every accepted
local alignment remains in `alignments`, including overlapping and alternative
support. The primary fragment mosaic splits query coordinates at alignment
boundaries, assigns each atomic interval a deterministic primary alignment,
and keeps other alignments as secondary evidence. It can combine different
query regions from different contigs, but it never treats those contigs as
physically joined.

The result distinguishes:

- `base_covered_bases`: query bases paired to a metagenome base by CIGAR `=` or
  `X`; the union counts each query base once.
- `aligned_span_bases`: query bases inside accepted alignment spans, including
  query bases crossed by a CIGAR deletion.
- `unsupported_gaps`: query intervals with no accepted base support.
- `alignment_deletions`: query intervals crossed by a deletion inside an
  accepted alignment.
- `common_sequence_intervals` and `repeat_only_intervals`: support associated
  with common or repetitive seed evidence.

Unsupported gaps are not alignment deletions. Overlapping or duplicate contigs
cannot inflate unique query coverage. Common transposons, resistance cassettes,
terminal repeats, and chromosome-shared regions may create candidate evidence;
an alignment supported only by common/repetitive sequence is labelled in the
alignment role and requires extra interpretation.

The example workflow uses `alert` for a sketch candidate, `supported` for
candidate plus positional evidence, and `confirmed` only for an independently
reviewed downstream result. Jam-rs itself does not emit `confirmed` as a
biological presence call. Read mapping, breadth/depth, marker checks,
long-read linkage, assembly-graph context, or another independent method is
required before assigning that label. Sequence support also does not prove an
autonomous element or a host association.

## Resource forms and limits

Catalog resources accept local paths, `file://`, `http://`, `https://`, and
`s3://bucket/key` locators. Local paths in a catalog are resolved relative to
the catalog file. HTTP(S) and S3 access uses the shared range/stream resource
layer, retries bounded requests, and caches checked blocks within the memory
budget. Signed query parameters and user-info credentials are removed from
metadata, errors, and cache identities shown to users. See
[REMOTE_RESOURCES.md](REMOTE_RESOURCES.md).

The candidate `.jam` is always memory mapped locally. A local path or `file://`
locator is opened directly. HTTP(S) and S3-compatible locators require
`--cache-dir`; jam-rs downloads the complete database, validates its size and
SHA-256, stores it under a redacted locator/version identity, and then mmaps
the cached file. Objects without an ETag or Last-Modified token are downloaded
again before a content-addressed cache entry is selected, rather than being
silently reused by locator and size alone.

With a `jma_index` catalog resource, the reader validates the JMA header,
directory, contig table, and checksum-bound sidecar, then fetches only seed
buckets matching query seeds and sequence blocks intersecting accepted chains.
Gap rescue may request additional buckets and blocks only for unresolved query
regions. Raw resources enter through the streaming FASTA/FASTQ parser, but the
current raw candidate path retains parsed contigs and can materialize the
complete assembly in memory. A JMA row without an index
uses eager full-section access only when full-download fallback is enabled and
records that fallback explicitly.

`--threads` bounds candidate processing and `--io-concurrency` bounds candidate
resource tasks. Global `--memory-target` sets the resource/cache budget used by
the command. Seed occurrences, anchors, chains, alignment windows, retained
alignments, and output records are bounded. The memory target is not a hard RSS
ceiling: allocator, thread-stack, parser, mmap, and operating-system overhead
remain.

## Reproducibility and compatibility

Record the catalog, `.jam` sidecar, JMA source checksums, raw-resource versions,
trace JSONL, and command line together. JMA archives are format version 1;
existing `.jam` files remain version 3 and are not rewritten by `jam trace`.
The trace JSON schema is version `2.0.0`. A legacy `.jam` without a sidecar may
still be readable, but its original catalog provenance cannot be reconstructed.

Bias-assisted `.jam` retrieval is optional and catalog-specific. In bias mode,
retained and weighted evidence is labelled separately, and the uniform-hash
E-value is `null`/`NA`; it is not calibrated significance. The uniform index
remains the baseline. Bias tables are tied to their positive set, chromosome
background, k-mer size, and sampling configuration.
