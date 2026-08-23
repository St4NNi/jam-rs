# One-plasmid-to-many-metagenome trace search

The trace command answers a bounded screening question: does one cataloged
plasmid have sketch and positional sequence evidence in any of a collection of
metagenomic assemblies? It is a candidate finder for known and near-known
plasmid traces. It does not call plasmid presence, assign a host, assemble a
plasmid, or replace read mapping and assembly-graph review.

The workflow is:

```text
one plasmid FASTA
    -> existing .jam candidate index
    -> catalog lookup
    -> optional JMA positional archive (or raw assembly fallback)
    -> bounded seed/chaining and alignment
    -> JSONL evidence for downstream confirmation
```

The direction is intentional: the input is exactly one plasmid record and the
database contains one sample per metagenome. A catalog row identifies the
same sample ID stored in the `.jam` index and points to its JMA archive, raw
assembly, or both.

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
  --memory 4
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

Create a catalog beside the archive directory. Relative resource paths are
resolved relative to the catalog file:

```text
metagenome_id	jma	raw
sample_001.fasta	archives/sample_001.jma	../assemblies/sample_001.fasta
sample_002.fasta	archives/sample_002.jma	../assemblies/sample_002.fasta
```

Run one plasmid query across all catalog candidates:

```bash
jam trace \
  --plasmid plasmids/p_001.fasta \
  --database metagenomes.jam \
  --catalog metagenomes.tsv \
  --output p_001.trace.jsonl.zst \
  --sensitivity balanced \
  --min-shared 3 \
  --min-plasmid-containment 0.0 \
  --min-metagenome-containment 0.0 \
  --threads 8 \
  --memory 4
```

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

## Candidate and positional stages

The existing `.jam` index is the first, fast retrieval stage. It reports the
number of shared retained hashes and both containment directions:

```text
plasmid_containment     = shared_hashes / plasmid_hashes
metagenome_containment  = shared_hashes / metagenome_hashes
```

These denominators answer different questions and are not interchangeable.
`--min-shared`, `--min-plasmid-containment`, and
`--min-metagenome-containment` are independent filters applied before the
candidate limit. Candidates sort deterministically by score, shared hashes,
metagenome containment, and metagenome identifier; `rank` is one-based within
the retained candidate list. The `.jam` stage remains sketch evidence.

For each retained candidate, the runner prefers JMA when a catalog row has
one. It validates the JMA header, section directory, checksums, seed identity,
and sequence ranges, then performs bounded positional seeding, chaining, and
banded alignment. If JMA fails and a raw resource is present, the runner keeps
the JMA failure as a structured record and retries the candidate through the
raw resource. A catalog row with only raw data uses that path directly.

JMA v1 stores k=31 primary and optional k=21 rescue positional seeds. The
primary k=31 value is a specificity-oriented workflow contract, not a claim
that k=31 is universally optimal. The k=21 rescue level helps short or more
divergent represented traces. Profiles and their bounds are described in
[SENSITIVITY.md](SENSITIVITY.md); no new below-21 collision analysis is part
of this release.

## Reading evidence

An `alert` means a `.jam` candidate passed the configured sketch thresholds.
`supported` means the retained candidate also has positional alignment and
coverage evidence in the selected resource. Coverage is nonredundant query
coverage: overlapping contigs and secondary alignments cannot count the same
plasmid base twice. A split trace may therefore have a strong assembly-level
coverage value without proving that its fragments are linked on one molecule.

Use `supported_bases`, `supported_fraction`, primary intervals, secondary
overlap, gaps, alignment identity, and both containment directions together.
For a circular plasmid, an origin-crossing alignment has two ordered query
segments and `origin_crossing=true`; reported target coordinates remain
forward stored-contig coordinates even on the reverse strand.

`confirmed` is not emitted by jam-rs as a biological conclusion. A downstream
workflow may assign it only after independent evidence such as read mapping,
coverage breadth/depth, marker consistency, long-read linkage, or
assembly-graph context has been reviewed. The example's confirmation command
is a hand-off point, not an implementation of mapping.

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

The JMA reader validates its header, directory, and payload sections when
opened; the current v1 implementation decodes those sections before serving
sequence-range requests. Raw resources are streamed through the FASTA/FASTQ
parser, so a raw fallback can read the complete assembly resource. A JMA
archive still provides checked sequence-range semantics; this release does not
claim remote window-only transfer.

`--threads` bounds candidate processing, while global `--memory` sets the
resource/cache budget used by the command. Seed occurrences, anchors, chains,
alignment windows, retained alignments, and output records are bounded. The
memory value is an operational budget, not a hard RSS ceiling: allocator,
thread-stack, parser, mmap, and operating-system overhead remain.

## Reproducibility and compatibility

Record the catalog, `.jam` sidecar, JMA source checksums, raw-resource versions,
trace JSONL, and command line together. JMA archives are format version 1;
existing `.jam` files remain version 3 and are not rewritten by `jam trace`.
The trace JSON schema is version `1.0.0`. A legacy `.jam` without a sidecar may
still be readable, but its original catalog provenance cannot be reconstructed.

Bias-assisted `.jam` retrieval is optional and catalog-specific. In bias mode,
retained and weighted evidence is labelled separately, and the uniform-hash
E-value is `null`/`NA`; it is not calibrated significance. The uniform index
remains the baseline. Bias tables are tied to their positive set, chromosome
background, k-mer size, and sampling configuration.
