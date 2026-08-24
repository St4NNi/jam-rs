# Remote trace resources

The trace catalog can point to local files or remote assembly resources. The
resource layer gives JMA and raw FASTA/FASTQ readers one checked interface for
metadata, byte ranges, streams, retries, and cache accounting.

This document applies to the existing `.jam` plus JMA/raw trace path. Jam
Index build, append, and trace operation is local-only: `jam trace --index`
accepts a local directory, memory maps its parts, and does not use remote
ranges, persistent page caches, object uploads, JMA resources, or full-object
fallback. See [JAM_INDEX.md](JAM_INDEX.md).

Jam Index parts reference local external assemblies rather than embedding
sequence. Plain FASTA uses FAI ranges, BGZF FASTA uses FAI/GZI ranges, and
normal gzip uses one candidate-only sequential stream. No assembly is opened
for a noncandidate metagenome. These reads do not use the remote resource
contract described below.

## Locator syntax

Each catalog resource is one of:

```text
assemblies/sample_001.fasta
file:///data/assemblies/sample_001.fasta
https://objects.example.org/metagenomes/sample_001.jma
s3://surveillance-catalog/metagenomes/sample_001.jma
```

Production JMA rows also provide the checksum-bound `jma_index` object. The
index identifies exact archive ranges and is not guessed from a signed archive
URL.

Relative local paths are resolved relative to the catalog file. Absolute local
paths and `file://` paths are accepted as written. HTTP, HTTPS, and S3
locators require a host/bucket and object key as appropriate.

The transport implementation uses `curl` for HTTP(S). An `s3://bucket/key`
locator is mapped to an HTTPS endpoint (or `JAM_S3_ENDPOINT` for an
S3-compatible service); jam-rs does not embed an AWS SDK or perform SigV4
signing. Use a public object, a pre-authorized endpoint, or an HTTPS signed URL
supplied by the surrounding pipeline. Credentials should not be embedded in
catalogs or logs. Signed URL query and fragment values are removed from the
locator recorded in run output and errors.

## Remote candidate database

The `.jam` candidate index must be a local file before `QueryEngine` can mmap
it. `jam trace` bridges a remote locator through an explicit disk cache:

```bash
jam trace \
  --query p.fasta \
  --query-kind plasmid \
  --database https://objects.example.org/metagenomes.jam \
  --cache-dir cache/jam \
  --catalog metagenomes.tsv \
  --output p.trace.jsonl
```

The cache key includes the redacted locator, ETag or Last-Modified token, and
declared size. A cache sidecar records the downloaded SHA-256, and every reuse
checks size and checksum before mmap. An incomplete, identity-mismatched, or
corrupt cache entry fails explicitly and is not overwritten automatically.
When the server supplies no version token, jam-rs downloads the object before
selecting a content-addressed entry; it does not assume that a same-sized
object is unchanged.

Trace output remains a local writable path. Upload or object-store publication
belongs to the surrounding pipeline so an interrupted JSONL stream cannot be
mistaken for an atomically completed remote object.

## Catalog rows

The minimum TSV form is:

```text
metagenome_id	jma	jma_index	raw
sample_001	https://objects.example.org/m/sample_001.jma	https://objects.example.org/m/sample_001.jma.idx.json	https://objects.example.org/m/sample_001.fasta.gz
sample_002	s3://surveillance-catalog/m/sample_002.jma	s3://surveillance-catalog/m/sample_002.jma.idx.json	s3://surveillance-catalog/m/sample_002.fasta.gz
sample_003			file:///data/m/sample_003.fasta
```

The equivalent JSON is either an array or an object with an `entries` array:

```json
[
  {
    "metagenome_id": "sample_001",
    "jma": "https://objects.example.org/m/sample_001.jma",
    "jma_index": "https://objects.example.org/m/sample_001.jma.idx.json",
    "raw": "https://objects.example.org/m/sample_001.fasta.gz"
  },
  {
    "metagenome_id": "sample_002",
    "jma": "s3://surveillance-catalog/m/sample_002.jma",
    "jma_index": "s3://surveillance-catalog/m/sample_002.jma.idx.json"
  }
]
```

`metagenome_id` must be unique and must match the sample name in the `.jam`
candidate database. At least one of `jma` or `raw` is required. When both are
present, JMA is preferred. A permitted resource-level JMA failure is retained
as a structured failure before an explicit raw fallback is attempted.
Structural, checksum, format, sidecar-identity, and seed-level compatibility
failures do not trigger unanchored raw fallback. A raw-only row is valid but
does not provide indexed archive access.

## Range and streaming behavior

The resource contract uses checked, zero-based, half-open byte ranges. Before
allocation or transport, offset-plus-length overflow and out-of-bounds ranges
are rejected. Metadata records object size, version information (ETag or
last-modified when available), and whether ranges are advertised.

HTTP(S) and S3 readers retry a bounded number of failed requests. A ranged
response must be `206 Partial Content` with the exact requested
`Content-Range`. A `200` response is accepted only through the configured,
size-bounded full-object fallback. When range support is unavailable and
fallback is disabled, the candidate is a structured failure rather than an
empty result. The fallback is explicit in the result counters and does not
make a complete-object read look like selective access.

Indexed JMA opening uses this request sequence:

```text
fixed header -> section directory -> contig table -> sidecar
    -> query-intersecting seed buckets
    -> accepted chains
    -> chain-implicated sequence blocks
    -> optional unresolved-gap rescue buckets and blocks
```

The local implementation performs positional reads from the archive and
sidecar files. HTTP and S3 use the same byte-range contract. They must not read
the complete seed payload, complete sequence payload, unrelated contigs, or a
noncandidate resource in indexed mode. Resource metrics report
`metadata_requests`, `head_requests`, `get_requests`, `range_requests`,
`stream_requests`, `requested_bytes`, `returned_bytes`, `decoded_bytes`,
`cache_bytes`, `retries`, `full_object_fallbacks`, `seed_buckets_read`, and
`sequence_blocks_read`. Gap rescue adds per-round counts for seed buckets,
tested keys, anchors, chains, sequence blocks, alignment windows, and newly
supported query bases.

Raw FASTA/FASTQ resources are streamed into the parser and may require the
complete input stream; the current raw candidate path retains parsed contigs,
so it may also materialize the complete assembly in memory. Semantic
seed-level, format, checksum, and sidecar-identity failures are never eligible
for raw fallback. A resource transport failure may use a catalog raw fallback
only where the configured failure policy permits it; the failure is retained
in the result.

## Cache and redaction rules

Cache identities combine the redacted locator, resource version token, and
object size. If a remote object changes identity, old blocks are invalidated
instead of being mixed with the new object. Indexed JMA blocks are additionally
bound to their exact byte range and sidecar checksum. Cache block size, total
cache size, timeout, retry count, and full-download fallback are part of the
resource opening contract. `--cache-block-bytes`,
`--request-timeout-seconds`, `--max-retries`, and
`--no-full-download-fallback` expose those controls; global
`--memory-target` derives the in-memory cache target. `--cache-dir` is the disk
location used only when a remote candidate `.jam` must be materialized.

Indexed archive blocks are cached by redacted locator, object version, size,
and exact byte range. A sidecar whose archive identity or block checksum does
not match is rejected before its content is used.

Only the redacted locator is safe for logs and JSON. The redaction removes URL
user-info and query/fragment suffixes, for example:

```text
https://user:secret@example.org/a.jma?X-Amz-Signature=...
    -> https://example.org/a.jma
```

Do not use the redacted value as an access credential. Keep the original
catalog and deployment secret configuration under the pipeline's access
controls, and archive the trace output without reconstructing signed URLs.

## Operational guidance

For reproducibility, store the catalog, JMA/raw object version or checksum, the
`.jam` sidecar, the trace command, and the JSONL output together. Test remote
access with a small archive before a large run. A transport failure, checksum
failure, cache identity change, or malformed resource should remain visible in
the `failures` array and must not be interpreted as a negative biological
result. See [TRACE_JSON.md](TRACE_JSON.md) for consumer behavior.

## Current measurement limits

The current snapshot includes local and mock range-read evidence, but no actual
S3 cold-cache or warm-cache measurement. S3 deployment should preflight range
support, object size/version metadata, endpoint behavior, credentials, and the
chosen cache block size before a production run. A block-size study on a small
fixture does not establish the best value for a large object collection.

The indexed path is selective only when the JMA sidecar is present and the
object honors ranges. Without a sidecar, strict mode stops with
`jma_index_required`; explicit full-download fallback is bounded and recorded.
Raw FASTA/FASTQ fallback is a separate complete-stream path and does not provide
indexed byte selectivity. None of these fallbacks turns sketch or alignment
evidence into a biological presence call.
