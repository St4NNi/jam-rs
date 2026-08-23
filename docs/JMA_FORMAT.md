# JMA v1 archive format

JMA is the positional archive used by the trace workflow. It is separate from
the existing `.jam` sketch database: `.jam` remains format version 3 and is
never parsed as JMA.

## Layout

All integers are unsigned little-endian values. The writer never serializes a
Rust structure directly, and the reader checks every addition, conversion,
offset, length, and slice before using it.

```text
fixed header (160 bytes)
section directory (64 bytes per entry)
contig table and UTF-8 names
packed sequence blocks
k=31 seed section (optional only for an explicitly incomplete archive)
k=21 seed section (optional only for an explicitly incomplete archive)
```

The fixed header contains the magic `JMA\0`, format version `1`, flags, contig
count, total bases, source SHA-256, directory offset and length, a SHA-256
header checksum, and up to two seed-level identities. New archives also store
the existing `jam-seed-chain-align-v1` tag, algorithm version `1`, and the
optional finite entropy threshold used for seed construction. That binary tag
is retained for JMA v1 compatibility; the generic query workflow identifiers
are recorded in the sidecar described below. Header checksum coverage includes
those metadata bytes with the checksum field zeroed.

Legacy JMA v1 headers with zeroed metadata bytes can still be decoded, but
production `jam trace` requires the algorithm tag and version so an archive
with unknown positional semantics is not silently interpreted as current
evidence. An explicit entropy threshold of `0.0` is distinct from an absent
threshold.

Each directory entry contains a section kind, flags, byte offset, stored
length, uncompressed length, and a SHA-256 checksum. Version 1 stores sections
uncompressed, so stored and uncompressed lengths must match. Section ranges may
not overlap the header, directory, or another section.

## Contigs and sequences

The contig section stores strictly increasing IDs, base lengths, and checked
UTF-8 names. Sequence data is divided into sorted blocks. Each block stores its
contig ID, zero-based start, base length, two-bit packed bases, and a one-bit
ambiguity mask. A set ambiguity bit decodes as `N`; lowercase bases are
uppercased and `U` is treated as `T`. A reader rejects an uncovered range or
conflicting overlapping blocks.

Ranges use the workflow-wide zero-based half-open convention `[start, end)`.

## Seeds

Seed sections contain a k-mer size, one or more explicit density levels, and
fixed-size records with:

```text
hash, exact packed canonical k-mer, contig ID, orientation, position
```

The exact packed k-mer is compared after the hash. A hash collision therefore
cannot by itself create a seed occurrence. The frozen production identities
are k=31 primary seeds and optional k=21 rescue seeds. Density levels are
ordered from densest (lowest scale) to sparsest and can be selected explicitly
by scale. The default lookup uses the densest level represented by the
archive.

## Integrity and compatibility

Header and every loaded section have SHA-256 checksums. Corrupt magic,
versions, lengths, offsets, section overlap, names, sequence blocks, seed
records, or checksums are reported as JMA errors; malformed data is never
reported as an empty archive. The source checksum identifies the input used to
build the archive and is retained as provenance; it is not a checksum of the
JMA container itself.

JMA v1 is a new format and does not alter or migrate existing `.jam` files.
Future incompatible JMA changes must use a new format version.

## Range index

`jam archive` also writes `archive.jma.idx.json`. The sidecar does not change
JMA v1 bytes. It records the complete archive SHA-256 and is bound at runtime
to the archive's size, fixed-header SHA-256, source SHA-256, format version,
and algorithm identity. It contains:

- absolute byte ranges and SHA-256 values for each encoded sequence block;
- absolute fixed-record ranges and SHA-256 values for nonempty 16-bit hash
  prefix buckets at each `(k, scale)` seed level;
- the contig, seed-level, and algorithm identities needed to validate those
  ranges against the JMA header and section directory;
- the additive `workflow_identifiers` object:

  ```json
  {
    "screen_algorithm": "jam-fracminhash-screen-v1",
    "local_alignment_algorithm": "jam-exact-seed-chain-banded-v1",
    "mosaic_algorithm": "jam-fragment-mosaic-v1",
    "trace_workflow": "jam-trace-v1"
  }
  ```

These identifiers describe the screening, local alignment, fragment mosaic,
and end-to-end workflow layers. They do not alter the JMA header or payload.

Indexed opening follows this bounded sequence:

```text
fixed JMA header
    -> section directory
    -> contig table
    -> checksum-bound sidecar
    -> seed buckets for the current query seeds
    -> accepted seed chains
    -> sequence blocks intersecting those chains
    -> optional buckets and blocks for unresolved-gap rescue
```

A seed query reads only its matching prefix bucket and still compares the
complete hash and packed canonical k-mer. Sequence reads fetch only blocks
intersecting an accepted chain window. Every fetched range is checked against
its sidecar checksum before decoding. Local JMA resources use positional file
reads through this same range contract; the local `.jam` candidate index is
materialized for memory-mapped sketch search.

The run and candidate records expose the corresponding resource counters:
`metadata_requests`, `head_requests`, `get_requests`, `range_requests`,
`stream_requests`, `requested_bytes`, `returned_bytes`, `decoded_bytes`,
`cache_bytes`, `retries`, `full_object_fallbacks`, `seed_buckets_read`, and
`sequence_blocks_read`. Gap-rescue records additionally report buckets
requested, keys tested, anchors, chains, blocks fetched, alignment windows, and
newly supported query bases per round.

The sidecar is bound to the archive size, header checksum, source checksum,
format version, and binary algorithm identity. A legacy sidecar that lacks
`workflow_identifiers` remains JSON-decodable for inspection, but indexed
opening rejects it with a compatibility error because its workflow provenance
is incomplete. Rebuild the sidecar with the current archive builder. A sidecar
with incompatible identifiers is rejected before any payload range is used.

A JMA v1 archive without a sidecar remains decodable through the eager reader.
`jam trace` uses that path only when full-download fallback is enabled and
records `jma_full_download_fallback`; strict mode emits
`jma_index_required` instead.

Seed density is part of compatibility. A requested profile may use an equal or
denser stored seed level, but a sparser archive produces
`seed_level_mismatch` and no raw positional evidence.
