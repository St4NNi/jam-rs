# JMA format 1

JMA format 1 is a deterministic, self-contained archive for one metagenome.
One catalog row points at one `.jma` object and its checksum. Query-time trace
does not use a companion index or a text file.

Jam Index is a separate local format. It contains neither JMA objects nor
nucleotide sequence and does not open JMA objects.
Its complete selected contigs generate positions in memory. See
[JAM_INDEX.md](JAM_INDEX.md). This document applies only to the existing
`.jam` plus JMA trace path.

## Object envelope

The first 256 bytes are the fixed little-endian superblock. It contains the
format identifier (`JMAF1\0\0\0`), format version `1`, layout identifier
`0x31414d4a`, flags, section count, directory offset and length, object size,
contig count, total bases, the source-assembly SHA-256, an object checksum,
the section-directory checksum, and the superblock checksum. All offsets and
lengths are unsigned 64-bit half-open byte coordinates. The checksum fields
are finalized without recursively including themselves, so archive output is
deterministic.

The section directory contains 64-byte records:

```text
kind u32, flags u32, offset u64, stored_length u64,
decoded_length u64, sha256[32]
```

Required records are unique. A required unknown kind, overlapping range,
range outside the declared object, integer overflow, or duplicate required
kind is rejected. Optional unknown kinds are skipped. Every known section is
checksummed before it is decoded. The superblock's section-directory checksum
covers the complete directory bytes, so descriptors cannot be changed without
invalidating the superblock. The seed-section descriptor checksum intentionally
covers only the seed directory prefix (scheme records and page records); the
key and occurrence payloads are independently checked by each selected page.
Other section descriptor checksums cover their complete section payload.

## Required sections

The writer emits these sections in deterministic order:

1. metadata: `JMA`, version/layout identifiers, source checksum, builder
   version, source commit, `jamhash_u64_v1`, and the optional entropy filter
   setting used to create positional seeds.
2. contigs: fixed-width contig records followed by one UTF-8 name string
   table. Names are not parsed by seed lookup.
3. embedded sketch: sorted unique retained hash values.
4. seeds: an open list of scheme descriptors, fixed page directories, fixed
   key records, and fixed occurrence records.
5. sequence directory: one fixed record per independent block.
6. sequence payload: complete two-bit assembly blocks and IUPAC ambiguity
   records.

An optional gear directory can be added without changing the required section
set. Unknown optional sections are safely ignored.

## Seed pages

The seed section begins with its own directory. Each scheme descriptor records
its numeric scheme identity, algorithm identity, span, informative-base count,
density parameter, bucket width, key encoding, occurrence encoding, and flags.
There is no two-scheme header limit. Every seed page contains a separate key
range and occurrence range. Keys store the digest, exact packed canonical
selected bases, and an occurrence index. Occurrences store contig, position,
span, and strand.

A lookup selects pages by the high hash prefix, reads only those key and
occurrence ranges, verifies the page checksum, then compares the exact packed
key. A digest collision can increase work but cannot create an anchor.

## Sequence blocks

Sequence directory records contain contig ID, base start/count, absolute object
payload offset, stored and decoded byte lengths, codec, ambiguity range, and checksum.
The shared sequence codec writes `raw2bit` or `zstd2bit` blocks: A/C/G/T use
two bits and every non-ACGT base remains an IUPAC ambiguity record. Blocks are
independently addressable and cannot cross a contig boundary. The directory is
read when the archive opens; payload and adjacent ambiguity bytes are fetched
only for requested sequence ranges.

One archive may mix `raw2bit` and `zstd2bit` blocks. Fixed-size and Gear
content-dependent boundary policies use the same directory and per-block
checksum rules; their boundary parameters are independent of Gear seed
fragment parameters.

## Access model

Opening a range resource reads the superblock, section directory, metadata,
contig/name table, seed directory, and sequence directory. It does not read
seed occurrences or sequence payloads. Local readers map the object once and
borrow selected section/page slices; remote readers use bounded byte ranges.
Selected key/occurrence pages and selected sequence blocks are fetched on
demand, and counters expose mapped, resident, decoded, and transferred bytes
plus request counts.

JMA stores complete contigs but does not imply physical linkage between
separate contigs. Trace evidence is not proof of an autonomous plasmid or
phage.
