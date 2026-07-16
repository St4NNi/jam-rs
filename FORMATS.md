<!-- The byte layouts and JSON fields of the files of a shard and of an index. -->
<!-- Copyright (c) 2023-2026 Sebastian Beyvers -->
<!-- SPDX-License-Identifier: MIT -->

# jam file formats

A shard is one directory with `manifest.json`, `sketch`, `positions` and `packed/`. An index of
several shards adds `index.json`. Binary fields are little endian.

## manifest.json

The manifest of a shard: `manifest.json` with the parameters, the entries, and the size and
SHA-256 digest of every file.

| Field | Meaning |
| --- | --- |
| `format_version` | 1 |
| `parameters` | `k`, `w`, `min_entropy`, `bin_bits`, `bucket_bits` |
| `entries` | in entry id order: `name`, `contigs`, `bases`, `packed` (file) |
| `sketch`, `positions` | files |
| `seeds`, `color_sets`, `seed_locations` | totals of the sketch and the positions |

A file is `path` (relative to the shard directory; packed files to the packed directory),
`bytes` and `sha256` (lower case hex).

## index.json

The index file of a resharded index: `index.json` in the index directory, listing its shards
in order with the entries each one holds, so a search or a server can open the whole index.

| Field | Meaning |
| --- | --- |
| `format_version` | 1 |
| `shards` | in order: `directory` (relative to the index directory), `first_entry`, `entries`, `seed_locations`, `manifest_sha256` |

The entries of a shard are `first_entry..first_entry + entries` in index order; entry `i` of
the shard is entry `first_entry + i` of the index. `manifest_sha256` is the SHA-256 digest of
the shard's `manifest.json`, so a replaced shard is noticed.

## sketch

The sketch file (stage 1): for every stored seed its color set, found from the seed's hash
with a few memory reads, and every distinct color set once. The file is memory-mapped and kept
in RAM.

A seed's hash selects a bucket by its top `bucket_bits` bits; the next 16 bits are its
fingerprint. A build chooses the fewest bucket bits (8 to 28) that hold at most 256 seeds per
bucket on average, counting seed locations as an upper bound of the seeds. A query k-mer then
matches a stored seed by chance about 0.4% of the time, and the directory costs about 0.03 bytes
per seed. Seeds with the same bucket and fingerprint are one stored seed, whose color set is the
union of theirs. Color set ids are ranked by frequency: the set of the most seeds is 0.

Layout, little endian; every section starts at a multiple of 8 bytes:

| Bytes | Field |
| --- | --- |
| 0..16 | magic `jam sketch` padded with zero bytes |
| 16..20 | version (u32, 1) |
| 20 | bucket bits b (u8, 8 to 28) |
| 21 | fingerprint bits (u8, 16) |
| 22..24 | zero |
| 24..32 | entries (u64) |
| 32..40 | stored seeds (u64) |
| 40..48 | color sets s (u64) |
| 48..56 | bucket block bytes (u64) |
| 56..64 | color set bytes (u64) |
| 64..80 | CRC-32 of the four sections below, in order (u32 each) |
| 80..124 | zero |
| 124..128 | CRC-32 of bytes 0..124 (u32) |
| directory | 2^b + 1 offsets (u64) of the bucket blocks, relative to the first block |
| bucket blocks | one block per non-empty bucket, then zero bytes to a multiple of 8, at least 8 |
| color set offsets | s + 1 offsets (u64), relative to the first color set |
| color sets | per set: entry count, first entry, then each gap minus one, all varints |

A bucket block of n seeds:

| Bytes | Field |
| --- | --- |
| 0..4 | n (u32) |
| 4..8 | escape count x (low 24 bits) and id width w (high 8 bits, 0 to 32) |
| 2n | the fingerprints (u16), ascending |
| 4x | escaped color set ids (u32) |
| ceil(n w / 8) | one w-bit value per seed, bit-packed, first value in the low bits |

A value v below 2^w - x is the color set id itself; a larger v names escaped id
v - (2^w - x). Each bucket chooses w and x for the fewest bits, so a few rare color sets do not
widen every value.

## positions

The positions file (stage 2): for every entry a section with its seed locations, so a seed of
one candidate entry is found with two small reads, and the contig starts of the entry, so
alignment windows can be cut at contig ends without the packed sequence. It works
memory-mapped from local disk and as byte ranges from S3.

A seed is stored by its key: the top 44 bits of its hash (`POSITION_KEY_BITS`), at least as
many as the sketch keeps (bucket bits and fingerprint), so the sketch of any size is built
from the positions alone. A seed location is stored as the bin within the entry, shifted left
by one, with the low bit set when the canonical k-mer is the reverse complement of the contig
strand.

Layout, little endian:

| Bytes | Field |
| --- | --- |
| 0..16 | magic `jam positions` padded with zero bytes |
| 16..20 | version (u32, 2) |
| 20 | key bits (u8, 44) |
| 21 | bin bits (u8) |
| 22..24 | zero |
| 24..32 | entries e (u64) |
| 32..36 | CRC-32 of the section table (u32) |
| 36..124 | zero |
| 124..128 | CRC-32 of bytes 0..124 (u32) |
| 128.. | section table: e rows of 40 bytes, see below |
| sections | one per entry, in entry order, each at a multiple of 8 |

A section table row:

| Bytes | Field |
| --- | --- |
| 0..8 | section offset in the file (u64) |
| 8..16 | section length (u64) |
| 16..24 | seed locations (u64) |
| 24..28 | CRC-32 of the section (u32) |
| 28 | slot bits d (u8) |
| 29 | suffix bits s = key bits - d (u8) |
| 30 | location bits p (u8) |
| 31 | zero |
| 32..40 | contigs c (u64) |

A section with c contigs and r seed locations, sorted by key and then seed location:

| Bytes | Field |
| --- | --- |
| 8 (c + 1) | contig starts in entry coordinates, then the entry's bases (u64 each) |
| 4 (2^d + 1) | slot directory: index of the first seed location of each slot (u32); a key's slot is its top d bits |
| ceil(r (s + p) / 8) | one record per seed location, bit-packed, first record in the low bits: the low s bits of the key, then the seed location in p bits |
| 8 | zero |

d is chosen so a slot holds about 16 seed locations; a lookup reads 8 bytes of the directory and
the records of one slot. A seed that occurs at several bins has one record per bin.

## packed sequence

The packed sequence of one entry (stage 3): its contigs at 2 bits per base, with contig names
and the runs of bytes other than A, C, G and T. Lower case a, c, g and t are stored as upper
case; other bytes keep their value, in upper case.

Layout, little endian ("entry coordinates" count bases over all contigs in order):

| Bytes | Field |
| --- | --- |
| 0..16 | magic `jam-rs 2bit seq` and a zero byte |
| 16..20 | version (u32, 1) |
| 20..24 | zero |
| 24..32 | contigs c (u64) |
| 32..40 | bases n (u64) |
| 40..48 | runs r of other bytes (u64) |
| 48..56 | name bytes m (u64) |
| 56..64 | offset of the packed bases (u64) |
| 64..68 | CRC-32 of the tables, from byte 128 to the packed bases (u32) |
| 68..72 | CRC-32 of bytes 0..68 (u32) |
| 72..128 | zero |
| 8 (c + 1) | contig starts in entry coordinates, then n (u64 each) |
| 8 c | end of each contig name in the name bytes (u64 each) |
| m, padded to 8 | contig names, the first word of each FASTA header |
| 16 r | runs: start (u64), length (u32), byte (u8), three zero bytes |
| 4 ceil(n / 16384) | CRC-32 of every block of 4,096 packed bytes |
| ceil(n / 4) | the bases, four per byte, the first base in the low bits |
