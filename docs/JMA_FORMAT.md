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
header checksum, and up to two seed-level identities. Header checksum coverage
includes reserved bytes with the checksum field zeroed.

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
