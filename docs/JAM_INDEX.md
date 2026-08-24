# Jam Index format 1

A Jam Index is one logical local dataset containing a root manifest and one or
more independent, immutable parts. It is designed for query-centered trace
search across metagenomic assemblies without collection-wide nucleotide
positions.

Every part contains:

1. `screen.jam`: a pure `.jam` format-3 metagenome screening shard;
2. `part.bin`: a compact metagenome/contig directory, position-free contig
   signatures, and complete two-bit contig sequences with IUPAC exceptions.

The public format version is `1`. Ordinary `.jam` format 3 and JMA format 1
are independent formats and retain their existing behavior.

## Dataset layout

```text
metagenomes-index/
    manifest.json
    parts/
        part-000000/
            screen.jam
            part.bin
        part-000001/
            screen.jam
            part.bin
        ...
```

Part IDs and directories are consecutive. Append adds new directories after
the current final part; it never changes the bytes of an existing part. One
metagenome cannot span parts.

## Manifest

`manifest.json` has schema identity `jam-index-manifest-v1` and records:

- format version;
- the complete screen-selection policy;
- source-manifest SHA-256;
- total metagenomes, contigs, bases, and estimated signatures;
- ordered part descriptors.

Each part descriptor records its ID, relative directory, screen/data file
names, metagenome/contig/base counts, estimated signature count, logical
screen/signature/packed-sequence bytes, and SHA-256 for both files. Relative
paths cannot escape the index root. Totals must equal the sum of all part
descriptors.

The selected screen policy is:

```text
policy_id               = contig-minhash-v1
k                       = 21
hash_id                 = jamhash_u64_v1
zero_excluded           = true
contig minimum          = 16 hashes
contig maximum          = 256 hashes
bases per contig hash   = 1024
whole-metagenome hashes = 512
query window            = 256 bases
```

## Screening shard

`screen.jam` is an ordinary `.jam` format-3 file. Every metagenome is one
sample. Its stored hash set is the union of:

- one bottom-k signature for every contig, using the length-dependent budget;
- one fixed 512-hash bottom-k sketch across the complete metagenome.

A hash is stored once per sample even when selected by both rules. Hash zero
is excluded. The shard contains no contig IDs, nucleotide positions, or
sequence payload.

Every part can be opened and searched independently. Collection search merges
bounded top candidates from all parts and performs an exact recount before
retaining shared-hash detail.

## Part data envelope

`part.bin` starts with a 512-byte little-endian header identified by
`JAMIDX1P`. The header records version 1, collection counts, section offsets
and lengths, five section SHA-256 values, and a header checksum.

The sections are:

1. sequence payload;
2. sorted signature records;
3. fixed-width metagenome records;
4. fixed-width contig records;
5. UTF-8 string table.

All offsets and lengths are checked half-open ranges. Section overlap,
overflow, checksum mismatch, invalid directory totals, invalid contig ranges,
and malformed strings fail closed.

## Directories and signatures

A fixed 96-byte metagenome record binds:

- metagenome ID;
- contiguous range of contig IDs;
- total bases;
- screen hash count;
- source sequence SHA-256.

A fixed 128-byte contig record binds:

- global part-local contig ID and metagenome ID;
- UTF-8 contig name;
- base and signature counts;
- packed and ambiguity payload ranges;
- sequence checksum.

Each 16-byte signature record contains:

```text
jamhash u64, contig_id u32, flags u32
```

Flags distinguish contig-selected and whole-metagenome-selected hashes. The
table is sorted by hash for binary lookup. It contains contig IDs but no
nucleotide positions. A Stage-2 lookup maps only shared Stage-1 hashes to
bounded contig candidates.

Hash equality is never an alignment anchor. Stage 3 regenerates canonical
k=31 and k=21 positions from selected complete contigs and compares the exact
packed canonical k-mer before creating an anchor.

## Sequence payload

A/C/G/T use two bits per base. Every other supported IUPAC symbol is stored as
an exact `(position, symbol)` exception. Decoding a complete selected contig
checks its packed-plus-ambiguity checksum before returning bases.

The sequence payload is contig-addressable but not block-addressed. Jam Index
trace reads complete selected contigs by design. Noncandidate parts may be
searched through `screen.jam`, but their sequence payloads are not decoded.

## Local access and integrity

Both `screen.jam` and `part.bin` are memory mapped. Jam Index build, append,
and trace accept local paths only. No remote-object fallback, companion
positional index, or JMA lookup participates in this path.

The source catalog checksum binds each build or append operation. A stale,
missing, duplicate, nonlocal, or checksum-mismatched source fails before a
part is published. A new part is staged independently and the root manifest is
updated only after every new part passes construction.

## Evidence boundary

Contig signatures route work; they are not alignments. Base alignments are
reported only after selected sequence is decoded and exact packed seeds or
complete-contig base equality are verified. Separate contigs remain separate
sequence observations and are never reported as physically joined. Sequence
support does not prove that a plasmid, phage, or other element is autonomous.
