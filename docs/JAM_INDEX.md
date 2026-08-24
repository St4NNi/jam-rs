# Jam Index format 1

A Jam Index is one logical local dataset containing a root manifest and
independent, immutable parts. It routes one query into selected metagenomes and
contigs without storing nucleotide positions or nucleotide sequence.

Every part contains:

1. `screen.jam`, a pure `.jam` format-3 metagenome screening shard;
2. `part.bin`, containing hash-to-contig postings, metagenome and contig
   metadata, and checksummed external assembly references.

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

Twenty parts are the default large-collection target; 20 to 30 parts are the
intended range. Parts are balanced by source bases and estimated signature
count. One metagenome always remains wholly within one part. Each part builds
and searches independently, and append creates only new part directories.

## Manifest

`manifest.json` has schema identity `jam-index-manifest-v1` and records:

- format version;
- the complete screen-selection policy;
- source-manifest SHA-256;
- total metagenomes, contigs, bases, and estimated signatures;
- ordered part descriptors.

Each part descriptor records its ID, relative directory, screen/data file
names, metagenome/contig/base counts, estimated signature count, logical
screen/posting/source-reference bytes, and SHA-256 for both files. Relative
paths cannot escape the index root. Totals must equal the sum of all parts.

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

Every part can be searched independently. Collection search merges bounded
top candidates from all parts and performs an exact recount before retaining
shared-hash detail for globally selected metagenomes.

## Part data envelope

`part.bin` starts with a 512-byte little-endian header identified by
`JAMIDX1P`. The header records format version 1, collection counts, section
offsets and lengths, five section SHA-256 values, and a header checksum.

The sections are:

1. fixed-width external source-reference records;
2. posting headers followed by delta-coded contig-ID payloads;
3. fixed-width metagenome records;
4. fixed-width contig records;
5. a UTF-8 string table.

All offsets and lengths are checked half-open ranges. Section overlap,
overflow, checksum mismatch, invalid directory totals, invalid contig ranges,
and malformed strings fail closed.

## External references and contigs

Each source-reference record binds one metagenome to:

- assembly path, size, and catalog SHA-256;
- access mode;
- optional FAI and GZI paths and checksums.

Each contig record binds its part-local contig and metagenome IDs, name, base
count, source ordinal, FAI byte geometry when available, and sequence SHA-256.
It contains no packed bases.

Supported access modes are:

- plain FASTA with FAI: exact selected-contig byte ranges;
- BGZF FASTA with FAI/GZI: selected BGZF blocks and contig ranges;
- normal gzip: one sequential pass for a selected candidate assembly.

Normal gzip is a candidate-only fallback. No assembly file is opened for a
metagenome that did not survive collection screening.

## Aligned contig postings

The posting table does not repeat hashes. Its headers have the same logical
bucket/sample/hash order as the entries in `screen.jam`. A Stage-1 hit carries
that stable entry ordinal directly into Stage 2.

Each 16-byte posting header stores only:

```text
payload offset u64, contig count u32, payload length u32
```

The payload contains sorted, delta-coded part-local contig IDs. Posting count
must equal the number of logical `.jam` entries, and malformed, unsorted, or
out-of-range postings fail closed. Stage 2 reads postings only for shared hashes
belonging to globally selected metagenomes.

Hash equality is never an alignment anchor. Stage 3 reads complete selected
contigs from their assembly files, generates dense canonical k=31 and k=21
positions in memory, and verifies the packed canonical k-mer before creating
an anchor.

## Local access and integrity

Both `screen.jam` and `part.bin` are memory mapped. External assemblies remain
ordinary user-managed files. Query-time validation checks source path, size,
access mode, sidecar identity, and every selected contig checksum. A stale or
corrupt selected source fails closed.

Build and append stage each new part independently. The root manifest is
published only after all new parts pass construction. Append never rewrites an
existing part.

## Evidence boundary

Contig signatures route work; they are not alignments. Base alignments are
reported only after exact packed seed verification or complete-contig base
equality. Separate contigs remain separate sequence observations and are never
reported as physically joined. Sequence support does not prove that a plasmid,
phage, or other element is autonomous.
