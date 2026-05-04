# jam architecture

jam (just another minimizer) finds small sequences (genes, plasmids, contigs; 1 kb to a few Mb)
in very large sets of assemblies and aligns them, base by base, with CIGAR strings. One index
answers three questions, one per stage:

1. **Sketch**: which entries share enough seeds with the query to hold it, and roughly how much
   of it?
2. **Locate**: where in those entries do the shared seeds lie?
3. **Align**: how exactly does the query align there, and does the entry qualify?

A search can stop after any stage. Each stage reads more data than the one before, so each one
only works on what the previous stage kept.

## 1. Terms

Every module, type, option, comment and document uses these terms with exactly this meaning.

| Term | Meaning |
| --- | --- |
| entry | One unit of the index: one sequence (a plasmid) or one set of sequences (an assembly). Results are reported per entry. |
| contig | One sequence (one FASTA record) of an entry. |
| query | The sequence searched for. It is linear or circular. |
| k-mer | A substring of length k; jam always uses the canonical form (the smaller of the k-mer and its reverse complement). |
| seed | A k-mer of an entry that the index stores: the smallest k-mer (by hash) of each minimizer window, and the smallest of every clean run too short for a full window. |
| minimizer window | w consecutive k-mers of a contig; each holds at least one seed. |
| query k-mer | Any k-mer of the query. The query side uses all its k-mers; only the entry side selects seeds. |
| seed hit | A query k-mer that equals a seed of an entry. |
| sketch | The seeds of an entry. Stage 1 compares the query k-mers with the sketches of all entries. |
| color set | The set of entries that hold a seed. |
| containment estimate | The share of the query that an entry is estimated to hold, from its seed hits alone. |
| coverage bound | The highest query coverage an entry can still reach, given where its seed hits lie on the query. It is wrong at most with the chance set by the loss budget. |
| loss budget | The largest accepted chance of missing a qualifying entry (for example 0.01%). Every shortcut of a search stays within it. |
| candidate entry | An entry that stage 1 keeps. |
| bin | A block of 2^b bases (1,024 by default) of an entry's sequence. Positions in the positional index are bins. |
| positional index | For each seed, the bins (with strand) where it occurs in each entry. |
| posting | One stored bin and strand of one seed occurrence. |
| fingerprint | The few bits of a seed's hash that the index stores to tell seeds apart; two seeds can share one (a false seed hit). |
| candidate bin | A bin with seed hits of a candidate entry. |
| alignment window | The stretch of a contig that stage 3 reads around a candidate bin. |
| packed sequence | The contigs of an entry, 2 bits per base, with a table of contig names, lengths and runs of N. |
| anchor | An exact match of length a (15 by default) between the query and an alignment window. |
| chain | Anchors in order along nearly the same diagonal. |
| alignment | A gapped, base-level alignment of a query interval to a contig interval, with score, CIGAR and counts. |
| fragment | An alignment placed in query and contig coordinates; a fragment over the end of a circular query has two query segments. |
| mosaic | The fragments of one entry chosen to cover the query without counting a query base twice. |
| query coverage | The share of query bases covered by an entry's mosaic. |
| identity | Matches divided by alignment columns, over the fragments of the mosaic together. |
| decision | The rule that reports an entry: query coverage, identity and the contig rule (each counted contig is aligned at least halfway at the identity limit). |
| shard | A part of an index (a set of entries) that one process opens. An index has one or more shards. |
| coordinator | The server that sends a query to every shard and merges their results. |
| screen | The reverse direction: the small sequences form the index and large assemblies stream past it once. Results use the same decision. |
| manifest | The JSON file of a shard: parameters, entries, files and their SHA-256 digests. |
| block cache | Fixed-size blocks of packed sequences kept in RAM (and optionally on local disk) between reads. |

Words not used: sample (say entry), catalog (say index), locator (say positional index), trace,
a bare "window" (say minimizer window or alignment window), a bare "hit" (say seed hit), and
"sketch" in the MinHash sense.

All hashing uses jamhash (`jamhash_u64` of the 2-bit canonical k-mer, k at most 32).

## 2. One search, step by step

```
query
  -> query k-mers, hashed
  -> stage 1, sketch     (RAM)      seed hits per entry -> containment estimate, coverage bound
                                    -> candidate entries
  -> stage 2, locate     (RAM/disk) postings of the candidate entries -> candidate bins
  -> stage 3, align      (S3/disk)  alignment windows from packed sequences -> anchors -> chains
                                    -> alignments (CIGAR) -> fragments -> mosaic -> decision
  -> reported entries, streamed as each one is decided
```

Each entry moves through three states and is streamed to the caller as it changes: possible
(after stage 1, with its containment estimate and coverage bound), located (after stage 2) and
confirmed (after stage 3, with alignments). Entries whose coverage bound falls below the
requested query coverage are dropped at once. Rare seeds are looked up first, so the bounds
tighten early.

## 3. Index on disk

A shard is one directory. All files are laid out for memory mapping: fixed-width little-endian
fields, aligned arrays, offsets instead of pointers, a version number and a checksum per section.

| File | Stage | Holds | Read |
| --- | --- | --- | --- |
| `manifest.json` | all | parameters (k, w, bin size, bucket bits, fingerprint bits), entries (name, contigs, bases), file sizes and SHA-256 digests | at open |
| `sketch` | 1 | a directory of hash buckets; per bucket the sorted fingerprints of its seeds and the color set id of each; then every distinct color set once, entries as sorted ids | memory-mapped and kept in RAM |
| `positions` | 2 | per entry a section with its seeds' postings, grouped by seed hash, so that one candidate entry needs few reads | memory-mapped from local disk, or byte ranges from S3 |
| `packed/` | 3 | one packed sequence file per entry | byte ranges from local disk or S3 |

The interfaces the search uses:

- `Sketch::color_set(hash) -> Option<ColorSetId>` and `ColorSets::entries(id)`: stage 1.
- `Positions::postings(entry, hash)`: stage 2, the bins and strands of one seed in one entry.
- `PackedSequence::read(entry, contig, range)`: stage 3, bases of one contig.

Joining shards keeps every entry's positions and packed sequence as they are and rebuilds only
the sketch, whose color sets gain the entries of all joined shards.

## 4. Source layout

A module uses only modules listed above it.

| Module | Job |
| --- | --- |
| `seeds` | canonical k-mers, jamhash, minimizer windows, the low-complexity filter |
| `loss` | the chance that a present query stretch shows no seed hit; the loss budget |
| `storage` | reading bytes by offset: local files, S3 byte ranges, the block cache |
| `index` | the on-disk index: manifest, packed sequences, positional index, building and joining shards |
| `align` | alignment types, CIGAR, the banded alignment kernel, anchors and chains |
| `decide` | fragments, mosaic, query coverage, identity and the decision |
| `search` | the three stages for one query against one shard, with memory limits and streaming |
| `screen` | the reverse direction over streamed assemblies |
| `service` | `serve` and remote search over TCP: protocol, server, client, coordinator |
| `cli` | the commands: `index`, `search`, `screen`, `serve`, `stats` |

## 5. Code rules

- Use the terms of section 1 exactly. No abbreviations beyond k, w, CIGAR, S3, TCP and RAM.
- One job per module. Files stay under about 800 lines and functions under about 80; split by
  step, not by helper type.
- Every public item has a doc comment of one or two plain sentences.
- One error type per module; bad input never panics; `unwrap` only in tests.
- No dead code, no experiment flags, no counters that nothing reads.
- Tests sit next to the code they test. Every file format has a round-trip test; the search has
  end-to-end tests on small generated entries.
- Moves and renames are separate commits from changes of behaviour or format.

## 6. Plan

jam 0.11 is built new in this layout, from the first commit with jamhash and the formats of
section 3. Code of jam 0.10 is reused only where it fits (the alignment kernel, the S3 reader,
the packed sequence format), and is cleaned while it moves. Results are compared with jam 0.10 on
test sets: the reported entries must agree, except where the different seeds explain it.

1. **Foundations, in parallel**: `seeds` and `index`; `storage`; `align`; `decide` and `loss`.
2. **Search**: the three stages, the coverage bound with a loss budget, rare seeds first, entry
   states streamed as they change.
3. **Screen, service and commands**; the coordinator spreads stage 3 work over idle shards
   (every shard can read every packed sequence from S3).
4. **Speed**: faster than LexicMap on every query length, on local disk and on network storage.
