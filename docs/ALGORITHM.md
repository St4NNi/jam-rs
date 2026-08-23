# Trace algorithm

## Identifier and scope

The positional trace algorithm is identified as:

```text
jam-seed-chain-align-v1
```

Algorithm version `1` combines candidate retrieval from an existing `.jam`
index with exact positional seed lookup, bounded anchor chaining, fixed-band
local affine-gap alignment, and nonredundant plasmid-coordinate coverage.
It searches one plasmid query against many metagenome assemblies. It is not a
general sequence aligner, read mapper, plasmid assembler, or plasmid-origin
classifier.

The implementation is distinct from BLAST and minimap2. It does not use their
seed-selection, chaining, extension, indexing, or reporting algorithms, and no
novelty or general performance claim follows from that distinction.

## Coordinates

All reported intervals are zero-based and half-open: `[start,end)`. Query
coordinates use the forward plasmid sequence. Target coordinates use the
forward stored contig even when the accepted alignment is reverse-strand.

A circular query span may cross the chosen sequence origin once. Such a span
is represented by two ordered, non-wrapping query intervals and
`origin_crossing=true`. Unwrapped query coordinates used internally by
chaining are not substituted for the normalized output intervals.

## Pipeline

```text
read exactly one plasmid FASTA record
build its candidate-query FracMinHash sketch
query the metagenome .jam index
apply independent shared-hash and containment filters
rank deterministically and retain the candidate limit

for each retained metagenome:
    resolve its JMA archive or raw assembly resource
    if JMA is available:
        extract retained k=31 query seeds
        look up exact packed-kmer occurrences
        discard over-repetitive seed groups
        build bounded monotone chains per contig and strand
        if no primary chain exists and rescue is enabled:
            repeat with retained k=21 seeds
        retrieve target windows for retained chains
        run fixed-band local affine-gap alignment on the hinted strand
    if the catalog intentionally supplies only a raw assembly:
        tile raw contigs into bounded overlapping windows
        run alignment on both strands
    rank and deduplicate alignments
    select alignments that contribute novel plasmid support
    union supported plasmid intervals and compute complement gaps
    emit one candidate result or structured failure
```

## Candidate sketch

Candidate retrieval uses every canonical k-mer window accepted by the FASTA
parser. It does not use minimizers. Each packed canonical k-mer is hashed with
`jamhash_u64_v1`; hash zero is excluded. The `.jam` database threshold,
entropy setting, and optional catalog-specific bias table determine retained
hashes. Duplicate hashes are removed before querying the memory-mapped `.jam`
version-3 index.

The candidate filters have independent denominators:

```text
plasmid_containment    = shared_hashes / plasmid_hashes
metagenome_containment = shared_hashes / metagenome_hashes
```

Candidate retrieval is sketch evidence only.

## Positional seeds

The positional seed rule also examines every canonical k-mer window; it is not
a minimizer scheme. For k-mer size `k` and FracMinHash scale `s`, a seed is
retained exactly when:

```text
h = jamhash_u64_v1(canonical_packed_kmer)
h != 0
h < floor(u64::MAX / s)
```

Smaller scale values retain denser nested sets. Ambiguous windows are skipped.
The primary seed size is k=31. The optional rescue size is k=21; no smaller
k-mer size is accepted.

The JMA lookup key is the complete tuple:

```text
(k, hash, canonical_packed_kmer)
```

The stored occurrence contains contig ID, forward k-mer start, and the target
canonicalization orientation. Comparing the packed k-mer as well as the hash
ensures that a hash collision between unequal k-mers cannot become an anchor.

## Anchors and repetitive seeds

An anchor joins a query seed position to one exact JMA occurrence. Relative
strand is the exclusive-or of the query and target canonicalization flags.

If an exact seed has more occurrences than the profile's
`max_occurrences`, the complete occurrence group is discarded. Groups exactly
at the limit are retained. Anchors are sorted by contig, strand, query
position, target position, hash, packed k-mer, k, and orientation flags;
duplicates are removed and the profile anchor limit is applied
deterministically.

## Chaining

Anchors are grouped by contig, strand, and k. k=31 and k=21 anchors are not
mixed in one chain. Query anchors are duplicated at one plasmid-length offset
to permit one circular-origin crossing; a source anchor cannot be used twice.

For anchor `i`, the dynamic program initializes `score[i]=100` and considers
at most `P` earlier anchors, where `P=min(max_anchors,256)`. A predecessor `j`
must have a positive query gap and a positive orientation-correct target gap,
both within their configured limits. The recurrence is:

```text
extension = 100 - gap_penalty * abs(query_gap - target_gap)
candidate = score[j] + extension
score[i]  = max(score[i], candidate)
```

All shipped profiles use `gap_penalty=1`; this is not an affine chain gap
model. Ties prefer more anchors and then deterministic coordinate/hash keys.
Chain extraction requires the profile's minimum anchor count and is bounded by
the profile chain limit.

For `m` anchors in one group, one pass costs
`O(m log m + mP)`. Up to `C` retained-chain passes give a worst-case bound of
`O(C(m log m + mP))` for that group, plus initial sorting. The anchor and chain
limits are correctness-affecting resource bounds: a repeat-rich candidate can
lose later evidence when a limit is reached.

## Alignment

The alignment kernel is local alignment with three affine-gap states. It is
neither global, semiglobal, nor overlap alignment. Non-positive scores restart
at zero, and traceback begins at the highest positive cell and stops at a
restart cell.

For query base `q` and target base `t`:

```text
M(i,j) = best(M,I,D at i-1,j-1) + (match if q==t else mismatch)
I(i,j) = max(I(i,j-1) + extend,
             M(i,j-1) + open + extend,
             D(i,j-1) + open + extend)
D(i,j) = max(D(i-1,j) + extend,
             M(i-1,j) + open + extend,
             I(i-1,j) + open + extend)
```

Only positive state values are retained. The first base of a gap therefore
costs `gap_open + gap_extend`. Shipped scores are match `+2`, mismatch `-3`,
gap open `-5`, and gap extend `-1`.

The dynamic-programming matrix is restricted to a fixed half-band around the
chain-derived diagonal. The band is clipped to target bounds, never widens,
and has no X-drop. A traceback that needs a cell outside the band is an error;
there is no automatic wider retry in algorithm version 1.

Unaligned query and target prefixes/suffixes are represented by the reported
local intervals, not by soft-clipping operations in the generated CIGAR. The
CIGAR uses `=`, `X`, `I`, and `D`. Reverse-strand targets are reverse
complemented for alignment and mapped back to forward target coordinates.

## Coverage and gaps

`=` and `X` consume query and target and count as supported plasmid bases.
`D` consumes query and contributes to the aligned span but not supported
coverage. `I` consumes target only. Consequently a coverage gap is the
complement of supported plasmid coordinates and is not the same as an
alignment deletion.

Alignments are ranked by supported bases, alignment score, chain score,
matches, fewer edits, contig ID, target interval, and a deterministic final
tie-breaker. In that order, an alignment is primary if it contributes at least
one previously unsupported plasmid base. Fully redundant alignments remain
visible as secondary evidence but never increase `supported_bases`.

Primary supported intervals are validated, sorted, and unioned. Coverage is
therefore bounded by plasmid length, and overlapping or duplicate contigs do
not count a plasmid base twice. Gaps are the circular complement of this union.
Separate supporting contigs are not assumed to be physically connected.

## Default profiles

| Parameter | fast | balanced | sensitive |
| --- | ---: | ---: | ---: |
| primary k / scale | 31 / 500 | 31 / 200 | 31 / 100 |
| rescue k / scale | disabled | 21 / 500 | 21 / 200 |
| max occurrences per seed | 64 | 128 | 256 |
| max anchors per candidate | 25,000 | 100,000 | 250,000 |
| max chains per candidate | 4 | 8 | 16 |
| minimum chain anchors | 3 | 3 | 2 |
| maximum alignment window | 250,000 | 1,000,000 | 4,000,000 |
| alignment half-band | 64 | 128 | 256 |
| candidate concurrency target | 4 | 4 | 2 |

The profile names are execution settings, not calibrated biological
sensitivity levels.

## Known failure boundaries

- A true trace can be missed when retained seeds are absent, suppressed as
  repetitive, or truncated by an anchor/chain/window limit.
- Rescue seeds are used only when primary chaining yields no chain; primary
  and rescue evidence are not merged.
- Rearrangements require separate monotone chains and may be truncated by the
  chain limit.
- A true alignment outside the fixed band can fail or be missed because the
  band does not widen.
- Shared mobile elements, chromosomal integrations, repeats, and contamination
  can produce genuine sequence evidence without an autonomous plasmid.
- Local alignment omits unmatched prefixes and suffixes; reported coverage
  must be interpreted from query intervals rather than CIGAR alone.
- Resource and memory limits can turn a pathological candidate into a
  structured failure rather than an exhaustive search.
- A JMA seed level sparser than the requested profile is a structured
  `seed_level_mismatch` failure. It does not fall back to unanchored raw
  alignment, because that would change the requested algorithm semantics.

An accepted alignment is sequence evidence. It does not prove plasmid
autonomy, physical linkage between contigs, or biological presence without
independent confirmation.
