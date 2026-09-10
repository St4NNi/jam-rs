#!/usr/bin/env python3
"""Audit bounded JIDX evidence and a prefix/suffix reference-only routing hypothesis."""
import argparse
import collections
import hashlib
import json
import mmap
import struct
from pathlib import Path

TAG = 1 << 63
ROUTE = struct.Struct("<IQBB")


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").digest()


def save(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def canonical(key, k):
    reverse, value = 0, key
    for _ in range(k):
        reverse = (reverse << 2) | ((value & 3) ^ 3)
        value >>= 2
    return (reverse, True) if reverse < key else (key, False)


def words(sequence, k):
    key = valid = 0
    for end, base in enumerate(sequence.upper()):
        value = "ACGT".find(base)
        if value < 0:
            key = valid = 0
            continue
        key = ((key << 2) | value) & ((1 << (2 * k)) - 1)
        valid += 1
        if valid >= k:
            yield end + 1 - k, *canonical(key, k)


class Index:
    def __init__(self, path):
        self.path = Path(path)
        with self.path.open("rb") as stream:
            self.data = mmap.mmap(stream.fileno(), 0, access=mmap.ACCESS_READ)
        assert self.data[:8] == b"JIDX\0\0\0\0" and self.u16(8) == 4
        self.sections = {}
        for at in range(144, 360, 24):
            self.sections[self.u16(at)] = (self.u64(at + 8), self.u64(at + 16))
        hasher = hashlib.sha256()
        for at in range(512, len(self.data), 1024 * 1024):
            hasher.update(self.data[at:at + 1024 * 1024])
        assert hasher.digest() == self.data[112:144]
        self.documents = []
        for i in range(self.u32(24)):
            at = self.sections[2][0] + 80 * i
            self.documents.append((self.string(at), self.u32(at + 24)))
        self.contigs = []
        for i in range(self.u32(28)):
            at = self.sections[3][0] + 40 * i
            self.contigs.append((self.u32(at), self.string(at + 4), self.u64(at + 16)))

    def u16(self, at):
        return struct.unpack_from("<H", self.data, at)[0]

    def u32(self, at):
        return struct.unpack_from("<I", self.data, at)[0]

    def u64(self, at):
        return struct.unpack_from("<Q", self.data, at)[0]

    def string(self, at):
        start = self.sections[1][0] + self.u32(at)
        return self.data[start:start + self.u32(at + 4)].decode()

    def seed(self, ordinal):
        assert 0 <= ordinal < self.u64(32)
        return struct.unpack_from("<QQII", self.data, self.sections[4][0] + ordinal * 24)

    def varint(self, at):
        value = shift = 0
        while True:
            byte = self.data[at]
            at += 1
            value |= (byte & 127) << shift
            if byte < 128:
                return value, at
            shift += 7
            assert shift < 70

    def occurrences(self, ordinal):
        key, relative, members, zero = self.seed(ordinal)
        assert zero == 0
        at = self.sections[5][0] + relative
        first_at, position_bytes = at, 0
        k = 15 if key & TAG else 21
        previous_member = -1
        for _ in range(members):
            member, local, payload = struct.unpack_from("<IIQ", self.data, at)
            assert member > previous_member
            previous_member = member
            first_contig = self.documents[member][1]
            if local != (1 << 32) - 1:
                rows = [(local, payload >> 1, bool(payload & 1))]
                at += 16
            else:
                count, length = struct.unpack_from("<QQ", self.data, at + 16)
                position_bytes += length
                cursor = self.sections[6][0] + payload
                end = cursor + length
                rows, local, position = [], 0, 0
                for i in range(count):
                    packed, cursor = self.varint(cursor)
                    delta, cursor = self.varint(cursor)
                    new_local = local + (packed >> 1)
                    position = position + delta if i and new_local == local else delta
                    local = new_local
                    rows.append((local, position, bool(packed & 1)))
                assert cursor == end
                at += 32
            for local, position, reverse in rows:
                contig = first_contig + local
                assert self.contigs[contig][0] == member
                assert position + k <= self.contigs[contig][2]
                yield member, contig, position, reverse
        self.last_membership_bytes = at - first_at
        self.last_position_bytes = position_bytes


def audit(args):
    index = Index(args.index)
    counts = [collections.Counter(), collections.Counter()]
    starts = [[set() for _ in index.contigs] for _ in range(2)]
    routes = []
    for ordinal in range(index.u64(32)):
        key, _, members, _ = index.seed(ordinal)
        scheme = int(bool(key & TAG))
        rows = list(index.occurrences(ordinal))
        unique = set(rows)
        count = counts[scheme]
        count.update(keys=1, memberships=members, occurrences=len(rows),
                     dictionary_bytes=24, membership_bytes=index.last_membership_bytes,
                     position_bytes=index.last_position_bytes,
                     duplicate_physical_records=len(rows) - len(unique),
                     repeated_key_placements=max(0, len(unique) - 1),
                     keys_with_multiple_placements=int(len(unique) > 1))
        for _, contig, position, _ in unique:
            starts[scheme][contig].add(position)
        if scheme == 0 and args.routes:
            for offset in (0, 6):
                short, reverse = canonical((key >> (2 * (6 - offset))) & ((1 << 30) - 1), 15)
                routes.append((short, ordinal, offset, reverse))
    assert sum(c["occurrences"] for c in counts) == index.u64(40)
    gap_rows = []
    for scheme, k in ((0, 21), (1, 15)):
        longest = uncovered = empty = short_contigs = 0
        for (_, _, length), positions in zip(index.contigs, starts[scheme]):
            ordered = sorted(positions)
            empty += not ordered
            start_domain = max(0, length - k + 1)
            short_contigs += length < k
            longest = max(longest, max((b - a - 1 for a, b in zip([-1] + ordered, ordered + [start_domain])), default=start_domain))
            end = 0
            for position in ordered:
                uncovered = max(uncovered, position - end)
                end = max(end, position + k)
            uncovered = max(uncovered, length - end)
        gap_rows.append({"k": k, "longest_interval_without_selected_start": longest,
                         "longest_base_interval_uncovered_by_seed": uncovered,
                         "contigs_without_seed": empty, "contigs_shorter_than_k": short_contigs,
                         "start_domain": "0 <= start <= length-k; ambiguous windows included as unseeded"})
    overlap = sum(len(a & b) for a, b in zip(*starts))
    union = sum(len(a | b) for a, b in zip(*starts))
    result = {"index": str(args.index), "bytes": len(index.data), "sha256": digest(args.index).hex(),
              "window": index.u16(360), "source_bases": sum(c[2] for c in index.contigs),
              "contigs": len(index.contigs), "schemes_k21_k15": counts, "gaps": gap_rows,
              "unique_locations_ignoring_scheme_and_orientation": union,
              "locations_in_both_schemes": overlap,
              "section_bytes": {str(k): v[1] for k, v in index.sections.items()},
              "header_padding_bytes": len(index.data) - sum(v[1] for v in index.sections.values())}
    if args.routes:
        routes.sort()
        body = hashlib.sha256()
        with args.routes.open("xb") as stream:
            stream.write(b"JEBRTE01" + digest(args.index) + bytes(32) + struct.pack("<Q", len(routes)))
            for row in routes:
                encoded = ROUTE.pack(*row)
                stream.write(encoded)
                body.update(encoded)
            stream.seek(40)
            stream.write(body.digest())
        fanouts = collections.Counter(row[0] for row in routes)
        result["shared_anchor_hypothesis"] = {
            "context_offsets_in_canonical_k21": [0, 6], "context_length": 15,
            "references": len(routes), "complete_flat_route_bytes": args.routes.stat().st_size,
            "distinct_short_keys": len(fanouts), "maximum_reference_fanout": max(fanouts.values(), default=0),
            "short_keys_with_multiple_references": sum(v > 1 for v in fanouts.values()),
            "sha256": digest(args.routes).hex(), "independent_position_bytes": 0,
            "semantics": "sorted exact short key plus anchor ordinal, offset and orientation; no full-key absence filter",
            "not_a_selected_product_route": True,
            "unpriced": ["production block authentication directory", "production direct alignment integration"],
        }
    save(args.output, result)


def fasta(path):
    name, parts = None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(parts)
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line)
    if name is not None:
        yield name, "".join(parts)


def route(args):
    index = Index(args.index)
    with args.routes.open("rb") as stream:
        data = mmap.mmap(stream.fileno(), 0, access=mmap.ACCESS_READ)
    assert data[:8] == b"JEBRTE01" and data[8:40] == digest(args.index)
    assert hashlib.sha256(data[80:]).digest() == data[40:72]
    count = struct.unpack_from("<Q", data, 72)[0]
    assert len(data) == 80 + count * ROUTE.size
    with args.output.open("x") as output:
        for name, sequence in fasta(args.query):
            query = sequence + sequence[:14] if args.circular else sequence
            matches, references, decoded, comparisons = [], 0, 0, 0
            query_keys = collections.defaultdict(list)
            for position, key, reverse in words(query, 15):
                if position < len(sequence):
                    query_keys[key].append((position, reverse))
            for key, query_positions in sorted(query_keys.items()):
                lo, hi = 0, count
                while lo < hi:
                    middle = (lo + hi) // 2
                    comparisons += 1
                    if ROUTE.unpack_from(data, 80 + middle * ROUTE.size)[0] < key:
                        lo = middle + 1
                    else:
                        hi = middle
                while lo < count:
                    short, ordinal, offset, sub_reverse = ROUTE.unpack_from(data, 80 + lo * ROUTE.size)
                    if short != key:
                        break
                    references += 1
                    for member, contig, position, anchor_reverse in index.occurrences(ordinal):
                        decoded += 1
                        target_position = position + (6 - offset if anchor_reverse else offset)
                        matches.append({"k": 15, "sample": index.documents[member][0],
                                        "contig": index.contigs[contig][1], "target_position": target_position,
                                        "target_canonical": bool(anchor_reverse) ^ bool(sub_reverse),
                                        "query_positions": query_positions})
                    lo += 1
            unique = {(h["sample"], h["contig"], h["target_position"], h["target_canonical"]) for h in matches}
            output.write(json.dumps({"query_id": name, "short_dictionary_probes": len(query_keys),
                                     "dictionary_comparisons": comparisons, "anchor_references_followed": references,
                                     "unique_derived_locations": len(unique), "duplicate_derived_records": len(matches) - len(unique),
                                     "anchor_positions_decoded": decoded, "placements": matches}) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["audit", "route"])
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--routes", type=Path)
    parser.add_argument("--query", type=Path)
    parser.add_argument("--circular", action="store_true")
    args = parser.parse_args()
    if args.mode == "audit":
        audit(args)
    else:
        assert args.routes and args.query
        route(args)


if __name__ == "__main__":
    main()
