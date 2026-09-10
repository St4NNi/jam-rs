#!/usr/bin/env python3
"""Reproduce the six frozen sparse 50 to 250 kb development queries on BCF.

The rotated LexicMap file is a preserved, unused preparation artifact.
Retained comparisons use queries.jam.fasta with lexicmap_circular.py instead.
"""

import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path


TARGET_SHA256 = "0282d1afef69eafd4a0d42875bf50f85a892672ec128923eeadbdc95f36ecd9b"
SEED_SHA256 = "1169a7ce2bbe0be3b22af8e3a49c358da05ecbfec0ad9696230f2c78e03b26c2"
SAMTOOLS_SHA256 = "a75a28d2551e2073257748c9057c18f7d250b092303a40cebf40a5a45b6a002a"
BASES = b"ACGT"
SPECS = (
    ("long_050k_exact", 50_000, 400, "exact", "linear"),
    ("long_075k_exact", 75_000, 800, "exact", "linear"),
    ("long_100k_substitution", 100_000, 800, "substitution_dispersed", "linear"),
    ("long_150k_indel", 150_000, 800, "balanced_indel_5", "linear"),
    ("long_200k_exact", 200_000, 1_600, "exact", "linear"),
    ("long_250k_origin", 250_000, 400, "exact", "circular"),
)


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def dna(key: bytes, label: str, length: int) -> bytes:
    result = bytearray()
    counter = 0
    while len(result) < length:
        block = hashlib.sha256(key + b"\0" + label.encode() + counter.to_bytes(8, "little")).digest()
        result.extend(BASES[value & 3] for value in block)
        counter += 1
    return bytes(result[:length])


def reverse_complement(sequence: bytes) -> bytes:
    return sequence.translate(bytes.maketrans(b"ACGT", b"TGCA"))[::-1]


def extract(samtools: Path, target: dict, contig: dict, start: int, width: int) -> bytes:
    region = f"{contig['name']}:{start + 1}-{start + width}"
    run = subprocess.run([str(samtools), "faidx", target["bgzf"], region],
                         check=False, capture_output=True)
    if run.returncode:
        raise ValueError(f"target extraction failed for {target['source_id']}")
    sequence = b"".join(line.strip() for line in run.stdout.splitlines()
                        if not line.startswith(b">" )).upper().replace(b"U", b"T")
    if len(sequence) != width or any(base not in BASES for base in sequence):
        raise ValueError("target interval is not unambiguous ACGT")
    return sequence


def write_fasta(path: Path, records: list[tuple[str, bytes, str]]) -> None:
    with path.open("xb") as stream:
        for name, sequence, topology in records:
            stream.write(f">{name} topology={topology}\n".encode())
            for start in range(0, len(sequence), 80):
                stream.write(sequence[start:start + 80] + b"\n")
        stream.flush()
        os.fsync(stream.fileno())


def query_ranges(length: int, topology: str, width: int) -> list[list[list[int]]]:
    if topology == "circular":
        rest = []
        gap = (length - 7 * width) // 8
        for index in range(7):
            start = gap + index * (gap + width)
            rest.append([[start, start + width]])
        return [[[length - width // 2, length], [0, width // 2]], *rest]
    gap = (length - 8 * width) // 9
    return [[[gap + index * (gap + width), gap + index * (gap + width) + width]]
            for index in range(8)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--seed", type=Path, required=True)
    parser.add_argument("--samtools", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise SystemExit("new output required")
    for path, expected in ((args.targets, TARGET_SHA256), (args.seed, SEED_SHA256),
                           (args.samtools, SAMTOOLS_SHA256)):
        if not path.is_file() or path.is_symlink() or digest(path) != expected:
            raise SystemExit(f"input identity differs: {path}")
    key = args.seed.read_bytes()
    frozen = json.loads(args.targets.read_text(encoding="utf-8"))
    targets = frozen["targets"]
    by_member = {}
    for target in targets:
        eligible = [contig for contig in target["selected_contigs"]
                    if int(contig["length"]) >= 1_600]
        ranked = sorted(eligible, key=lambda contig: hashlib.sha256(
            key + b"\0long-query-contig\0" + target["source_id"].encode()
            + b"\0" + contig["name"].encode()).digest())
        by_member[target["source_id"]] = (target, ranked)
    members = sorted(by_member, key=lambda name: hashlib.sha256(
        key + b"\0long-query-member\0" + name.encode()).digest())
    used = set()
    records = []
    lex_records = []
    truth_rows = []
    for query_index, (name, length, width, transform, topology) in enumerate(SPECS):
        query = bytearray(dna(key, f"long-background:{name}", length))
        ranges = query_ranges(length, topology, width)
        components = []
        selected_members = [members[(query_index * 3 + index) % len(members)] for index in range(8)]
        for component_index, (member, segments) in enumerate(zip(selected_members, ranges, strict=True)):
            target, contigs = by_member[member]
            available = [contig for contig in contigs if (member, contig["name"]) not in used
                         and int(contig["length"]) >= width]
            if not available:
                raise ValueError("insufficient independent eligible contigs")
            contig = available[0]
            used.add((member, contig["name"]))
            room = int(contig["length"]) - width + 1
            initial = int.from_bytes(hashlib.sha256(
                key + b"\0long-query-start\0" + name.encode()
                + component_index.to_bytes(4, "little")).digest()[:8], "little") % room
            target_sequence = None
            for attempt in range(min(room, 64)):
                start = (initial + attempt * max(1, room // 64)) % room
                try:
                    target_sequence = extract(args.samtools, target, contig, start, width)
                    break
                except ValueError as error:
                    if "unambiguous" not in str(error):
                        raise
            if target_sequence is None:
                raise ValueError("no unambiguous target interval found")
            strand = "reverse" if component_index % 2 else "forward"
            source_query = reverse_complement(target_sequence) if strand == "reverse" else target_sequence
            if transform == "substitution_dispersed":
                changed = bytearray(source_query)
                for position in range(15, len(changed), 31):
                    changed[position] = BASES[(BASES.index(changed[position]) + 1) % 4]
                source_query = bytes(changed)
            elif transform == "balanced_indel_5":
                source_query = (source_query[:200] + source_query[205:500]
                                + dna(key, f"long-insertion:{name}:{component_index}", 5)
                                + source_query[500:])
            offset = 0
            for left, right in segments:
                count = right - left
                query[left:right] = source_query[offset:offset + count]
                offset += count
            components.append({
                "source_id": target["source_id"],
                "contig": contig["name"],
                "strand": strand,
                "target_interval": [start, start + width],
                "query_segments": segments,
                "transform": transform,
                "fragment_sha256": hashlib.sha256(target_sequence).hexdigest(),
            })
        sequence = bytes(query)
        records.append((name, sequence, topology))
        rotation = length // 2 if topology == "circular" else 0
        lex_sequence = sequence[rotation:] + sequence[:rotation]
        lex_records.append((name, lex_sequence, "linear"))
        truth_rows.append({
            "query_id": name,
            "split": "development",
            "query_bases": length,
            "retained_target_bases": 8 * width,
            "retained_fraction": 8 * width / length,
            "transform": transform,
            "topology": topology,
            "expected_source_components": 8,
            "expected_distinct_source_pairs": 8,
            "independent_target_fragments": True,
            "lexicmap_rotation_left": rotation,
            "components": components,
        })
    args.output.mkdir(mode=0o700, parents=True)
    write_fasta(args.output / "queries.jam.fasta", records)
    write_fasta(args.output / "queries.lexicmap.fasta", lex_records)
    with (args.output / "truth.jsonl").open("x", encoding="utf-8") as stream:
        for row in truth_rows:
            stream.write(json.dumps(row, sort_keys=True) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    files = {}
    for filename in ("queries.jam.fasta", "queries.lexicmap.fasta", "truth.jsonl"):
        path = args.output / filename
        files[filename] = {"bytes": path.stat().st_size, "sha256": digest(path)}
    manifest = {
        "format": "jam-shared-core-long-development-v1",
        "status": "frozen_before_candidate_reader_settings_or_timing",
        "target_sha256": TARGET_SHA256,
        "seed_sha256": SEED_SHA256,
        "samtools_sha256": SAMTOOLS_SHA256,
        "queries": len(SPECS),
        "length_range_inclusive": [50_000, 250_000],
        "independent_fragment_constructions": len(SPECS),
        "fragment_reuse_across_queries": False,
        "jam_topology": "per FASTA header token",
        "lexicmap_circular_representation": "single deterministic rotation, same sequence length, no duplicated bases",
        "final_split": "sealed and unused",
        "files": files,
        "rows": [{key: value for key, value in row.items() if key != "components"}
                 for row in truth_rows],
    }
    with (args.output / "manifest.json").open("x", encoding="utf-8") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())


if __name__ == "__main__":
    main()
