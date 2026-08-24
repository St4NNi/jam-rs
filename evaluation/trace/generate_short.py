#!/usr/bin/env python3
"""Generate an anonymous deterministic short-trace benchmark collection."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
from pathlib import Path


LENGTHS = (80, 100, 160, 250, 500, 1000)
IDENTITIES = (100, 99, 97, 95, 90)


def dna(rng: random.Random, length: int) -> str:
    return "".join(rng.choices("ACGT", k=length))


def revcomp(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def mutate(sequence: str, identity: int, seed: int) -> str:
    if identity == 100:
        return sequence
    rng = random.Random(seed)
    result = list(sequence)
    changes = max(1, round(len(sequence) * (100 - identity) / 100))
    for position in rng.sample(range(len(result)), min(changes, len(result))):
        choices = [base for base in "ACGT" if base != result[position]]
        result[position] = choices[rng.randrange(len(choices))]
    return "".join(result)


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="ascii") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset : offset + 80] + "\n")


def checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"refusing to reuse output directory: {args.output}")
    assemblies = args.output / "assemblies"
    assemblies.mkdir(parents=True)

    rng = random.Random(0x4A414D53484F5254)
    query = list(dna(rng, 10_000))
    repeated = "".join(query[1800:2048])
    query[6200:6448] = repeated
    query = "".join(query)
    background = dna(rng, 5_000)
    rows: list[dict[str, str]] = []
    sources: list[tuple[str, Path, str]] = []

    def add_case(
        case_id: str,
        records: list[tuple[str, str]],
        intervals: list[list[int]],
        length: int,
        identity: int,
        orientation: str,
        truth_class: str,
        expected: bool,
        standalone: bool = False,
    ) -> None:
        filename = f"{case_id}.fasta"
        path = assemblies / filename
        write_fasta(path, records)
        sources.append((filename, path.resolve(), checksum(path)))
        rows.append(
            {
                "metagenome_id": filename,
                "case_group": "short_grid" if truth_class == "edited_fragment" else "special",
                "truth_class": truth_class,
                "expected_candidate": str(expected).lower(),
                "trace_length": str(length),
                "target_identity": str(identity),
                "orientation": orientation,
                "standalone": str(standalone).lower(),
                "query_intervals_json": json.dumps(intervals, separators=(",", ":")),
            }
        )

    for length in LENGTHS:
        for identity in IDENTITIES:
            start = (length * 17 + identity * 31) % (len(query) - length)
            source = query[start : start + length]
            for orientation in ("forward", "reverse"):
                suffix = "f" if orientation == "forward" else "r"
                case_id = f"short_l{length:04d}_i{identity:03d}_{suffix}"
                edited = mutate(source, identity, length * 1000 + identity)
                if orientation == "reverse":
                    edited = revcomp(edited)
                add_case(
                    case_id,
                    [("background", background), ("trace", edited)],
                    [[start, start + length]],
                    length,
                    identity,
                    orientation,
                    "edited_fragment",
                    True,
                )

    add_case(
        "origin_l0160",
        [("background", background), ("trace", query[-80:] + query[:80])],
        [[len(query) - 80, len(query)], [0, 80]],
        160,
        100,
        "forward",
        "origin_crossing",
        True,
        True,
    )
    add_case(
        "multi_fragments",
        [
            ("background", background),
            ("trace_a", query[900:1060]),
            ("trace_b", query[4200:4360]),
            ("trace_c", revcomp(query[8100:8260])),
        ],
        [[900, 1060], [4200, 4360], [8100, 8260]],
        480,
        100,
        "mixed",
        "separate_fragments",
        True,
    )
    add_case(
        "overlap_fragments",
        [
            ("background", background),
            ("trace_a", query[3000:3300]),
            ("trace_b", query[3200:3500]),
        ],
        [[3000, 3500]],
        500,
        100,
        "forward",
        "overlap",
        True,
    )
    add_case(
        "rare_l0160",
        [("rare_trace", query[7400:7560])],
        [[7400, 7560]],
        160,
        100,
        "forward",
        "rare_standalone",
        True,
        True,
    )
    add_case(
        "repeat_shared",
        [("chromosome_repeat", repeated)],
        [[1800, 2048], [6200, 6448]],
        len(repeated),
        100,
        "forward",
        "chromosome_shared",
        False,
        True,
    )
    add_case(
        "integrated_trace",
        [("chromosome", background[:2000] + query[5200:5700] + background[2000:])],
        [[5200, 5700]],
        500,
        100,
        "forward",
        "integrated_element",
        False,
    )
    add_case(
        "unrelated",
        [("chromosome", dna(rng, 8_000))],
        [],
        0,
        0,
        "none",
        "unrelated",
        False,
    )

    query_path = args.output / "query.fasta"
    write_fasta(query_path, [("anonymous_query", query)])
    with (args.output / "sources.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metagenome_id", "resource_uri", "sha256"))
        writer.writerows(sources)
    fields = list(rows[0])
    with (args.output / "truth.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    metadata = {
        "schema_version": "1.0.0",
        "seed": "0x4A414D53484F5254",
        "query_length": len(query),
        "positive_cases": sum(row["expected_candidate"] == "true" for row in rows),
        "negative_cases": sum(row["expected_candidate"] == "false" for row in rows),
        "case_count": len(rows),
        "lengths": list(LENGTHS),
        "identities": list(IDENTITIES),
        "query_sha256": checksum(query_path),
    }
    (args.output / "dataset.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(metadata, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
