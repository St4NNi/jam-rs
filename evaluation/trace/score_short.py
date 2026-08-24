#!/usr/bin/env python3
"""Score candidate and alignment recovery for the anonymous short-trace matrix."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def union(intervals: list[list[int]]) -> list[list[int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if start >= end:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return merged


def overlap(left: list[list[int]], right: list[list[int]]) -> int:
    total = 0
    left_index = 0
    right_index = 0
    while left_index < len(left) and right_index < len(right):
        start = max(left[left_index][0], right[right_index][0])
        end = min(left[left_index][1], right[right_index][1])
        total += max(0, end - start)
        if left[left_index][1] <= right[right_index][1]:
            left_index += 1
        else:
            right_index += 1
    return total


def load_trace(path: Path) -> dict[str, dict]:
    results: dict[str, dict] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            record = json.loads(line)
            if record.get("record_type") == "metagenome_result":
                results[record["metagenome_id"]] = record
    return results


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trace", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.truth.open(encoding="utf-8", newline="") as handle:
        truth = list(csv.DictReader(handle, delimiter="\t"))
    results = load_trace(args.trace)
    positive = [row for row in truth if row["expected_candidate"] == "true"]
    negative = [row for row in truth if row["expected_candidate"] == "false"]
    scored: list[dict] = []
    for row in truth:
        result = results.get(row["metagenome_id"])
        expected = union(json.loads(row["query_intervals_json"]))
        observed = union(
            [
                [int(interval["start"]), int(interval["end"])]
                for interval in (result or {}).get("coverage", {}).get("primary_intervals", [])
            ]
        )
        truth_bases = sum(end - start for start, end in expected)
        recovered = overlap(expected, observed)
        scored.append(
            {
                **row,
                "candidate": result is not None and result.get("candidate") is not None,
                "aligned": bool((result or {}).get("alignments")),
                "truth_bases": truth_bases,
                "recovered_bases": recovered,
                "truth_base_recall": recovered / truth_bases if truth_bases else None,
                "supporting_contigs": len(
                    ((result or {}).get("primary_fragment_mosaic") or {}).get(
                        "supporting_contigs", []
                    )
                ),
            }
        )

    positive_scored = [row for row in scored if row["expected_candidate"] == "true"]
    negative_scored = [row for row in scored if row["expected_candidate"] == "false"]
    groups: dict[tuple[int, int, str], list[dict]] = {}
    for row in positive_scored:
        if row["case_group"] != "short_grid":
            continue
        key = (int(row["trace_length"]), int(row["target_identity"]), row["orientation"])
        groups.setdefault(key, []).append(row)
    grid = []
    for (length, identity, orientation), rows in sorted(groups.items()):
        truth_bases = sum(row["truth_bases"] for row in rows)
        recovered = sum(row["recovered_bases"] for row in rows)
        grid.append(
            {
                "trace_length": length,
                "target_identity": identity,
                "orientation": orientation,
                "cases": len(rows),
                "candidates": sum(row["candidate"] for row in rows),
                "aligned": sum(row["aligned"] for row in rows),
                "truth_base_recall": recovered / truth_bases if truth_bases else None,
            }
        )
    exact_160 = next(row for row in scored if row["truth_class"] == "rare_standalone")
    special = {
        row["truth_class"]: {
            "candidate": row["candidate"],
            "aligned": row["aligned"],
            "truth_base_recall": row["truth_base_recall"],
            "supporting_contigs": row["supporting_contigs"],
        }
        for row in scored
        if row["case_group"] == "special"
    }
    total_truth = sum(row["truth_bases"] for row in positive_scored)
    total_recovered = sum(row["recovered_bases"] for row in positive_scored)
    output = {
        "schema_version": "1.0.0",
        "positive_cases": len(positive),
        "candidate_cases": sum(row["candidate"] for row in positive_scored),
        "aligned_cases": sum(row["aligned"] for row in positive_scored),
        "candidate_recall": sum(row["candidate"] for row in positive_scored) / len(positive),
        "alignment_recall": sum(row["aligned"] for row in positive_scored) / len(positive),
        "truth_base_recall": total_recovered / total_truth,
        "exact_160_standalone": {
            "candidate": exact_160["candidate"],
            "aligned": exact_160["aligned"],
            "truth_base_recall": exact_160["truth_base_recall"],
        },
        "biological_negatives": {
            "cases": len(negative),
            "selected": sum(row["candidate"] for row in negative_scored),
            "aligned": sum(row["aligned"] for row in negative_scored),
        },
        "special_cases": special,
        "grid": grid,
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in output.items() if key != "grid"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
