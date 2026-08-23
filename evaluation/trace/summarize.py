#!/usr/bin/env python3
"""Convert raw trace measurement JSON files into a stable TSV summary."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


FIELDS = [
    "stage",
    "label",
    "profile",
    "cache_state",
    "threads",
    "exit_code",
    "wall_seconds",
    "cpu_seconds",
    "peak_rss_kib",
    "candidate_count",
    "alignment_count",
    "coverage_mean",
    "gap_bases",
    "largest_gap",
    "remote_bytes",
    "cache_bytes",
    "cache_hits",
    "metadata_requests",
    "range_requests",
    "stream_requests",
    "server_requests",
    "server_bytes",
    "exact_candidate_recovery",
    "exact_alignment_recovery",
    "partial_candidate_recovery",
    "partial_alignment_recovery",
    "split_overlap_candidate_recovery",
    "split_overlap_alignment_recovery",
    "short_fragment_candidate_recovery",
    "short_fragment_alignment_recovery",
]


def value(data: dict, key: str) -> object:
    if key == "cpu_seconds":
        return float(data.get("user_seconds", 0.0)) + float(data.get("system_seconds", 0.0))
    if key in {"remote_bytes", "cache_bytes", "cache_hits", "metadata_requests", "range_requests", "stream_requests"}:
        return data.get("resource_metrics", {}).get(key, 0)
    if key == "server_requests":
        return data.get("server_metrics", {}).get("requests", 0)
    if key == "server_bytes":
        return data.get("server_metrics", {}).get("bytes_served", 0)
    raw = data.get(key)
    return "NA" if raw is None else raw


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--measurements", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for path in sorted(args.measurements.glob("*.json")):
        data = json.loads(path.read_text(encoding="utf-8"))
        row = {field: value(data, field) for field in FIELDS}
        row["measurement_file"] = path.name
        rows.append(row)
    fields = ["measurement_file", *FIELDS]
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
