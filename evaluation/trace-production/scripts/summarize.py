#!/usr/bin/env python3
"""Regenerate production benchmark summaries from raw measurement JSON files.

Only completed raw cases contribute to timing quantiles.  Accuracy metrics are
weighted from truth counts rather than averaged as percentages.  Bootstrap
intervals resample query-level means with a stable seed derived from the group
key, making summary regeneration deterministic.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "1.0.0"
NUMERIC_FIELDS = (
    "wall_seconds",
    "cpu_seconds",
    "user_seconds",
    "system_seconds",
    "peak_rss_kib",
    "candidate_count",
    "alignment_count",
    "candidate_recall",
    "final_recall",
    "false_candidate_count",
    "interval_precision",
    "interval_recall",
    "gap_boundary_error_bases",
)
IO_FIELDS = (
    "metadata_requests",
    "head_requests",
    "get_requests",
    "range_requests",
    "requested_bytes",
    "returned_bytes",
    "decoded_bytes",
    "cache_bytes",
    "cache_hits",
    "cache_misses",
    "stale_cache_rejections",
    "retries",
)


def number(value: Any) -> float | None:
    if isinstance(value, bool):
        return float(value)
    if isinstance(value, (int, float)) and math.isfinite(float(value)):
        return float(value)
    return None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def percentile(values: list[float], percentile_value: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * percentile_value
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def bootstrap(values: list[float], key: str, iterations: int = 2000) -> tuple[float | None, float | None]:
    if not values:
        return None, None
    if len(values) == 1:
        return values[0], values[0]
    seed = int.from_bytes(hashlib.sha256(key.encode("utf-8")).digest()[:8], "big")
    rng = random.Random(seed)
    estimates: list[float] = []
    for _ in range(iterations):
        sample = [values[rng.randrange(len(values))] for _ in values]
        estimates.append(sum(sample) / len(sample))
    return percentile(estimates, 0.025), percentile(estimates, 0.975)


def load_raw(raw_dir: Path) -> tuple[list[dict[str, Any]], list[str]]:
    rows: list[dict[str, Any]] = []
    paths: list[str] = []
    for path in sorted(raw_dir.rglob("measurement.json")):
        try:
            value = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if isinstance(value, dict):
            rows.append(value)
            paths.append(str(path.relative_to(raw_dir)))
    return rows, paths


def key_for(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        row.get("track"),
        row.get("resource_mode"),
        row.get("cache_state"),
        row.get("cpu_threads"),
        row.get("io_tasks"),
        row.get("sweep_dimension"),
        str(row.get("sweep_value")),
    )


def aggregate(group_key: tuple[Any, ...], rows: list[dict[str, Any]]) -> dict[str, Any]:
    complete = [row for row in rows if row.get("status") == "complete" and row.get("exit_code") == 0]
    unavailable = sum(row.get("status") == "unavailable" for row in rows)
    failed = sum(row.get("status") == "failed" or (row.get("status") == "complete" and row.get("exit_code") not in (None, 0)) for row in rows)
    result: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "track": group_key[0],
        "resource_mode": group_key[1],
        "cache_state": group_key[2],
        "cpu_threads": group_key[3],
        "io_tasks": group_key[4],
        "sweep_dimension": group_key[5],
        "sweep_value": group_key[6],
        "raw_count": len(rows),
        "complete_count": len(complete),
        "failed_count": failed,
        "unavailable_count": unavailable,
        "status": "complete" if complete else ("unavailable" if unavailable and not failed else "failed"),
    }
    for field in NUMERIC_FIELDS:
        values = [item for item in (number(row.get(field)) for row in complete) if item is not None]
        query_values: dict[str, list[float]] = {}
        for row in complete:
            value = number(row.get(field))
            if value is not None:
                query_values.setdefault(str(row.get("query_id", "")), []).append(value)
        query_means = [sum(values) / len(values) for values in query_values.values() if values]
        result[field] = {
            "median": percentile(values, 0.5),
            "p90": percentile(values, 0.9),
            "p95": percentile(values, 0.95),
            "p99": percentile(values, 0.99),
            "bootstrap_95": bootstrap(query_means, "|".join(map(str, group_key)) + ":" + field),
            "query_count": len(query_means),
        }
    # Truth metrics are recomputed from counts where available.  This avoids
    # overweighting a query with more contigs or repeated order runs.
    accuracy_by_query: dict[str, dict[str, Any]] = {}
    for row in sorted(complete, key=lambda item: (str(item.get("query_id", "")), int(item.get("order", 0) or 0), int(item.get("repeat", 0) or 0), str(item.get("cache_state", "")))):
        accuracy_by_query.setdefault(str(row.get("query_id", "")), row)
    accuracy_rows = list(accuracy_by_query.values())
    candidate_rows = [row for row in accuracy_rows if number(row.get("candidate_recall")) is not None]
    final_rows = [row for row in accuracy_rows if number(row.get("final_recall")) is not None]
    false_counts = [value for value in (number(row.get("false_candidate_count")) for row in accuracy_rows) if value is not None]
    negative_counts = [value for value in (number(row.get("negative_pair_count")) for row in accuracy_rows) if value is not None]
    result["candidate_recall"] = {
        "median": percentile([number(row["candidate_recall"]) for row in candidate_rows if number(row.get("candidate_recall")) is not None], 0.5),
        "query_count": len(candidate_rows),
    }
    result["final_recall"] = {
        "median": percentile([number(row["final_recall"]) for row in final_rows if number(row.get("final_recall")) is not None], 0.5),
        "query_count": len(final_rows),
    }
    false_total = sum(value for value in false_counts if value is not None)
    negative_total = sum(value for value in negative_counts if value is not None)
    result["false_candidate_count"] = false_total if false_counts else None
    result["negative_pair_count"] = negative_total if negative_counts else None
    result["false_candidate_rate"] = false_total / negative_total if negative_total else None
    for field in ("interval_precision", "interval_recall", "gap_boundary_error_bases"):
        values = [number(row.get(field)) for row in accuracy_rows]
        values = [value for value in values if value is not None]
        result[field + "_mean"] = sum(values) / len(values) if values else None
    io_totals: dict[str, int | float] = {}
    for row in complete:
        metrics = row.get("resource_metrics", {})
        if isinstance(metrics, dict):
            for field in IO_FIELDS:
                value = number(metrics.get(field))
                if value is not None:
                    io_totals[field] = io_totals.get(field, 0) + value
    result["io_totals"] = io_totals
    phase_values: dict[str, list[float]] = {}
    for row in complete:
        phases = row.get("phase_timings", {})
        if not isinstance(phases, dict):
            continue
        for phase, details in phases.items():
            if isinstance(details, dict):
                value = number(details.get("wall_seconds"))
                if value is not None:
                    phase_values.setdefault(str(phase), []).append(value)
    result["phase_timing_medians"] = {
        phase: percentile(values, 0.5) for phase, values in sorted(phase_values.items())
    }
    archive_rows = [row.get("archive_economics", {}) for row in complete if isinstance(row.get("archive_economics"), dict)]
    break_even_wall = [number(item.get("break_even_queries_wall")) for item in archive_rows]
    break_even_wall = [item for item in break_even_wall if item is not None]
    break_even_cpu = [number(item.get("break_even_queries_cpu")) for item in archive_rows]
    break_even_cpu = [item for item in break_even_cpu if item is not None]
    result["break_even_queries_wall"] = sum(break_even_wall) / len(break_even_wall) if break_even_wall else None
    result["break_even_queries_cpu"] = sum(break_even_cpu) / len(break_even_cpu) if break_even_cpu else None
    return result


def write_tsv(path: Path, groups: list[dict[str, Any]]) -> None:
    fields = [
        "schema_version", "track", "resource_mode", "cache_state", "cpu_threads", "io_tasks",
        "sweep_dimension", "sweep_value", "status", "raw_count", "complete_count", "failed_count",
        "unavailable_count", "wall_median", "wall_p90", "wall_p95", "wall_p99", "wall_bootstrap_low",
        "wall_bootstrap_high", "cpu_median", "rss_median", "candidate_recall", "final_recall",
        "false_candidate_count", "false_candidate_rate", "interval_precision_mean", "interval_recall_mean",
        "gap_boundary_error_bases_mean", "break_even_queries_wall", "break_even_queries_cpu",
        "metadata_requests", "head_requests", "get_requests", "range_requests",
        "requested_bytes", "returned_bytes", "decoded_bytes", "cache_bytes", "cache_hits",
        "cache_misses", "stale_cache_rejections",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for group in groups:
            wall = group.get("wall_seconds", {})
            cpu = group.get("cpu_seconds", {})
            rss = group.get("peak_rss_kib", {})
            io_totals = group.get("io_totals", {})
            row = {
                **group,
                "wall_median": wall.get("median"),
                "wall_p90": wall.get("p90"),
                "wall_p95": wall.get("p95"),
                "wall_p99": wall.get("p99"),
                "wall_bootstrap_low": (wall.get("bootstrap_95") or [None, None])[0],
                "wall_bootstrap_high": (wall.get("bootstrap_95") or [None, None])[1],
                "cpu_median": cpu.get("median"),
                "rss_median": rss.get("median"),
                "candidate_recall": (group.get("candidate_recall") or {}).get("median"),
                "final_recall": (group.get("final_recall") or {}).get("median"),
                **{field: io_totals.get(field) for field in IO_FIELDS},
            }
            writer.writerow({key: "NA" if value is None else value for key, value in row.items()})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    raw_dir = args.raw_dir.resolve()
    output_dir = args.output_dir.resolve()
    if not raw_dir.is_dir():
        raise SystemExit(f"raw directory is not a directory: {raw_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    rows, raw_paths = load_raw(raw_dir)
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(key_for(row), []).append(row)
    groups = [aggregate(key, grouped[key]) for key in sorted(grouped, key=lambda item: tuple(str(part) for part in item))]
    raw_hashes = {path: sha256(raw_dir / path) for path in raw_paths}
    summary = {
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "raw_directory": str(raw_dir),
        "raw_measurement_count": len(rows),
        "raw_measurements": raw_hashes,
        "groups": groups,
        "statistics": {"quantiles": ["median", "p90", "p95", "p99"], "bootstrap_iterations": 2000, "bootstrap_seed": "sha256(group-key)"},
        "accuracy_note": "final_recall means aligned trace recovery where expected_final is present; it is not an autonomous plasmid confirmation call",
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_tsv(output_dir / "summary.tsv", groups)
    (output_dir / "raw-manifest.json").write_text(json.dumps({"schema_version": SCHEMA_VERSION, "files": raw_hashes}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
