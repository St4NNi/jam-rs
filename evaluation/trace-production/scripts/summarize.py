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
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "1.0.0"
ABLATION_NAMES = (
    "local_aligner",
    "mosaic",
    "dense_k31_rescue",
    "k21_rescue",
    "common_sequence_evidence",
    "selective_reads",
)
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
    if isinstance(value, str):
        try:
            parsed = float(value)
        except ValueError:
            return None
        return parsed if math.isfinite(parsed) else None
    return None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_trace(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    try:
        if path.suffix.lower() in {".zst", ".zstd"}:
            completed = subprocess.run(
                ["zstd", "--decompress", "--stdout", str(path)],
                check=True,
                stdout=subprocess.PIPE,
            )
            text = completed.stdout.decode("utf-8")
        else:
            text = path.read_text(encoding="utf-8")
    except (OSError, subprocess.SubprocessError, UnicodeError):
        return []
    records: list[dict[str, Any]] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        try:
            value = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(value, dict):
            records.append(value)
    return records


def interval_union(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted((item for item in intervals if item[0] <= item[1])):
        if not merged or start > merged[-1][1]:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged


def interval_length(intervals: list[tuple[int, int]]) -> int:
    return sum(end - start for start, end in interval_union(intervals))


def interval_overlap(left: list[tuple[int, int]], right: list[tuple[int, int]]) -> int:
    return sum(
        max(0, min(left_end, right_end) - max(left_start, right_start))
        for left_start, left_end in interval_union(left)
        for right_start, right_end in interval_union(right)
    )


def parse_intervals(value: Any) -> list[tuple[int, int]]:
    if value is None or value == "":
        return []
    if isinstance(value, str):
        try:
            value = json.loads(value)
        except json.JSONDecodeError:
            return []
    if isinstance(value, dict):
        value = value.get("intervals", [value])
    if not isinstance(value, list):
        return []
    intervals: list[tuple[int, int]] = []
    for item in value:
        if isinstance(item, dict):
            start, end = item.get("start"), item.get("end")
        elif isinstance(item, (list, tuple)) and len(item) == 2:
            start, end = item
        else:
            continue
        try:
            start, end = int(start), int(end)
        except (TypeError, ValueError):
            continue
        if 0 <= start <= end:
            intervals.append((start, end))
    return intervals


def parse_bool(value: Any) -> bool | None:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        value = value.strip().lower()
        if value in {"true", "1", "yes", "positive", "supported"}:
            return True
        if value in {"false", "0", "no", "negative", "absent"}:
            return False
    return None


def enrich_truth_metrics(row: dict[str, Any]) -> None:
    truth_value = row.get("truth_path")
    trace_value = row.get("trace_output")
    if not isinstance(truth_value, str) or not truth_value:
        return
    if not isinstance(trace_value, str) or not trace_value:
        return
    truth_path = Path(truth_value)
    trace_path = Path(trace_value)
    if not truth_path.is_file() or not trace_path.is_file():
        parent = row.get("measurement_parent")
        if isinstance(parent, str):
            if not truth_path.is_file():
                truth_path = Path(parent) / truth_path.name
            if not trace_path.is_file():
                trace_path = Path(parent) / trace_path.name
    if not truth_path.is_file() or not trace_path.is_file():
        return
    records = read_trace(trace_path)
    results = {
        str(record.get("metagenome_id")): record
        for record in records
        if record.get("record_type") == "metagenome_result"
    }
    if not results:
        return
    try:
        with truth_path.open(encoding="utf-8", newline="") as handle:
            truth_rows = list(csv.DictReader(handle, delimiter="\t"))
    except (OSError, csv.Error):
        return
    query_id = str(row.get("query_id", ""))
    by_metagenome: dict[str, list[dict[str, str]]] = {}
    for truth in truth_rows:
        if truth.get("query_id") and truth.get("query_id") != query_id:
            continue
        metagenome_id = truth.get("metagenome_id") or truth.get("assembly_id") or truth.get("sample_id")
        if metagenome_id:
            by_metagenome.setdefault(metagenome_id, []).append(truth)
    candidate_total = candidate_hits = final_total = final_hits = false_candidates = negatives = 0
    predicted_bases = expected_bases = overlap_bases = 0
    for metagenome_id, truth_group in by_metagenome.items():
        result = results.get(metagenome_id)
        has_candidate = bool(result and result.get("candidate") is not None)
        has_alignment = bool(result and result.get("alignments"))
        expected_candidates = [parse_bool(item.get("expected_candidate")) for item in truth_group]
        expected_finals = [parse_bool(item.get("expected_final")) for item in truth_group]
        if True in expected_candidates:
            candidate_total += 1
            candidate_hits += int(has_candidate)
        if False in expected_candidates:
            negatives += 1
            false_candidates += int(has_candidate)
        if True in expected_finals:
            final_total += 1
            final_hits += int(has_alignment)
        expected: list[tuple[int, int]] = []
        for truth in truth_group:
            expected.extend(parse_intervals(
                truth.get("expected_intervals")
                or truth.get("truth_intervals")
                or truth.get("query_intervals_json")
            ))
        expected = interval_union(expected)
        if expected:
            expected_bases += interval_length(expected)
            coverage = (result or {}).get("coverage") or {}
            actual = parse_intervals(coverage.get("primary_intervals") or coverage.get("supported_intervals"))
            predicted_bases += interval_length(actual)
            overlap_bases += interval_overlap(actual, expected)
    if row.get("candidate_recall") is None and candidate_total:
        row["candidate_recall"] = candidate_hits / candidate_total
    if row.get("final_recall") is None and final_total:
        row["final_recall"] = final_hits / final_total
    if row.get("false_candidate_count") is None and negatives:
        row["false_candidate_count"] = false_candidates
    if row.get("negative_pair_count") is None and negatives:
        row["negative_pair_count"] = negatives
    if row.get("interval_precision") is None and predicted_bases:
        row["interval_precision"] = overlap_bases / predicted_bases
    if row.get("interval_recall") is None and expected_bases:
        row["interval_recall"] = overlap_bases / expected_bases


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
            value["measurement_parent"] = str(path.parent)
            rows.append(normalize_raw(value))
            paths.append(str(path.relative_to(raw_dir)))
    return rows, paths


def enrich_from_trace(row: dict[str, Any]) -> dict[str, Any]:
    trace_value = row.get("trace_output")
    if not isinstance(trace_value, str) or not trace_value:
        return row
    trace_path = Path(trace_value)
    if not trace_path.is_file():
        measurement_parent = row.get("measurement_parent")
        if isinstance(measurement_parent, str):
            trace_path = Path(measurement_parent) / trace_path.name
    records = read_trace(trace_path)
    if not records:
        return row
    header = next((record for record in records if record.get("record_type") == "run_header"), {})
    results = [record for record in records if record.get("record_type") == "metagenome_result"]
    row["query_id"] = row.get("query_id") or header.get("query_id") or header.get("plasmid_id") or ""
    row["query_kind"] = row.get("query_kind") or header.get("query_kind") or "unknown"
    row["topology"] = row.get("topology") or header.get("topology_requested") or "auto"
    row["candidate_count"] = sum(record.get("candidate") is not None for record in results)
    row["alignment_count"] = sum(len(record.get("alignments", [])) for record in results if isinstance(record.get("alignments", []), list))
    footer = next((record for record in records if record.get("record_type") == "run_footer"), {})
    if not row.get("resource_metrics") and isinstance(footer.get("resource_metrics"), dict):
        row["resource_metrics"] = footer["resource_metrics"]
    row["trace_result_count"] = len(results)
    row["trace_query_id"] = header.get("query_id") or header.get("plasmid_id")
    return row


def normalize_raw(value: dict[str, Any]) -> dict[str, Any]:
    """Normalize v1 and v2 measurement aliases without changing raw values."""

    row = enrich_from_trace(dict(value))
    row["query_id"] = row.get("query_id") or row.get("plasmid_id") or row.get("query") or ""
    row["query_kind"] = row.get("query_kind") or row.get("kind") or "unknown"
    row["topology"] = row.get("topology") or row.get("topology_requested") or "auto"
    row["query_class"] = row.get("query_class") or row.get("class") or row["query_kind"]
    row["stratum"] = row.get("stratum") or row.get("query_stratum") or row["query_kind"]
    ablation = row.get("ablation")
    if not isinstance(ablation, dict):
        ablation = {}
    row["ablation"] = {
        name: str(ablation.get(name, "unavailable")) for name in ABLATION_NAMES
    }
    for field in ("candidate_recall", "final_recall", "interval_precision", "interval_recall", "gap_boundary_error_bases"):
        if row.get(field) is None and isinstance(row.get("metrics"), dict):
            row[field] = row["metrics"].get(field)
    enrich_truth_metrics(row)
    return row


def key_for(row: dict[str, Any]) -> tuple[Any, ...]:
    ablation = row.get("ablation", {})
    return (
        row.get("track"),
        row.get("query_kind", "unknown"),
        row.get("topology", "auto"),
        row.get("stratum", row.get("query_kind", "unknown")),
        row.get("resource_mode"),
        row.get("cache_state"),
        row.get("cpu_threads"),
        row.get("io_tasks"),
        row.get("sweep_dimension"),
        str(row.get("sweep_value")),
        *(str(ablation.get(name, "unavailable")) for name in ABLATION_NAMES),
        row.get("matched_recall_target"),
    )


def _mean(rows: list[dict[str, Any]], field: str) -> float | None:
    values = [number(row.get(field)) for row in rows]
    values = [value for value in values if value is not None]
    return sum(values) / len(values) if values else None


def aggregate(group_key: tuple[Any, ...], rows: list[dict[str, Any]]) -> dict[str, Any]:
    complete = [row for row in rows if row.get("status") == "complete" and row.get("exit_code") == 0]
    unavailable = sum(row.get("status") == "unavailable" for row in rows)
    failed = sum(row.get("status") == "failed" or (row.get("status") == "complete" and row.get("exit_code") not in (None, 0)) for row in rows)
    result: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "track": group_key[0],
        "query_kind": group_key[1],
        "topology": group_key[2],
        "stratum": group_key[3],
        "resource_mode": group_key[4],
        "cache_state": group_key[5],
        "cpu_threads": group_key[6],
        "io_tasks": group_key[7],
        "sweep_dimension": group_key[8],
        "sweep_value": group_key[9],
        "ablation": {
            name: group_key[10 + index] for index, name in enumerate(ABLATION_NAMES)
        },
        "matched_recall_target": group_key[16],
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
    target = number(group_key[16])
    if target is None:
        result["interval_metrics_at_matched_recall"] = {
            "status": "unavailable",
            "reason": "matched_recall_target was not supplied",
        }
    else:
        matched = [
            row for row in accuracy_rows
            if (candidate_recall := number(row.get("candidate_recall"))) is not None
            and candidate_recall >= target
        ]
        result["interval_metrics_at_matched_recall"] = {
            "status": "available" if matched else "unavailable",
            "target": target,
            "query_count": len(matched),
            "interval_precision": _mean(matched, "interval_precision"),
            "interval_recall": _mean(matched, "interval_recall"),
            "gap_boundary_error_bases": _mean(matched, "gap_boundary_error_bases"),
        }
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
        "schema_version", "track", "query_kind", "query_class", "topology", "stratum",
        "resource_mode", "cache_state", "cpu_threads", "io_tasks", "sweep_dimension", "sweep_value",
        "status", "raw_count", "complete_count", "failed_count",
        "unavailable_count", "wall_median", "wall_p90", "wall_p95", "wall_p99", "wall_bootstrap_low",
        "wall_bootstrap_high", "cpu_median", "rss_median", "candidate_recall", "final_recall",
        "false_candidate_count", "false_candidate_rate", "interval_precision_mean", "interval_recall_mean",
        "gap_boundary_error_bases_mean", "break_even_queries_wall", "break_even_queries_cpu",
        "matched_recall_target", "matched_interval_precision", "matched_interval_recall",
        "matched_gap_boundary_error_bases", "matched_recall_status",
        *("ablation_" + name for name in ABLATION_NAMES),
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
            matched = group.get("interval_metrics_at_matched_recall", {})
            row = {
                **group,
                "query_class": group.get("query_kind", "unknown"),
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
                "matched_interval_precision": matched.get("interval_precision"),
                "matched_interval_recall": matched.get("interval_recall"),
                "matched_gap_boundary_error_bases": matched.get("gap_boundary_error_bases"),
                "matched_recall_status": matched.get("status"),
                **{
                    "ablation_" + name: (group.get("ablation") or {}).get(name)
                    for name in ABLATION_NAMES
                },
                **{field: io_totals.get(field) for field in IO_FIELDS},
            }
            writer.writerow({key: "NA" if value is None else value for key, value in row.items()})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--matched-recall",
        type=float,
        help="report interval metrics for complete query rows meeting this candidate recall",
    )
    args = parser.parse_args()
    raw_dir = args.raw_dir.resolve()
    output_dir = args.output_dir.resolve()
    if not raw_dir.is_dir():
        raise SystemExit(f"raw directory is not a directory: {raw_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    rows, raw_paths = load_raw(raw_dir)
    if args.matched_recall is not None:
        if not 0.0 <= args.matched_recall <= 1.0:
            raise SystemExit("--matched-recall must be between 0 and 1")
        for row in rows:
            row["matched_recall_target"] = args.matched_recall
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
        "query_strata": sorted({str(row.get("stratum", "unknown")) for row in rows}),
        "query_kinds": sorted({str(row.get("query_kind", "unknown")) for row in rows}),
        "topologies": sorted({str(row.get("topology", "auto")) for row in rows}),
        "ablation_dimensions": list(ABLATION_NAMES),
        "statistics": {"quantiles": ["median", "p90", "p95", "p99"], "bootstrap_iterations": 2000, "bootstrap_seed": "sha256(group-key)"},
        "matched_recall_policy": "interval metrics are reported at a supplied matched_recall_target; without one they remain unavailable",
        "accuracy_note": "final_recall means aligned trace recovery where expected_final is present; it is not an autonomous plasmid confirmation call",
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_tsv(output_dir / "summary.tsv", groups)
    (output_dir / "raw-manifest.json").write_text(json.dumps({"schema_version": SCHEMA_VERSION, "files": raw_hashes}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
