#!/usr/bin/env python3
"""Run one production trace case and record raw, reproducible measurements.

This module deliberately does not build software, install comparators, or
interpret missing inputs as negative results.  It starts one child process so
``RUSAGE_CHILDREN`` gives an uncontaminated wall/CPU/RSS record.  The JSONL
trace is parsed only for counters and truth comparison; the raw trace remains
the source of biological evidence.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import resource
import subprocess
import time
import urllib.parse
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "1.0.0"
SECRET_KEYS = {"x-amz-credential", "x-amz-signature", "x-amz-security-token", "signature", "token", "password"}
METRIC_KEYS = (
    "metadata_requests",
    "head_requests",
    "get_requests",
    "range_requests",
    "stream_requests",
    "requested_bytes",
    "returned_bytes",
    "decoded_bytes",
    "remote_bytes",
    "cache_bytes",
    "cache_hits",
    "cache_misses",
    "cache_evictions",
    "stale_cache_rejections",
    "retries",
    "full_object_fallbacks",
    "seed_buckets_read",
    "sequence_blocks_read",
)


def redact_text(value: str) -> str:
    """Remove credentials from command/URI text before it enters an artifact."""

    value = re.sub(r"(https?://)([^/@\s]+):([^/@\s]+)@", r"\1<redacted>@", value)
    parsed = urllib.parse.urlsplit(value)
    if parsed.scheme and parsed.netloc and parsed.query:
        pairs = urllib.parse.parse_qsl(parsed.query, keep_blank_values=True)
        safe = []
        for key, item in pairs:
            safe.append((key, "<redacted>" if key.lower() in SECRET_KEYS else item))
        value = urllib.parse.urlunsplit(
            (parsed.scheme, parsed.netloc, parsed.path, urllib.parse.urlencode(safe), parsed.fragment)
        )
    return value


def redacted_command(command: list[str]) -> list[str]:
    return [redact_text(item) for item in command]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    if path.suffix.lower() in {".zst", ".zstd"}:
        completed = subprocess.run(
            ["zstd", "--decompress", "--stdout", str(path)],
            check=True,
            stdout=subprocess.PIPE,
        )
        text = completed.stdout.decode("utf-8")
    else:
        text = path.read_text(encoding="utf-8")
    records: list[dict[str, Any]] = []
    for line in text.splitlines():
        if line.strip():
            value = json.loads(line)
            if isinstance(value, dict):
                records.append(value)
    return records


def number(value: Any) -> int | float | None:
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, (int, float)):
        return value
    return None


def nested_metrics(value: Any) -> dict[str, int | float]:
    """Collect only known metric names from headers/results/footers."""

    wanted = set(METRIC_KEYS)
    result: dict[str, int | float] = {}
    if isinstance(value, dict):
        for key, child in value.items():
            normalized = str(key).lower().replace("-", "_")
            if normalized in wanted:
                item = number(child)
                if item is not None:
                    result[normalized] = result.get(normalized, 0) + item
            for child_key, child_value in nested_metrics(child).items():
                result[child_key] = result.get(child_key, 0) + child_value
    elif isinstance(value, list):
        for child in value:
            for child_key, child_value in nested_metrics(child).items():
                result[child_key] = result.get(child_key, 0) + child_value
    return result


def parse_bool(value: Any) -> bool | None:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "yes", "1", "positive", "present"}:
            return True
        if normalized in {"false", "no", "0", "negative", "absent"}:
            return False
    return None


def parse_intervals(value: Any) -> list[tuple[int, int]]:
    if value is None or value == "":
        return []
    if isinstance(value, str):
        try:
            value = json.loads(value)
        except json.JSONDecodeError:
            intervals = []
            for item in value.split(","):
                parts = item.strip().replace(":", "-").split("-")
                if len(parts) == 2 and all(part.strip().isdigit() for part in parts):
                    intervals.append((int(parts[0]), int(parts[1])))
            return intervals
    if isinstance(value, dict):
        start, end = value.get("start"), value.get("end")
        if isinstance(start, int) and isinstance(end, int) and 0 <= start <= end:
            return [(start, end)]
        return []
    if not isinstance(value, list):
        return []
    intervals: list[tuple[int, int]] = []
    for item in value:
        if isinstance(item, dict):
            start, end = item.get("start"), item.get("end")
            if isinstance(start, int) and isinstance(end, int) and 0 <= start <= end:
                intervals.append((start, end))
        elif isinstance(item, (list, tuple)) and len(item) == 2:
            start, end = item
            if isinstance(start, int) and isinstance(end, int) and 0 <= start <= end:
                intervals.append((start, end))
    return intervals


def union_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted((item for item in intervals if item[0] <= item[1])):
        if not merged or start > merged[-1][1]:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged


def interval_length(intervals: list[tuple[int, int]]) -> int:
    return sum(end - start for start, end in union_intervals(intervals))


def interval_intersection(left: list[tuple[int, int]], right: list[tuple[int, int]]) -> int:
    total = 0
    for left_start, left_end in union_intervals(left):
        for right_start, right_end in union_intervals(right):
            total += max(0, min(left_end, right_end) - max(left_start, right_start))
    return total


def boundary_error(actual: list[tuple[int, int]], expected: list[tuple[int, int]]) -> int | None:
    if not actual or not expected:
        return None
    actual = union_intervals(actual)
    expected = union_intervals(expected)
    if len(actual) != len(expected):
        return abs(interval_length(actual) - interval_length(expected)) + abs(len(actual) - len(expected))
    return sum(
        abs(actual_start - expected_start) + abs(actual_end - expected_end)
        for (actual_start, actual_end), (expected_start, expected_end) in zip(actual, expected)
    )


def load_truth(path: Path | None) -> dict[tuple[str, str], dict[str, str]]:
    if path is None or not path.is_file():
        return {}
    with path.open(encoding="utf-8", newline="") as handle:
        return {
            (str(row.get("query_id", "")), str(row.get("metagenome_id", ""))): row
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("query_id") and row.get("metagenome_id")
        }


def result_metrics(records: list[dict[str, Any]], truth_path: Path | None) -> dict[str, Any]:
    results = {
        (str(record.get("plasmid_id", "")), str(record.get("metagenome_id", ""))): record
        for record in records
        if record.get("record_type") == "metagenome_result"
    }
    footer = next((record for record in records if record.get("record_type") == "run_footer"), {})
    metrics = nested_metrics(footer)
    alignments = sum(
        len(record.get("alignments", []))
        for record in results.values()
        if isinstance(record.get("alignments", []), list)
    )
    candidates = sum(record.get("candidate") is not None for record in results.values())
    aligned = sum(bool(record.get("alignments")) for record in results.values())
    truth = load_truth(truth_path)
    candidate_tp = candidate_total = final_tp = final_total = false_candidates = negatives = 0
    overlap_bases = predicted_bases = expected_bases = 0
    gap_errors: list[int] = []
    for key, row in truth.items():
        record = results.get(key)
        has_candidate = bool(record and record.get("candidate") is not None)
        has_alignment = bool(record and record.get("alignments"))
        expected_candidate = parse_bool(row.get("expected_candidate"))
        expected_final = parse_bool(row.get("expected_final"))
        if expected_candidate is True:
            candidate_total += 1
            candidate_tp += int(has_candidate)
        elif expected_candidate is False:
            negatives += 1
            false_candidates += int(has_candidate)
        if expected_final is True:
            final_total += 1
            final_tp += int(has_alignment)
        expected = parse_intervals(row.get("expected_intervals") or row.get("truth_intervals"))
        if expected:
            expected_bases += interval_length(expected)
            if record:
                coverage = record.get("coverage") or {}
                actual = parse_intervals(coverage.get("primary_intervals"))
                predicted_bases += interval_length(actual)
                overlap_bases += interval_intersection(actual, expected)
                gaps = parse_intervals(row.get("expected_gaps"))
                actual_gaps = parse_intervals(coverage.get("gaps"))
                error = boundary_error(actual_gaps, gaps)
                if error is not None:
                    gap_errors.append(error)
    if not truth:
        accuracy = {
            "truth_status": "unavailable",
            "candidate_recall": None,
            "final_recall": None,
            "false_candidate_count": None,
            "negative_pair_count": None,
            "interval_precision": None,
            "interval_recall": None,
            "gap_boundary_error_bases": None,
        }
    else:
        accuracy = {
            "truth_status": "available",
            "candidate_recall": candidate_tp / candidate_total if candidate_total else None,
            "final_recall": final_tp / final_total if final_total else None,
            "false_candidate_count": false_candidates,
            "negative_pair_count": negatives,
            "interval_precision": overlap_bases / predicted_bases if predicted_bases else None,
            "interval_recall": overlap_bases / expected_bases if expected_bases else None,
            "gap_boundary_error_bases": sum(gap_errors) / len(gap_errors) if gap_errors else None,
        }
    return {
        "metagenome_count": len(results),
        "candidate_count": candidates,
        "alignment_count": alignments,
        "resource_metrics": {key: metrics.get(key, 0) for key in METRIC_KEYS if key in metrics},
        **accuracy,
    }


def unavailable_record(args: argparse.Namespace, reason: str) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "unavailable",
        "unavailable_reason": reason,
        "command": redacted_command(args.command),
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
        "exit_code": None,
        "wall_seconds": None,
        "user_seconds": None,
        "system_seconds": None,
        "peak_rss_kib": None,
        "resource_metrics": {},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--trace-output", type=Path)
    parser.add_argument("--truth", type=Path)
    parser.add_argument("--stdout", type=Path)
    parser.add_argument("--stderr", type=Path)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--phase-json", type=Path)
    parser.add_argument("--track", default="")
    parser.add_argument("--query-id", default="")
    parser.add_argument("--metagenome-set", default="")
    parser.add_argument("--resource-mode", default="local")
    parser.add_argument("--cache-state", default="local-process-cold")
    parser.add_argument("--cpu-threads", type=int, default=1)
    parser.add_argument("--io-tasks", type=int, default=1)
    parser.add_argument("--order", type=int, default=0)
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument("--sweep-dimension", default="baseline")
    parser.add_argument("--sweep-value", default="baseline")
    parser.add_argument("--archive-json", type=Path)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = list(args.command)
    if command and command[0] == "--":
        command = command[1:]
    args.command = command
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if not command:
        parser.error("a command is required after --")

    stdout_handle = args.stdout.open("w", encoding="utf-8") if args.stdout else None
    stderr_handle = args.stderr.open("w", encoding="utf-8") if args.stderr else None
    started = datetime.now(timezone.utc).isoformat()
    start = time.perf_counter()
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    try:
        completed = subprocess.run(command, check=False, stdout=stdout_handle, stderr=stderr_handle)
    finally:
        if stdout_handle:
            stdout_handle.close()
        if stderr_handle:
            stderr_handle.close()
    elapsed = time.perf_counter() - start
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    trace_exists = bool(args.trace_output and args.trace_output.is_file())
    records: list[dict[str, Any]] = []
    trace_error: str | None = None
    if trace_exists:
        try:
            records = read_jsonl(args.trace_output)
        except (OSError, ValueError, subprocess.SubprocessError) as error:
            trace_error = f"trace output could not be parsed: {error}"
        if not records:
            trace_error = "trace output contained no JSON records"
    elif completed.returncode == 0:
        trace_error = "trace output was not created"
    result = result_metrics(records, args.truth) if records else result_metrics([], args.truth)
    metadata: dict[str, Any] = {}
    if args.metadata and args.metadata.is_file():
        try:
            metadata = json.loads(args.metadata.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            metadata = {"status": "unavailable", "reason": "invalid metadata JSON"}
    phase_timings: dict[str, Any] = {}
    if args.phase_json and args.phase_json.is_file():
        try:
            phase_timings = json.loads(args.phase_json.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            phase_timings = {"status": "unavailable", "reason": "invalid phase JSON"}
    if not phase_timings:
        phase_timings = {
            "trace_total": {
                "wall_seconds": elapsed,
                "cpu_seconds": (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime),
                "source": "process measurement; internal phases not exposed by this executable",
            }
        }
    archive: dict[str, Any] = {}
    if args.archive_json and args.archive_json.is_file():
        try:
            archive = json.loads(args.archive_json.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            archive = {"status": "unavailable", "reason": "invalid archive JSON"}
    status = "complete" if completed.returncode == 0 and trace_error is None else "failed"
    output = {
        "schema_version": SCHEMA_VERSION,
        "status": status,
        "track": args.track,
        "query_id": args.query_id,
        "metagenome_set": args.metagenome_set,
        "resource_mode": args.resource_mode,
        "cache_state": args.cache_state,
        "cpu_threads": args.cpu_threads,
        "io_tasks": args.io_tasks,
        "order": args.order,
        "repeat": args.repeat,
        "sweep_dimension": args.sweep_dimension,
        "sweep_value": args.sweep_value,
        "command": redacted_command(command),
        "started_at_utc": started,
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
        "exit_code": completed.returncode,
        "wall_seconds": elapsed,
        "user_seconds": after.ru_utime - before.ru_utime,
        "system_seconds": after.ru_stime - before.ru_stime,
        "cpu_seconds": (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime),
        "peak_rss_kib": after.ru_maxrss,
        "trace_output": str(args.trace_output) if args.trace_output else None,
        "trace_sha256": sha256(args.trace_output) if trace_exists else None,
        "trace_parse_error": trace_error,
        "phase_timings": phase_timings,
        "archive_economics": archive,
        "metadata": metadata,
        **result,
    }
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
