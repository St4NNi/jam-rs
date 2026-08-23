#!/usr/bin/env python3
"""Measure one trace command and extract its machine-readable evidence.

The process is intentionally launched by a fresh Python interpreter for each
measurement, so Linux child peak RSS is not contaminated by earlier runs.
Environment creation and compilation are never performed here.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import resource
import subprocess
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


def read_jsonl(path: Path) -> list[dict]:
    if path.suffix.lower() in {".zst", ".zstd"}:
        completed = subprocess.run(
            ["zstd", "--decompress", "--stdout", str(path)],
            check=True,
            stdout=subprocess.PIPE,
        )
        text = completed.stdout.decode("utf-8")
    else:
        text = path.read_text(encoding="utf-8")
    return [json.loads(line) for line in text.splitlines() if line.strip()]


def load_truth(path: Path | None) -> list[dict[str, str]]:
    if path is None:
        return []
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fetch_json(url: str | None) -> dict:
    if not url:
        return {}
    try:
        with urllib.request.urlopen(url, timeout=10) as response:
            return json.loads(response.read().decode("utf-8"))
    except (OSError, urllib.error.URLError, ValueError):
        return {}


def trace_metrics(path: Path | None, truth_path: Path | None) -> dict:
    if path is None or not path.exists():
        return {
            "trace_records": 0,
            "metagenome_count": 0,
            "candidate_count": 0,
            "alignment_count": 0,
            "coverage_mean": None,
            "gap_bases": 0,
            "largest_gap": 0,
            "resource_metrics": {},
        }
    try:
        records = read_jsonl(path)
    except (OSError, ValueError, subprocess.SubprocessError):
        records = []
    results = {
        record.get("metagenome_id"): record
        for record in records
        if record.get("record_type") == "metagenome_result"
    }
    footer = next(
        (record for record in records if record.get("record_type") == "run_footer"),
        {},
    )
    coverages = [
        result["coverage"].get("supported_fraction")
        for result in results.values()
        if result.get("coverage") is not None
        and result["coverage"].get("supported_fraction") is not None
    ]
    gap_bases = sum(
        gap.get("length", 0)
        for result in results.values()
        for gap in (result.get("coverage") or {}).get("gaps", [])
    )
    largest_gap = max(
        [
            (result.get("coverage") or {}).get("largest_gap", 0)
            for result in results.values()
        ]
        or [0]
    )
    metrics = {
        "trace_records": len(records),
        "metagenome_count": len(results),
        "candidate_count": sum(result.get("candidate") is not None for result in results.values()),
        "alignment_count": sum(len(result.get("alignments", [])) for result in results.values()),
        "coverage_mean": sum(coverages) / len(coverages) if coverages else None,
        "gap_bases": gap_bases,
        "largest_gap": largest_gap,
        "resource_metrics": footer.get("resource_metrics", {}),
    }
    truth = load_truth(truth_path)
    for class_name in ("exact", "reverse_complement", "partial", "split_overlap"):
        rows = [row for row in truth if row.get("truth_class") == class_name]
        expected = [row for row in rows if row.get("expected_candidate", "false").lower() == "true"]
        candidates = sum(results.get(row.get("metagenome_id"), {}).get("candidate") is not None for row in expected)
        aligned = sum(bool(results.get(row.get("metagenome_id"), {}).get("alignments")) for row in expected)
        metrics[f"{class_name}_candidate_recovery"] = f"{candidates}/{len(expected)}"
        metrics[f"{class_name}_alignment_recovery"] = f"{aligned}/{len(expected)}"
    short_rows = [row for row in truth if row.get("short_fragment", "false").lower() == "true"]
    short_expected = [row for row in short_rows if row.get("expected_candidate", "false").lower() == "true"]
    short_candidates = sum(results.get(row.get("metagenome_id"), {}).get("candidate") is not None for row in short_expected)
    short_aligned = sum(bool(results.get(row.get("metagenome_id"), {}).get("alignments")) for row in short_expected)
    metrics["short_fragment_candidate_recovery"] = f"{short_candidates}/{len(short_expected)}"
    metrics["short_fragment_alignment_recovery"] = f"{short_aligned}/{len(short_expected)}"
    return metrics


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--stdout", type=Path)
    parser.add_argument("--stderr", type=Path)
    parser.add_argument("--trace-output", type=Path)
    parser.add_argument("--truth", type=Path)
    parser.add_argument("--server-url")
    parser.add_argument("--stage", default="trace_search")
    parser.add_argument("--label", default="")
    parser.add_argument("--profile", default="")
    parser.add_argument("--cache-state", default="unspecified")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--extra-json", type=Path)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = list(args.command)
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        parser.error("a command is required after --")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    stdout_handle = args.stdout.open("w", encoding="utf-8") if args.stdout else None
    stderr_handle = args.stderr.open("w", encoding="utf-8") if args.stderr else None
    started_at = datetime.now(timezone.utc).isoformat()
    server_before = fetch_json(args.server_url)
    start = time.perf_counter()
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    try:
        completed = subprocess.run(command, stdout=stdout_handle, stderr=stderr_handle, check=False)
    finally:
        if stdout_handle:
            stdout_handle.close()
        if stderr_handle:
            stderr_handle.close()
    elapsed = time.perf_counter() - start
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    trace = trace_metrics(args.trace_output, args.truth)
    resource_metrics = trace.pop("resource_metrics", {})
    server_after = fetch_json(args.server_url)
    server_metrics = {
        key: value - server_before.get(key, 0)
        if isinstance(value, (int, float)) and isinstance(server_before.get(key, 0), (int, float))
        else value
        for key, value in server_after.items()
    }
    extra = {}
    if args.extra_json and args.extra_json.exists():
        extra = json.loads(args.extra_json.read_text(encoding="utf-8"))
    output = {
        "schema_version": "1.0.0",
        "stage": args.stage,
        "label": args.label,
        "profile": args.profile,
        "cache_state": args.cache_state,
        "threads": args.threads,
        "command": command,
        "cwd": os.getcwd(),
        "started_at_utc": started_at,
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
        "exit_code": completed.returncode,
        "wall_seconds": elapsed,
        "user_seconds": after.ru_utime - before.ru_utime,
        "system_seconds": after.ru_stime - before.ru_stime,
        "peak_rss_kib": after.ru_maxrss,
        "resource_metrics": resource_metrics,
        "server_metrics": server_metrics,
        **trace,
        "extra": extra,
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
