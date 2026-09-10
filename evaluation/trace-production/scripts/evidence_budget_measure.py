#!/usr/bin/env python3
"""Run a frozen, rotated small-search plan and record complete invocation cost."""

import argparse
import hashlib
import json
import math
import os
import re
import runpy
import subprocess
import time
from pathlib import Path


FORMAT = "jam-evidence-budget-measure-v1"
TIME_FIELDS = {
    "User time (seconds)": ("user_seconds", float),
    "System time (seconds)": ("system_seconds", float),
    "Percent of CPU this job got": ("cpu_percent", lambda value: float(value.rstrip("%"))),
    "Maximum resident set size (kbytes)": ("max_rss_kib", int),
    "Major (requiring I/O) page faults": ("major_faults", int),
    "Minor (reclaiming a frame) page faults": ("minor_faults", int),
    "File system inputs": ("filesystem_inputs", int),
    "File system outputs": ("filesystem_outputs", int),
    "Elapsed (wall clock) time (h:mm:ss or m:ss)": ("gnu_elapsed_seconds", None),
}


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def save(path: Path, value) -> None:
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())


def parse_time(path: Path) -> dict:
    values = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        for label, (name, cast) in TIME_FIELDS.items():
            prefix = label + ": "
            if stripped.startswith(prefix):
                raw = stripped.removeprefix(prefix)
                if name == "gnu_elapsed_seconds":
                    parts = [float(value) for value in raw.split(":")]
                    values[name] = sum(value * 60**index for index, value in enumerate(reversed(parts)))
                else:
                    values[name] = cast(raw)
    if set(values) != {item[0] for item in TIME_FIELDS.values()}:
        raise ValueError("GNU time output is incomplete")
    return values


def validate_plan(path: Path, expected: str, timeout: int) -> dict:
    if digest(path) != expected:
        raise ValueError("measurement plan identity differs")
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan.get("format") != FORMAT or not plan.get("runs"):
        raise ValueError("invalid measurement plan")
    labels = [row.get("label") for row in plan["runs"]]
    measurements = [row.get("measurement") for row in plan["runs"]]
    if (len(labels) != len(set(labels))
            or any(not re.fullmatch(r"[A-Za-z0-9_.-]+", value or "")
                   for value in labels + measurements)):
        raise ValueError("measurement labels must be unique safe names")
    required = ("label", "measurement", "method", "variant", "phase", "repetition",
                "topology", "binary", "binary_sha256", "output_kind",
                "output_suffix", "expected_query_records", "command")
    for row in plan["runs"]:
        if any(field not in row for field in required):
            raise ValueError("incomplete measurement row")
        if (row["topology"] not in ("linear", "circular", "native")
                or row["output_kind"] not in ("jam_jsonl", "opaque_native")
                or not re.fullmatch(r"\.[a-z0-9.]+", row["output_suffix"])
                or row["command"].count("{output}") != 1):
            raise ValueError("each command needs one output placeholder and explicit topology")
        binary = Path(row["binary"])
        if not binary.is_file() or binary.is_symlink() or digest(binary) != row["binary_sha256"]:
            raise ValueError(f"binary identity differs: {row['label']}")
        if not row["command"] or Path(row["command"][0]) != binary:
            raise ValueError(f"command does not invoke its bound binary: {row['label']}")
        if row["phase"] not in ("first-observed", "warmup", "measured"):
            raise ValueError("invalid phase or expected query count")
        expected_queries = row["expected_query_records"]
        if ((row["output_kind"] == "jam_jsonl"
             and (not isinstance(expected_queries, int) or expected_queries < 1))
                or (row["output_kind"] == "opaque_native" and expected_queries is not None)):
            raise ValueError("output kind and expected query count differ")
        if timeout > 900:
            raise ValueError("exploratory run timeout cannot exceed 15 minutes")
    for binding in plan.get("bindings", []):
        path = Path(binding["path"])
        if not path.is_file() or path.is_symlink() or digest(path) != binding["sha256"]:
            raise ValueError(f"input binding differs: {path}")
    return plan


def merge_outputs(output: Path, plan: dict, records: list[dict]) -> list[dict]:
    merged_dir = output / "merged"
    merged_dir.mkdir()
    summaries = []
    groups = {}
    by_label = {row["label"]: row for row in records}
    for planned in plan["runs"]:
        groups.setdefault(planned["measurement"], []).append(planned)
    for measurement, rows in groups.items():
        identity = {(row["method"], row["variant"], row["phase"], row["repetition"]) for row in rows}
        if (len(rows) not in (1, 2) or len(identity) != 1
                or len({row["output_kind"] for row in rows}) != 1
                or (len(rows) == 2 and {row["topology"] for row in rows} != {"linear", "circular"})
                or (len(rows) == 1 and rows[0]["output_kind"] == "jam_jsonl"
                    and rows[0]["expected_query_records"] != 1)):
            raise ValueError(f"measurement topology grouping differs: {measurement}")
        seen = set()
        merged = None
        if rows[0]["output_kind"] == "jam_jsonl":
            merged = merged_dir / f"{measurement}.jsonl"
            with merged.open("x", encoding="utf-8") as stream:
                for row in rows:
                    source = output / "detail" / row["label"] / f"results{row['output_suffix']}"
                    for line in source.read_text(encoding="utf-8").splitlines():
                        item = json.loads(line)
                        if item["query_id"] in seen:
                            raise ValueError(f"duplicate merged query: {item['query_id']}")
                        seen.add(item["query_id"])
                        stream.write(json.dumps(item, sort_keys=True) + "\n")
        observed = [by_label[row["label"]] for row in rows]
        summaries.append({"measurement": measurement, "method": rows[0]["method"],
                          "variant": rows[0]["variant"], "phase": rows[0]["phase"],
                          "repetition": rows[0]["repetition"], "output_kind": rows[0]["output_kind"],
                          "query_records": len(seen) if merged else None,
                          "complete_batch_wall_seconds": sum(item["wall_seconds"] for item in observed),
                          "complete_batch_gnu_elapsed_seconds": sum(
                              item["gnu_elapsed_seconds"] for item in observed),
                          "user_seconds": sum(item["resources"]["user_seconds"] for item in observed),
                          "system_seconds": sum(item["resources"]["system_seconds"] for item in observed),
                          "major_faults": sum(item["resources"]["major_faults"] for item in observed),
                          "process_rss_max_kib": max(item["resources"]["max_rss_kib"] for item in observed),
                          "process_memory_headroom_bytes": plan["process_memory_limit_bytes"]
                          - 1024 * max(item["resources"]["max_rss_kib"] for item in observed),
                          "native_outputs": [item["output"] for item in observed],
                          "merged_results": str(merged) if merged else None,
                          "merged_results_sha256": digest(merged) if merged else None})
        if summaries[-1]["complete_batch_wall_seconds"] >= 300:
            raise ValueError(f"small complete measurement exceeded five minutes: {measurement}")
    return summaries


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--expected-plan-sha256", required=True)
    parser.add_argument("--gnu-time", type=Path, default=Path("/usr/bin/time"))
    parser.add_argument("--timeout-command", type=Path, default=Path("/usr/bin/timeout"))
    parser.add_argument("--cgroup-helper", type=Path)
    parser.add_argument("--expected-cgroup-helper-sha256")
    parser.add_argument("--timeout-seconds", type=int, default=270)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (args.output.exists() or not args.gnu_time.is_file() or not args.timeout_command.is_file()
            or not 1 <= args.timeout_seconds <= 900):
        raise SystemExit("new output, GNU time, and a 1..900 second timeout are required")
    plan = validate_plan(args.plan, args.expected_plan_sha256, args.timeout_seconds)
    process_limit = plan.get("process_memory_limit_bytes")
    if not isinstance(process_limit, int) or process_limit <= 0:
        raise SystemExit("plan lacks a positive process memory limit")
    if bool(args.cgroup_helper) != bool(args.expected_cgroup_helper_sha256):
        raise SystemExit("cgroup helper path and identity must be supplied together")
    cgroup = None
    before = None
    if args.cgroup_helper:
        if digest(args.cgroup_helper) != args.expected_cgroup_helper_sha256:
            raise SystemExit("cgroup helper identity differs")
        cgroup = runpy.run_path(str(args.cgroup_helper))
        before = cgroup["snapshot"](os.getpid())
        expected_memory = plan.get("job_memory_limit_bytes")
        if (not before["job_cgroup_found"] or expected_memory is None
                or before["groups"][-1].get("memory.max") != expected_memory):
            raise SystemExit("Slurm job cgroup or declared memory limit differs")
    elif os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Slurm measurements require bound cgroup accounting")
    args.output.mkdir(mode=0o700)
    detail = args.output / "detail"
    detail.mkdir()
    records = []
    all_started = time.perf_counter()
    for row in plan["runs"]:
        elapsed = time.perf_counter() - all_started
        effective_timeout = min(args.timeout_seconds, math.floor(900 - elapsed - 5))
        if effective_timeout < 1:
            raise SystemExit("measurement plan exceeded the fifteen-minute exploration boundary")
        run_dir = detail / row["label"]
        run_dir.mkdir()
        result_path = run_dir / f"results{row['output_suffix']}"
        command = [str(result_path) if value == "{output}" else value for value in row["command"]]
        time_path, stdout_path, stderr_path = run_dir / "time.txt", run_dir / "stdout", run_dir / "stderr"
        timed = [str(args.timeout_command), "--signal=TERM", "--kill-after=5",
                 f"{effective_timeout}s", str(args.gnu_time), "--verbose", "--output",
                 str(time_path), "--"] + command
        started = time.perf_counter()
        with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr:
            result = subprocess.run(timed, stdout=stdout, stderr=stderr, check=False)
            returncode = result.returncode
            timed_out = returncode in (124, 137)
        wall = time.perf_counter() - started
        if timed_out or returncode or not result_path.is_file():
            save(args.output / "FAILURE.json", {"status": "failed", "label": row["label"],
                                                "timed_out": timed_out,
                                                "returncode": returncode,
                                                "wall_seconds": wall})
            raise SystemExit(f"measurement failed: {row['label']}")
        resources = parse_time(time_path)
        record = {field: row[field] for field in ("label", "measurement", "method", "variant",
                                                   "phase", "repetition", "topology")}
        record.update({"command": command, "wall_seconds": wall,
                       "effective_timeout_seconds": effective_timeout,
                       "gnu_elapsed_seconds": resources["gnu_elapsed_seconds"], "resources": resources,
                       "output": str(result_path), "output_bytes": result_path.stat().st_size,
                       "results_sha256": digest(result_path),
                       "result_records": (sum(1 for line in result_path.read_bytes().splitlines() if line)
                                          if row["output_kind"] == "jam_jsonl" else None)})
        if not result_path.stat().st_size:
            raise ValueError(f"native output is empty: {row['label']}")
        if (row["output_kind"] == "jam_jsonl"
                and record["result_records"] != row["expected_query_records"]):
            raise ValueError(f"result query count differs: {row['label']}")
        if resources["max_rss_kib"] * 1024 >= process_limit:
            raise ValueError(f"process memory limit crossed: {row['label']}")
        save(run_dir / "record.json", record)
        records.append(record)
    summaries = merge_outputs(args.output, plan, records)
    cgroup_summary = None
    if cgroup is not None:
        after = cgroup["snapshot"](os.getpid())
        delta = cgroup["subtract"](before, after)
        events = delta["groups"][-1].get("memory.events_delta", {})
        if events.get("status") != "captured" or any(
                events.get("values", {}).get(name, 0) for name in ("max", "oom", "oom_kill", "oom_group_kill")):
            raise ValueError("job memory boundary was crossed or unavailable")
        save(args.output / "cgroup-before.json", before)
        save(args.output / "cgroup-after.json", after)
        save(args.output / "cgroup-delta.json", delta)
        cgroup_summary = {"helper_sha256": digest(args.cgroup_helper),
                          "job_memory_limit_bytes": after["groups"][-1]["memory.max"],
                          "job_memory_peak_bytes": after["groups"][-1].get("memory.peak"),
                          "memory_events_delta": events["values"]}
    save(args.output / "runs.json", records)
    summary = {"status": "complete", "plan_sha256": digest(args.plan),
               "cache_policy": plan.get("cache_policy"),
               "schedule_policy": plan.get("schedule_policy"),
               "measurements": summaries,
               "job_cgroup": cgroup_summary,
               "timing_scope": "sum of complete linear and circular native process invocations",
               "physical_storage_io_status": "GNU time filesystem counters are kernel accounting, not storage bytes"}
    save(args.output / "summary.json", summary)
    save(args.output / "COMPLETE.json", {"status": "complete", "summary_sha256": digest(args.output / "summary.json"),
                                         "runs_sha256": digest(args.output / "runs.json")})


if __name__ == "__main__":
    main()
