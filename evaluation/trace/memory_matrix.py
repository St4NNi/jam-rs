#!/usr/bin/env python3
"""Run a bounded trace resource and memory matrix.

This is a measurement harness, not a workload model or a result
interpreter.  It launches an already-built ``jam trace`` binary through the
existing ``measure.py`` wrapper and stores the raw process and resource
counters alongside a machine-readable matrix record.  The input manifest is
deliberately explicit so a run can be reproduced without downloading data.

Example manifest::

    {
      "database": "tools/out/dataset/metagenomes.jam",
      "catalog": "tools/out/dataset/catalog.tsv",
      "queries": [
        {"id": "small", "class": "small", "path": "queries/small.fa"},
        {"id": "large", "class": "large", "path": "queries/large.fa"},
        {"id": "repeat-rich", "class": "repeat-rich", "path": "queries/repeat-rich.fa"}
      ],
      "archives": [
        {"path": "tools/out/dataset/assemblies/sample.fasta.jma"}
      ]
    }

Paths in the manifest are relative to the manifest file.  The default
output location is intentionally not implicit: callers must provide a new
output directory.  No existing output directory is reused or removed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[2]
PROFILES = ("fast", "balanced", "sensitive")
THREADS = (1, 4, 8, 16)
CANDIDATE_LIMITS = (1, 100)
CACHE_BLOCK_BYTES = (16 * 1024, 64 * 1024, 256 * 1024, 1024 * 1024)
QUERY_CLASSES = ("small", "large", "repeat-rich")

SUMMARY_FIELDS = (
    "schema_version",
    "query_id",
    "query_class",
    "profile",
    "threads",
    "candidate_limit",
    "cache_block_bytes",
    "repeat",
    "status",
    "exit_code",
    "wall_seconds",
    "user_seconds",
    "system_seconds",
    "peak_rss_kib",
    "candidate_count",
    "alignment_count",
    "output_bytes",
    "archive_count",
    "archive_size_bytes",
    "archive_build_wall_seconds",
    "archive_build_cpu_seconds",
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
    "seed_count_proxy",
    "bucket_count_proxy",
    "anchor_count_proxy",
    "chain_count_proxy",
    "sequence_count_proxy",
    "alignment_count_proxy",
    "output_bytes_proxy",
    "thread_count_proxy",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def command_output(command: list[str]) -> str | None:
    try:
        completed = subprocess.run(
            command,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
    except OSError:
        return None
    text = completed.stdout.strip()
    return text if completed.returncode == 0 else text or None


def resolve_path(value: str | os.PathLike[str], base: Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (base / path).resolve()


def load_manifest(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise SystemExit(f"matrix manifest is not a file: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise SystemExit(f"cannot read matrix manifest {path}: {error}") from error
    if not isinstance(value, dict):
        raise SystemExit("matrix manifest must contain a JSON object")
    return value


def _query_entries(manifest: dict[str, Any], base: Path) -> list[dict[str, Any]]:
    raw = manifest.get("queries")
    if raw is None and manifest.get("query") is not None:
        raw = [manifest["query"]]
    if not isinstance(raw, list) or not raw:
        raise SystemExit("matrix manifest requires a non-empty queries list")
    queries: list[dict[str, Any]] = []
    for index, item in enumerate(raw):
        if not isinstance(item, dict):
            raise SystemExit(f"query entry {index} must be an object")
        query_id = str(item.get("id") or item.get("name") or f"query-{index + 1}")
        query_class = str(item.get("class") or query_id).lower().replace("_", "-")
        query_path = item.get("path") or item.get("plasmid") or item.get("query")
        if not isinstance(query_path, str):
            raise SystemExit(f"query {query_id!r} has no path")
        path = resolve_path(query_path, base)
        if not path.is_file():
            raise SystemExit(f"query {query_id!r} is not a file: {path}")
        queries.append({"id": query_id, "class": query_class, "path": path})
    return queries


def _manifest_resources(manifest: dict[str, Any], base: Path) -> tuple[Path, Path]:
    database = manifest.get("database")
    catalog = manifest.get("catalog")
    if not isinstance(database, str) or not isinstance(catalog, str):
        raise SystemExit("matrix manifest requires database and catalog paths")
    database_path = resolve_path(database, base)
    catalog_path = resolve_path(catalog, base)
    if not database_path.is_file():
        raise SystemExit(f"candidate database is not a file: {database_path}")
    if not catalog_path.is_file():
        raise SystemExit(f"trace catalog is not a file: {catalog_path}")
    return database_path, catalog_path


def _int_value(value: Any) -> int | None:
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int):
        return value
    if isinstance(value, float) and value.is_integer():
        return int(value)
    return None


def _number(value: Any) -> int | float | None:
    if isinstance(value, bool):
        return int(value)
    return value if isinstance(value, (int, float)) else None


def _nested_values(value: Any, wanted: set[str]) -> list[int | float]:
    values: list[int | float] = []
    if isinstance(value, dict):
        for key, child in value.items():
            normalized = str(key).lower().replace("-", "_")
            if normalized in wanted:
                number = _number(child)
                if number is not None:
                    values.append(number)
            values.extend(_nested_values(child, wanted))
    elif isinstance(value, list):
        for child in value:
            values.extend(_nested_values(child, wanted))
    return values


def _proxy(records: list[dict[str, Any]], aliases: Iterable[str]) -> int | float | None:
    wanted = {alias.lower().replace("-", "_") for alias in aliases}
    values = _nested_values(records, wanted)
    if not values:
        return None
    return sum(values)


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
    values: list[dict[str, Any]] = []
    for line in text.splitlines():
        if line.strip():
            decoded = json.loads(line)
            if isinstance(decoded, dict):
                values.append(decoded)
    return values


def trace_summary(trace_path: Path) -> dict[str, Any]:
    records = read_jsonl(trace_path)
    results = [record for record in records if record.get("record_type") == "metagenome_result"]
    footer = next(
        (record for record in records if record.get("record_type") == "run_footer"),
        {},
    )
    resource_metrics = footer.get("resource_metrics")
    if not isinstance(resource_metrics, dict):
        resource_metrics = {}
    alignments = sum(
        len(record.get("alignments", []))
        for record in results
        if isinstance(record.get("alignments", []), list)
    )
    return {
        "candidate_count": sum(record.get("candidate") is not None for record in results),
        "alignment_count": alignments,
        "resource_metrics": resource_metrics,
        "records": records,
    }


def _catalog_archive_paths(catalog: Path) -> list[Path]:
    paths: list[Path] = []
    try:
        with catalog.open(encoding="utf-8", newline="") as handle:
            rows = csv.DictReader(handle, delimiter="\t")
            for row in rows:
                value = row.get("jma") or row.get("resource") or row.get("archive")
                if not value or "://" in value:
                    continue
                candidate = Path(value)
                if not candidate.is_absolute():
                    candidate = (catalog.parent / candidate).resolve()
                if candidate.is_file():
                    paths.append(candidate)
    except (OSError, csv.Error):
        return []
    return paths


def archive_summary(manifest: dict[str, Any], base: Path, catalog: Path) -> dict[str, Any]:
    entries = manifest.get("archives")
    paths: list[Path] = []
    measurement_paths: list[Path] = []
    if isinstance(entries, list):
        for entry in entries:
            if isinstance(entry, str):
                paths.append(resolve_path(entry, base))
            elif isinstance(entry, dict):
                value = entry.get("path") or entry.get("jma") or entry.get("archive")
                if isinstance(value, str):
                    paths.append(resolve_path(value, base))
                measurement = entry.get("measurement") or entry.get("build_measurement")
                if isinstance(measurement, str):
                    measurement_paths.append(resolve_path(measurement, base))
    if not paths:
        paths = _catalog_archive_paths(catalog)
    sizes = [path.stat().st_size for path in paths if path.is_file()]
    wall_values: list[float] = []
    cpu_values: list[float] = []
    measurements = manifest.get("archive_build_measurements")
    if isinstance(measurements, list):
        for item in measurements:
            candidate = resolve_path(item, base) if isinstance(item, str) else None
            if candidate:
                measurement_paths.append(candidate)
    for candidate in measurement_paths:
        if candidate.is_file():
            try:
                value = json.loads(candidate.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            if isinstance(value.get("wall_seconds"), (int, float)):
                wall_values.append(float(value["wall_seconds"]))
            cpu = float(value.get("user_seconds", 0)) + float(value.get("system_seconds", 0))
            cpu_values.append(cpu)
    return {
        "archive_count": len(sizes),
        "archive_size_bytes": sum(sizes) if sizes else None,
        "archive_build_wall_seconds": sum(wall_values) if wall_values else None,
        "archive_build_cpu_seconds": sum(cpu_values) if cpu_values else None,
    }


def collect_environment(jam: Path, manifest: Path) -> dict[str, Any]:
    stat = os.statvfs(manifest.parent)
    return {
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "python": platform.python_version(),
        "cpu_count": os.cpu_count(),
        "rustc": command_output(["rustc", "--version", "--verbose"]),
        "cargo": command_output(["cargo", "--version"]),
        "jam": command_output([str(jam), "--version"]),
        "source_commit": command_output(["git", "rev-parse", "HEAD"]),
        "storage": {
            "filesystem_block_bytes": stat.f_bsize,
            "available_bytes_at_start": stat.f_bavail * stat.f_frsize,
        },
        "timing_policy": (
            "compiler, installation, environment setup, and comparator setup are outside timed commands"
        ),
    }


def _resource_value(metrics: dict[str, Any], key: str) -> int | float | None:
    value = metrics.get(key)
    return _number(value)


def enrich_measurement(
    measurement: dict[str, Any],
    trace_path: Path,
    query: dict[str, Any],
    profile: str,
    threads: int,
    candidate_limit: int,
    block_bytes: int,
    repeat: int,
    archive: dict[str, Any],
) -> dict[str, Any]:
    trace = trace_summary(trace_path)
    metrics = trace["resource_metrics"]
    records = trace["records"]
    component_aliases = {
        "seed_count_proxy": ("seed_count", "seeds_examined", "query_seeds", "seed_queries"),
        "bucket_count_proxy": ("bucket_count", "buckets_examined", "seed_buckets_read"),
        "anchor_count_proxy": ("anchor_count", "anchors_generated", "anchors_considered"),
        "chain_count_proxy": ("chain_count", "chains_generated", "chains_considered"),
        "sequence_count_proxy": ("sequence_count", "sequence_windows", "windows_retrieved", "sequence_blocks_read"),
        "alignment_count_proxy": ("alignment_count", "alignments_attempted", "alignments_succeeded"),
    }
    proxies = {name: _proxy(records, aliases) for name, aliases in component_aliases.items()}
    output_bytes = trace_path.stat().st_size if trace_path.is_file() else None
    proxies["output_bytes_proxy"] = output_bytes
    proxies["thread_count_proxy"] = threads
    row: dict[str, Any] = {
        "schema_version": "1.0.0",
        "query_id": query["id"],
        "query_class": query["class"],
        "profile": profile,
        "threads": threads,
        "candidate_limit": candidate_limit,
        "cache_block_bytes": block_bytes,
        "repeat": repeat,
        "status": "complete" if measurement.get("exit_code") == 0 else "failed",
        "exit_code": measurement.get("exit_code"),
        "wall_seconds": measurement.get("wall_seconds"),
        "user_seconds": measurement.get("user_seconds"),
        "system_seconds": measurement.get("system_seconds"),
        "peak_rss_kib": measurement.get("peak_rss_kib"),
        "candidate_count": trace["candidate_count"],
        "alignment_count": trace["alignment_count"],
        "output_bytes": output_bytes,
        **archive,
    }
    for key in (
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
    ):
        row[key] = _resource_value(metrics, key)
    row.update(proxies)
    row["raw_measurement"] = measurement
    row["resource_metrics"] = metrics
    return row


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(SUMMARY_FIELDS) + "\n")
        for row in rows:
            handle.write(
                "\t".join(
                    "" if row.get(field) is None else str(row.get(field))
                    for field in SUMMARY_FIELDS
                )
                + "\n"
            )


def parse_int_list(value: str, name: str) -> tuple[int, ...]:
    try:
        values = tuple(int(part) for part in value.split(",") if part.strip())
    except ValueError as error:
        raise SystemExit(f"{name} must be a comma-separated integer list") from error
    if not values or any(item <= 0 for item in values):
        raise SystemExit(f"{name} must contain positive integers")
    return values


def parse_profile_list(value: str) -> tuple[str, ...]:
    values = tuple(part.strip() for part in value.split(",") if part.strip())
    unknown = sorted(set(values) - set(PROFILES))
    if not values or unknown:
        raise SystemExit(f"profiles must be drawn from {', '.join(PROFILES)}")
    return values


def run_matrix(args: argparse.Namespace) -> int:
    manifest_path = args.manifest.resolve()
    manifest = load_manifest(manifest_path)
    base = manifest_path.parent
    database, catalog = _manifest_resources(manifest, base)
    queries = _query_entries(manifest, base)
    if not args.allow_missing_query_classes:
        observed = {query["class"] for query in queries}
        missing = [query_class for query_class in QUERY_CLASSES if query_class not in observed]
        if missing:
            raise SystemExit(
                "matrix manifest is missing query classes: "
                + ", ".join(missing)
                + "; pass --allow-missing-query-classes only for a partial measurement"
            )

    output = args.output_dir.resolve()
    if output.exists():
        raise SystemExit(f"refusing to reuse existing matrix output directory: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir()
    raw_dir = output / "raw"
    records_dir = output / "records"
    raw_dir.mkdir()
    records_dir.mkdir()

    jam = args.jam.resolve()
    if not jam.is_file() or not os.access(jam, os.X_OK):
        raise SystemExit(f"jam binary is not executable: {jam}")
    environment = collect_environment(jam, manifest_path)
    (output / "environment.json").write_text(
        json.dumps(environment, indent=2) + "\n", encoding="utf-8"
    )
    (output / "input_manifest.json").write_text(
        json.dumps(
            {
                "manifest_path": str(manifest_path),
                "manifest_sha256": sha256(manifest_path),
                "database": str(database),
                "catalog": str(catalog),
                "queries": [
                    {"id": item["id"], "class": item["class"], "path": str(item["path"])}
                    for item in queries
                ],
                "threads": list(args.threads),
                "profiles": list(args.profiles),
                "candidate_limits": list(args.candidates),
                "cache_block_bytes": list(args.cache_blocks),
                "repeats": args.repeats,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    archive = archive_summary(manifest, base, catalog)
    measure_script = Path(__file__).with_name("measure.py")
    rows: list[dict[str, Any]] = []
    failed = 0
    case_count = 0
    for query in queries:
        for profile in args.profiles:
            for threads in args.threads:
                for candidate_limit in args.candidates:
                    for block_bytes in args.cache_blocks:
                        for repeat in range(1, args.repeats + 1):
                            case_count += 1
                            stem = (
                                f"{query['id']}__{query['class']}__{profile}__t{threads}"
                                f"__c{candidate_limit}__b{block_bytes}__r{repeat}"
                            )
                            case = raw_dir / stem
                            case.mkdir()
                            trace_output = case / "trace.jsonl"
                            measurement_path = case / "measurement.json"
                            stdout_path = case / "stdout.log"
                            stderr_path = case / "stderr.log"
                            cache_dir = case / "cache"
                            cache_dir.mkdir()
                            command = [
                                str(jam),
                                "--threads",
                                str(threads),
                                "--memory-target",
                                str(args.memory_target),
                                "--silent",
                                "trace",
                                "--plasmid",
                                str(query["path"]),
                                "--database",
                                str(database),
                                "--catalog",
                                str(catalog),
                                "--output",
                                str(trace_output),
                                "--plasmid-id",
                                str(query["id"]),
                                "--sensitivity",
                                profile,
                                "--top-candidates",
                                str(candidate_limit),
                                "--cache-dir",
                                str(cache_dir),
                                "--cache-block-bytes",
                                str(block_bytes),
                            ]
                            measured = subprocess.run(
                                [
                                    sys.executable,
                                    str(measure_script),
                                    "--output",
                                    str(measurement_path),
                                    "--stdout",
                                    str(stdout_path),
                                    "--stderr",
                                    str(stderr_path),
                                    "--trace-output",
                                    str(trace_output),
                                    "--stage",
                                    "trace_memory_matrix",
                                    "--label",
                                    stem,
                                    "--profile",
                                    profile,
                                    "--cache-state",
                                    "matrix-new-cache",
                                    "--threads",
                                    str(threads),
                                    "--",
                                    *command,
                                ],
                                check=False,
                            )
                            if not measurement_path.is_file():
                                failed += 1
                                continue
                            measurement = json.loads(measurement_path.read_text(encoding="utf-8"))
                            row = enrich_measurement(
                                measurement,
                                trace_output,
                                query,
                                profile,
                                threads,
                                candidate_limit,
                                block_bytes,
                                repeat,
                                archive,
                            )
                            row["command"] = command
                            row["case_directory"] = str(case)
                            (records_dir / f"{stem}.json").write_text(
                                json.dumps(row, indent=2) + "\n", encoding="utf-8"
                            )
                            rows.append(row)
                            failed += int(measured.returncode != 0)
                            del measured
    (output / "summary.json").write_text(
        json.dumps(
            {
                "schema_version": "1.0.0",
                "case_count": case_count,
                "record_count": len(rows),
                "failed_count": failed,
                "rows": rows,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    write_tsv(output / "summary.tsv", rows)
    checksums: dict[str, str] = {}
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name != "checksums.json":
            checksums[str(path.relative_to(output))] = sha256(path)
    (output / "checksums.json").write_text(
        json.dumps(checksums, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "output": str(output),
                "case_count": case_count,
                "record_count": len(rows),
                "failed_count": failed,
            }
        )
    )
    return 1 if failed else 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--jam", type=Path, default=ROOT / "target" / "release" / "jam")
    parser.add_argument("--threads", type=lambda value: parse_int_list(value, "threads"), default=THREADS)
    parser.add_argument("--profiles", type=parse_profile_list, default=PROFILES)
    parser.add_argument("--candidates", type=lambda value: parse_int_list(value, "candidates"), default=CANDIDATE_LIMITS)
    parser.add_argument(
        "--cache-blocks",
        type=lambda value: parse_int_list(value, "cache-blocks"),
        default=CACHE_BLOCK_BYTES,
        help="comma-separated cache block sizes in bytes",
    )
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--memory-target", type=int, default=4)
    parser.add_argument(
        "--allow-missing-query-classes",
        action="store_true",
        help="permit a partial manifest without small, large, and repeat-rich queries",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.repeats <= 0 or args.memory_target <= 0:
        raise SystemExit("repeats and memory-target must be positive")
    return run_matrix(args)


if __name__ == "__main__":
    raise SystemExit(main())
