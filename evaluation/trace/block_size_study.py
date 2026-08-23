#!/usr/bin/env python3
"""Measure JMA archive and resource-cache block-size combinations.

The study keeps archive sequence-block construction and resource cache block
reads as separate dimensions.  It builds fresh archives for each requested
archive block size, then runs the same trace query against each archive set
with each requested cache block size.  Every output directory must be new;
the script never deletes or overwrites a previous study.

The manifest format is shared with ``memory_matrix.py`` and additionally
accepts an ``assemblies`` list.  Each assembly entry has ``id`` and ``path``
fields.  If it is omitted, local ``raw`` paths in the catalog are used.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from memory_matrix import (  # noqa: E402
    CACHE_BLOCK_BYTES,
    ROOT,
    _manifest_resources,
    _query_entries,
    collect_environment,
    enrich_measurement,
    load_manifest,
    resolve_path,
    sha256,
)


SUMMARY_FIELDS = (
    "schema_version",
    "archive_block_bases",
    "cache_block_bytes",
    "query_id",
    "query_class",
    "profile",
    "threads",
    "status",
    "build_exit_code",
    "trace_exit_code",
    "archive_count",
    "archive_size_bytes",
    "archive_build_wall_seconds",
    "archive_build_cpu_seconds",
    "archive_build_peak_rss_kib",
    "wall_seconds",
    "user_seconds",
    "system_seconds",
    "peak_rss_kib",
    "candidate_count",
    "alignment_count",
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


def _assembly_entries(manifest: dict[str, Any], base: Path, catalog: Path) -> list[dict[str, Any]]:
    entries = manifest.get("assemblies")
    if isinstance(entries, list) and entries:
        result: list[dict[str, Any]] = []
        for index, item in enumerate(entries):
            if not isinstance(item, dict):
                raise SystemExit(f"assembly entry {index} must be an object")
            assembly_id = str(item.get("id") or item.get("name") or f"assembly-{index + 1}")
            value = item.get("path") or item.get("assembly") or item.get("raw")
            if not isinstance(value, str):
                raise SystemExit(f"assembly {assembly_id!r} has no path")
            path = resolve_path(value, base)
            if not path.is_file():
                raise SystemExit(f"assembly {assembly_id!r} is not a file: {path}")
            result.append({"id": assembly_id, "path": path})
        return result

    result = []
    try:
        import csv

        with catalog.open(encoding="utf-8", newline="") as handle:
            for index, row in enumerate(csv.DictReader(handle, delimiter="\t")):
                value = row.get("raw") or row.get("assembly")
                if not value or "://" in value:
                    continue
                path = Path(value)
                if not path.is_absolute():
                    path = (catalog.parent / path).resolve()
                if not path.is_file():
                    raise SystemExit(f"catalog raw resource is not a file: {path}")
                result.append(
                    {
                        "id": row.get("metagenome_id") or row.get("id") or f"assembly-{index + 1}",
                        "path": path,
                    }
                )
    except (OSError, csv.Error) as error:
        raise SystemExit(f"cannot read local assembly paths from {catalog}: {error}") from error
    if not result:
        raise SystemExit("manifest requires assemblies or a catalog with local raw paths")
    return result


def _run_measure(
    measure_script: Path,
    output: Path,
    command: list[str],
    *,
    stage: str,
    label: str,
    profile: str,
    threads: int,
    trace_output: Path | None = None,
) -> dict[str, Any] | None:
    output.mkdir(parents=True, exist_ok=True)
    measurement = output / "measurement.json"
    stdout = output / "stdout.log"
    stderr = output / "stderr.log"
    args = [
        sys.executable,
        str(measure_script),
        "--output",
        str(measurement),
        "--stdout",
        str(stdout),
        "--stderr",
        str(stderr),
        "--stage",
        stage,
        "--label",
        label,
        "--profile",
        profile,
        "--threads",
        str(threads),
    ]
    if trace_output is not None:
        args.extend(["--trace-output", str(trace_output)])
    args.extend(["--", *command])
    subprocess.run(args, check=False)
    if not measurement.is_file():
        return None
    return json.loads(measurement.read_text(encoding="utf-8"))


def _write_catalog(path: Path, archives: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("metagenome_id\tjma\traw\n")
        for entry in archives:
            handle.write(
                f"{entry['id']}\t{entry['jma']}\t{entry['raw']}\n"
            )


def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(SUMMARY_FIELDS) + "\n")
        for row in rows:
            handle.write(
                "\t".join(
                    "" if row.get(field) is None else str(row.get(field))
                    for field in SUMMARY_FIELDS
                )
                + "\n"
            )


def _positive_list(value: str, name: str) -> tuple[int, ...]:
    try:
        values = tuple(int(item) for item in value.split(",") if item.strip())
    except ValueError as error:
        raise SystemExit(f"{name} must be a comma-separated integer list") from error
    if not values or any(item <= 0 for item in values):
        raise SystemExit(f"{name} must contain positive integers")
    return values


def run_study(args: argparse.Namespace) -> int:
    manifest_path = args.manifest.resolve()
    manifest = load_manifest(manifest_path)
    base = manifest_path.parent
    database, source_catalog = _manifest_resources(manifest, base)
    assemblies = _assembly_entries(manifest, base, source_catalog)
    queries = _query_entries(manifest, base)
    output = args.output_dir.resolve()
    if output.exists():
        raise SystemExit(f"refusing to reuse existing block study output directory: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir()
    raw_dir = output / "raw"
    raw_dir.mkdir()
    archive_root = output / "archives"
    archive_root.mkdir()
    records_dir = output / "records"
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
                "source_catalog": str(source_catalog),
                "assemblies": [
                    {"id": item["id"], "path": str(item["path"])} for item in assemblies
                ],
                "queries": [
                    {"id": item["id"], "class": item["class"], "path": str(item["path"])}
                    for item in queries
                ],
                "archive_block_bases": list(args.archive_blocks),
                "cache_block_bytes": list(args.cache_blocks),
                "threads": args.threads,
                "profile": args.profile,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    measure_script = SCRIPT_DIR / "measure.py"
    rows: list[dict[str, Any]] = []
    failed = 0
    for archive_block in args.archive_blocks:
        archive_dir = archive_root / f"archive_b{archive_block}"
        archive_dir.mkdir()
        built: list[dict[str, Any]] = []
        build_measurements: list[dict[str, Any]] = []
        for assembly in assemblies:
            archive_path = archive_dir / f"{assembly['id']}.jma"
            build_dir = raw_dir / f"archive_b{archive_block}" / f"build__{assembly['id']}"
            command = [
                str(jam),
                "--memory-target",
                str(args.memory_target),
                "--silent",
                "archive",
                "--input",
                str(assembly["path"]),
                "--output",
                str(archive_path),
                "--block-bases",
                str(archive_block),
            ]
            measurement = _run_measure(
                measure_script,
                build_dir,
                command,
                stage="jma_archive_build",
                label=f"archive_b{archive_block}__{assembly['id']}",
                profile="archive",
                threads=args.threads,
            )
            if measurement is None:
                failed += 1
                continue
            build_measurements.append(measurement)
            failed += int(measurement.get("exit_code") != 0)
            if measurement.get("exit_code") == 0 and archive_path.is_file():
                built.append(
                    {
                        "id": assembly["id"],
                        "jma": archive_path.resolve(),
                        "raw": assembly["path"].resolve(),
                    }
                )
        if not built:
            continue
        catalog = archive_dir / "catalog.tsv"
        _write_catalog(catalog, built)
        archive_metrics = {
            "archive_count": len(built),
            "archive_size_bytes": sum(item["jma"].stat().st_size for item in built),
            "archive_build_wall_seconds": sum(
                float(item.get("wall_seconds", 0)) for item in build_measurements
            ),
            "archive_build_cpu_seconds": sum(
                float(item.get("user_seconds", 0)) + float(item.get("system_seconds", 0))
                for item in build_measurements
            ),
            "archive_build_peak_rss_kib": max(
                (item.get("peak_rss_kib") for item in build_measurements),
                default=None,
            ),
        }
        for cache_block in args.cache_blocks:
            for query in queries:
                stem = f"archive_b{archive_block}__cache_b{cache_block}__{query['id']}"
                case = raw_dir / stem
                trace_output = case / "trace.jsonl"
                cache_dir = case / "cache"
                command = [
                    str(jam),
                    "--threads",
                    str(args.threads),
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
                    args.profile,
                    "--top-candidates",
                    str(args.candidates),
                    "--cache-dir",
                    str(cache_dir),
                    "--cache-block-bytes",
                    str(cache_block),
                ]
                measurement = _run_measure(
                    measure_script,
                    case,
                    command,
                    stage="trace_cache_block_study",
                    label=stem,
                    profile=args.profile,
                    threads=args.threads,
                    trace_output=trace_output,
                )
                if measurement is None:
                    failed += 1
                    continue
                failed += int(measurement.get("exit_code") != 0)
                row = enrich_measurement(
                    measurement,
                    trace_output,
                    query,
                    args.profile,
                    args.threads,
                    args.candidates,
                    cache_block,
                    1,
                    archive_metrics,
                )
                row.update(
                    {
                        "archive_block_bases": archive_block,
                        "build_exit_code": max(
                            (item.get("exit_code", 1) for item in build_measurements),
                            default=1,
                        ),
                        "trace_exit_code": measurement.get("exit_code"),
                        "archive_build_peak_rss_kib": archive_metrics["archive_build_peak_rss_kib"],
                        "command": command,
                        "case_directory": str(case),
                    }
                )
                (records_dir / f"{stem}.json").write_text(
                    json.dumps(row, indent=2) + "\n", encoding="utf-8"
                )
                rows.append(row)
    summary = {
        "schema_version": "1.0.0",
        "record_count": len(rows),
        "failed_count": failed,
        "rows": rows,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    _write_tsv(output / "summary.tsv", rows)
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
    parser.add_argument("--archive-blocks", type=lambda value: _positive_list(value, "archive-blocks"), default=CACHE_BLOCK_BYTES)
    parser.add_argument("--cache-blocks", type=lambda value: _positive_list(value, "cache-blocks"), default=CACHE_BLOCK_BYTES)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--profile", choices=("fast", "balanced", "sensitive"), default="balanced")
    parser.add_argument("--candidates", type=int, default=100)
    parser.add_argument("--memory-target", type=int, default=4)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.threads <= 0 or args.candidates <= 0 or args.memory_target <= 0:
        raise SystemExit("threads, candidates, and memory-target must be positive")
    return run_study(args)


if __name__ == "__main__":
    raise SystemExit(main())
