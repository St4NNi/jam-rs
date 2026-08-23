#!/usr/bin/env python3
"""Run reproducible production trace benchmark cases.

The benchmark is manifest-driven.  No data are downloaded or synthesized by
this script.  Missing track resources, optional remote resources, and
unsupported command parameters become explicit ``unavailable`` raw records;
they are never turned into negative biological results.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import random
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import parse_qsl, urlencode, urlsplit, urlunsplit


ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
SCHEMA_VERSION = "1.0.0"
CPU_VALUES = (1, 4, 8, 16)
IO_VALUES = (1, 8, 32)
DEFAULT_SWEEP: dict[str, list[Any]] = {
    "candidate_cutoff": [0.0, 0.01, 0.05],
    "shared_hashes": [3, 5, 10],
    "seed_density": ["fast", "balanced", "sensitive"],
    "rescue": ["profile-default"],
    "chain": ["profile-default"],
    "alignment": ["profile-default"],
    "identity": ["not-filtered"],
    "band": ["profile-default"],
}
SUPPORTED_DEFAULT = {"candidate_cutoff", "shared_hashes", "seed_density"}
PLACEHOLDERS = {
    "jam",
    "cpu_threads",
    "io_tasks",
    "query",
    "database",
    "catalog",
    "output",
    "query_id",
    "sensitivity",
    "candidate_cutoff",
    "shared_hashes",
    "seed_density",
    "rescue",
    "chain",
    "alignment",
    "identity",
    "band",
    "cache_dir",
    "cache_block_bytes",
    "top_candidates",
    "max_alignments",
}
SECRET_KEYS = {"authorization", "password", "secret", "signature", "token", "x-amz-credential", "x-amz-security-token", "x-amz-signature"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def command_output(command: list[str]) -> str | None:
    try:
        completed = subprocess.run(command, check=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    except OSError:
        return None
    text = completed.stdout.strip()
    return text if completed.returncode == 0 else text or None


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise SystemExit(f"cannot read JSON manifest {path}: {error}") from error
    if not isinstance(value, dict):
        raise SystemExit(f"JSON manifest must contain an object: {path}")
    return value


def resolve(value: str | os.PathLike[str], base: Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (base / path).resolve()


def slug(value: Any) -> str:
    text = str(value).strip().lower()
    return "".join(char if char.isalnum() or char in "._-" else "_" for char in text)[:80] or "na"


def redact_text(value: str) -> str:
    value = re.sub(r"(https?://)([^/@\s]+):([^/@\s]+)@", r"\1<redacted>@", value)
    parsed = urlsplit(value)
    if parsed.scheme and parsed.netloc and parsed.query:
        pairs = [(key, "<redacted>" if key.lower() in SECRET_KEYS else item) for key, item in parse_qsl(parsed.query, keep_blank_values=True)]
        value = urlunsplit((parsed.scheme, parsed.netloc, parsed.path, urlencode(pairs), parsed.fragment))
    return value


def redact_value(value: Any) -> Any:
    if isinstance(value, str):
        return redact_text(value)
    if isinstance(value, list):
        return [redact_value(item) for item in value]
    if isinstance(value, dict):
        return {key: ("<redacted>" if str(key).lower() in SECRET_KEYS else redact_value(item)) for key, item in value.items()}
    return value


def query_entries(manifest: dict[str, Any], base: Path) -> dict[str, dict[str, Any]]:
    raw = manifest.get("queries", [])
    if not isinstance(raw, list):
        raise SystemExit("manifest queries must be a list")
    result = {}
    for index, item in enumerate(raw):
        if not isinstance(item, dict):
            raise SystemExit(f"query {index} must be an object")
        query_id = str(item.get("id") or item.get("name") or f"query-{index + 1}")
        path_value = item.get("path") or item.get("plasmid") or item.get("query")
        result[query_id] = {
            **item,
            "id": query_id,
            "path": resolve(path_value, base) if isinstance(path_value, str) else None,
            "class": str(item.get("class") or "unspecified"),
        }
    return result


def track_entries(manifest: dict[str, Any], queries: dict[str, dict[str, Any]], base: Path) -> list[dict[str, Any]]:
    raw_tracks = manifest.get("tracks")
    if raw_tracks is None:
        raw_tracks = [{"id": "A", "kind": "controlled_truth", "queries": list(queries)}]
    if not isinstance(raw_tracks, list) or not raw_tracks:
        raise SystemExit("manifest tracks must be a non-empty list")
    tracks: list[dict[str, Any]] = []
    for index, item in enumerate(raw_tracks):
        if not isinstance(item, dict):
            raise SystemExit(f"track {index} must be an object")
        track = dict(item)
        track["id"] = str(track.get("id") or chr(ord("A") + index))
        track["kind"] = str(track.get("kind") or "controlled_truth")
        selected = track.get("queries", list(queries))
        if not isinstance(selected, list):
            raise SystemExit(f"track {track['id']} queries must be a list")
        resolved_queries: list[dict[str, Any]] = []
        for entry in selected:
            query_id = str(entry.get("id")) if isinstance(entry, dict) else str(entry)
            if query_id in queries:
                resolved_queries.append(queries[query_id])
            elif isinstance(entry, dict):
                resolved_queries.append(
                    {
                        **entry,
                        "id": query_id,
                        "path": resolve(entry["path"], base) if isinstance(entry.get("path"), str) else None,
                        "class": str(entry.get("class") or "unspecified"),
                    }
                )
        track["queries"] = resolved_queries
        if isinstance(track.get("truth"), str):
            track["truth"] = resolve(track["truth"], base)
        tracks.append(track)
    return tracks


def resources_for(track: dict[str, Any], manifest: dict[str, Any], base: Path, mode: str) -> tuple[str | None, str | None, str | None]:
    root = manifest.get("resources", {})
    if not isinstance(root, dict):
        root = {}
    resource = track.get(mode) if isinstance(track.get(mode), dict) else {}
    if mode == "local":
        resource = {**root, **resource}
    else:
        remote = root.get("remote", {}) if isinstance(root.get("remote"), dict) else {}
        resource = {**remote, **resource}
    database = resource.get("database")
    catalog = resource.get("catalog")
    reason = None
    if not isinstance(database, str) or not isinstance(catalog, str):
        reason = f"{mode} database/catalog is not configured"
    elif not ("://" in database) and not resolve(database, base).is_file():
        reason = f"{mode} database is missing: {database}"
    elif not ("://" in catalog) and not resolve(catalog, base).is_file():
        reason = f"{mode} catalog is missing: {catalog}"
    database_value = database if isinstance(database, str) else None
    catalog_value = catalog if isinstance(catalog, str) else None
    if database_value and "://" not in database_value:
        database_value = str(resolve(database_value, base))
    if catalog_value and "://" not in catalog_value:
        catalog_value = str(resolve(catalog_value, base))
    return database_value, catalog_value, reason


def archive_economics(manifest: dict[str, Any], base: Path) -> dict[str, Any]:
    archive = manifest.get("archive_economics", manifest.get("archive", {}))
    if not isinstance(archive, dict):
        return {"status": "unavailable", "reason": "archive economics not configured"}
    result: dict[str, Any] = {"status": "available"}
    for key in (
        "archive_count",
        "archive_size_bytes",
        "archive_build_wall_seconds",
        "archive_build_cpu_seconds",
        "baseline_query_wall_seconds",
        "baseline_query_cpu_seconds",
    ):
        if isinstance(archive.get(key), (int, float)):
            result[key] = archive[key]
    paths = archive.get("paths", archive.get("archives", []))
    if isinstance(paths, list):
        sizes = []
        for item in paths:
            if isinstance(item, str):
                path = resolve(item, base)
                if path.is_file():
                    sizes.append(path.stat().st_size)
        if sizes:
            result["archive_count"] = len(sizes)
            result["archive_size_bytes"] = sum(sizes)
    if len(result) == 1:
        return {"status": "unavailable", "reason": "archive economics has no measured values"}
    if result.get("baseline_query_wall_seconds") and result.get("archive_build_wall_seconds") is not None:
        result["break_even_queries_wall"] = result["archive_build_wall_seconds"] / result["baseline_query_wall_seconds"]
    else:
        result["break_even_queries_wall"] = None
    if result.get("baseline_query_cpu_seconds") and result.get("archive_build_cpu_seconds") is not None:
        result["break_even_queries_cpu"] = result["archive_build_cpu_seconds"] / result["baseline_query_cpu_seconds"]
    else:
        result["break_even_queries_cpu"] = None
    return result


def sweep_cases(manifest: dict[str, Any]) -> list[tuple[str, Any]]:
    configured = manifest.get("sweep", {})
    configured = configured if isinstance(configured, dict) else {}
    sweep: dict[str, list[Any]] = {}
    for dimension, defaults in DEFAULT_SWEEP.items():
        values = configured.get(dimension, defaults)
        if not isinstance(values, list) or not values:
            values = defaults
        sweep[dimension] = values
    cases: list[tuple[str, Any]] = [("baseline", "baseline")]
    # One dimension at a time keeps the matrix interpretable and prevents a
    # combinatorial explosion.  A custom command template may opt in to all
    # dimensions; the raw record still names each requested value.
    for dimension, values in sweep.items():
        for value in values:
            cases.append((dimension, value))
    return cases


def command_template(manifest: dict[str, Any]) -> list[str] | None:
    value = manifest.get("command_template")
    if value is None:
        return None
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise SystemExit("command_template must be a list of strings")
    unknown = set()
    for item in value:
        for chunk in item.split("{")[1:]:
            if "}" in chunk:
                name = chunk.split("}", 1)[0]
                if name not in PLACEHOLDERS:
                    unknown.add(name)
    if unknown:
        raise SystemExit(f"command_template contains unknown placeholders: {sorted(unknown)}")
    return value


def command_for(template: list[str] | None, jam: Path, query: dict[str, Any], database: str, catalog: str, output: Path, cache: Path, cpu: int, io_tasks: int, params: dict[str, Any], memory_target: int, block_bytes: int, top_candidates: int, max_alignments: int) -> tuple[list[str], set[str]]:
    values = {
        "jam": str(jam),
        "cpu_threads": cpu,
        "io_tasks": io_tasks,
        "query": str(query["path"]),
        "database": database,
        "catalog": catalog,
        "output": str(output),
        "query_id": query["id"],
        "sensitivity": params.get("seed_density", "balanced") if isinstance(params.get("seed_density"), str) and params.get("seed_density") in {"fast", "balanced", "sensitive"} else "balanced",
        "candidate_cutoff": params.get("candidate_cutoff", 0.0),
        "shared_hashes": params.get("shared_hashes", 3),
        "seed_density": params.get("seed_density", "balanced"),
        "rescue": params.get("rescue", "profile-default"),
        "chain": params.get("chain", "profile-default"),
        "alignment": params.get("alignment", "profile-default"),
        "identity": params.get("identity", "not-filtered"),
        "band": params.get("band", "profile-default"),
        "cache_dir": str(cache),
        "cache_block_bytes": block_bytes,
        "top_candidates": top_candidates,
        "max_alignments": max_alignments,
    }
    if template is not None:
        rendered = [item.format(**values) for item in template]
        used = {name for name in PLACEHOLDERS if any("{" + name + "}" in item for item in template)}
        return rendered, used
    if query["path"] is None:
        return [], {"query"}
    return [
        str(jam),
        "--threads", str(cpu),
        "--memory-target", str(memory_target),
        "--silent", "trace",
        "--plasmid", str(query["path"]),
        "--database", database,
        "--catalog", catalog,
        "--output", str(output),
        "--plasmid-id", str(query["id"]),
        "--sensitivity", str(values["sensitivity"]),
        "--min-shared", str(values["shared_hashes"]),
        "--min-plasmid-containment", str(values["candidate_cutoff"]),
        "--min-metagenome-containment", str(values["candidate_cutoff"]),
        "--top-candidates", str(top_candidates),
        "--max-alignments", str(max_alignments),
        "--cache-dir", str(cache),
        "--cache-block-bytes", str(block_bytes),
    ], {"cpu_threads", "query", "database", "catalog", "output", "query_id", "sensitivity", "shared_hashes", "candidate_cutoff", "cache_dir"}


def environment(jam: Path, manifest: Path) -> dict[str, Any]:
    stat = os.statvfs(manifest.parent)
    return {
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "cpu_count": os.cpu_count(),
        "python": platform.python_version(),
        "rustc": command_output(["rustc", "--version", "--verbose"]),
        "cargo": command_output(["cargo", "--version"]),
        "jam": command_output([str(jam), "--version"]),
        "source_commit": command_output(["git", "rev-parse", "HEAD"]),
        "source_dirty": bool(command_output(["git", "status", "--porcelain"])),
        "release_flags": {"profile": "release", "opt_level": 3, "lto": True, "codegen_units": 1, "panic": "abort"},
        "storage": {"filesystem_block_bytes": stat.f_bsize, "available_bytes_at_start": stat.f_bavail * stat.f_frsize},
        "cache_policy": "cold means a fresh process and case cache without harness warmup; kernel cache eviction is not attempted; warm means one unmeasured invocation followed by timed repeats",
        "timing_policy": "build, installation, environment creation, and comparator setup are outside timed commands",
    }


def write_record_unavailable(path: Path, **fields: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    record = {
        "schema_version": SCHEMA_VERSION,
        "status": "unavailable",
        "unavailable_reason": fields.pop("reason", "prerequisite unavailable"),
        "exit_code": None,
        "wall_seconds": None,
        "user_seconds": None,
        "system_seconds": None,
        "peak_rss_kib": None,
        "candidate_count": None,
        "alignment_count": None,
        "candidate_recall": None,
        "final_recall": None,
        "false_candidate_count": None,
        "interval_precision": None,
        "interval_recall": None,
        "gap_boundary_error_bases": None,
        "resource_metrics": {},
        **redact_value(fields),
    }
    path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jam", type=Path, default=ROOT / "target" / "release" / "jam")
    parser.add_argument("--tracks", default="", help="comma-separated track IDs; default is all manifest tracks")
    parser.add_argument("--modes", default="local", help="comma-separated resource modes")
    parser.add_argument("--cpu", default=",".join(str(value) for value in CPU_VALUES))
    parser.add_argument("--io", default=",".join(str(value) for value in IO_VALUES))
    parser.add_argument("--orders", type=int, default=3)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--memory-target", type=int, default=4)
    parser.add_argument("--cache-block-bytes", type=int, default=1024 * 1024)
    parser.add_argument("--top-candidates", type=int, default=250)
    parser.add_argument("--max-alignments", type=int, default=256)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    manifest_path = args.manifest.resolve()
    manifest = load_json(manifest_path)
    if manifest.get("schema_version") not in {None, SCHEMA_VERSION}:
        raise SystemExit(f"unsupported production manifest schema: {manifest.get('schema_version')}")
    output = args.output.resolve()
    if output.exists():
        raise SystemExit(f"refusing to reuse existing output directory: {output}")
    if args.orders <= 0 or args.repeats <= 0:
        raise SystemExit("orders and repeats must be positive")
    try:
        cpus = tuple(int(value) for value in args.cpu.split(",") if value)
        ios = tuple(int(value) for value in args.io.split(",") if value)
    except ValueError as error:
        raise SystemExit("--cpu and --io require comma-separated integers") from error
    if not cpus or any(value <= 0 for value in cpus) or not ios or any(value <= 0 for value in ios):
        raise SystemExit("--cpu and --io require positive values")
    base = manifest_path.parent
    queries = query_entries(manifest, base)
    tracks = track_entries(manifest, queries, base)
    selected_tracks = set(filter(None, (item.strip() for item in args.tracks.split(","))))
    if selected_tracks:
        tracks = [track for track in tracks if track["id"] in selected_tracks]
    modes = tuple(item.strip() for item in args.modes.split(",") if item.strip())
    if not modes:
        raise SystemExit("--modes must not be empty")
    jam = args.jam.resolve()
    template = command_template(manifest)
    output.mkdir(parents=True)
    raw_dir = output / "raw"
    raw_dir.mkdir()
    (output / "manifest.json").write_text(json.dumps(redact_value({**manifest, "source_manifest": str(manifest_path), "source_manifest_sha256": sha256(manifest_path)}), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output / "environment.json").write_text(json.dumps(environment(jam, manifest_path), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    archive = archive_economics(manifest, base)
    case_count = 0
    unavailable_count = 0
    executed_count = 0
    failed_count = 0
    sweep = sweep_cases(manifest)
    for track in tracks:
        truth = track.get("truth") if isinstance(track.get("truth"), Path) else None
        for mode in modes:
            database, catalog, resource_reason = resources_for(track, manifest, base, mode)
            for order in range(args.orders):
                selected = list(track.get("queries", []))
                if order:
                    random.Random(0x4A414D + order).shuffle(selected)
                if not selected:
                    selected = [{"id": "", "path": None, "class": "unavailable"}]
                for query in selected:
                    for cpu in cpus:
                        for io_tasks in ios:
                            for dimension, value in sweep:
                                params = {dimension: value} if dimension != "baseline" else {}
                                # The default CLI has no independent knobs for
                                # several algorithm thresholds.  Mark those
                                # cases unavailable unless a custom command
                                # template explicitly consumes the placeholder.
                                command_for(template, jam, query, database or "", catalog or "", Path("output"), Path("cache"), cpu, io_tasks, params, args.memory_target, args.cache_block_bytes, args.top_candidates, args.max_alignments)
                                unsupported = dimension not in {"baseline"} and dimension not in SUPPORTED_DEFAULT
                                if template is not None and dimension != "baseline":
                                    placeholder = dimension
                                    unsupported = "{" + placeholder + "}" not in " ".join(template)
                                if template is None and io_tasks != 1:
                                    unsupported = True
                                if template is not None and io_tasks != 1 and "{io_tasks}" not in " ".join(template):
                                    unsupported = True
                                if template is None and catalog and "://" in catalog:
                                    unsupported = True
                                for cache_state, repeats in ((f"{mode}-process-cold", (0,)), (f"{mode}-warm", tuple(range(1, args.repeats + 1)))):
                                    for repeat in repeats:
                                        case_count += 1
                                        stem = f"{track['id']}__{mode}__o{order}__{slug(query.get('id'))}__c{cpu}__i{io_tasks}__{slug(dimension)}-{slug(value)}__{cache_state}__r{repeat}"
                                        case = raw_dir / stem
                                        trace = case / "trace.jsonl"
                                        measurement = case / "measurement.json"
                                        cache = case / "cache"
                                        cache.mkdir(parents=True, exist_ok=True)
                                        case_fields = {
                                            "track": track["id"],
                                            "track_kind": track["kind"],
                                            "query_id": query.get("id", ""),
                                            "resource_mode": mode,
                                            "cache_state": cache_state,
                                            "cpu_threads": cpu,
                                            "io_tasks": io_tasks,
                                            "order": order,
                                            "repeat": repeat,
                                            "sweep_dimension": dimension,
                                            "sweep_value": value,
                                            "parameter_status": "unsupported" if unsupported else "available",
                                            "archive_economics": archive,
                                        }
                                        reason = resource_reason
                                        if query.get("path") is None:
                                            reason = "query path is not configured"
                                        elif not query["path"].is_file():
                                            reason = f"query path is missing: {query['path']}"
                                        elif truth is not None and not truth.is_file():
                                            reason = f"truth manifest is missing: {truth}"
                                        elif unsupported:
                                            reason = f"parameter {dimension} is not exposed by the selected command"
                                        if reason or database is None or catalog is None:
                                            write_record_unavailable(measurement, reason=reason or "resource is unavailable", **case_fields)
                                            unavailable_count += 1
                                            continue
                                        command, _ = command_for(template, jam, query, database, catalog, trace, cache, cpu, io_tasks, params, args.memory_target, args.cache_block_bytes, args.top_candidates, args.max_alignments)
                                        if args.dry_run:
                                            write_record_unavailable(measurement, reason="dry run", command=command, **case_fields)
                                            unavailable_count += 1
                                            continue
                                        case.mkdir(parents=True, exist_ok=True)
                                        archive_path = case / "archive.json"
                                        archive_path.write_text(json.dumps(archive, indent=2, sort_keys=True) + "\n", encoding="utf-8")
                                        if cache_state.endswith("-warm"):
                                            warmup = case / "warmup.jsonl"
                                            warmup_command, _ = command_for(template, jam, query, database, catalog, warmup, cache, cpu, io_tasks, params, args.memory_target, args.cache_block_bytes, args.top_candidates, args.max_alignments)
                                            subprocess.run(warmup_command, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                                        measure_command = [
                                            sys.executable, str(SCRIPT_DIR / "measure.py"),
                                            "--output", str(measurement), "--trace-output", str(trace),
                                            "--stdout", str(case / "stdout.log"), "--stderr", str(case / "stderr.log"),
                                            "--track", str(track["id"]), "--query-id", str(query.get("id", "")),
                                            "--metagenome-set", str(track.get("metagenome_set", "")),
                                            "--resource-mode", mode, "--cache-state", cache_state,
                                            "--cpu-threads", str(cpu), "--io-tasks", str(io_tasks),
                                            "--order", str(order), "--repeat", str(repeat),
                                            "--sweep-dimension", dimension, "--sweep-value", value,
                                            "--archive-json", str(archive_path),
                                        ]
                                        if truth is not None:
                                            measure_command.extend(["--truth", str(truth)])
                                        measure_command.extend(["--", *command])
                                        completed = subprocess.run(measure_command, check=False)
                                        executed_count += 1
                                        if completed.returncode != 0:
                                            failed_count += 1
    checksums = {}
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name != "checksums.json":
            checksums[str(path.relative_to(output))] = sha256(path)
    (output / "checksums.json").write_text(json.dumps(checksums, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    run_record = {
        "schema_version": SCHEMA_VERSION,
        "status": "failed" if failed_count else ("complete" if executed_count or unavailable_count else "unavailable"),
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "manifest_sha256": sha256(manifest_path),
        "tracks": [track["id"] for track in tracks],
        "resource_modes": list(modes),
        "cpu_threads": list(cpus),
        "io_tasks": list(ios),
        "orders": args.orders,
        "repeats": args.repeats,
        "order_policy": "order 0 follows manifest order; later orders use deterministic Fisher-Yates shuffles seeded by 0x4A414D + order",
        "case_count": case_count,
        "executed_count": executed_count,
        "failed_count": failed_count,
        "unavailable_count": unavailable_count,
        "sweep_dimensions": [dimension for dimension, _ in sweep],
        "cache_policy": "local-process-cold is a fresh process/cache; local-warm has one unmeasured warmup; kernel cache eviction is not attempted",
        "output_policy": "raw records are append-free per case and summaries are regenerated from them",
    }
    (output / "run.json").write_text(json.dumps(run_record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(run_record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
