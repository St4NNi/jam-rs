#!/usr/bin/env python3
"""Run the read-derived assembly benchmark track.

Track B is deliberately manifest driven.  The harness does not simulate reads,
guess an assembler, or turn an absent output into a negative result.  A run is
usable only when the read manifest, abundance grid, pinned tools, checksums,
and coordinate extractor are all present.  ``--self-check`` validates those
contracts without executing any external command.

The coordinate extractor is an explicit part of the manifest.  It must report
the portion of the query plasmid that was actually assembled.  Jam coverage is
intersected with those reported plasmid intervals, so sequence that was not
assembled cannot improve the Track B recovery value.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "1.0.0"
UNAVAILABLE_EXIT = 3
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")
IMAGE_DIGEST_RE = re.compile(r"^sha256:[0-9a-fA-F]{64}$")
PLACEHOLDER_RE = re.compile(r"(?:<[^>]+>|TODO|CHANGEME|REPLACE[-_ ]WITH)", re.I)
QUERY_KINDS = {"plasmid", "phage"}
TOPOLOGIES = {"circular", "linear", "unknown"}
SECRET_RE = re.compile(
    r"(?i)((?:[?&]|\b)(?:x-amz-signature|x-amz-credential|x-amz-security-token|authorization|password|token|secret)=)[^&\s]+"
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_load(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read JSON manifest {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"manifest must contain a JSON object: {path}")
    return value


def issue_missing(missing: list[str], text: str) -> None:
    if text not in missing:
        missing.append(text)


def issue_error(errors: list[str], text: str) -> None:
    if text not in errors:
        errors.append(text)


def valid_text(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip()) and not PLACEHOLDER_RE.search(value)


def redact(value: Any) -> Any:
    if isinstance(value, str):
        return SECRET_RE.sub(r"\1REDACTED", value)
    if isinstance(value, list):
        return [redact(item) for item in value]
    if isinstance(value, dict):
        return {key: redact(item) for key, item in value.items()}
    return value


def validate_license(value: Any, label: str, errors: list[str]) -> None:
    if not isinstance(value, dict):
        issue_error(errors, f"{label}: license object is required")
        return
    for field in ("spdx_id", "source", "redistribution"):
        if not valid_text(value.get(field)):
            issue_error(errors, f"{label}: license.{field} must be a non-placeholder string")


def validate_sha(value: Any, label: str, errors: list[str]) -> None:
    if not isinstance(value, str) or not SHA256_RE.fullmatch(value):
        issue_error(errors, f"{label}: sha256 must be 64 hexadecimal characters")


def validate_tool(value: Any, label: str, missing: list[str], errors: list[str]) -> None:
    if not isinstance(value, dict):
        issue_missing(missing, f"{label}: pinned command, image, and version")
        return
    for field in ("name", "version", "image", "image_digest"):
        if not valid_text(value.get(field)):
            issue_missing(missing, f"{label}.{field}")
    digest = value.get("image_digest")
    if digest is not None and not IMAGE_DIGEST_RE.fullmatch(str(digest)):
        issue_error(errors, f"{label}.image_digest must be a sha256 image digest")
    version = value.get("version")
    if isinstance(version, str) and version.strip().lower() in {
        "latest",
        "unknown",
        "unpinned",
        "not supplied",
        "not available",
        "n/a",
        "na",
    }:
        issue_error(errors, f"{label}.version must be pinned")
    command = value.get("command")
    if not isinstance(command, list) or not command or not all(isinstance(part, str) and part for part in command):
        issue_missing(missing, f"{label}.command")


def resolve_path(raw: Any, base: Path) -> Path | None:
    if not isinstance(raw, str) or not raw.strip():
        return None
    path = Path(raw)
    return path if path.is_absolute() else (base / path).resolve()


def validate_file(
    raw_path: Any,
    declared_sha: Any,
    license_value: Any,
    label: str,
    base: Path,
    missing: list[str],
    errors: list[str],
    check_content: bool,
) -> None:
    validate_license(license_value, label, errors)
    if not isinstance(declared_sha, str) or not SHA256_RE.fullmatch(declared_sha):
        issue_error(errors, f"{label}: sha256 must be 64 hexadecimal characters")
    path = resolve_path(raw_path, base)
    if path is None:
        issue_missing(missing, f"{label}: input path")
        return
    if not path.is_file():
        issue_missing(missing, f"{label}: file not found ({path})")
        return
    if check_content and isinstance(declared_sha, str) and SHA256_RE.fullmatch(declared_sha):
        actual = sha256(path)
        if actual.lower() != declared_sha.lower():
            issue_missing(missing, f"{label}: checksum mismatch (declared {declared_sha}, actual {actual})")


def query_entries(
    manifest: dict[str, Any],
    manifest_path: Path,
    missing: list[str],
    errors: list[str],
    check_content: bool,
) -> list[dict[str, Any]]:
    """Normalize new query_elements and legacy plasmids into one representation."""
    new_entries = manifest.get("query_elements")
    legacy_entries = manifest.get("plasmids")
    if new_entries is not None and legacy_entries is not None:
        issue_error(errors, "provide query_elements or legacy plasmids, not both")
        return []
    legacy = new_entries is None
    entries = legacy_entries if legacy else new_entries
    if not isinstance(entries, list) or not entries:
        issue_missing(missing, "at least one query element")
        return []
    valid: list[dict[str, Any]] = []
    seen: set[str] = set()
    for index, entry in enumerate(entries):
        label = f"{'plasmids' if legacy else 'query_elements'}[{index}]"
        if not isinstance(entry, dict):
            issue_error(errors, f"{label}: query element object is required")
            continue
        query_id = entry.get("plasmid_id") if legacy else entry.get("query_id")
        if not valid_text(query_id):
            issue_error(errors, f"{label}: {'plasmid_id' if legacy else 'query_id'} is required")
            continue
        query_kind = entry.get("query_kind", "plasmid" if legacy else None)
        topology = entry.get("topology", "circular" if legacy else "unknown")
        if query_kind not in QUERY_KINDS:
            issue_error(errors, f"{label}.query_kind must be plasmid or phage")
        if topology not in TOPOLOGIES:
            issue_error(errors, f"{label}.topology must be circular, linear, or unknown")
        if query_id in seen:
            issue_error(errors, f"duplicate query_id: {query_id}")
            continue
        seen.add(query_id)
        validate_file(
            entry.get("fasta"),
            entry.get("sha256"),
            entry.get("license"),
            label,
            manifest_path.parent,
            missing,
            errors,
            check_content,
        )
        path = resolve_path(entry.get("fasta"), manifest_path.parent)
        if path is not None and path.is_file():
            valid.append({**entry, "query_id": query_id, "query_kind": query_kind, "topology": topology, "_fasta_path": path})
    return valid


def validate_read_manifest(
    read_manifest: dict[str, Any],
    read_path: Path,
    parent: dict[str, Any],
    missing: list[str],
    errors: list[str],
    check_content: bool,
) -> list[dict[str, Any]]:
    if read_manifest.get("schema_version") != SCHEMA_VERSION:
        issue_error(errors, "read manifest schema_version must be 1.0.0")
    if not valid_text(read_manifest.get("manifest_id")):
        issue_error(errors, "read manifest manifest_id is required")
    if read_manifest.get("manifest_id") != parent.get("read_manifest_id"):
        issue_error(errors, "read manifest id does not match the benchmark manifest")
    validate_license(read_manifest.get("license"), "read manifest", errors)
    declared_manifest_sha = parent.get("read_manifest_sha256")
    validate_sha(declared_manifest_sha, "read_manifest_sha256", errors)
    if read_path.is_file() and isinstance(declared_manifest_sha, str) and SHA256_RE.fullmatch(declared_manifest_sha):
        actual = sha256(read_path)
        if actual.lower() != declared_manifest_sha.lower():
            issue_missing(missing, f"read manifest checksum mismatch (declared {declared_manifest_sha}, actual {actual})")
    else:
        issue_missing(missing, f"read manifest file not found ({read_path})")

    read_sets = read_manifest.get("read_sets")
    if not isinstance(read_sets, list) or not read_sets:
        issue_missing(missing, "read manifest must contain at least one read_set")
        return []
    seen: set[str] = set()
    valid_sets: list[dict[str, Any]] = []
    for index, read_set in enumerate(read_sets):
        label = f"read_sets[{index}]"
        if not isinstance(read_set, dict):
            issue_error(errors, f"{label}: object is required")
            continue
        read_set_id = read_set.get("read_set_id")
        if not valid_text(read_set_id):
            issue_error(errors, f"{label}.read_set_id is required")
            continue
        if read_set_id in seen:
            issue_error(errors, f"duplicate read_set_id: {read_set_id}")
            continue
        seen.add(read_set_id)
        validate_license(read_set.get("license"), label, errors)
        r1_path = resolve_path(read_set.get("r1"), read_path.parent)
        validate_file(
            read_set.get("r1"),
            read_set.get("r1_sha256"),
            read_set.get("license"),
            f"{label}.r1",
            read_path.parent,
            missing,
            errors,
            check_content,
        )
        r2_path = resolve_path(read_set.get("r2"), read_path.parent)
        if read_set.get("r2") is not None:
            validate_file(
                read_set.get("r2"),
                read_set.get("r2_sha256"),
                read_set.get("license"),
                f"{label}.r2",
                read_path.parent,
                missing,
                errors,
                check_content,
            )
        if r1_path is not None and r1_path.is_file():
            valid_sets.append({**read_set, "_r1_path": r1_path, "_r2_path": r2_path})
    return valid_sets


def validate_manifest(
    manifest: dict[str, Any],
    manifest_path: Path,
    read_manifest: dict[str, Any],
    read_manifest_path: Path,
    check_content: bool,
) -> tuple[list[str], list[str], list[dict[str, Any]], list[dict[str, Any]]]:
    missing: list[str] = []
    errors: list[str] = []
    if manifest.get("schema_version") != SCHEMA_VERSION:
        issue_error(errors, "benchmark schema_version must be 1.0.0")
    if manifest.get("track") != "B":
        issue_error(errors, "benchmark track must be B")
    if not valid_text(manifest.get("dataset_id")):
        issue_error(errors, "dataset_id is required")
    validate_license(manifest.get("license"), "benchmark", errors)
    for field in ("read_manifest_id", "read_manifest_sha256", "abundance_grid", "mixer", "assembler", "coordinate_extractor", "jam"):
        if field not in manifest:
            issue_missing(missing, f"benchmark field: {field}")
    if "query_elements" not in manifest and "plasmids" not in manifest:
        issue_missing(missing, "query_elements (or legacy plasmids)")
    validate_read_manifest(read_manifest, read_manifest_path, manifest, missing, errors, check_content)
    abundance_grid = manifest.get("abundance_grid")
    labels: set[str] = set()
    valid_abundances: list[dict[str, Any]] = []
    if not isinstance(abundance_grid, list) or not abundance_grid:
        issue_missing(missing, "abundance_grid with at least one level")
    else:
        for index, level in enumerate(abundance_grid):
            label = f"abundance_grid[{index}]"
            if not isinstance(level, dict) or not valid_text(level.get("label")):
                issue_error(errors, f"{label}: label and fraction are required")
                continue
            fraction = level.get("fraction")
            if not isinstance(fraction, (int, float)) or isinstance(fraction, bool) or not 0 < fraction <= 1:
                issue_error(errors, f"{label}.fraction must be in (0, 1]")
            if level["label"] in labels:
                issue_error(errors, f"duplicate abundance label: {level['label']}")
            labels.add(level["label"])
            valid_abundances.append(level)

    valid_queries = query_entries(manifest, manifest_path, missing, errors, check_content)

    for field in ("mixer", "assembler", "coordinate_extractor", "jam"):
        validate_tool(manifest.get(field), field, missing, errors)
    return errors, missing, valid_abundances, valid_queries


def prepare_output(path: Path) -> None:
    if path.exists():
        raise ValueError(f"refusing to reuse existing output directory: {path}")
    path.mkdir(parents=True)


def render_command(command: list[str], values: dict[str, Any]) -> list[str]:
    rendered: list[str] = []
    for part in command:
        try:
            value = part.format(**{key: str(item) for key, item in values.items()})
        except KeyError as exc:
            raise ValueError(f"command contains unknown placeholder: {{{exc.args[0]}}}") from exc
        rendered.append(value)
    return rendered


def run_command(command: list[str], cwd: Path, log_path: Path, tool: dict[str, Any], stage: str) -> dict[str, Any]:
    started = utc_now()
    start = time.perf_counter()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with log_path.open("w", encoding="utf-8") as log:
            completed = subprocess.run(command, cwd=cwd, stdout=log, stderr=subprocess.STDOUT, check=False)
        error = None
    except OSError as exc:
        completed = None
        error = str(exc)
    return {
        "stage": stage,
        "tool": {"name": tool.get("name"), "version": tool.get("version"), "image": tool.get("image"), "image_digest": tool.get("image_digest")},
        "command": command,
        "started_at_utc": started,
        "finished_at_utc": utc_now(),
        "wall_seconds": time.perf_counter() - start,
        "exit_code": completed.returncode if completed is not None else None,
        "log": str(log_path),
        "error": error,
    }


def fasta_sequence(path: Path) -> str:
    sequence: list[str] = []
    with path.open(encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            sequence.append(line.upper())
    value = "".join(sequence)
    if not value or set(value) - set("ACGTN"):
        raise ValueError(f"FASTA has no valid assembled sequence: {path}")
    return value


def union_intervals(intervals: list[tuple[int, int]], length: int | None = None) -> list[tuple[int, int]]:
    cleaned: list[tuple[int, int]] = []
    for start, end in intervals:
        if length is not None:
            start = max(0, min(length, start))
            end = max(0, min(length, end))
        if end > start:
            cleaned.append((start, end))
    cleaned.sort()
    result: list[tuple[int, int]] = []
    for start, end in cleaned:
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(result[-1][1], end))
        else:
            result.append((start, end))
    return result


def intersection_length(left: list[tuple[int, int]], right: list[tuple[int, int]]) -> int:
    total = 0
    for left_start, left_end in left:
        for right_start, right_end in right:
            total += max(0, min(left_end, right_end) - max(left_start, right_start))
    return total


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        raise ValueError(f"trace output was not created: {path}")
    if path.suffix.lower() in {".zst", ".zstd"}:
        completed = subprocess.run(["zstd", "--decompress", "--stdout", str(path)], stdout=subprocess.PIPE, check=False)
        if completed.returncode != 0:
            raise ValueError(f"cannot decompress trace output: {path}")
        text = completed.stdout.decode("utf-8")
    else:
        text = path.read_text(encoding="utf-8")
    records: list[dict[str, Any]] = []
    for line_number, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        value = json.loads(line)
        if not isinstance(value, dict):
            raise ValueError(f"trace record {line_number} is not an object")
        records.append(value)
    return records


def coordinate_truth(
    path: Path,
    query_id: str,
    query_kind: str,
    topology: str,
    query_length: int,
    job_id: str,
) -> tuple[str, list[dict[str, Any]], list[tuple[int, int]]]:
    try:
        value = json_load(path)
    except ValueError as exc:
        raise ValueError(f"invalid coordinate extractor output: {exc}") from exc
    if value.get("schema_version") != SCHEMA_VERSION or value.get("job_id") != job_id:
        raise ValueError("coordinate output schema_version or job_id does not match")
    output_id = value.get("query_id", value.get("plasmid_id"))
    output_length = value.get("query_length", value.get("plasmid_length"))
    if output_id != query_id or output_length != query_length:
        raise ValueError("coordinate output query identity or length does not match manifest")
    if value.get("query_kind", query_kind) != query_kind or value.get("topology", topology) != topology:
        raise ValueError("coordinate output query_kind or topology does not match manifest")
    sequence = value.get("assembled_query_sequence", value.get("assembled_plasmid_sequence"))
    if not isinstance(sequence, str) or not sequence:
        raise ValueError("coordinate output must contain assembled_query_sequence")
    sequence = sequence.upper()
    declared_sha = value.get("assembled_query_sequence_sha256", value.get("assembled_plasmid_sequence_sha256"))
    actual_sha = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    if declared_sha != actual_sha:
        raise ValueError("coordinate output assembled query sequence checksum does not match sequence")
    raw_intervals = value.get("intervals")
    if not isinstance(raw_intervals, list) or not raw_intervals:
        raise ValueError("coordinate output must contain at least one query interval")
    intervals: list[tuple[int, int]] = []
    normalized: list[dict[str, Any]] = []
    for index, interval in enumerate(raw_intervals):
        if not isinstance(interval, dict):
            raise ValueError(f"coordinate interval {index} is not an object")
        start = interval.get("query_start", interval.get("plasmid_start"))
        end = interval.get("query_end", interval.get("plasmid_end"))
        if not isinstance(start, int) or not isinstance(end, int) or not 0 <= start < end <= query_length:
            raise ValueError(f"coordinate interval {index} is outside the query element")
        if not valid_text(interval.get("contig_id")):
            raise ValueError(f"coordinate interval {index} lacks contig_id")
        intervals.append((start, end))
        normalized.append(interval)
    union = union_intervals(intervals, query_length)
    if not union:
        raise ValueError("coordinate output has zero assembled query bases")
    return sequence, normalized, union


def trace_support(path: Path, query_length: int, assembled: list[tuple[int, int]]) -> dict[str, Any]:
    records = read_jsonl(path)
    all_intervals: list[tuple[int, int]] = []
    for record in records:
        if record.get("record_type") != "metagenome_result":
            continue
        coverage = record.get("coverage") or {}
        for interval in coverage.get("primary_intervals", []):
            if not isinstance(interval, dict):
                continue
            start, end = interval.get("start"), interval.get("end")
            if isinstance(start, int) and isinstance(end, int) and 0 <= start < end <= query_length:
                all_intervals.append((start, end))
    if not any(record.get("record_type") == "metagenome_result" for record in records):
        raise ValueError("trace output contains no metagenome_result record")
    supported = union_intervals(all_intervals, query_length)
    assembled_bases = sum(end - start for start, end in assembled)
    observed_bases = intersection_length(supported, assembled)
    return {
        "trace_records": len(records),
        "supported_intervals": [{"start": start, "end": end} for start, end in supported],
        "assembled_intervals": [{"start": start, "end": end} for start, end in assembled],
        "supported_bases_in_assembled_sequence": observed_bases,
        "supported_fraction_of_assembled_sequence": observed_bases / assembled_bases if assembled_bases else None,
        "supported_fraction_of_query": observed_bases / query_length if query_length else None,
        "assembled_query_bases": assembled_bases,
    }


def unavailable_files(output: Path, run: dict[str, Any], missing: list[str]) -> None:
    run["status"] = "unavailable"
    run["missing_prerequisites"] = sorted(set(missing))
    (output / "run.json").write_text(json.dumps(run, indent=2) + "\n", encoding="utf-8")
    result = {
        "schema_version": SCHEMA_VERSION,
        "track": "B",
        "record_type": "run_status",
        "run_id": run["run_id"],
        "status": "unavailable",
        "missing_prerequisites": sorted(set(missing)),
        "note": "No recovery or accuracy value is reported until all prerequisites are available.",
    }
    (output / "results.jsonl").write_text(json.dumps(result) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> int:
    manifest_path = args.manifest.resolve()
    read_manifest_path = args.read_manifest.resolve()
    try:
        manifest = json_load(manifest_path)
        read_manifest = json_load(read_manifest_path)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    errors, missing, abundances, queries = validate_manifest(
        manifest, manifest_path, read_manifest, read_manifest_path, check_content=True
    )
    if errors:
        print("invalid Track B manifest:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 2
    run_id = f"track-b-{int(time.time())}"
    if args.self_check:
        if args.output_dir.exists():
            issue_missing(missing, f"output directory already exists; choose a new path ({args.output_dir})")
        status = "ready" if not missing else "unavailable"
        print(json.dumps({"schema_version": SCHEMA_VERSION, "track": "B", "status": status, "missing_prerequisites": sorted(set(missing))}, indent=2))
        return 0
    try:
        prepare_output(args.output_dir)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    run_record: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "track": "B",
        "record_type": "run_status",
        "run_id": run_id,
        "started_at_utc": utc_now(),
        "manifest": str(manifest_path),
        "read_manifest": str(read_manifest_path),
        "abundance_grid": abundances,
        "commands": redact({name: manifest.get(name) for name in ("mixer", "assembler", "coordinate_extractor", "jam")}),
    }
    if missing:
        unavailable_files(args.output_dir, run_record, missing)
        return UNAVAILABLE_EXIT

    jobs_dir = args.output_dir / "jobs"
    jobs_dir.mkdir()
    records: list[dict[str, Any]] = []
    failed = False
    for read_set in read_manifest["read_sets"]:
        read_set_id = read_set["read_set_id"]
        r1_path = resolve_path(read_set["r1"], read_manifest_path.parent)
        r2_path = resolve_path(read_set.get("r2"), read_manifest_path.parent)
        if r1_path is None:
            continue
        for query in queries:
            query_id = query["query_id"]
            query_kind = query["query_kind"]
            topology = query["topology"]
            query_path = query["_fasta_path"]
            query_sequence = fasta_sequence(query_path)
            for abundance in abundances:
                label = abundance["label"]
                job_id = f"{read_set_id}__{query_id}__{label}"
                job_dir = jobs_dir / job_id
                job_dir.mkdir()
                mixed_r1 = job_dir / "mixed_R1.fastq.gz"
                mixed_r2 = job_dir / "mixed_R2.fastq.gz"
                assembly = job_dir / "assembly.fasta"
                coordinates = job_dir / "coordinates.json"
                trace_output = job_dir / "trace.jsonl"
                values = {
                    "job_id": job_id,
                    "read_manifest": str(read_manifest_path),
                    "read_set_id": read_set_id,
                    "input_r1": str(r1_path),
                    "input_r2": str(r2_path) if r2_path is not None else "",
                    "query_fasta": str(query_path),
                    "query_id": query_id,
                    "query_kind": query_kind,
                    "topology": topology,
                    "query_length": len(query_sequence),
                    # Legacy aliases remain available to old command templates.
                    "plasmid": str(query_path),
                    "plasmid_id": query_id,
                    "plasmid_length": len(query_sequence),
                    "abundance": abundance["fraction"],
                    "abundance_label": label,
                    "output_dir": str(job_dir),
                    "mixed_r1": str(mixed_r1),
                    "mixed_r2": str(mixed_r2),
                    "assembly": str(assembly),
                    "coordinates": str(coordinates),
                    "trace_output": str(trace_output),
                }
                (job_dir / "job.json").write_text(json.dumps(values, indent=2) + "\n", encoding="utf-8")
                commands: list[dict[str, Any]] = []
                status = "complete"
                failure: str | None = None
                stages = (
                    ("mix", manifest["mixer"], [mixed_r1, mixed_r2] if r2_path is not None else [mixed_r1]),
                    ("assemble", manifest["assembler"], [assembly]),
                    ("coordinates", manifest["coordinate_extractor"], [coordinates]),
                    ("jam_trace", manifest["jam"], [trace_output]),
                )
                for stage, tool, expected_outputs in stages:
                    try:
                        command = render_command(tool["command"], values)
                    except ValueError as exc:
                        status, failure = "failed", str(exc)
                        break
                    record = run_command(command, job_dir, job_dir / f"{stage}.log", tool, stage)
                    record["command"] = redact(record["command"])
                    commands.append(record)
                    if record["exit_code"] != 0:
                        status = "failed"
                        failure = f"{stage} command failed with exit code {record['exit_code']}"
                        break
                    missing_outputs = [str(expected) for expected in expected_outputs if not expected.is_file()]
                    if missing_outputs:
                        status = "failed"
                        failure = f"{stage} command did not create {', '.join(missing_outputs)}"
                        break
                output_record: dict[str, Any] = {
                    "schema_version": SCHEMA_VERSION,
                    "track": "B",
                    "record_type": "job_result",
                    "run_id": run_id,
                    "job_id": job_id,
                    "read_set_id": read_set_id,
                    "query_id": query_id,
                    "query_kind": query_kind,
                    "topology": topology,
                    "query_length": len(query_sequence),
                    "abundance_label": label,
                    "abundance_fraction": abundance["fraction"],
                    "status": status,
                    "commands": commands,
                    "failure": failure,
                }
                if query_kind == "plasmid":
                    output_record["plasmid_id"] = query_id
                if status == "complete":
                    try:
                        sequence, coordinates_value, assembled_intervals = coordinate_truth(
                            coordinates, query_id, query_kind, topology, len(query_sequence), job_id
                        )
                        sequence_sha = hashlib.sha256(sequence.encode("ascii")).hexdigest()
                        sequence_path = job_dir / "assembled_query.fasta"
                        with sequence_path.open("w", encoding="ascii") as handle:
                            handle.write(f">{query_id}|{job_id}|assembled\n")
                            for offset in range(0, len(sequence), 80):
                                handle.write(sequence[offset : offset + 80] + "\n")
                        support = trace_support(trace_output, len(query_sequence), assembled_intervals)
                        if query_kind == "plasmid":
                            support["assembled_plasmid_bases"] = support["assembled_query_bases"]
                            support["supported_fraction_of_plasmid"] = support["supported_fraction_of_query"]
                        sequence_record = {
                            "path": str(sequence_path),
                            "sha256": sequence_sha,
                            "fasta_sha256": sha256(sequence_path),
                            "length": len(sequence),
                        }
                        output_record["assembled_query_sequence"] = sequence_record
                        if query_kind == "plasmid":
                            output_record["plasmid_id"] = query_id
                            output_record["assembled_plasmid_sequence"] = sequence_record
                        output_record.update(
                            {
                                "assembled_coordinates": coordinates_value,
                                "assembled_intervals": support["assembled_intervals"],
                                "trace_evidence": support,
                                "measurement_scope": "intersection of Jam-supported intervals and independently extracted assembled query-element intervals",
                            }
                        )
                    except (OSError, ValueError, UnicodeError, json.JSONDecodeError, subprocess.SubprocessError) as exc:
                        output_record["status"] = "failed"
                        output_record["failure"] = f"cannot validate assembled sequence or trace evidence: {exc}"
                if output_record["status"] != "complete":
                    failed = True
                records.append(output_record)
    run_record.update(
        {
            "finished_at_utc": utc_now(),
            "status": "failed" if failed else "complete",
            "jobs_total": len(records),
            "jobs_complete": sum(record["status"] == "complete" for record in records),
            "jobs_failed": sum(record["status"] == "failed" for record in records),
        }
    )
    (args.output_dir / "results.jsonl").write_text(
        "".join(json.dumps(record, sort_keys=True) + "\n" for record in records), encoding="utf-8"
    )
    (args.output_dir / "run.json").write_text(json.dumps(run_record, indent=2) + "\n", encoding="utf-8")
    return 1 if failed else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--read-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--self-check", action="store_true", help="validate manifests without executing tools")
    args = parser.parse_args()
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
