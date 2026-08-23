#!/usr/bin/env python3
"""Run the independently supported natural-positive benchmark track.

Track C keeps independent support separate from Jam evidence.  A metagenome
without an external support record is ``unresolved``; it is never called
negative merely because it is absent from the support manifest.  The harness
requires explicit checksums and licensing metadata, records the pinned Jam
command, and emits ``status=unavailable`` when the required collection or
support evidence is not present.  ``--self-check`` performs validation only.
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
    validate_sha(declared_sha, f"{label}.sha256", errors)
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
    """Normalize query_elements and legacy query_plasmids entries."""
    new_entries = manifest.get("query_elements")
    legacy_entries = manifest.get("query_plasmids")
    if new_entries is not None and legacy_entries is not None:
        issue_error(errors, "provide query_elements or legacy query_plasmids, not both")
        return []
    legacy = new_entries is None
    entries = legacy_entries if legacy else new_entries
    if not isinstance(entries, list) or not entries:
        issue_missing(missing, "at least one query element")
        return []
    valid: list[dict[str, Any]] = []
    seen: set[str] = set()
    for index, entry in enumerate(entries):
        label = f"{'query_plasmids' if legacy else 'query_elements'}[{index}]"
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
        validate_file(entry.get("fasta"), entry.get("sha256"), entry.get("license"), label, manifest_path.parent, missing, errors, check_content)
        path = resolve_path(entry.get("fasta"), manifest_path.parent)
        if path is not None and path.is_file():
            valid.append({**entry, "query_id": query_id, "query_kind": query_kind, "topology": topology, "_fasta_path": path})
    return valid


def validate_support_manifest(
    support: dict[str, Any], support_path: Path, parent: dict[str, Any], missing: list[str], errors: list[str], check_content: bool
) -> list[dict[str, Any]]:
    if support.get("schema_version") != SCHEMA_VERSION:
        issue_error(errors, "support manifest schema_version must be 1.0.0")
    if not valid_text(support.get("support_manifest_id")):
        issue_error(errors, "support_manifest_id is required")
    if support.get("support_manifest_id") != parent.get("support_manifest_id"):
        issue_error(errors, "support manifest id does not match the natural benchmark manifest")
    validate_license(support.get("license"), "support manifest", errors)
    declared_sha = parent.get("support_manifest_sha256")
    validate_sha(declared_sha, "support_manifest_sha256", errors)
    if support_path.is_file() and isinstance(declared_sha, str) and SHA256_RE.fullmatch(declared_sha):
        actual = sha256(support_path)
        if actual.lower() != declared_sha.lower():
            issue_missing(missing, f"support manifest checksum mismatch (declared {declared_sha}, actual {actual})")
    else:
        issue_missing(missing, f"support manifest file not found ({support_path})")

    records = support.get("support_records")
    if not isinstance(records, list) or not records:
        issue_missing(missing, "at least one independent support record")
        return []
    valid: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for index, record in enumerate(records):
        label = f"support_records[{index}]"
        if not isinstance(record, dict):
            issue_error(errors, f"{label}: object is required")
            continue
        for field in ("metagenome_id", "support_type", "evidence_id", "evidence_source"):
            if not valid_text(record.get(field)):
                issue_error(errors, f"{label}.{field} is required")
        query_id = record.get("query_id", record.get("plasmid_id"))
        query_kind = record.get("query_kind", "plasmid" if record.get("plasmid_id") is not None else None)
        topology = record.get("topology", "circular" if query_kind == "plasmid" else "unknown")
        if not valid_text(query_id):
            issue_error(errors, f"{label}.query_id is required")
        if query_kind not in QUERY_KINDS:
            issue_error(errors, f"{label}.query_kind must be plasmid or phage")
        if topology not in TOPOLOGIES:
            issue_error(errors, f"{label}.topology must be circular, linear, or unknown")
        validate_sha(record.get("evidence_sha256"), f"{label}.evidence_sha256", errors)
        validate_license(record.get("license"), label, errors)
        support_type = record.get("support_type")
        if isinstance(support_type, str) and support_type.lower() in {"truth", "none", "unknown", "absence"}:
            issue_error(errors, f"{label}.support_type must identify independent evidence")
        key = (str(record.get("metagenome_id")), str(query_id))
        if key in seen:
            issue_error(errors, f"duplicate support record for {key[0]} and {key[1]}")
            continue
        seen.add(key)
        evidence_path = resolve_path(record.get("evidence_path"), support_path.parent)
        if record.get("evidence_path") is not None:
            if evidence_path is None or not evidence_path.is_file():
                issue_missing(missing, f"{label}: evidence file not found ({evidence_path})")
            elif check_content:
                actual = sha256(evidence_path)
                if actual.lower() != str(record.get("evidence_sha256")).lower():
                    issue_missing(missing, f"{label}: evidence checksum mismatch (declared {record.get('evidence_sha256')}, actual {actual})")
        elif not valid_text(record.get("evidence_uri")):
            issue_missing(missing, f"{label}: evidence_path or evidence_uri")
        valid.append({**record, "query_id": query_id, "query_kind": query_kind, "topology": topology})
    return valid


def validate_manifest(
    manifest: dict[str, Any], manifest_path: Path, support: dict[str, Any], support_path: Path, check_content: bool
) -> tuple[list[str], list[str], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    missing: list[str] = []
    errors: list[str] = []
    if manifest.get("schema_version") != SCHEMA_VERSION:
        issue_error(errors, "benchmark schema_version must be 1.0.0")
    if manifest.get("track") != "C":
        issue_error(errors, "benchmark track must be C")
    if not valid_text(manifest.get("dataset_id")):
        issue_error(errors, "dataset_id is required")
    validate_license(manifest.get("license"), "benchmark", errors)
    for field in ("support_manifest_id", "support_manifest_sha256", "metagenomes", "jam"):
        if field not in manifest:
            issue_missing(missing, f"benchmark field: {field}")
    if "query_elements" not in manifest and "query_plasmids" not in manifest:
        issue_missing(missing, "query_elements (or legacy query_plasmids)")
    valid_support = validate_support_manifest(support, support_path, manifest, missing, errors, check_content)

    valid_queries = query_entries(manifest, manifest_path, missing, errors, check_content)

    metagenomes = manifest.get("metagenomes")
    valid_metagenomes: list[dict[str, Any]] = []
    if not isinstance(metagenomes, list) or not metagenomes:
        issue_missing(missing, "at least one natural metagenome assembly")
    else:
        seen = set()
        for index, metagenome in enumerate(metagenomes):
            label = f"metagenomes[{index}]"
            if not isinstance(metagenome, dict) or not valid_text(metagenome.get("metagenome_id")):
                issue_error(errors, f"{label}: metagenome_id and assembly metadata are required")
                continue
            metagenome_id = metagenome["metagenome_id"]
            if metagenome_id in seen:
                issue_error(errors, f"duplicate metagenome_id: {metagenome_id}")
                continue
            seen.add(metagenome_id)
            validate_file(metagenome.get("assembly"), metagenome.get("sha256"), metagenome.get("license"), label, manifest_path.parent, missing, errors, check_content)
            path = resolve_path(metagenome.get("assembly"), manifest_path.parent)
            if path is not None and path.is_file():
                valid_metagenomes.append({**metagenome, "_assembly_path": path})
    validate_tool(manifest.get("jam"), "jam", missing, errors)
    query_metadata = {query["query_id"]: query for query in valid_queries}
    for record in valid_support:
        query = query_metadata.get(record["query_id"])
        if query is None:
            issue_error(errors, f"support record references unknown query_id: {record['query_id']}")
        elif record["query_kind"] != query["query_kind"] or record["topology"] != query["topology"]:
            issue_error(errors, f"support record query metadata does not match query_id: {record['query_id']}")
    return errors, missing, valid_queries, valid_metagenomes, valid_support


def prepare_output(path: Path) -> None:
    if path.exists():
        raise ValueError(f"refusing to reuse existing output directory: {path}")
    path.mkdir(parents=True)


def render_command(command: list[str], values: dict[str, Any]) -> list[str]:
    rendered: list[str] = []
    for part in command:
        try:
            rendered.append(part.format(**{key: str(item) for key, item in values.items()}))
        except KeyError as exc:
            raise ValueError(f"command contains unknown placeholder: {{{exc.args[0]}}}") from exc
    return rendered


def run_command(command: list[str], cwd: Path, log_path: Path, tool: dict[str, Any]) -> dict[str, Any]:
    started = utc_now()
    start = time.perf_counter()
    with log_path.open("w", encoding="utf-8") as log:
        try:
            completed = subprocess.run(command, cwd=cwd, stdout=log, stderr=subprocess.STDOUT, check=False)
            error = None
        except OSError as exc:
            completed = None
            error = str(exc)
    return {
        "tool": {"name": tool.get("name"), "version": tool.get("version"), "image": tool.get("image"), "image_digest": tool.get("image_digest")},
        "command": command,
        "started_at_utc": started,
        "finished_at_utc": utc_now(),
        "wall_seconds": time.perf_counter() - start,
        "exit_code": completed.returncode if completed is not None else None,
        "log": str(log_path),
        "error": error,
    }


def trace_summary(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ValueError(f"Jam did not create trace output: {path}")
    if path.suffix.lower() in {".zst", ".zstd"}:
        completed = subprocess.run(["zstd", "--decompress", "--stdout", str(path)], stdout=subprocess.PIPE, check=False)
        if completed.returncode != 0:
            raise ValueError(f"cannot decompress trace output: {path}")
        text = completed.stdout.decode("utf-8")
    else:
        text = path.read_text(encoding="utf-8")
    candidates = 0
    alignments = 0
    coverages: list[float] = []
    records = 0
    for line_number, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        value = json.loads(line)
        if not isinstance(value, dict):
            raise ValueError(f"trace record {line_number} is not an object")
        records += 1
        if value.get("record_type") != "metagenome_result":
            continue
        if value.get("candidate") is not None:
            candidates += 1
        alignments += len(value.get("alignments", []))
        coverage = value.get("coverage") or {}
        if isinstance(coverage.get("supported_fraction"), (int, float)):
            coverages.append(float(coverage["supported_fraction"]))
    return {
        "records": records,
        "candidate_records": candidates,
        "alignment_count": alignments,
        "mean_supported_fraction": sum(coverages) / len(coverages) if coverages else None,
    }


def unavailable_files(output: Path, run: dict[str, Any], missing: list[str]) -> None:
    run["status"] = "unavailable"
    run["missing_prerequisites"] = sorted(set(missing))
    (output / "run.json").write_text(json.dumps(run, indent=2) + "\n", encoding="utf-8")
    result = {
        "schema_version": SCHEMA_VERSION,
        "track": "C",
        "record_type": "run_status",
        "run_id": run["run_id"],
        "status": "unavailable",
        "missing_prerequisites": sorted(set(missing)),
        "note": "No biological absence is inferred from unavailable or un-supported natural records.",
    }
    (output / "results.jsonl").write_text(json.dumps(result) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> int:
    manifest_path = args.manifest.resolve()
    support_path = args.support_manifest.resolve()
    try:
        manifest = json_load(manifest_path)
        support = json_load(support_path)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    errors, missing, queries, metagenomes, supports = validate_manifest(
        manifest, manifest_path, support, support_path, check_content=True
    )
    if errors:
        print("invalid Track C manifest:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 2
    run_id = f"track-c-{int(time.time())}"
    if args.self_check:
        if args.output_dir.exists():
            issue_missing(missing, f"output directory already exists; choose a new path ({args.output_dir})")
        status = "ready" if not missing else "unavailable"
        print(json.dumps({"schema_version": SCHEMA_VERSION, "track": "C", "status": status, "missing_prerequisites": sorted(set(missing))}, indent=2))
        return 0
    try:
        prepare_output(args.output_dir)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    run_record: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "track": "C",
        "record_type": "run_status",
        "run_id": run_id,
        "started_at_utc": utc_now(),
        "manifest": str(manifest_path),
        "support_manifest": str(support_path),
        "jam": redact(manifest.get("jam")),
        "interpretation_rule": "independent support is reported separately; missing support is unresolved, never absence",
    }
    if missing:
        unavailable_files(args.output_dir, run_record, missing)
        return UNAVAILABLE_EXIT

    support_by_pair = {(record["metagenome_id"], record["query_id"]): record for record in supports}
    jobs_dir = args.output_dir / "jobs"
    jobs_dir.mkdir()
    records: list[dict[str, Any]] = []
    failed = False
    for query in queries:
        for metagenome in metagenomes:
            query_id = query["query_id"]
            query_kind = query["query_kind"]
            topology = query["topology"]
            metagenome_id = metagenome["metagenome_id"]
            job_id = f"{query_id}__{metagenome_id}"
            job_dir = jobs_dir / job_id
            job_dir.mkdir()
            trace_output = job_dir / "trace.jsonl"
            query_path = query["_fasta_path"]
            assembly_path = metagenome["_assembly_path"]
            values = {
                "job_id": job_id,
                "query_id": query_id,
                "query_kind": query_kind,
                "topology": topology,
                "query_fasta": str(query_path),
                # Legacy command templates remain supported for plasmid runs.
                "plasmid_id": query_id,
                "plasmid": str(query_path),
                "metagenome_id": metagenome_id,
                "metagenome": str(assembly_path),
                "output_dir": str(job_dir),
                "trace_output": str(trace_output),
            }
            support_record = support_by_pair.get((metagenome_id, query_id))
            commands: list[dict[str, Any]] = []
            failure: str | None = None
            status = "complete"
            try:
                command = render_command(manifest["jam"]["command"], values)
                command_record = run_command(command, job_dir, job_dir / "jam.log", manifest["jam"])
                command_record["command"] = redact(command_record["command"])
                commands.append(command_record)
                if command_record["exit_code"] != 0:
                    status = "failed"
                    failure = f"Jam command failed with exit code {command_record['exit_code']}"
                else:
                    summary = trace_summary(trace_output)
            except (OSError, ValueError, UnicodeError, json.JSONDecodeError, subprocess.SubprocessError) as exc:
                status = "failed"
                failure = str(exc)
                summary = None
            if status == "failed":
                failed = True
            records.append(
                {
                    "schema_version": SCHEMA_VERSION,
                    "track": "C",
                    "record_type": "natural_result",
                    "run_id": run_id,
                    "job_id": job_id,
                    "query_id": query_id,
                    "query_kind": query_kind,
                    "topology": topology,
                    "metagenome_id": metagenome_id,
                    "status": status,
                    "trace_output": str(trace_output) if trace_output.is_file() else None,
                    "trace_summary": summary,
                    "commands": commands,
                    "independent_support": redact(support_record),
                    "support_status": "supported" if support_record is not None else "not_provided",
                    "biological_interpretation": "supported_by_independent_evidence" if support_record is not None else "unresolved_no_independent_support",
                    "failure": failure,
                }
            )
            if query_kind == "plasmid":
                records[-1]["plasmid_id"] = query_id
    run_record.update(
        {
            "finished_at_utc": utc_now(),
            "status": "failed" if failed else "complete",
            "jobs_total": len(records),
            "jobs_failed": sum(record["status"] == "failed" for record in records),
            "jobs_with_independent_support": sum(record["support_status"] == "supported" for record in records),
            "jobs_unresolved": sum(record["support_status"] != "supported" for record in records),
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
    parser.add_argument("--support-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--self-check", action="store_true", help="validate manifests without executing Jam")
    args = parser.parse_args()
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
