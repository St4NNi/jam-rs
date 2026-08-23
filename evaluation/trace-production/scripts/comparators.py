#!/usr/bin/env python3
"""Run matched candidate and direct comparator tracks.

Each selected assembly is supplied unchanged to the requested comparator.  The
BLAST tabular and minimap2 PAF output is normalized by the shared
``scripts/coverage_normalize.py`` implementation.  A precomputed ``jam
trace`` JSONL result can be normalized with the same command, or the runner
can invoke ``jam trace`` when its database and catalog are supplied.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlsplit, urlunsplit

SCRIPT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(SCRIPT_ROOT / "scripts"))

from coverage_normalize import (  # noqa: E402
    NormalizationError,
    _truth_intervals,
    normalize_records,
    parse_blast6,
    parse_paf,
    read_fasta_length,
)


SCHEMA_VERSION = "1.0.0"
FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".fastq", ".fq"}
_URI_RE = re.compile(r"\b[a-z][a-z0-9+.-]*://[^\s]+", re.IGNORECASE)


@dataclass(frozen=True)
class Assembly:
    metagenome_id: str
    path: Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def redact(value: str) -> str:
    """Remove credentials and query fragments from URI-like command values."""

    if "://" not in value:
        return value
    parsed = urlsplit(value)
    if not parsed.scheme or not parsed.netloc:
        return value
    host = parsed.hostname or ""
    if parsed.port:
        host = f"{host}:{parsed.port}"
    return urlunsplit((parsed.scheme, host, parsed.path, "", ""))


def command_version(executable: str) -> str | None:
    flag = "-version" if Path(executable).name == "blastn" else "--version"
    try:
        completed = subprocess.run(
            [executable, flag],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
    except OSError:
        return None
    text = completed.stdout.strip()
    return text.splitlines()[0] if text else None


def command_text(command: list[str]) -> list[str]:
    return [redact(str(part)) for part in command]


def redact_message(message: str) -> str:
    return _URI_RE.sub(lambda match: redact(match.group(0)), message)


def _truth_checksum(path: Path | None) -> str | None:
    return sha256(path) if path and path.is_file() else None


def _field(row: dict[str, str], *names: str) -> str | None:
    for name in names:
        value = row.get(name)
        if value is not None and value.strip():
            return value.strip()
    return None


def _truthy(value: str | None) -> bool:
    return value is not None and value.strip().lower() in {"1", "true", "yes", "selected"}


def _manifest_rows(path: Path) -> list[dict[str, str]]:
    if path.suffix.lower() == ".jsonl":
        rows: list[dict[str, str]] = []
        seen: set[str] = set()
        with path.open(encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, 1):
                if not raw.strip():
                    continue
                try:
                    record = json.loads(raw)
                except json.JSONDecodeError as exc:
                    raise ValueError(
                        f"invalid JSONL candidate manifest at line {line_number}: {exc.msg}"
                    ) from exc
                if not isinstance(record, dict):
                    raise ValueError(
                        f"candidate manifest line {line_number} is not a JSON object"
                    )
                if record.get("record_type") != "metagenome_result":
                    continue
                metagenome_id = record.get("metagenome_id")
                if not isinstance(metagenome_id, str) or not metagenome_id.strip():
                    raise ValueError(
                        f"candidate manifest line {line_number} has no valid metagenome_id"
                    )
                metagenome_id = metagenome_id.strip()
                if metagenome_id in seen:
                    raise ValueError(
                        f"candidate manifest repeats metagenome_id {metagenome_id!r}"
                    )
                seen.add(metagenome_id)
                rows.append(
                    {
                        "metagenome_id": metagenome_id,
                        "candidate": "true" if record.get("candidate") is not None else "false",
                    }
                )
        if not rows:
            raise ValueError("candidate JSONL contains no metagenome_result records")
        return rows
    with path.open(encoding="utf-8", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = "\t" if "\t" in sample else ","
        return [dict(row) for row in csv.DictReader(handle, delimiter=delimiter)]


def discover_assemblies(assemblies_dir: Path) -> list[Assembly]:
    if not assemblies_dir.is_dir():
        raise ValueError(f"assemblies directory does not exist: {assemblies_dir}")
    assemblies = [
        Assembly(path.name, path)
        for path in assemblies_dir.iterdir()
        if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES
    ]
    return sorted(assemblies, key=lambda item: (item.metagenome_id, str(item.path)))


def select_assemblies(
    assemblies_dir: Path,
    candidate_manifest: Path | None,
    scope: str,
    direct_limit: int | None,
) -> tuple[list[Assembly], dict]:
    available = {item.metagenome_id: item for item in discover_assemblies(assemblies_dir)}
    if candidate_manifest is None:
        if scope == "candidate":
            selected = list(available.values())
            source = "assemblies_dir_without_candidate_manifest"
        else:
            selected = list(available.values())
            source = "assemblies_dir"
    else:
        manifest_base = candidate_manifest.parent
        selected = []
        seen: set[str] = set()
        for row in _manifest_rows(candidate_manifest):
            candidate_flag = _field(row, "candidate", "selected", "include")
            if candidate_flag is not None and not _truthy(candidate_flag):
                continue
            metagenome_id = _field(row, "metagenome_id", "assembly_id", "sample_id", "id")
            assembly_value = _field(
                row, "assembly", "assembly_path", "raw", "path", "file", "resource"
            )
            if metagenome_id is None and assembly_value:
                metagenome_id = Path(assembly_value).name
            if not metagenome_id or metagenome_id in seen:
                continue
            assembly_path = Path(assembly_value) if assembly_value else available.get(metagenome_id, Assembly(metagenome_id, Path())).path
            if not assembly_path.is_absolute():
                assembly_path = (manifest_base / assembly_path).resolve()
            if not assembly_path.is_file():
                fallback = available.get(metagenome_id)
                if fallback is not None:
                    assembly_path = fallback.path
            if not assembly_path.is_file():
                raise ValueError(
                    f"candidate manifest row {metagenome_id!r} has no readable assembly: {assembly_path}"
                )
            selected.append(Assembly(metagenome_id, assembly_path))
            seen.add(metagenome_id)
        selected.sort(key=lambda item: (item.metagenome_id, str(item.path)))
        source = str(candidate_manifest)
    total_count = len(selected) if candidate_manifest is not None else len(available)
    if scope == "direct" and candidate_manifest is not None:
        # A direct comparison always uses the complete collection represented by
        # the manifest, regardless of candidate flags.
        selected = []
        seen = set()
        for row in _manifest_rows(candidate_manifest):
            metagenome_id = _field(row, "metagenome_id", "assembly_id", "sample_id", "id")
            assembly_value = _field(
                row, "assembly", "assembly_path", "raw", "path", "file", "resource"
            )
            if metagenome_id is None and assembly_value:
                metagenome_id = Path(assembly_value).name
            if not metagenome_id or metagenome_id in seen:
                continue
            assembly_path = Path(assembly_value) if assembly_value else available.get(metagenome_id, Assembly(metagenome_id, Path())).path
            if not assembly_path.is_absolute():
                assembly_path = (candidate_manifest.parent / assembly_path).resolve()
            if not assembly_path.is_file():
                fallback = available.get(metagenome_id)
                if fallback is not None:
                    assembly_path = fallback.path
            if not assembly_path.is_file():
                raise ValueError(
                    f"direct manifest row {metagenome_id!r} has no readable assembly: {assembly_path}"
                )
            selected.append(Assembly(metagenome_id, assembly_path))
            seen.add(metagenome_id)
        selected.sort(key=lambda item: (item.metagenome_id, str(item.path)))
        total_count = len(selected)
    if scope == "direct" and direct_limit is not None:
        if direct_limit < 1:
            raise ValueError("--direct-limit must be positive")
        selected = selected[:direct_limit]
    effective_scope = (
        "candidate"
        if scope == "candidate"
        else "direct_subset"
        if direct_limit is not None and len(selected) < total_count
        else "direct_full"
    )
    return selected, {
        "source": source,
        "total_count": total_count,
        "selected_count": len(selected),
        "limit": direct_limit if scope == "direct" else None,
        "ids": [item.metagenome_id for item in selected],
        "effective_scope": effective_scope,
    }


def read_text(path: Path) -> str:
    if path.suffix.lower() in {".zst", ".zstd"}:
        completed = subprocess.run(
            ["zstd", "--decompress", "--stdout", str(path)],
            check=True,
            stdout=subprocess.PIPE,
        )
        return completed.stdout.decode("utf-8")
    return path.read_text(encoding="utf-8")


def _truth_rows(path: Path) -> list[dict[str, str]]:
    if path.suffix.lower() == ".json":
        value = json.loads(path.read_text(encoding="utf-8"))
        rows = value.get("intervals", value) if isinstance(value, dict) else value
        return [row for row in rows if isinstance(row, dict)] if isinstance(rows, list) else []
    csv.field_size_limit(1024 * 1024 * 1024)
    with path.open(encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def _truth_query_filter(path: Path | None, query_id: str | None) -> str | None:
    """Select an exact truth query ID without rejecting suffixed controlled IDs."""

    if path is None or query_id is None:
        return None
    try:
        rows = _truth_rows(path)
    except (OSError, json.JSONDecodeError, UnicodeError):
        return None
    exact = {
        row.get("query_id") or row.get("plasmid_id") or row.get("reference_id")
        for row in rows
    }
    if query_id in exact:
        return query_id
    aliases = {
        row.get("query_id")
        for row in rows
        if row.get("query_id")
        and (
            row.get("source_query_id") == query_id
            or query_id.startswith(f"{row.get('source_query_id', '')}__")
        )
    }
    return next(iter(aliases)) if len(aliases) == 1 else None


def load_truth(
    path: Path | None,
    metagenome_id: str,
    query_length: int,
    query_id: str | None = None,
) -> list[dict[str, int]] | None:
    if path is None:
        return None
    return _truth_intervals(
        path,
        metagenome_id,
        query_length,
        _truth_query_filter(path, query_id),
    )


def result_for_records(
    assembly: Assembly,
    records: list[dict],
    query_length: int,
    truth_path: Path | None,
    command: list[str],
    version: str | None,
    extra: dict | None = None,
    query_id: str | None = None,
) -> dict:
    truth = load_truth(truth_path, assembly.metagenome_id, query_length, query_id)
    coverage = normalize_records(records, query_length, truth)
    output = {
        "metagenome_id": assembly.metagenome_id,
        "assembly": str(assembly.path),
        "status": "measured",
        "command": command_text(command),
        "version": version,
        "records": records,
        "coverage": coverage,
    }
    if extra:
        output.update(extra)
    return output


def failed_result(
    assembly: Assembly,
    query_length: int,
    error: str,
    command: list[str],
    version: str | None,
) -> dict:
    return {
        "metagenome_id": assembly.metagenome_id,
        "assembly": str(assembly.path),
        "status": "failed",
        "command": command_text(command),
        "version": version,
        "records": [],
        "coverage": normalize_records([], query_length, None),
        "error": redact_message(error),
    }


def aggregate_metrics(results: list[dict], query_length: int) -> dict:
    measured = [item for item in results if item.get("status") == "measured"]
    observed = sum(item["coverage"]["supported_bases"] for item in measured)
    truth_bases = sum(
        (item.get("coverage", {}).get("metrics") or {}).get("truth_bases", 0)
        for item in measured
    )
    intersection = sum(
        (item.get("coverage", {}).get("metrics") or {}).get("intersection_bases", 0)
        for item in measured
    )
    gap_error = sum(
        (item.get("coverage", {}).get("metrics") or {}).get("gap_error_bases", 0)
        for item in measured
    )
    precision_values = [
        item["coverage"]["metrics"]["interval_precision"]
        for item in measured
        if item.get("coverage", {}).get("metrics", {}).get("interval_precision") is not None
    ]
    recall_values = [
        item["coverage"]["metrics"]["interval_recall"]
        for item in measured
        if item.get("coverage", {}).get("metrics", {}).get("interval_recall") is not None
    ]
    return {
        "assemblies_total": len(results),
        "assemblies_with_records": sum(bool(item.get("records")) for item in measured),
        "observed_bases": observed,
        "truth_bases": truth_bases or None,
        "intersection_bases": intersection if truth_bases else None,
        "base_precision": intersection / observed if observed and truth_bases else None,
        "base_recall": intersection / truth_bases if truth_bases else None,
        "interval_precision": sum(precision_values) / len(precision_values) if precision_values else None,
        "interval_recall": sum(recall_values) / len(recall_values) if recall_values else None,
        "gap_error_bases": gap_error if truth_bases else None,
        "gap_error_fraction": gap_error / (query_length * len(measured)) if truth_bases and measured else None,
        "truth_status": "measured" if truth_bases else "unavailable",
    }


def run_alignment_tool(
    tool: str,
    executable: str,
    query: Path,
    query_id: str,
    query_length: int,
    assemblies: list[Assembly],
    truth_path: Path | None,
    args: argparse.Namespace,
) -> tuple[list[dict], list[list[str]]]:
    version = command_version(executable)
    results: list[dict] = []
    commands: list[list[str]] = []
    for assembly in assemblies:
        if tool == "blastn":
            command = [
                executable,
                "-query",
                str(query),
                "-subject",
                str(assembly.path),
                "-task",
                args.blast_task,
                "-outfmt",
                "6 qseqid sseqid pident length qlen qstart qend sstart send bitscore",
                "-evalue",
                str(args.evalue),
                "-max_hsps",
                str(args.max_hsps),
                "-num_threads",
                str(args.threads),
            ]
        else:
            command = [
                executable,
                "-x",
                args.minimap_preset,
                "--secondary=no",
                "-t",
                str(args.threads),
                "-c",
                "--eqx",
                str(assembly.path),
                str(query),
            ]
        commands.append(command)
        try:
            completed = subprocess.run(
                command,
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            if completed.returncode != 0:
                results.append(failed_result(assembly, query_length, completed.stderr.strip() or f"exit {completed.returncode}", command, version))
                continue
            if tool == "blastn":
                records = parse_blast6(completed.stdout, query_id, query_length, assembly.metagenome_id, args.min_pident)
            else:
                records = parse_paf(completed.stdout, query_id, query_length, assembly.metagenome_id, args.min_mapq)
            results.append(
                result_for_records(
                    assembly,
                    records,
                    query_length,
                    truth_path,
                    command,
                    version,
                    query_id=query_id,
                )
            )
        except (OSError, NormalizationError, UnicodeError) as exc:
            results.append(failed_result(assembly, query_length, str(exc), command, version))
    return results, commands


def jam_result_records(
    text: str,
    query_id: str,
    query_length: int,
    assemblies: list[Assembly],
    truth_path: Path | None,
    command: list[str],
    version: str | None,
) -> list[dict]:
    records = [json.loads(line) for line in text.splitlines() if line.strip()]
    header = next((record for record in records if record.get("record_type") == "run_header"), None)
    if header is None:
        raise NormalizationError("jam JSONL has no run_header record")
    if header.get("plasmid_id") != query_id or int(header.get("plasmid_length", 0)) != query_length:
        raise NormalizationError("jam JSONL query identity or length disagrees with the selected query")
    by_id = {
        record.get("metagenome_id"): record
        for record in records
        if record.get("record_type") == "metagenome_result"
    }
    results: list[dict] = []
    for assembly in assemblies:
        source = by_id.get(assembly.metagenome_id)
        if source is None:
            results.append(
                result_for_records(
                    assembly,
                    [],
                    query_length,
                    truth_path,
                    command,
                    version,
                    {"jam_status": "not_selected"},
                    query_id=query_id,
                )
            )
            continue
        intervals = []
        for interval in (source.get("coverage") or {}).get("primary_intervals", []):
            intervals.append(
                {
                    "metagenome_id": assembly.metagenome_id,
                    "contig_id": "jam-primary",
                    "strand": "unknown",
                    "intervals": [{"start": int(interval["start"]), "end": int(interval["end"])}],
                }
            )
        result = result_for_records(
            assembly,
            intervals,
            query_length,
            truth_path,
            command,
            version,
            {"jam_status": source.get("status")},
            query_id,
        )
        results.append(result)
    return results


def run_jam(
    executable: str,
    query: Path,
    query_id: str,
    query_length: int,
    assemblies: list[Assembly],
    truth_path: Path | None,
    args: argparse.Namespace,
) -> tuple[list[dict], list[list[str]]]:
    version = command_version(executable)
    if args.jam_output:
        command = [executable, "trace", "(precomputed output)", str(args.jam_output)]
        try:
            return jam_result_records(read_text(args.jam_output), query_id, query_length, assemblies, truth_path, command, version), [command]
        except (OSError, subprocess.SubprocessError, json.JSONDecodeError, ValueError, NormalizationError) as exc:
            return [failed_result(item, query_length, str(exc), command, version) for item in assemblies], [command]
    required = (args.jam_database, args.jam_catalog)
    if any(value is None for value in required):
        raise ValueError("jam requires --jam-output or both --jam-database and --jam-catalog")
    with tempfile.TemporaryDirectory(prefix="jam-trace-comparator-") as directory:
        output = Path(directory) / "trace.jsonl"
        command = [
            executable,
            "--threads",
            str(args.threads),
            "--silent",
            "trace",
            "--plasmid",
            str(query),
            "--database",
            str(args.jam_database),
            "--catalog",
            str(args.jam_catalog),
            "--output",
            str(output),
            "--plasmid-id",
            query_id,
            "--sensitivity",
            args.sensitivity,
            "--min-shared",
            str(args.min_shared),
            "--top-candidates",
            str(args.top_candidates),
            "--max-alignments",
            str(args.max_alignments),
        ]
        completed = subprocess.run(command, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if completed.returncode != 0:
            error = completed.stderr.strip() or f"exit {completed.returncode}"
            return [failed_result(item, query_length, error, command, version) for item in assemblies], [command]
        try:
            return jam_result_records(output.read_text(encoding="utf-8"), query_id, query_length, assemblies, truth_path, command, version), [command]
        except (OSError, json.JSONDecodeError, ValueError, NormalizationError) as exc:
            return [failed_result(item, query_length, str(exc), command, version) for item in assemblies], [command]


def unavailable_result(
    tool: str, scope: str, query_id: str, query_length: int, selection: dict, error: str
) -> dict:
    return {
        "schema_version": SCHEMA_VERSION,
        "tool": tool,
        "status": "unavailable",
        "scope": scope,
        "query_id": query_id,
        "query_length": query_length,
        "selection": selection,
        "results": [],
        "metrics": {"truth_status": "unavailable"},
        "error": error,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tool", choices=("jam", "blastn", "minimap2"), required=True)
    parser.add_argument("--query", type=Path, required=True)
    parser.add_argument("--assemblies-dir", type=Path, required=True)
    parser.add_argument("--candidate-manifest", type=Path)
    parser.add_argument("--scope", choices=("candidate", "direct"), default="candidate")
    parser.add_argument("--direct-limit", type=int)
    parser.add_argument("--truth", type=Path)
    parser.add_argument("--query-id")
    parser.add_argument("--query-length", type=int)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--executable")
    parser.add_argument("--jam-output", type=Path)
    parser.add_argument("--jam-database")
    parser.add_argument("--jam-catalog", type=Path)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--sensitivity", choices=("fast", "balanced", "sensitive"), default="balanced")
    parser.add_argument("--min-shared", type=int, default=3)
    parser.add_argument("--top-candidates", type=int, default=250)
    parser.add_argument("--max-alignments", type=int, default=256)
    parser.add_argument("--blast-task", default="megablast")
    parser.add_argument("--evalue", default="1e-10")
    parser.add_argument("--max-hsps", type=int, default=10)
    parser.add_argument("--min-pident", type=float, default=0.0)
    parser.add_argument("--minimap-preset", default="asm5")
    parser.add_argument("--min-mapq", type=int, default=0)
    args = parser.parse_args()
    try:
        query_id, query_length = read_fasta_length(args.query, args.query_id)
        if args.query_id is not None:
            query_id = args.query_id
        if args.query_length is not None and args.query_length != query_length:
            raise ValueError("--query-length disagrees with query FASTA")
        assemblies, selection = select_assemblies(
            args.assemblies_dir, args.candidate_manifest, args.scope, args.direct_limit
        )
        scope = selection["effective_scope"]
        executable = args.executable or shutil.which(args.tool)
        if executable is None:
            result = unavailable_result(
                args.tool,
                scope,
                query_id,
                query_length,
                selection,
                f"executable not found: {args.executable or args.tool}",
            )
            result["input_manifest_sha256"] = _truth_checksum(args.candidate_manifest)
            result["truth_manifest_sha256"] = _truth_checksum(args.truth)
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            return 0
        if args.tool == "jam":
            results, commands = run_jam(executable, args.query, query_id, query_length, assemblies, args.truth, args)
        else:
            results, commands = run_alignment_tool(args.tool, executable, args.query, query_id, query_length, assemblies, args.truth, args)
        failed = sum(item.get("status") == "failed" for item in results)
        measured = sum(item.get("status") == "measured" for item in results)
        status = "measured" if measured else "failed"
        result = {
            "schema_version": SCHEMA_VERSION,
            "tool": args.tool,
            "status": status,
            "scope": scope,
            "query_id": query_id,
            "query_length": query_length,
            "version": command_version(executable),
            "commands": [command_text(command) for command in commands],
            "selection": selection,
            "input_manifest_sha256": _truth_checksum(args.candidate_manifest),
            "truth_manifest_sha256": _truth_checksum(args.truth),
            "results": results,
            "metrics": aggregate_metrics(results, query_length),
            "failed_count": failed,
        }
    except (OSError, ValueError, NormalizationError, subprocess.SubprocessError) as exc:
        parser.error(str(exc))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0 if result["status"] in {"measured", "unavailable"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
