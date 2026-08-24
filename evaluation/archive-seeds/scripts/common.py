#!/usr/bin/env python3
"""Shared safety and measurement helpers for archive-seed experiments.

The runners in this directory deliberately keep the experiment protocol
outside the trace implementation.  A manifest is a list of argument arrays,
not a shell script.  Planning is the default; execution is opt-in and every
case is bounded, writes only below the checkout, and refuses an existing
output.  Child JSON may contain stage/request counters produced by the trace
runner.  The helpers preserve those records while also exposing a stable
measurement envelope for archive, seed, local, remote, and comparator runs.
"""

from __future__ import annotations

import hashlib
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Iterable

ROOT = Path(__file__).resolve().parents[3]
TMP_ROOT = ROOT / "tools" / "out" / "archive-seeds" / "tmp"
MANIFEST_SCHEMA = "archive-seeds-experiment-manifest-v1"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


class ExperimentError(ValueError):
    """Raised when a manifest, path, or result cannot be handled safely."""


def workspace_path(value: str | Path, *, field: str) -> Path:
    """Resolve a path and reject paths outside the jam-rs checkout."""

    path = Path(value)
    if not path.is_absolute():
        path = ROOT / path
    try:
        resolved = path.resolve()
    except OSError as exc:
        raise ExperimentError(f"{field} cannot be resolved: {path}: {exc}") from exc
    if resolved != ROOT and ROOT not in resolved.parents:
        raise ExperimentError(f"{field} must remain below workspace {ROOT}: {value}")
    return resolved


def ensure_new_output(path: str | Path) -> Path:
    """Validate a workspace output without creating or replacing it."""

    output = workspace_path(path, field="output")
    if output.exists():
        raise ExperimentError(f"refusing to overwrite existing output: {output}")
    return output


def load_manifest(path: str | Path, expected_kind: str) -> dict[str, Any]:
    manifest_path = workspace_path(path, field="manifest")
    try:
        value = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ExperimentError(f"cannot read manifest {manifest_path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ExperimentError("experiment manifest must be an object")
    if value.get("schema_version") != MANIFEST_SCHEMA:
        raise ExperimentError(
            f"manifest schema_version must be {MANIFEST_SCHEMA!r}"
        )
    if value.get("kind") != expected_kind:
        raise ExperimentError(
            f"manifest kind {value.get('kind')!r} does not match {expected_kind!r}"
        )
    cases = value.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ExperimentError("experiment manifest requires a non-empty cases list")
    for index, case in enumerate(cases):
        _validate_case(case, index, expected_kind)
    return value


def _validate_case(case: Any, index: int, kind: str) -> None:
    if not isinstance(case, dict):
        raise ExperimentError(f"cases[{index}] must be an object")
    command = case.get("command")
    if (
        not isinstance(command, list)
        or not command
        or not all(isinstance(item, str) and item for item in command)
    ):
        raise ExperimentError(f"cases[{index}].command must be a non-empty string list")
    metadata = case.get("metadata", {})
    if not isinstance(metadata, dict):
        raise ExperimentError(f"cases[{index}].metadata must be an object")
    for key in ("output_paths", "input_paths"):
        paths = metadata.get(key, [])
        if not isinstance(paths, list) or not all(isinstance(item, str) for item in paths):
            raise ExperimentError(f"cases[{index}].metadata.{key} must be a string list")
        for item in paths:
            workspace_path(item, field=f"cases[{index}].metadata.{key}")
    cwd = case.get("cwd")
    if cwd is not None:
        workspace_path(cwd, field=f"cases[{index}].cwd")
    environment = case.get("env", {})
    if not isinstance(environment, dict) or not all(
        isinstance(key, str) and isinstance(value, str)
        for key, value in environment.items()
    ):
        raise ExperimentError(f"cases[{index}].env must map strings to strings")
    supplied_tmpdir = environment.get("TMPDIR")
    if supplied_tmpdir is not None:
        workspace_path(supplied_tmpdir, field=f"cases[{index}].env.TMPDIR")
    timeout = metadata.get("timeout_seconds", 60)
    if isinstance(timeout, bool) or not isinstance(timeout, (int, float)):
        raise ExperimentError(f"cases[{index}].metadata.timeout_seconds must be numeric")
    if timeout <= 0 or timeout > 600:
        raise ExperimentError(
            f"cases[{index}].metadata.timeout_seconds must be in (0, 600]"
        )
    cache_state = metadata.get("cache_state")
    if cache_state is not None and not isinstance(cache_state, str):
        raise ExperimentError(f"cases[{index}].metadata.cache_state must be a string")
    if kind == "local_matrix" and cache_state not in {"process-cold", "process-warm"}:
        raise ExperimentError(
            f"cases[{index}] local matrix requires cache_state process-cold or process-warm"
        )
    if kind == "remote_matrix" and cache_state not in {"cold", "warm"}:
        raise ExperimentError(
            f"cases[{index}] remote matrix requires cache_state cold or warm"
        )


def _command(raw: Any, index: int) -> list[str]:
    if not isinstance(raw, list) or not raw or not all(
        isinstance(item, str) and item for item in raw
    ):
        raise ExperimentError(f"cases[{index}].command must be a non-empty string list")
    return list(raw)


def _parse_json_output(stdout: bytes) -> Any:
    text = stdout.decode("utf-8", errors="replace").strip()
    if not text:
        return None
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        records: list[Any] = []
        for line in text.splitlines():
            if not line.strip():
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError:
                return None
        return records or None


def _as_nonnegative_int(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value if value >= 0 else None
    if isinstance(value, float) and value.is_integer() and value >= 0:
        return int(value)
    return None


def _walk_dicts(value: Any) -> Iterable[dict[str, Any]]:
    if isinstance(value, dict):
        yield value
        for nested in value.values():
            yield from _walk_dicts(nested)
    elif isinstance(value, list):
        for nested in value:
            yield from _walk_dicts(nested)


def _first_counter(value: Any, names: tuple[str, ...]) -> int:
    for item in _walk_dicts(value):
        for name in names:
            number = _as_nonnegative_int(item.get(name))
            if number is not None:
                return number
    return 0


def _records(value: Any, names: tuple[str, ...]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for item in _walk_dicts(value):
        for name in names:
            nested = item.get(name)
            if isinstance(nested, list):
                records.extend(record for record in nested if isinstance(record, dict))
    return records


def _normal_stage(value: dict[str, Any], index: int) -> dict[str, Any]:
    aliases = {
        "stage": ("stage", "round", "stage_id"),
        "name": ("name", "stage_name"),
        "seed_pages_read": ("seed_pages_read", "seed_pages_requested"),
        "keys_tested": ("keys_tested", "seed_keys_tested"),
        "occurrences_decoded": ("occurrences_decoded", "occurrences_after_limits"),
        "anchors_retained": ("anchors_retained", "anchors_after_limits"),
        "chains_retained": ("chains_retained", "chains_accepted"),
        "sequence_blocks_read": ("sequence_blocks_read", "sequence_blocks_fetched"),
        "alignment_attempts": ("alignment_attempts", "alignment_windows_attempted"),
        "new_query_bases_supported": ("new_query_bases_supported",),
        "wall_micros": ("wall_micros",),
        "cpu_micros": ("cpu_micros",),
        "bytes_read": ("bytes_read", "bytes_returned"),
        "bytes_requested": ("bytes_requested",),
        "range_requests": ("range_requests", "coalesced_range_requests"),
        "rss_max_kib": ("rss_max_kib", "peak_rss_kib"),
    }
    result: dict[str, Any] = {}
    for field, field_aliases in aliases.items():
        if field == "name":
            result[field] = next(
                (str(value[name]) for name in field_aliases if name in value),
                f"stage-{index}",
            )
        else:
            result[field] = next(
                (
                    number
                    for name in field_aliases
                    for number in [_as_nonnegative_int(value.get(name))]
                    if number is not None
                ),
                0,
            )
    return result


def _normal_request(value: dict[str, Any], index: int) -> dict[str, Any]:
    aliases = {
        "request_id": ("request_id", "id"),
        "kind": ("kind", "type", "resource"),
        "offset": ("offset", "start"),
        "length": ("length", "bytes_requested"),
        "bytes_requested": ("bytes_requested", "requested_bytes", "length"),
        "bytes_returned": ("bytes_returned", "returned_bytes", "bytes_read"),
        "wall_micros": ("wall_micros",),
        "retry_count": ("retry_count", "retries"),
    }
    result: dict[str, Any] = {"request_id": index}
    for field, field_aliases in aliases.items():
        if field in {"request_id", "kind"}:
            result[field] = next(
                (value[name] for name in field_aliases if name in value),
                result[field],
            ) if field == "request_id" else next(
                (str(value[name]) for name in field_aliases if name in value),
                "unknown",
            )
        else:
            result[field] = next(
                (
                    number
                    for name in field_aliases
                    for number in [_as_nonnegative_int(value.get(name))]
                    if number is not None
                ),
                0,
            )
    result["cache_hit"] = bool(value.get("cache_hit", False))
    return result


def _measurement_envelope(result: Any, *, wall_micros: int, cpu_micros: int, rss_max_kib: int) -> dict[str, Any]:
    stages = [_normal_stage(item, index) for index, item in enumerate(_records(result, ("stages", "stage_metrics")))]
    requests = [_normal_request(item, index) for index, item in enumerate(_records(result, ("requests", "request_log", "range_requests")))]
    metric_aliases = {
        "archive_bytes": ("archive_bytes", "object_bytes", "archive_size_bytes"),
        "metadata_bytes_read": ("metadata_bytes_read",),
        "metadata_requests": ("metadata_requests",),
        "seed_bytes_read": ("seed_bytes_read",),
        "sequence_bytes_read": ("sequence_bytes_read",),
        "bytes_read": ("bytes_read", "returned_bytes", "bytes_returned"),
        "bytes_requested": ("bytes_requested", "requested_bytes"),
        "decoded_bytes": ("decoded_bytes",),
        "decoded_bases": ("decoded_bases", "decoded_sequence_bases"),
        "mapped_bytes": ("mapped_bytes",),
        "resident_bytes": ("resident_bytes",),
        "remote_bytes": ("remote_bytes",),
        "range_requests": ("range_requests", "request_count", "requests"),
        "stream_requests": ("stream_requests",),
        "coalesced_range_requests": ("coalesced_range_requests",),
        "full_object_requests": (
            "full_object_requests",
            "full_object_fallbacks",
            "full_downloads",
        ),
        "head_requests": ("head_requests",),
        "get_requests": ("get_requests",),
        "retries": ("retries", "retry_count"),
        "cache_bytes": ("cache_bytes",),
        "cache_hits": ("cache_hits",),
        "cache_misses": ("cache_misses",),
        "seed_pages_read": ("seed_pages_read", "seed_buckets_read"),
        "sequence_blocks_read": ("sequence_blocks_read",),
    }
    metrics = {
        name: _first_counter(result, aliases) for name, aliases in metric_aliases.items()
    }
    if not metrics["range_requests"] and requests:
        metrics["range_requests"] = len(requests)
    if not metrics["bytes_requested"] and requests:
        metrics["bytes_requested"] = sum(item["bytes_requested"] for item in requests)
    if not metrics["bytes_read"] and requests:
        metrics["bytes_read"] = sum(item["bytes_returned"] for item in requests)
    return {
        "timing": {
            "wall_micros": wall_micros,
            "cpu_micros": cpu_micros,
        },
        "rss_max_kib": rss_max_kib,
        "metrics": metrics,
        "stages": stages,
        "requests": requests,
    }


def _filesystem_metrics(metadata: dict[str, Any]) -> dict[str, Any]:
    """Add source/output and JMA section sizes when an executed case produced them."""

    input_bytes = 0
    output_bytes = 0
    archive_bytes = 0
    seed_index_bytes = 0
    sequence_directory_bytes = 0
    sequence_payload_bytes = 0
    contig_table_bytes = 0
    seed_index_bytes_by_scheme: dict[str, int] = {}
    for raw_path in metadata.get("input_paths", []):
        path = workspace_path(raw_path, field="measurement input")
        if path.is_file():
            input_bytes += path.stat().st_size
    for raw_path in metadata.get("output_paths", []):
        path = workspace_path(raw_path, field="measurement output")
        if not path.is_file():
            continue
        size = path.stat().st_size
        output_bytes += size
        if path.suffix != ".jma":
            continue
        archive_bytes += size
        try:
            seed_section: tuple[int, int] | None = None
            with path.open("rb") as handle:
                header = handle.read(256)
                if len(header) != 256 or header[:8] != b"JMAF1\0\0\0":
                    continue
                section_count = int.from_bytes(header[20:24], "little")
                directory_offset = int.from_bytes(header[24:32], "little")
                directory_length = int.from_bytes(header[32:40], "little")
                if section_count > 1_000_000 or directory_length != section_count * 64:
                    continue
                if directory_offset + directory_length > size:
                    continue
                handle.seek(directory_offset)
                directory = handle.read(directory_length)
            for offset in range(0, len(directory), 64):
                kind = int.from_bytes(directory[offset : offset + 4], "little")
                section_offset = int.from_bytes(directory[offset + 8 : offset + 16], "little")
                section_length = int.from_bytes(directory[offset + 16 : offset + 24], "little")
                if kind == 2:
                    contig_table_bytes += section_length
                elif kind == 4:
                    seed_index_bytes += section_length
                    seed_section = (section_offset, section_length)
                elif kind == 5:
                    sequence_directory_bytes += section_length
                elif kind == 6:
                    sequence_payload_bytes += section_length
            if seed_section is not None:
                section_offset, section_length = seed_section
                with path.open("rb") as handle:
                    if section_length < 40 or section_offset + section_length > size:
                        continue
                    handle.seek(section_offset)
                    seed_header = handle.read(40)
                    scheme_count = int.from_bytes(seed_header[4:8], "little")
                    page_count = int.from_bytes(seed_header[8:12], "little")
                    scheme_offset = int.from_bytes(seed_header[16:24], "little")
                    page_offset = int.from_bytes(seed_header[24:32], "little")
                    data_offset = int.from_bytes(seed_header[32:40], "little")
                    if scheme_count > 1_000_000 or page_count > 1_000_000:
                        continue
                    if scheme_offset + scheme_count * 64 > section_length:
                        continue
                    if page_offset + page_count * 104 > section_length:
                        continue
                    if data_offset > section_length:
                        continue
                    handle.seek(section_offset + scheme_offset)
                    schemes = handle.read(scheme_count * 64)
                    handle.seek(section_offset + page_offset)
                    pages = handle.read(page_count * 104)
                for scheme_index in range(scheme_count):
                    scheme_record = schemes[scheme_index * 64 : (scheme_index + 1) * 64]
                    scheme_id = int.from_bytes(scheme_record[0:4], "little")
                    first_page = int.from_bytes(scheme_record[24:28], "little")
                    scheme_page_count = int.from_bytes(scheme_record[28:32], "little")
                    if first_page + scheme_page_count > page_count:
                        continue
                    scheme_size = 64 + scheme_page_count * 104
                    for page_index in range(first_page, first_page + scheme_page_count):
                        page_record = pages[page_index * 104 : (page_index + 1) * 104]
                        key_offset = int.from_bytes(page_record[24:32], "little")
                        key_length = int.from_bytes(page_record[32:40], "little")
                        occurrence_offset = int.from_bytes(page_record[40:48], "little")
                        occurrence_length = int.from_bytes(page_record[48:56], "little")
                        if (
                            key_offset + key_length > section_length
                            or occurrence_offset + occurrence_length > section_length
                        ):
                            continue
                        scheme_size += key_length + occurrence_length
                    key = str(scheme_id)
                    seed_index_bytes_by_scheme[key] = (
                        seed_index_bytes_by_scheme.get(key, 0) + scheme_size
                    )
        except OSError:
            continue
    return {
        "source_input_bytes": input_bytes,
        "output_bytes": output_bytes,
        "archive_bytes": archive_bytes,
        "seed_index_bytes": seed_index_bytes,
        "seed_index_bytes_by_scheme": seed_index_bytes_by_scheme,
        "sequence_directory_bytes": sequence_directory_bytes,
        "sequence_payload_bytes": sequence_payload_bytes,
        "contig_table_bytes": contig_table_bytes,
    }


def _blast_json_alignments(value: Any) -> Any:
    """Extract interval-only records from BLAST JSON2 output."""

    alignments: list[dict[str, Any]] = []
    for item in _walk_dicts(value):
        hsps = item.get("hsps")
        if not isinstance(hsps, list):
            continue
        for hsp in hsps:
            if not isinstance(hsp, dict):
                continue
            start = hsp.get("query_from")
            end = hsp.get("query_to")
            if isinstance(start, int) and isinstance(end, int):
                low, high = sorted((start, end))
                record: dict[str, Any] = {
                    "query_segments": [{"start": low - 1, "end": high}]
                }
                query_sequence = hsp.get("qseq")
                target_sequence = hsp.get("hseq")
                if isinstance(query_sequence, str) and isinstance(target_sequence, str):
                    runs: list[dict[str, Any]] = []
                    for query_base, target_base in zip(query_sequence, target_sequence):
                        if query_base == "-":
                            operation = "deletion"
                        elif target_base == "-":
                            operation = "insertion"
                        elif query_base.upper() == target_base.upper():
                            operation = "equal"
                        else:
                            operation = "substitution"
                        if runs and runs[-1]["operation"] == operation:
                            runs[-1]["length"] += 1
                        else:
                            runs.append({"operation": operation, "length": 1})
                    record["edit_script"] = runs
                alignments.append(record)
    return {"alignments": alignments}


def _paf_alignments(text: str) -> Any:
    alignments: list[dict[str, Any]] = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        try:
            start, end = int(fields[2]), int(fields[3])
        except ValueError:
            continue
        record: dict[str, Any] = {"query_segments": [{"start": start, "end": end}]}
        for tag in fields[12:]:
            if tag.startswith("cg:Z:"):
                record["cigar"] = tag[5:]
                break
        alignments.append(record)
    return {"alignments": alignments}


def _normalize_comparator_result(result: Any, metadata: dict[str, Any]) -> Any:
    normalization = metadata.get("normalization")
    if not isinstance(normalization, dict):
        return result
    from tools.trace_failure_analysis.comparators import ComparatorError, normalize_comparator
    from tools.trace_failure_analysis.intervals import IntervalError, parse_intervals

    query_length = normalization.get("query_length")
    if isinstance(query_length, bool) or not isinstance(query_length, int) or query_length <= 0:
        raise ExperimentError("comparator normalization requires a positive query_length")
    format_name = str(normalization.get("format", "json"))
    if format_name == "blast-json" and not isinstance(result, (dict, list)):
        result = _blast_json_alignments(result)
    elif format_name == "blast-json":
        blast_records = _blast_json_alignments(result)
        if blast_records["alignments"]:
            result = blast_records
    elif format_name == "paf":
        result = _paf_alignments(str(result or ""))
    truth_value: Any = normalization.get("truth_intervals", [])
    truth_path = normalization.get("truth_path")
    if truth_path is not None:
        truth_file = workspace_path(truth_path, field="normalization.truth_path")
        try:
            truth_value = json.loads(truth_file.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ExperimentError(f"cannot read comparator truth {truth_file}: {exc}") from exc
    try:
        truth = parse_intervals(truth_value, field="truth_intervals")
        normalized = normalize_comparator(
            result,
            query_length,
            truth_intervals=truth,
            query_deletion_operation=str(
                normalization.get("query_deletion_operation", "I")
            ),
        )
    except (ComparatorError, IntervalError) as exc:
        raise ExperimentError(f"comparator normalization failed: {exc}") from exc
    target = metadata.get("target_recall")
    if isinstance(target, (int, float)) and not isinstance(target, bool):
        recall = normalized["metrics"].get("base_recall")
        normalized["matched_recall"] = {
            "target": float(target),
            "observed": recall,
            "within_tolerance": recall is not None and abs(float(recall) - float(target)) <= float(metadata.get("recall_tolerance", 0.01)),
        }
    return normalized


def _prepare_case_outputs(case: dict[str, Any], index: int) -> None:
    metadata = case.get("metadata", {})
    for item in metadata.get("output_paths", []):
        output = workspace_path(item, field=f"cases[{index}].metadata.output_paths")
        if output.exists():
            raise ExperimentError(f"refusing to overwrite existing case output: {output}")


def _attach_measurement(row: dict[str, Any], measurement: dict[str, Any]) -> None:
    """Expose stable nested counters and convenient row-level aliases."""

    row["measurement"] = measurement
    row["stage_metrics"] = measurement["stages"]
    row["request_metrics"] = measurement["requests"]
    row["rss_max_kib"] = measurement["rss_max_kib"]
    row["wall_micros"] = measurement["timing"]["wall_micros"]
    row["cpu_micros"] = measurement["timing"]["cpu_micros"]


def run_manifest(
    manifest: dict[str, Any],
    *,
    execute: bool,
    max_case_seconds: float = 60.0,
) -> dict[str, Any]:
    """Plan or execute bounded cases and return the stable measurement envelope."""

    if max_case_seconds <= 0 or max_case_seconds > 600:
        raise ExperimentError("max_case_seconds must be in (0, 600]")
    rows: list[dict[str, Any]] = []
    manifest_name = str(manifest.get("name", "unnamed"))
    for index, case in enumerate(manifest["cases"]):
        command = _command(case.get("command"), index)
        metadata = case.get("metadata", {})
        _prepare_case_outputs(case, index)
        label = str(case.get("name", f"case-{index}"))
        cache_state = metadata.get("cache_state")
        row: dict[str, Any] = {
            "name": label,
            "command": command,
            "metadata": metadata,
            "cache_state": cache_state,
            "cache_policy": metadata.get(
                "cache_policy",
                "page-cache-state-uncontrolled" if cache_state else None,
            ),
            "status": "planned",
            "measurement": _measurement_envelope(
                None, wall_micros=0, cpu_micros=0, rss_max_kib=0
            ),
        }
        _attach_measurement(row, row["measurement"])
        if not execute:
            rows.append(row)
            continue
        cwd_value = case.get("cwd")
        cwd = workspace_path(cwd_value, field=f"cases[{index}].cwd") if cwd_value else ROOT
        environment = os.environ.copy()
        supplied_environment = case.get("env", {})
        environment.update(supplied_environment)
        case_tmp = TMP_ROOT / manifest_name / f"case-{index}-{label}"
        case_tmp = workspace_path(case_tmp, field="case TMPDIR")
        case_tmp.mkdir(parents=True, exist_ok=True)
        environment["TMPDIR"] = str(case_tmp)
        timeout = min(float(metadata.get("timeout_seconds", max_case_seconds)), max_case_seconds)
        before_cpu = _child_cpu_seconds()
        before_rss = _child_rss_kib()
        started = time.perf_counter_ns()
        try:
            completed = subprocess.run(
                command,
                cwd=cwd,
                env=environment,
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired as exc:
            elapsed = (time.perf_counter_ns() - started) // 1_000
            row.update(
                {
                    "status": "timeout",
                    "returncode": None,
                    "error": f"case exceeded {timeout:g}s",
                    "measurement": _measurement_envelope(
                        None,
                        wall_micros=elapsed,
                        cpu_micros=round(max(0.0, _child_cpu_seconds() - before_cpu) * 1_000_000),
                        rss_max_kib=max(0, _child_rss_kib(), before_rss),
                    ),
                }
            )
            _attach_measurement(row, row["measurement"])
            rows.append(row)
            continue
        except OSError as exc:
            row.update({"status": "unavailable", "error": str(exc), "returncode": None})
            rows.append(row)
            continue
        elapsed = (time.perf_counter_ns() - started) // 1_000
        cpu_micros = round(max(0.0, _child_cpu_seconds() - before_cpu) * 1_000_000)
        rss_max_kib = max(0, _child_rss_kib(), before_rss)
        parsed = _parse_json_output(completed.stdout)
        if parsed is None and metadata.get("comparator_format") == "paf":
            parsed = completed.stdout.decode("utf-8", errors="replace")
        try:
            parsed = _normalize_comparator_result(parsed, metadata)
        except ExperimentError as exc:
            row.update({"status": "failed", "error": str(exc), "returncode": completed.returncode})
            rows.append(row)
            continue
        row.update(
            {
                "status": "complete" if completed.returncode == 0 else "failed",
                "returncode": completed.returncode,
                "wall_micros": elapsed,
                "cpu_micros": cpu_micros,
                "stdout_bytes": len(completed.stdout),
                "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(),
                "stderr_bytes": len(completed.stderr),
                "stderr_tail": completed.stderr.decode("utf-8", errors="replace")[-2000:],
                "result": parsed,
                "measurement": _measurement_envelope(
                    parsed,
                    wall_micros=elapsed,
                    cpu_micros=cpu_micros,
                    rss_max_kib=rss_max_kib,
                ),
            }
        )
        _attach_measurement(row, row["measurement"])
        row["measurement"]["metrics"].update(_filesystem_metrics(metadata))
        if metadata.get("range_only_required"):
            full_reads = row["measurement"]["metrics"].get("full_object_requests", 0)
            has_io_observation = any(
                row["measurement"]["metrics"].get(name, 0) > 0
                for name in ("range_requests", "bytes_requested", "bytes_read")
            )
            row["selective_io_check"] = {
                "required": True,
                "full_object_requests": full_reads,
                "status": "fail" if full_reads else ("pass" if has_io_observation else "not-recorded"),
            }
        rows.append(row)
    return {
        "schema_version": "archive-seeds-experiment-v1",
        "kind": manifest["kind"],
        "mode": "execute" if execute else "plan",
        "manifest": manifest_name,
        "workspace": str(ROOT),
        "tmp_policy": "workspace-only",
        "max_case_seconds": max_case_seconds,
        "cases": rows,
    }


def _child_cpu_seconds() -> float:
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    return float(usage.ru_utime + usage.ru_stime)


def _child_rss_kib() -> int:
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    # Linux reports KiB; macOS reports bytes.  jam-rs experiments run on
    # Linux, but retaining this branch keeps the field honest elsewhere.
    value = int(usage.ru_maxrss)
    return value if os.name == "posix" else value // 1024


def write_new(path: str | Path, value: Any, *, pretty: bool = False) -> None:
    output = ensure_new_output(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(value, indent=2 if pretty else None, sort_keys=True) + "\n"
    try:
        with output.open("x", encoding="utf-8") as handle:
            handle.write(text)
    except FileExistsError as exc:
        raise ExperimentError(f"refusing to overwrite existing output: {output}") from exc
