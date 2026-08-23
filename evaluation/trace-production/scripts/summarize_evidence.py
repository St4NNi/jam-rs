#!/usr/bin/env python3
"""Build a deterministic evidence summary from the production raw tree.

The raw tree is treated as an immutable measurement manifest.  This script
validates every JSON, TSV, and JSONL input before deriving any rows, hashes all
raw files, and writes only to a new output directory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import statistics
import sys
from pathlib import Path
from typing import Any, Iterable


SCHEMA_VERSION = "evidence-summary-v1"
TSV_FIELDS = ("section", "case", "metric", "value", "unit", "status", "note")
EXPECTED_FRAGMENT_CASES = (
    "jam_fragments_fast_s1_t8.json",
    "jam_fragments_fast_s3_t8.json",
    "jam_fragments_balanced_s1_t8.json",
    "jam_fragments_balanced_s3_t8.json",
    "jam_fragments_sensitive_s1_t8.json",
    "jam_fragments_sensitive_s3_t8.json",
)
ACCESSION_QUERIES = ("AP000342.1", "J01749.1", "NC_008272.1")
ACCESSION_FILES = (
    "environment.json",
    "generation.json",
    "local_balanced_t1_cold.json",
    "local_balanced_t8_cold.json",
    "local_balanced_t8_cold.jsonl",
    "summary.tsv",
)
CONTROLLED_NORMALIZED = (
    "jam_fragments_fast_s1_t8.json",
    "jam_fragments_fast_s3_t8.json",
    "jam_fragments_balanced_s1_t8.json",
    "jam_fragments_balanced_s3_t8.json",
    "jam_fragments_sensitive_s1_t8.json",
    "jam_fragments_sensitive_s3_t8.json",
    "jam-q2000-identity.json",
    "jam-q10000-identity.json",
    "jam-q50000-identity.json",
    "jam-q10000-origin.json",
    "jam-q10000-fragments.json",
    "blastn-q10000-candidates-t8.json",
    "blastn-q10000-fragments-direct-t8.json",
    "minimap2-q10000-candidates-t8.json",
    "minimap2-q10000-fragments-direct-t8.json",
)
REQUIRED_STATIC = (
    "controlled-v1/input/controlled-source.json",
    "controlled-v1/input/fragments-indels-manifest.json",
    "controlled-v1/input/identity-coverage-manifest.json",
    "controlled-v1/input/origin-manifest.json",
    "controlled-v1/truth/fragments-indels.tsv",
    "controlled-v1/truth/identity-coverage.tsv",
    "controlled-v1/truth/origin.tsv",
    "resource-scaling-v1/block/environment.json",
    "resource-scaling-v1/block/input_manifest.json",
    "resource-scaling-v1/block/summary.json",
    "resource-scaling-v1/block/summary.tsv",
    "resource-scaling-v1/memory/environment.json",
    "resource-scaling-v1/memory/input_manifest.json",
    "resource-scaling-v1/memory/summary.json",
    "resource-scaling-v1/memory/summary.tsv",
    "resource-scaling-v1/http/cold.measurement.json",
    "resource-scaling-v1/http/warm.measurement.json",
    "resource-scaling-v1/http/http.metrics.json",
    "accession-spike-v1/source/dataset.json",
    "accession-spike-v1/source/source_verification.tsv",
    "accession-spike-v1/source/truth.tsv",
)


JsonValue = dict[str, Any] | list[Any]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite_number(value: Any) -> int | float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    if not math.isfinite(float(value)):
        return None
    return value


def median(values: Iterable[Any]) -> int | float | None:
    numbers = [finite_number(value) for value in values]
    numbers = [value for value in numbers if value is not None]
    return statistics.median(numbers) if numbers else None


def maximum(values: Iterable[Any]) -> int | float | None:
    numbers = [finite_number(value) for value in values]
    numbers = [value for value in numbers if value is not None]
    return max(numbers) if numbers else None


def relative(path: Path, root: Path) -> str:
    return path.relative_to(root).as_posix()


def read_json(path: Path) -> JsonValue:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{path}: invalid JSON: {error}") from error
    if not isinstance(value, (dict, list)):
        raise ValueError(f"{path}: JSON root must be an object or array")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    csv.field_size_limit(1024 * 1024 * 1024)
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            try:
                header = next(reader)
            except StopIteration as error:
                raise ValueError(f"{path}: empty TSV") from error
            if not header or any(not field for field in header):
                raise ValueError(f"{path}: TSV header contains an empty field")
            if len(set(header)) != len(header):
                raise ValueError(f"{path}: TSV header contains duplicate fields")
            rows: list[dict[str, str]] = []
            for line_number, row in enumerate(reader, 2):
                if len(row) != len(header):
                    raise ValueError(
                        f"{path}:{line_number}: expected {len(header)} fields, got {len(row)}"
                    )
                rows.append(dict(zip(header, row)))
            return rows
    except csv.Error as error:
        raise ValueError(f"{path}: invalid TSV: {error}") from error
    except OSError as error:
        raise ValueError(f"{path}: cannot read TSV: {error}") from error


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    try:
        with path.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip():
                    raise ValueError(f"{path}:{line_number}: blank JSONL line")
                try:
                    value = json.loads(line)
                except json.JSONDecodeError as error:
                    raise ValueError(f"{path}:{line_number}: invalid JSON: {error}") from error
                if not isinstance(value, dict):
                    raise ValueError(f"{path}:{line_number}: JSONL record is not an object")
                records.append(value)
    except OSError as error:
        raise ValueError(f"{path}: cannot read JSONL: {error}") from error
    if not records:
        raise ValueError(f"{path}: empty JSONL")
    if records[0].get("record_type") != "run_header":
        raise ValueError(f"{path}: first record is not run_header")
    if records[-1].get("record_type") != "run_footer":
        raise ValueError(f"{path}: last record is not run_footer")
    return records


def validate_inputs(raw_root: Path) -> tuple[dict[str, JsonValue], dict[str, list[dict[str, str]]], dict[str, list[dict[str, Any]]], list[str]]:
    json_values: dict[str, JsonValue] = {}
    tsv_values: dict[str, list[dict[str, str]]] = {}
    jsonl_values: dict[str, list[dict[str, Any]]] = {}
    errors: list[str] = []

    required = set(REQUIRED_STATIC)
    required.update(
        f"controlled-v1/measurements/{name}" for name in EXPECTED_FRAGMENT_CASES
    )
    required.update(
        f"controlled-v1/normalized/{name}" for name in CONTROLLED_NORMALIZED
    )
    for thread in (1, 4, 8, 16):
        required.add(f"resource-scaling-v1/large-catalog/jam_q50000_sensitive_c100_t{thread}.json")
    for query in ACCESSION_QUERIES:
        required.update(f"accession-spike-v1/{query}/{name}" for name in ACCESSION_FILES)

    for rel in sorted(required):
        if not (raw_root / rel).is_file():
            errors.append(f"missing expected raw input: {rel}")

    for path in sorted(raw_root.rglob("*")):
        if not path.is_file():
            continue
        rel = relative(path, raw_root)
        try:
            if path.suffix == ".json":
                json_values[rel] = read_json(path)
            elif path.suffix == ".tsv":
                tsv_values[rel] = read_tsv(path)
            elif path.suffix == ".jsonl":
                jsonl_values[rel] = read_jsonl(path)
        except ValueError as error:
            errors.append(str(error))

    return json_values, tsv_values, jsonl_values, errors


def raw_status(value: dict[str, Any]) -> str:
    explicit = value.get("status")
    if explicit in {"failed", "unavailable"}:
        return str(explicit)
    if value.get("exit_code") not in (None, 0):
        return "failed"
    return "measured"


def metric_value(value: Any) -> Any:
    if value is None:
        return None
    number = finite_number(value)
    if number is not None:
        return number
    if isinstance(value, (str, bool)):
        return value
    return None


def ratio(value: Any) -> tuple[Any, str]:
    if isinstance(value, str):
        match = re.fullmatch(r"\s*(\d+)\s*/\s*(\d+)\s*", value)
        if match:
            numerator, denominator = (int(part) for part in match.groups())
            if denominator == 0:
                return "NA", f"reported recovery {value}; denominator is zero"
            return numerator / denominator, f"reported recovery {value}"
    return metric_value(value), "reported by the raw measurement"


def source_note(note: str, source_files: Iterable[str]) -> str:
    sources = ", ".join(sorted(set(source_files)))
    return f"{note}; source: {sources}" if sources else note


def add_row(
    rows: list[dict[str, Any]],
    section: str,
    case: str,
    metric: str,
    value: Any,
    unit: str,
    status: str,
    note: str,
    source_files: Iterable[str] = (),
) -> None:
    normalized = metric_value(value)
    if status != "measured":
        normalized = "NA"
    elif normalized is None or normalized == "NA":
        status = "unavailable"
        normalized = "NA"
    rows.append(
        {
            "section": section,
            "case": case,
            "metric": metric,
            "value": normalized,
            "unit": unit,
            "status": status,
            "note": source_note(note, source_files),
            "source_files": sorted(set(source_files)),
        }
    )


def add_measurement_fields(
    rows: list[dict[str, Any]],
    section: str,
    case: str,
    value: dict[str, Any],
    source: str,
    fields: Iterable[tuple[str, str]],
    note: str,
) -> None:
    status = raw_status(value)
    for field, unit in fields:
        add_row(rows, section, case, field, value.get(field), unit, status, note, [source])


def nested(value: dict[str, Any], parent: str, child: str) -> Any:
    container = value.get(parent)
    return container.get(child) if isinstance(container, dict) else None


def load_object(values: dict[str, JsonValue], rel: str) -> dict[str, Any]:
    value = values.get(rel)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {rel}")
    return value


def summarize_controlled(
    rows: list[dict[str, Any]],
    json_values: dict[str, JsonValue],
    tsv_values: dict[str, list[dict[str, str]]],
) -> None:
    normalized_root = "controlled-v1/normalized"
    for name in CONTROLLED_NORMALIZED:
        rel = f"{normalized_root}/{name}"
        value = load_object(json_values, rel)
        stem = Path(name).stem
        if "identity" in stem:
            section = "controlled_identity"
        elif "origin" in stem:
            section = "controlled_origin"
        elif "fragments" in stem:
            section = "controlled_fragments"
        else:
            section = "controlled_comparator"
        case = stem
        status = str(value.get("status", "measured"))
        if status not in {"measured", "failed", "unavailable"}:
            status = "measured"
        metrics = value.get("metrics") if isinstance(value.get("metrics"), dict) else {}
        for field, unit in (
            ("assemblies_total", "assemblies"),
            ("assemblies_with_records", "assemblies"),
            ("observed_bases", "bases"),
            ("truth_bases", "bases"),
            ("intersection_bases", "bases"),
            ("base_precision", "fraction"),
            ("base_recall", "fraction"),
            ("interval_precision", "fraction"),
            ("interval_recall", "fraction"),
            ("gap_error_bases", "bases"),
            ("gap_error_fraction", "fraction"),
        ):
            add_row(
                rows,
                section,
                case,
                field,
                metrics.get(field),
                unit,
                status,
                "normalized truth comparison metric",
                [rel],
            )
        if value.get("tool"):
            add_row(rows, section, case, "tool", value["tool"], "identifier", status, "normalized result tool", [rel])
        if value.get("version"):
            add_row(rows, section, case, "version", value["version"], "identifier", status, "normalized result version", [rel])

    measurement_root = "controlled-v1/measurements"
    for name in EXPECTED_FRAGMENT_CASES:
        rel = f"{measurement_root}/{name}"
        value = load_object(json_values, rel)
        match = re.fullmatch(r"jam_fragments_(fast|balanced|sensitive)_s(\d+)_t8", Path(name).stem)
        if not match:
            continue
        profile, min_shared = match.groups()
        case = f"profile={profile};min_shared={min_shared};threads=8"
        add_measurement_fields(
            rows,
            "controlled_fragment_measurements",
            case,
            value,
            rel,
            (
                ("wall_seconds", "seconds"),
                ("user_seconds", "seconds"),
                ("system_seconds", "seconds"),
                ("peak_rss_kib", "KiB"),
                ("candidate_count", "candidates"),
                ("alignment_count", "alignments"),
                ("coverage_mean", "fraction"),
                ("gap_bases", "bases"),
            ),
            "Jam fragment profile measurement",
        )

    comparator_raw = {
        ("blastn", "candidate"): "controlled-v1/measurements/blastn_q10000_candidates_t8.json",
        ("blastn", "direct_full"): "controlled-v1/measurements/blastn_q10000_fragments_direct_t8.json",
        ("minimap2", "candidate"): "controlled-v1/measurements/minimap2_q10000_candidates_t8.json",
        ("minimap2", "direct_full"): "controlled-v1/measurements/minimap2_q10000_fragments_direct_t8.json",
    }
    comparator_normalized = {
        ("blastn", "candidate"): "controlled-v1/normalized/blastn-q10000-candidates-t8.json",
        ("blastn", "direct_full"): "controlled-v1/normalized/blastn-q10000-fragments-direct-t8.json",
        ("minimap2", "candidate"): "controlled-v1/normalized/minimap2-q10000-candidates-t8.json",
        ("minimap2", "direct_full"): "controlled-v1/normalized/minimap2-q10000-fragments-direct-t8.json",
    }
    for (tool, scope), normalized_rel in sorted(comparator_normalized.items()):
        normalized = load_object(json_values, normalized_rel)
        raw_rel = comparator_raw[(tool, scope)]
        raw = load_object(json_values, raw_rel)
        case = f"{tool}/{scope}"
        status = str(normalized.get("status", raw_status(raw)))
        if status not in {"measured", "failed", "unavailable"}:
            status = raw_status(raw)
        metrics = normalized.get("metrics") if isinstance(normalized.get("metrics"), dict) else {}
        for field, unit in (
            ("assemblies_total", "assemblies"),
            ("assemblies_with_records", "assemblies"),
            ("observed_bases", "bases"),
            ("truth_bases", "bases"),
            ("intersection_bases", "bases"),
            ("base_precision", "fraction"),
            ("base_recall", "fraction"),
            ("interval_precision", "fraction"),
            ("interval_recall", "fraction"),
            ("gap_error_bases", "bases"),
            ("gap_error_fraction", "fraction"),
        ):
            add_row(rows, "comparator_accuracy", case, field, metrics.get(field), unit, status, "normalized comparator metric", [normalized_rel])
        add_measurement_fields(
            rows,
            "comparator_timing",
            case,
            raw,
            raw_rel,
            (("wall_seconds", "seconds"), ("user_seconds", "seconds"), ("system_seconds", "seconds"), ("peak_rss_kib", "KiB")),
            "comparator measurement timing and resident memory",
        )
        if normalized.get("version"):
            add_row(rows, "comparator_accuracy", case, "version", normalized["version"], "identifier", status, "comparator version from normalized result", [normalized_rel])


def summarize_accession(
    rows: list[dict[str, Any]],
    json_values: dict[str, JsonValue],
    jsonl_values: dict[str, list[dict[str, Any]]],
    tsv_values: dict[str, list[dict[str, str]]],
) -> None:
    truth_rel = "accession-spike-v1/source/truth.tsv"
    truth_rows = tsv_values.get(truth_rel, [])
    for query in ACCESSION_QUERIES:
        query_rows = [item for item in truth_rows if item.get("expected_reference") == query]
        query_rel = f"accession-spike-v1/{query}"
        query_truth_classes = sorted({item.get("truth_class", "") for item in query_rows if item.get("truth_class")})
        add_row(
            rows,
            "accession_spike",
            query,
            "truth_rows",
            len(query_rows),
            "rows",
            "measured" if query_rows else "unavailable",
            f"truth classes: {', '.join(query_truth_classes) if query_truth_classes else 'not reported'}",
            [truth_rel],
        )
        for measurement_name in ("local_balanced_t1_cold.json", "local_balanced_t8_cold.json"):
            rel = f"{query_rel}/{measurement_name}"
            value = load_object(json_values, rel)
            case = f"{query}/threads={value.get('threads', 'NA')};cache={value.get('cache_state', 'NA')}"
            add_measurement_fields(
                rows,
                "accession_spike",
                case,
                value,
                rel,
                (
                    ("wall_seconds", "seconds"),
                    ("user_seconds", "seconds"),
                    ("system_seconds", "seconds"),
                    ("peak_rss_kib", "KiB"),
                    ("candidate_count", "candidates"),
                    ("alignment_count", "alignments"),
                    ("coverage_mean", "fraction"),
                    ("gap_bases", "bases"),
                    ("largest_gap", "bases"),
                ),
                "accession-spike Jam screening measurement",
            )
            for field in (
                "exact_candidate_recovery",
                "exact_alignment_recovery",
                "reverse_complement_candidate_recovery",
                "reverse_complement_alignment_recovery",
                "partial_candidate_recovery",
                "partial_alignment_recovery",
                "split_overlap_candidate_recovery",
                "split_overlap_alignment_recovery",
                "short_fragment_candidate_recovery",
                "short_fragment_alignment_recovery",
            ):
                value_out, note = ratio(value.get(field))
                add_row(rows, "accession_spike", case, field, value_out, "fraction", raw_status(value), note, [rel])
            resource = value.get("resource_metrics") if isinstance(value.get("resource_metrics"), dict) else {}
            for field, unit in (
                ("metadata_requests", "requests"),
                ("range_requests", "requests"),
                ("requested_bytes", "bytes"),
                ("returned_bytes", "bytes"),
                ("decoded_bytes", "bytes"),
                ("cache_bytes", "bytes"),
                ("cache_hits", "requests"),
                ("cache_misses", "requests"),
            ):
                add_row(rows, "accession_spike", case, field, resource.get(field), unit, raw_status(value), "resource metric", [rel])

        jsonl_rel = f"{query_rel}/local_balanced_t8_cold.jsonl"
        records = jsonl_values.get(jsonl_rel, [])
        add_row(rows, "accession_spike_output", query, "jsonl_records", len(records), "records", "measured" if records else "unavailable", "validated run JSONL records", [jsonl_rel])
        add_row(rows, "accession_spike_output", query, "jsonl_result_records", sum(item.get("record_type") == "metagenome_result" for item in records), "records", "measured" if records else "unavailable", "validated metagenome result records", [jsonl_rel])

        summary_rel = f"{query_rel}/summary.tsv"
        for measurement_file, metric, note in (
            ("jam_index.json", "jam_index_wall_seconds", "Jam sketch index timing"),
            ("jma_index.json", "jma_index_wall_seconds", "JMA archive index timing"),
        ):
            summary_rows = [item for item in tsv_values.get(summary_rel, []) if item.get("measurement_file") == measurement_file]
            value = summary_rows[0].get("wall_seconds") if summary_rows else None
            value = float(value) if value not in (None, "", "NA") else None
            add_row(rows, "accession_spike_index", query, metric, value, "seconds", "measured" if summary_rows else "unavailable", note, [summary_rel])


def summarize_block_and_memory(
    rows: list[dict[str, Any]],
    json_values: dict[str, JsonValue],
) -> None:
    block_rel = "resource-scaling-v1/block/summary.json"
    block_obj = load_object(json_values, block_rel)
    block_rows = block_obj.get("rows") if isinstance(block_obj.get("rows"), list) else []
    groups: dict[tuple[Any, Any], list[dict[str, Any]]] = {}
    for item in block_rows:
        if isinstance(item, dict):
            groups.setdefault((item.get("archive_block_bases"), item.get("cache_block_bytes")), []).append(item)
    for (archive_block, cache_block), group in sorted(groups.items(), key=lambda pair: (str(pair[0][0]), str(pair[0][1]))):
        case = f"archive_block={archive_block};cache_block={cache_block}"
        status = "measured" if any(raw_status(item) == "measured" for item in group) else "failed"
        for field, unit in (
            ("wall_seconds", "seconds"),
            ("archive_build_wall_seconds", "seconds"),
            ("archive_build_cpu_seconds", "seconds"),
            ("peak_rss_kib", "KiB"),
            ("archive_build_peak_rss_kib", "KiB"),
            ("archive_size_bytes", "bytes"),
            ("cache_bytes", "bytes"),
            ("range_requests", "requests"),
            ("requested_bytes", "bytes"),
            ("returned_bytes", "bytes"),
            ("decoded_bytes", "bytes"),
        ):
            add_row(rows, "block_size", case, f"{field}_median", median(item.get(field) for item in group), unit, status, "median across query classes", [block_rel])
        add_row(rows, "block_size", case, "failed_count", sum(raw_status(item) == "failed" for item in group), "cases", "measured", "explicit failed case count", [block_rel])
        add_row(rows, "block_size", case, "full_object_fallback_count", sum(int(item.get("full_object_fallbacks") or 0) for item in group), "fallbacks", "measured", "explicit fallback count", [block_rel])

    memory_rel = "resource-scaling-v1/memory/summary.json"
    memory_obj = load_object(json_values, memory_rel)
    memory_rows = memory_obj.get("rows") if isinstance(memory_obj.get("rows"), list) else []
    complete = [item for item in memory_rows if isinstance(item, dict) and raw_status(item) == "measured"]
    case = "all-288-cases" if len(memory_rows) == 288 else f"all-{len(memory_rows)}-cases"
    add_row(rows, "memory", case, "case_count", len(memory_rows), "cases", "measured" if memory_rows else "unavailable", "raw memory matrix case count", [memory_rel])
    for field, unit in (("wall_seconds", "seconds"), ("peak_rss_kib", "KiB"), ("candidate_count", "candidates"), ("alignment_count", "alignments")):
        add_row(rows, "memory", case, f"{field}_median", median(item.get(field) for item in complete), unit, "measured" if complete else "unavailable", "median over complete memory cases", [memory_rel])
    add_row(rows, "memory", case, "peak_rss_kib_max", maximum(item.get("peak_rss_kib") for item in complete), "KiB", "measured" if complete else "unavailable", "maximum over complete memory cases", [memory_rel])
    add_row(rows, "memory", case, "complete_count", len(complete), "cases", "measured", "explicit status count", [memory_rel])
    add_row(rows, "memory", case, "failed_count", sum(raw_status(item) == "failed" for item in memory_rows), "cases", "measured", "explicit status count", [memory_rel])
    add_row(rows, "memory", case, "unavailable_count", sum(raw_status(item) == "unavailable" for item in memory_rows), "cases", "measured", "explicit status count", [memory_rel])
    for field, unit in (("full_object_fallbacks", "fallbacks"), ("stale_cache_rejections", "rejections"), ("retries", "retries")):
        add_row(rows, "memory", case, f"{field}_total", sum(int(item.get(field) or 0) for item in memory_rows), unit, "measured", "sum of explicit raw counters", [memory_rel])


def summarize_http_and_catalog(
    rows: list[dict[str, Any]],
    json_values: dict[str, JsonValue],
) -> None:
    for cache_state in ("cold", "warm"):
        rel = f"resource-scaling-v1/http/{cache_state}.measurement.json"
        value = load_object(json_values, rel)
        case = f"http/{cache_state}"
        add_measurement_fields(rows, "http", case, value, rel, (("wall_seconds", "seconds"), ("user_seconds", "seconds"), ("system_seconds", "seconds"), ("peak_rss_kib", "KiB")), "HTTP range measurement")
        server = value.get("server_metrics") if isinstance(value.get("server_metrics"), dict) else {}
        resource = value.get("resource_metrics") if isinstance(value.get("resource_metrics"), dict) else {}
        for field, unit in (
            ("responses_206", "responses"),
            ("responses_200", "responses"),
            ("requests", "requests"),
            ("range_requests", "requests"),
            ("head_requests", "requests"),
            ("get_requests", "requests"),
            ("bytes_served", "bytes"),
        ):
            source = server if field in server else resource
            add_row(rows, "http", case, field, source.get(field), unit, raw_status(value), "HTTP server or client metric", [rel])
    metrics_rel = "resource-scaling-v1/http/http.metrics.json"
    metrics = load_object(json_values, metrics_rel)
    for field, unit in (("responses_206", "responses"), ("responses_200", "responses"), ("requests", "requests"), ("bytes_served", "bytes")):
        add_row(rows, "http", "combined-cold-warm", field, metrics.get(field), unit, "measured", "combined HTTP fixture metric", [metrics_rel])

    for thread in (1, 4, 8, 16):
        rel = f"resource-scaling-v1/large-catalog/jam_q50000_sensitive_c100_t{thread}.json"
        value = load_object(json_values, rel)
        case = f"threads={thread}"
        add_measurement_fields(rows, "large_catalog", case, value, rel, (("wall_seconds", "seconds"), ("user_seconds", "seconds"), ("system_seconds", "seconds"), ("peak_rss_kib", "KiB"), ("candidate_count", "candidates"), ("alignment_count", "alignments")), "large-catalog Jam measurement")
        resource = value.get("resource_metrics") if isinstance(value.get("resource_metrics"), dict) else {}
        for field, unit in (("range_requests", "requests"), ("requested_bytes", "bytes"), ("returned_bytes", "bytes"), ("decoded_bytes", "bytes"), ("full_object_fallbacks", "fallbacks")):
            add_row(rows, "large_catalog", case, field, resource.get(field), unit, raw_status(value), "large-catalog resource metric", [rel])


def add_unavailable(rows: list[dict[str, Any]], case: str, note: str) -> None:
    add_row(rows, "availability", case, "status", "NA", "status", "unavailable", note)


def format_value(value: Any) -> str:
    if value is None or value == "NA":
        return "NA"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return format(value, ".12g")
    return str(value).replace("\t", " ").replace("\n", " ")


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_value(row.get(field)) for field in TSV_FIELDS})


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    raw_root = args.raw_root.resolve()
    output_dir = args.output_dir.resolve()
    if not raw_root.is_dir():
        raise SystemExit(f"raw root is not a directory: {raw_root}")
    if output_dir.exists():
        raise SystemExit(f"refusing existing output directory: {output_dir}")
    if output_dir == raw_root or raw_root in output_dir.parents:
        raise SystemExit("output directory must not be inside raw root")

    files = [path for path in sorted(raw_root.rglob("*")) if path.is_file()]
    checksums = {relative(path, raw_root): sha256(path) for path in files}
    json_values, tsv_values, jsonl_values, errors = validate_inputs(raw_root)
    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        raise SystemExit(f"raw input validation failed with {len(errors)} error(s)")

    rows: list[dict[str, Any]] = []
    summarize_controlled(rows, json_values, tsv_values)
    summarize_accession(rows, json_values, jsonl_values, tsv_values)
    summarize_block_and_memory(rows, json_values)
    summarize_http_and_catalog(rows, json_values)
    add_unavailable(rows, "actual_s3", "No actual S3 object-storage run is present; the HTTP range fixture is not an S3 run.")
    add_unavailable(rows, "read_derived_assemblies", "No read-derived assembly track is present in the raw tree.")
    add_unavailable(rows, "natural_supported_positives", "No natural or independently supported positive track is present in the raw tree.")
    add_unavailable(rows, "1000_real_assemblies", "The 1,000-real-assembly release scale was not run.")
    add_unavailable(rows, "100_plasmid_queries", "The 100-query release scale was not run.")
    rows.sort(key=lambda row: (row["section"], row["case"], row["metric"], row["unit"], row["status"]))

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir()
    checksum_value = {
        "schema_version": SCHEMA_VERSION,
        "raw_root": str(raw_root),
        "file_count": len(checksums),
        "files": checksums,
    }
    summary_value = {
        "schema_version": SCHEMA_VERSION,
        "raw_root": str(raw_root),
        "raw_file_count": len(checksums),
        "raw_checksums_file": "raw-checksums.json",
        "validation": {
            "status": "valid",
            "json_file_count": len(json_values),
            "tsv_file_count": len(tsv_values),
            "jsonl_file_count": len(jsonl_values),
            "errors": [],
        },
        "row_count": len(rows),
        "rows": rows,
        "availability_note": "unavailable means no source artifact was present; it is not a negative biological result",
    }
    (output_dir / "summary.json").write_text(json.dumps(summary_value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_tsv(output_dir / "summary.tsv", rows)
    (output_dir / "raw-checksums.json").write_text(json.dumps(checksum_value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
