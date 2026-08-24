"""Load and attribute ``TraceDiagnosticReport`` JSON values.

The Rust models are serialized with snake-case field names.  This module keeps
the parser tolerant of a JSONL run containing a nested ``diagnostics`` member
and of the older rescue-round counter names, while the emitted analysis stays
on one explicit schema.  No base is silently counted twice: coordinate-aware
reports use interval union, and count-only reports reject overlapping truth
intervals.
"""

from __future__ import annotations

import json
import subprocess
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

from .intervals import (
    Interval,
    IntervalError,
    interval_length,
    interval_metrics,
    parse_interval,
    parse_intervals,
    subtract_intervals,
    union_intervals,
)

FAILURE_CATEGORIES = (
    "candidate_miss",
    "no_retained_seed",
    "no_matching_seed",
    "repetitive_seed_suppression",
    "anchor_cap",
    "no_valid_chain",
    "chain_limit",
    "sequence_budget",
    "alignment_band_failure",
    "alignment_rejection",
    "alignment_cap",
    "mosaic_selection",
    "other",
)

STAGE_FIELDS = (
    "stage",
    "name",
    "seed_pages_read",
    "keys_tested",
    "occurrences_decoded",
    "anchors_retained",
    "chains_retained",
    "sequence_blocks_read",
    "alignment_attempts",
    "new_query_bases_supported",
    "wall_micros",
    "cpu_micros",
    "bytes_read",
)

_STAGE_ALIASES: dict[str, tuple[str, ...]] = {
    "seed_pages_read": ("seed_pages_read", "seed_buckets_requested", "seed_pages_requested"),
    "keys_tested": ("keys_tested", "seed_keys_tested"),
    "occurrences_decoded": ("occurrences_decoded", "occurrences_after_limits", "occurrences_decoded"),
    "anchors_retained": ("anchors_retained", "anchors_created", "anchors_after_limits"),
    "chains_retained": ("chains_retained", "chains_accepted"),
    "sequence_blocks_read": ("sequence_blocks_read", "sequence_blocks_fetched"),
    "alignment_attempts": ("alignment_attempts", "alignment_windows_attempted"),
    "new_query_bases_supported": ("new_query_bases_supported",),
    "wall_micros": ("wall_micros",),
    "cpu_micros": ("cpu_micros",),
    "bytes_read": ("bytes_read",),
}

_SEED_FIELDS = (
    "scheme_id",
    "query_seeds_emitted",
    "matching_target_keys",
    "matching_occurrences_before_limits",
    "matching_occurrences_after_limits",
)

_SEED_ALIASES: dict[str, tuple[str, ...]] = {
    "query_seeds_emitted": ("query_seeds_emitted", "seeds_emitted", "query_seed_count"),
    "matching_target_keys": ("matching_target_keys", "target_keys_matched", "matching_keys"),
    "matching_occurrences_before_limits": (
        "matching_occurrences_before_limits",
        "occurrences_before_limits",
        "matching_occurrences",
    ),
    "matching_occurrences_after_limits": (
        "matching_occurrences_after_limits",
        "occurrences_after_limits",
        "retained_occurrences",
    ),
}


class DiagnosticError(ValueError):
    """Raised when a diagnostic report cannot support exact accounting."""


def _integer(value: Any, field: str, *, default: int = 0) -> int:
    if value is None:
        return default
    if isinstance(value, bool):
        raise DiagnosticError(f"{field} must be an integer")
    try:
        number = int(value)
    except (TypeError, ValueError) as exc:
        raise DiagnosticError(f"{field} must be an integer") from exc
    if isinstance(value, float) and value != number:
        raise DiagnosticError(f"{field} must be an integer")
    if number < 0:
        raise DiagnosticError(f"{field} cannot be negative")
    return number


def _bool(value: Any, *, default: bool = False) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and value in {0, 1}:
        return bool(value)
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1", "yes", "selected", "accepted"}:
            return True
        if normalized in {"false", "0", "no", "rejected", "not_selected"}:
            return False
    raise DiagnosticError(f"boolean value is malformed: {value!r}")


def _first(value: dict[str, Any], names: Iterable[str], default: Any = None) -> Any:
    for name in names:
        if name in value:
            return value[name]
    return default


def _category(value: Any) -> str | None:
    if value is None or value == "":
        return None
    normalized = str(value).strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "candidate_miss": "candidate_miss",
        "candidate_missing": "candidate_miss",
        "no_seed": "no_retained_seed",
        "no_retained_seed": "no_retained_seed",
        "no_matching_seed": "no_matching_seed",
        "repetitive_seed_suppression": "repetitive_seed_suppression",
        "repetitive_suppression": "repetitive_seed_suppression",
        "anchor_cap": "anchor_cap",
        "no_valid_chain": "no_valid_chain",
        "chain_limit": "chain_limit",
        "sequence_budget": "sequence_budget",
        "alignment_band_failure": "alignment_band_failure",
        "band_failure": "alignment_band_failure",
        "alignment_rejection": "alignment_rejection",
        "alignment_cap": "alignment_cap",
        "mosaic_selection": "mosaic_selection",
        "other": "other",
    }
    return aliases.get(normalized, "other")


def _interval_or_none(value: Any, *, field: str) -> dict[str, int] | None:
    if value is None:
        return None
    try:
        start, end = parse_interval(value, field=field)
    except IntervalError as exc:
        raise DiagnosticError(str(exc)) from exc
    return {"start": start, "end": end}


def _normal_seed(raw: Any, index: int) -> dict[str, int]:
    if not isinstance(raw, dict):
        raise DiagnosticError(f"seed_schemes[{index}] must be an object")
    scheme_id = _first(raw, ("scheme_id", "id"), 0)
    result: dict[str, int] = {"scheme_id": _integer(scheme_id, f"seed_schemes[{index}].scheme_id")}
    for field in _SEED_FIELDS[1:]:
        result[field] = _integer(
            _first(raw, _SEED_ALIASES[field]),
            f"seed_schemes[{index}].{field}",
        )
    return result


def _normal_stage(raw: Any, index: int) -> dict[str, Any]:
    if not isinstance(raw, dict):
        raise DiagnosticError(f"stages[{index}] must be an object")
    stage = _integer(_first(raw, ("stage", "round"), index), f"stages[{index}].stage")
    name = str(_first(raw, ("name", "stage_name"), f"stage-{stage}"))
    result: dict[str, Any] = {"stage": stage, "name": name}
    for field in STAGE_FIELDS[2:]:
        result[field] = _integer(
            _first(raw, _STAGE_ALIASES[field]),
            f"stages[{index}].{field}",
        )
    return result


def _truth_interval(raw: Any, index: int, *, candidate_default: bool = False) -> dict[str, Any]:
    if not isinstance(raw, dict):
        raise DiagnosticError(f"truth_intervals[{index}] must be an object")
    interval_value = _first(raw, ("truth_interval", "interval"), raw)
    try:
        start, end = parse_interval(interval_value, field=f"truth_intervals[{index}].truth_interval")
    except IntervalError as exc:
        raise DiagnosticError(str(exc)) from exc
    if start == end:
        raise DiagnosticError(f"truth_intervals[{index}] must be non-empty")
    result: dict[str, Any] = {
        "truth_interval": {"start": start, "end": end},
        "candidate_selected": _bool(
            _first(raw, ("candidate_selected", "selected")),
            default=candidate_default,
        ),
        "seed_schemes": [
            _normal_seed(seed, seed_index)
            for seed_index, seed in enumerate(
                _first(raw, ("seed_schemes", "schemes"), []) or []
            )
        ],
        "anchors_before_limits": _integer(
            _first(raw, ("anchors_before_limits", "anchors_before")),
            f"truth_intervals[{index}].anchors_before_limits",
        ),
        "anchors_after_limits": _integer(
            _first(raw, ("anchors_after_limits", "anchors_after")),
            f"truth_intervals[{index}].anchors_after_limits",
        ),
        "best_chain_anchor_count": _integer(
            _first(raw, ("best_chain_anchor_count", "chain_anchor_count")),
            f"truth_intervals[{index}].best_chain_anchor_count",
        ),
        "best_chain_query_span": _integer(
            _first(raw, ("best_chain_query_span", "chain_query_span")),
            f"truth_intervals[{index}].best_chain_query_span",
        ),
        "best_chain_target_span": _integer(
            _first(raw, ("best_chain_target_span", "chain_target_span")),
            f"truth_intervals[{index}].best_chain_target_span",
        ),
        "sequence_blocks_requested": _integer(
            _first(raw, ("sequence_blocks_requested", "sequence_blocks")),
            f"truth_intervals[{index}].sequence_blocks_requested",
        ),
        "sequence_window_decoded": _interval_or_none(
            _first(raw, ("sequence_window_decoded", "decoded_window")),
            field=f"truth_intervals[{index}].sequence_window_decoded",
        ),
        "alignment_attempted": _bool(_first(raw, ("alignment_attempted", "alignment"))),
        "alignment_band_widths": [
            _integer(width, f"truth_intervals[{index}].alignment_band_widths[{width_index}]")
            for width_index, width in enumerate(
                _first(raw, ("alignment_band_widths", "band_widths"), []) or []
            )
        ],
        "alignment_band_edge_touched": _bool(
            _first(raw, ("alignment_band_edge_touched", "band_edge_touched"))
        ),
        "alignment_accepted": _bool(_first(raw, ("alignment_accepted", "accepted"))),
        "alignment_query_span": _integer(
            _first(raw, ("alignment_query_span", "aligned_query_span")),
            f"truth_intervals[{index}].alignment_query_span",
        ),
        "mosaic_supported_bases": _integer(
            _first(raw, ("mosaic_supported_bases", "supported_bases")),
            f"truth_intervals[{index}].mosaic_supported_bases",
        ),
        "alignment_or_output_cap_reached": _bool(
            _first(raw, ("alignment_or_output_cap_reached", "output_cap_reached"))
        ),
        "failure_category": _category(_first(raw, ("failure_category", "failure"))),
    }
    return result


def _report_supported_intervals(raw: dict[str, Any]) -> list[Interval]:
    candidates = [
        raw.get("supported_intervals"),
        raw.get("mosaic_supported_intervals"),
    ]
    coverage = raw.get("coverage")
    if isinstance(coverage, dict):
        candidates.extend((coverage.get("supported_intervals"), coverage.get("covered_intervals")))
    mosaic = raw.get("mosaic")
    if isinstance(mosaic, dict):
        candidates.append(mosaic.get("covered_intervals"))
    for value in candidates:
        if value:
            try:
                return parse_intervals(value, field="supported_intervals")
            except IntervalError as exc:
                raise DiagnosticError(str(exc)) from exc
    return []


def normalize_report(raw: dict[str, Any]) -> dict[str, Any]:
    """Normalize one serialized ``TraceDiagnosticReport``."""

    if not isinstance(raw, dict):
        raise DiagnosticError("diagnostic report must be an object")
    nested = raw.get("diagnostics", raw.get("trace_diagnostics"))
    if isinstance(nested, dict):
        raw = nested
    truth_raw = raw.get("truth_intervals", raw.get("truth", [])) or []
    if not isinstance(truth_raw, list):
        raise DiagnosticError("truth_intervals must be a list")
    candidate_default = _bool(raw.get("candidate_selected"), default=False)
    truth = [_truth_interval(item, index, candidate_default=candidate_default) for index, item in enumerate(truth_raw)]
    stages_raw = raw.get("stages", raw.get("stage_metrics", [])) or []
    if not isinstance(stages_raw, list):
        raise DiagnosticError("stages must be a list")
    stages = [_normal_stage(item, index) for index, item in enumerate(stages_raw)]
    metagenome_id = str(raw.get("metagenome_id", raw.get("assembly", "unknown")))
    return {
        "schema_version": str(raw.get("schema_version", "trace-diagnostics-v1")),
        "metagenome_id": metagenome_id,
        "truth_intervals": truth,
        "stages": stages,
        "supported_intervals": _report_supported_intervals(raw),
        "source": raw,
    }


def _iter_report_values(value: Any) -> Iterable[dict[str, Any]]:
    if isinstance(value, list):
        for item in value:
            yield from _iter_report_values(item)
        return
    if not isinstance(value, dict):
        return
    for key in ("diagnostics", "trace_diagnostics", "diagnostic_report"):
        nested = value.get(key)
        if isinstance(nested, (dict, list)):
            yield from _iter_report_values(nested)
            return
    records = value.get("records")
    if isinstance(records, list) and "truth_intervals" not in value:
        yielded = False
        for item in records:
            before = yielded
            found = list(_iter_report_values(item))
            if found:
                yielded = True
                yield from found
            elif isinstance(item, dict) and "diagnostics" in item:
                yielded = True
            if before and not found:
                continue
        if yielded:
            return
    if "truth_intervals" in value or "stage_metrics" in value or "stages" in value:
        yield value


def _read_text(path: Path) -> str:
    if path.suffix.lower() not in {".zst", ".zstd"}:
        return path.read_text(encoding="utf-8")
    try:
        completed = subprocess.run(
            ["zstd", "--decompress", "--stdout", str(path)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        raise DiagnosticError(f"cannot decompress {path}: {exc}") from exc
    try:
        return completed.stdout.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise DiagnosticError(f"{path} is not UTF-8 JSON") from exc


def load_reports(path: str | Path) -> list[dict[str, Any]]:
    """Load direct reports, JSON arrays, or JSONL trace output."""

    source = Path(path)
    text = _read_text(source)
    try:
        value = json.loads(text)
    except json.JSONDecodeError:
        values: list[Any] = []
        for line_number, line in enumerate(text.splitlines(), 1):
            if not line.strip():
                continue
            try:
                values.append(json.loads(line))
            except json.JSONDecodeError as exc:
                raise DiagnosticError(f"invalid JSON at {source}:{line_number}") from exc
        value = values
    reports = [normalize_report(item) for item in _iter_report_values(value)]
    if not reports:
        raise DiagnosticError(f"{source} contains no TraceDiagnosticReport values")
    return reports


def _seed_totals(truth: dict[str, Any]) -> dict[str, int]:
    totals = {
        "query_seeds_emitted": 0,
        "matching_target_keys": 0,
        "matching_occurrences_before_limits": 0,
        "matching_occurrences_after_limits": 0,
    }
    for scheme in truth["seed_schemes"]:
        for field in totals:
            totals[field] += scheme[field]
    return totals


def _choose_category(truth: dict[str, Any], report: dict[str, Any]) -> str:
    explicit = truth.get("failure_category")
    if explicit in FAILURE_CATEGORIES:
        return explicit
    if not truth["candidate_selected"]:
        return "candidate_miss"
    seeds = _seed_totals(truth)
    if seeds["query_seeds_emitted"] == 0:
        return "no_retained_seed"
    if seeds["matching_target_keys"] == 0:
        return "no_matching_seed"
    if (
        seeds["matching_occurrences_before_limits"] > 0
        and seeds["matching_occurrences_after_limits"] == 0
    ):
        return "repetitive_seed_suppression"
    if truth["anchors_before_limits"] > truth["anchors_after_limits"]:
        return "anchor_cap"
    if truth["anchors_after_limits"] > 0 and truth["best_chain_anchor_count"] == 0:
        return "no_valid_chain"
    if truth["best_chain_anchor_count"] > 0 and truth["sequence_blocks_requested"] == 0:
        return "sequence_budget"
    if truth["alignment_or_output_cap_reached"]:
        return "alignment_cap"
    if truth["alignment_band_edge_touched"] and not truth["alignment_accepted"]:
        return "alignment_band_failure"
    if truth["alignment_attempted"] and not truth["alignment_accepted"]:
        return "alignment_rejection"
    if truth["alignment_accepted"] and truth["mosaic_supported_bases"] < interval_length(
        [(truth["truth_interval"]["start"], truth["truth_interval"]["end"])]
    ):
        return "mosaic_selection"
    return "other"


def _truth_interval_tuple(truth: dict[str, Any]) -> Interval:
    interval = truth["truth_interval"]
    return interval["start"], interval["end"]


def _split_at_boundaries(
    intervals: Iterable[Interval], boundaries: Iterable[int]
) -> list[Interval]:
    """Split intervals at every truth boundary without merging the pieces."""

    boundary_set = set(boundaries)
    split: list[Interval] = []
    for start, end in intervals:
        cuts = sorted({start, end} | {boundary for boundary in boundary_set if start < boundary < end})
        split.extend((left, right) for left, right in zip(cuts, cuts[1:]) if left < right)
    return split


def attribute_report(report: dict[str, Any]) -> dict[str, Any]:
    """Assign every missing truth base in one report to one failure category."""

    normalized = (
        report
        if {"truth_intervals", "stages", "supported_intervals", "source"}.issubset(report)
        else normalize_report(report)
    )
    truth = normalized["truth_intervals"]
    truth_intervals = [_truth_interval_tuple(item) for item in truth]
    supported_intervals = normalized["supported_intervals"]
    coordinate_mode = bool(supported_intervals)
    if coordinate_mode:
        truth_union = union_intervals(truth_intervals)
        supported_union = union_intervals(supported_intervals)
        truth_boundaries = [coordinate for interval in truth_intervals for coordinate in interval]
        missing_intervals = _split_at_boundaries(
            subtract_intervals(truth_union, supported_union), truth_boundaries
        )
        metrics = interval_metrics(supported_union, truth_union)
        attributions: list[dict[str, Any]] = []
        for start, end in missing_intervals:
            owners = [
                item
                for item in truth
                if _truth_interval_tuple(item)[0] < end and start < _truth_interval_tuple(item)[1]
            ]
            categories = {_choose_category(item, normalized) for item in owners}
            category = next(iter(categories)) if len(categories) == 1 else "other"
            attributions.append(
                {
                    "interval": {"start": start, "end": end},
                    "missing_bases": end - start,
                    "category": category,
                }
            )
        truth_bases = metrics["truth_bases"]
        supported_bases = metrics["intersection_bases"]
        accounting_mode = "coordinate"
    else:
        # Count-only reports cannot prove a unique answer for overlapping truth
        # intervals.  Fail closed instead of inflating the missing-base count.
        if interval_length(truth_intervals) != sum(end - start for start, end in truth_intervals):
            raise DiagnosticError(
                f"{normalized['metagenome_id']}: overlapping truth intervals require supported_intervals"
            )
        attributions = []
        truth_bases = 0
        supported_bases = 0
        for item in truth:
            start, end = _truth_interval_tuple(item)
            length = end - start
            supported = item["mosaic_supported_bases"]
            if supported > length:
                raise DiagnosticError(
                    f"{normalized['metagenome_id']}: mosaic_supported_bases exceeds truth interval"
                )
            truth_bases += length
            supported_bases += supported
            missing = length - supported
            if missing:
                attributions.append(
                    {
                        "interval": {"start": start, "end": end},
                        "missing_bases": missing,
                        "category": _choose_category(item, normalized),
                    }
                )
        metrics = {
            "observed_bases": supported_bases,
            "truth_bases": truth_bases,
            "intersection_bases": supported_bases,
            "base_precision": None,
            "base_recall": supported_bases / truth_bases if truth_bases else None,
            "observed_intervals": sum(1 for item in truth if item["mosaic_supported_bases"]),
            "truth_intervals": len(truth),
            "interval_precision": None,
            "interval_recall": (
                sum(1 for item in truth if item["mosaic_supported_bases"]) / len(truth)
                if truth
                else None
            ),
        }
        accounting_mode = "count_only"
    missing_total = sum(item["missing_bases"] for item in attributions)
    expected_missing = truth_bases - supported_bases
    if missing_total != expected_missing:
        raise DiagnosticError(
            f"{normalized['metagenome_id']}: attribution totals {missing_total}, expected {expected_missing}"
        )
    category_counts = Counter()
    for item in attributions:
        category_counts[item["category"]] += item["missing_bases"]
    if set(category_counts) - set(FAILURE_CATEGORIES):
        raise DiagnosticError("attribution contains an unknown failure category")
    stage_totals = {
        field: sum(stage[field] for stage in normalized["stages"])
        for field in STAGE_FIELDS[2:]
    }
    return {
        "schema_version": "trace-failure-attribution-v1",
        "metagenome_id": normalized["metagenome_id"],
        "accounting_mode": accounting_mode,
        "truth_bases": truth_bases,
        "supported_bases": supported_bases,
        "missing_truth_bases": missing_total,
        "metrics": metrics,
        "attributions": attributions,
        "failure_counts": {
            category: category_counts.get(category, 0) for category in FAILURE_CATEGORIES
        },
        "stages": normalized["stages"],
        "stage_totals": stage_totals,
    }


def summarize_reports(reports: Iterable[dict[str, Any]]) -> dict[str, Any]:
    """Aggregate attributions and every stage counter across reports."""

    attributed = [attribute_report(report) for report in reports]
    truth_bases = sum(item["truth_bases"] for item in attributed)
    supported_bases = sum(item["supported_bases"] for item in attributed)
    observed_bases = sum(int(item["metrics"]["observed_bases"]) for item in attributed)
    intersection_bases = sum(int(item["metrics"]["intersection_bases"]) for item in attributed)
    missing = sum(item["missing_truth_bases"] for item in attributed)
    failures = Counter()
    for item in attributed:
        failures.update(item["failure_counts"])
    stage_totals = {
        field: sum(item["stage_totals"][field] for item in attributed)
        for field in STAGE_FIELDS[2:]
    }
    precision_available = bool(attributed) and all(
        item["accounting_mode"] == "coordinate" for item in attributed
    )
    base_precision = (
        intersection_bases / observed_bases if precision_available and observed_bases else None
    )
    observed_intervals = sum(int(item["metrics"]["observed_intervals"]) for item in attributed)
    interval_precision = None
    if precision_available and observed_intervals:
        interval_precision = sum(
            float(item["metrics"]["interval_precision"] or 0.0)
            * int(item["metrics"]["observed_intervals"])
            for item in attributed
        ) / observed_intervals
    return {
        "schema_version": "trace-failure-attribution-summary-v1",
        "report_count": len(attributed),
        "truth_bases": truth_bases,
        "supported_bases": supported_bases,
        "observed_bases": observed_bases,
        "intersection_bases": intersection_bases,
        "missing_truth_bases": missing,
        "base_precision": base_precision,
        "base_recall": intersection_bases / truth_bases if truth_bases else None,
        "interval_precision": interval_precision,
        "failure_counts": {
            category: failures.get(category, 0) for category in FAILURE_CATEGORIES
        },
        "stage_totals": stage_totals,
        "reports": attributed,
    }
