"""Assertions for the committed normalized 40-case baseline."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any


class BaselineError(ValueError):
    """Raised when a normalized baseline does not match the fixture."""


def _number(value: Any, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise BaselineError(f"{field} is not numeric")
    if not math.isfinite(float(value)):
        raise BaselineError(f"{field} is not finite")
    return float(value)


def _metrics(value: dict[str, Any]) -> dict[str, Any]:
    if isinstance(value.get("metrics"), dict):
        return value["metrics"]
    if isinstance(value.get("coverage"), dict) and isinstance(value["coverage"].get("metrics"), dict):
        return value["coverage"]["metrics"]
    raise BaselineError("normalized artifact has no metrics object")


def assert_baseline(
    normalized_path: str | Path,
    fixture_path: str | Path,
) -> dict[str, Any]:
    """Validate aggregate counts and precision/recall against the fixture.

    The normalized production artifact stores one result per assembly under
    ``results``.  The helper also accepts a single normalized record so it can
    be used by focused fixture tests.
    """

    normalized_file = Path(normalized_path)
    fixture_file = Path(fixture_path)
    try:
        normalized = json.loads(normalized_file.read_text(encoding="utf-8"))
        fixture = json.loads(fixture_file.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise BaselineError(f"cannot read baseline input: {exc}") from exc
    if not isinstance(normalized, dict) or not isinstance(fixture, dict):
        raise BaselineError("baseline inputs must be JSON objects")
    metrics = _metrics(normalized)
    results = normalized.get("results")
    if isinstance(results, list) and results:
        case_count = len(results)
        selected_candidates = sum(
            1
            for result in results
            if isinstance(result, dict)
            and result.get("status") not in {"unavailable", "failed"}
            and result.get("jam_status") not in {"no_candidate", "failed"}
        )
    else:
        case_count = int(metrics.get("assemblies_total", 1))
        selected_candidates = int(metrics.get("assemblies_with_records", case_count))
    actual: dict[str, Any] = {
        "case_count": case_count,
        "selected_candidates": selected_candidates,
        "assemblies_with_records": int(metrics.get("assemblies_with_records", selected_candidates)),
        "truth_bases": int(metrics.get("truth_bases", 0)),
        "observed_bases": int(metrics.get("observed_bases", 0)),
        "intersection_bases": int(metrics.get("intersection_bases", 0)),
        "missing_truth_bases": int(metrics.get("gap_error_bases", metrics.get("missing_truth_bases", 0))),
        "base_precision": _number(metrics.get("base_precision"), "base_precision"),
        "base_recall": _number(metrics.get("base_recall"), "base_recall"),
        "interval_precision": _number(metrics.get("interval_precision"), "interval_precision"),
        "interval_recall": _number(metrics.get("interval_recall"), "interval_recall"),
    }
    tolerance = float(fixture.get("tolerance", 1e-12))
    integer_fields = (
        "case_count",
        "selected_candidates",
        "assemblies_with_records",
        "truth_bases",
        "observed_bases",
        "intersection_bases",
        "missing_truth_bases",
    )
    float_fields = (
        "base_precision",
        "base_recall",
        "interval_precision",
        "interval_recall",
    )
    for field in integer_fields:
        if int(fixture.get(field, -1)) != actual[field]:
            raise BaselineError(f"{field}: expected {fixture.get(field)!r}, got {actual[field]!r}")
    for field in float_fields:
        expected = _number(fixture.get(field), field)
        if abs(expected - actual[field]) > tolerance:
            raise BaselineError(f"{field}: expected {expected!r}, got {actual[field]!r}")
    return {"status": "pass", "actual": actual, "fixture": fixture_file.as_posix()}
