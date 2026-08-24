"""Normalize comparator alignments onto the shared query coordinate system."""

from __future__ import annotations

from typing import Any, Iterable

from .intervals import Interval, IntervalError, interval_length, interval_metrics, project_alignment, union_intervals


class ComparatorError(ValueError):
    """Raised when comparator evidence cannot be projected safely."""


def _records(value: Any) -> list[dict[str, Any]]:
    if isinstance(value, list):
        records: list[dict[str, Any]] = []
        for item in value:
            records.extend(_records(item))
        return records
    if not isinstance(value, dict):
        return []
    if isinstance(value.get("alignments"), list):
        return [item for item in value["alignments"] if isinstance(item, dict)]
    if isinstance(value.get("records"), list):
        records = []
        for item in value["records"]:
            records.extend(_records(item))
        return records
    if any(key in value for key in ("query_segments", "segments", "query_start", "query_end", "cigar")):
        return [value]
    return []


def _query_segments(record: dict[str, Any]) -> Any:
    if "query_segments" in record:
        return record["query_segments"]
    if "segments" in record:
        return record["segments"]
    if "interval" in record:
        return record["interval"]
    if "query_start" in record and "query_end" in record:
        return {"start": record["query_start"], "end": record["query_end"]}
    raise ComparatorError("alignment has no query coordinates")


def normalize_comparator(
    value: Any,
    query_length: int,
    *,
    truth_intervals: Iterable[Interval] = (),
    query_deletion_operation: str = "I",
) -> dict[str, Any]:
    """Return normalized supported intervals and shared overlap metrics."""

    normalized_alignments: list[dict[str, Any]] = []
    supported: list[Interval] = []
    aligned: list[Interval] = []
    deletions: list[Interval] = []
    for index, record in enumerate(_records(value)):
        try:
            projected = project_alignment(
                _query_segments(record),
                query_length,
                cigar=record.get("cigar"),
                edit_script=record.get("edit_script"),
                query_deletion_operation=query_deletion_operation,
            )
        except IntervalError as exc:
            raise ComparatorError(f"alignment[{index}] cannot be normalized: {exc}") from exc
        supported.extend(projected["supported_intervals"])
        aligned.extend(projected["aligned_intervals"])
        deletions.extend(projected["alignment_deletions"])
        normalized_alignments.append(
            {
                "alignment_index": index,
                "supported_intervals": [
                    {"start": start, "end": end}
                    for start, end in projected["supported_intervals"]
                ],
                "aligned_intervals": [
                    {"start": start, "end": end}
                    for start, end in projected["aligned_intervals"]
                ],
                "alignment_deletions": [
                    {"start": start, "end": end}
                    for start, end in projected["alignment_deletions"]
                ],
            }
        )
    supported_union = union_intervals(supported)
    aligned_union = union_intervals(aligned)
    deletion_union = union_intervals(deletions)
    truth_union = union_intervals(truth_intervals)
    metrics = interval_metrics(supported_union, truth_union) if truth_union else {
        "observed_bases": interval_length(supported_union),
        "truth_bases": 0,
        "intersection_bases": 0,
        "gap_error_bases": interval_length(supported_union),
        "base_precision": None,
        "base_recall": None,
        "observed_intervals": len(supported_union),
        "truth_intervals": 0,
        "interval_precision": None,
        "interval_recall": None,
    }
    return {
        "schema_version": "normalized-comparator-v1",
        "query_length": query_length,
        "alignments": normalized_alignments,
        "supported_intervals": [
            {"start": start, "end": end} for start, end in supported_union
        ],
        "aligned_intervals": [{"start": start, "end": end} for start, end in aligned_union],
        "alignment_deletions": [
            {"start": start, "end": end} for start, end in deletion_union
        ],
        "metrics": metrics,
    }
