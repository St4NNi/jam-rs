"""Small, shared interval utilities for trace and comparator experiments.

All intervals are zero-based and half-open.  Empty intervals are ignored by
union and subtraction, while reversed or non-integral coordinates are errors.
The projection helper accepts both Jam edit-script names and SAM/PAF CIGAR
operations so every benchmark path can use one coordinate normalization.
"""

from __future__ import annotations

import re
from typing import Any, Iterable

Interval = tuple[int, int]
_CIGAR_RE = re.compile(r"([0-9]+)([MIDNSHP=X])")


class IntervalError(ValueError):
    """Raised when an interval or alignment projection is malformed."""


def _integer(value: Any, field: str) -> int:
    if isinstance(value, bool):
        raise IntervalError(f"{field} must be an integer")
    try:
        converted = int(value)
    except (TypeError, ValueError) as exc:
        raise IntervalError(f"{field} must be an integer") from exc
    if isinstance(value, float) and value != converted:
        raise IntervalError(f"{field} must be an integer")
    return converted


def parse_interval(value: Any, *, field: str = "interval") -> Interval:
    """Parse one ``{start, end}`` or ``[start, end]`` interval."""

    if isinstance(value, dict):
        if "interval" in value and not {"start", "end"}.issubset(value):
            return parse_interval(value["interval"], field=field)
        if "start" not in value or "end" not in value:
            raise IntervalError(f"{field} requires start and end")
        start = _integer(value["start"], f"{field}.start")
        end = _integer(value["end"], f"{field}.end")
    elif isinstance(value, (list, tuple)) and len(value) == 2:
        start = _integer(value[0], f"{field}.start")
        end = _integer(value[1], f"{field}.end")
    else:
        raise IntervalError(f"{field} must be an object or two-item array")
    if start < 0 or end < 0 or start > end:
        raise IntervalError(f"{field} must satisfy 0 <= start <= end")
    return start, end


def parse_intervals(value: Any, *, field: str = "intervals") -> list[Interval]:
    """Parse a list of intervals while retaining deterministic input order."""

    if value is None or value == "":
        return []
    if isinstance(value, str):
        import json

        try:
            value = json.loads(value)
        except json.JSONDecodeError as exc:
            raise IntervalError(f"{field} is not valid JSON") from exc
    if isinstance(value, dict):
        value = value.get("intervals", [value])
    if not isinstance(value, (list, tuple)):
        raise IntervalError(f"{field} must be a list")
    return [parse_interval(item, field=f"{field}[{index}]") for index, item in enumerate(value)]


def interval_length(intervals: Iterable[Interval]) -> int:
    """Return the length of the union of ``intervals``."""

    return sum(end - start for start, end in union_intervals(intervals))


def union_intervals(intervals: Iterable[Interval]) -> list[Interval]:
    """Merge sorted, overlapping, or directly adjacent intervals."""

    ordered = sorted((start, end) for start, end in intervals if start < end)
    merged: list[Interval] = []
    for start, end in ordered:
        if start < 0 or end < 0 or start > end:
            raise IntervalError("interval must satisfy 0 <= start <= end")
        if not merged or start > merged[-1][1]:
            merged.append((start, end))
        else:
            merged[-1] = merged[-1][0], max(merged[-1][1], end)
    return merged


def subtract_intervals(base: Iterable[Interval], covered: Iterable[Interval]) -> list[Interval]:
    """Return ``base - covered`` as disjoint, sorted intervals."""

    remaining: list[Interval] = []
    covered_union = union_intervals(covered)
    for start, end in union_intervals(base):
        cursor = start
        for covered_start, covered_end in covered_union:
            if covered_end <= cursor:
                continue
            if covered_start >= end:
                break
            if covered_start > cursor:
                remaining.append((cursor, min(covered_start, end)))
            cursor = max(cursor, covered_end)
            if cursor >= end:
                break
        if cursor < end:
            remaining.append((cursor, end))
    return remaining


def intersection_intervals(left: Iterable[Interval], right: Iterable[Interval]) -> list[Interval]:
    """Return the disjoint intersection of two interval collections."""

    intersections: list[Interval] = []
    for left_start, left_end in union_intervals(left):
        for right_start, right_end in union_intervals(right):
            start = max(left_start, right_start)
            end = min(left_end, right_end)
            if start < end:
                intersections.append((start, end))
    return union_intervals(intersections)


def intersection_length(left: Iterable[Interval], right: Iterable[Interval]) -> int:
    return interval_length(intersection_intervals(left, right))


def interval_metrics(observed: Iterable[Interval], truth: Iterable[Interval]) -> dict[str, int | float | None]:
    """Compute nonredundant base and interval overlap metrics.

    ``gap_error_bases`` is the symmetric difference: unsupported truth bases
    plus observed bases outside truth support.
    """

    observed_union = union_intervals(observed)
    truth_union = union_intervals(truth)
    observed_bases = interval_length(observed_union)
    truth_bases = interval_length(truth_union)
    intersection_bases = intersection_length(observed_union, truth_union)
    base_precision = intersection_bases / observed_bases if observed_bases else None
    base_recall = intersection_bases / truth_bases if truth_bases else None
    # An interval is recovered when it has any observed overlap.  This is the
    # same conservative interval-level interpretation used by the production
    # normalized MGE summaries.
    recovered = sum(
        1
        for truth_start, truth_end in truth_union
        if intersection_length([(truth_start, truth_end)], observed_union) > 0
    )
    observed_recovered = sum(
        1
        for observed_start, observed_end in observed_union
        if intersection_length([(observed_start, observed_end)], truth_union) > 0
    )
    interval_recall = recovered / len(truth_union) if truth_union else None
    interval_precision = observed_recovered / len(observed_union) if observed_union else None
    return {
        "observed_bases": observed_bases,
        "truth_bases": truth_bases,
        "intersection_bases": intersection_bases,
        "gap_error_bases": observed_bases + truth_bases - 2 * intersection_bases,
        "base_precision": base_precision,
        "base_recall": base_recall,
        "observed_intervals": len(observed_union),
        "truth_intervals": len(truth_union),
        "interval_precision": interval_precision,
        "interval_recall": interval_recall,
    }


def _cigar_runs(cigar: str) -> list[tuple[str, int]]:
    if not cigar or cigar == "*":
        return []
    runs: list[tuple[str, int]] = []
    offset = 0
    for match in _CIGAR_RE.finditer(cigar):
        if match.start() != offset:
            raise IntervalError(f"invalid CIGAR at byte {offset}: {cigar!r}")
        length = int(match.group(1))
        if length <= 0:
            raise IntervalError("CIGAR contains a zero-length operation")
        runs.append((match.group(2), length))
        offset = match.end()
    if offset != len(cigar):
        raise IntervalError(f"invalid CIGAR at byte {offset}: {cigar!r}")
    return runs


def _edit_runs(edit_script: Any) -> list[tuple[str, int]]:
    if edit_script in (None, []):
        return []
    if not isinstance(edit_script, list):
        raise IntervalError("edit_script must be a list")
    names = {
        "equal": "=",
        "match": "=",
        "m": "M",
        "substitution": "X",
        "mismatch": "X",
        "x": "X",
        "insertion": "I",
        "i": "I",
        "deletion": "D",
        "d": "D",
        "soft_clip": "S",
        "softclip": "S",
        "s": "S",
    }
    runs: list[tuple[str, int]] = []
    for index, item in enumerate(edit_script):
        if not isinstance(item, dict):
            raise IntervalError(f"edit_script[{index}] must be an object")
        operation = str(item.get("operation", "")).strip().lower()
        if operation not in names:
            raise IntervalError(f"unsupported edit operation: {operation!r}")
        length = _integer(item.get("length"), f"edit_script[{index}].length")
        if length <= 0:
            raise IntervalError("edit_script contains a zero-length operation")
        runs.append((names[operation], length))
    return runs


def _split_segments(value: Any, query_length: int) -> list[Interval]:
    raw = value
    if isinstance(raw, dict):
        raw = raw.get("segments", raw.get("query_segments", raw.get("interval", raw)))
    if isinstance(raw, dict):
        raw = [raw]
    if not isinstance(raw, (list, tuple)):
        raise IntervalError("query_segments must be a list")
    if query_length <= 0:
        raise IntervalError("query_length must be positive")
    result: list[Interval] = []
    for index, item in enumerate(raw):
        if isinstance(item, dict) and "interval" in item and not {"start", "end"}.issubset(item):
            item = item["interval"]
        if isinstance(item, dict):
            if "start" not in item or "end" not in item:
                raise IntervalError(f"query_segments[{index}] requires start and end")
            try:
                start = int(item["start"])
                end = int(item["end"])
            except (TypeError, ValueError) as exc:
                raise IntervalError(f"query_segments[{index}] requires integer coordinates") from exc
        elif isinstance(item, (list, tuple)) and len(item) == 2:
            try:
                start = int(item[0])
                end = int(item[1])
            except (TypeError, ValueError) as exc:
                raise IntervalError(f"query_segments[{index}] requires integer coordinates") from exc
        else:
            raise IntervalError(f"query_segments[{index}] must be an interval")
        if start < 0 or end < 0 or start > query_length or end > query_length:
            raise IntervalError("query segment exceeds query length")
        if start <= end:
            if start < end:
                result.append((start, end))
        else:
            # Accept a compact circular segment and make its two ordinary
            # pieces explicit before checking overlap.
            result.extend(((start, query_length), (0, end)))
    if interval_length(result) != sum(end - start for start, end in result):
        raise IntervalError("query segments overlap")
    return result


def project_alignment(
    query_segments: Any,
    query_length: int,
    *,
    cigar: str | None = None,
    edit_script: Any = None,
    query_deletion_operation: str = "D",
) -> dict[str, Any]:
    """Project alignment support onto query coordinates.

    ``query_deletion_operation`` is ``D`` for the Jam edit-script convention
    and ``I`` for SAM/PAF CIGARs.  In either case that operation consumes query
    bases which are aligned but not supported coverage.
    """

    if query_deletion_operation not in {"I", "D"}:
        raise IntervalError("query_deletion_operation must be I or D")
    segments = _split_segments(query_segments, query_length)
    runs = _edit_runs(edit_script) if edit_script not in (None, []) else _cigar_runs(cigar or "")
    if not runs:
        supported = list(segments)
        return {
            "query_segments": segments,
            "supported_intervals": supported,
            "aligned_intervals": list(segments),
            "alignment_deletions": [],
            "clipped_intervals": [],
            "clipped_bases": 0,
            "query_consumed": interval_length(segments),
        }
    available = sum(end - start for start, end in segments)
    cursor = 0
    supported: list[Interval] = []
    aligned: list[Interval] = []
    deletions: list[Interval] = []
    clipped: list[Interval] = []
    clipped_bases = 0

    def map_span(offset: int, length: int) -> list[Interval]:
        if offset < 0 or length < 0 or offset + length > available:
            raise IntervalError("edit script consumes more query sequence than segments")
        pieces: list[Interval] = []
        remaining_offset = offset
        remaining = length
        for start, end in segments:
            segment_length = end - start
            if remaining_offset >= segment_length:
                remaining_offset -= segment_length
                continue
            piece_start = start + remaining_offset
            take = min(remaining, end - piece_start)
            if take:
                pieces.append((piece_start, piece_start + take))
            remaining -= take
            remaining_offset = 0
            if not remaining:
                break
        if remaining:
            raise IntervalError("edit span could not be mapped to query segments")
        return pieces

    for operation, length in runs:
        consumes_query = operation in {"M", "=", "X", query_deletion_operation, "S"}
        if not consumes_query:
            continue
        if cursor + length > available:
            # Terminal soft clips are allowed outside the selected segment;
            # all other over-consumption is malformed.
            if operation == "S" and cursor >= available:
                clipped_bases += length
                continue
            raise IntervalError("edit script consumes more query sequence than segments")
        pieces = map_span(cursor, length)
        cursor += length
        if operation == "S":
            clipped.extend(pieces)
            clipped_bases += length
        elif operation == query_deletion_operation:
            aligned.extend(pieces)
            deletions.extend(pieces)
        else:
            aligned.extend(pieces)
            supported.extend(pieces)
    return {
        "query_segments": segments,
        "supported_intervals": union_intervals(supported),
        "aligned_intervals": union_intervals(aligned),
        "alignment_deletions": union_intervals(deletions),
        "clipped_intervals": union_intervals(clipped),
        "clipped_bases": clipped_bases,
        "query_consumed": cursor,
    }
