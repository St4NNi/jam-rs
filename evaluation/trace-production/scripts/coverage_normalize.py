#!/usr/bin/env python3
"""Normalize alignment evidence onto one checked query coordinate system.

All normalized support is represented as zero-based, half-open query
intervals.  BLASTn and LexicMap report one-based inclusive query coordinates;
PAF reports zero-based half-open coordinates; and jam trace JSON uses explicit
query segments.  A single projection routine handles direct spans and edit
scripts so supported bases, aligned deletions, and unsupported gaps remain
separate evidence components.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any, Iterable, Sequence


SCHEMA_VERSION = "1.0.0"
NORMALIZATION_SCHEMA_VERSION = "normalized-mge-v1"
_CIGAR_RE = re.compile(r"([0-9]+)([MIDNSHP=X])")


class NormalizationError(ValueError):
    """Raised when comparator coordinates or evidence are malformed."""


def read_fasta_length(path: Path, query_id: str | None = None) -> tuple[str, int]:
    """Read one FASTA record and return its identifier and sequence length."""

    name: str | None = None
    length = 0
    with path.open(encoding="ascii") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split(None, 1)[0]
                if name is None or (query_id is not None and current == query_id):
                    name = current
                    length = 0
                elif query_id is None:
                    raise NormalizationError(
                        f"query FASTA contains more than one record: {path}"
                    )
            elif name is not None and (query_id is None or name == query_id):
                length += len("".join(line.split()))
    if name is None:
        raise NormalizationError(f"query FASTA has no selected record: {path}")
    if query_id is not None and name != query_id:
        raise NormalizationError(f"query id {query_id!r} is absent from {path}")
    if length <= 0:
        raise NormalizationError(f"query FASTA record {name!r} is empty")
    return name, length


def _checked_interval(start: int, end: int, query_length: int) -> list[dict[str, int]]:
    """Validate a circular interval and return ordinary half-open pieces."""

    if query_length <= 0:
        raise NormalizationError("query length must be positive")
    if start < 0 or end < 0 or start > query_length or end > query_length:
        raise NormalizationError(
            f"interval [{start}, {end}) lies outside query length {query_length}"
        )
    if start <= end:
        return [{"start": start, "end": end}] if start != end else []
    return [
        {"start": start, "end": query_length},
        {"start": 0, "end": end},
    ]


def _query_segments(value: Any, query_length: int) -> list[dict[str, int]]:
    """Decode one or more query segments without assigning topology.

    The caller may provide a single interval, a list of intervals, or records
    containing an ``interval`` member.  A wrapping interval is split at the
    query origin, but the function does not conclude that the query is
    circular; that fact must come from an explicit input field.
    """

    if isinstance(value, dict):
        value = value.get("segments", value.get("query_segments", value.get("interval", value)))
    if isinstance(value, dict):
        value = [value]
    if not isinstance(value, (list, tuple)):
        raise NormalizationError("query segments must be an interval or list of intervals")
    segments: list[dict[str, int]] = []
    for item in value:
        if isinstance(item, dict) and "interval" in item:
            item = item["interval"]
        if not isinstance(item, dict):
            raise NormalizationError("query segment must be an object")
        try:
            start = int(item["start"])
            end = int(item["end"])
        except (KeyError, TypeError, ValueError) as exc:
            raise NormalizationError("query segment requires integer start and end") from exc
        segments.extend(_checked_interval(start, end, query_length))
    if not segments:
        raise NormalizationError("query alignment has no non-empty query segment")
    # Query segments may be ordered around an origin.  They must not overlap,
    # because a CIGAR cursor cannot unambiguously map to two query bases.
    if _length(union_intervals(segments, query_length)) != _length(segments):
        raise NormalizationError("query alignment segments overlap")
    return segments


def _map_query_span(
    segments: Sequence[dict[str, int]], offset: int, length: int, query_length: int
) -> list[dict[str, int]]:
    """Map a linear edit-script cursor onto ordered query segments."""

    if offset < 0 or length < 0:
        raise NormalizationError("query edit span cannot be negative")
    available = _length(segments)
    if offset > available or length > available - offset:
        raise NormalizationError(
            f"query edit span [{offset}, {offset + length}) exceeds {available} query bases"
        )
    pieces: list[dict[str, int]] = []
    remaining_offset = offset
    remaining = length
    for segment in segments:
        segment_length = segment["end"] - segment["start"]
        if remaining_offset >= segment_length:
            remaining_offset -= segment_length
            continue
        start = segment["start"] + remaining_offset
        take = min(remaining, segment["end"] - start)
        if take:
            pieces.append({"start": start, "end": start + take})
        remaining = remaining - take
        remaining_offset = 0
        if not remaining:
            break
    if remaining:
        raise NormalizationError("query edit span could not be mapped to query segments")
    return [piece for piece in pieces if piece["start"] < piece["end"]]


def _cigar_runs(cigar: str) -> list[tuple[str, int]]:
    if not cigar or cigar == "*":
        return []
    runs: list[tuple[str, int]] = []
    offset = 0
    for match in _CIGAR_RE.finditer(cigar):
        if match.start() != offset:
            raise NormalizationError(f"invalid CIGAR at byte {offset}: {cigar!r}")
        length = int(match.group(1))
        if length <= 0:
            raise NormalizationError(f"CIGAR contains a zero-length operation: {cigar!r}")
        runs.append((match.group(2), length))
        offset = match.end()
    if offset != len(cigar):
        raise NormalizationError(f"invalid CIGAR at byte {offset}: {cigar!r}")
    return runs


def _edit_script_runs(edit_script: Any) -> list[tuple[str, int]]:
    if edit_script is None:
        return []
    if not isinstance(edit_script, list):
        raise NormalizationError("edit_script must be a list")
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
    for item in edit_script:
        if not isinstance(item, dict):
            raise NormalizationError("edit_script entries must be objects")
        operation = str(item.get("operation", "")).strip().lower()
        try:
            length = int(item["length"])
        except (KeyError, TypeError, ValueError) as exc:
            raise NormalizationError("edit_script entries require integer length") from exc
        if operation not in names:
            raise NormalizationError(f"unsupported edit operation: {operation!r}")
        if length <= 0:
            raise NormalizationError("edit_script contains a zero-length operation")
        runs.append((names[operation], length))
    return runs


def calculate_query_intervals(
    query_segments: Any,
    query_length: int,
    *,
    cigar: str | None = None,
    edit_script: Any = None,
    query_deletion_operation: str = "D",
) -> dict[str, Any]:
    """Project accepted alignment evidence onto query coordinates.

    ``query_deletion_operation`` identifies the operation that consumes query
    bases while aligning them to a target gap.  Jam's edit script uses ``D``
    for this operation.  SAM/PAF/LexicMap CIGAR conventions use ``I`` because
    the query is inserted relative to the target.  In either convention such
    bases are aligned deletions: they are part of the aligned span but are not
    supported coverage.  Unsupported gaps are calculated later as the
    complement of supported coverage.
    """

    if query_deletion_operation not in {"I", "D"}:
        raise NormalizationError("query_deletion_operation must be I or D")
    segments = _query_segments(query_segments, query_length)
    runs = (
        _edit_script_runs(edit_script)
        if edit_script not in (None, [])
        else _cigar_runs(cigar or "")
    )
    if not runs:
        supported = list(segments)
        return {
            "query_segments": segments,
            "supported_intervals": supported,
            "aligned_intervals": list(segments),
            "alignment_deletions": [],
            "clipped_intervals": [],
            "clipped_bases": 0,
            "query_consumed": _length(segments),
        }

    supported: list[dict[str, int]] = []
    aligned: list[dict[str, int]] = []
    deletions: list[dict[str, int]] = []
    clipped: list[dict[str, int]] = []
    cursor = 0
    clipped_bases = 0
    support_operations = {"M", "=", "X"}
    for operation, length in runs:
        consumes_query = operation in support_operations or operation in {
            query_deletion_operation,
            "S",
        }
        if not consumes_query:
            # Target-only insertion/deletion/skip operations do not consume
            # any query coordinate and therefore cannot create query support.
            continue
        if operation == "S" and cursor + length > _length(segments):
            # Some external formats include terminal soft clips outside the
            # declared aligned query span.  Preserve their count without
            # inventing coordinates; the complement remains unsupported.
            clipped_bases += length
            cursor += length
            continue
        pieces = _map_query_span(segments, cursor, length, query_length)
        cursor += length
        if operation in support_operations:
            supported.extend(pieces)
            aligned.extend(pieces)
        elif operation == query_deletion_operation:
            aligned.extend(pieces)
            deletions.extend(pieces)
        else:
            clipped.extend(pieces)
            clipped_bases += length
    available = _length(segments)
    if cursor != available and cursor < available:
        raise NormalizationError(
            f"edit script consumes {cursor} query bases but segments contain {available}"
        )
    return {
        "query_segments": segments,
        "supported_intervals": supported,
        "aligned_intervals": aligned,
        "alignment_deletions": deletions,
        "clipped_intervals": clipped,
        "clipped_bases": clipped_bases,
        "query_consumed": cursor,
    }


def union_intervals(intervals: Iterable[dict[str, int]], query_length: int) -> list[dict[str, int]]:
    """Checked union of non-wrapping query intervals."""

    checked: list[tuple[int, int]] = []
    for item in intervals:
        if not isinstance(item, dict):
            raise NormalizationError("interval must be an object")
        try:
            start = int(item["start"])
            end = int(item["end"])
        except (KeyError, TypeError, ValueError) as exc:
            raise NormalizationError("interval requires integer start and end") from exc
        if start < 0 or end < start or end > query_length:
            raise NormalizationError(
                f"non-wrapping interval [{start}, {end}) exceeds query length {query_length}"
            )
        if start != end:
            checked.append((start, end))
    checked.sort()
    merged: list[list[int]] = []
    for start, end in checked:
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [{"start": start, "end": end} for start, end in merged]


def complement(intervals: list[dict[str, int]], query_length: int) -> list[dict[str, int]]:
    """Return unsupported linear query intervals."""

    gaps: list[dict[str, int]] = []
    cursor = 0
    for item in intervals:
        if item["start"] > cursor:
            gaps.append({"start": cursor, "end": item["start"]})
        cursor = max(cursor, item["end"])
    if cursor < query_length:
        gaps.append({"start": cursor, "end": query_length})
    return gaps


def _length(intervals: Iterable[dict[str, int]]) -> int:
    return sum(item["end"] - item["start"] for item in intervals)


def _intersection_length(
    left: list[dict[str, int]], right: list[dict[str, int]]
) -> int:
    total = 0
    i = 0
    j = 0
    while i < len(left) and j < len(right):
        start = max(left[i]["start"], right[j]["start"])
        end = min(left[i]["end"], right[j]["end"])
        if start < end:
            total += end - start
        if left[i]["end"] < right[j]["end"]:
            i += 1
        else:
            j += 1
    return total


def _overlap_count(
    observed: list[dict[str, int]], truth: list[dict[str, int]]
) -> tuple[int, int]:
    observed_hits = sum(
        any(
            max(item["start"], expected["start"])
            < min(item["end"], expected["end"])
            for expected in truth
        )
        for item in observed
    )
    truth_hits = sum(
        any(
            max(item["start"], expected["start"])
            < min(item["end"], expected["end"])
            for item in observed
        )
        for expected in truth
    )
    return observed_hits, truth_hits


def interval_metrics(
    observed: Iterable[dict[str, int]],
    truth: Iterable[dict[str, int]],
    query_length: int,
) -> dict[str, float | int | None]:
    """Calculate base and interval metrics using checked half-open intervals.

    ``gap_error_bases`` is the symmetric difference between truth support and
    observed support.  It therefore counts both unsupported truth bases and
    false supported bases, rather than treating an alignment gap as a base
    mismatch.  Interval precision and recall count intervals with any positive
    overlap after each side has been unioned.
    """

    observed_union = union_intervals(observed, query_length)
    truth_union = union_intervals(truth, query_length)
    observed_bases = _length(observed_union)
    truth_bases = _length(truth_union)
    intersection = _intersection_length(observed_union, truth_union)
    overlap_observed, overlap_truth = _overlap_count(observed_union, truth_union)
    union_bases = observed_bases + truth_bases - intersection
    return {
        "observed_bases": observed_bases,
        "truth_bases": truth_bases,
        "intersection_bases": intersection,
        "base_precision": intersection / observed_bases if observed_bases else None,
        "base_recall": intersection / truth_bases if truth_bases else None,
        "interval_precision": (
            overlap_observed / len(observed_union) if observed_union else None
        ),
        "interval_recall": overlap_truth / len(truth_union) if truth_union else None,
        "gap_error_bases": union_bases - intersection,
        "gap_error_fraction": (
            (union_bases - intersection) / query_length if query_length else None
        ),
    }


def parse_cigar_query_intervals(
    cigar: str, qstart: int, qend: int, query_length: int
) -> list[dict[str, int]]:
    """Return supported PAF query intervals while retaining deletion evidence."""

    projection = calculate_query_intervals(
        [{"start": qstart, "end": qend}],
        query_length,
        cigar=cigar,
        query_deletion_operation="I",
    )
    return projection["supported_intervals"]


def parse_blast6(
    text: str,
    query_id: str,
    query_length: int,
    metagenome_id: str,
    min_pident: float = 0.0,
    *,
    query_kind: str = "unknown",
    topology_requested: str = "unknown",
) -> list[dict]:
    """Parse BLAST outfmt 6 with qseqid/sseqid/pident/length/qlen/qstart/qend."""

    rows: list[dict] = []
    for line_number, raw in enumerate(text.splitlines(), 1):
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 9:
            raise NormalizationError(
                f"BLAST tabular line {line_number} has {len(fields)} fields; expected at least 9"
            )
        try:
            qseqid, sseqid = fields[0], fields[1]
            pident = float(fields[2])
            qlen = int(fields[4])
            qstart = int(fields[5])
            qend = int(fields[6])
            sstart = int(fields[7])
            send = int(fields[8])
        except ValueError as exc:
            raise NormalizationError(f"invalid BLAST line {line_number}: {raw!r}") from exc
        if qseqid != query_id:
            continue
        if qlen != query_length:
            raise NormalizationError(
                f"BLAST qlen {qlen} disagrees with query length {query_length}"
            )
        if pident < min_pident:
            continue
        if qstart < 1 or qend < 1 or qstart > query_length or qend > query_length:
            raise NormalizationError(f"BLAST query coordinates are out of range: {raw!r}")
        query_segments = _checked_interval(
            min(qstart, qend) - 1, max(qstart, qend), query_length
        )
        rows.append(
            {
                "metagenome_id": metagenome_id,
                "contig_id": sseqid,
                "strand": "reverse" if qend < qstart else "forward",
                "query_segments": query_segments,
                "intervals": query_segments,
                "aligned_intervals": query_segments,
                "alignment_deletions": [],
                "unsupported_intervals": [],
                "query_kind": query_kind,
                "topology_requested": topology_requested,
                "pident": pident,
                "query_start_1based": qstart,
                "query_end_1based": qend,
                "target_start_1based": sstart,
                "target_end_1based": send,
            }
        )
    return rows


def parse_paf(
    text: str,
    query_id: str,
    query_length: int,
    metagenome_id: str,
    min_mapq: int = 0,
    *,
    query_kind: str = "unknown",
    topology_requested: str = "unknown",
) -> list[dict]:
    """Parse PAF, projecting an optional ``cg:Z:`` tag onto query intervals."""

    rows: list[dict] = []
    for line_number, raw in enumerate(text.splitlines(), 1):
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 12:
            raise NormalizationError(
                f"PAF line {line_number} has {len(fields)} fields; expected at least 12"
            )
        try:
            qname = fields[0]
            qlen, qstart, qend = (int(value) for value in fields[1:4])
            strand = fields[4]
            tname = fields[5]
            tlen, tstart, tend = (int(value) for value in fields[6:9])
            matches, block_length, mapq = (int(value) for value in fields[9:12])
        except ValueError as exc:
            raise NormalizationError(f"invalid PAF line {line_number}: {raw!r}") from exc
        if qname != query_id:
            continue
        if qlen != query_length:
            raise NormalizationError(
                f"PAF qlen {qlen} disagrees with query length {query_length}"
            )
        if strand not in {"+", "-"}:
            raise NormalizationError(f"PAF strand is not + or -: {raw!r}")
        if not 0 <= qstart <= qend <= qlen:
            raise NormalizationError(f"PAF query coordinates are out of range: {raw!r}")
        if not 0 <= tstart <= tend <= tlen:
            raise NormalizationError(f"PAF target coordinates are out of range: {raw!r}")
        if mapq < min_mapq:
            continue
        cigar = None
        for tag in fields[12:]:
            if tag.startswith("cg:Z:"):
                cigar = tag[5:]
                break
        projection = calculate_query_intervals(
            [{"start": qstart, "end": qend}],
            query_length,
            cigar=cigar,
            query_deletion_operation="I",
        )
        rows.append(
            {
                "metagenome_id": metagenome_id,
                "contig_id": tname,
                "strand": "reverse" if strand == "-" else "forward",
                "query_segments": projection["query_segments"],
                "intervals": projection["supported_intervals"],
                "aligned_intervals": projection["aligned_intervals"],
                "alignment_deletions": projection["alignment_deletions"],
                "unsupported_intervals": projection["clipped_intervals"],
                "query_kind": query_kind,
                "topology_requested": topology_requested,
                "matches": matches,
                "block_length": block_length,
                "mapq": mapq,
                "query_start": qstart,
                "query_end": qend,
                "target_start": tstart,
                "target_end": tend,
                "cigar": cigar,
            }
        )
    return rows


def parse_blastn(
    text: str,
    query_id: str,
    query_length: int,
    metagenome_id: str,
    min_pident: float = 0.0,
    *,
    query_kind: str = "unknown",
    topology_requested: str = "unknown",
) -> list[dict]:
    """Parse BLASTn tabular output with either positional or named columns.

    The production comparator uses the positional outfmt 6 emitted by
    :func:`parse_blast6`.  BLAST's ``# Fields:`` header is also accepted so a
    caller can retain optional CIGAR/aligned-sequence columns without
    changing the coordinate calculator.
    """

    header: list[str] | None = None
    data: list[str] = []
    for raw in text.splitlines():
        line = raw.rstrip("\n")
        if not line.strip():
            continue
        if line.startswith("#"):
            fields_text = line.split(":", 1)[1] if ":" in line else ""
            if line.lower().startswith("# fields:"):
                header = [
                    re.sub(r"[^a-z0-9]+", "", item.lower())
                    for item in fields_text.split(",")
                ]
            continue
        data.append(line)
    if header is None:
        return parse_blast6(
            "\n".join(data),
            query_id,
            query_length,
            metagenome_id,
            min_pident,
            query_kind=query_kind,
            topology_requested=topology_requested,
        )

    aliases = {
        "query": ("qseqid", "query", "qname", "queryid", "queryaccver"),
        "target": ("sseqid", "subject", "sname", "subjectid", "subjectaccver"),
        "pident": ("pident", "identity"),
        "qlen": ("qlen", "querylength"),
        "qstart": ("qstart", "querystart"),
        "qend": ("qend", "queryend"),
        "sstart": ("sstart", "subjectstart"),
        "send": ("send", "subjectend"),
    }

    def index(name: str, required: bool = True) -> int | None:
        for candidate in aliases[name]:
            if candidate in header:
                return header.index(candidate)
        if required:
            raise NormalizationError(f"BLAST header has no {name} column")
        return None

    query_index = index("query")
    target_index = index("target")
    pident_index = index("pident")
    qlen_index = index("qlen")
    qstart_index = index("qstart")
    qend_index = index("qend")
    sstart_index = index("sstart")
    send_index = index("send")
    cigar_index = next(
        (header.index(name) for name in ("cigar", "cg", "cigarstring") if name in header),
        None,
    )
    rows: list[dict] = []
    for line_number, raw in enumerate(data, 1):
        fields = raw.split("\t")
        try:
            values = {
                "query": fields[query_index],
                "target": fields[target_index],
                "pident": float(fields[pident_index]),
                "qlen": int(fields[qlen_index]),
                "qstart": int(fields[qstart_index]),
                "qend": int(fields[qend_index]),
                "sstart": int(fields[sstart_index]),
                "send": int(fields[send_index]),
            }
            cigar = fields[cigar_index] if cigar_index is not None else None
        except (IndexError, TypeError, ValueError) as exc:
            raise NormalizationError(f"invalid BLAST header row {line_number}: {raw!r}") from exc
        if values["query"] != query_id:
            continue
        if values["qlen"] != query_length:
            raise NormalizationError(
                f"BLAST qlen {values['qlen']} disagrees with query length {query_length}"
            )
        if values["pident"] < min_pident:
            continue
        qstart, qend = values["qstart"], values["qend"]
        if qstart < 1 or qend < 1 or qstart > query_length or qend > query_length:
            raise NormalizationError(f"BLAST query coordinates are out of range: {raw!r}")
        projection = calculate_query_intervals(
            [{"start": min(qstart, qend) - 1, "end": max(qstart, qend)}],
            query_length,
            cigar=cigar,
            query_deletion_operation="I",
        )
        rows.append(
            {
                "metagenome_id": metagenome_id,
                "contig_id": values["target"],
                "strand": "reverse" if qend < qstart else "forward",
                "query_segments": projection["query_segments"],
                "intervals": projection["supported_intervals"],
                "aligned_intervals": projection["aligned_intervals"],
                "alignment_deletions": projection["alignment_deletions"],
                "unsupported_intervals": projection["clipped_intervals"],
                "query_kind": query_kind,
                "topology_requested": topology_requested,
                "pident": values["pident"],
                "query_start_1based": qstart,
                "query_end_1based": qend,
                "target_start_1based": values["sstart"],
                "target_end_1based": values["send"],
                "cigar": cigar,
            }
        )
    return rows


_LEXICMAP_COLUMNS = (
    "query",
    "qlen",
    "hits",
    "sgenome",
    "sseqid",
    "qcovgnm",
    "cls",
    "hsp",
    "qcovhsp",
    "alenhsp",
    "pident",
    "gaps",
    "qstart",
    "qend",
    "sstart",
    "send",
    "sstr",
    "slen",
    "evalue",
    "bitscore",
)


def _normalized_column(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def parse_lexicmap(
    text: str,
    query_id: str,
    query_length: int,
    metagenome_id: str,
    min_pident: float = 0.0,
    *,
    query_kind: str = "unknown",
    topology_requested: str = "unknown",
) -> list[dict]:
    """Parse LexicMap's tabular search output.

    LexicMap's documented default columns use one-based inclusive ``qstart``
    and ``qend``.  The parser accepts the header emitted by LexicMap and the
    same fixed column order when a header was removed.  Optional ``cigar`` or
    ``cg`` columns are projected with external-tool (PAF/SAM) semantics.
    """

    rows: list[dict] = []
    header: list[str] | None = None
    data: list[str] = []
    for raw in text.splitlines():
        line = raw.rstrip("\n")
        if not line.strip() or set(line.replace("\t", "").replace(" ", "")) <= {"-"}:
            continue
        fields = line.split("\t") if "\t" in line else line.split()
        normalized = [_normalized_column(item) for item in fields]
        if "query" in normalized and "qstart" in normalized and "qend" in normalized:
            header = normalized
            continue
        if line.startswith("#"):
            continue
        data.append(line)
    if header is None:
        header = list(_LEXICMAP_COLUMNS)

    def index(*names: str, required: bool = True) -> int | None:
        for name in names:
            normalized = _normalized_column(name)
            if normalized in header:
                return header.index(normalized)
        if required:
            raise NormalizationError(f"LexicMap output has no required column {names[0]!r}")
        return None

    query_index = index("query", "qseqid")
    qlen_index = index("qlen", "query_length")
    qstart_index = index("qstart")
    qend_index = index("qend")
    target_index = index("sseqid", "subject", "target")
    genome_index = index("sgenome", "genome", "metagenome_id", required=False)
    strand_index = index("sstr", "strand", required=False)
    pident_index = index("pident", "identity", required=False)
    cigar_index = index("cigar", "cg", "cigar_string", required=False)
    optional_names = {
        name: index(name, required=False)
        for name in ("qcovgnm", "qcovhsp", "alenhsp", "gaps", "evalue", "bitscore")
    }

    for line_number, raw in enumerate(data, 1):
        fields = raw.split("\t") if "\t" in raw else raw.split()
        try:
            current_query = fields[query_index]
            qlen = int(fields[qlen_index])
            qstart = int(fields[qstart_index])
            qend = int(fields[qend_index])
            target = fields[target_index]
        except (IndexError, TypeError, ValueError) as exc:
            raise NormalizationError(f"invalid LexicMap row {line_number}: {raw!r}") from exc
        if current_query != query_id:
            continue
        if qlen != query_length:
            raise NormalizationError(
                f"LexicMap qlen {qlen} disagrees with query length {query_length}"
            )
        if qstart < 1 or qend < 1 or qstart > query_length or qend > query_length:
            raise NormalizationError(f"LexicMap query coordinates are out of range: {raw!r}")
        pident = None
        if pident_index is not None:
            try:
                pident = float(fields[pident_index])
            except (IndexError, TypeError, ValueError) as exc:
                raise NormalizationError(f"invalid LexicMap identity: {raw!r}") from exc
            if pident < min_pident:
                continue
        cigar = fields[cigar_index] if cigar_index is not None and cigar_index < len(fields) else None
        projection = calculate_query_intervals(
            [{"start": min(qstart, qend) - 1, "end": max(qstart, qend)}],
            query_length,
            cigar=cigar,
            query_deletion_operation="I",
        )
        record: dict[str, Any] = {
            # The caller's assembly/resource identifier is the canonical
            # grouping key.  LexicMap's optional sgenome is retained as
            # source metadata rather than silently replacing that key.
            "metagenome_id": metagenome_id,
            "contig_id": target,
            "strand": (
                "reverse"
                if strand_index is not None
                and strand_index < len(fields)
                and fields[strand_index].strip() in {"-", "reverse"}
                else "forward"
            ),
            "query_segments": projection["query_segments"],
            "intervals": projection["supported_intervals"],
            "aligned_intervals": projection["aligned_intervals"],
            "alignment_deletions": projection["alignment_deletions"],
            "unsupported_intervals": projection["clipped_intervals"],
            "query_kind": query_kind,
            "topology_requested": topology_requested,
            "query_start_1based": qstart,
            "query_end_1based": qend,
            "cigar": cigar,
        }
        if pident is not None:
            record["pident"] = pident
        if genome_index is not None and genome_index < len(fields):
            record["source_genome_id"] = fields[genome_index]
        for name, column_index in optional_names.items():
            if column_index is not None and column_index < len(fields):
                record[name] = fields[column_index]
        rows.append(record)
    return rows


def _truth_intervals(
    path: Path | None,
    metagenome_id: str | None,
    query_length: int,
    query_id: str | None = None,
) -> list[dict[str, int]]:
    if path is None:
        return []
    if path.suffix.lower() == ".json":
        value = json.loads(path.read_text(encoding="utf-8"))
        rows = value.get("intervals", value) if isinstance(value, dict) else value
        if not isinstance(rows, list):
            raise NormalizationError("truth JSON must be a list or an object with intervals")
    else:
        csv.field_size_limit(1024 * 1024 * 1024)
        with path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    intervals: list[dict[str, int]] = []
    for row in rows:
        if not isinstance(row, dict):
            raise NormalizationError("truth records must be objects")
        row_id = row.get("metagenome_id") or row.get("assembly_id") or row.get("sample_id")
        if row_id is None and row.get("assembly_path"):
            row_id = Path(str(row["assembly_path"])).name
        if metagenome_id is not None:
            exact_id = row_id == metagenome_id
            stem_id = any(
                metagenome_id.endswith(suffix)
                and metagenome_id[: -len(suffix)] == row_id
                for suffix in (".fa", ".fasta", ".fna", ".fastq", ".fq")
            )
            if not exact_id and not stem_id:
                continue
        row_query_id = (
            row.get("query_id")
            or row.get("query_accession")
            or row.get("plasmid_id")
            or row.get("reference_id")
        )
        if query_id is not None and row_query_id not in {None, query_id}:
            continue
        declared_length = row.get("query_length")
        if declared_length is not None and str(declared_length).strip():
            try:
                if isinstance(declared_length, bool):
                    raise TypeError
                if isinstance(declared_length, int):
                    parsed_length = declared_length
                elif isinstance(declared_length, str):
                    parsed_length = int(declared_length.strip())
                else:
                    raise TypeError
            except (TypeError, ValueError) as exc:
                raise NormalizationError(f"invalid truth query_length: {row!r}") from exc
            if parsed_length != query_length:
                continue
        nested = next(
            (
                row.get(name)
                for name in ("intervals", "query_intervals", "truth_intervals")
                if isinstance(row.get(name), list)
            ),
            None,
        )
        if isinstance(nested, list):
            for interval in nested:
                if not isinstance(interval, dict):
                    raise NormalizationError(f"invalid truth interval: {row!r}")
                try:
                    intervals.extend(
                        _checked_interval(
                            int(interval["start"]), int(interval["end"]), query_length
                        )
                    )
                except (KeyError, TypeError, ValueError) as exc:
                    raise NormalizationError(f"invalid truth interval: {row!r}") from exc
            continue
        encoded_intervals = row.get("query_intervals_json")
        if encoded_intervals is not None:
            if isinstance(encoded_intervals, str):
                try:
                    encoded_intervals = json.loads(encoded_intervals)
                except json.JSONDecodeError as exc:
                    raise NormalizationError(f"invalid query_intervals_json: {row!r}") from exc
            if not isinstance(encoded_intervals, list):
                raise NormalizationError(f"query_intervals_json must be a list: {row!r}")
            for interval in encoded_intervals:
                if isinstance(interval, dict):
                    try:
                        intervals.extend(
                            _checked_interval(
                                int(interval["start"]), int(interval["end"]), query_length
                            )
                        )
                    except (KeyError, TypeError, ValueError) as exc:
                        raise NormalizationError(f"invalid query interval: {row!r}") from exc
                    continue
                if not isinstance(interval, (list, tuple)) or len(interval) != 2:
                    raise NormalizationError(f"invalid query interval: {row!r}")
                try:
                    start, end = (int(value) for value in interval)
                    intervals.extend(_checked_interval(start, end, query_length))
                except (TypeError, ValueError) as exc:
                    raise NormalizationError(f"invalid query interval: {row!r}") from exc
            continue
        start_value = next(
            (
                row.get(name)
                for name in (
                    "truth_start",
                    "interval_start",
                    "plasmid_start",
                    "query_start",
                    "start",
                )
            ),
            None,
        )
        end_value = next(
            (
                row.get(name)
                for name in (
                    "truth_end",
                    "interval_end",
                    "plasmid_end",
                    "query_end",
                    "end",
                )
            ),
            None,
        )
        if start_value is None or end_value is None:
            continue
        try:
            intervals.extend(
                {"start": item["start"], "end": item["end"]}
                for item in _checked_interval(int(start_value), int(end_value), query_length)
            )
        except (TypeError, ValueError) as exc:
            raise NormalizationError(f"invalid truth interval: {row!r}") from exc
    return intervals


def normalize_records(
    records: list[dict],
    query_length: int,
    truth: list[dict[str, int]] | None = None,
    *,
    query_kind: str = "unknown",
    topology_requested: str = "unknown",
    coordinate_model: str = "undetermined",
    topology_evidence: str = "undetermined",
    source_schema_version: str | None = None,
) -> dict:
    # Secondary/alternative alignments are retained in ``records`` but do not
    # contribute to the nonredundant coverage union.  External comparator rows
    # have no ``primary`` field and therefore retain the historical behavior.
    primary_records = [record for record in records if record.get("primary", True)]
    intervals = [piece for record in primary_records for piece in record.get("intervals", [])]
    supported = union_intervals(intervals, query_length)
    aligned = union_intervals(
        [
            piece
            for record in primary_records
            for piece in record.get("aligned_intervals", record.get("intervals", []))
        ],
        query_length,
    )
    alignment_deletions = union_intervals(
        [piece for record in primary_records for piece in record.get("alignment_deletions", [])],
        query_length,
    )
    unsupported = complement(supported, query_length)
    result: dict = {
        "query_kind": query_kind,
        "topology_requested": topology_requested,
        "coordinate_model": coordinate_model,
        "topology_evidence": topology_evidence,
        "supported_intervals": supported,
        "supported_bases": _length(supported),
        "supported_fraction": _length(supported) / query_length,
        "aligned_intervals": aligned,
        "aligned_bases": _length(aligned),
        "aligned_fraction": _length(aligned) / query_length,
        "alignment_deletions": alignment_deletions,
        "alignment_deletion_bases": _length(alignment_deletions),
        "unsupported_gaps": unsupported,
        # ``gaps`` is retained as the v1 normalized-output spelling.
        "gaps": unsupported,
    }
    if source_schema_version is not None:
        result["source_schema_version"] = source_schema_version
    if truth is None:
        result["truth_status"] = "unavailable"
    else:
        truth_union = union_intervals(truth, query_length)
        result["truth_intervals"] = truth_union
        result["metrics"] = interval_metrics(supported, truth_union, query_length)
        result["truth_status"] = "measured"
    return result


def _normalize_jam_jsonl(
    text: str, query_id: str | None = None, metagenome_id: str | None = None
) -> tuple[str, int, list[dict], dict[str, Any]]:
    records = [json.loads(line) for line in text.splitlines() if line.strip()]
    header = next((record for record in records if record.get("record_type") == "run_header"), None)
    if header is None:
        raise NormalizationError("jam JSONL has no run_header record")
    selected_query = query_id or header.get("query_id") or header.get("plasmid_id")
    query_length = int(header.get("query_length", header.get("plasmid_length", 0)))
    if not selected_query or query_length <= 0:
        raise NormalizationError("jam JSONL header lacks plasmid id or length")
    explicit_topology = header.get("topology_requested")
    if not isinstance(explicit_topology, str):
        candidate_topology = header.get("topology")
        explicit_topology = candidate_topology if isinstance(candidate_topology, str) else "unknown"
    context = {
        "source_schema_version": header.get("schema_version", "unknown"),
        "query_kind": header.get("query_kind", "unknown"),
        "topology_requested": explicit_topology,
        "coordinate_model": header.get("coordinate_model", "undetermined"),
        "topology_evidence": header.get("topology_evidence", "undetermined"),
    }
    output: list[dict] = []
    for record in records:
        if record.get("record_type") != "metagenome_result":
            continue
        current_id = record.get("metagenome_id")
        if metagenome_id is not None and current_id != metagenome_id:
            continue
        result_topology = record.get("topology_requested")
        if not isinstance(result_topology, str):
            candidate_topology = record.get("topology")
            result_topology = candidate_topology if isinstance(candidate_topology, str) else context["topology_requested"]
        result_context = {
            key: record.get(key, value)
            for key, value in context.items()
        }
        result_context["topology_requested"] = result_topology
        for alignment in record.get("alignments", []) or []:
            if not isinstance(alignment, dict):
                raise NormalizationError("jam alignment records must be objects")
            query_value = alignment.get("query_segments", alignment.get("query_interval"))
            if query_value is None:
                # Keep malformed accepted records visible to downstream
                # review, but they cannot contribute query-coordinate support.
                output.append(
                    {
                        "metagenome_id": current_id,
                        "contig_id": alignment.get("contig_id", "unknown"),
                        "strand": alignment.get("strand", "unknown"),
                        "intervals": [],
                        "aligned_intervals": [],
                        "alignment_deletions": [],
                        "unsupported_intervals": [],
                        "source": "jam",
                        "accepted": True,
                        "primary": bool(alignment.get("primary", True)),
                        "source_record": alignment,
                        **result_context,
                    }
                )
                continue
            projection = calculate_query_intervals(
                query_value,
                query_length,
                cigar=alignment.get("cigar"),
                edit_script=alignment.get("edit_script"),
                query_deletion_operation="D",
            )
            output.append(
                {
                    "metagenome_id": current_id,
                    "contig_id": alignment.get("contig_id", "unknown"),
                    "strand": alignment.get("strand", "unknown"),
                    "intervals": projection["supported_intervals"],
                    "aligned_intervals": projection["aligned_intervals"],
                    "alignment_deletions": projection["alignment_deletions"],
                    "unsupported_intervals": projection["clipped_intervals"],
                    "query_segments": projection["query_segments"],
                    "source": "jam",
                    "accepted": True,
                    "primary": bool(alignment.get("primary", True)),
                    "source_record": alignment,
                    **result_context,
                }
            )
        if record.get("alignments"):
            continue
        # v1 streams and no-alignment v2 records still carry the selected
        # primary mosaic/coverage.  Preserve that established fallback.
        coverage = record.get("coverage") or {}
        mosaic = record.get("primary_fragment_mosaic") or {}
        intervals = coverage.get("primary_intervals") or mosaic.get("covered_intervals") or []
        for interval in intervals:
            for piece in _checked_interval(
                int(interval["start"]), int(interval["end"]), query_length
            ):
                output.append(
                    {
                        "metagenome_id": current_id,
                        "contig_id": "jam-primary",
                        "strand": "unknown",
                        "query_segments": [piece],
                        "intervals": [piece],
                        "aligned_intervals": [piece],
                        "alignment_deletions": [],
                        "unsupported_intervals": [],
                        "source": "jam",
                        "accepted": True,
                        "primary": True,
                        **result_context,
                    }
                )
    return selected_query, query_length, output, context


def normalize_jam_jsonl(
    text: str, query_id: str | None = None, metagenome_id: str | None = None
) -> tuple[str, int, list[dict]]:
    """Compatibility wrapper for callers using the original three-tuple API."""

    selected_query, query_length, output, _context = _normalize_jam_jsonl(
        text, query_id, metagenome_id
    )
    return selected_query, query_length, output


def normalize_file(
    input_path: Path,
    input_format: str,
    query_id: str,
    query_length: int,
    metagenome_id: str,
    truth_path: Path | None = None,
    min_pident: float = 0.0,
    min_mapq: int = 0,
    *,
    query_kind: str = "unknown",
    topology_requested: str = "unknown",
) -> dict:
    text = input_path.read_text(encoding="utf-8")
    context: dict[str, Any] = {
        "query_kind": query_kind,
        "topology_requested": topology_requested,
        "coordinate_model": "undetermined",
        "topology_evidence": "undetermined",
    }
    if input_format in {"blast6", "blastn"}:
        parser = parse_blastn if input_format == "blastn" else parse_blast6
        records = parser(
            text,
            query_id,
            query_length,
            metagenome_id,
            min_pident,
            query_kind=query_kind,
            topology_requested=topology_requested,
        )
    elif input_format == "paf":
        records = parse_paf(
            text,
            query_id,
            query_length,
            metagenome_id,
            min_mapq,
            query_kind=query_kind,
            topology_requested=topology_requested,
        )
    elif input_format == "lexicmap":
        records = parse_lexicmap(
            text,
            query_id,
            query_length,
            metagenome_id,
            min_pident,
            query_kind=query_kind,
            topology_requested=topology_requested,
        )
    elif input_format == "jam":
        selected, discovered_length, records, jam_context = _normalize_jam_jsonl(
            text, query_id, metagenome_id
        )
        if selected != query_id or discovered_length != query_length:
            raise NormalizationError("jam JSONL query identity or length disagrees with arguments")
        context.update(jam_context)
    else:
        raise NormalizationError(f"unsupported input format: {input_format}")
    truth = (
        _truth_intervals(truth_path, metagenome_id, query_length, query_id)
        if truth_path
        else None
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "format": input_format,
        "normalization_schema": NORMALIZATION_SCHEMA_VERSION,
        "query_id": query_id,
        "query_length": query_length,
        "metagenome_id": metagenome_id,
        **context,
        "records": records,
        "coverage": normalize_records(records, query_length, truth, **context),
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument(
        "--format", choices=("blast6", "blastn", "paf", "lexicmap", "jam"), required=True
    )
    parser.add_argument("--query-id")
    parser.add_argument("--query-length", type=int)
    parser.add_argument("--query-fasta", type=Path)
    parser.add_argument("--metagenome-id")
    parser.add_argument("--truth", type=Path)
    parser.add_argument("--min-pident", type=float, default=0.0)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--query-kind", default="unknown")
    parser.add_argument("--topology-requested", default="unknown")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        if args.format == "jam":
            selected, query_length, records, context = _normalize_jam_jsonl(
                args.input.read_text(encoding="utf-8"), args.query_id, args.metagenome_id
            )
            query_id = args.query_id or selected
            if args.query_length is not None and args.query_length != query_length:
                raise NormalizationError("--query-length disagrees with jam JSONL header")
            if args.query_kind != "unknown":
                context["query_kind"] = args.query_kind
            if args.topology_requested != "unknown":
                context["topology_requested"] = args.topology_requested
            truth = (
                _truth_intervals(args.truth, args.metagenome_id, query_length, query_id)
                if args.truth
                else None
            )
            result = {
                "schema_version": SCHEMA_VERSION,
                "format": args.format,
                "normalization_schema": NORMALIZATION_SCHEMA_VERSION,
                "query_id": query_id,
                "query_length": query_length,
                "metagenome_id": args.metagenome_id or "all",
                **context,
                "records": records,
                "coverage": normalize_records(records, query_length, truth, **context),
                "input_sha256": _sha256(args.input),
            }
        else:
            if args.query_fasta:
                fasta_id, fasta_length = read_fasta_length(args.query_fasta, args.query_id)
                query_id = args.query_id or fasta_id
                query_length = args.query_length or fasta_length
                if query_length != fasta_length:
                    raise NormalizationError("--query-length disagrees with query FASTA")
            elif args.query_id and args.query_length:
                query_id, query_length = args.query_id, args.query_length
            else:
                raise NormalizationError(
                    "BLAST/PAF normalization requires --query-id and --query-length or --query-fasta"
                )
            result = normalize_file(
                args.input,
                args.format,
                query_id,
                query_length,
                args.metagenome_id,
                args.truth,
                args.min_pident,
                args.min_mapq,
                query_kind=args.query_kind,
                topology_requested=args.topology_requested,
            )
            result["input_sha256"] = _sha256(args.input)
    except (OSError, json.JSONDecodeError, NormalizationError) as exc:
        parser.error(str(exc))
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    else:
        sys.stdout.write(encoded)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
