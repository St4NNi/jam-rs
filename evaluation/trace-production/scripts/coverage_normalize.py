#!/usr/bin/env python3
"""Normalize sequence-comparison intervals onto one checked query coordinate system.

The comparator runners use query coordinates as zero-based, half-open plasmid
coordinates.  BLAST tabular output is converted from its one-based inclusive
coordinates; PAF output is already zero-based half-open and may additionally
carry a ``cg:Z:`` CIGAR tag.  The same interval union and truth comparison is
used for both formats and for jam trace JSONL output.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Iterable


SCHEMA_VERSION = "1.0.0"
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
    """Project supported query bases from a PAF CIGAR.

    PAF's query is the plasmid.  Match, equal, and substitution operations
    support query bases.  Query insertions and soft clips advance the query
    cursor but are not supported by the target; target-only deletion/skip
    operations do not advance it.
    """

    if not cigar or cigar == "*":
        return _checked_interval(qstart, qend, query_length)
    cursor = qstart
    consumed = 0
    supported: list[dict[str, int]] = []
    offset = 0
    for match in _CIGAR_RE.finditer(cigar):
        if match.start() != offset:
            raise NormalizationError(f"invalid CIGAR at byte {offset}: {cigar!r}")
        length = int(match.group(1))
        operation = match.group(2)
        if length <= 0:
            raise NormalizationError(f"CIGAR contains a zero-length operation: {cigar!r}")
        if operation in {"M", "=", "X", "I", "S"}:
            consumed += length
            if operation in {"M", "=", "X"}:
                supported.extend(_checked_interval(cursor, cursor + length, query_length))
            cursor += length
        elif operation in {"D", "N", "H", "P"}:
            pass
        offset = match.end()
    if offset != len(cigar):
        raise NormalizationError(f"invalid CIGAR at byte {offset}: {cigar!r}")
    if cursor != qend:
        raise NormalizationError(
            f"CIGAR query span ends at {cursor}, but PAF declares {qend}: {cigar!r}"
        )
    if consumed < qend - qstart:
        raise NormalizationError("CIGAR query consumption is inconsistent with PAF span")
    return supported


def parse_blast6(
    text: str,
    query_id: str,
    query_length: int,
    metagenome_id: str,
    min_pident: float = 0.0,
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
        intervals = _checked_interval(min(qstart, qend) - 1, max(qstart, qend), query_length)
        rows.append(
            {
                "metagenome_id": metagenome_id,
                "contig_id": sseqid,
                "strand": "reverse" if qend < qstart else "forward",
                "intervals": intervals,
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
        intervals = parse_cigar_query_intervals(cigar, qstart, qend, query_length)
        rows.append(
            {
                "metagenome_id": metagenome_id,
                "contig_id": tname,
                "strand": "reverse" if strand == "-" else "forward",
                "intervals": intervals,
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


def _truth_intervals(
    path: Path | None, metagenome_id: str | None, query_length: int
) -> list[dict[str, int]]:
    if path is None:
        return []
    if path.suffix.lower() == ".json":
        value = json.loads(path.read_text(encoding="utf-8"))
        rows = value.get("intervals", value) if isinstance(value, dict) else value
        if not isinstance(rows, list):
            raise NormalizationError("truth JSON must be a list or an object with intervals")
    else:
        with path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    intervals: list[dict[str, int]] = []
    for row in rows:
        if not isinstance(row, dict):
            raise NormalizationError("truth records must be objects")
        row_id = row.get("metagenome_id") or row.get("assembly_id") or row.get("sample_id")
        if metagenome_id is not None and row_id not in {None, metagenome_id}:
            continue
        nested = row.get("intervals")
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
    records: list[dict], query_length: int, truth: list[dict[str, int]] | None = None
) -> dict:
    intervals = [piece for record in records for piece in record.get("intervals", [])]
    supported = union_intervals(intervals, query_length)
    result: dict = {
        "supported_intervals": supported,
        "supported_bases": _length(supported),
        "supported_fraction": _length(supported) / query_length,
        "gaps": complement(supported, query_length),
    }
    if truth is None:
        result["truth_status"] = "unavailable"
    else:
        truth_union = union_intervals(truth, query_length)
        result["truth_intervals"] = truth_union
        result["metrics"] = interval_metrics(supported, truth_union, query_length)
        result["truth_status"] = "measured"
    return result


def normalize_jam_jsonl(
    text: str, query_id: str | None = None, metagenome_id: str | None = None
) -> tuple[str, int, list[dict]]:
    records = [json.loads(line) for line in text.splitlines() if line.strip()]
    header = next((record for record in records if record.get("record_type") == "run_header"), None)
    if header is None:
        raise NormalizationError("jam JSONL has no run_header record")
    selected_query = query_id or header.get("plasmid_id")
    query_length = int(header.get("plasmid_length", 0))
    if not selected_query or query_length <= 0:
        raise NormalizationError("jam JSONL header lacks plasmid id or length")
    output: list[dict] = []
    for record in records:
        if record.get("record_type") != "metagenome_result":
            continue
        current_id = record.get("metagenome_id")
        if metagenome_id is not None and current_id != metagenome_id:
            continue
        for interval in (record.get("coverage") or {}).get("primary_intervals", []):
            for piece in _checked_interval(
                int(interval["start"]), int(interval["end"]), query_length
            ):
                output.append(
                    {
                        "metagenome_id": current_id,
                        "contig_id": "jam-primary",
                        "strand": "unknown",
                        "intervals": [piece],
                        "source": "jam",
                    }
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
) -> dict:
    text = input_path.read_text(encoding="utf-8")
    if input_format == "blast6":
        records = parse_blast6(text, query_id, query_length, metagenome_id, min_pident)
    elif input_format == "paf":
        records = parse_paf(text, query_id, query_length, metagenome_id, min_mapq)
    elif input_format == "jam":
        selected, discovered_length, records = normalize_jam_jsonl(text, query_id, metagenome_id)
        if selected != query_id or discovered_length != query_length:
            raise NormalizationError("jam JSONL query identity or length disagrees with arguments")
    else:
        raise NormalizationError(f"unsupported input format: {input_format}")
    truth = _truth_intervals(truth_path, metagenome_id, query_length) if truth_path else None
    return {
        "schema_version": SCHEMA_VERSION,
        "format": input_format,
        "query_id": query_id,
        "query_length": query_length,
        "metagenome_id": metagenome_id,
        "records": records,
        "coverage": normalize_records(records, query_length, truth),
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
    parser.add_argument("--format", choices=("blast6", "paf", "jam"), required=True)
    parser.add_argument("--query-id")
    parser.add_argument("--query-length", type=int)
    parser.add_argument("--query-fasta", type=Path)
    parser.add_argument("--metagenome-id")
    parser.add_argument("--truth", type=Path)
    parser.add_argument("--min-pident", type=float, default=0.0)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        if args.format == "jam":
            selected, query_length, records = normalize_jam_jsonl(
                args.input.read_text(encoding="utf-8"), args.query_id, args.metagenome_id
            )
            query_id = args.query_id or selected
            if args.query_length is not None and args.query_length != query_length:
                raise NormalizationError("--query-length disagrees with jam JSONL header")
            truth = _truth_intervals(args.truth, args.metagenome_id, query_length) if args.truth else None
            result = {
                "schema_version": SCHEMA_VERSION,
                "format": args.format,
                "query_id": query_id,
                "query_length": query_length,
                "metagenome_id": args.metagenome_id or "all",
                "records": records,
                "coverage": normalize_records(records, query_length, truth),
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
