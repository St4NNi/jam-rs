#!/usr/bin/env python3
"""Materialize coordinate-truth Track A trace benchmark cases.

This generator deliberately consumes a verified source manifest.  It never
creates a background sequence and never treats a copied background as a new
independent metagenome.  Controlled cases are recorded as derivatives of a
named source background, with the source checksum repeated in every status
and truth record.

The source manifest is intentionally small and JSON based.  The supported
shape is documented by ``manifests/controlled-source-v1.schema.json``.  A
catalog is a list of FASTA records (or one catalog FASTA with a record list),
and backgrounds are real assembly FASTA files with required SHA-256 values.

Examples::

    evaluation/trace-production/scripts/controlled.py \
        --source-root /data/verified-plasmids \
        --output tools/out/trace-production/track-a-j01749 \
        --query-length 2000 --identity 100 --coverage 100 \
        --fragments 1 --orientation forward --indel none --limit 1

    evaluation/trace-production/scripts/controlled.py --self-check

The default matrix is intentionally broader than a normal CI run.  Use
filters or ``--limit`` to materialize a bounded release slice; the output
manifest retains planned cells and explicit unsupported statuses for every
cell considered by that invocation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import re
import shutil
import sys
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence


SCRIPT_VERSION = "controlled-track-a-v1"
MATRIX_SCHEMA = "jam-trace-controlled-matrix-v1"
TRUTH_SCHEMA = "jam-trace-controlled-truth-v1"
SOURCE_SCHEMA = "jam-trace-controlled-source-v1"
SEED = 0x4A414D5452414345
DNA = "ACGT"
VALID_BASES = frozenset("ACGTN")
MIN_SEED = 21
ANNOTATED_LABELS = frozenset({"shared_mobile_element", "chromosome_shared", "repeat_rich"})
CSV_FIELD_LIMIT = 1 << 30

DEFAULT_QUERY_LENGTHS = (2_000, 10_000, 50_000, 100_000, 200_000)
DEFAULT_IDENTITIES = (100, 99, 97, 95, 90, 85)
DEFAULT_COVERAGES = (5, 10, 20, 50, 80, 100)
DEFAULT_FRAGMENT_COUNTS = (1, 2, 5, 20)
DEFAULT_ORIENTATIONS = ("forward", "reverse")
DEFAULT_INDELS = (
    "none",
    "short_insertion",
    "short_deletion",
    "long_insertion",
    "long_deletion",
)
DEFAULT_QUERY_KINDS = ("plasmid", "phage")
DEFAULT_TOPOLOGIES = ("linear", "circular", "unknown")
DEFAULT_TERMINAL_REPEATS = ("none", "direct", "inverted")
VALID_QUERY_KINDS = frozenset({"plasmid", "phage", "other", "unknown"})
VALID_TOPOLOGIES = frozenset({"linear", "circular", "unknown"})
VALID_TERMINAL_REPEATS = frozenset({"none", "direct", "inverted"})


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_int(seed: int, *parts: object) -> int:
    material = "|".join([str(seed), *(str(part) for part in parts)]).encode("utf-8")
    return int.from_bytes(hashlib.sha256(material).digest()[:8], "big")


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._") or "unnamed"


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(table)[::-1].upper()


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    seen: set[str] = set()
    with path.open(encoding="ascii") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    sequence = "".join(chunks).upper()
                    _validate_sequence(path, name, sequence)
                    records.append((name, sequence))
                name = line[1:].split()[0]
                if not name:
                    raise ValueError(f"empty FASTA identifier in {path}:{line_number}")
                if name in seen:
                    raise ValueError(f"duplicate FASTA identifier {name!r} in {path}")
                seen.add(name)
                chunks = []
            elif name is None:
                raise ValueError(f"FASTA sequence precedes header in {path}:{line_number}")
            else:
                chunks.append(line)
    if name is not None:
        sequence = "".join(chunks).upper()
        _validate_sequence(path, name, sequence)
        records.append((name, sequence))
    if not records:
        raise ValueError(f"FASTA is empty: {path}")
    return records


def _validate_sequence(path: Path, name: str, sequence: str) -> None:
    if not sequence:
        raise ValueError(f"empty sequence {name!r} in {path}")
    invalid = sorted(set(sequence) - VALID_BASES)
    if invalid:
        raise ValueError(f"unsupported bases {invalid!r} in {path}:{name}")


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset : offset + 80] + "\n")


@dataclass(frozen=True)
class SequenceSource:
    identifier: str
    path: Path
    sha256: str
    length: int
    labels: tuple[str, ...] = ()
    role: str = ""


@dataclass(frozen=True)
class SequenceRecord:
    identifier: str
    sequence: str
    source: SequenceSource
    source_record_index: int
    labels: tuple[str, ...] = ()
    query_kind: str = "unknown"
    topology: str = "unknown"
    terminal_repeat: str = "none"


@dataclass(frozen=True)
class SourceData:
    root: Path
    manifest_path: Path
    manifest_sha256: str
    manifest: dict[str, Any]
    catalog_files: tuple[SequenceSource, ...]
    catalog_records: tuple[SequenceRecord, ...]
    background_files: tuple[SequenceSource, ...]
    background_records: dict[str, tuple[tuple[str, str], ...]]


@dataclass(frozen=True)
class Cell:
    case_id: str
    query: SequenceRecord | None
    requested_query_length: int
    identity: int
    coverage: int
    fragments: int
    orientation: str
    indel: str
    background: SequenceSource | None
    label: str
    query_kind: str
    topology: str
    terminal_repeat: str


def _manifest_file(root: Path, explicit: Path | None) -> Path:
    if explicit is not None:
        candidate = explicit if explicit.is_absolute() else root / explicit
        if not candidate.is_file():
            raise SystemExit(f"source manifest does not exist: {candidate}")
        return candidate.resolve()
    for name in ("manifest.json", "dataset.json"):
        candidate = root / name
        if candidate.is_file():
            return candidate.resolve()
    raise SystemExit(f"source root has no manifest.json or dataset.json: {root}")


def _entry_path(root: Path, value: Any, field: str) -> Path:
    if not isinstance(value, str) or not value:
        raise SystemExit(f"source {field} path must be a non-empty string")
    path = Path(value)
    if not path.is_absolute():
        path = root / path
    path = path.resolve()
    if not path.is_file():
        raise SystemExit(f"source {field} file does not exist: {path}")
    return path


def _entry_list(manifest: dict[str, Any], singular: str, plural: str) -> list[dict[str, Any]]:
    value = manifest.get(plural)
    if value is None:
        value = manifest.get(singular)
    if value is None:
        return []
    if isinstance(value, str):
        return [{"path": value}]
    if isinstance(value, dict):
        return [value]
    if not isinstance(value, list) or not all(isinstance(item, dict) for item in value):
        raise SystemExit(f"source {plural} must be a list of objects")
    return list(value)


def _metadata_records(manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    values = manifest.get("source_records") or manifest.get("records") or []
    result: dict[str, dict[str, Any]] = {}
    if isinstance(values, list):
        for item in values:
            if isinstance(item, dict):
                identifier = item.get("id") or item.get("accession") or item.get("identifier")
                if isinstance(identifier, str) and identifier:
                    result[identifier] = item
    return result


def _metadata_value(item: dict[str, Any], key: str, allowed: frozenset[str], default: str) -> str:
    value = item.get(key, default)
    if not isinstance(value, str):
        return default
    value = value.lower()
    return value if value in allowed else default


def _make_source(root: Path, item: dict[str, Any], fallback_id: str, field: str) -> SequenceSource:
    path = _entry_path(root, item.get("path") or item.get("file"), field)
    identifier = str(item.get("id") or item.get("accession") or fallback_id)
    expected = item.get("sha256") or item.get("checksum")
    if not isinstance(expected, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", expected):
        raise SystemExit(f"source {field} {identifier!r} lacks a required SHA-256 checksum")
    actual = sha256_file(path)
    if actual.lower() != expected.lower():
        raise SystemExit(f"checksum mismatch for {field} {identifier!r}: {path}")
    labels = item.get("labels", ())
    if isinstance(labels, str):
        labels = (labels,)
    if not isinstance(labels, (list, tuple)):
        raise SystemExit(f"source {field} {identifier!r} labels must be a list")
    return SequenceSource(
        identifier=identifier,
        path=path,
        sha256=actual,
        length=int(item["length"]) if item.get("length") is not None else 0,
        labels=tuple(str(label) for label in labels),
        role=str(item.get("role", "")),
    )


def load_source(root: Path, explicit_manifest: Path | None = None) -> SourceData:
    root = root.resolve()
    if not root.is_dir():
        raise SystemExit(f"source root is not a directory: {root}")
    manifest_path = _manifest_file(root, explicit_manifest)
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise SystemExit(f"cannot read source manifest {manifest_path}: {error}") from error
    if not isinstance(manifest, dict):
        raise SystemExit(f"source manifest must be a JSON object: {manifest_path}")
    declared_schema = manifest.get("schema_version")
    if declared_schema not in (None, SOURCE_SCHEMA, "1.0.0"):
        raise SystemExit(f"unsupported source manifest schema: {declared_schema}")

    metadata = _metadata_records(manifest)
    catalog_items = _entry_list(manifest, "catalog", "catalog_files")
    if not catalog_items:
        catalog_path_value = manifest.get("catalog_path") or "catalog.fasta"
        catalog_items = [{"path": catalog_path_value}]
    catalog_files: list[SequenceSource] = []
    catalog_records: list[SequenceRecord] = []
    for index, item in enumerate(catalog_items):
        fallback = Path(str(item.get("path") or item.get("file") or f"catalog-{index}")).stem
        if not item.get("sha256"):
            metadata_item = metadata.get(str(item.get("id") or item.get("accession") or fallback), {})
            expected_file_hash = item.get("file_sha256") or metadata_item.get("file_sha256")
            if expected_file_hash:
                item = {**item, "sha256": expected_file_hash}
            elif len(catalog_items) == 1 and manifest.get("catalog_sha256"):
                item = {**item, "sha256": manifest["catalog_sha256"]}
        source = _make_source(root, item, fallback, "catalog")
        records = read_fasta(source.path)
        catalog_files.append(source)
        for record_index, (identifier, sequence) in enumerate(records):
            metadata_item = metadata.get(identifier, {})
            catalog_flag = metadata_item.get("catalog")
            if catalog_flag is False:
                continue
            expected_length = metadata_item.get("length")
            if expected_length is not None and int(expected_length) != len(sequence):
                raise SystemExit(f"length mismatch for catalog record {identifier!r}")
            expected_record_hash = metadata_item.get("sequence_sha256")
            if expected_record_hash and sha256_bytes(sequence.encode("ascii")) != expected_record_hash:
                raise SystemExit(f"sequence checksum mismatch for catalog record {identifier!r}")
            labels = metadata_item.get("labels", source.labels)
            if isinstance(labels, str):
                labels = (labels,)
            query_kind = _metadata_value(metadata_item, "query_kind", VALID_QUERY_KINDS, "unknown")
            topology = _metadata_value(metadata_item, "topology", VALID_TOPOLOGIES, "unknown")
            terminal_repeat = _metadata_value(
                metadata_item, "terminal_repeat", VALID_TERMINAL_REPEATS, "none"
            )
            catalog_records.append(
                SequenceRecord(
                    identifier=identifier,
                    sequence=sequence,
                    source=SequenceSource(
                        identifier=source.identifier,
                        path=source.path,
                        sha256=source.sha256,
                        length=len(sequence),
                        labels=tuple(str(label) for label in labels),
                        role=str(metadata_item.get("role", source.role)),
                    ),
                    source_record_index=record_index,
                    labels=tuple(str(label) for label in labels),
                    query_kind=query_kind,
                    topology=topology,
                    terminal_repeat=terminal_repeat,
                )
            )
    if not catalog_records:
        raise SystemExit("verified catalog contains no usable records")

    # Explicit ``backgrounds`` is required.  A generic ``assemblies`` list is
    # often already a spike-in result and cannot be assumed to be a clean real
    # metagenome background.
    background_items = _entry_list(manifest, "background", "backgrounds")
    if not background_items:
        raise SystemExit(
            "source manifest must provide explicit checksum-verified backgrounds; "
            "an assemblies list is not accepted as a clean background source"
        )
    background_files: list[SequenceSource] = []
    background_records: dict[str, tuple[tuple[str, str], ...]] = {}
    for index, item in enumerate(background_items):
        fallback = Path(str(item.get("path") or item.get("file") or f"background-{index}")).stem
        source = _make_source(root, item, fallback, "background")
        records = tuple(read_fasta(source.path))
        if source.length and source.length != sum(len(sequence) for _, sequence in records):
            raise SystemExit(f"length mismatch for background {source.identifier!r}")
        background_files.append(source)
        background_records[source.identifier] = records
    if len({item.identifier for item in catalog_files}) != len(catalog_files):
        raise SystemExit("catalog file identifiers are not unique")
    if len({item.identifier for item in background_files}) != len(background_files):
        raise SystemExit("background identifiers are not unique")
    return SourceData(
        root=root,
        manifest_path=manifest_path,
        manifest_sha256=sha256_file(manifest_path),
        manifest=manifest,
        catalog_files=tuple(catalog_files),
        catalog_records=tuple(catalog_records),
        background_files=tuple(background_files),
        background_records=background_records,
    )


def matrix_path() -> Path:
    return Path(__file__).resolve().parents[1] / "manifests" / "controlled-matrix-v1.json"


def load_matrix(path: Path | None = None) -> dict[str, Any]:
    path = path or matrix_path()
    try:
        matrix = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise SystemExit(f"cannot read matrix manifest {path}: {error}") from error
    if not isinstance(matrix, dict) or matrix.get("schema_version") != MATRIX_SCHEMA:
        raise SystemExit(f"unsupported controlled matrix schema: {path}")
    required = {
        "query_lengths",
        "identities",
        "coverages",
        "fragment_counts",
        "orientations",
        "indel_profiles",
    }
    missing = sorted(required - matrix.keys())
    if missing:
        raise SystemExit(f"controlled matrix misses axes: {', '.join(missing)}")
    for stratum in matrix.get("query_strata", []):
        if not isinstance(stratum, dict):
            raise SystemExit("query_strata entries must be objects")
        if stratum.get("query_kind") not in VALID_QUERY_KINDS:
            raise SystemExit("query_strata has an unsupported query_kind")
        if stratum.get("topology") not in VALID_TOPOLOGIES:
            raise SystemExit("query_strata has an unsupported topology")
        if stratum.get("terminal_repeat") not in VALID_TERMINAL_REPEATS:
            raise SystemExit("query_strata has an unsupported terminal_repeat")
    return matrix


def choose_query(
    records: Sequence[SequenceRecord],
    requested_length: int,
    query_kind: str = "plasmid",
    topology: str = "unknown",
    terminal_repeat: str = "none",
) -> SequenceRecord | None:
    eligible = [
        record
        for record in records
        if len(record.sequence) >= requested_length
        and (
            query_kind == "unknown"
            or record.query_kind == query_kind
            or (query_kind == "plasmid" and record.query_kind == "unknown")
        )
        and (topology == "unknown" or record.topology in {topology, "unknown"})
        and (
            terminal_repeat == "none"
            or record.terminal_repeat in {terminal_repeat, "none"}
        )
    ]
    if not eligible:
        return None
    return min(eligible, key=lambda record: (len(record.sequence), record.identifier))


def query_slice(record: SequenceRecord, requested_length: int, seed: int) -> tuple[str, int]:
    if requested_length > len(record.sequence):
        raise ValueError("requested query length exceeds source record")
    if requested_length == len(record.sequence):
        return record.sequence, 0
    start = stable_int(seed, "query-window", record.identifier, requested_length) % (
        len(record.sequence) - requested_length + 1
    )
    return record.sequence[start : start + requested_length], start


def circular_positions(sequence_length: int, start: int, length: int) -> list[int]:
    return [((start + offset) % sequence_length) for offset in range(length)]


def circular_slice(sequence: str, start: int, length: int) -> str:
    return "".join(sequence[position] for position in circular_positions(len(sequence), start, length))


def intervals_for_positions(positions: Sequence[int], sequence_length: int) -> list[list[int]]:
    if not positions:
        return []
    intervals: list[list[int]] = []
    begin = positions[0]
    previous = positions[0]
    for position in positions[1:]:
        if position != previous + 1:
            intervals.append([begin, previous + 1])
            begin = position
        previous = position
    intervals.append([begin, previous + 1])
    if intervals and intervals[0][0] == 0 and intervals[-1][1] == sequence_length:
        # Keep the two pieces separate: this explicitly records an origin
        # crossing instead of silently flattening a circular coordinate.
        return intervals
    return intervals


def requested_overlap(fragment_count: int, coverage_bases: int, query_length: int) -> int:
    if fragment_count <= 1:
        return 0
    minimum_piece = coverage_bases // fragment_count
    return min(1000, max(MIN_SEED, coverage_bases // 100), max(0, minimum_piece // 2))


def split_intervals(
    query_length: int,
    start: int,
    coverage_bases: int,
    fragment_count: int,
) -> tuple[list[tuple[int, int]], int]:
    if fragment_count < 1 or coverage_bases < fragment_count:
        return [], 0
    overlap = requested_overlap(fragment_count, coverage_bases, query_length)
    edges = [round(index * coverage_bases / fragment_count) for index in range(fragment_count + 1)]
    intervals: list[tuple[int, int]] = []
    for index in range(fragment_count):
        begin = edges[index] - (overlap if index else 0)
        end = edges[index + 1] + (overlap if index + 1 < fragment_count else 0)
        if end <= begin:
            return [], overlap
        intervals.append(((start + begin) % query_length, end - begin))
    return intervals, overlap


def split_linear_intervals(
    query_length: int,
    start: int,
    coverage_bases: int,
    fragment_count: int,
) -> tuple[list[tuple[int, int]], int]:
    """Partition a linear query without modulo or origin crossing."""
    if fragment_count < 1 or coverage_bases < fragment_count:
        return [], 0
    overlap = requested_overlap(fragment_count, coverage_bases, query_length)
    edges = [round(index * coverage_bases / fragment_count) for index in range(fragment_count + 1)]
    intervals: list[tuple[int, int]] = []
    for index in range(fragment_count):
        begin = edges[index] - (overlap if index else 0)
        end = edges[index + 1] + (overlap if index + 1 < fragment_count else 0)
        if end <= begin:
            return [], overlap
        fragment_start = start + begin
        fragment_end = start + end
        if fragment_start < 0 or fragment_end > query_length:
            return [], overlap
        intervals.append((fragment_start, end - begin))
    return intervals, overlap


def _substitution_count(length: int, identity: int) -> int:
    return min(length, max(0, int(round(length * (100 - identity) / 100))))


def mutate_fragment(
    sequence: str,
    source_positions: list[int],
    identity: int,
    indel: str,
    rng: random.Random,
) -> tuple[str, dict[str, Any]]:
    if len(sequence) != len(source_positions):
        raise ValueError("sequence and coordinate list differ")
    chars = list(sequence)
    substitutions = _substitution_count(len(chars), identity)
    candidates = list(range(len(chars)))
    rng.shuffle(candidates)
    mutation_events: list[dict[str, Any]] = []
    for offset in sorted(candidates[:substitutions]):
        reference = chars[offset]
        choices = [base for base in DNA if base != reference]
        alternate = choices[rng.randrange(len(choices))]
        chars[offset] = alternate
        mutation_events.append(
            {
                "type": "substitution",
                "oriented_offset": offset,
                "query_position": source_positions[offset],
                "reference_base": reference,
                "alternate_base": alternate,
            }
        )

    if indel == "none":
        return "".join(chars), {
            "events": mutation_events,
            "substitution_count": substitutions,
            "substitution_identity_percent": 100.0 * (len(chars) - substitutions) / len(chars),
        }
    if indel not in DEFAULT_INDELS:
        raise ValueError(f"unsupported indel profile: {indel}")
    event_length = 3 if indel.startswith("short_") else max(32, min(512, len(chars) // 10))
    if indel.endswith("deletion"):
        if event_length >= len(chars):
            raise ValueError("deletion is longer than the fragment")
        offset = stable_int(rng.getrandbits(64), "deletion-offset", len(chars), indel) % (
            len(chars) - event_length + 1
        )
        deleted_positions = source_positions[offset : offset + event_length]
        del chars[offset : offset + event_length]
        mutation_events.append(
            {
                "type": "deletion",
                "oriented_offset": offset,
                "length": event_length,
                "query_positions": deleted_positions,
            }
        )
    else:
        offset = stable_int(rng.getrandbits(64), "insertion-offset", len(chars), indel) % (len(chars) + 1)
        inserted = "".join(rng.choice(DNA) for _ in range(event_length))
        chars[offset:offset] = list(inserted)
        mutation_events.append(
            {
                "type": "insertion",
                "oriented_offset": offset,
                "length": event_length,
                "inserted_sequence": inserted,
            }
        )
    return "".join(chars), {
        "events": mutation_events,
        "substitution_count": substitutions,
        "substitution_identity_percent": 100.0 * (len(sequence) - substitutions) / len(sequence),
    }


def _axis_values(
    matrix: dict[str, Any],
    key: str,
    values: Sequence[Any] | None,
    default: Sequence[Any] = (),
) -> tuple[Any, ...]:
    if values is None:
        return tuple(matrix.get(key, default))
    return tuple(values)


def _case_id(
    query: SequenceRecord | None,
    query_length: int,
    identity: int,
    coverage: int,
    fragments: int,
    orientation: str,
    indel: str,
    background: SequenceSource | None,
    query_kind: str = "plasmid",
    topology: str = "unknown",
    terminal_repeat: str = "none",
) -> str:
    query_id = safe_name(query.identifier if query else "unavailable")
    background_id = safe_name(background.identifier if background else "unavailable")
    base = (
        f"q{query_length}__{query_id}__i{identity}__c{coverage}__n{fragments}__"
        f"{orientation}__{indel}__bg_{background_id}"
    )
    if query_kind == "plasmid" and topology == "unknown" and terminal_repeat == "none":
        return base
    return f"{base}__kind_{safe_name(query_kind)}__top_{safe_name(topology)}__repeat_{safe_name(terminal_repeat)}"


def iter_cells(
    source: SourceData,
    matrix: dict[str, Any],
    *,
    query_lengths: Sequence[int] | None = None,
    identities: Sequence[int] | None = None,
    coverages: Sequence[int] | None = None,
    fragments: Sequence[int] | None = None,
    orientations: Sequence[str] | None = None,
    indels: Sequence[str] | None = None,
    background_id: str | None = None,
    label: str = "standard",
    query_kinds: Sequence[str] | None = None,
    topologies: Sequence[str] | None = None,
    terminal_repeats: Sequence[str] | None = None,
) -> Iterator[Cell]:
    backgrounds = [item for item in source.background_files if background_id is None or item.identifier == background_id]
    if background_id is not None and not backgrounds:
        raise SystemExit(f"requested background is absent: {background_id}")
    if not backgrounds:
        return
    for requested_length, identity, coverage, fragment_count, orientation, indel, background, query_kind, topology, terminal_repeat in product(
        _axis_values(matrix, "query_lengths", query_lengths),
        _axis_values(matrix, "identities", identities),
        _axis_values(matrix, "coverages", coverages),
        _axis_values(matrix, "fragment_counts", fragments),
        _axis_values(matrix, "orientations", orientations),
        _axis_values(matrix, "indel_profiles", indels),
        backgrounds,
        _axis_values(matrix, "query_kinds", query_kinds, DEFAULT_QUERY_KINDS[:1]),
        _axis_values(matrix, "topologies", topologies, DEFAULT_TOPOLOGIES[-1:]),
        _axis_values(matrix, "terminal_repeats", terminal_repeats, DEFAULT_TERMINAL_REPEATS[:1]),
    ):
        query = choose_query(
            source.catalog_records,
            int(requested_length),
            str(query_kind),
            str(topology),
            str(terminal_repeat),
        )
        yield Cell(
            case_id=_case_id(
                query,
                int(requested_length),
                int(identity),
                int(coverage),
                int(fragment_count),
                str(orientation),
                str(indel),
                background,
                str(query_kind),
                str(topology),
                str(terminal_repeat),
            ),
            query=query,
            requested_query_length=int(requested_length),
            identity=int(identity),
            coverage=int(coverage),
            fragments=int(fragment_count),
            orientation=str(orientation),
            indel=str(indel),
            background=background,
            label=label,
            query_kind=str(query_kind),
            topology=str(topology),
            terminal_repeat=str(terminal_repeat),
        )


def cell_status(cell: Cell, query_length: int, overlap: int | None = None) -> tuple[str, str]:
    if cell.query is None:
        return "unsupported", "no_catalog_record_reaches_requested_query_length"
    if cell.query_kind not in VALID_QUERY_KINDS:
        return "unsupported", "unsupported_query_kind"
    if cell.topology not in VALID_TOPOLOGIES:
        return "unsupported", "unsupported_topology"
    if cell.terminal_repeat not in VALID_TERMINAL_REPEATS:
        return "unsupported", "unsupported_terminal_repeat"
    if cell.label == "origin_crossing" and cell.topology != "circular":
        return "unsupported", "origin_crossing_requires_circular_coordinate_model"
    if cell.label in ANNOTATED_LABELS:
        available_labels = set(cell.query.labels)
        if cell.background is not None:
            available_labels.update(cell.background.labels)
        if cell.label not in available_labels:
            return "unsupported", "requested_label_has_no_source_annotation"
    if cell.identity not in DEFAULT_IDENTITIES or not 0 < cell.coverage <= 100:
        return "unsupported", "axis_value_outside_declared_matrix"
    if cell.fragments < 1:
        return "unsupported", "fragment_count_must_be_positive"
    coverage_bases = math.ceil(query_length * cell.coverage / 100)
    if coverage_bases < cell.fragments:
        return "unsupported", "coverage_has_fewer_bases_than_fragments"
    intervals, requested = split_intervals(query_length, 0, coverage_bases, cell.fragments)
    if not intervals:
        return "unsupported", "fragment_intervals_could_not_be_partitioned"
    if cell.fragments > 1 and requested == 0:
        return "unsupported", "requested_overlap_would_exceed_fragment_size"
    minimum_fragment = min(length for _, length in intervals)
    if minimum_fragment < MIN_SEED:
        return "unsupported", "fragment_shorter_than_k21_seed"
    if cell.indel in {"long_insertion", "long_deletion"} and minimum_fragment < 64:
        return "unsupported", "long_indel_requires_at_least_64_query_bases_per_fragment"
    if cell.orientation not in DEFAULT_ORIENTATIONS:
        return "unsupported", "unsupported_orientation"
    if cell.indel not in DEFAULT_INDELS:
        return "unsupported", "unsupported_indel_profile"
    return "materialize", ""


def _truth_fields() -> list[str]:
    return [
        "schema_version",
        "case_id",
        "metagenome_id",
        "query_id",
        "query_kind",
        "topology_requested",
        "coordinate_model",
        "topology_evidence",
        "terminal_repeat_type",
        "terminal_repeat_length",
        "source_query_id",
        "source_query_sha256",
        "reference_id",
        "background_id",
        "background_source_sha256",
        "source_manifest_sha256",
        "truth_class",
        "query_length",
        "source_query_start",
        "source_query_end",
        "identity_target_percent",
        "substitution_identity_percent",
        "coverage_target_percent",
        "coverage_target_bases",
        "coverage_union_bases",
        "fragment_count",
        "fragment_id",
        "query_fragment_start",
        "query_fragment_end",
        "orientation",
        "origin_crossing",
        "indel_profile",
        "mutation_events_json",
        "query_intervals_json",
        "contig_id",
        "contig_start",
        "contig_end",
        "query_sha256",
        "contig_sha256",
        "background_reused",
        "label",
    ]


def _status_fields() -> list[str]:
    return [
        "schema_version",
        "case_id",
        "status",
        "reason",
        "metagenome_id",
        "query_id",
        "query_kind",
        "topology_requested",
        "terminal_repeat_type",
        "source_query_id",
        "source_query_sha256",
        "background_id",
        "background_source_sha256",
        "source_manifest_sha256",
        "requested_query_length",
        "available_query_length",
        "identity_target_percent",
        "coverage_target_percent",
        "fragment_count",
        "orientation",
        "indel_profile",
        "label",
        "materialized_assembly",
    ]


def _write_tsv(path: Path, fields: Sequence[str], rows: Iterable[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            normalized = {
                field: str(row.get(field, "")).lower() if isinstance(row.get(field, ""), bool) else row.get(field, "")
                for field in fields
            }
            writer.writerow(normalized)


def _background_records_for_output(
    background: Sequence[tuple[str, str]],
) -> Iterator[tuple[str, str]]:
    for name, sequence in background:
        yield name, sequence


def query_output_id(cell: Cell, query_length: int) -> str:
    base = f"{safe_name(cell.query.identifier if cell.query else 'unavailable')}__q{query_length}"
    if cell.query_kind == "plasmid" and cell.topology == "unknown" and cell.terminal_repeat == "none":
        return base
    return (
        f"{base}__kind_{safe_name(cell.query_kind)}__top_{safe_name(cell.topology)}"
        f"__repeat_{safe_name(cell.terminal_repeat)}"
    )


def apply_terminal_repeat(sequence: str, repeat_type: str) -> tuple[str, dict[str, Any]]:
    """Construct a deterministic terminal-repeat query without topology claims."""
    if repeat_type == "none":
        return sequence, {"type": "none", "length": 0, "intervals": []}
    if repeat_type not in {"direct", "inverted"}:
        raise ValueError(f"unsupported terminal repeat type: {repeat_type}")
    repeat_length = min(100, max(MIN_SEED, len(sequence) // 10))
    if repeat_length * 2 > len(sequence):
        repeat_length = len(sequence) // 2
    if repeat_length < MIN_SEED:
        raise ValueError("terminal repeat requires at least two k21-sized regions")
    tail = sequence[-repeat_length:]
    prefix = tail if repeat_type == "direct" else reverse_complement(tail)
    transformed = prefix + sequence[repeat_length:]
    return transformed, {
        "type": repeat_type,
        "length": repeat_length,
        "intervals": [[0, repeat_length], [len(sequence) - repeat_length, len(sequence)]],
    }


def materialize_cell(
    output: Path,
    source: SourceData,
    cell: Cell,
    query_sequence: str,
    source_query_start: int,
    background_records: Sequence[tuple[str, str]],
    background_reused: bool,
    seed: int,
    repeat_metadata: dict[str, Any] | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    query_length = len(query_sequence)
    coverage_bases = math.ceil(query_length * cell.coverage / 100)
    # A reverse case reverses each fragment, not the source query record.  The
    # recorded intervals always remain in forward query coordinates.
    if cell.topology == "circular":
        origin_start = query_length - max(1, coverage_bases // 3) if cell.coverage else 0
        selected_start = (
            origin_start
            if cell.label == "origin_crossing"
            else stable_int(seed, cell.case_id, "segment") % query_length
        )
        intervals, overlap = split_intervals(query_length, selected_start, coverage_bases, cell.fragments)
    else:
        max_start = query_length - coverage_bases
        selected_start = stable_int(seed, cell.case_id, "segment") % (max_start + 1)
        intervals, overlap = split_linear_intervals(
            query_length, selected_start, coverage_bases, cell.fragments
        )
    assembly_records = list(_background_records_for_output(background_records))
    truth_rows: list[dict[str, Any]] = []
    case_dir = output / "cases" / safe_name(cell.case_id)
    case_dir.mkdir(parents=True, exist_ok=False)
    for index, (fragment_start, fragment_length) in enumerate(intervals, 1):
        origin_crossing = cell.topology == "circular" and fragment_start + fragment_length > query_length
        if cell.topology == "circular":
            positions = circular_positions(query_length, fragment_start, fragment_length)
            exact = circular_slice(query_sequence, fragment_start, fragment_length)
        else:
            positions = list(range(fragment_start, fragment_start + fragment_length))
            exact = query_sequence[fragment_start : fragment_start + fragment_length]
        if cell.orientation == "reverse":
            exact = reverse_complement(exact)
            positions = list(reversed(positions))
        rng = random.Random(stable_int(seed, cell.case_id, "fragment", index))
        observed, mutations = mutate_fragment(exact, positions, cell.identity, cell.indel, rng)
        contig_id = f"trace_{safe_name(cell.case_id)}_fragment_{index:02d}"
        assembly_records.append((contig_id, observed))
        query_intervals = intervals_for_positions(sorted(set(positions)), query_length)
        repeat_length = int((repeat_metadata or {}).get("length", 0))
        truth_rows.append(
            {
                "schema_version": TRUTH_SCHEMA,
                "case_id": cell.case_id,
                "metagenome_id": safe_name(cell.case_id),
                "query_id": query_output_id(cell, query_length),
                "query_kind": cell.query_kind,
                "topology_requested": cell.topology,
                "coordinate_model": "wrap" if cell.topology == "circular" else "linear",
                # Controlled construction does not infer biological topology.
                "topology_evidence": "undetermined",
                "terminal_repeat_type": cell.terminal_repeat,
                "terminal_repeat_length": repeat_length,
                "source_query_id": cell.query.identifier if cell.query else "",
                "source_query_sha256": sha256_bytes(cell.query.sequence.encode("ascii")) if cell.query else "",
                "reference_id": cell.query.identifier if cell.query else "",
                "background_id": cell.background.identifier if cell.background else "",
                "background_source_sha256": cell.background.sha256 if cell.background else "",
                "source_manifest_sha256": source.manifest_sha256,
                "truth_class": "controlled_spike_in",
                "query_length": query_length,
                "source_query_start": source_query_start,
                "source_query_end": source_query_start + query_length,
                "identity_target_percent": cell.identity,
                "substitution_identity_percent": f"{mutations['substitution_identity_percent']:.6f}",
                "coverage_target_percent": cell.coverage,
                "coverage_target_bases": coverage_bases,
                "coverage_union_bases": coverage_bases,
                "fragment_count": cell.fragments,
                "fragment_id": index,
                "query_fragment_start": fragment_start,
                "query_fragment_end": fragment_start + fragment_length,
                "orientation": cell.orientation,
                "origin_crossing": origin_crossing,
                "indel_profile": cell.indel,
                "mutation_events_json": json.dumps(mutations["events"], sort_keys=True, separators=(",", ":")),
                "query_intervals_json": json.dumps(query_intervals, separators=(",", ":")),
                "contig_id": contig_id,
                "contig_start": 0,
                "contig_end": len(observed),
                "query_sha256": sha256_bytes(query_sequence.encode("ascii")),
                "contig_sha256": sha256_bytes(observed.encode("ascii")),
                "background_reused": background_reused,
                "label": cell.label,
            }
        )
    assembly_path = output / "assemblies" / f"{safe_name(cell.case_id)}.fasta"
    write_fasta(assembly_path, assembly_records)
    case_metadata = {
        "schema_version": TRUTH_SCHEMA,
        "case_id": cell.case_id,
        "metagenome_id": safe_name(cell.case_id),
        "query_id": query_output_id(cell, query_length),
        "query_length": query_length,
        "query_sha256": sha256_bytes(query_sequence.encode("ascii")),
        "query_source_id": cell.query.identifier if cell.query else None,
        "query_source_interval": [source_query_start, source_query_start + query_length],
        "background_id": cell.background.identifier if cell.background else None,
        "background_sha256": cell.background.sha256 if cell.background else None,
        "background_reused": background_reused,
        "assembly_path": str(assembly_path.relative_to(output)),
        "assembly_sha256": sha256_file(assembly_path),
        "coverage_target_percent": cell.coverage,
        "coverage_target_bases": coverage_bases,
        "fragment_count": cell.fragments,
        "overlap_bases": overlap,
        "orientation": cell.orientation,
        "identity_target_percent": cell.identity,
        "indel_profile": cell.indel,
        "label": cell.label,
        "query_kind": cell.query_kind,
        "topology_requested": cell.topology,
        "coordinate_model": "wrap" if cell.topology == "circular" else "linear",
        # Construction metadata is deliberately not a biological topology call.
        "topology_evidence": "undetermined",
        "terminal_repeat_type": cell.terminal_repeat,
        "terminal_repeat": repeat_metadata or {"type": "none", "length": 0, "intervals": []},
        "truth_rows": len(truth_rows),
        "source_manifest_sha256": source.manifest_sha256,
    }
    (case_dir / "case.json").write_text(json.dumps(case_metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (case_dir / "truth.json").write_text(json.dumps(truth_rows, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return truth_rows, case_metadata


def materialize(
    source: SourceData,
    matrix: dict[str, Any],
    matrix_manifest_path: Path,
    output: Path,
    *,
    seed: int,
    query_lengths: Sequence[int] | None,
    identities: Sequence[int] | None,
    coverages: Sequence[int] | None,
    fragments: Sequence[int] | None,
    orientations: Sequence[str] | None,
    indels: Sequence[str] | None,
    background_id: str | None,
    label: str,
    query_kinds: Sequence[str] | None,
    topologies: Sequence[str] | None,
    terminal_repeats: Sequence[str] | None,
    case_id: str | None,
    limit: int | None,
) -> dict[str, Any]:
    output = output.resolve()
    if label not in matrix.get("labels", []):
        raise SystemExit(f"unsupported controlled case label: {label}")
    if output.exists():
        raise SystemExit(f"refusing to reuse existing output directory: {output}")
    output.mkdir(parents=True)
    (output / "source").mkdir()
    (output / "queries").mkdir()
    (output / "assemblies").mkdir()
    (output / "cases").mkdir()

    # Copy each verified source object at most once.  Case assemblies are
    # explicitly controlled derivatives and carry the original background ID.
    copied_catalog: list[dict[str, Any]] = []
    for index, item in enumerate(source.catalog_files):
        target = output / "source" / f"catalog_{index:02d}_{safe_name(item.identifier)}{item.path.suffix or '.fasta'}"
        shutil.copy2(item.path, target)
        copied_catalog.append({"id": item.identifier, "path": str(target.relative_to(output)), "sha256": item.sha256})
    copied_backgrounds: list[dict[str, Any]] = []
    for item in source.background_files:
        target = output / "source" / f"background_{safe_name(item.identifier)}{item.path.suffix or '.fasta'}"
        shutil.copy2(item.path, target)
        copied_backgrounds.append({"id": item.identifier, "path": str(target.relative_to(output)), "sha256": item.sha256})

    all_cells = list(
        iter_cells(
            source,
            matrix,
            query_lengths=query_lengths,
            identities=identities,
            coverages=coverages,
            fragments=fragments,
            orientations=orientations,
            indels=indels,
            background_id=background_id,
            label=label,
            query_kinds=query_kinds,
            topologies=topologies,
            terminal_repeats=terminal_repeats,
        )
    )
    all_cells.sort(key=lambda cell: cell.case_id)
    if case_id is not None:
        all_cells = [cell for cell in all_cells if cell.case_id == case_id]
        if not all_cells:
            raise SystemExit(f"requested case ID is not in the selected matrix: {case_id}")
    if limit is not None:
        if limit < 1:
            raise SystemExit("--limit must be positive")
        all_cells = all_cells[:limit]

    query_cache: dict[tuple[str, int, str], tuple[str, int, dict[str, Any]]] = {}
    query_files: dict[str, dict[str, Any]] = {}
    status_rows: list[dict[str, Any]] = []
    truth_rows: list[dict[str, Any]] = []
    case_metadata: list[dict[str, Any]] = []
    background_use: dict[str, int] = {}
    for cell in all_cells:
        available_length = len(cell.query.sequence) if cell.query else ""
        status, reason = cell_status(cell, cell.requested_query_length)
        metagenome_id = safe_name(cell.case_id)
        if status == "materialize":
            assert cell.query is not None and cell.background is not None
            key = (cell.query.identifier, cell.requested_query_length, cell.terminal_repeat)
            if key not in query_cache:
                source_sequence, source_query_start = query_slice(
                    cell.query, cell.requested_query_length, seed
                )
                query_sequence, repeat_metadata = apply_terminal_repeat(
                    source_sequence, cell.terminal_repeat
                )
                query_cache[key] = (query_sequence, source_query_start, repeat_metadata)
            query_sequence, source_query_start, repeat_metadata = query_cache[key]
            query_key = query_output_id(cell, cell.requested_query_length)
            query_path = output / "queries" / f"{query_key}.fasta"
            if query_key not in query_files:
                write_fasta(query_path, [(query_key, query_sequence)])
                query_files[query_key] = {
                    "id": query_key,
                    "source_id": cell.query.identifier,
                    "source_sha256": sha256_bytes(cell.query.sequence.encode("ascii")),
                    "length": len(query_sequence),
                    "source_interval": [source_query_start, source_query_start + len(query_sequence)],
                    "sha256": sha256_bytes(query_sequence.encode("ascii")),
                    "path": str(query_path.relative_to(output)),
                    "query_kind": cell.query_kind,
                    "topology_requested": cell.topology,
                    "coordinate_model": "wrap" if cell.topology == "circular" else "linear",
                    "topology_evidence": "undetermined",
                    "terminal_repeat": repeat_metadata,
                }
            background_records = source.background_records[cell.background.identifier]
            background_use[cell.background.identifier] = background_use.get(cell.background.identifier, 0) + 1
            rows, metadata = materialize_cell(
                output,
                source,
                cell,
                query_sequence,
                source_query_start,
                background_records,
                background_reused=True,
                seed=seed,
                repeat_metadata=repeat_metadata,
            )
            truth_rows.extend(rows)
            case_metadata.append(metadata)
            assembly_path = metadata["assembly_path"]
        else:
            assembly_path = ""
        status_rows.append(
            {
                "schema_version": TRUTH_SCHEMA,
                "case_id": cell.case_id,
                "status": status,
                "reason": reason,
                "metagenome_id": metagenome_id,
                "query_id": query_output_id(cell, cell.requested_query_length) if cell.query else "",
                "query_kind": cell.query_kind,
                "topology_requested": cell.topology,
                "terminal_repeat_type": cell.terminal_repeat,
                "source_query_id": cell.query.identifier if cell.query else "",
                "source_query_sha256": sha256_bytes(cell.query.sequence.encode("ascii")) if cell.query else "",
                "background_id": cell.background.identifier if cell.background else "",
                "background_source_sha256": cell.background.sha256 if cell.background else "",
                "source_manifest_sha256": source.manifest_sha256,
                "requested_query_length": cell.requested_query_length,
                "available_query_length": available_length,
                "identity_target_percent": cell.identity,
                "coverage_target_percent": cell.coverage,
                "fragment_count": cell.fragments,
                "orientation": cell.orientation,
                "indel_profile": cell.indel,
                "label": cell.label,
                "materialized_assembly": assembly_path,
            }
        )

    _write_tsv(output / "truth.tsv", _truth_fields(), truth_rows)
    _write_tsv(output / "status.tsv", _status_fields(), status_rows)
    (output / "truth.jsonl").write_text(
        "".join(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n" for row in truth_rows),
        encoding="utf-8",
    )
    matrix_checksum = sha256_file(matrix_manifest_path)
    output_manifest = {
        "schema_version": TRUTH_SCHEMA,
        "generator": SCRIPT_VERSION,
        "seed": seed,
        "seed_hex": hex(seed),
        "matrix_schema_version": matrix.get("schema_version"),
        "matrix_manifest": str(matrix_manifest_path.resolve()),
        "matrix_manifest_sha256": matrix_checksum,
        "source_root": str(source.root),
        "source_manifest": str(source.manifest_path),
        "source_manifest_sha256": source.manifest_sha256,
        "source_dataset_id": source.manifest.get("dataset_id") or source.manifest.get("dataset_kind"),
        "source_catalog_files": copied_catalog,
        "source_background_files": copied_backgrounds,
        "background_policy": {
            "backgrounds_are_source_objects": True,
            "controlled_derivatives_are_not_independent_metagenomes": True,
            "background_use_counts": background_use,
        },
        "planned_cells": len(all_cells),
        "materialized_cells": sum(1 for row in status_rows if row["status"] == "materialize"),
        "unsupported_cells": sum(1 for row in status_rows if row["status"] != "materialize"),
        "truth_rows": len(truth_rows),
        "queries": sorted(query_files.values(), key=lambda value: value["path"]),
        "cases": sorted(case_metadata, key=lambda value: value["case_id"]),
        "status_path": "status.tsv",
        "truth_path": "truth.tsv",
        "generation_parameters": {
            "query_lengths": list(_axis_values(matrix, "query_lengths", query_lengths)),
            "identities": list(_axis_values(matrix, "identities", identities)),
            "coverages": list(_axis_values(matrix, "coverages", coverages)),
            "fragment_counts": list(_axis_values(matrix, "fragment_counts", fragments)),
            "orientations": list(_axis_values(matrix, "orientations", orientations)),
            "indel_profiles": list(_axis_values(matrix, "indel_profiles", indels)),
            "query_kinds": list(_axis_values(matrix, "query_kinds", query_kinds, DEFAULT_QUERY_KINDS[:1])),
            "topologies": list(_axis_values(matrix, "topologies", topologies, DEFAULT_TOPOLOGIES[-1:])),
            "terminal_repeats": list(
                _axis_values(matrix, "terminal_repeats", terminal_repeats, DEFAULT_TERMINAL_REPEATS[:1])
            ),
            "background_id": background_id,
            "label": label,
            "case_id": case_id,
            "limit": limit,
        },
        "limitations": [
            "controlled cases reuse named source backgrounds and are not independent real metagenomes",
            "identity is applied by deterministic substitutions before optional indels",
            "unsupported matrix cells remain in status.tsv instead of being silently omitted",
            "read-derived assemblies and natural independently supported positives are separate benchmark tracks",
            "query_kind and topology are controlled construction strata; topology_evidence remains undetermined",
            "terminal repeat cases alter the query sequence deterministically and preserve unique query-coordinate coverage",
        ],
    }
    (output / "manifest.json").write_text(json.dumps(output_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    output_checksums: list[dict[str, str]] = []
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name != "checksums.tsv":
            output_checksums.append({"path": str(path.relative_to(output)), "sha256": sha256_file(path)})
    _write_tsv(output / "checksums.tsv", ["path", "sha256"], output_checksums)
    return output_manifest


def validate_output(output: Path) -> dict[str, Any]:
    output = output.resolve()
    csv.field_size_limit(CSV_FIELD_LIMIT)
    manifest_path = output / "manifest.json"
    truth_path = output / "truth.tsv"
    status_path = output / "status.tsv"
    for path in (manifest_path, truth_path, status_path):
        if not path.is_file():
            raise SystemExit(f"controlled output is missing {path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != TRUTH_SCHEMA:
        raise SystemExit("controlled output has an unsupported schema")
    with truth_path.open(encoding="utf-8", newline="") as handle:
        truth_rows = list(csv.DictReader(handle, delimiter="\t"))
    with status_path.open(encoding="utf-8", newline="") as handle:
        status_rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(status_rows) != int(manifest.get("planned_cells", -1)):
        raise SystemExit("status row count does not match manifest")
    materialized = {row["case_id"] for row in status_rows if row["status"] == "materialize"}
    truth_cases = {row["case_id"] for row in truth_rows}
    if truth_cases != materialized:
        raise SystemExit("truth cases do not match materialized status cases")
    for row in truth_rows:
        query_sha = row["query_sha256"]
        contig_sha = row["contig_sha256"]
        if len(query_sha) != 64 or len(contig_sha) != 64:
            raise SystemExit(f"invalid sequence checksum in truth row {row['case_id']}")
        if int(row["contig_start"]) < 0 or int(row["contig_end"]) < int(row["contig_start"]):
            raise SystemExit(f"invalid contig coordinate in truth row {row['case_id']}")
        intervals = json.loads(row["query_intervals_json"])
        if any(not (0 <= int(begin) < int(end) <= int(row["query_length"])) for begin, end in intervals):
            raise SystemExit(f"query interval outside bounds in truth row {row['case_id']}")
    return {
        "output": str(output),
        "planned_cells": len(status_rows),
        "materialized_cells": len(materialized),
        "truth_rows": len(truth_rows),
        "unsupported_cells": len(status_rows) - len(materialized),
    }


def self_check() -> dict[str, Any]:
    """Run pure generator invariants without creating benchmark files."""
    sequence = "ACGTN" * 100
    if reverse_complement(reverse_complement(sequence)) != sequence:
        raise SystemExit("reverse-complement invariant failed")
    positions = circular_positions(len(sequence), len(sequence) - 7, 20)
    if positions[:7] != list(range(len(sequence) - 7, len(sequence))) or positions[7:] != list(range(13)):
        raise SystemExit("circular coordinate invariant failed")
    intervals, overlap = split_intervals(2_000, 1_900, 1_000, 5)
    if len(intervals) != 5 or overlap <= 0:
        raise SystemExit("overlapping fragment invariant failed")
    rng = random.Random(7)
    observed, evidence = mutate_fragment("A" * 200, list(range(200)), 90, "long_insertion", rng)
    if len(observed) <= 200 or evidence["substitution_count"] != 20:
        raise SystemExit("mutation invariant failed")
    linear_intervals, _ = split_linear_intervals(2_000, 100, 1_000, 5)
    if any(start < 0 or start + length > 2_000 for start, length in linear_intervals):
        raise SystemExit("linear coordinate invariant failed")
    repeat_input = "ACGT" * 80 + "AAAA" * 20 + "TGCA" * 20
    direct, direct_meta = apply_terminal_repeat(repeat_input, "direct")
    inverted, inverted_meta = apply_terminal_repeat(repeat_input, "inverted")
    if direct == repeat_input or inverted == repeat_input:
        raise SystemExit("terminal repeat construction invariant failed")
    if direct_meta["type"] != "direct" or inverted_meta["type"] != "inverted":
        raise SystemExit("terminal repeat metadata invariant failed")
    matrix = load_matrix()
    for key, expected in {
        "query_lengths": DEFAULT_QUERY_LENGTHS,
        "identities": DEFAULT_IDENTITIES,
        "coverages": DEFAULT_COVERAGES,
        "fragment_counts": DEFAULT_FRAGMENT_COUNTS,
        "orientations": DEFAULT_ORIENTATIONS,
    }.items():
        if tuple(matrix[key]) != expected:
            raise SystemExit(f"matrix axis changed unexpectedly: {key}")
    strata = matrix.get("query_strata", [])
    if not any(item.get("query_kind") == "plasmid" for item in strata):
        raise SystemExit("matrix lacks plasmid query stratum")
    if not any(item.get("query_kind") == "phage" for item in strata):
        raise SystemExit("matrix lacks phage query stratum")
    if not {item.get("topology") for item in strata}.issuperset({"linear", "circular", "unknown"}):
        raise SystemExit("matrix lacks topology strata")
    return {"status": "ok", "script_version": SCRIPT_VERSION, "matrix_schema": MATRIX_SCHEMA}


def _parse_int_list(value: str, name: str) -> list[int]:
    try:
        values = [int(item) for item in value.split(",") if item]
    except ValueError as error:
        raise argparse.ArgumentTypeError(f"{name} must be comma-separated integers") from error
    if not values:
        raise argparse.ArgumentTypeError(f"{name} must not be empty")
    return values


def _parse_string_list(value: str, name: str) -> list[str]:
    values = [item for item in value.split(",") if item]
    if not values:
        raise argparse.ArgumentTypeError(f"{name} must not be empty")
    return values


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-check", action="store_true", help="check pure generator invariants without writing files")
    parser.add_argument("--validate-output", type=Path, help="validate an existing generated output directory")
    parser.add_argument("--source-root", type=Path, help="checksum-verified source dataset root")
    parser.add_argument("--source-manifest", type=Path, help="source manifest path, relative to --source-root unless absolute")
    parser.add_argument("--output", type=Path, help="new output directory")
    parser.add_argument("--matrix", type=Path, help="controlled matrix manifest")
    parser.add_argument("--seed", type=lambda value: int(value, 0), default=SEED)
    parser.add_argument("--query-length", type=lambda value: _parse_int_list(value, "query-length"))
    parser.add_argument("--identity", type=lambda value: _parse_int_list(value, "identity"))
    parser.add_argument("--coverage", type=lambda value: _parse_int_list(value, "coverage"))
    parser.add_argument("--fragments", type=lambda value: _parse_int_list(value, "fragments"))
    parser.add_argument("--orientation", type=lambda value: _parse_string_list(value, "orientation"))
    parser.add_argument("--indel", type=lambda value: _parse_string_list(value, "indel"))
    parser.add_argument("--query-kind", type=lambda value: _parse_string_list(value, "query-kind"))
    parser.add_argument("--topology", type=lambda value: _parse_string_list(value, "topology"))
    parser.add_argument(
        "--terminal-repeat", type=lambda value: _parse_string_list(value, "terminal-repeat")
    )
    parser.add_argument("--background-id")
    parser.add_argument("--label", default="standard")
    parser.add_argument("--case-id")
    parser.add_argument("--limit", type=int)
    args = parser.parse_args(argv)

    if args.self_check:
        print(json.dumps(self_check(), sort_keys=True))
        return 0
    if args.validate_output is not None:
        print(json.dumps(validate_output(args.validate_output), sort_keys=True))
        return 0
    if args.source_root is None or args.output is None:
        parser.error("--source-root and --output are required unless --self-check or --validate-output is used")
    source = load_source(args.source_root, args.source_manifest)
    matrix = load_matrix(args.matrix)
    result = materialize(
        source,
        matrix,
        (args.matrix or matrix_path()).resolve(),
        args.output,
        seed=args.seed,
        query_lengths=args.query_length,
        identities=args.identity,
        coverages=args.coverage,
        fragments=args.fragments,
        orientations=args.orientation,
        indels=args.indel,
        background_id=args.background_id,
        label=args.label,
        query_kinds=args.query_kind,
        topologies=args.topology,
        terminal_repeats=args.terminal_repeat,
        case_id=args.case_id,
        limit=args.limit,
    )
    print(
        json.dumps(
            {
                "output": str(args.output.resolve()),
                "planned_cells": result["planned_cells"],
                "materialized_cells": result["materialized_cells"],
                "unsupported_cells": result["unsupported_cells"],
                "truth_rows": result["truth_rows"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"controlled Track A generation failed: {error}", file=sys.stderr)
        raise SystemExit(2) from error
