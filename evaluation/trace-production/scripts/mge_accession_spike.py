#!/usr/bin/env python3
"""Generate accession-pinned plasmid/phage trace cases on real backgrounds.

The input FASTA is never downloaded or changed.  It must contain exactly the
accession-versioned records selected from the accompanying manifest, and every
selected record must pass its declared length and sequence SHA-256 check.  The
generator appends declared reference-derived contigs to copies of the supplied
real chromosome or metagenome-background records.  It does not infer topology,
terminal repeats, shared regions, or biological origin: those facts must be
declared in the manifest.

Rows whose checksums or topology annotations still need live verification are
kept in the public manifest.  Generation refuses them by default.  The
explicit ``--verified-only`` mode can materialize the verified subset and
records the omitted rows in ``dataset.json``; it never treats them as negative
evidence.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO


SCHEMA_VERSION = "jam-trace-mge-accession-spike-v1"
DATASET_IDS = ("forward", "reverse", "origin_crossing", "terminal_repeat", "shared_region")
SOURCE_FIELDS = (
    "accession",
    "display_name",
    "query_kind",
    "topology",
    "role",
    "catalog",
    "background_group",
    "length",
    "sequence_sha256",
    "record_url",
    "verification_status",
    "terminal_repeat_start",
    "terminal_repeat_end",
    "shared_group",
    "shared_region_start",
    "shared_region_end",
)
TRUTH_FIELDS = (
    "case_id",
    "query_accession",
    "query_kind",
    "topology",
    "expected_reference",
    "truth_class",
    "construction",
    "background_group",
    "background_accessions_json",
    "source_length",
    "source_sha256",
    "source_start",
    "source_end",
    "query_intervals_json",
    "contig_id",
    "contig_start",
    "contig_end",
    "orientation",
    "origin_crossing",
    "terminal_repeat",
    "shared_group",
    "coverage_bases",
    "coverage_fraction",
    "assembly_path",
    "assembly_sha256",
)
STATUS_FIELDS = (
    "case_id",
    "query_accession",
    "query_kind",
    "topology",
    "case_type",
    "status",
    "reason",
    "source_sha256",
    "background_group",
)
VERIFICATION_FIELDS = SOURCE_FIELDS + ("observed_length", "observed_sha256", "verified")


@dataclass(frozen=True)
class Source:
    accession: str
    display_name: str
    query_kind: str
    topology: str
    role: str
    catalog: bool
    background_group: str | None
    length: int
    sequence_sha256: str
    record_url: str
    verification_status: str
    terminal_repeat_start: int | None
    terminal_repeat_end: int | None
    shared_group: str | None
    shared_region_start: int | None
    shared_region_end: int | None


@dataclass(frozen=True)
class Fragment:
    identifier: str
    sequence: str
    source_start: int | None
    source_end: int | None
    query_intervals: tuple[tuple[int, int], ...]
    orientation: str
    origin_crossing: bool
    terminal_repeat: bool
    shared_group: str | None


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path) -> TextIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="ascii", newline="")
    return path.open(encoding="ascii", newline="")


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier: str | None = None
    sequence: list[str] = []
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records[identifier] = finish_sequence(path, identifier, sequence)
                identifier = line[1:].split()[0]
                if not identifier:
                    raise ValueError(f"empty FASTA identifier in {path}:{line_number}")
                if identifier in records:
                    raise ValueError(f"duplicate FASTA identifier {identifier!r} in {path}")
                sequence = []
            elif identifier is None:
                raise ValueError(f"sequence before first FASTA header in {path}:{line_number}")
            else:
                sequence.append(line)
    if identifier is not None:
        records[identifier] = finish_sequence(path, identifier, sequence)
    if not records:
        raise ValueError(f"source FASTA is empty: {path}")
    return records


def finish_sequence(path: Path, identifier: str, chunks: list[str]) -> str:
    sequence = "".join(chunks).upper()
    if not sequence:
        raise ValueError(f"empty sequence {identifier!r} in {path}")
    invalid = sorted(set(sequence) - set("ACGTRYSWKMBDHVN"))
    if invalid:
        raise ValueError(f"unsupported bases {invalid!r} in {path}:{identifier}")
    return sequence


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset : offset + 80] + "\n")


def optional_int(value: str, field: str, accession: str) -> int | None:
    if value in {"", "NA", "na", "null"}:
        return None
    try:
        result = int(value)
    except ValueError as error:
        raise ValueError(f"invalid {field} for {accession}: {value}") from error
    if result < 0:
        raise ValueError(f"negative {field} for {accession}: {value}")
    return result


def read_sources(path: Path) -> list[Source]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != SOURCE_FIELDS:
            raise ValueError(f"unexpected MGE source manifest columns in {path}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"MGE source manifest is empty: {path}")
    sources: list[Source] = []
    seen: set[str] = set()
    for row in rows:
        if any(not row.get(field) for field in SOURCE_FIELDS):
            raise ValueError(f"MGE source manifest has an empty field: {row}")
        accession = row["accession"]
        if accession in seen:
            raise ValueError(f"duplicate MGE accession: {accession}")
        seen.add(accession)
        digest = row["sequence_sha256"].lower()
        if digest != "unverified" and (len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest)):
            raise ValueError(f"invalid sequence SHA-256 for {accession}")
        try:
            length = int(row["length"])
        except ValueError as error:
            raise ValueError(f"invalid sequence length for {accession}: {row['length']}") from error
        if length < 1:
            raise ValueError(f"sequence length must be positive for {accession}")
        catalog_value = row["catalog"].lower()
        if catalog_value not in {"true", "false"}:
            raise ValueError(f"catalog must be true or false for {accession}")
        status = row["verification_status"].lower()
        if status not in {"verified", "needs_live_verification"}:
            raise ValueError(f"unsupported verification status for {accession}: {status}")
        sources.append(
            Source(
                accession=accession,
                display_name=row["display_name"],
                query_kind=row["query_kind"],
                topology=row["topology"],
                role=row["role"],
                catalog=catalog_value == "true",
                background_group=None if row["background_group"] == "NA" else row["background_group"],
                length=length,
                sequence_sha256=digest,
                record_url=row["record_url"],
                verification_status=status,
                terminal_repeat_start=optional_int(row["terminal_repeat_start"], "terminal_repeat_start", accession),
                terminal_repeat_end=optional_int(row["terminal_repeat_end"], "terminal_repeat_end", accession),
                shared_group=None if row["shared_group"] == "NA" else row["shared_group"],
                shared_region_start=optional_int(row["shared_region_start"], "shared_region_start", accession),
                shared_region_end=optional_int(row["shared_region_end"], "shared_region_end", accession),
            )
        )
    return sources


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]


def source_intervals(start: int, end: int) -> tuple[tuple[int, int], ...]:
    return ((start, end),)


def fragments_for(source: Source, sequence: str) -> dict[str, Fragment | None]:
    length = len(sequence)
    reverse = Fragment(
        identifier=f"{source.accession}__reverse",
        sequence=reverse_complement(sequence),
        source_start=0,
        source_end=length,
        query_intervals=source_intervals(0, length),
        orientation="reverse",
        origin_crossing=False,
        terminal_repeat=False,
        shared_group=None,
    )
    forward = Fragment(
        identifier=f"{source.accession}__forward",
        sequence=sequence,
        source_start=0,
        source_end=length,
        query_intervals=source_intervals(0, length),
        orientation="forward",
        origin_crossing=False,
        terminal_repeat=False,
        shared_group=None,
    )
    origin: Fragment | None = None
    if source.topology == "circular":
        offset = max(1, length // 3)
        rotated = sequence[offset:] + sequence[:offset]
        origin = Fragment(
            identifier=f"{source.accession}__origin",
            sequence=rotated,
            source_start=None,
            source_end=None,
            query_intervals=((offset, length), (0, offset)),
            orientation="forward",
            origin_crossing=True,
            terminal_repeat=False,
            shared_group=None,
        )
    terminal: Fragment | None = None
    if source.terminal_repeat_start is not None and source.terminal_repeat_end is not None:
        start = source.terminal_repeat_start
        end = source.terminal_repeat_end
        if end <= start or end > length:
            raise ValueError(f"terminal repeat coordinates are outside {source.accession}")
        terminal = Fragment(
            identifier=f"{source.accession}__terminal_repeat",
            sequence=sequence[start:end],
            source_start=start,
            source_end=end,
            query_intervals=source_intervals(start, end),
            orientation="forward",
            origin_crossing=False,
            terminal_repeat=True,
            shared_group=None,
        )
    shared: Fragment | None = None
    if source.shared_group is not None and source.shared_region_start is not None and source.shared_region_end is not None:
        start = source.shared_region_start
        end = source.shared_region_end
        if end <= start or end > length:
            raise ValueError(f"shared-region coordinates are outside {source.accession}")
        shared = Fragment(
            identifier=f"{source.accession}__shared_region",
            sequence=sequence[start:end],
            source_start=start,
            source_end=end,
            query_intervals=source_intervals(start, end),
            orientation="forward",
            origin_crossing=False,
            terminal_repeat=False,
            shared_group=source.shared_group,
        )
    return {"forward": forward, "reverse": reverse, "origin_crossing": origin, "terminal_repeat": terminal, "shared_region": shared}


def _write_tsv(path: Path, fields: tuple[str, ...], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(
            {
                field: str(row.get(field, "")).lower() if isinstance(row.get(field, ""), bool) else row.get(field, "")
                for field in fields
            }
            for row in rows
        )


def generate(sources_fasta: Path, source_manifest: Path, output: Path, verified_only: bool) -> dict[str, object]:
    sources_fasta = sources_fasta.resolve()
    source_manifest = source_manifest.resolve()
    output = output.resolve()
    if not sources_fasta.is_file():
        raise SystemExit(f"source FASTA does not exist: {sources_fasta}")
    if not source_manifest.is_file():
        raise SystemExit(f"source manifest does not exist: {source_manifest}")
    if output.exists():
        raise SystemExit(f"refusing to reuse existing output directory: {output}")
    all_sources = read_sources(source_manifest)
    unverified = [source.accession for source in all_sources if source.verification_status != "verified"]
    if unverified and not verified_only:
        raise SystemExit(
            "source manifest contains records requiring live verification: "
            + ", ".join(unverified)
            + "; update their lengths/checksums or pass --verified-only"
        )
    sources = [source for source in all_sources if source.verification_status == "verified"]
    sequences = read_fasta(sources_fasta)
    expected = {source.accession for source in sources}
    actual = set(sequences)
    if actual != expected:
        raise SystemExit(f"source FASTA accession mismatch; missing={sorted(expected - actual)}, extra={sorted(actual - expected)}")
    for source in sources:
        observed = sequences[source.accession]
        digest = sha256_bytes(observed.encode("ascii"))
        if len(observed) != source.length or digest != source.sequence_sha256:
            raise SystemExit(f"checksum or length verification failed for {source.accession}")

    backgrounds = [source for source in sources if source.query_kind in {"chromosome", "metagenome_assembly"} and not source.catalog]
    if not backgrounds:
        raise SystemExit("verified source manifest has no real chromosome or metagenome background records")
    background_group_names = {source.background_group for source in backgrounds}
    if len(background_group_names) != 1 or None in background_group_names:
        raise SystemExit("all selected backgrounds must belong to one declared background_group")
    background_group = next(iter(background_group_names))
    queries = [source for source in sources if source.query_kind in {"plasmid", "phage"}]
    if not queries:
        raise SystemExit("verified source manifest has no plasmid or phage query records")

    output.mkdir(parents=True, exist_ok=False)
    assemblies_dir = output / "assemblies"
    assemblies_dir.mkdir()
    queries_dir = output / "queries"
    queries_dir.mkdir()
    query_files: list[dict[str, object]] = []
    for source in queries:
        query_path = queries_dir / f"{source.accession}.fasta"
        write_fasta(query_path, [(source.accession, sequences[source.accession])])
        query_files.append(
            {
                "accession": source.accession,
                "query_kind": source.query_kind,
                "topology": source.topology,
                "file": str(query_path.relative_to(output)),
                "sha256": sha256_file(query_path),
            }
        )
    catalog_sources = [source for source in sources if source.catalog]
    write_fasta(output / "catalog.fasta", [(source.accession, sequences[source.accession]) for source in catalog_sources])
    background_records = [(f"background_{source.accession}", sequences[source.accession]) for source in backgrounds]
    background_accessions = [source.accession for source in backgrounds]

    verification_rows: list[dict[str, object]] = []
    selected_ids = {source.accession for source in sources}
    for source in all_sources:
        selected = source.accession in selected_ids
        verification_rows.append(
            {
                **{field: getattr(source, field) if field not in {"catalog", "background_group"} else source.catalog if field == "catalog" else source.background_group or "NA" for field in SOURCE_FIELDS},
                "observed_length": len(sequences[source.accession]) if selected else "",
                "observed_sha256": sha256_bytes(sequences[source.accession].encode("ascii")) if selected else "",
                "verified": "true" if selected else "false",
            }
        )
    _write_tsv(output / "source_verification.tsv", VERIFICATION_FIELDS, verification_rows)

    truth_rows: list[dict[str, object]] = []
    status_rows: list[dict[str, object]] = []
    assembly_metadata: list[dict[str, object]] = []
    for source in queries:
        fragments = fragments_for(source, sequences[source.accession])
        for case_type in DATASET_IDS:
            fragment = fragments[case_type]
            case_id = f"{source.accession}__{case_type}"
            if fragment is None:
                reason = "case requires declared source topology or coordinate annotation"
                status_rows.append(
                    {
                        "case_id": case_id,
                        "query_accession": source.accession,
                        "query_kind": source.query_kind,
                        "topology": source.topology,
                        "case_type": case_type,
                        "status": "unsupported",
                        "reason": reason,
                        "source_sha256": source.sequence_sha256,
                        "background_group": background_group,
                    }
                )
                continue
            assembly_records = list(background_records) + [(fragment.identifier, fragment.sequence)]
            assembly_path = assemblies_dir / f"{case_id}.fasta"
            write_fasta(assembly_path, assembly_records)
            assembly_sha = sha256_file(assembly_path)
            query_intervals_json = json.dumps(fragment.query_intervals, separators=(",", ":"))
            status_rows.append(
                {
                    "case_id": case_id,
                    "query_accession": source.accession,
                    "query_kind": source.query_kind,
                    "topology": source.topology,
                    "case_type": case_type,
                    "status": "materialized",
                    "reason": "",
                    "source_sha256": source.sequence_sha256,
                    "background_group": background_group,
                }
            )
            truth_rows.append(
                {
                    "case_id": case_id,
                    "query_accession": source.accession,
                    "query_kind": source.query_kind,
                    "topology": source.topology,
                    "expected_reference": source.accession,
                    "truth_class": f"controlled_{case_type}",
                    "construction": "declared source sequence transformation; not a biological presence call",
                    "background_group": background_group,
                    "background_accessions_json": json.dumps(background_accessions, separators=(",", ":")),
                    "source_length": len(sequences[source.accession]),
                    "source_sha256": source.sequence_sha256,
                    "source_start": fragment.source_start if fragment.source_start is not None else "multiple",
                    "source_end": fragment.source_end if fragment.source_end is not None else "multiple",
                    "query_intervals_json": query_intervals_json,
                    "contig_id": fragment.identifier,
                    "contig_start": 0,
                    "contig_end": len(fragment.sequence),
                    "orientation": fragment.orientation,
                    "origin_crossing": str(fragment.origin_crossing).lower(),
                    "terminal_repeat": str(fragment.terminal_repeat).lower(),
                    "shared_group": fragment.shared_group or "NA",
                    "coverage_bases": len(set(position for begin, end in fragment.query_intervals for position in range(begin, end))),
                    "coverage_fraction": f"{len(set(position for begin, end in fragment.query_intervals for position in range(begin, end))) / len(sequences[source.accession]):.6f}",
                    "assembly_path": str(assembly_path.relative_to(output)),
                    "assembly_sha256": assembly_sha,
                }
            )
            assembly_metadata.append(
                {
                    "case_id": case_id,
                    "file": str(assembly_path.relative_to(output)),
                    "sha256": assembly_sha,
                    "background_group": background_group,
                    "background_accessions": background_accessions,
                    "controlled_derivative": True,
                }
            )

    _write_tsv(output / "truth.tsv", TRUTH_FIELDS, truth_rows)
    _write_tsv(output / "status.tsv", STATUS_FIELDS, status_rows)
    metadata = {
        "schema_version": SCHEMA_VERSION,
        "generator": "mge_accession_spike.py",
        "source_manifest": source_manifest.name,
        "source_manifest_sha256": sha256_file(source_manifest),
        "source_fasta": sources_fasta.name,
        "source_fasta_sha256": sha256_file(sources_fasta),
        "selected_source_accessions": sorted(selected_ids),
        "omitted_unverified_accessions": sorted(set(unverified) - selected_ids),
        "background_group": background_group,
        "background_accessions": background_accessions,
        "query_accessions": [source.accession for source in queries],
        "queries": query_files,
        "catalog_accessions": [source.accession for source in catalog_sources],
        "source_records": [
            {
                "accession": source.accession,
                "display_name": source.display_name,
                "query_kind": source.query_kind,
                "topology": source.topology,
                "role": source.role,
                "catalog": source.catalog,
                "background_group": source.background_group,
                "length": source.length,
                "sequence_sha256": source.sequence_sha256,
                "record_url": source.record_url,
                "verification_status": source.verification_status,
                "terminal_repeat": [source.terminal_repeat_start, source.terminal_repeat_end],
                "shared_region": [source.shared_region_start, source.shared_region_end],
                "shared_group": source.shared_group,
            }
            for source in all_sources
        ],
        "assemblies": assembly_metadata,
        "truth_rows": len(truth_rows),
        "status_rows": len(status_rows),
        "limitations": [
            "Records marked needs_live_verification are omitted only with explicit --verified-only and are not negative evidence.",
            "Background records are supplied real chromosome or metagenome sequences reused across controlled derivatives; derivatives are not independent natural metagenomes.",
            "Topology, terminal-repeat coordinates, and shared-region coordinates are used only when declared in the source manifest.",
            "The forward, reverse, origin-crossing, terminal-repeat, and shared-region labels describe sequence constructions, not biological origin or autonomous plasmid status.",
            "This output does not establish sensitivity, specificity, metagenome absence, or confirmed plasmid/phage presence.",
        ],
    }
    (output / "dataset.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources-fasta", type=Path, required=True)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--verified-only",
        action="store_true",
        help="materialize only records marked verified; omitted records are listed in dataset.json",
    )
    args = parser.parse_args()
    try:
        metadata = generate(args.sources_fasta, args.source_manifest, args.output, args.verified_only)
    except (OSError, ValueError) as error:
        print(f"MGE accession generation failed: {error}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "output": str(args.output.resolve()),
                "schema_version": metadata["schema_version"],
                "truth_rows": metadata["truth_rows"],
                "status_rows": metadata["status_rows"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
