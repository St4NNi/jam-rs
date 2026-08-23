#!/usr/bin/env python3
"""Generate deterministic accession-pinned plasmid spike-in assemblies.

The input FASTA is expected to contain the accession-versioned records listed
in ``manifests/accession-sources-v1.tsv``.  The importer verifies every record
before writing anything.  It then fragments the three supplied chromosome
records with an accession-derived seed and appends exact plasmid-derived
contigs to four controlled assemblies.  No sequences are downloaded, invented,
or treated as independent natural metagenomes by this script.

Only generation is supported.  Existing output directories are refused so a
previous materialization cannot be silently overwritten.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import random
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO


SCHEMA_VERSION = "jam-trace-accession-spike-v1"
DATASET_IDS = ("exact_full", "fragmented_trace", "near_known", "background_only")
SOURCE_FIELDS = (
    "accession",
    "name",
    "role",
    "catalog",
    "length",
    "sequence_sha256",
    "record_url",
)
TRUTH_FIELDS = (
    "dataset",
    "spike_accession",
    "expected_reference",
    "truth_class",
    "construction",
    "source_length",
    "source_query_sha256",
    "source_start",
    "source_end",
    "orientation",
    "coverage_bases",
    "coverage_percent",
    "overlap_bases",
    "coordinates_json",
    "injected_contigs",
    "injected_bases",
)
VERIFICATION_FIELDS = (
    "accession",
    "name",
    "role",
    "catalog",
    "length",
    "sequence_sha256",
    "record_url",
    "verified",
)


@dataclass(frozen=True)
class Source:
    accession: str
    name: str
    role: str
    catalog: bool
    length: int
    sequence_sha256: str
    record_url: str


@dataclass(frozen=True)
class InjectedContig:
    identifier: str
    sequence: str
    source_start: int
    source_end: int
    orientation: str
    source_positions: tuple[int, ...]


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
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
                    records[identifier] = _finish_sequence(path, identifier, sequence)
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
        records[identifier] = _finish_sequence(path, identifier, sequence)
    if not records:
        raise ValueError(f"source FASTA is empty: {path}")
    return records


def _finish_sequence(path: Path, identifier: str, sequence: list[str]) -> str:
    value = "".join(sequence).upper()
    if not value:
        raise ValueError(f"empty sequence {identifier!r} in {path}")
    invalid = sorted(set(value) - set("ACGTN"))
    if invalid:
        raise ValueError(f"unsupported bases {invalid!r} in {path}:{identifier}")
    return value


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset : offset + 80] + "\n")


def read_sources(path: Path) -> list[Source]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != SOURCE_FIELDS:
            raise ValueError(f"unexpected source manifest columns in {path}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"source manifest is empty: {path}")
    sources: list[Source] = []
    seen: set[str] = set()
    for row in rows:
        if any(not row.get(field) for field in SOURCE_FIELDS):
            raise ValueError(f"source manifest has an empty field: {row}")
        accession = row["accession"]
        if accession in seen:
            raise ValueError(f"duplicate accession in source manifest: {accession}")
        seen.add(accession)
        digest = row["sequence_sha256"].lower()
        if len(digest) != 64 or any(character not in "0123456789abcdef" for character in digest):
            raise ValueError(f"invalid sequence SHA-256 for {accession}")
        sources.append(
            Source(
                accession=accession,
                name=row["name"],
                role=row["role"],
                catalog=row["catalog"].lower() == "true",
                length=int(row["length"]),
                sequence_sha256=digest,
                record_url=row["record_url"],
            )
        )
    return sources


def source_seed(accession: str) -> int:
    return int.from_bytes(hashlib.sha256(accession.encode("ascii")).digest()[:8], "big")


def fragment_background(accession: str, sequence: str) -> list[tuple[str, str]]:
    """Split a real chromosome record without changing its sequence."""
    rng = random.Random(source_seed(accession))
    records: list[tuple[str, str]] = []
    offset = 0
    index = 1
    while offset < len(sequence):
        end = min(len(sequence), offset + rng.randint(8_000, 20_000))
        if len(sequence) - end < 1_000:
            end = len(sequence)
        records.append((f"background_{accession}_{index:04d}", sequence[offset:end]))
        offset = end
        index += 1
    return records


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def split_with_overlap(identifier: str, sequence: str, parts: int, overlap: int) -> list[InjectedContig]:
    records: list[InjectedContig] = []
    for index in range(parts):
        start = max(0, index * len(sequence) // parts - (overlap if index else 0))
        end = min(len(sequence), (index + 1) * len(sequence) // parts + (overlap if index + 1 < parts else 0))
        fragment = sequence[start:end]
        records.append(
            InjectedContig(
                identifier=f"{identifier}_{index + 1}",
                sequence=fragment,
                source_start=start,
                source_end=end,
                orientation="forward",
                source_positions=tuple(range(start, end)),
            )
        )
    return records


def whole_contig(identifier: str, sequence: str, orientation: str = "forward") -> InjectedContig:
    if orientation == "reverse":
        observed = reverse_complement(sequence)
        positions = tuple(reversed(range(len(sequence))))
    else:
        observed = sequence
        positions = tuple(range(len(sequence)))
    return InjectedContig(
        identifier=identifier,
        sequence=observed,
        source_start=0,
        source_end=len(sequence),
        orientation=orientation,
        source_positions=positions,
    )


def interval_json(contigs: list[InjectedContig]) -> str:
    return json.dumps(
        [
            {
                "contig_id": contig.identifier,
                "contig_start": 0,
                "contig_end": len(contig.sequence),
                "source_start": contig.source_start,
                "source_end": contig.source_end,
                "orientation": contig.orientation,
                "source_positions_sha256": sha256_bytes(
                    ",".join(str(position) for position in contig.source_positions).encode("ascii")
                ),
            }
            for contig in contigs
        ],
        sort_keys=True,
        separators=(",", ":"),
    )


def _write_tsv(path: Path, fields: tuple[str, ...], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def make_truth_row(
    dataset: str,
    spike_accession: str,
    expected_reference: str,
    truth_class: str,
    construction: str,
    source_sequence: str,
    contigs: list[InjectedContig],
    source_start: int | str,
    source_end: int | str,
    coverage_bases: int,
    coverage_percent: float,
    overlap_bases: int,
) -> dict[str, object]:
    return {
        "dataset": dataset,
        "spike_accession": spike_accession,
        "expected_reference": expected_reference,
        "truth_class": truth_class,
        "construction": construction,
        "source_length": len(source_sequence),
        "source_query_sha256": sha256_bytes(source_sequence.encode("ascii")),
        "source_start": source_start,
        "source_end": source_end,
        "orientation": contigs[0].orientation if len({item.orientation for item in contigs}) == 1 else "mixed",
        "coverage_bases": coverage_bases,
        "coverage_percent": f"{coverage_percent:.6f}",
        "overlap_bases": overlap_bases,
        "coordinates_json": interval_json(contigs),
        "injected_contigs": len(contigs),
        "injected_bases": sum(len(item.sequence) for item in contigs),
    }


def generate(sources_fasta: Path, source_manifest: Path, output: Path) -> dict[str, object]:
    sources_fasta = sources_fasta.resolve()
    source_manifest = source_manifest.resolve()
    output = output.resolve()
    if not sources_fasta.is_file():
        raise SystemExit(f"source FASTA does not exist: {sources_fasta}")
    if not source_manifest.is_file():
        raise SystemExit(f"source manifest does not exist: {source_manifest}")
    if output.exists():
        raise SystemExit(f"refusing to reuse existing output directory: {output}")

    sources = read_sources(source_manifest)
    sequences = read_fasta(sources_fasta)
    expected_accessions = {source.accession for source in sources}
    actual_accessions = set(sequences)
    if actual_accessions != expected_accessions:
        missing = sorted(expected_accessions - actual_accessions)
        extra = sorted(actual_accessions - expected_accessions)
        raise SystemExit(f"source FASTA accession mismatch; missing={missing}, extra={extra}")
    for source in sources:
        sequence = sequences[source.accession]
        digest = sha256_bytes(sequence.encode("ascii"))
        if len(sequence) != source.length or digest != source.sequence_sha256:
            raise SystemExit(f"source sequence verification failed for {source.accession}")

    output.mkdir(parents=True, exist_ok=False)
    assemblies_dir = output / "assemblies"
    assemblies_dir.mkdir()
    verification_rows = [
        {
            "accession": source.accession,
            "name": source.name,
            "role": source.role,
            "catalog": str(source.catalog).lower(),
            "length": source.length,
            "sequence_sha256": source.sequence_sha256,
            "record_url": source.record_url,
            "verified": "true",
        }
        for source in sources
    ]
    _write_tsv(output / "source_verification.tsv", VERIFICATION_FIELDS, verification_rows)

    catalog_accessions = [source.accession for source in sources if source.catalog]
    catalog_path = output / "catalog.fasta"
    write_fasta(catalog_path, [(accession, sequences[accession]) for accession in catalog_accessions])

    background_accessions = [source.accession for source in sources if source.role == "chromosome background"]
    if len(background_accessions) != 3:
        raise SystemExit(f"expected exactly three chromosome backgrounds, got {len(background_accessions)}")
    background_records = [
        record
        for accession in background_accessions
        for record in fragment_background(accession, sequences[accession])
    ]
    assemblies: dict[str, list[tuple[str, str]]] = {dataset_id: list(background_records) for dataset_id in DATASET_IDS}
    truth_rows: list[dict[str, object]] = []

    def add_truth(
        dataset: str,
        spike_accession: str,
        expected_reference: str,
        truth_class: str,
        construction: str,
        contigs: list[InjectedContig],
        source_start: int | str,
        source_end: int | str,
        coverage_bases: int,
        overlap_bases: int = 0,
    ) -> None:
        assemblies[dataset].extend((contig.identifier, contig.sequence) for contig in contigs)
        source_sequence = sequences[spike_accession]
        truth_rows.append(
            make_truth_row(
                dataset,
                spike_accession,
                expected_reference,
                truth_class,
                construction,
                source_sequence,
                contigs,
                source_start,
                source_end,
                coverage_bases,
                coverage_bases * 100 / len(source_sequence),
                overlap_bases,
            )
        )

    pbr322 = sequences["J01749.1"]
    add_truth(
        "exact_full",
        "J01749.1",
        "J01749.1",
        "exact",
        "complete accession",
        [whole_contig("spike_pBR322", pbr322)],
        0,
        len(pbr322),
        len(pbr322),
    )
    r100 = sequences["AP000342.1"]
    add_truth(
        "exact_full",
        "AP000342.1",
        "AP000342.1",
        "exact",
        "complete accession",
        [whole_contig("spike_R100", r100)],
        0,
        len(r100),
        len(r100),
    )
    pkjk5 = sequences["NC_008272.1"]
    add_truth(
        "exact_full",
        "NC_008272.1",
        "NC_008272.1",
        "exact",
        "complete accession",
        [whole_contig("spike_pKJK5", pkjk5)],
        0,
        len(pkjk5),
        len(pkjk5),
    )
    add_truth(
        "fragmented_trace",
        "J01749.1",
        "J01749.1",
        "reverse_complement",
        "complete reverse complement",
        [whole_contig("spike_pBR322_reverse", pbr322, "reverse")],
        0,
        len(pbr322),
        len(pbr322),
    )
    r100_fragments = split_with_overlap("spike_R100_split", r100, 5, 1_000)
    add_truth(
        "fragmented_trace",
        "AP000342.1",
        "AP000342.1",
        "split_overlap",
        "five fragments with 1000 bp boundary overlap",
        r100_fragments,
        "multiple",
        "multiple",
        len(r100),
        4 * 1_000,
    )
    partial_start = 2 * len(pkjk5) // 5
    partial_end = 3 * len(pkjk5) // 5
    partial = InjectedContig(
        identifier="spike_pKJK5_partial",
        sequence=pkjk5[partial_start:partial_end],
        source_start=partial_start,
        source_end=partial_end,
        orientation="forward",
        source_positions=tuple(range(partial_start, partial_end)),
    )
    add_truth(
        "fragmented_trace",
        "NC_008272.1",
        "NC_008272.1",
        "partial",
        "central 20 percent",
        [partial],
        partial_start,
        partial_end,
        partial_end - partial_start,
    )
    puc19 = sequences["M77789.2"]
    add_truth(
        "near_known",
        "M77789.2",
        "J01749.1",
        "near_known_accession",
        "complete pUC19 screened against catalog pBR322",
        [whole_contig("spike_pUC19", puc19)],
        0,
        len(puc19),
        len(puc19),
    )
    r1 = sequences["KY749247.1"]
    r1_fragments = split_with_overlap("spike_R1_split", r1, 6, 1_000)
    add_truth(
        "near_known",
        "KY749247.1",
        "AP000342.1",
        "near_known_accession",
        "six R1 fragments with 1000 bp boundary overlap screened against R100",
        r1_fragments,
        "multiple",
        "multiple",
        len(r1),
        5 * 1_000,
    )

    assembly_metadata: list[dict[str, object]] = []
    for dataset_id in DATASET_IDS:
        path = assemblies_dir / f"{dataset_id}.fasta"
        write_fasta(path, assemblies[dataset_id])
        assembly_metadata.append(
            {
                "dataset": dataset_id,
                "file": str(path.relative_to(output)),
                "contigs": len(assemblies[dataset_id]),
                "bases": sum(len(sequence) for _, sequence in assemblies[dataset_id]),
                "sha256": sha256_file(path),
                "background_source_accessions": background_accessions,
                "controlled_derivative": True,
            }
        )
    _write_tsv(output / "truth.tsv", TRUTH_FIELDS, truth_rows)

    metadata = {
        "schema_version": SCHEMA_VERSION,
        "generator": "accession_spike.py",
        "source_manifest": source_manifest.name,
        "source_manifest_sha256": sha256_file(source_manifest),
        "source_fasta": sources_fasta.name,
        "source_fasta_sha256": sha256_file(sources_fasta),
        "source_endpoint": "NCBI nuccore accession-version records; URLs are recorded in the source manifest",
        "construction_seed": "sha256(accession) first 64 bits for chromosome fragmentation",
        "catalog_accessions": catalog_accessions,
        "catalog_bases": sum(len(sequences[accession]) for accession in catalog_accessions),
        "catalog_sha256": sha256_file(catalog_path),
        "background_accessions": background_accessions,
        "source_records": [
            {
                "accession": source.accession,
                "name": source.name,
                "role": source.role,
                "catalog": source.catalog,
                "length": source.length,
                "sequence_sha256": source.sequence_sha256,
                "record_url": source.record_url,
            }
            for source in sources
        ],
        "assemblies": assembly_metadata,
        "truth_targets": len(truth_rows),
        "limitations": [
            "Assembly-level accession spike-in; reads were not simulated and no assembler was run.",
            "Chromosome contigs are deterministic perfect fragments without sequencing or assembly errors.",
            "Every controlled assembly reuses the same three supplied chromosome sources; these are not independent natural metagenomes.",
            "Near-known rows are candidate-retrieval evidence against the declared catalog, not confirmed plasmid presence.",
            "This dataset does not establish biological sensitivity, specificity, or absence from a metagenome.",
        ],
    }
    (output / "dataset.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources-fasta", type=Path, required=True)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    try:
        metadata = generate(arguments.sources_fasta, arguments.source_manifest, arguments.output)
    except (OSError, ValueError) as error:
        print(f"accession spike generation failed: {error}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "output": str(arguments.output.resolve()),
                "schema_version": metadata["schema_version"],
                "truth_targets": metadata["truth_targets"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
