#!/usr/bin/env python3
"""Prepare deterministic trace benchmark inputs.

The default input is intentionally synthetic.  ``--dataset-root`` imports an
already materialized dataset (for example an accession-backed spike-in set)
without downloading or modifying it.  Only the selected query record and the
assembly files are copied into the benchmark work directory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import shutil
from pathlib import Path


PROFILES = {
    "smoke": {"query_length": 4_000, "background_length": 4_000},
    "standard": {"query_length": 10_000, "background_length": 8_000},
    "large": {"query_length": 25_000, "background_length": 12_000},
}


def random_dna(rng: random.Random, length: int) -> str:
    return "".join(rng.choices("ACGT", k=length))


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="ascii") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset : offset + 80] + "\n")


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    sequence: list[str] = []
    with path.open(encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(sequence).upper()))
                name = line[1:].split()[0]
                sequence = []
            elif name is None:
                raise ValueError(f"sequence data precedes a FASTA header in {path}")
            else:
                sequence.append(line)
    if name is not None:
        records.append((name, "".join(sequence).upper()))
    return records


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_truth(path: Path, rows: list[dict[str, str]]) -> None:
    fields = ["metagenome_id", "truth_class", "expected_candidate", "short_fragment"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def synthetic(output: Path, profile: str) -> str:
    settings = PROFILES[profile]
    rng = random.Random(0x4A414D5452414345)
    query = random_dna(rng, settings["query_length"])
    query_length = len(query)
    repeat = query[settings["query_length"] // 3 : settings["query_length"] // 3 + 600]
    assemblies = output / "assemblies"
    assemblies.mkdir(parents=True, exist_ok=False)
    records: list[tuple[str, list[tuple[str, str]], str, bool, bool]] = [
        ("sample_exact.fasta", [("exact_contig", query)], "exact", True, False),
        (
            "sample_reverse.fasta",
            [("reverse_contig", reverse_complement(query))],
            "reverse_complement",
            True,
            False,
        ),
        (
            "sample_partial.fasta",
            [("partial_contig", query[query_length // 4 : 3 * query_length // 4])],
            "partial",
            True,
            True,
        ),
        (
            "sample_split.fasta",
            [
                ("split_1", query[: query_length // 3 + 100]),
                ("split_2", query[query_length // 3 - 100 : 2 * query_length // 3 + 100]),
                ("split_3", query[2 * query_length // 3 - 100 :]),
            ],
            "split_overlap",
            True,
            True,
        ),
        (
            "sample_ambiguous.fasta",
            [("mobile_element_only", repeat)],
            "ambiguous_repeat",
            False,
            True,
        ),
        (
            "sample_background.fasta",
            [("chromosome_like", random_dna(rng, settings["background_length"]))],
            "background",
            False,
            False,
        ),
    ]
    truth: list[dict[str, str]] = []
    for filename, contigs, truth_class, expected, short_fragment in records:
        path = assemblies / filename
        write_fasta(path, contigs)
        truth.append(
            {
                "metagenome_id": filename,
                "truth_class": truth_class,
                "expected_candidate": "true" if expected else "false",
                "short_fragment": "true" if short_fragment else "false",
            }
        )

    query_path = output / "query.fasta"
    write_fasta(query_path, [("synthetic_query", query)])
    write_truth(output / "truth.tsv", truth)
    metadata = {
        "schema_version": "1.0.0",
        "dataset_kind": "synthetic",
        "profile": profile,
        "seed": "0x4A414D5452414345",
        "query_id": "synthetic_query",
        "query_length": len(query),
        "assembly_count": len(records),
        "contig_count": sum(len(contigs) for _, contigs, *_ in records),
        "assembly_bases": sum(len(sequence) for _, contigs, *_ in records for _, sequence in contigs),
        "query_sha256": sha256(query_path),
    }
    (output / "dataset.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return "synthetic_query"


def imported(output: Path, dataset_root: Path, query_id: str | None) -> str:
    catalog_path = dataset_root / "catalog.fasta"
    assembly_root = dataset_root / "assemblies"
    if not catalog_path.is_file():
        raise SystemExit(f"dataset root has no catalog.fasta: {dataset_root}")
    if not assembly_root.is_dir():
        raise SystemExit(f"dataset root has no assemblies directory: {dataset_root}")
    catalog = read_fasta(catalog_path)
    if not catalog:
        raise SystemExit(f"catalog is empty: {catalog_path}")
    selected_id = query_id or catalog[0][0]
    selected = next((sequence for name, sequence in catalog if name == selected_id), None)
    if selected is None:
        raise SystemExit(f"query id {selected_id!r} is absent from {catalog_path}")

    assemblies = output / "assemblies"
    assemblies.mkdir(parents=True, exist_ok=False)
    source_files = sorted(
        path
        for path in assembly_root.iterdir()
        if path.is_file() and path.suffix.lower() in {".fa", ".fasta", ".fna", ".fastq", ".fq"}
    )
    if not source_files:
        raise SystemExit(f"dataset has no FASTA/FASTQ assemblies: {assembly_root}")
    for source in source_files:
        shutil.copy2(source, assemblies / source.name)

    query_path = output / "query.fasta"
    write_fasta(query_path, [(selected_id, selected)])
    truth_rows: list[dict[str, str]] = []
    source_truth = dataset_root / "truth.tsv"
    by_dataset: dict[str, list[dict[str, str]]] = {}
    if source_truth.is_file():
        with source_truth.open(encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                by_dataset.setdefault(row.get("dataset", ""), []).append(row)
    for source in source_files:
        sample_id = source.name
        sample_key = source.stem
        matching = [
            row
            for key in (sample_id, sample_key)
            for row in by_dataset.get(key, [])
            if row.get("expected_reference") == selected_id
        ]
        if matching:
            truth_class = matching[0].get("truth_class", "external_truth")
            expected = "true"
            short_fragment = "true" if truth_class in {"partial", "split_overlap", "reverse_complement"} else "false"
        else:
            truth_class = "background_or_not_in_truth"
            expected = "false"
            short_fragment = "false"
        truth_rows.append(
            {
                "metagenome_id": sample_id,
                "truth_class": truth_class,
                "expected_candidate": expected,
                "short_fragment": short_fragment,
            }
        )
    write_truth(output / "truth.tsv", truth_rows)
    metadata = {
        "schema_version": "1.0.0",
        "dataset_kind": "imported",
        "source_root": str(dataset_root.resolve()),
        "source_dataset_sha256": sha256(dataset_root / "dataset.json") if (dataset_root / "dataset.json").is_file() else None,
        "catalog_path": str(catalog_path.resolve()),
        "catalog_sha256": sha256(catalog_path),
        "query_id": selected_id,
        "query_length": len(selected),
        "assembly_count": len(source_files),
        "source_truth_path": str(source_truth.resolve()) if source_truth.is_file() else None,
    }
    (output / "dataset.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return selected_id


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--profile", choices=sorted(PROFILES), default="smoke")
    parser.add_argument("--dataset-root", type=Path)
    parser.add_argument("--query-id")
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"refusing to reuse existing dataset work directory: {args.output}")
    args.output.mkdir(parents=True)
    if args.dataset_root:
        query_id = imported(args.output, args.dataset_root.resolve(), args.query_id)
    else:
        query_id = synthetic(args.output, args.profile)
    print(json.dumps({"query_id": query_id, "dataset": str(args.output.resolve())}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
