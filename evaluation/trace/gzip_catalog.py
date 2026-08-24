#!/usr/bin/env python3
"""Materialize deterministic ordinary-gzip sources for trace benchmarks."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import shutil
from pathlib import Path


def checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def compress(source: Path, output: Path) -> None:
    with source.open("rb") as reader, output.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=1, mtime=0) as writer:
            shutil.copyfileobj(reader, writer, length=1024 * 1024)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"refusing to reuse output directory: {args.output}")
    assemblies = args.output / "assemblies"
    assemblies.mkdir(parents=True)
    with args.catalog.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    output_rows = []
    for index, row in enumerate(rows):
        source = Path(row["resource_uri"])
        output = assemblies / f"source-{index:06d}.fasta.gz"
        compress(source, output)
        output_rows.append(
            {
                "metagenome_id": row["metagenome_id"],
                "resource_uri": str(output.resolve()),
                "sha256": checksum(output),
            }
        )
    with (args.output / "sources.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("metagenome_id", "resource_uri", "sha256"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(output_rows)
    print(args.output / "sources.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
