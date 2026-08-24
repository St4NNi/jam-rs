#!/usr/bin/env python3
"""Write a trace resource catalog for local or mock-HTTP JMA resources."""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assemblies-dir", type=Path, required=True)
    parser.add_argument("--jma-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mode", choices=("local", "http"), default="local")
    parser.add_argument("--base-url")
    args = parser.parse_args()
    if args.mode == "http" and not args.base_url:
        parser.error("--base-url is required in HTTP mode")
    rows: list[tuple[str, str, str]] = []
    assemblies = sorted(
        path
        for path in args.assemblies_dir.iterdir()
        if path.is_file() and path.suffix.lower() in {".fa", ".fasta", ".fna", ".fastq", ".fq"}
    )
    for assembly in assemblies:
        archive_name = f"{assembly.name}.jma"
        archive = args.jma_dir / archive_name
        if not archive.is_file():
            raise SystemExit(f"missing JMA archive for {assembly.name}: {archive}")
        if args.mode == "local":
            resource = os.path.relpath(archive, args.output.parent)
        else:
            resource = f"{args.base_url.rstrip('/')}/{archive_name}"
        rows.append((assembly.name, resource, sha256(archive)))
    if not rows:
        raise SystemExit(f"no assemblies found in {args.assemblies_dir}")
    with args.output.open("w", encoding="utf-8") as handle:
        handle.write("metagenome_id\tresource_uri\tsha256\n")
        for metagenome_id, resource_uri, checksum in rows:
            handle.write(f"{metagenome_id}\t{resource_uri}\t{checksum}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
