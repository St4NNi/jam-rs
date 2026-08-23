#!/usr/bin/env python3
"""Write a trace resource catalog for local or mock-HTTP JMA resources."""

from __future__ import annotations

import argparse
from pathlib import Path


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
    rows: list[tuple[str, str, str, str]] = []
    assemblies = sorted(
        path
        for path in args.assemblies_dir.iterdir()
        if path.is_file() and path.suffix.lower() in {".fa", ".fasta", ".fna", ".fastq", ".fq"}
    )
    for assembly in assemblies:
        archive_name = f"{assembly.name}.jma"
        archive = args.jma_dir / archive_name
        index_name = f"{archive_name}.idx.json"
        index = args.jma_dir / index_name
        if not archive.is_file():
            raise SystemExit(f"missing JMA archive for {assembly.name}: {archive}")
        if not index.is_file():
            raise SystemExit(f"missing JMA range index for {assembly.name}: {index}")
        if args.mode == "local":
            resource = str(archive.relative_to(args.output.parent))
            index_resource = str(index.relative_to(args.output.parent))
        else:
            resource = f"{args.base_url.rstrip('/')}/{archive_name}"
            index_resource = f"{args.base_url.rstrip('/')}/{index_name}"
        rows.append(
            (
                assembly.name,
                resource,
                index_resource,
                str(assembly.relative_to(args.output.parent)),
            )
        )
    if not rows:
        raise SystemExit(f"no assemblies found in {args.assemblies_dir}")
    with args.output.open("w", encoding="utf-8") as handle:
        handle.write("metagenome_id\tjma\tjma_index\traw\n")
        for metagenome_id, jma, jma_index, raw in rows:
            handle.write(f"{metagenome_id}\t{jma}\t{jma_index}\t{raw}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
