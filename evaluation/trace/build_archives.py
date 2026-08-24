#!/usr/bin/env python3
"""Build one JMA archive per assembly for a trace benchmark.

This helper is intentionally separate from timing.  The benchmark wrapper
times this complete archive-construction stage through ``measure.py`` and
does not include compiler or environment setup.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jam", type=Path, required=True)
    parser.add_argument("--assemblies-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--primary-scale", type=int, default=100)
    parser.add_argument("--rescue-scale", type=int, default=100)
    parser.add_argument("--sequence-block-policy", choices=("fixed", "gear"), default="fixed")
    parser.add_argument("--sequence-block-codec", choices=("raw2bit", "zstd2bit"), default="raw2bit")
    parser.add_argument("--block-bases", type=int, default=1024 * 1024)
    parser.add_argument("--gear-min-bases", type=int, default=16 * 1024)
    parser.add_argument("--gear-target-bases", type=int, default=64 * 1024)
    parser.add_argument("--gear-max-bases", type=int, default=256 * 1024)
    parser.add_argument(
        "--gear-table",
        choices=("single-base", "dinucleotide", "packed-four-base"),
        default="single-base",
    )
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--memory-target", type=int, default=4)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise SystemExit(f"refusing to reuse archive directory: {args.output_dir}")
    args.output_dir.mkdir(parents=True)
    assemblies = sorted(
        path
        for path in args.assemblies_dir.iterdir()
        if path.is_file() and path.suffix.lower() in {".fa", ".fasta", ".fna", ".fastq", ".fq"}
    )
    if not assemblies:
        raise SystemExit(f"no FASTA/FASTQ assemblies found in {args.assemblies_dir}")
    built: list[dict[str, str | int]] = []
    for assembly in assemblies:
        output = args.output_dir / f"{assembly.name}.jma"
        command = [
            str(args.jam),
            "--threads",
            str(args.threads),
            "--memory-target",
            str(args.memory_target),
            "--silent",
            "archive",
            "--input",
            str(assembly),
            "--output",
            str(output),
            "--primary-scale",
            str(args.primary_scale),
            "--rescue-scale",
            str(args.rescue_scale),
            "--sequence-block-policy",
            args.sequence_block_policy,
            "--sequence-block-codec",
            args.sequence_block_codec,
            "--block-bases",
            str(args.block_bases),
            "--gear-min-bases",
            str(args.gear_min_bases),
            "--gear-target-bases",
            str(args.gear_target_bases),
            "--gear-max-bases",
            str(args.gear_max_bases),
            "--gear-table",
            args.gear_table,
        ]
        subprocess.run(command, check=True)
        if not output.is_file():
            raise SystemExit(f"archive command did not create JMA object: {output}")
        built.append(
            {
                "metagenome_id": assembly.name,
                "source_assembly": str(assembly.resolve()),
                "source_assembly_sha256": sha256(assembly),
                "resource_uri": str(output.resolve()),
                "archive_bytes": output.stat().st_size,
                "sha256": sha256(output),
            }
        )
    (args.output_dir / "build_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": "1.0.0",
                "primary_k": 31,
                "rescue_k": 21,
                "primary_scale": args.primary_scale,
                "rescue_scale": args.rescue_scale,
                "sequence_block_policy": args.sequence_block_policy,
                "sequence_block_codec": args.sequence_block_codec,
                "archives": built,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
