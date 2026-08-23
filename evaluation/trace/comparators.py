#!/usr/bin/env python3
"""Optional candidate comparisons against locally installed BLASTn/minimap2."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path


def run_blastn(tool: str, query: Path, assemblies: list[Path], output: Path) -> dict:
    rows = []
    for assembly in assemblies:
        completed = subprocess.run(
            [
                tool,
                "-query",
                str(query),
                "-subject",
                str(assembly),
                "-task",
                "megablast",
                "-outfmt",
                "6 qseqid sseqid pident length qlen qstart qend sstart send bitscore",
                "-evalue",
                "1e-10",
                "-max_target_seqs",
                "5",
            ],
            check=True,
            stdout=subprocess.PIPE,
            text=True,
        )
        hits = [line for line in completed.stdout.splitlines() if line.strip()]
        rows.append({"metagenome_id": assembly.name, "hits": len(hits)})
    return {"tool": "blastn", "version": version(tool), "rows": rows}


def run_minimap2(tool: str, query: Path, assemblies: list[Path]) -> dict:
    rows = []
    for assembly in assemblies:
        completed = subprocess.run(
            [tool, "-x", "asm5", "--secondary=no", str(assembly), str(query)],
            check=True,
            stdout=subprocess.PIPE,
            text=True,
        )
        hits = [line for line in completed.stdout.splitlines() if line.strip() and not line.startswith("@")]
        rows.append({"metagenome_id": assembly.name, "hits": len(hits)})
    return {"tool": "minimap2", "version": version(tool), "rows": rows}


def version(tool: str) -> str:
    completed = subprocess.run([tool, "--version"], check=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return completed.stdout.strip().splitlines()[0] if completed.stdout.strip() else "unknown"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tool", choices=("blastn", "minimap2"), required=True)
    parser.add_argument("--query", type=Path, required=True)
    parser.add_argument("--assemblies-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    executable = shutil.which(args.tool)
    if executable is None:
        result = {"schema_version": "1.0.0", "tool": args.tool, "status": "unavailable"}
    else:
        assemblies = sorted(
            path
            for path in args.assemblies_dir.iterdir()
            if path.is_file() and path.suffix.lower() in {".fa", ".fasta", ".fna", ".fastq", ".fq"}
        )
        if args.tool == "blastn":
            result = run_blastn(executable, args.query, assemblies, args.output)
        else:
            result = run_minimap2(executable, args.query, assemblies)
        result["schema_version"] = "1.0.0"
        result["status"] = "measured"
        result["candidate_count"] = sum(row["hits"] > 0 for row in result["rows"])
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if result["status"] == "unavailable":
        return 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
