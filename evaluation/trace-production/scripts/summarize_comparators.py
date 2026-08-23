#!/usr/bin/env python3
"""Summarize normalized comparator results and their separate timing records."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "1.0.0"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_object(path: Path) -> dict[str, Any] | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return value if isinstance(value, dict) else None


def timing_path(result_path: Path) -> Path:
    direct = result_path.with_name(f"{result_path.stem}-measurement.json")
    if direct.is_file():
        return direct
    if result_path.stem == "lexicmap":
        alternate = result_path.with_name("lexicmap-search-measurement.json")
        if alternate.is_file():
            return alternate
    return direct


def observed_version(result: dict[str, Any]) -> str | None:
    for item in result.get("results", []):
        if isinstance(item, dict) and item.get("version"):
            return str(item["version"])
    metadata = result.get("tool_metadata")
    if isinstance(metadata, dict) and metadata.get("version"):
        return str(metadata["version"])
    return None


def rows(raw_dir: Path) -> tuple[list[dict[str, Any]], dict[str, str]]:
    output: list[dict[str, Any]] = []
    hashes: dict[str, str] = {}
    for path in sorted(raw_dir.rglob("*.json")):
        if path.name.endswith("-measurement.json"):
            continue
        result = read_object(path)
        if result is None or not result.get("tool") or not isinstance(result.get("metrics"), dict):
            continue
        timing_file = timing_path(path)
        timing = read_object(timing_file) if timing_file.is_file() else {}
        index_file = path.with_name("lexicmap-index-measurement.json")
        index = read_object(index_file) if result.get("tool") == "lexicmap" and index_file.is_file() else {}
        relative = path.relative_to(raw_dir)
        metrics = result["metrics"]
        selection = result.get("selection") if isinstance(result.get("selection"), dict) else {}
        output.append(
            {
                "schema_version": SCHEMA_VERSION,
                "case": relative.parent.as_posix(),
                "configuration": path.stem,
                "tool": result.get("tool"),
                "version": observed_version(result),
                "scope": result.get("scope"),
                "status": result.get("status"),
                "selected_assemblies": selection.get("selected_count"),
                "assemblies_total": metrics.get("assemblies_total"),
                "assemblies_with_records": metrics.get("assemblies_with_records"),
                "base_precision": metrics.get("base_precision"),
                "base_recall": metrics.get("base_recall"),
                "interval_precision": metrics.get("interval_precision"),
                "interval_recall": metrics.get("interval_recall"),
                "gap_error_bases": metrics.get("gap_error_bases"),
                "search_wall_seconds": timing.get("wall_seconds"),
                "search_cpu_seconds": timing.get("cpu_seconds"),
                "search_peak_rss_kib": timing.get("peak_rss_kib"),
                "index_wall_seconds": index.get("wall_seconds"),
                "index_cpu_seconds": index.get("cpu_seconds"),
                "index_peak_rss_kib": index.get("peak_rss_kib"),
            }
        )
        for used in (path, timing_file, index_file):
            if used.is_file():
                hashes[used.relative_to(raw_dir).as_posix()] = sha256(used)
    return output, dict(sorted(hashes.items()))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    raw_dir = args.raw_dir.resolve()
    output_dir = args.output_dir.resolve()
    if not raw_dir.is_dir():
        raise SystemExit(f"raw comparator directory does not exist: {raw_dir}")
    if output_dir.exists():
        raise SystemExit(f"refusing to reuse output directory: {output_dir}")
    output_dir.mkdir(parents=True)
    summary_rows, hashes = rows(raw_dir)
    fields = list(summary_rows[0]) if summary_rows else ["schema_version", "case", "configuration", "tool", "status"]
    with (output_dir / "summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(
            {field: "NA" if row.get(field) is None else row.get(field) for field in fields}
            for row in summary_rows
        )
    summary = {
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "raw_directory": str(raw_dir),
        "rows": summary_rows,
        "raw_files": hashes,
        "limits": [
            "Candidate rows use the exact assemblies selected by the corresponding Jam run.",
            "Index construction and search are recorded separately for LexicMap.",
            "Construction-truth precision is not biological specificity in the reused chromosome backgrounds.",
        ],
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
