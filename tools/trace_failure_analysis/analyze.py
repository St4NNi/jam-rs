#!/usr/bin/env python3
"""Attribute missing truth bases and summarize stage-level trace counters."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.trace_failure_analysis.diagnostics import DiagnosticError, load_reports, summarize_reports


def _write_new(path: Path, value: object, *, pretty: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    serialized = json.dumps(value, indent=2 if pretty else None, sort_keys=True) + "\n"
    try:
        with path.open("x", encoding="utf-8") as handle:
            handle.write(serialized)
    except FileExistsError as exc:
        raise DiagnosticError(f"refusing to overwrite existing output: {path}") from exc


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="TraceDiagnosticReport JSON or JSONL")
    parser.add_argument("--output", required=True, type=Path, help="new attribution summary JSON")
    parser.add_argument("--pretty", action="store_true", help="indent the JSON output")
    args = parser.parse_args()
    try:
        result = summarize_reports(load_reports(args.input))
        _write_new(args.output, result, pretty=args.pretty)
    except (DiagnosticError, OSError) as exc:
        parser.exit(1, f"diagnostic analysis failed: {exc}\n")
    print(json.dumps({"status": "pass", "reports": result["report_count"], "output": args.output.as_posix()}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
