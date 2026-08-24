#!/usr/bin/env python3
"""Normalize BLASTn/minimap2/JAM alignment JSON through one interval path."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from common import ExperimentError, ensure_new_output, workspace_path
from tools.trace_failure_analysis.comparators import ComparatorError, normalize_comparator
from tools.trace_failure_analysis.intervals import IntervalError, parse_intervals


def _read(path: Path) -> object:
    text = path.read_text(encoding="utf-8")
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        values = []
        for line_number, line in enumerate(text.splitlines(), 1):
            if not line.strip():
                continue
            try:
                values.append(json.loads(line))
            except json.JSONDecodeError as exc:
                raise ComparatorError(f"invalid JSON at {path}:{line_number}") from exc
        return values


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--query-length", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--truth", type=Path, help="JSON interval list, optional")
    parser.add_argument("--query-deletion-operation", choices=("I", "D"), default="I")
    parser.add_argument("--pretty", action="store_true")
    args = parser.parse_args()
    try:
        input_path = workspace_path(args.input, field="input")
        output_path = ensure_new_output(args.output)
        truth = []
        if args.truth:
            truth_path = workspace_path(args.truth, field="truth")
            truth = parse_intervals(json.loads(truth_path.read_text(encoding="utf-8")), field="truth")
        result = normalize_comparator(
            _read(input_path),
            args.query_length,
            truth_intervals=truth,
            query_deletion_operation=args.query_deletion_operation,
        )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        serialized = json.dumps(result, indent=2 if args.pretty else None, sort_keys=True) + "\n"
        with output_path.open("x", encoding="utf-8") as handle:
            handle.write(serialized)
    except (ComparatorError, ExperimentError, IntervalError, OSError, json.JSONDecodeError) as exc:
        parser.exit(1, f"comparator normalization failed: {exc}\n")
    print(json.dumps({"status": "pass", "output": output_path.as_posix()}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
