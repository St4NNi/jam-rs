#!/usr/bin/env python3
"""Plan or run bounded remote range-query measurements."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from tools.trace_failure_analysis.experiment import ExperimentError, load_manifest, run_manifest, write_new


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--pretty", action="store_true")
    args = parser.parse_args()
    try:
        result = run_manifest(load_manifest(args.manifest, "remote_matrix"), execute=args.execute)
        write_new(args.output, result, pretty=args.pretty)
    except (ExperimentError, OSError) as exc:
        parser.exit(1, f"remote matrix failed: {exc}\n")
    print(json.dumps({"status": "pass", "mode": result["mode"], "cases": len(result["cases"])}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
