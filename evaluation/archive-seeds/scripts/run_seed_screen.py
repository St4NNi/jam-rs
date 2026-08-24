#!/usr/bin/env python3
"""Plan or run the staged seed-method correctness screen."""

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
    parser.add_argument("--execute", action="store_true", help="execute commands instead of planning")
    parser.add_argument("--pretty", action="store_true")
    args = parser.parse_args()
    try:
        result = run_manifest(load_manifest(args.manifest, "seed_screen"), execute=args.execute)
        write_new(args.output, result, pretty=args.pretty)
    except (ExperimentError, OSError) as exc:
        parser.exit(1, f"seed screen failed: {exc}\n")
    print(json.dumps({"status": "pass", "mode": result["mode"], "cases": len(result["cases"])}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
