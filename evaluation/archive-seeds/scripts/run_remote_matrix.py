#!/usr/bin/env python3
"""Plan or run bounded remote range-query measurements."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from common import ExperimentError, ensure_new_output, load_manifest, run_manifest, write_new


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--pretty", action="store_true")
    parser.add_argument("--max-case-seconds", type=float, default=60.0)
    args = parser.parse_args()
    try:
        ensure_new_output(args.output)
        result = run_manifest(
            load_manifest(args.manifest, "remote_matrix"),
            execute=args.execute,
            max_case_seconds=args.max_case_seconds,
        )
        write_new(args.output, result, pretty=args.pretty)
    except (ExperimentError, OSError) as exc:
        parser.exit(1, f"remote matrix failed: {exc}\n")
    print(json.dumps({"status": "pass", "mode": result["mode"], "cases": len(result["cases"])}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
