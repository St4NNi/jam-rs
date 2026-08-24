#!/usr/bin/env python3
"""Assert the normalized 40-case baseline against its checked-in fixture."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.trace_failure_analysis.baseline import BaselineError, assert_baseline


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--normalized", required=True, type=Path)
    parser.add_argument(
        "--fixture",
        type=Path,
        default=ROOT / "tests/data/archive_seeds/baseline_40_case.json",
    )
    args = parser.parse_args()
    try:
        result = assert_baseline(args.normalized, args.fixture)
    except BaselineError as exc:
        parser.exit(1, f"baseline assertion failed: {exc}\n")
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
