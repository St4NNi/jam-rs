#!/usr/bin/env python3
"""Validate production benchmark artifacts without interpreting absence as negative."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "1.0.0"
SECRET_PATTERN = re.compile(r"(?i)(x-amz-signature|x-amz-credential|x-amz-security-token|password|token)=([^&\s]+)|://[^/\s:@]+:[^/\s@]+@")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fail(errors: list[str], message: str) -> None:
    errors.append(message)


def validate_measurement(path: Path, errors: list[str]) -> dict[str, Any] | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        fail(errors, f"{path}: invalid JSON: {error}")
        return None
    if not isinstance(value, dict):
        fail(errors, f"{path}: measurement is not an object")
        return None
    if value.get("schema_version") != SCHEMA_VERSION:
        fail(errors, f"{path}: unsupported schema version")
    if value.get("status") not in {"complete", "failed", "unavailable"}:
        fail(errors, f"{path}: invalid status")
    command = value.get("command", [])
    if not isinstance(command, list) or not all(isinstance(item, str) for item in command):
        fail(errors, f"{path}: command must be a string list")
    if any(SECRET_PATTERN.search(item) for item in command if isinstance(item, str)):
        fail(errors, f"{path}: command contains a credential-like query parameter")
    if value.get("status") == "complete":
        for field in ("wall_seconds", "user_seconds", "system_seconds", "peak_rss_kib", "resource_metrics"):
            if field not in value:
                fail(errors, f"{path}: complete measurement missing {field}")
        trace_path = value.get("trace_output")
        if trace_path:
            candidate = Path(trace_path)
            if not candidate.is_file():
                candidate = path.parent / candidate.name
            if not candidate.is_file():
                fail(errors, f"{path}: trace output is missing: {trace_path}")
    if value.get("status") == "unavailable" and not value.get("unavailable_reason"):
        fail(errors, f"{path}: unavailable measurement has no reason")
    return value


def validate_checksums(output: Path, errors: list[str]) -> None:
    path = output / "checksums.json"
    if not path.is_file():
        return
    try:
        checksums = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        fail(errors, f"{path}: invalid checksum JSON: {error}")
        return
    if not isinstance(checksums, dict):
        fail(errors, f"{path}: checksums must be an object")
        return
    for relative, expected in checksums.items():
        candidate = output / relative
        if not candidate.is_file():
            fail(errors, f"checksum target is missing: {relative}")
        elif sha256(candidate) != expected:
            fail(errors, f"checksum mismatch: {relative}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--regenerate-summary", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    raw = output / "raw"
    errors: list[str] = []
    if not raw.is_dir():
        fail(errors, f"raw directory is missing: {raw}")
    records: list[dict[str, Any]] = []
    seen: set[str] = set()
    if raw.is_dir():
        for path in sorted(raw.rglob("measurement.json")):
            record = validate_measurement(path, errors)
            if record is None:
                continue
            key = "|".join(str(record.get(field, "")) for field in ("track", "resource_mode", "query_id", "order", "cpu_threads", "io_tasks", "sweep_dimension", "sweep_value", "cache_state", "repeat"))
            if key in seen:
                fail(errors, f"duplicate raw measurement key: {key}")
            seen.add(key)
            records.append(record)
    validate_checksums(output, errors)
    if args.regenerate_summary and raw.is_dir():
        with tempfile.TemporaryDirectory(prefix="trace-summary-", dir=output.parent) as temporary:
            summary_dir = Path(temporary)
            completed = subprocess.run([sys.executable, str(Path(__file__).with_name("summarize.py")), "--raw-dir", str(raw), "--output-dir", str(summary_dir)], check=False)
            if completed.returncode != 0:
                fail(errors, "summary regeneration failed")
            summary = summary_dir / "summary.json"
            if not summary.is_file():
                fail(errors, "summary regeneration did not write summary.json")
            else:
                existing = output / "summary.json"
                if not existing.is_file():
                    existing = output / "summary" / "summary.json"
                if existing.is_file():
                    try:
                        generated_value = json.loads(summary.read_text(encoding="utf-8"))
                        existing_value = json.loads(existing.read_text(encoding="utf-8"))
                        for field in ("raw_measurement_count", "raw_measurements", "groups"):
                            if generated_value.get(field) != existing_value.get(field):
                                fail(errors, f"stored summary differs from raw regeneration in {field}")
                    except (OSError, ValueError) as error:
                        fail(errors, f"stored summary is not valid JSON: {error}")
    result = {
        "schema_version": SCHEMA_VERSION,
        "status": "valid" if not errors else "invalid",
        "measurement_count": len(records),
        "errors": errors,
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
