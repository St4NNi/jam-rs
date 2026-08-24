"""Bounded, manifest-driven experiment runner shared by benchmark scripts."""

from __future__ import annotations

import hashlib
import json
import os
import resource
import subprocess
import time
from pathlib import Path
from typing import Any


class ExperimentError(ValueError):
    """Raised when a benchmark manifest or output is unsafe or malformed."""


def load_manifest(path: str | Path, expected_kind: str) -> dict[str, Any]:
    manifest_path = Path(path)
    try:
        value = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ExperimentError(f"cannot read manifest {manifest_path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ExperimentError("experiment manifest must be an object")
    kind = str(value.get("kind", ""))
    if kind != expected_kind:
        raise ExperimentError(f"manifest kind {kind!r} does not match {expected_kind!r}")
    cases = value.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ExperimentError("experiment manifest requires a non-empty cases list")
    return value


def _command(raw: Any, index: int) -> list[str]:
    if not isinstance(raw, list) or not raw or not all(isinstance(item, str) and item for item in raw):
        raise ExperimentError(f"cases[{index}].command must be a non-empty string list")
    return list(raw)


def _parse_json_output(stdout: bytes) -> Any:
    text = stdout.decode("utf-8", errors="replace").strip()
    if not text:
        return None
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        records: list[Any] = []
        for line in text.splitlines():
            if not line.strip():
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError:
                return None
        return records or None


def _child_cpu_seconds() -> float:
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    return float(usage.ru_utime + usage.ru_stime)


def run_manifest(manifest: dict[str, Any], *, execute: bool) -> dict[str, Any]:
    """Plan or execute manifest cases without shell interpolation.

    Dry-run is the default for all stage scripts.  Execution is explicit and
    bounded by the commands and environment supplied in the manifest.
    """

    rows: list[dict[str, Any]] = []
    for index, case in enumerate(manifest["cases"]):
        if not isinstance(case, dict):
            raise ExperimentError(f"cases[{index}] must be an object")
        command = _command(case.get("command"), index)
        label = str(case.get("name", f"case-{index}"))
        row: dict[str, Any] = {
            "name": label,
            "command": command,
            "metadata": case.get("metadata", {}),
            "status": "planned",
        }
        if not execute:
            rows.append(row)
            continue
        cwd = case.get("cwd")
        if cwd is not None and not isinstance(cwd, str):
            raise ExperimentError(f"cases[{index}].cwd must be a string")
        environment = os.environ.copy()
        supplied_environment = case.get("env", {})
        if not isinstance(supplied_environment, dict) or not all(
            isinstance(key, str) and isinstance(value, str)
            for key, value in supplied_environment.items()
        ):
            raise ExperimentError(f"cases[{index}].env must map strings to strings")
        environment.update(supplied_environment)
        before_cpu = _child_cpu_seconds()
        started = time.perf_counter_ns()
        try:
            completed = subprocess.run(
                command,
                cwd=cwd,
                env=environment,
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        except OSError as exc:
            row.update({"status": "failed", "error": str(exc), "returncode": None})
            rows.append(row)
            continue
        elapsed = time.perf_counter_ns() - started
        after_cpu = _child_cpu_seconds()
        row.update(
            {
                "status": "complete" if completed.returncode == 0 else "failed",
                "returncode": completed.returncode,
                "wall_micros": elapsed // 1_000,
                "cpu_micros": round(max(0.0, after_cpu - before_cpu) * 1_000_000),
                "stdout_bytes": len(completed.stdout),
                "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(),
                "stderr_bytes": len(completed.stderr),
                "stderr_tail": completed.stderr.decode("utf-8", errors="replace")[-2000:],
                "result": _parse_json_output(completed.stdout),
            }
        )
        rows.append(row)
    return {
        "schema_version": "archive-seeds-experiment-v1",
        "kind": manifest["kind"],
        "mode": "execute" if execute else "plan",
        "manifest": manifest.get("name", "unnamed"),
        "cases": rows,
    }


def write_new(path: str | Path, value: Any, *, pretty: bool = False) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(value, indent=2 if pretty else None, sort_keys=True) + "\n"
    try:
        with output.open("x", encoding="utf-8") as handle:
            handle.write(text)
    except FileExistsError as exc:
        raise ExperimentError(f"refusing to overwrite existing output: {output}") from exc
