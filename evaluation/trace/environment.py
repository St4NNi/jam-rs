#!/usr/bin/env python3
"""Record benchmark environment and tool versions without running setup."""

from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def command_output(command: list[str]) -> str | None:
    try:
        completed = subprocess.run(command, check=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    except OSError:
        return None
    text = completed.stdout.strip()
    return text if completed.returncode == 0 else text or None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jam", type=Path, required=True)
    parser.add_argument("--dataset", type=Path, required=True)
    args = parser.parse_args()
    cpu = {}
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.exists():
        for line in cpuinfo.read_text(encoding="utf-8").splitlines():
            if ":" in line:
                key, value = (part.strip() for part in line.split(":", 1))
                if key in {"model name", "cpu cores", "siblings"} and key not in cpu:
                    cpu[key] = value
    tools = {}
    for tool in ("blastn", "minimap2", "curl", "zstd"):
        executable = shutil.which(tool)
        tools[tool] = command_output([executable, "--version"]) if executable else None
    data = {
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "cpu_count": os.cpu_count(),
        "cpu": cpu,
        "python": platform.python_version(),
        "rustc": command_output(["rustc", "--version", "--verbose"]),
        "cargo": command_output(["cargo", "--version"]),
        "jam": command_output([str(args.jam), "--version"]),
        "source_commit": command_output(["git", "rev-parse", "HEAD"]),
        "source_dirty": bool(command_output(["git", "status", "--porcelain"])),
        "release_flags": {
            "profile": "release",
            "opt_level": 3,
            "lto": True,
            "codegen_units": 1,
            "panic": "abort",
        },
        "tools": tools,
        "dataset_metadata": json.loads((args.dataset / "dataset.json").read_text(encoding="utf-8")),
        "cache_policy": "process-cold is the first timed invocation without a harness warm-up; warm is after one unmeasured invocation; kernel page-cache eviction is not attempted",
        "timing_policy": "compiler, installation, environment creation, and comparator setup are outside timed commands",
    }
    args.output.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
