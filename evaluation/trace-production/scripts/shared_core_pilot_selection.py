#!/usr/bin/env python3
"""Freeze a query-independent whole-contig prefix selection for the pilot."""

import argparse
import hashlib
import json
import os
from pathlib import Path


SOURCE_SHA256 = "0282d1afef69eafd4a0d42875bf50f85a892672ec128923eeadbdc95f36ecd9b"
QUOTA = 27_000_000
LIMIT = 250_000_000


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-freeze", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise SystemExit("new output required")
    if digest(args.source_freeze) != SOURCE_SHA256:
        raise SystemExit("source freeze identity differs")
    source = json.loads(args.source_freeze.read_text(encoding="utf-8"))
    if (source.get("status") != "frozen"
            or source.get("query_dependent_selection") is not False
            or len(source.get("targets", [])) != 12):
        raise ValueError("source freeze is not the expected query-independent collection")
    targets = []
    total_bases = 0
    total_contigs = 0
    for target in source["targets"]:
        fai = Path(target["source_fai"])
        bgzf = Path(target["source_bgzf"])
        gzi = Path(target["source_gzi"])
        for path in (fai, bgzf, gzi):
            if not path.is_file() or path.is_symlink():
                raise ValueError(f"source resource differs: {path}")
        selected = []
        selected_bases = 0
        with fai.open(encoding="utf-8") as stream:
            for ordinal, line in enumerate(stream):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    raise ValueError(f"invalid FAI row: {fai}:{ordinal + 1}")
                name, length = fields[0], int(fields[1])
                if not name or length < 1:
                    raise ValueError("invalid FAI contig")
                selected.append({"name": name, "length": length, "fai_ordinal": ordinal})
                selected_bases += length
                if selected_bases >= QUOTA:
                    break
        original = [(item["name"], int(item["length"]), int(item["fai_ordinal"]))
                    for item in target["selected_contigs"]]
        expanded = [(item["name"], item["length"], item["fai_ordinal"])
                    for item in selected[:len(original)]]
        if original != expanded:
            raise ValueError("pilot does not preserve the complete frozen contig prefix")
        row = {
            "name": target["name"],
            "source_id": target["source_id"],
            "source_name": target["source_name"],
            "source_bgzf": str(bgzf),
            "source_fai": str(fai),
            "source_gzi": str(gzi),
            "source_identity": {
                "bgzf_bytes": bgzf.stat().st_size,
                "fai_bytes": fai.stat().st_size,
                "gzi_bytes": gzi.stat().st_size,
                "fai_sha256": digest(fai),
                "gzi_sha256": digest(gzi),
            },
            "selected_bases": selected_bases,
            "selected_contigs": selected,
            "quota_reached": selected_bases >= QUOTA,
            "original_fixture_bases": int(target["selected_bases"]),
            "original_fixture_contigs": len(original),
        }
        targets.append(row)
        total_bases += selected_bases
        total_contigs += len(selected)
    if total_bases > LIMIT:
        raise ValueError("pilot selection exceeds declared biological base limit")
    result = {
        "format": "jam-shared-core-pilot-selection-v1",
        "status": "frozen_before_candidate_timing_or_query_results",
        "selection": "whole contigs in recorded source FAI order through 27000000 cumulative bases per member, or source exhaustion",
        "query_dependent_selection": False,
        "source_freeze": str(args.source_freeze),
        "source_freeze_sha256": SOURCE_SHA256,
        "quota_bases_per_member": QUOTA,
        "limit_bases": LIMIT,
        "selected_bases": total_bases,
        "selected_contigs": total_contigs,
        "members": len(targets),
        "partial_assemblies": True,
        "final_split": "sealed and unused",
        "targets": targets,
    }
    args.output.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())


if __name__ == "__main__":
    main()
