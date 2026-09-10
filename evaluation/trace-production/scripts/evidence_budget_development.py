#!/usr/bin/env python3
"""Create query-only development fixtures after a target collection is frozen.

The target manifest and indexes must exist before this program is run. This
program reads target sequence but never writes or changes target material.
"""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


FORMAT = "jam-evidence-budget-development-v1"
BASES = b"ACGT"
CASES = (
    "exact80", "exact160", "reverse", "origin_crossing",
    "retained_50", "retained_10", "retained_05",
    "fragment_160", "fragment_400", "fragment_1600",
    "contigs_1", "contigs_2", "contigs_8",
    "substitution_dispersed", "substitution_clustered", "periodic_k21",
    "indel_short", "indel_long", "fragmented", "common", "shared", "negative",
)


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def dna(key: bytes, label: str, length: int) -> bytes:
    output = bytearray()
    counter = 0
    while len(output) < length:
        block = hashlib.sha256(key + b"\0" + label.encode() + counter.to_bytes(8, "little")).digest()
        output.extend(BASES[value & 3] for value in block)
        counter += 1
    return bytes(output[:length])


def reverse_complement(sequence: bytes) -> bytes:
    return sequence.translate(bytes.maketrans(b"ACGT", b"TGCA"))[::-1]


def canonical_kmers(sequence: bytes, k: int) -> set[bytes]:
    return {min(word, reverse_complement(word))
            for word in (sequence[start:start + k] for start in range(len(sequence) - k + 1))}


def intervals(length: int, count: int, width: int) -> list[list[int]]:
    gap = (length - count * width) // (count + 1)
    if gap < 1:
        raise ValueError("query is too short for the requested components")
    return [[gap + index * (gap + width), gap + index * (gap + width) + width]
            for index in range(count)]


def case_spec(case: str) -> tuple[int, list[list[int]], str, str]:
    if case == "exact80":
        return 480, [[191, 271]], "exact", "forward"
    if case == "exact160":
        return 720, [[211, 371]], "exact", "forward"
    if case == "reverse":
        return 2_000, [[500, 900]], "exact", "reverse"
    if case == "origin_crossing":
        return 2_000, [[1_800, 2_000], [0, 200]], "exact", "forward"
    if case.startswith("retained_"):
        length = {"retained_50": 6_400, "retained_10": 32_000, "retained_05": 64_000}[case]
        return length, intervals(length, 8, 400), "exact", "forward"
    if case.startswith("fragment_"):
        width = int(case.removeprefix("fragment_"))
        return 32_000, intervals(32_000, 3_200 // width, width), "exact", "forward"
    if case.startswith("contigs_"):
        return 32_000, intervals(32_000, 8, 400), "exact", "forward"
    if case in ("substitution_dispersed", "substitution_clustered", "periodic_k21"):
        return 2_000, [[400, 1_200]], case, "forward"
    if case in ("indel_short", "indel_long"):
        return 2_000, [[400, 1_200]], case, "forward"
    if case == "fragmented":
        return 4_000, [[400, 800], [2_400, 2_800]], "exact", "forward"
    if case in ("common", "shared"):
        return 2_000, [[600, 1_000]], "exact", "forward"
    if case == "negative":
        return 2_000, [], "none", "forward"
    raise ValueError(f"unknown case: {case}")


def write_fasta(path: Path, records: list[tuple[str, bytes]]) -> None:
    with path.open("xb") as stream:
        for name, sequence in records:
            stream.write(b">" + name.encode() + b"\n")
            for start in range(0, len(sequence), 80):
                stream.write(sequence[start:start + 80] + b"\n")


def target_slice(samtools: Path, target: dict, contig: str, start: int, end: int) -> bytes:
    if ":" in contig or start < 0 or end <= start:
        raise ValueError("unsafe target interval")
    command = [str(samtools), "faidx", target["bgzf"], f"{contig}:{start + 1}-{end}"]
    result = subprocess.run(command, check=False, capture_output=True)
    if result.returncode:
        raise ValueError(f"target extraction failed: {target['name']}/{contig}")
    sequence = b"".join(line.strip() for line in result.stdout.splitlines() if not line.startswith(b">"))
    sequence = sequence.upper().replace(b"U", b"T")
    if len(sequence) != end - start or any(base not in BASES for base in sequence):
        raise ValueError("target interval is not unambiguous ACGT")
    return sequence


def validate_manifest(path: Path, expected: str) -> tuple[dict, list[dict]]:
    if digest(path) != expected:
        raise ValueError("target manifest identity differs")
    manifest = json.loads(path.read_text(encoding="utf-8"))
    targets = manifest.get("targets", [])
    if (len(targets) != 12 or len({item["name"] for item in targets}) != 12
            or len({item.get("source_id") for item in targets}) != 12):
        raise ValueError("development target manifest must contain 12 unique members")
    for target in targets:
        if any(not target.get(field) for field in ("name", "source_id", "bgzf", "fai", "gzi",
                                                    "fasta", "selected_contigs", "selected_bases")):
            raise ValueError("incomplete target manifest row")
        for field in ("bgzf", "fai", "gzi", "fasta"):
            file = Path(target[field])
            expected_file = target.get("sha256", {}).get(field)
            if not file.is_file() or file.is_symlink() or not expected_file or digest(file) != expected_file:
                raise ValueError(f"target file identity differs: {target['name']}/{field}")
        contigs = target["selected_contigs"]
        if not contigs or len({item["name"] for item in contigs}) != len(contigs):
            raise ValueError("each target must bind a nonempty unique selected-contig list")
        selected_bases = sum(int(item["length"]) for item in contigs)
        if selected_bases != int(target["selected_bases"]):
            raise ValueError("selected target base accounting differs")
        if [int(item["fai_ordinal"]) for item in contigs] != list(range(len(contigs))):
            raise ValueError("selected target contigs are not the contiguous FAI prefix")
        if selected_bases < 1_000_000:
            if target.get("quota_reached") is not False or not target.get("source_fai"):
                raise ValueError("sub-quota target lacks explicit source exhaustion")
            source_rows = []
            with Path(target["source_fai"]).open(encoding="utf-8") as stream:
                for line in stream:
                    fields = line.rstrip("\n").split("\t")
                    source_rows.append((fields[0], int(fields[1])))
            expected_rows = [(item["name"], int(item["length"])) for item in contigs]
            if source_rows != expected_rows:
                raise ValueError("sub-quota target did not exhaust its source FAI")
        elif target.get("quota_reached") is not True:
            raise ValueError("target at or above quota lacks the recorded quota boundary")
    return manifest, targets


def placement_pool(targets: list[dict], minimum_length: int = 1_600) -> list[tuple[dict, dict]]:
    return [(target, contig) for target in targets for contig in target["selected_contigs"]
            if int(contig["length"]) >= minimum_length]


def choose_group(key: bytes, label: str, targets: list[dict], count: int) -> list[tuple[dict, dict]]:
    eligible = [(target, [contig for contig in target["selected_contigs"]
                          if int(contig["length"]) >= 1_600]) for target in targets]
    eligible = [(target, contigs) for target, contigs in eligible if len(contigs) >= count]
    if not eligible:
        raise ValueError(f"no target member has {count} eligible contigs for {label}")
    target, contigs = min(eligible, key=lambda item: hashlib.sha256(
        key + b"\0member\0" + label.encode() + b"\0" + item[0]["name"].encode()).digest())
    ranked = sorted(contigs, key=lambda contig: hashlib.sha256(
        key + b"\0contig\0" + label.encode() + b"\0" + contig["name"].encode()).digest())
    return [(target, contig) for contig in ranked[:count]]


def target_interval(key: bytes, label: str, target: dict, contig: dict, width: int) -> tuple[int, int]:
    length = int(contig["length"])
    if length < width:
        raise ValueError("selected contig is shorter than the component")
    room = length - width + 1
    value = int.from_bytes(hashlib.sha256(key + b"\0start\0" + label.encode()).digest()[:8], "little")
    start = value % room
    return start, start + width


def choose_slice(samtools: Path, key: bytes, label: str, target: dict,
                 contig: dict, width: int) -> tuple[int, int, bytes]:
    initial, _ = target_interval(key, label, target, contig, width)
    room = int(contig["length"]) - width + 1
    step = max(1, room // 32)
    for attempt in range(min(32, room)):
        start = (initial + attempt * step) % room
        try:
            sequence = target_slice(samtools, target, contig["name"], start, start + width)
            return start, start + width, sequence
        except ValueError as error:
            if "not unambiguous ACGT" not in str(error):
                raise
    raise ValueError(f"no unambiguous bounded slice for {target['name']}/{contig['name']}")


def mutate(sequence: bytes, transform: str, key: bytes, label: str) -> tuple[bytes, list[list[int]]]:
    changed = bytearray(sequence)
    available = [[0, len(sequence)]]
    if transform == "substitution_dispersed":
        count = max(1, len(changed) // 31)
        positions = sorted(range(len(changed)), key=lambda position: hashlib.sha256(
            key + b"\0substitution\0" + label.encode() + position.to_bytes(8, "little")).digest())[:count]
    elif transform == "substitution_clustered":
        positions = range(320, min(380, len(changed)), 3)
    elif transform == "periodic_k21":
        positions = range(0, len(changed), 20)
    else:
        positions = ()
    for position in positions:
        changed[position] = BASES[(BASES.index(changed[position]) + 1) % 4]
    if transform == "periodic_k21" and canonical_kmers(sequence, 21) & canonical_kmers(bytes(changed), 21):
        raise ValueError("periodic development component retained an exact canonical k21")
    if transform in ("indel_short", "indel_long"):
        gap = 5 if transform == "indel_short" else 30
        query = sequence[:200] + sequence[200 + gap:500 + gap]
        query += dna(key, f"deleted:{label}", gap) + sequence[500 + gap:]
        available = [[0, 500], [500 + gap, len(sequence)]]
        return query, available
    return bytes(changed), available


def shared_group(audit: dict, case: str, replicate: int) -> dict:
    groups = audit.get("groups", {}).get(case, [])
    if len(groups) <= replicate:
        raise ValueError(f"frozen manifest lacks query-independent {case} group r{replicate:02d}")
    group = groups[replicate]
    placements = group.get("placements", [])
    by_member = {}
    for item in placements:
        by_member[item["name"]] = by_member.get(item["name"], 0) + 1
    if case == "common" and max(by_member.values(), default=0) < 2:
        raise ValueError("common group lacks two physical placements in one member")
    if case == "shared" and len(by_member) < 2:
        raise ValueError("shared group has fewer than two members")
    if not isinstance(group.get("sequence_sha256"), str) or len(group["sequence_sha256"]) != 64:
        raise ValueError("shared group lacks a sequence identity")
    return group


def map_available(local: list[list[int]], segments: list[list[int]]) -> list[list[int]]:
    output = []
    component_offset = 0
    for query_start, query_end in segments:
        segment_length = query_end - query_start
        for left, right in local:
            overlap_start = max(left, component_offset)
            overlap_end = min(right, component_offset + segment_length)
            if overlap_start < overlap_end:
                output.append([query_start + overlap_start - component_offset,
                               query_start + overlap_end - component_offset])
        component_offset += segment_length
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-manifest", type=Path, required=True)
    parser.add_argument("--expected-target-sha256", required=True)
    parser.add_argument("--shared-groups", type=Path, required=True)
    parser.add_argument("--expected-shared-groups-sha256", required=True)
    parser.add_argument("--family-seed-file", type=Path, required=True)
    parser.add_argument("--samtools", type=Path, required=True)
    parser.add_argument("--expected-samtools-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (args.output.exists() or not args.samtools.is_file() or args.samtools.is_symlink()
            or digest(args.samtools) != args.expected_samtools_sha256
            or not args.family_seed_file.is_file() or args.family_seed_file.is_symlink()):
        raise SystemExit("new output and existing samtools are required")
    key = args.family_seed_file.read_bytes()
    if len(key) < 32:
        raise SystemExit("development family seed must contain at least 32 bytes")
    manifest, targets = validate_manifest(args.target_manifest, args.expected_target_sha256)
    if digest(args.shared_groups) != args.expected_shared_groups_sha256:
        raise ValueError("shared-group audit identity differs")
    shared_audit = json.loads(args.shared_groups.read_text(encoding="utf-8"))
    if (shared_audit.get("format") != "jam-target-repeat-audit-v1"
            or shared_audit.get("target_manifest_sha256") != digest(args.target_manifest)):
        raise ValueError("shared-group audit is not bound to the frozen targets")
    by_name = {item["name"]: item for item in targets}
    pool = placement_pool(targets)
    if not pool:
        raise ValueError("frozen target corpus has no contig of at least 1,600 bases")
    query_records, linear_records, circular_records, truth = [], [], [], []
    unavailable = []
    available_cases = []
    for case in CASES:
        if case == "common" and not shared_audit.get("groups", {}).get(case):
            unavailable.append({"case": case, "queries_omitted": 2,
                                "reason": "no frozen target-only common group"})
            continue
        available_cases.append(case)
        for replicate in range(2):
            query_id = f"budget_development_{case}_r{replicate:02d}"
            length, query_ranges, transform, strand = case_spec(case)
            query = bytearray(dna(key, f"background:{query_id}", length))
            components = []
            if case in ("common", "shared"):
                group = shared_group(shared_audit, case, replicate)
                placements = group["placements"]
                shared_widths = {int(item["end"]) - int(item["start"]) for item in placements}
                if len(shared_widths) != 1 or next(iter(shared_widths)) not in (160, 400):
                    raise ValueError("shared group width must be exactly 160 or 400 bases")
                shared_width = next(iter(shared_widths))
                query_ranges = [[600, 600 + shared_width]]
                groups = [placements]
            elif case == "negative":
                groups = []
            else:
                distinct = (int(case.removeprefix("contigs_")) if case.startswith("contigs_")
                            else 2 if case == "fragmented" or case.startswith(("retained_", "fragment_")) else 1)
                selected = choose_group(key, query_id, targets, distinct)
                groups = [[{"name": selected[index % distinct][0]["name"],
                            "contig": selected[index % distinct][1]["name"]}]
                          for index in range(len(query_ranges) if case != "origin_crossing" else 1)]
            for index, placements in enumerate(groups):
                ranges = query_ranges if case == "origin_crossing" else [query_ranges[index]]
                width = sum(end - start for start, end in ranges)
                first = placements[0]
                target = by_name[first["name"]]
                contig = next(item for item in target["selected_contigs"] if item["name"] == first["contig"])
                start = first.get("start")
                end = first.get("end")
                if start is None or end is None:
                    start, end, target_bases = choose_slice(
                        args.samtools, key, f"{query_id}/{index}", target, contig, width
                    )
                else:
                    target_bases = target_slice(args.samtools, target, first["contig"], start, end)
                    if hashlib.sha256(target_bases).hexdigest() != group["sequence_sha256"]:
                        raise ValueError("shared group sequence identity differs")
                if end - start != width:
                    raise ValueError("shared placement width differs from query component")
                source_query = reverse_complement(target_bases) if strand == "reverse" else target_bases
                source_query, available_local = mutate(source_query, transform, key, f"{query_id}/{index}")
                offset = 0
                for left, right in ranges:
                    query[left:right] = source_query[offset:offset + right - left]
                    offset += right - left
                available = map_available(available_local, ranges)
                for placement in placements:
                    ptarget = by_name[placement["name"]]
                    pstart, pend = int(placement.get("start", start)), int(placement.get("end", end))
                    observed = target_slice(args.samtools, ptarget, placement["contig"], pstart, pend)
                    if observed != target_bases:
                        raise ValueError("shared group placements are not byte-identical")
                    components.append({"source_id": ptarget.get("source_id", ptarget["name"]),
                                       "contig": placement["contig"], "strand": strand,
                                       "target_interval": [pstart, pend], "query_segments": ranges,
                                       "available_query_intervals": available})
            record = (query_id, bytes(query))
            query_records.append(record)
            (circular_records if case == "origin_crossing" else linear_records).append(record)
            truth.append({"query_id": query_id, "split": "development", "case": case,
                          "replicate": replicate, "components": components})
    args.output.mkdir(mode=0o700)
    write_fasta(args.output / "queries.fasta", query_records)
    write_fasta(args.output / "queries.linear.fasta", linear_records)
    write_fasta(args.output / "queries.circular.fasta", circular_records)
    calibration = next(record for record in query_records
                       if record[0] == "budget_development_reverse_r00")
    write_fasta(args.output / "queries.single-positive.fasta", [calibration])
    truth_path = args.output / "truth.jsonl"
    truth_path.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in truth), encoding="utf-8")
    files = {name: {"bytes": (args.output / name).stat().st_size, "sha256": digest(args.output / name)}
             for name in ("queries.fasta", "queries.linear.fasta", "queries.circular.fasta",
                          "queries.single-positive.fasta", "truth.jsonl")}
    seal = {"format": FORMAT, "status": "complete", "split": "development",
            "queries": len(query_records), "requested_classes": list(CASES),
            "available_classes": available_cases, "unavailable_strata": unavailable,
            "replicates_per_class": 2,
            "target_manifest_sha256": digest(args.target_manifest),
            "shared_groups_sha256": digest(args.shared_groups),
            "family_seed_sha256": hashlib.sha256(key).hexdigest(), "files": files,
            "samtools_sha256": digest(args.samtools),
            "final_split_status": "not_created_or_revealed"}
    (args.output / "fixture-manifest.json").write_text(
        json.dumps(seal, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
