#!/usr/bin/env python3
"""Run LexicMap on mixed topology queries and project doubled circular evidence."""

import argparse
import copy
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path


DETAIL_COLUMNS = (
    "query", "qlen", "hits", "sgenome", "sseqid", "qcovGnm", "cls", "hsp",
    "qcovHSP", "alenHSP", "pident", "gaps", "qstart", "qend", "sstart", "send",
    "sstr", "slen", "evalue", "bitscore", "cigar", "qseq", "sseq", "align",
)
CIGAR = re.compile(r"([1-9][0-9]*)([M=XID])")
HEX64 = re.compile(r"[0-9a-f]{64}")
LEXICMAP_NORMALIZE = bytes.maketrans(b"ACGTUMRWSYKVHDBN", b"ACGTTAAACCGAAACA")

def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()

def sequence_digest(sequence: bytes) -> str:
    return hashlib.sha256(sequence).hexdigest()

def read_fasta(path: Path) -> list[dict]:
    records, header, sequence = [], None, bytearray()
    for line in path.read_bytes().splitlines():
        if line.startswith(b">"):
            if header is not None:
                records.append(fasta_record(header, bytes(sequence)))
            header, sequence = line[1:].decode(), bytearray()
        elif header is None:
            raise ValueError("sequence before query header")
        else:
            sequence.extend(line.strip().upper())
    if header is not None:
        records.append(fasta_record(header, bytes(sequence)))
    identifiers = [record["query_id"] for record in records]
    if not records or len(set(identifiers)) != len(identifiers):
        raise ValueError("empty query FASTA or duplicate query ID")
    return records

def fasta_record(header: str, sequence: bytes) -> dict:
    tokens = header.split()
    topology = [token.removeprefix("topology=") for token in tokens[1:]
                if token.startswith("topology=")]
    if not tokens or not sequence or topology not in (["linear"], ["circular"]):
        raise ValueError("each query requires one topology=linear or topology=circular token")
    return {"query_id": tokens[0], "header": header, "topology": topology[0],
            "sequence": sequence}

def write_transformed(path: Path, records: list[dict]) -> list[dict]:
    manifest = []
    with path.open("x", encoding="ascii") as stream:
        for record in records:
            biological = record["sequence"]
            physical = transformed_sequence(record)
            stream.write(f">{record['header']}\n")
            for start in range(0, len(physical), 80):
                stream.write(physical[start:start + 80].decode("ascii") + "\n")
            manifest.append({
                "query_id": record["query_id"], "topology": record["topology"],
                "biological_bases": len(biological), "physical_bases": len(physical),
                "source_sequence_sha256": sequence_digest(biological),
                "transformed_sequence_sha256": sequence_digest(physical),
            })
    return manifest

def transformed_sequence(record: dict) -> bytes:
    return record["sequence"] * (2 if record["topology"] == "circular" else 1)

def read_native(path: Path, records: list[dict]) -> dict:
    expected = {record["query_id"]: record for record in records}
    counts = {"native_hsps": 0, "over_one_traversal_hsps": 0}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        rows = csv.reader((line for line in stream if not line.startswith("#")), delimiter="\t")
        if tuple(next(rows, ())) != DETAIL_COLUMNS:
            raise ValueError("unexpected LexicMap detail columns")
        for number, values in enumerate(rows, 2):
            if len(values) != len(DETAIL_COLUMNS):
                raise ValueError(f"wrong LexicMap column count at row {number}")
            row = dict(zip(DETAIL_COLUMNS, values, strict=True))
            record = expected.get(row["query"])
            if record is None:
                raise ValueError("unknown LexicMap query")
            physical = transformed_sequence(record)
            qstart, qend = int(row["qstart"]) - 1, int(row["qend"])
            if int(row["qlen"]) != len(physical) or not (0 <= qstart < qend <= len(physical)):
                raise ValueError("LexicMap transformed query coordinates differ")
            query = row["qseq"].replace("-", "").encode().upper().translate(LEXICMAP_NORMALIZE)
            expected_query = physical[qstart:qend].translate(LEXICMAP_NORMALIZE)
            if query != expected_query:
                raise ValueError("LexicMap query alignment differs from transformed input")
            runs = [(int(length), operation) for length, operation in CIGAR.findall(row["cigar"])]
            if not runs or "".join(f"{length}{operation}" for length, operation in runs) != row["cigar"]:
                raise ValueError("invalid LexicMap CIGAR")
            counts["native_hsps"] += 1
            counts["over_one_traversal_hsps"] += int(qend - qstart > len(record["sequence"]))
    return counts

def atomic_json(path: Path, value: dict) -> None:
    if path.exists() or path.is_symlink():
        raise ValueError(f"output exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent,
                                     prefix=f".{path.name}.", delete=False) as stream:
        temporary = Path(stream.name)
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    try:
        os.link(temporary, path)
        sync_directory(path.parent)
    finally:
        temporary.unlink(missing_ok=True)


def sync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def union(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    output = []
    for start, end in sorted(intervals):
        if start >= end:
            raise ValueError("invalid query interval")
        if output and start <= output[-1][1]:
            output[-1] = (output[-1][0], max(output[-1][1], end))
        else:
            output.append((start, end))
    return output


def complement(intervals: list[tuple[int, int]], length: int) -> list[tuple[int, int]]:
    output, start = [], 0
    for left, right in intervals:
        if start < left:
            output.append((start, left))
        start = right
    if start < length:
        output.append((start, length))
    return output


def projected_segments(start: int, span: int, length: int) -> list[tuple[int, int]]:
    if not (0 <= start < 2 * length and 0 < span <= length and start + span <= 2 * length):
        raise ValueError("invalid transformed query interval")
    left = start % length
    if left + span <= length:
        return [(left, left + span)]
    return [(left, length), (0, left + span - length)]


def support(fragment: dict) -> list[tuple[int, int]]:
    positions = [position for segment in fragment["query_segments"]
                 for position in range(int(segment["start"]), int(segment["end"]))]
    runs = [(int(item["length"]), item["operation"])
            for item in fragment["alignment"]["edit_script"]]
    if "".join(f"{length}{operation}" for length, operation in runs) \
            != fragment["alignment"]["cigar"]:
        raise ValueError("CIGAR and edit script differ")
    selected, offset = [], 0
    for length, operation in runs:
        if operation in ("=", "X"):
            selected.extend(positions[offset:offset + length])
        if operation in ("=", "X", "D"):
            offset += length
    if offset != len(positions):
        raise ValueError("CIGAR query consumption differs")
    return union([(position, position + 1) for position in selected])


def fragment_identity(fragment: dict) -> str:
    alignment = fragment["alignment"]
    value = {
        "contig_id": fragment["contig_id"], "query_segments": fragment["query_segments"],
        "alignment": {key: alignment[key] for key in (
            "score", "strand", "target_interval", "matches", "substitutions", "insertions",
            "deletions", "cigar", "edit_script", "identity")},
    }
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def project_records(rows: dict[str, dict], transform: dict, transform_sha256: str) -> dict:
    bindings = {record["query_id"]: record for record in transform["records"]}
    if set(rows) != set(bindings):
        raise ValueError("projected query IDs differ from transform manifest")
    totals = {"input_fragments": 0, "excluded_over_one_traversal": 0,
              "deduplicated_copy_fragments": 0, "retained_fragments": 0,
              "origin_spanning_fragments": 0}
    for query_id, row in rows.items():
        binding = bindings[query_id]
        length = int(binding["biological_bases"])
        if row.get("query_length") != int(binding["physical_bases"]):
            raise ValueError("normalized result query length differs from transform")
        row["query_length"] = length
        if binding["topology"] == "linear":
            fragments = sum(
                len(wrapped.get("trace", wrapped)["mosaic"]["primary"])
                + len(wrapped.get("trace", wrapped)["mosaic"]["alternatives"])
                for wrapped in row["metagenomes"])
            totals["input_fragments"] += fragments
            totals["retained_fragments"] += fragments
            row["lexicmap_circular_projection"] = {
                "format": "jam-lexicmap-circular-projection-v1",
                "transform_sha256": transform_sha256,
            }
            continue
        retained_metagenomes = []
        for wrapped in row["metagenomes"]:
            trace = wrapped.get("trace", wrapped)
            mosaic = trace["mosaic"]
            candidates = [entry["fragment"] for entry in mosaic["primary"]] + mosaic["alternatives"]
            totals["input_fragments"] += len(candidates)
            prepared, identities = [], set()
            for fragment in candidates:
                alignment = fragment["alignment"]
                segments = fragment["query_segments"]
                if len(segments) != 1:
                    raise ValueError("normalized LexicMap fragment is not linear")
                start, end = int(segments[0]["start"]), int(segments[0]["end"])
                span = end - start
                if span > length:
                    totals["excluded_over_one_traversal"] += 1
                    continue
                if binding["topology"] == "circular":
                    projected = projected_segments(start, span, length)
                    totals["origin_spanning_fragments"] += int(len(projected) == 2)
                    fragment["query_segments"] = [
                        {"start": left, "end": right} for left, right in projected]
                    alignment["query_interval"] = {"start": 0, "end": span}
                    alignment.setdefault("lexicmap", {}).update({
                        "physical_query_interval": {"start": start, "end": end},
                        "circular_transform_sha256": transform_sha256,
                    })
                identity = fragment_identity(fragment)
                if identity in identities:
                    totals["deduplicated_copy_fragments"] += 1
                    continue
                identities.add(identity)
                prepared.append((fragment, support(fragment)))
            covered, primary, alternatives = [], [], []
            for fragment, fragment_support in prepared:
                before = sum(end - start for start, end in covered)
                combined = union(covered + fragment_support)
                added = sum(end - start for start, end in combined) - before
                if added:
                    covered = combined
                    primary.append({"newly_supported_bases": added, "fragment": fragment})
                else:
                    alternatives.append(fragment)
            if not primary and not alternatives:
                continue
            totals["retained_fragments"] += len(primary) + len(alternatives)
            mosaic.update({
                "query_length": length, "covered_bases": sum(end - start for start, end in covered),
                "covered_intervals": [{"start": start, "end": end} for start, end in covered],
                "gaps": [{"start": start, "end": end} for start, end in complement(covered, length)],
                "primary": primary, "alternatives": alternatives,
            })
            retained_metagenomes.append(wrapped)
        row["metagenomes"] = retained_metagenomes
        row["lexicmap_circular_projection"] = {
            "format": "jam-lexicmap-circular-projection-v1",
            "transform_sha256": transform_sha256,
        }
    return totals


def load_transform(path: Path, expected_sha256: str, original_query: Path,
                   transformed_query: Path) -> dict:
    if not HEX64.fullmatch(expected_sha256) or digest(path) != expected_sha256:
        raise ValueError("circular transform manifest identity differs")
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("format") != "jam-lexicmap-doubled-query-v1" \
            or value.get("source_query_sha256") != digest(original_query):
        raise ValueError("circular transform manifest binding differs")
    records = value.get("records")
    if not isinstance(records, list) or len({item.get("query_id") for item in records}) != len(records):
        raise ValueError("invalid circular transform records")
    source = {record["query_id"]: record for record in read_fasta(original_query)}
    if {item.get("query_id") for item in records} != set(source):
        raise ValueError("circular transform query IDs differ")
    for item in records:
        biological, physical = item.get("biological_bases"), item.get("physical_bases")
        expected = biological * (2 if item.get("topology") == "circular" else 1)
        record = source[item["query_id"]]
        if (not isinstance(biological, int) or biological <= 0 or physical != expected
                or record["topology"] != item.get("topology")
                or len(record["sequence"]) != biological
                or sequence_digest(record["sequence"]) != item.get("source_sequence_sha256")
                or sequence_digest(transformed_sequence(record))
                != item.get("transformed_sequence_sha256")):
            raise ValueError("invalid circular transform record")
    transformed = value.get("transformed_query", {})
    binary_path = Path(value.get("binary", ""))
    native_path = path.parent / value.get("native_output_name", "")
    if (value.get("adapter_sha256") != digest(Path(__file__))
            or not binary_path.is_file() or binary_path.is_symlink()
            or digest(binary_path) != value.get("binary_sha256")
            or not transformed_query.is_file() or transformed_query.is_symlink()
            or digest(transformed_query) != transformed.get("sha256")
            or not native_path.is_file() or native_path.is_symlink()
            or digest(native_path) != value.get("native_output_sha256")):
        raise ValueError("circular transform retained files differ")
    value["offline_native_validation"] = read_native(native_path, list(source.values()))
    return value


def self_test() -> None:
    linear = {"topology": "linear", "sequence": b"ACGT"}
    circular = {"topology": "circular", "sequence": b"ACGT"}
    assert transformed_sequence(linear) == b"ACGT"
    assert transformed_sequence(circular) == b"ACGTACGT"
    transform = {"records": [
        {"query_id": "q", "topology": "circular", "biological_bases": 10,
         "physical_bases": 20},
    ]}
    alignment = {"score": -2, "strand": "reverse", "query_interval": {"start": 8, "end": 14},
                 "target_interval": {"start": 20, "end": 26}, "matches": 5,
                 "substitutions": 0, "insertions": 1, "deletions": 1,
                 "cigar": "2=1D1I3=", "edit_script": [
                     {"length": 2, "operation": "="}, {"length": 1, "operation": "D"},
                     {"length": 1, "operation": "I"}, {"length": 3, "operation": "="}],
                 "identity": 5 / 7}
    fragment = {"contig_id": 3, "query_segments": [{"start": 8, "end": 14}],
                "alignment": alignment}
    copy_fragment = copy.deepcopy(fragment)
    copy_fragment["query_segments"] = [{"start": 2, "end": 5}]
    copy_fragment["alignment"].update({
        "query_interval": {"start": 2, "end": 5}, "target_interval": {"start": 30, "end": 33},
        "score": 6, "matches": 3, "insertions": 0, "deletions": 0, "cigar": "3=",
        "edit_script": [{"length": 3, "operation": "="}], "identity": 1.0,
    })
    second_copy = copy.deepcopy(copy_fragment)
    second_copy["query_segments"] = [{"start": 12, "end": 15}]
    second_copy["alignment"]["query_interval"] = {"start": 12, "end": 15}
    over = copy.deepcopy(fragment)
    over["query_segments"] = [{"start": 1, "end": 12}]
    over["alignment"]["query_interval"] = {"start": 1, "end": 12}
    row = {"query_id": "q", "query_length": 20, "metagenomes": [{"trace": {"mosaic": {
        "primary": [{"fragment": fragment}, {"fragment": copy_fragment},
                    {"fragment": second_copy}],
        "alternatives": [over], "covered_intervals": [], "gaps": [], "covered_bases": 0,
        "query_length": 20}}}]}
    totals = project_records({"q": row}, transform, "a" * 64)
    projected = row["metagenomes"][0]["trace"]["mosaic"]["primary"][0]["fragment"]
    assert projected["query_segments"] == [{"start": 8, "end": 10}, {"start": 0, "end": 4}]
    assert projected["alignment"]["strand"] == "reverse"
    assert projected["alignment"]["cigar"] == "2=1D1I3="
    assert totals == {"input_fragments": 4, "excluded_over_one_traversal": 1,
                      "deduplicated_copy_fragments": 1, "retained_fragments": 2,
                      "origin_spanning_fragments": 1}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--prevalidated-binary-sha256")
    parser.add_argument("--prevalidated-adapter-sha256")
    parser.add_argument("--index", type=Path)
    parser.add_argument("--query", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--threads", type=int)
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if any(value is None for value in (
            args.binary, args.prevalidated_binary_sha256, args.prevalidated_adapter_sha256,
            args.index, args.query, args.output, args.work_dir, args.threads)):
        raise SystemExit("all adapter arguments are required")
    manifest_output = Path(str(args.output) + ".transform.json")
    if (not HEX64.fullmatch(args.prevalidated_binary_sha256)
            or not HEX64.fullmatch(args.prevalidated_adapter_sha256)
            or not args.binary.is_file()
            or args.binary.is_symlink()
            or not args.index.exists() or not args.query.is_file() or args.query.is_symlink()
            or args.output.exists() or manifest_output.exists() or args.work_dir.exists()
            or args.threads <= 0):
        raise SystemExit("new outputs and exact existing inputs are required")
    args.work_dir.mkdir(mode=0o700, parents=True)
    temporary = args.work_dir / "temp"
    temporary.mkdir(mode=0o700)
    transformed = args.work_dir / "queries.lexicmap.fasta"
    records = read_fasta(args.query)
    record_manifest = write_transformed(transformed, records)
    native = args.work_dir / "native.tsv.gz"
    command = [str(args.binary), "--quiet", "--threads", str(args.threads), "search",
               "--index", str(args.index), "--align-min-match-pident", "80",
               "--align-min-match-len", "40", "--min-qcov-per-hsp", "0",
               "--min-qcov-per-genome", "0", "--max-evalue", "1e300",
               "--top-n-chains", "0", "--top-n-genomes", "0", "--max-query-conc", "8",
               "--out-file", str(native), "--all", str(transformed)]
    environment = os.environ.copy()
    environment.update({"TMPDIR": str(temporary), "TMP": str(temporary), "TEMP": str(temporary)})
    with (args.work_dir / "stdout").open("x") as stdout, \
            (args.work_dir / "stderr").open("x") as stderr:
        run = subprocess.run(command, stdout=stdout, stderr=stderr, env=environment, check=False)
    if run.returncode or not native.is_file():
        raise SystemExit(f"LexicMap failed with exit code {run.returncode}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with native.open("rb") as source, tempfile.NamedTemporaryFile(
            "wb", dir=args.output.parent, prefix=f".{args.output.name}.", delete=False) as stream:
        temporary_output = Path(stream.name)
        shutil.copyfileobj(source, stream)
        stream.flush()
        os.fsync(stream.fileno())
    try:
        os.link(temporary_output, args.output)
        sync_directory(args.output.parent)
    finally:
        temporary_output.unlink(missing_ok=True)
    manifest = {
        "format": "jam-lexicmap-doubled-query-v1", "status": "complete",
        "adapter_sha256": args.prevalidated_adapter_sha256,
        "binary": str(args.binary), "binary_sha256": args.prevalidated_binary_sha256,
        "binary_identity_mode": "prevalidated outside timed adapter",
        "binary_stat": {key: value for key, value in zip(
            ("device", "inode", "bytes", "mtime_ns"),
            (args.binary.stat().st_dev, args.binary.stat().st_ino, args.binary.stat().st_size,
             args.binary.stat().st_mtime_ns), strict=True)},
        "source_query_sha256": digest(args.query),
        "transformed_query": {"name": transformed.name, "sha256": digest(transformed)},
        "native_output_name": args.output.name, "native_output_sha256": digest(args.output),
        "records": record_manifest, "biological_query_bases": sum(
            row["biological_bases"] for row in record_manifest),
        "physical_query_bases": sum(row["physical_bases"] for row in record_manifest),
        "native_validation": "pending symmetric offline scoring audit", "command": command,
        "projection": "offline modulo biological length before scientific scoring",
        "limitation": "doubled circular input may change LexicMap chaining; HSPs over one biological traversal are excluded and counted",
    }
    atomic_json(manifest_output, manifest)


if __name__ == "__main__":
    main()
