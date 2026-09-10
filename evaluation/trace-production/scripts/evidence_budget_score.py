#!/usr/bin/env python3
"""Score one frozen method and attribute every missed component's first loss."""

import argparse
import hashlib
import importlib.util
import json
import math
import subprocess
import sys
from pathlib import Path


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def load_jsonl(path: Path, key: str) -> dict[str, dict]:
    output = {}
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line:
            raise ValueError(f"blank JSONL line {number}: {path}")
        row = json.loads(line)
        identity = row[key]
        if identity in output:
            raise ValueError(f"duplicate {key}: {identity}")
        output[identity] = row
    return output


def load_scorer(path: Path):
    spec = importlib.util.spec_from_file_location("bound_profile_score", path)
    if spec is None or spec.loader is None:
        raise ValueError("cannot load base scorer")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def normalize_results(args, output: Path) -> Path:
    rows = load_jsonl(args.results, "query_id")
    membership = {}
    with args.membership.open(encoding="utf-8") as stream:
        header = stream.readline().rstrip("\n").split("\t")
        for line in stream:
            values = line.rstrip("\n").split("\t")
            row = dict(zip(header, values, strict=True))
            membership[row["name"]] = int(row["shard_ordinal"])
    changed = False
    for row in rows.values():
        if args.expected_manifest_sha256:
            index = row.get("index", {})
            expected = {"manifest_sha256": args.expected_manifest_sha256,
                        "body_sha256": args.expected_body_sha256,
                        "seed_k": args.expected_seed_k,
                        "rescue_k15": args.expected_rescue_k15}
            if any(index.get(field) != value for field, value in expected.items()):
                raise ValueError("single-index result identity differs")
        if "root_sha256" not in row:
            if not args.expected_manifest_sha256:
                raise ValueError("collection result lacks its serialized root identity")
            row["root_sha256"] = args.expected_root_sha256
            changed = True
        converted = []
        for item in row["metagenomes"]:
            if "trace" in item:
                converted.append(item)
            else:
                converted.append({"shard_ordinal": membership[item["name"]], "trace": item})
                changed = True
        row["metagenomes"] = converted
    if not changed:
        return args.results
    with output.open("x", encoding="utf-8") as stream:
        for row in rows.values():
            stream.write(json.dumps(row, sort_keys=True) + "\n")
    return output


def scoring_truth(args, output: Path) -> tuple[Path, str]:
    rows = load_jsonl(args.truth, "query_id")
    classes = {}
    for row in rows.values():
        classes[row["case"]] = classes.get(row["case"], 0) + 1
    if len(rows) == 42 and len(classes) == 21 and set(classes.values()) == {2}:
        return args.truth, "development"
    with output.open("x", encoding="utf-8") as stream:
        for row in rows.values():
            copied = dict(row)
            copied["split"] = "profile"
            stream.write(json.dumps(copied, sort_keys=True) + "\n")
    return output, "profile"


def score_command(args, results: Path, truth: Path, split: str, output: Path) -> list[str]:
    return [
        sys.executable, str(args.base_scorer), "--results", str(results),
        "--queries", str(args.queries), "--truth", str(truth),
        "--membership", str(args.membership), "--samtools", str(args.samtools),
        "--expected-root-sha256", args.expected_root_sha256,
        "--method", args.method, "--split", split,
        "--expected-queries", str(args.expected_queries),
        "--expected-members", str(args.expected_members), "--output", str(output),
    ]


def probe_rows(path: Path | None, truth: dict[str, dict]) -> dict[str, dict]:
    if path is None:
        return {}
    rows = load_jsonl(path, "query_id")
    if set(rows) != set(truth):
        raise ValueError("probe result query IDs differ from truth")
    for row in rows.values():
        for field in ("query_words_k21_k15", "present_keys_k21_k15",
                      "decoded_memberships_k21_k15", "decoded_positions_k21_k15"):
            if (not isinstance(row.get(field), list) or len(row[field]) != 2
                    or any(not isinstance(value, int) or value < 0 for value in row[field])):
                raise ValueError("invalid independent reader probe counters")
        if not isinstance(row.get("placements"), list):
            raise ValueError("reader probe placements are missing")
        if (not isinstance(row.get("dictionary_probes"), int) or row["dictionary_probes"] < 0
                or sum(row["present_keys_k21_k15"]) > row["dictionary_probes"]
                or row["dictionary_probes"] > sum(row["query_words_k21_k15"])):
            raise ValueError("reader probe dictionary counters differ")
        if (not isinstance(row.get("probe_seconds"), (int, float))
                or not math.isfinite(row["probe_seconds"]) or row["probe_seconds"] < 0):
            raise ValueError("reader probe timing is invalid")
        if len(row["placements"]) != sum(row["decoded_positions_k21_k15"]):
            raise ValueError("reader probe placement count differs from decoded positions")
        for placement in row["placements"]:
            if (placement.get("k") not in (15, 21) or not isinstance(placement.get("sample"), str)
                    or not isinstance(placement.get("contig"), str)
                    or not isinstance(placement.get("target_position"), int)
                    or placement["target_position"] < 0
                    or not isinstance(placement.get("target_canonical"), bool)
                    or not isinstance(placement.get("query_positions"), list)
                    or any(not isinstance(pair, list) or len(pair) != 2
                           or not isinstance(pair[0], int) or pair[0] < 0
                           or not isinstance(pair[1], bool)
                           for pair in placement.get("query_positions", []))):
                raise ValueError("reader probe placement is malformed")
        expected_circular = truth[row["query_id"]]["case"] == "origin_crossing"
        if row.get("circular") is not expected_circular:
            raise ValueError("reader probe query topology differs")
    return rows


def probe_totals(rows: dict[str, dict]) -> dict | None:
    if not rows:
        return None
    output = {"query_words_k21_k15": [0, 0], "present_keys_k21_k15": [0, 0],
              "decoded_memberships_k21_k15": [0, 0], "decoded_positions_k21_k15": [0, 0],
              "dictionary_probes": 0, "probe_seconds": 0.0}
    for row in rows.values():
        for field in ("query_words_k21_k15", "present_keys_k21_k15",
                      "decoded_memberships_k21_k15", "decoded_positions_k21_k15"):
            for scheme in range(2):
                output[field][scheme] += row[field][scheme]
        output["dictionary_probes"] += row["dictionary_probes"]
        output["probe_seconds"] += float(row["probe_seconds"])
    output["semantics"] = "exhaustive independent direct-reader audit, excluded from native timing"
    return output


def max_region_hits(pairs: set[tuple[int, int]], max_diagonal_drift: int = 64) -> int:
    regions = []
    for query, target in sorted(pairs):
        diagonal = target - query
        selected = None
        for region in reversed(regions):
            if (query > region["query_end"] and target > region["target_end"]
                    and min(region["diagonal_min"], diagonal) + max_diagonal_drift
                    >= max(region["diagonal_max"], diagonal)):
                selected = region
                break
        if selected is None:
            regions.append({"query_end": query, "target_end": target,
                            "diagonal_min": diagonal, "diagonal_max": diagonal, "hits": 1})
        else:
            selected["query_end"] = query
            selected["target_end"] = target
            selected["diagonal_min"] = min(selected["diagonal_min"], diagonal)
            selected["diagonal_max"] = max(selected["diagonal_max"], diagonal)
            selected["hits"] += 1
    return max((region["hits"] for region in regions), default=0)


def query_offset(position: int, k: int, segments: list[list[int]], query_length: int) -> int | None:
    selected = []
    for delta in range(k):
        query_position = (position + delta) % query_length
        offset = 0
        found = None
        for start, end in segments:
            if start <= query_position < end:
                found = offset + query_position - start
                break
            offset += end - start
        selected.append(found)
    if any(value is None for value in selected):
        return None
    return selected[0] if selected == list(range(selected[0], selected[0] + k)) else None


def probe_component(row: dict | None, item: dict, sample: str, query_length: int) -> dict | None:
    if row is None:
        return None
    pairs = {21: set(), 15: set()}
    reverse = item["strand"] == "reverse"
    target_start, target_end = item["target_interval"]
    for placement in row["placements"]:
        k = int(placement["k"])
        position = int(placement["target_position"])
        if (k not in pairs or placement["sample"] != sample or placement["contig"] != item["contig"]
                or position < target_start or position + k > target_end):
            continue
        for query_position, query_canonical in placement["query_positions"]:
            offset = query_offset(int(query_position), k, item["query_segments"], query_length)
            if offset is not None and (bool(placement["target_canonical"]) ^ bool(query_canonical)) == reverse:
                oriented_target = -position if reverse else position
                pairs[k].add((int(query_position), oriented_target))
    return {"primary_exact_pairs": len(pairs[21]), "rescue_exact_pairs": len(pairs[15]),
            "primary_max_product_region_hits": max_region_hits(pairs[21]),
            "rescue_max_product_region_hits": max_region_hits(pairs[15]),
            "product_diagonal_drift_bases": 64}


def classify_first_loss(evidence: dict | None) -> tuple[str, str]:
    if evidence is None:
        return "exact_evidence_audit_unavailable", "independent_reader_probe"
    primary, rescue = evidence["primary_exact_pairs"], evidence["rescue_exact_pairs"]
    primary_chain = evidence["primary_max_product_region_hits"]
    rescue_chain = evidence["rescue_max_product_region_hits"]
    if primary + rescue == 0:
        return "no_usable_exact_evidence", "independent_reader_probe"
    if primary_chain < 2 and rescue_chain < 3:
        return "insufficient_oriented_monotonic_evidence", "independent_reader_probe"
    return "after_exact_evidence_unknown_native_admission_or_alignment", "query_level_boundary"


def component_ledger(args, scorer, truth: dict[str, dict], results: dict[str, dict], probes: dict[str, dict]) -> list[dict]:
    queries = scorer.read_fasta(args.queries)
    by_name, by_source = scorer.membership(args.membership, args.expected_members)
    cache = {}
    output = []
    for query_id, expected in truth.items():
        query = queries[query_id]
        recovered = {index: [] for index in range(len(expected["components"]))}
        for wrapped in results[query_id]["metagenomes"]:
            trace = wrapped["trace"]
            member = by_name[trace["name"]]
            contigs = {int(item["id"]): item["name"] for item in trace["contigs"]}
            mosaic = trace["mosaic"]
            fragments = [entry["fragment"] for entry in mosaic["primary"]] + mosaic["alternatives"]
            for fragment in fragments:
                alignment = fragment["alignment"]
                start = int(alignment["target_interval"]["start"])
                end = int(alignment["target_interval"]["end"])
                contig = contigs[int(fragment["contig_id"])]
                target = scorer.target_slice(args.samtools, member, contig, start, end, cache)
                _, pairs = scorer.validate_fragment(fragment, query, target, start, end)
                for index, item in enumerate(expected["components"]):
                    if (by_source[item["source_id"]]["name"] == trace["name"]
                            and item["contig"] == contig and item["strand"] == alignment["strand"]):
                        recovered[index].extend(scorer.credited_support(pairs, tuple(item["target_interval"])))
        for index, item in enumerate(expected["components"]):
            available = [tuple(value) for value in item["available_query_intervals"]]
            available_bases = scorer.bases(available)
            recovered_bases = scorer.overlap(recovered[index], available)
            strict = recovered_bases == available_bases
            threshold = scorer.component_pass(recovered_bases, available_bases, expected["case"])
            evidence = probe_component(
                probes.get(query_id), item, by_source[item["source_id"]]["name"], len(query)
            )
            first_loss = None
            if 0 < recovered_bases < available_bases:
                first_loss = "direct_alignment_or_mosaic_selection"
                loss_scope = "reported_direct_alignment"
            elif not strict:
                first_loss, loss_scope = classify_first_loss(evidence)
            else:
                loss_scope = None
            output.append({"query_id": query_id, "case": expected["case"],
                           "component": index, "source_id": item["source_id"],
                           "contig": item["contig"], "target_interval": item["target_interval"],
                           "available_bases": available_bases, "recovered_bases": recovered_bases,
                           "lost_bases": available_bases - recovered_bases,
                           "strict_recovered": strict, "threshold_recovered": threshold,
                           "exact_evidence": evidence,
                           "first_loss_stage": first_loss,
                           "first_loss_scope": loss_scope,
                           "exact_evidence_limitation": (
                               "component-local audit; nearby anchors outside truth may route a wider alignment"
                               if first_loss in ("no_usable_exact_evidence",
                                                 "insufficient_oriented_monotonic_evidence") else None)})
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-scorer", type=Path, required=True)
    parser.add_argument("--expected-base-scorer-sha256", required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--membership", type=Path, required=True)
    parser.add_argument("--samtools", type=Path, required=True)
    parser.add_argument("--expected-samtools-sha256", required=True)
    parser.add_argument("--expected-root-sha256", required=True)
    parser.add_argument("--method", choices=("jam", "lexicmap"), required=True)
    parser.add_argument("--expected-queries", type=int, default=42)
    parser.add_argument("--expected-members", type=int, default=12)
    parser.add_argument("--expected-manifest-sha256")
    parser.add_argument("--expected-body-sha256")
    parser.add_argument("--expected-seed-k", type=int)
    parser.add_argument("--expected-rescue-k15", action=argparse.BooleanOptionalAction)
    parser.add_argument("--probe-results", type=Path)
    parser.add_argument("--probe-index", type=Path)
    parser.add_argument("--expected-probe-index-sha256")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (args.output.exists() or digest(args.base_scorer) != args.expected_base_scorer_sha256
            or not args.samtools.is_file() or args.samtools.is_symlink()
            or digest(args.samtools) != args.expected_samtools_sha256):
        raise SystemExit("new output and exact base scorer identity are required")
    index_bindings = (args.expected_manifest_sha256, args.expected_body_sha256,
                      args.expected_seed_k, args.expected_rescue_k15)
    if any(value is not None for value in index_bindings) != all(value is not None for value in index_bindings):
        raise SystemExit("all single-index result bindings must be supplied together")
    if (bool(args.probe_results) != bool(args.probe_index)
            or bool(args.probe_index) != bool(args.expected_probe_index_sha256)):
        raise SystemExit("probe results, index path, and index identity must be supplied together")
    if args.probe_index and (not args.probe_index.is_file() or args.probe_index.is_symlink()
                             or digest(args.probe_index) != args.expected_probe_index_sha256):
        raise SystemExit("probe index identity differs")
    args.output.mkdir(mode=0o700)
    normalized_results = normalize_results(args, args.output / "normalized-results.jsonl")
    bound_truth, base_split = scoring_truth(args, args.output / "scoring-truth.jsonl")
    base_output = args.output / "score.json"
    command = score_command(args, normalized_results, bound_truth, base_split, base_output)
    run = subprocess.run(command, check=False, capture_output=True, text=True)
    (args.output / "score.stdout").write_text(run.stdout, encoding="utf-8")
    (args.output / "score.stderr").write_text(run.stderr, encoding="utf-8")
    if run.returncode or not base_output.is_file():
        raise SystemExit("base scorer did not produce an assessment")
    base = json.loads(base_output.read_text(encoding="utf-8"))
    if base.get("status") != "evaluated":
        assessment = {"status": "unevaluable", "base_score": base,
                      "base_scorer_sha256": digest(args.base_scorer), "command": command}
    else:
        truth = load_jsonl(args.truth, "query_id")
        results = load_jsonl(normalized_results, "query_id")
        probes = probe_rows(args.probe_results, truth)
        scorer = load_scorer(args.base_scorer)
        components = component_ledger(args, scorer, truth, results, probes)
        first_losses = {}
        for row in components:
            stage = row["first_loss_stage"]
            if stage:
                first_losses[stage] = first_losses.get(stage, 0) + 1
        assessment = {"status": "evaluated", "method": args.method,
                      "base_scorer_sha256": digest(args.base_scorer),
                      "results_sha256": digest(args.results),
                      "normalized_results_sha256": digest(normalized_results),
                      "normalization_applied": normalized_results != args.results,
                      "queries_sha256": digest(args.queries),
                      "truth_sha256": digest(args.truth),
                      "scoring_truth_sha256": digest(bound_truth),
                      "base_scorer_split_adapter": base_split if bound_truth != args.truth else None,
                      "probe_results_sha256": digest(args.probe_results) if args.probe_results else None,
                      "probe_index_sha256": digest(args.probe_index) if args.probe_index else None,
                      "probe_totals": probe_totals(probes),
                      "samtools_sha256": digest(args.samtools),
                      "identity_score_and_cigar_validity": "passed_by_bound_base_scorer",
                      "metrics": base, "component_ledger": components,
                      "first_loss_components": first_losses,
                      "first_loss_attribution_status": "available" if probes else "probe_unavailable"}
    (args.output / "assessment.json").write_text(
        json.dumps(assessment, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
