#!/usr/bin/env python3
"""Build one frozen, query-independent density variant with the installed JIDX path."""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def save(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--targets-sha256", required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--binary-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--window", type=int, choices=[16, 32, 64], required=True)
    args = parser.parse_args()
    job = os.environ.get("SLURM_JOB_ID", "")
    if not job.isdigit() or os.environ.get("SLURM_CPUS_PER_TASK") != "4":
        raise SystemExit("four CPU Slurm job required")
    scratch = Path("/var/scratch/sbeyvers") / f"jam-evidence-build-{job}"
    stage = scratch / "candidate"
    if scratch.exists() or args.output.exists():
        raise SystemExit("output or owned scratch exists")
    if args.output.parent != Path("/vol/plasmidhunter/data/tests/jampub/v0.10/18_generic_trace/.claude/20260910/evidence-budget"):
        raise SystemExit("output boundary")
    if digest(args.binary) != args.binary_sha256:
        raise SystemExit("binary identity")
    frozen_path = args.targets / "FROZEN.json"
    if digest(frozen_path) != args.targets_sha256:
        raise SystemExit("frozen target identity")
    frozen = json.loads(frozen_path.read_text())
    if frozen.get("status") != "frozen" or frozen.get("query_dependent_selection") is not False:
        raise SystemExit("target selection is not frozen and query independent")
    manifest = args.targets / "jidx-manifest.json"
    if digest(manifest) != frozen["jidx_manifest_sha256"]:
        raise SystemExit("target manifest identity")
    if digest(args.targets / "lexicmap-inputs.txt") != frozen["lexicmap_inputs_sha256"]:
        raise SystemExit("target FASTA list identity")
    targets = json.loads(manifest.read_text())["metagenomes"]
    if len(targets) != 12 or len({t["name"] for t in targets}) != 12:
        raise SystemExit("twelve distinct target members required")
    sources = [Path(line) for line in (args.targets / "lexicmap-inputs.txt").read_text().splitlines()]
    if len(sources) != 12:
        raise SystemExit("twelve whole-contig FASTA sources required")
    for target in frozen["targets"]:
        for kind in ("bgzf", "fai", "gzi", "fasta"):
            if digest(target[kind]) != target["sha256"][kind]:
                raise SystemExit(f"target {kind} identity")
    scratch.mkdir(mode=0o700)
    stage.mkdir(mode=0o700)
    os.environ.update(TMPDIR=str(scratch), TEMP=str(scratch), TMP=str(scratch))
    if shutil.disk_usage(scratch).free < 4 * 1024**3:
        raise SystemExit("scratch reserve below 4 GiB")
    aliases = []
    for target, source in zip(frozen["targets"], sources, strict=True):
        if Path(target["name"]).name != target["name"] or not target["name"]:
            raise SystemExit("sample name is not a filename")
        alias = scratch / target["name"]
        alias.symlink_to(source)
        aliases.append(alias)
    commands = [
        [str(args.binary), "--threads", "4", "--memory", "4", "--silent", "sketch",
         *map(str, aliases), "--output", str(stage / "database.jam"), "--kmer-size", "21",
         "--fscale", "100", "--complexity", "0", "--temp-dir", str(scratch)],
        [str(args.binary), "--threads", "4", "--silent", "jidx",
         "--database", str(stage / "database.jam"), "--manifest", str(manifest),
         "--output", str(stage / "index.jidx"), "--kmer-size", "21",
         "--minimizer-window", str(args.window), "--rescue-k15"],
    ]
    try:
        for ordinal, command in enumerate(commands):
            with (stage / f"build-{ordinal}.stdout").open("xb") as stdout, (stage / f"build-{ordinal}.stderr").open("xb") as stderr:
                subprocess.run(["/usr/bin/time", "-v", "-o", str(stage / f"build-{ordinal}.time"), *command],
                               check=True, timeout=280, stdout=stdout, stderr=stderr)
        index = stage / "index.jidx"
        with index.open("rb") as stream:
            header = stream.read(512)
            body = hashlib.file_digest(stream, "sha256").digest()
        if header[112:144] != body:
            raise ValueError("JIDX body digest")
        files = {name: {"bytes": (stage / name).stat().st_size, "sha256": digest(stage / name)}
                 for name in ("database.jam", "index.jidx")}
        sequence_files = {str(Path(t["bgzf"]).resolve()) for t in targets}
        sequence_bytes = sum(Path(path).stat().st_size for path in sequence_files)
        complete = sum(v["bytes"] for v in files.values()) + sequence_bytes + manifest.stat().st_size
        save(stage / "COMPLETE.json", {
            "status": "complete", "window": args.window, "k": 21, "rescue_k15": True,
            "target_frozen_sha256": digest(frozen_path), "target_manifest_sha256": digest(manifest),
            "binary": str(args.binary), "binary_sha256": args.binary_sha256, "commands": commands,
            "job_id": job, "hostname": os.uname().nodename, "cpu_affinity": sorted(os.sched_getaffinity(0)),
            "files": files, "bgzf_bytes_once": sequence_bytes,
            "query_ready_bytes": complete, "manifest_bytes": manifest.stat().st_size,
            "accounting": "JAM plus JIDX plus BGZF once plus runtime manifest; GZI and contig metadata embedded in JIDX",
            "excluded_build_evaluation_inputs": ["uncompressed FASTA", "external FAI and GZI", "FROZEN.json", "build logs"],
            "targets": frozen_path.as_posix(), "frozen_status": frozen.get("status"),
        })
        publication = args.output.with_name(args.output.name + f".publish-{job}")
        publication.mkdir(mode=0o700)
        for source in stage.iterdir():
            destination = publication / source.name
            with source.open("rb") as input_stream, destination.open("xb") as output_stream:
                shutil.copyfileobj(input_stream, output_stream, 1024 * 1024)
                output_stream.flush()
                os.fsync(output_stream.fileno())
            if digest(source) != digest(destination):
                raise ValueError("published copy differs")
        directory = os.open(publication, os.O_RDONLY | os.O_DIRECTORY)
        os.fsync(directory)
        os.close(directory)
        os.rename(publication, args.output)
        parent = os.open(args.output.parent, os.O_RDONLY | os.O_DIRECTORY)
        os.fsync(parent)
        os.close(parent)
    except Exception as error:
        failure = args.output.with_name(args.output.name + f".failed-{job}")
        failure.mkdir(mode=0o700)
        for source in stage.iterdir():
            if source.suffix in (".stdout", ".stderr", ".time"):
                shutil.copyfile(source, failure / source.name)
        save(failure / "FAILURE.json", {"status": "failed", "job_id": job,
             "exception": repr(error), "commands": commands,
             "target_frozen_sha256": args.targets_sha256,
             "binary_sha256": args.binary_sha256})
        raise
    finally:
        if scratch.parent == Path("/var/scratch/sbeyvers") and scratch.name == f"jam-evidence-build-{job}" and not scratch.is_symlink():
            shutil.rmtree(scratch)


if __name__ == "__main__":
    main()
