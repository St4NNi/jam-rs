#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 IMAGE [OUTPUT_DIRECTORY]" >&2
    exit 2
fi

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMAGE="$1"
OUTPUT="${2:-$ROOT/tools/out/container-validation-$(date -u +%Y%m%dT%H%M%SZ)}"

if [[ -e "$OUTPUT" ]]; then
    echo "Refusing to reuse existing container-validation directory: $OUTPUT" >&2
    exit 1
fi
mkdir -p "$OUTPUT/tmp"

container_uid="$(docker run --rm --entrypoint /usr/bin/id "$IMAGE" -u)"
if [[ "$container_uid" -eq 0 ]]; then
    echo "Container image runs as root" >&2
    exit 1
fi

docker run --rm "$IMAGE" --version
docker run --rm --entrypoint /usr/bin/curl "$IMAGE" --version >/dev/null
docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env TMPDIR=/output/tmp \
    --volume "$ROOT:/workspace:ro" \
    --volume "$OUTPUT:/output" \
    "$IMAGE" sketch /workspace/tests/testfiles/short.fa \
    --output /output/catalog.jam \
    --singleton \
    --kmer-size 21 \
    --fscale 1 \
    --threads 2 \
    --memory 1 \
    --silent

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env TMPDIR=/output/tmp \
    --volume "$ROOT:/workspace:ro" \
    --volume "$OUTPUT:/output" \
    "$IMAGE" sketch /workspace/tests/testfiles/short.fa \
    --output /output/metagenomes.jam \
    --kmer-size 21 \
    --fscale 1 \
    --threads 2 \
    --memory 1 \
    --silent
docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env TMPDIR=/output/tmp \
    --volume "$ROOT:/workspace:ro" \
    --volume "$OUTPUT:/output" \
    "$IMAGE" archive \
    --input /workspace/tests/testfiles/short.fa \
    --output /output/short.jma \
    --primary-scale 1 \
    --rescue-scale 1 \
    --block-bases 256 \
    --silent
printf 'metagenome_id\tjma\tjma_index\traw\nshort.fa\t/output/short.jma\t/output/short.jma.idx.json\t/workspace/tests/testfiles/short.fa\n' > "$OUTPUT/metagenomes.tsv"
docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env TMPDIR=/output/tmp \
    --volume "$ROOT:/workspace:ro" \
    --volume "$OUTPUT:/output" \
    "$IMAGE" trace \
    --query /workspace/tests/testfiles/short.fa \
    --query-kind other \
    --topology linear \
    --database /output/metagenomes.jam \
    --metagenomes /output/metagenomes.tsv \
    --output /output/trace.jsonl.zst \
    --sensitivity sensitive \
    --min-shared 1 \
    --top-candidates 1 \
    --threads 2 \
    --memory 1 \
    --silent

docker run --rm \
    --user "$(id -u):$(id -g)" \
    --env TMPDIR=/output/tmp \
    --volume "$ROOT:/workspace:ro" \
    --volume "$OUTPUT:/output" \
    "$IMAGE" screen \
    --input /workspace/tests/testfiles/short.fa \
    --database /output/catalog.jam \
    --output /output/contigs.tsv \
    --summary /output/summary.tsv \
    --metadata /output/run.json \
    --min-shared 1 \
    --threads 2 \
    --memory 1 \
    --silent

grep -q '^schema_version.*query_contig.*query_containment.*reference_containment' "$OUTPUT/contigs.tsv"
grep -q '^schema_version.*shared_hashes_union.*aggregate_reference_containment' "$OUTPUT/summary.tsv"
test -s "$OUTPUT/catalog.jam.json"
test -s "$OUTPUT/run.json"
test -s "$OUTPUT/trace.jsonl.zst"

echo "Container validation passed: $OUTPUT"
