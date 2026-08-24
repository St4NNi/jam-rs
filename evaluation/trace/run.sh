#!/usr/bin/env bash
set -euo pipefail

# Trace-stage benchmark. Build the release binary before invoking this
# script; compilation is deliberately outside all timed commands.
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PYTHON="${PYTHON:-python3}"
JAM="${JAM:-$ROOT/target/release/jam}"
PROFILE="smoke"
DATASET_ROOT=""
QUERY_ID=""
WORK=""
REPEATS=3
THREADS=8
WITH_REMOTE=1
WITH_COMPARATORS=0

usage() {
    cat <<'EOF'
Usage: evaluation/trace/run.sh [options]

Options:
  --profile NAME       synthetic input size: smoke, standard, or large
  --dataset-root PATH  import an existing materialized dataset; no download
  --query-id ID        query record when --dataset-root has a multi-record catalog
  --work-dir PATH      new workspace output directory
  --repeats N          warm repeats per profile/thread (default: 3)
  --threads N          intended production thread count (default: 8)
  --no-remote          skip the local range-capable HTTP fixture
  --comparators        measure optional locally installed BLASTn and minimap2
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --profile)
            PROFILE="$2"
            shift 2
            ;;
        --dataset-root)
            DATASET_ROOT="$2"
            shift 2
            ;;
        --query-id)
            QUERY_ID="$2"
            shift 2
            ;;
        --work-dir)
            WORK="$2"
            shift 2
            ;;
        --repeats)
            REPEATS="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --no-remote)
            WITH_REMOTE=0
            shift
            ;;
        --comparators)
            WITH_COMPARATORS=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            usage >&2
            exit 2
            ;;
    esac
done

if [[ ! -x "$JAM" ]]; then
    echo "release binary not found or not executable: $JAM" >&2
    echo "build it before timing, for example: cargo build --release --locked" >&2
    exit 1
fi
if [[ "$REPEATS" -lt 1 || "$THREADS" -lt 1 ]]; then
    echo "--repeats and --threads must be positive integers" >&2
    exit 2
fi
if [[ -n "$DATASET_ROOT" && ! -d "$DATASET_ROOT" ]]; then
    echo "--dataset-root is not a directory: $DATASET_ROOT" >&2
    exit 2
fi
if [[ -z "$WORK" ]]; then
    WORK="$ROOT/tools/out/trace-benchmark/$(date -u +%Y%m%dT%H%M%SZ)"
fi
if [[ -e "$WORK" ]]; then
    echo "refusing to reuse existing benchmark work directory: $WORK" >&2
    exit 1
fi
mkdir -p "$WORK" "$WORK/indexes" "$WORK/outputs" "$WORK/measurements" "$WORK/logs"

GEN_ARGS=(--output "$WORK/dataset" --profile "$PROFILE")
if [[ -n "$DATASET_ROOT" ]]; then
    GEN_ARGS+=(--dataset-root "$DATASET_ROOT")
fi
if [[ -n "$QUERY_ID" ]]; then
    GEN_ARGS+=(--query-id "$QUERY_ID")
fi
"$PYTHON" "$ROOT/evaluation/trace/generate.py" "${GEN_ARGS[@]}" > "$WORK/generation.json"
QUERY_ID="$("$PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1], encoding="utf-8"))["query_id"])' "$WORK/generation.json")"
mapfile -t ASSEMBLIES < <(find "$WORK/dataset/assemblies" -maxdepth 1 -type f \( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' -o -name '*.fastq' -o -name '*.fq' \) -print | sort)
if [[ "${#ASSEMBLIES[@]}" -eq 0 ]]; then
    echo "dataset has no assembly files" >&2
    exit 1
fi

measure_command() {
    local name="$1"
    local stage="$2"
    local label="$3"
    local profile="$4"
    local cache_state="$5"
    local threads="$6"
    local trace_output="$7"
    local server_url="$8"
    shift 8
    local args=(
        --output "$WORK/measurements/${name}.json"
        --stdout "$WORK/logs/${name}.stdout"
        --stderr "$WORK/logs/${name}.stderr"
        --stage "$stage"
        --label "$label"
        --profile "$profile"
        --cache-state "$cache_state"
        --threads "$threads"
    )
    if [[ -n "$trace_output" ]]; then
        args+=(--trace-output "$trace_output" --truth "$WORK/dataset/truth.tsv")
    fi
    if [[ -n "$server_url" ]]; then
        args+=(--server-url "$server_url")
    fi
    "$PYTHON" "$ROOT/evaluation/trace/measure.py" "${args[@]}" -- "$@"
}

"$PYTHON" "$ROOT/evaluation/trace/environment.py" \
    --output "$WORK/environment.json" --jam "$JAM" --dataset "$WORK/dataset"

measure_command \
    jam_index catalog_index jam-sketch none none "$THREADS" "" "" \
    "$JAM" --threads "$THREADS" --memory-target 4 --silent sketch "${ASSEMBLIES[@]}" \
    --output "$WORK/indexes/metagenomes.jam" --kmer-size 31 --fscale 200

measure_command \
    jma_index catalog_index jma-archive none none "$THREADS" "" "" \
    "$PYTHON" "$ROOT/evaluation/trace/build_archives.py" --jam "$JAM" \
    --assemblies-dir "$WORK/dataset/assemblies" --output-dir "$WORK/indexes/jma" \
    --primary-scale 100 --rescue-scale 100 --threads "$THREADS" --memory-target 4

"$PYTHON" "$ROOT/evaluation/trace/make_catalog.py" \
    --assemblies-dir "$WORK/dataset/assemblies" --jma-dir "$WORK/indexes/jma" \
    --output "$WORK/catalog.local.tsv" --mode local

run_trace() {
    local name="$1"
    local profile="$2"
    local threads="$3"
    local cache_state="$4"
    local catalog="$5"
    local server_url="$6"
    local output="$WORK/outputs/${name}.jsonl"
    local command=(
        "$JAM" --threads "$threads" --memory-target 4 --silent trace
        --plasmid "$WORK/dataset/query.fasta"
        --database "$WORK/indexes/metagenomes.jam"
        --catalog "$catalog"
        --output "$output"
        --plasmid-id "$QUERY_ID"
        --sensitivity "$profile"
        --min-shared 3
        --min-plasmid-containment 0
        --min-metagenome-containment 0
        --top-candidates 250
        --max-alignments 256
    )
    measure_command "$name" trace_search jam_trace "$profile" "$cache_state" "$threads" "$output" "$server_url" "${command[@]}"
}

run_local_profile() {
    local profile="$1"
    local threads="$2"
    local stem="local_${profile}_t${threads}"
    # process-cold does not mean kernel-cache-cold; portable cache eviction is
    # privileged and intentionally not attempted by this harness.
    run_trace "${stem}_cold" "$profile" "$threads" process-cold "$WORK/catalog.local.tsv" ""
    local warmup="$WORK/outputs/${stem}_warmup.jsonl"
    "$JAM" --threads "$threads" --memory-target 4 --silent trace \
        --plasmid "$WORK/dataset/query.fasta" --database "$WORK/indexes/metagenomes.jam" \
        --catalog "$WORK/catalog.local.tsv" --output "$warmup" --plasmid-id "$QUERY_ID" \
        --sensitivity "$profile" --min-shared 3 --top-candidates 250 --max-alignments 256
    for repeat in $(seq 1 "$REPEATS"); do
        run_trace "${stem}_warm_r${repeat}" "$profile" "$threads" warm "$WORK/catalog.local.tsv" ""
    done
}

THREAD_VALUES=(1)
if [[ "$THREADS" -ne 1 ]]; then
    THREAD_VALUES+=("$THREADS")
fi
for profile in fast balanced sensitive; do
    for threads in "${THREAD_VALUES[@]}"; do
        run_local_profile "$profile" "$threads"
    done
done

if [[ "$WITH_REMOTE" -eq 1 ]]; then
    HTTP_PID=""
    cleanup_http() {
        if [[ -n "$HTTP_PID" ]]; then
            kill "$HTTP_PID" 2>/dev/null || true
            wait "$HTTP_PID" 2>/dev/null || true
        fi
    }
    trap cleanup_http EXIT
    "$PYTHON" "$ROOT/evaluation/trace/http_server.py" \
        --root "$WORK/indexes/jma" --metrics "$WORK/http.metrics.json" \
        > "$WORK/http.port" 2> "$WORK/http.stderr" &
    HTTP_PID=$!
    for _ in $(seq 1 100); do
        if grep -q '^PORT=' "$WORK/http.port"; then
            break
        fi
        sleep 0.1
    done
    if ! grep -q '^PORT=' "$WORK/http.port"; then
        echo "mock HTTP server did not publish a port" >&2
        exit 1
    fi
    PORT="$(sed -n 's/^PORT=//p' "$WORK/http.port" | head -n 1)"
    BASE_URL="http://127.0.0.1:${PORT}"
    "$PYTHON" "$ROOT/evaluation/trace/make_catalog.py" \
        --assemblies-dir "$WORK/dataset/assemblies" --jma-dir "$WORK/indexes/jma" \
        --output "$WORK/catalog.http.tsv" --mode http --base-url "$BASE_URL"
    # Balanced exercises both primary and rescue seeds; local runs cover all
    # three sensitivity profiles.
    run_trace remote_balanced_t1_cold balanced 1 process-cold "$WORK/catalog.http.tsv" "$BASE_URL/__metrics"
    remote_warmup="$WORK/outputs/remote_balanced_warmup.jsonl"
    "$JAM" --threads "$THREADS" --memory-target 4 --silent trace \
        --plasmid "$WORK/dataset/query.fasta" --database "$WORK/indexes/metagenomes.jam" \
        --catalog "$WORK/catalog.http.tsv" --output "$remote_warmup" --plasmid-id "$QUERY_ID" \
        --sensitivity balanced --min-shared 3 --top-candidates 250 --max-alignments 256
    for repeat in $(seq 1 "$REPEATS"); do
        run_trace "remote_balanced_t${THREADS}_warm_r${repeat}" balanced "$THREADS" warm \
            "$WORK/catalog.http.tsv" "$BASE_URL/__metrics"
    done
    cleanup_http
    HTTP_PID=""
    trap - EXIT
fi

if [[ "$WITH_COMPARATORS" -eq 1 ]]; then
    for tool in blastn minimap2; do
        measure_command "comparator_${tool}" comparator "$tool" none none 1 "" "" \
            "$PYTHON" "$ROOT/evaluation/trace/comparators.py" --tool "$tool" \
            --query "$WORK/dataset/query.fasta" --assemblies-dir "$WORK/dataset/assemblies" \
            --output "$WORK/outputs/comparator_${tool}.json"
    done
fi

"$PYTHON" "$ROOT/evaluation/trace/summarize.py" \
    --measurements "$WORK/measurements" --output "$WORK/summary.tsv"
echo "Trace benchmark complete: $WORK"
echo "Summary: $WORK/summary.tsv"
