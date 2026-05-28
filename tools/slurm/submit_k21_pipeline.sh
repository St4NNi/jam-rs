#!/usr/bin/env bash
set -euo pipefail

DEFAULT_OUTDIR="/vol/plasmidhunter/data/no_backup/hashes"
OUTDIR="${1:-${DEFAULT_OUTDIR}}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${HOME}/k21_pipeline"
SLURMOUT_DIR="${HOME}/slurmout"
ANALYSIS_DIR="${HOME}/k21_analysis"

MAP_TASKS=64
PARTITIONS=64
SORT_TASKS=$((4 * PARTITIONS))
HASHES=("jamhash" "murmur3" "xxhash3" "wang64")

SCRIPT_NAMES=(
    "k21_map.sbatch"
    "k21_prepare.sbatch"
    "k21_sort_reduce.sbatch"
    "k21_aggregate_hash.sbatch"
    "k21_aggregate_all.sbatch"
    "k21_linear_test.sbatch"
    "k21_linear_aggregate.sbatch"
    "k21_progress.sh"
    "submit_k21_pipeline.sh"
)

mkdir -p "${OUTDIR}" "${SLURMOUT_DIR}" "${ANALYSIS_DIR}"

for script_name in "${SCRIPT_NAMES[@]}"; do
    src="${SCRIPT_DIR}/${script_name}"
    dst="${ANALYSIS_DIR}/${script_name}"
    if [[ "${src}" != "${dst}" ]]; then
        install -m 755 "${src}" "${dst}"
    fi
done

if [[ "${SCRIPT_DIR}/README.md" != "${ANALYSIS_DIR}/README.md" ]]; then
    install -m 644 "${SCRIPT_DIR}/README.md" "${ANALYSIS_DIR}/README.md"
fi

if [[ ! -x "${BIN}" ]]; then
    echo "Missing executable: ${BIN}" >&2
    echo "Please place the compiled k21_pipeline binary at ${BIN}" >&2
    exit 1
fi

RUN_SCRIPT_DIR="${ANALYSIS_DIR}"

count_glob() {
    local pattern="$1"
    local files=()
    mapfile -t files < <(compgen -G "${pattern}" || true)
    printf '%s' "${#files[@]}"
}

apply_array_limit() {
    local spec="$1"
    local max_running="$2"
    if [[ -z "${spec}" || -z "${max_running}" ]]; then
        printf '%s' "${spec}"
        return
    fi
    if [[ "${spec}" == *%* ]]; then
        printf '%s' "${spec}"
        return
    fi
    printf '%s%%%s' "${spec}" "${max_running}"
}

map_is_complete() {
    local n
    n=$(count_glob "${OUTDIR}/_meta/map/map_task_*.json")
    [[ "${n}" -eq "${MAP_TASKS}" ]]
}

prepare_is_complete() {
    [[ -f "${OUTDIR}/_meta/partition_plan.json" ]]
}

missing_sort_array_spec() {
    local missing=()
    local task hash_idx partition hash_name stats_path
    for ((task = 0; task < SORT_TASKS; task++)); do
        hash_idx=$((task / PARTITIONS))
        partition=$((task % PARTITIONS))
        hash_name="${HASHES[hash_idx]}"
        printf -v stats_path '%s/%s/partition_stats/part_%02d.json' "${OUTDIR}" "${hash_name}" "${partition}"
        if [[ ! -f "${stats_path}" ]]; then
            missing+=("${task}")
        fi
    done
    if (( ${#missing[@]} == 0 )); then
        printf ''
    else
        local joined
        joined=$(IFS=, ; echo "${missing[*]}")
        printf '%s' "${joined}"
    fi
}

missing_agg_hash_array_spec() {
    local missing=()
    local i hash_name metrics
    for ((i = 0; i < 4; i++)); do
        hash_name="${HASHES[i]}"
        metrics="${OUTDIR}/${hash_name}/metrics.json"
        if [[ ! -f "${metrics}" ]]; then
            missing+=("${i}")
        fi
    done
    if (( ${#missing[@]} == 0 )); then
        printf ''
    else
        local joined
        joined=$(IFS=, ; echo "${missing[*]}")
        printf '%s' "${joined}"
    fi
}

build_dep_args() {
    local dep_job_id="$1"
    local -n _out_ref="$2"
    _out_ref=()
    if [[ -n "${dep_job_id}" ]]; then
        _out_ref+=(--kill-on-invalid-dep=yes --dependency="afterok:${dep_job_id}")
    fi
}

FORCE_MAP="${K21_FORCE_MAP:-0}"
FORCE_PREPARE="${K21_FORCE_PREPARE:-0}"
FORCE_SORT="${K21_FORCE_SORT:-0}"
FORCE_AGG_HASH="${K21_FORCE_AGG_HASH:-0}"
FORCE_AGG_ALL="${K21_FORCE_AGG_ALL:-0}"
FORCE_LINEAR="${K21_FORCE_LINEAR:-0}"

MAP_JOB_ID="SKIPPED"
PREPARE_JOB_ID="SKIPPED"
SORT_JOB_ID="SKIPPED"
AGG_HASH_JOB_ID="SKIPPED"
AGG_ALL_JOB_ID="SKIPPED"

MAP_ARRAY_SPEC="${K21_MAP_ARRAY_SPEC:-0-63}"
MAP_ARRAY_SPEC=$(apply_array_limit "${MAP_ARRAY_SPEC}" "${K21_MAP_ARRAY_MAX_RUNNING:-}")

SORT_MEM_ARGS=()
if [[ -n "${K21_SORT_MEM:-}" ]]; then
    SORT_MEM_ARGS+=(--mem="${K21_SORT_MEM}")
fi

SORT_ARRAY_SPEC=""
if [[ -n "${K21_SORT_ARRAY_SPEC:-}" ]]; then
    SORT_ARRAY_SPEC="${K21_SORT_ARRAY_SPEC}"
fi

echo "[submit] output dir: ${OUTDIR}"
echo "[submit] slurm logs: ${SLURMOUT_DIR}"
echo "[submit] binary: ${BIN}"
echo "[submit] script dir: ${RUN_SCRIPT_DIR}"

map_submitted=0
if [[ "${FORCE_MAP}" == "1" ]]; then
    echo "[submit] forcing map stage"
elif map_is_complete; then
    echo "[submit] map stage already complete, skipping"
else
    echo "[submit] map stage incomplete, submitting"
fi

if [[ "${FORCE_MAP}" == "1" || ! map_is_complete ]]; then
    MAP_JOB_ID=$(sbatch --parsable --array="${MAP_ARRAY_SPEC}" "${RUN_SCRIPT_DIR}/k21_map.sbatch" "${OUTDIR}")
    map_submitted=1
    echo "[submit] map job: ${MAP_JOB_ID} (array=${MAP_ARRAY_SPEC})"
fi

prepare_submitted=0
prepare_needed=0
if (( map_submitted == 1 )); then
    prepare_needed=1
elif [[ "${FORCE_PREPARE}" == "1" || ! prepare_is_complete ]]; then
    prepare_needed=1
fi

if (( prepare_needed == 1 )); then
    PREPARE_DEP_ARGS=()
    if (( map_submitted == 1 )); then
        build_dep_args "${MAP_JOB_ID}" PREPARE_DEP_ARGS
    fi
    PREPARE_JOB_ID=$(sbatch --parsable "${PREPARE_DEP_ARGS[@]}" "${RUN_SCRIPT_DIR}/k21_prepare.sbatch" "${OUTDIR}")
    prepare_submitted=1
    echo "[submit] prepare job: ${PREPARE_JOB_ID}"
else
    echo "[submit] prepare stage already complete, skipping"
fi

sort_submitted=0
if [[ "${FORCE_SORT}" == "1" ]]; then
    if [[ -z "${SORT_ARRAY_SPEC}" ]]; then
        SORT_ARRAY_SPEC="0-255"
    fi
elif (( map_submitted == 1 )); then
    if [[ -z "${SORT_ARRAY_SPEC}" ]]; then
        SORT_ARRAY_SPEC="0-255"
    fi
else
    if [[ -z "${SORT_ARRAY_SPEC}" ]]; then
        SORT_ARRAY_SPEC="$(missing_sort_array_spec)"
    fi
fi

SORT_ARRAY_SPEC=$(apply_array_limit "${SORT_ARRAY_SPEC}" "${K21_SORT_ARRAY_MAX_RUNNING:-}")

if [[ -n "${SORT_ARRAY_SPEC}" ]]; then
    SORT_DEP_ARGS=()
    if (( prepare_submitted == 1 )); then
        build_dep_args "${PREPARE_JOB_ID}" SORT_DEP_ARGS
    fi
    SORT_JOB_ID=$(sbatch --parsable "${SORT_DEP_ARGS[@]}" "${SORT_MEM_ARGS[@]}" --array="${SORT_ARRAY_SPEC}" "${RUN_SCRIPT_DIR}/k21_sort_reduce.sbatch" "${OUTDIR}")
    sort_submitted=1
    echo "[submit] sort-reduce job: ${SORT_JOB_ID} (array=${SORT_ARRAY_SPEC})"
else
    echo "[submit] sort-reduce already complete, skipping"
fi

agg_hash_submitted=0
AGG_HASH_ARRAY_SPEC=""
if [[ "${FORCE_AGG_HASH}" == "1" ]]; then
    AGG_HASH_ARRAY_SPEC="0-3"
elif (( sort_submitted == 1 )); then
    AGG_HASH_ARRAY_SPEC="0-3"
else
    AGG_HASH_ARRAY_SPEC="$(missing_agg_hash_array_spec)"
fi

AGG_HASH_ARRAY_SPEC=$(apply_array_limit "${AGG_HASH_ARRAY_SPEC}" "${K21_AGG_HASH_ARRAY_MAX_RUNNING:-}")

if [[ -n "${AGG_HASH_ARRAY_SPEC}" ]]; then
    AGG_HASH_DEP_ARGS=()
    if (( sort_submitted == 1 )); then
        build_dep_args "${SORT_JOB_ID}" AGG_HASH_DEP_ARGS
    fi
    AGG_HASH_JOB_ID=$(sbatch --parsable "${AGG_HASH_DEP_ARGS[@]}" --array="${AGG_HASH_ARRAY_SPEC}" "${RUN_SCRIPT_DIR}/k21_aggregate_hash.sbatch" "${OUTDIR}")
    agg_hash_submitted=1
    echo "[submit] aggregate-hash job: ${AGG_HASH_JOB_ID} (array=${AGG_HASH_ARRAY_SPEC})"
else
    echo "[submit] aggregate-hash already complete, skipping"
fi

agg_all_needed=0
if [[ "${FORCE_AGG_ALL}" == "1" ]]; then
    agg_all_needed=1
elif (( agg_hash_submitted == 1 )); then
    agg_all_needed=1
elif [[ ! -f "${OUTDIR}/combined_metrics.json" ]]; then
    agg_all_needed=1
fi

if (( agg_all_needed == 1 )); then
    AGG_ALL_DEP_ARGS=()
    if (( agg_hash_submitted == 1 )); then
        build_dep_args "${AGG_HASH_JOB_ID}" AGG_ALL_DEP_ARGS
    fi
    AGG_ALL_JOB_ID=$(sbatch --parsable "${AGG_ALL_DEP_ARGS[@]}" "${RUN_SCRIPT_DIR}/k21_aggregate_all.sbatch" "${OUTDIR}")
    echo "[submit] aggregate-all job: ${AGG_ALL_JOB_ID}"
else
    echo "[submit] aggregate-all already complete, skipping"
fi

LINEAR_MAP_JOB_ID="SKIPPED"
LINEAR_AGG_JOB_ID="SKIPPED"

linear_map_needed=0
if [[ "${FORCE_LINEAR}" == "1" ]]; then
    linear_map_needed=1
else
    linear_map_missing=0
    for ((lt = 0; lt < 32; lt++)); do
        printf -v lt_path '%s/linear_test/linear_task_%02d.json' "${OUTDIR}" "${lt}"
        if [[ ! -f "${lt_path}" ]]; then
            linear_map_missing=1
            break
        fi
    done
    if (( linear_map_missing == 1 )); then
        linear_map_needed=1
    fi
fi

if (( linear_map_needed == 1 )); then
    LINEAR_MAP_JOB_ID=$(sbatch --parsable "${RUN_SCRIPT_DIR}/k21_linear_test.sbatch" "${OUTDIR}")
    echo "[submit] linear-map job: ${LINEAR_MAP_JOB_ID} (array=0-31, independent)"
else
    echo "[submit] linear-map already complete, skipping"
fi

linear_agg_needed=0
if [[ "${FORCE_LINEAR}" == "1" ]]; then
    linear_agg_needed=1
elif (( linear_map_needed == 1 )); then
    linear_agg_needed=1
elif [[ ! -f "${OUTDIR}/linear_test/linear_test.json" ]]; then
    linear_agg_needed=1
fi

if (( linear_agg_needed == 1 )); then
    LINEAR_AGG_DEP_ARGS=()
    if [[ "${LINEAR_MAP_JOB_ID}" != "SKIPPED" ]]; then
        build_dep_args "${LINEAR_MAP_JOB_ID}" LINEAR_AGG_DEP_ARGS
    fi
    LINEAR_AGG_JOB_ID=$(sbatch --parsable "${LINEAR_AGG_DEP_ARGS[@]}" "${RUN_SCRIPT_DIR}/k21_linear_aggregate.sbatch" "${OUTDIR}")
    echo "[submit] linear-aggregate job: ${LINEAR_AGG_JOB_ID}"
else
    echo "[submit] linear-aggregate already complete, skipping"
fi

META_DIR="${OUTDIR}/_meta"
JOB_INFO_FILE="${META_DIR}/slurm_jobs.env"
mkdir -p "${META_DIR}"

{
    printf 'OUTDIR="%s"\n' "${OUTDIR}"
    printf 'MAP_JOB_ID="%s"\n' "${MAP_JOB_ID}"
    printf 'PREPARE_JOB_ID="%s"\n' "${PREPARE_JOB_ID}"
    printf 'SORT_JOB_ID="%s"\n' "${SORT_JOB_ID}"
    printf 'AGG_HASH_JOB_ID="%s"\n' "${AGG_HASH_JOB_ID}"
    printf 'AGG_ALL_JOB_ID="%s"\n' "${AGG_ALL_JOB_ID}"
    printf 'LINEAR_MAP_JOB_ID="%s"\n' "${LINEAR_MAP_JOB_ID}"
    printf 'LINEAR_AGG_JOB_ID="%s"\n' "${LINEAR_AGG_JOB_ID}"
    printf 'MAP_ARRAY_SPEC="%s"\n' "${MAP_ARRAY_SPEC}"
    printf 'SORT_ARRAY_SPEC="%s"\n' "${SORT_ARRAY_SPEC}"
    printf 'AGG_HASH_ARRAY_SPEC="%s"\n' "${AGG_HASH_ARRAY_SPEC}"
    printf 'SUBMITTED_AT="%s"\n' "$(date -Is)"
} >"${JOB_INFO_FILE}"

install -m 644 "${JOB_INFO_FILE}" "${RUN_SCRIPT_DIR}/last_jobs.env"

echo "[submit] done"
echo "  map:            ${MAP_JOB_ID}"
echo "  prepare:        ${PREPARE_JOB_ID}"
echo "  sort-reduce:    ${SORT_JOB_ID}"
echo "  aggregate-hash: ${AGG_HASH_JOB_ID}"
echo "  aggregate-all:  ${AGG_ALL_JOB_ID}"
echo "  linear-map:     ${LINEAR_MAP_JOB_ID}"
echo "  linear-agg:     ${LINEAR_AGG_JOB_ID}"
echo "  job info file:  ${JOB_INFO_FILE}"
echo "[submit] track progress with: ${RUN_SCRIPT_DIR}/k21_progress.sh ${OUTDIR}"
