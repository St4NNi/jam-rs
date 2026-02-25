#!/usr/bin/env bash
set -euo pipefail

DEFAULT_OUTDIR="/vol/plasmidhunter/data/no_backup/hashes"
OUTDIR="${1:-${DEFAULT_OUTDIR}}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${HOME}/k21_pipeline"
SLURMOUT_DIR="${HOME}/slurmout"
ANALYSIS_DIR="${HOME}/k21_analysis"

SCRIPT_NAMES=(
    "k21_map.sbatch"
    "k21_prepare.sbatch"
    "k21_sort_reduce.sbatch"
    "k21_aggregate_hash.sbatch"
    "k21_aggregate_all.sbatch"
    "k21_progress.sh"
    "submit_k21_pipeline.sh"
)

mkdir -p "${OUTDIR}"
mkdir -p "${SLURMOUT_DIR}"
mkdir -p "${ANALYSIS_DIR}"

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

echo "[submit] output dir: ${OUTDIR}"
echo "[submit] slurm logs: ${SLURMOUT_DIR}"
echo "[submit] binary: ${BIN}"
echo "[submit] script dir: ${RUN_SCRIPT_DIR}"

echo "[submit] submitting map stage"
MAP_JOB_ID=$(sbatch --parsable "${RUN_SCRIPT_DIR}/k21_map.sbatch" "${OUTDIR}")

echo "[submit] submitting prepare stage after ${MAP_JOB_ID}"
PREPARE_JOB_ID=$(sbatch --parsable --dependency="afterok:${MAP_JOB_ID}" "${RUN_SCRIPT_DIR}/k21_prepare.sbatch" "${OUTDIR}")

echo "[submit] submitting sort-reduce stage after ${PREPARE_JOB_ID}"
SORT_JOB_ID=$(sbatch --parsable --dependency="afterok:${PREPARE_JOB_ID}" "${RUN_SCRIPT_DIR}/k21_sort_reduce.sbatch" "${OUTDIR}")

echo "[submit] submitting aggregate-hash stage after ${SORT_JOB_ID}"
AGG_HASH_JOB_ID=$(sbatch --parsable --dependency="afterok:${SORT_JOB_ID}" "${RUN_SCRIPT_DIR}/k21_aggregate_hash.sbatch" "${OUTDIR}")

echo "[submit] submitting aggregate-all stage after ${AGG_HASH_JOB_ID}"
AGG_ALL_JOB_ID=$(sbatch --parsable --dependency="afterok:${AGG_HASH_JOB_ID}" "${RUN_SCRIPT_DIR}/k21_aggregate_all.sbatch" "${OUTDIR}")

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
    printf 'SUBMITTED_AT="%s"\n' "$(date -Is)"
} >"${JOB_INFO_FILE}"

install -m 644 "${JOB_INFO_FILE}" "${RUN_SCRIPT_DIR}/last_jobs.env"

echo "[submit] done"
echo "  map:           ${MAP_JOB_ID}"
echo "  prepare:       ${PREPARE_JOB_ID}"
echo "  sort-reduce:   ${SORT_JOB_ID}"
echo "  aggregate-hash:${AGG_HASH_JOB_ID}"
echo "  aggregate-all: ${AGG_ALL_JOB_ID}"
echo "  job info file: ${JOB_INFO_FILE}"
echo "[submit] track progress with: ${RUN_SCRIPT_DIR}/k21_progress.sh ${OUTDIR}"
