#!/usr/bin/env bash
#SBATCH --job-name=k21-sort-reduce
#SBATCH --partition=bcf
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --array=0-255
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

DEFAULT_OUTDIR="/vol/plasmidhunter/data/no_backup/hashes"
OUTDIR="${1:-${DEFAULT_OUTDIR}}"
BIN="${HOME}/k21_pipeline"
SLURMOUT_DIR="${HOME}/slurmout"

mkdir -p "${SLURMOUT_DIR}"
JOB_TAG="${SLURM_JOB_NAME:-k21-sort-reduce}-${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-0}"
exec >"${SLURMOUT_DIR}/${JOB_TAG}.out" 2>"${SLURMOUT_DIR}/${JOB_TAG}.err"

if [[ ! -x "${BIN}" ]]; then
    echo "Missing executable: ${BIN}" >&2
    exit 1
fi

export RAYON_NUM_THREADS="${SLURM_CPUS_PER_TASK:-128}"
export K21_THREADS="${SLURM_CPUS_PER_TASK:-128}"
export MALLOC_CONF="${MALLOC_CONF:-background_thread:true,dirty_decay_ms:5000,muzzy_decay_ms:5000,metadata_thp:auto}"

/usr/bin/time -v "${BIN}" sort-reduce "${OUTDIR}"
