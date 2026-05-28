#!/usr/bin/env bash
#SBATCH --job-name=k21-aggregate-all
#SBATCH --partition=bcf
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

DEFAULT_OUTDIR="/vol/plasmidhunter/data/no_backup/hashes"
OUTDIR="${1:-${DEFAULT_OUTDIR}}"
BIN="${HOME}/k21_pipeline"
SLURMOUT_DIR="${HOME}/slurmout"

mkdir -p "${SLURMOUT_DIR}"
JOB_TAG="${SLURM_JOB_NAME:-k21-aggregate-all}-${SLURM_JOB_ID:-manual}"
LOG_FILE="${SLURMOUT_DIR}/${JOB_TAG}.log"
exec >"${LOG_FILE}" 2>&1

if [[ ! -x "${BIN}" ]]; then
    echo "Missing executable: ${BIN}" >&2
    exit 1
fi

export MALLOC_CONF="${MALLOC_CONF:-background_thread:true,dirty_decay_ms:5000,muzzy_decay_ms:5000,metadata_thp:auto}"

/usr/bin/time -v "${BIN}" aggregate-all "${OUTDIR}"
