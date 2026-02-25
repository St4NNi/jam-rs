#!/usr/bin/env bash
#SBATCH --job-name=k21-map
#SBATCH --partition=bcf
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --array=0-63
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

DEFAULT_OUTDIR="/vol/plasmidhunter/data/no_backup/hashes"
OUTDIR="${1:-${DEFAULT_OUTDIR}}"
BIN="${HOME}/k21_pipeline"
SLURMOUT_DIR="${HOME}/slurmout"

mkdir -p "${SLURMOUT_DIR}"
JOB_TAG="${SLURM_JOB_NAME:-k21-map}-${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-0}"
exec >"${SLURMOUT_DIR}/${JOB_TAG}.out" 2>"${SLURMOUT_DIR}/${JOB_TAG}.err"

if [[ ! -x "${BIN}" ]]; then
    echo "Missing executable: ${BIN}" >&2
    exit 1
fi

export RAYON_NUM_THREADS="${SLURM_CPUS_PER_TASK:-128}"
export K21_THREADS="${SLURM_CPUS_PER_TASK:-128}"
export K21_IO_THREADS="${K21_IO_THREADS:-16}"
export K21_LOCAL_BUFFER_VALUES="${K21_LOCAL_BUFFER_VALUES:-32768}"
export K21_MAP_FILE_BUFFER_MB="${K21_MAP_FILE_BUFFER_MB:-64}"
export K21_IO_CHANNEL_ENTRIES="${K21_IO_CHANNEL_ENTRIES:-4000000}"
export K21_LOCAL_SCRATCH_BASE="${K21_LOCAL_SCRATCH_BASE:-/var/scratch/sbeyvers/k21}"
export MALLOC_CONF="${MALLOC_CONF:-background_thread:true,dirty_decay_ms:5000,muzzy_decay_ms:5000,metadata_thp:auto}"

/usr/bin/time -v "${BIN}" map "${OUTDIR}"
