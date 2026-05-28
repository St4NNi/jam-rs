#!/usr/bin/env bash
#SBATCH --job-name=k21-map
#SBATCH --partition=bcf
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --array=0-63
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -eEuo pipefail

cancel_array_on_error() {
    local rc=$?
    if [[ ${rc} -ne 0 && -n "${SLURM_ARRAY_JOB_ID:-}" ]]; then
        echo "[fail-fast] task failure, cancelling array ${SLURM_ARRAY_JOB_ID}" >&2
        scancel "${SLURM_ARRAY_JOB_ID}" >/dev/null 2>&1 || true
    fi
    exit "${rc}"
}

trap cancel_array_on_error ERR

DEFAULT_OUTDIR="/vol/plasmidhunter/data/no_backup/hashes"
OUTDIR="${1:-${DEFAULT_OUTDIR}}"
BIN="${HOME}/k21_pipeline"
SLURMOUT_DIR="${HOME}/slurmout"

mkdir -p "${SLURMOUT_DIR}"
JOB_TAG="${SLURM_JOB_NAME:-k21-map}-${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-0}"
LOG_FILE="${SLURMOUT_DIR}/${JOB_TAG}.log"
exec >"${LOG_FILE}" 2>&1

if [[ ! -x "${BIN}" ]]; then
    echo "Missing executable: ${BIN}" >&2
    exit 1
fi

if ! command -v cpb >/dev/null 2>&1; then
    echo "Missing command: cpb (required for local->/vol copy)" >&2
    exit 1
fi

if ! command -v flock >/dev/null 2>&1; then
    echo "Missing command: flock (required for per-node map lock)" >&2
    exit 1
fi

export RAYON_NUM_THREADS="${SLURM_CPUS_PER_TASK:-128}"
export K21_THREADS="${SLURM_CPUS_PER_TASK:-128}"
export K21_IO_THREADS="${K21_IO_THREADS:-16}"
export K21_LOCAL_BUFFER_VALUES="${K21_LOCAL_BUFFER_VALUES:-32768}"
export K21_MAP_FILE_BUFFER_MB="${K21_MAP_FILE_BUFFER_MB:-64}"
export K21_IO_CHANNEL_ENTRIES="${K21_IO_CHANNEL_ENTRIES:-4000000}"
export K21_IO_MAX_OPEN_WRITERS="${K21_IO_MAX_OPEN_WRITERS:-8}"
export K21_LOCAL_SCRATCH_BASE="${K21_LOCAL_SCRATCH_BASE:-/var/scratch/sbeyvers/k21}"
export K21_EXPECTED_TASK_SCRATCH_GIB="${K21_EXPECTED_TASK_SCRATCH_GIB:-1024}"
export K21_SCRATCH_HEADROOM_GIB="${K21_SCRATCH_HEADROOM_GIB:-128}"
export K21_ENFORCE_MIN_SCRATCH="${K21_ENFORCE_MIN_SCRATCH:-1}"
export MALLOC_CONF="${MALLOC_CONF:-background_thread:true,dirty_decay_ms:5000,muzzy_decay_ms:5000,metadata_thp:auto}"

mkdir -p "${K21_LOCAL_SCRATCH_BASE}" || true

node_name="$(hostname -s 2>/dev/null || hostname)"
node_lock_file="${K21_LOCAL_SCRATCH_BASE}/.map-node-${node_name}.lock"
exec 9>"${node_lock_file}"
if [[ "${K21_NODE_LOCK_TIMEOUT_SEC:-0}" == "0" ]]; then
    flock 9
else
    flock -w "${K21_NODE_LOCK_TIMEOUT_SEC}" 9
fi

avail_bytes=$(df -PB1 "${K21_LOCAL_SCRATCH_BASE}" 2>/dev/null | awk 'NR==2 {print $4}')
if [[ -n "${avail_bytes}" && "${avail_bytes}" =~ ^[0-9]+$ ]]; then
    avail_gib=$((avail_bytes / 1024 / 1024 / 1024))
    required_gib=$((K21_EXPECTED_TASK_SCRATCH_GIB + K21_SCRATCH_HEADROOM_GIB))
    echo "[meta] host=${node_name} task=${SLURM_ARRAY_TASK_ID:-na} scratch_free_gib=${avail_gib} required_gib=${required_gib}"
    if (( avail_gib < required_gib )); then
        echo "[fail] insufficient scratch on $(hostname): ${avail_gib} GiB free < ${required_gib} GiB required" >&2
        if [[ "${K21_ENFORCE_MIN_SCRATCH}" == "1" ]]; then
            df -h "${K21_LOCAL_SCRATCH_BASE}" >&2 || true
            df -i "${K21_LOCAL_SCRATCH_BASE}" >&2 || true
            exit 1
        else
            echo "[warn] continuing because K21_ENFORCE_MIN_SCRATCH=${K21_ENFORCE_MIN_SCRATCH}" >&2
        fi
    fi
else
    echo "[meta] host=$(hostname) task=${SLURM_ARRAY_TASK_ID:-na} scratch_check=unavailable" >&2
fi

/usr/bin/time -v "${BIN}" map "${OUTDIR}"
