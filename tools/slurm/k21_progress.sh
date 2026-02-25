#!/usr/bin/env bash
set -euo pipefail

DEFAULT_OUTDIR="/vol/plasmidhunter/data/no_backup/hashes"
OUTDIR="${1:-${DEFAULT_OUTDIR}}"
INTERVAL_SEC="${2:-5}"
JOB_INFO_FILE="${OUTDIR}/_meta/slurm_jobs.env"

HASHES=("jamhash" "murmur3" "xxhash3" "wang64")
MAP_TOTAL=64
SORT_TOTAL=256

if [[ -f "${JOB_INFO_FILE}" ]]; then
    # shellcheck source=/dev/null
    source "${JOB_INFO_FILE}"
fi

count_glob() {
    local pattern="$1"
    local files=()
    mapfile -t files < <(compgen -G "${pattern}" || true)
    printf '%s' "${#files[@]}"
}

progress_line() {
    local label="$1"
    local done="$2"
    local total="$3"
    local pct=0
    if (( total > 0 )); then
        pct=$(( done * 100 / total ))
    fi
    printf '%-18s %5d / %-5d (%3d%%)\n' "${label}" "${done}" "${total}" "${pct}"
}

squeue_state_summary() {
    if ! command -v squeue >/dev/null 2>&1; then
        return
    fi

    local job_ids=()
    local id
    for id in "${MAP_JOB_ID:-}" "${PREPARE_JOB_ID:-}" "${SORT_JOB_ID:-}" "${AGG_HASH_JOB_ID:-}" "${AGG_ALL_JOB_ID:-}"; do
        if [[ -n "${id}" ]]; then
            job_ids+=("${id}")
        fi
    done

    if (( ${#job_ids[@]} == 0 )); then
        return
    fi

    local ids_csv
    ids_csv=$(IFS=, ; echo "${job_ids[*]}")
    local queue
    queue=$(squeue -h -j "${ids_csv}" -o "%i %T" 2>/dev/null || true)

    local running=0
    local pending=0
    local completing=0
    local other=0
    local _jid state

    while read -r _jid state; do
        if [[ -z "${state:-}" ]]; then
            continue
        fi
        case "${state}" in
            RUNNING) running=$((running + 1)) ;;
            PENDING) pending=$((pending + 1)) ;;
            COMPLETING) completing=$((completing + 1)) ;;
            *) other=$((other + 1)) ;;
        esac
    done <<< "${queue}"

    printf 'Slurm queue         R:%d P:%d C:%d O:%d\n' "${running}" "${pending}" "${completing}" "${other}"
}

print_snapshot() {
    local map_done
    local prepare_done
    local sort_done=0
    local agg_hash_done=0
    local agg_all_done
    local hash

    map_done=$(count_glob "${OUTDIR}/_meta/map/map_task_*.json")

    if [[ -f "${OUTDIR}/_meta/partition_plan.json" ]]; then
        prepare_done=1
    else
        prepare_done=0
    fi

    for hash in "${HASHES[@]}"; do
        sort_done=$((sort_done + $(count_glob "${OUTDIR}/${hash}/partition_stats/part_*.json")))
        if [[ -f "${OUTDIR}/${hash}/metrics.json" ]]; then
            agg_hash_done=$((agg_hash_done + 1))
        fi
    done

    if [[ -f "${OUTDIR}/combined_metrics.json" ]]; then
        agg_all_done=1
    else
        agg_all_done=0
    fi

    printf '\n[%s] outdir=%s\n' "$(date '+%F %T')" "${OUTDIR}"
    progress_line "map tasks" "${map_done}" "${MAP_TOTAL}"
    progress_line "prepare" "${prepare_done}" "1"
    progress_line "sort-reduce" "${sort_done}" "${SORT_TOTAL}"
    progress_line "aggregate-hash" "${agg_hash_done}" "4"
    progress_line "aggregate-all" "${agg_all_done}" "1"
    squeue_state_summary
    printf 'Logs                %s\n' "${HOME}/slurmout"
}

if [[ "${ONCE:-0}" == "1" ]]; then
    print_snapshot
    exit 0
fi

while true; do
    clear
    print_snapshot
    sleep "${INTERVAL_SEC}"
done
