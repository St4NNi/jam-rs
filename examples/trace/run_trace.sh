#!/usr/bin/env bash
# Run one plasmid against one or more metagenomic assembly files.
#
# Usage:
#   bash examples/trace/run_trace.sh PLASMID_FASTA OUTPUT_DIR ASSEMBLY...
#
# Each assembly file is one metagenome sample in the candidate .jam index.
# The file may contain many contigs; do not pass --singleton when building the
# candidate index because singleton mode would make each contig a sample.

set -euo pipefail

if (( $# < 3 )); then
  printf 'usage: %s PLASMID_FASTA OUTPUT_DIR ASSEMBLY...\n' "$0" >&2
  exit 2
fi

JAM=${JAM:-jam}
THREADS=${JAM_THREADS:-1}
MEMORY_GB=${JAM_MEMORY_GB:-4}
SENSITIVITY=${JAM_TRACE_SENSITIVITY:-balanced}
PLASMID=$1
OUTPUT_DIR=$2
shift 2
ASSEMBLIES=("$@")

for input in "$PLASMID" "${ASSEMBLIES[@]}"; do
  if [[ ! -f "$input" ]]; then
    printf 'input is not a file: %s\n' "$input" >&2
    exit 1
  fi
done

if ! command -v "$JAM" >/dev/null 2>&1; then
  printf 'jam executable not found: %s\n' "$JAM" >&2
  exit 1
fi
if ! command -v jq >/dev/null 2>&1; then
  printf 'jq is required to extract the candidate report\n' >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR/archives"
DATABASE="$OUTPUT_DIR/metagenomes.jam"
CATALOG="$OUTPUT_DIR/metagenomes.tsv"
TRACE_OUTPUT=${TRACE_OUTPUT:-"$OUTPUT_DIR/trace.jsonl"}
REPORT="$OUTPUT_DIR/candidate_report.tsv"
RUN_MANIFEST="$OUTPUT_DIR/run_manifest.tsv"
CONFIRM_SCRIPT="$OUTPUT_DIR/confirm_with_mapping.sh"

case "$TRACE_OUTPUT" in
  *.zst|*.zstd)
    if command -v zstdcat >/dev/null 2>&1; then
      JSON_STREAM=(zstdcat)
    elif command -v unzstd >/dev/null 2>&1; then
      JSON_STREAM=(unzstd -c)
    else
      printf 'zstdcat or unzstd is required for compressed TRACE_OUTPUT: %s\n' \
        "$TRACE_OUTPUT" >&2
      exit 1
    fi
    ;;
  *)
    JSON_STREAM=(cat)
    ;;
esac

# Non-singleton mode gives one sample per input file. The resulting sample
# name is the input file basename, which is also the catalog metagenome_id.
"$JAM" --force --threads "$THREADS" --memory-target "$MEMORY_GB" sketch \
  "${ASSEMBLIES[@]}" \
  --output "$DATABASE" \
  --kmer-size 21 \
  --fscale 100

printf 'metagenome_id\tresource_uri\tsha256\n' > "$CATALOG"
declare -A seen_ids=()
for assembly in "${ASSEMBLIES[@]}"; do
  assembly_abs=$(realpath "$assembly")
  metagenome_id=$(basename "$assembly")
  if [[ -n "${seen_ids[$metagenome_id]+x}" ]]; then
    printf 'duplicate assembly basename, cannot create stable catalog ID: %s\n' \
      "$metagenome_id" >&2
    exit 1
  fi
  seen_ids["$metagenome_id"]=1
  archive="$OUTPUT_DIR/archives/$metagenome_id.jma"
  "$JAM" --force archive \
    --input "$assembly_abs" \
    --output "$archive" \
    --primary-scale 100 \
    --rescue-scale 100
  archive_sha256=$(sha256sum "$archive" | awk '{print $1}')
  printf '%s\t%s\t%s\n' \
    "$metagenome_id" \
    "archives/$metagenome_id.jma" \
    "$archive_sha256" >> "$CATALOG"
done

# The trace command records its own run header/footer and redacted input
# checksums. A plain .jsonl output keeps this example usable without a zstd
# command; set TRACE_OUTPUT to a .jsonl.zst path in a wrapper when desired.
"$JAM" --force --threads "$THREADS" --memory-target "$MEMORY_GB" trace \
  --plasmid "$PLASMID" \
  --database "$DATABASE" \
  --catalog "$CATALOG" \
  --output "$TRACE_OUTPUT" \
  --sensitivity "$SENSITIVITY" \
  --min-shared 3 \
  --min-plasmid-containment 0.0 \
  --min-metagenome-containment 0.0

# Extract stable candidate IDs and supporting contigs without loading a full
# result document into memory. These are evidence labels, not final calls:
# alert = sketch candidate, supported = positional alignment, confirmed is
# reserved for the independent mapping/context step below.
{
  printf 'metagenome_id\tevidence_label\tsupporting_contigs\tsupported_fraction\tcandidate_rank\n'
  "${JSON_STREAM[@]}" "$TRACE_OUTPUT" | jq -r '
    select(.record_type == "metagenome_result") |
    [
      .metagenome_id,
      (if .candidate == null then "no_candidate"
       elif (.alignments | length) == 0 then "alert"
       else "supported" end),
      ([.alignments[].contig_id] | unique | join(",")),
      (.coverage.supported_fraction // 0),
      (.candidate.rank // "")
    ] | @tsv
  '
} > "$REPORT"

"$JAM" stats --input "$DATABASE" --json > "$OUTPUT_DIR/database_stats.json"

{
  printf 'key\tvalue\n'
  printf 'jam_executable\t%s\n' "$JAM"
  printf 'threads\t%s\n' "$THREADS"
  printf 'memory_target_gib\t%s\n' "$MEMORY_GB"
  printf 'sensitivity\t%s\n' "$SENSITIVITY"
  printf 'candidate_database\t%s\n' "$DATABASE"
  printf 'catalog\t%s\n' "$CATALOG"
  printf 'trace_output\t%s\n' "$TRACE_OUTPUT"
  printf 'evidence_boundary\tjam-rs emits alert and supported; confirmed requires independent downstream evidence\n'
  sha256sum "$PLASMID" "$DATABASE" "$CATALOG" "$TRACE_OUTPUT"
} > "$RUN_MANIFEST"

# This file is intentionally a hand-off, not an automatic confirmation call.
# A reviewer supplies reads and decides whether an evidence row can be marked
# confirmed after inspecting mapping breadth/depth and relevant context.
cat > "$CONFIRM_SCRIPT" <<EOF
#!/usr/bin/env bash
set -euo pipefail

READS_R1=\${READS_R1:-reads_R1.fastq.gz}
READS_R2=\${READS_R2:-reads_R2.fastq.gz}
MAPPING_BAM=\"$OUTPUT_DIR/mapping.bam\"

# Independent confirmation example; jam-rs does not run this command.
minimap2 -ax sr -t \"\${THREADS:-8}\" \"$PLASMID\" \"\$READS_R1\" \"\$READS_R2\" | \\
  samtools sort -o \"\$MAPPING_BAM\"
samtools index \"\$MAPPING_BAM\"

# After reviewing mapping breadth/depth, marker consistency, and assembly
# context, write a separate report row with evidence_label=confirmed. Do not
# promote an alert or supported row to confirmed from sketch output alone.
printf 'evidence_label=confirmed is assigned by this downstream review only\\n'
EOF
chmod +x "$CONFIRM_SCRIPT"

printf 'Trace complete. Candidate report: %s\n' "$REPORT"
printf 'Downstream confirmation hand-off: %s\n' "$CONFIRM_SCRIPT"
