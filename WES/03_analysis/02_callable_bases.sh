#!/usr/bin/env bash
set -euo pipefail

# --- Args ---
SAMPLE=$1

# --- Constants ---
HF_SCC_ROOT="${HF_SCC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
HF_SCC_OUTPUT_ROOT="${HF_SCC_OUTPUT_ROOT:?Set HF_SCC_OUTPUT_ROOT to an absolute output directory}"
[[ "$HF_SCC_OUTPUT_ROOT" == /* ]] || { echo "ERROR: HF_SCC_OUTPUT_ROOT must be absolute" >&2; exit 2; }
BASE_DIR="${HF_SCC_OUTPUT_ROOT%/}/WES"
cd "$BASE_DIR"

OUTDIR="${BASE_DIR}/02_data/02_processed/callable_bases"
THREADS="${THREADS:-8}"
RUN_CONTRACT="${OUTDIR}/run_contract.tsv"

SEQ_DIR="${HF_SCC_WES_DATA_ROOT:?Set HF_SCC_WES_DATA_ROOT}"
BAM_PROJECT_PREFIX="${HF_SCC_WES_BAM_PROJECT_PREFIX:?Set HF_SCC_WES_BAM_PROJECT_PREFIX}"

# Thresholds are declared once per fresh run by the submission stage.
expected_header=$'assay\tmin_depth\tmapq\tbase_quality\texpected_sample_count\tcallable_column'
[[ -r $RUN_CONTRACT ]] || {
  echo "ERROR: missing WES callable run contract" >&2
  exit 6
}
mapfile -t contract_lines < "$RUN_CONTRACT"
[[ ${#contract_lines[@]} -eq 2 && ${contract_lines[0]} == "$expected_header" ]] || {
  echo "ERROR: malformed WES callable run contract" >&2
  exit 6
}
IFS=$'\t' read -r ASSAY DEPTH MAPQ BQ EXPECTED_COUNT CALLABLE_COLUMN \
  <<< "${contract_lines[1]}"
[[ $ASSAY == WES && $DEPTH =~ ^[0-9]+$ && $MAPQ =~ ^[0-9]+$ &&
   $BQ =~ ^[0-9]+$ && $EXPECTED_COUNT =~ ^[0-9]+$ &&
   -n $CALLABLE_COLUMN ]] || {
  echo "ERROR: malformed WES callable run contract" >&2
  exit 6
}

# --- Dirs ---
ROWS="$OUTDIR/rows"
mkdir -p "$ROWS"

# --- HELPERS ---

# Resolve sample BAM by wildcard search across all PIDs
resolve_sample_bam() {
  local sample="$1"
  local seq_dir="$2"
  local -a cands=()

  shopt -s nullglob
  cands=( "${seq_dir}/"*/"${sample}/paired/merged-alignment/${sample}_${BAM_PROJECT_PREFIX}"*_merged.mdup.bam )
  shopt -u nullglob

  if (( ${#cands[@]} != 1 )); then
    echo "ERROR: BAM resolution failed for sample '$sample'" >&2
    echo "  searched across all PIDs under: $seq_dir" >&2
    echo "  matches found: ${#cands[@]}" >&2
    if (( ${#cands[@]} > 0 )); then
      printf '  %s\n' "${cands[@]}" >&2
    fi
    return 5
  fi

  printf '%s\n' "${cands[0]}"
}

# --- Files ---
BAM="$(resolve_sample_bam "$SAMPLE" "$SEQ_DIR")"
BAI="${BAM}.bai"
ROW_TSV="$ROWS/${SAMPLE}.tsv"

# --- Overall callable bases >= 10x ---
CALLABLE10_ALL=$(
  samtools depth \
    -q "$MAPQ" \
    -Q "$BQ" \
    -s \
    -G 0xF04 \
    -J \
    "$BAM" |
  awk -v thr="$DEPTH" '{if($3>=thr) c10++} END{print (c10+0)}'
)

# --- Overall reads MAPQ20+ (not restricted to BED) ---
ALL_READS_MAPQ20=$(
  samtools view -@ "$THREADS" -c \
    -F 0xF04 -f 0x2 -q "$MAPQ" \
    "$BAM"
)

# Write row
printf "%s\t%s\t%s\t%s\t%s\n" \
  "$SAMPLE" \
  "$ALL_READS_MAPQ20" \
  "$CALLABLE10_ALL" \
  "$MAPQ" "$BQ" > "$ROW_TSV"

echo "[OK] ${SAMPLE}: allReadsMAPQ20=${ALL_READS_MAPQ20}, allCallable10x=${CALLABLE10_ALL}" >&2
