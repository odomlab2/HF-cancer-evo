#!/usr/bin/env bash
set -euo pipefail

# --- Args ---
SAMPLE=$1

# --- Constants ---
HF_SCC_ROOT="${HF_SCC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
HF_SCC_OUTPUT_ROOT="${HF_SCC_OUTPUT_ROOT:?Set HF_SCC_OUTPUT_ROOT to an absolute output directory}"
[[ "$HF_SCC_OUTPUT_ROOT" == /* ]] || { echo "ERROR: HF_SCC_OUTPUT_ROOT must be absolute" >&2; exit 2; }
BASE_DIR="${HF_SCC_OUTPUT_ROOT%/}/TES"
cd "$BASE_DIR"
OUTDIR="${BASE_DIR}/02_data/02_processed/callable_bases"
THREADS="${THREADS:-8}"
source "$HF_SCC_ROOT/TES/03_analysis/02_load_callable_contract.sh"
load_tes_callable_contract "$HF_SCC_ROOT" "$BASE_DIR"

DATA_ROOT="${HF_SCC_TES_DATA_ROOT:?Set HF_SCC_TES_DATA_ROOT}"
BAM_SUFFIX="${HF_SCC_TES_BAM_SUFFIX:?Set HF_SCC_TES_BAM_SUFFIX}"

# --- Dirs ---
ROWS="$OUTDIR/rows"
mkdir -p "$ROWS"

# --- Files ---
BAM="${DATA_ROOT}/${SAMPLE}/paired/merged-alignment/${SAMPLE}${BAM_SUFFIX}"
BAI="${BAM}.bai"
ROW_TSV="$ROWS/${SAMPLE}.tsv"

# --- Compute inclusive target size once per job ---
TARGET_SIZE=$(awk '{sum+=($3-$2+1)} END{print sum+0}' "$BED")

# --- Number of bases with coverage >= 10x ---
CALLABLE10_ON_TARGET=$(
  samtools depth -a -b "$BED" \
    -q "$MAPQ" \
    -Q "$BQ" \
    -s \
    -G 0xF04 \
    -J \
    "$BAM" |
  awk -v thr="$DEPTH" '{if($3>=thr) c10++} END{print (c10+0)}')

# --- On-target reads passing flags+MAPQ ---
ON_TARGET_READS=$(
  samtools view -@ "$THREADS" -c -F 0xF04 -f 0x2 -q "$MAPQ" -L "$BED" "$BAM")

# --- Overall callable bases >= 10x (not restricted to BED) ---
CALLABLE10_ALL=$(
  samtools depth \
    -q "$MAPQ" \
    -Q "$BQ" \
    -s \
    -G 0xF04 \
    -J \
    "$BAM" |
  awk -v thr="$DEPTH" '{if($3>=thr) c10++} END{print (c10+0)}')

# --- Overall reads MAPQ20+ (not restricted to BED) ---
ALL_READS_MAPQ20=$(
  samtools view -@ "$THREADS" -c \
    -F 0xF04 -f 0x2 -q "$MAPQ" \
    "$BAM")

# --- Percents ---
PCT10=$(awk -v c="$CALLABLE10_ON_TARGET" -v t="$TARGET_SIZE" 'BEGIN{if(t>0) printf "%.4f", 100.0*c/t; else print "NA"}')

# Write row
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
  "$SAMPLE" \
  "$TARGET_SIZE" \
  "$CALLABLE10_ON_TARGET" "$PCT10" \
  "$ON_TARGET_READS" \
  "$ALL_READS_MAPQ20" \
  "$CALLABLE10_ALL" \
  "$MAPQ" "$BQ" > "$ROW_TSV"

echo "[OK] ${SAMPLE}: onTarget10x=${CALLABLE10_ON_TARGET} (${PCT10}%%), onTargetReads=${ON_TARGET_READS}, allReadsMAPQ20=${ALL_READS_MAPQ20}, allCallable10x=${CALLABLE10_ALL}" >&2
