#!/usr/bin/env bash
set -euo pipefail

# --- CONFIG ---
HF_SCC_ROOT="${HF_SCC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
BASE_DIR="${HF_SCC_OUTPUT_ROOT:?Set HF_SCC_OUTPUT_ROOT to an absolute output directory}/TES"
cd "$BASE_DIR"
SAMPLE=$1 # iterate through the samples
REF="${HF_SCC_REFERENCE_FASTA:?Set HF_SCC_REFERENCE_FASTA}"
SEQ_DIR="${HF_SCC_TES_DATA_ROOT:?Set HF_SCC_TES_DATA_ROOT}"
TUMOUR_BAM_SUFFIX="${HF_SCC_TES_BAM_SUFFIX:?Set HF_SCC_TES_BAM_SUFFIX}"
THREADS="${THREADS:-24}"
NORMAL_MANIFEST="${HF_SCC_ROOT}/TES/01_source/tes_matched_normals.tsv"
source "${HF_SCC_ROOT}/TES/01_source/tes_matched_normal_resolver.sh"
NORMAL_ID=$(tes_resolve_matched_normal "$NORMAL_MANIFEST" "$SAMPLE")

# Resolve the same logical private inputs before creating any output.
TUMOUR_BAM="${SEQ_DIR}/${SAMPLE}/paired/merged-alignment/${SAMPLE}${TUMOUR_BAM_SUFFIX}"
GERM_DIR="${SEQ_DIR}/${NORMAL_ID}"
[[ -d "$GERM_DIR" ]] || {
  printf "ERROR: matched-normal directory is missing for '%s'\n" "$NORMAL_ID" >&2
  exit 3
}
GERM_BAM=$(find -L "$GERM_DIR" -type f -name '*_merged.mdup.bam' -print -quit)
[[ -n "$GERM_BAM" ]] || {
  printf "ERROR: no matched-normal BAM resolves for '%s'\n" "$NORMAL_ID" >&2
  exit 3
}
OUT_DIR="${BASE_DIR}/02_data/02_processed/mutation_calling"
mkdir -p "$OUT_DIR"
RUN_INFO_TXT="${OUT_DIR}/mutation_calling_run_info.txt"
RUN_INFO_LOCK="${RUN_INFO_TXT}.lock"
TMP_ECHO="$(mktemp -p "$OUT_DIR" ".${SAMPLE}.XXXX.echo")"

# --- ENV ---
: "${GATK_BIN:=gatk}"

# --- HELPERS ---

# Extract the unique sample name from a BAM header
bam_sm() {
  local bam="$1"
  mapfile -t names < <(
    samtools view -H "$bam" \
      | awk -F'\t' '$1=="@RG"{for(i=1;i<=NF;i++) if($i~/^SM:/){sub(/^SM:/,"",$i); print $i}}' \
      | sort -u
  )
  (( ${#names[@]} == 1 )) || return 2
  printf '%s\n' "${names[0]}"
}

# --- LOGGING ---
log() { printf '%s\n' "$*" | tee -a "$TMP_ECHO" >/dev/null; }
append_run_info() {
  local status=$?
  {
    flock -x 200
    {
      printf '=== %s | %s | status=%s ===\n' "$SAMPLE" "$(date -Is)" "$status"
      cat "$TMP_ECHO"
      printf '\n'
    } >> "$RUN_INFO_TXT"
  } 200>>"$RUN_INFO_LOCK"
  rm -f "$TMP_ECHO"
}
trap append_run_info EXIT

# --- PATHS ---
GERM_SM="$(bam_sm "$GERM_BAM")"

UNFILTERED="${OUT_DIR}/${SAMPLE}.unfiltered.vcf.gz"
FILTERED="${OUT_DIR}/${SAMPLE}.filtered.vcf.gz"

log "Tumour: $TUMOUR_BAM (SM=$SAMPLE)"
log "Normal: $GERM_BAM  (SM=$GERM_SM)"
log "Ref: $REF"
log "Threads: $THREADS"

# --- MUTECT2 ---
"$GATK_BIN" Mutect2 \
    -R "$REF" \
    -I "$TUMOUR_BAM" \
    -I "$GERM_BAM" \
    -normal "$GERM_SM" \
    --native-pair-hmm-threads "$THREADS" \
    -O "$UNFILTERED"

# --- FILTER ---
"$GATK_BIN" FilterMutectCalls \
    -R "$REF" \
    -V "$UNFILTERED" \
    -O "$FILTERED"
