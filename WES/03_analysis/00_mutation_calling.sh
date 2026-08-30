#!/usr/bin/env bash
set -euo pipefail

# --- CONFIG ---
HF_SCC_ROOT="${HF_SCC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
BASE_DIR="${HF_SCC_OUTPUT_ROOT:?Set HF_SCC_OUTPUT_ROOT to an absolute output directory}/WES"
cd "$BASE_DIR"
SAMPLE=$1 # iterate through the samples
REF="${HF_SCC_REFERENCE_FASTA:?Set HF_SCC_REFERENCE_FASTA}"
SEQ_DIR="${HF_SCC_WES_DATA_ROOT:?Set HF_SCC_WES_DATA_ROOT}"
BAM_PROJECT_PREFIX="${HF_SCC_WES_BAM_PROJECT_PREFIX:?Set HF_SCC_WES_BAM_PROJECT_PREFIX}"
THREADS="${THREADS:-24}"
OUT_DIR="${BASE_DIR}/02_data/02_processed/mutation_calling"
mkdir -p "$OUT_DIR"
RUN_INFO_TXT="${OUT_DIR}/mutation_calling_run_info.txt"
RUN_INFO_LOCK="${RUN_INFO_TXT}.lock"
TMP_ECHO="$(mktemp -p "$OUT_DIR" ".${SAMPLE}.XXXX.echo")"

# --- ENV ---
: "${GATK_BIN:=gatk}"

# --- HELPERS ---

# Get mouse number from sample name
extract_mouse_num() {
  local s="$(basename "$1")"
  if [[ "$s" =~ ^(hairfollicle|skin)([0-9]+) ]]; then
    printf '%s\n' "${BASH_REMATCH[2]}"
  elif [[ "$s" =~ -w([0-9]+) ]]; then
    printf '%s\n' "${BASH_REMATCH[1]}"
  else
    echo "ERROR: cannot parse mouse number from '$s'" >&2
    return 2
  fi
}

# Extract unique SM from BAM header
bam_sm() {
  local bam="$1"
  local -a names=()
  mapfile -t names < <(
    samtools view -H "$bam" \
      | awk -F'\t' '$1=="@RG"{
          for(i=1;i<=NF;i++){
            if($i~/^SM:/){sub(/^SM:/,"",$i); print $i}
          }
        }' \
      | sort -u
  )
  (( ${#names[@]} == 1 )) || {
    echo "ERROR: expected exactly 1 SM in BAM header, got ${#names[@]} for $bam" >&2
    return 2
  }
  printf '%s\n' "${names[0]}"
}

# Resolve tumour BAM by wildcard search across all PIDs
resolve_tumour_bam() {
  local sample="$1"
  local SEQ_DIR="$2"
  local -a cands=()

  shopt -s nullglob
  cands=( "${SEQ_DIR}/"*/"${sample}/paired/merged-alignment/${sample}_${BAM_PROJECT_PREFIX}"*_merged.mdup.bam )
  shopt -u nullglob

  if (( ${#cands[@]} != 1 )); then
    echo "ERROR: tumour BAM resolution failed for sample '$sample'" >&2
    echo "  searched across all PIDs under: $SEQ_DIR" >&2
    echo "  matches found: ${#cands[@]}" >&2
    if (( ${#cands[@]} > 0 )); then
      printf '  %s\n' "${cands[@]}" >&2
    fi
    return 5
  fi

  printf '%s\n' "${cands[0]}"
}

# Choose germline directory
choose_germline_bam_dir() {
  local num="$1"
  local seq_root="$2"
  local -a dirs=()

  shopt -s nullglob
  dirs=(
    "${seq_root}/"*/"liver${num}"
    "${seq_root}/"*/"liver${num}-w"*
  )
  shopt -u nullglob

  if (( ${#dirs[@]} != 1 )); then
    echo "ERROR: expected exactly 1 liver directory for mouse $num across all PIDs; found ${#dirs[@]}" >&2
    if (( ${#dirs[@]} > 0 )); then
      printf '  %s\n' "${dirs[@]}" >&2
    fi
    return 3
  fi

  printf '%s\n' "${dirs[0]}"
}

# Resolve unique merged.mdup.bam in a sample directory
resolve_unique_merged_bam_in_sample_dir() {
  local sample_dir="$1"
  local -a bams=()

  shopt -s nullglob
  bams=( "${sample_dir}/paired/merged-alignment/"*_merged.mdup.bam )
  shopt -u nullglob

  if (( ${#bams[@]} != 1 )); then
    echo "ERROR: expected exactly 1 *_merged.mdup.bam in $sample_dir/paired/merged-alignment; found ${#bams[@]}" >&2
    if (( ${#bams[@]} > 0 )); then
      printf '  %s\n' "${bams[@]}" >&2
    fi
    return 6
  fi

  printf '%s\n' "${bams[0]}"
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
TUMOUR_BAM="$(resolve_tumour_bam "$SAMPLE" "$SEQ_DIR")"

MOUSE_NUM="$(extract_mouse_num "$SAMPLE")"
GERM_DIR="$(choose_germline_bam_dir "$MOUSE_NUM" "$SEQ_DIR")"
GERM_BAM="$(resolve_unique_merged_bam_in_sample_dir "$GERM_DIR")"
GERM_SM="$(bam_sm "$GERM_BAM")"

UNFILTERED="${OUT_DIR}/${SAMPLE}.unfiltered.vcf.gz"
FILTERED="${OUT_DIR}/${SAMPLE}.filtered.vcf.gz"

log "Sample: $SAMPLE"
log "Tumour: $TUMOUR_BAM (SM=$SAMPLE)"
log "Mouse: $MOUSE_NUM"
log "Normal dir: $GERM_DIR"
log "Normal: $GERM_BAM (SM=$GERM_SM)"
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
