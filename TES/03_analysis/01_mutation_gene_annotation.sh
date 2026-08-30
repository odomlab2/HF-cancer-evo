#!/usr/bin/env bash
set -eo pipefail

# --- CONFIG ---
BASE_DIR="${HF_SCC_OUTPUT_ROOT:?Set HF_SCC_OUTPUT_ROOT to an absolute output directory}/TES"
VEP_BIN="${HF_SCC_VEP_BIN:?Set HF_SCC_VEP_BIN to the VEP executable}"
VEP_CACHE="${HF_SCC_VEP_CACHE:?Set HF_SCC_VEP_CACHE}"
SAMPLE=$1
REF="${HF_SCC_VEP_REFERENCE_FASTA:?Set HF_SCC_VEP_REFERENCE_FASTA}"
INPUT_DIR="${BASE_DIR}/02_data/02_processed/mutation_calling"
OUT_DIR="${BASE_DIR}/02_data/02_processed/mutation_annotation"
mkdir -p "$OUT_DIR"
THREADS="${THREADS:-12}"

echo "Processing $SAMPLE"
      
# full annotation options can be found in: https://grch37.ensembl.org/info/docs/tools/vep/script/vep_options.html
"$VEP_BIN" \
      --fork $THREADS \
      --dir_cache "$VEP_CACHE" \
      --fasta $REF \
      --vcf \
      --compress_output gzip \
      --species mus_musculus \
      --cache \
      --cache_version 102 \
      -i ${INPUT_DIR}/${SAMPLE}.filtered.vcf.gz \
      -o ${OUT_DIR}/${SAMPLE}_annot.vcf.gz \
      --symbol \
      --variant_class \
      --sift b \
      --nearest gene \
      --distance 5000 \
      --regulatory \
      --total_length \
      --numbers \
      --protein \
      --biotype \
      --domains \
      --offline \
      -a GRCm38 \
      --force_overwrite
