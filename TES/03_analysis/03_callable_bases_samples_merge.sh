#!/usr/bin/env bash
set -euo pipefail

# --- CONFIG ---
HF_SCC_ROOT="${HF_SCC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
HF_SCC_OUTPUT_ROOT="${HF_SCC_OUTPUT_ROOT:?Set HF_SCC_OUTPUT_ROOT to an absolute output directory}"
[[ "$HF_SCC_OUTPUT_ROOT" == /* ]] || { echo "ERROR: HF_SCC_OUTPUT_ROOT must be absolute" >&2; exit 2; }
BASE_DIR="${HF_SCC_OUTPUT_ROOT%/}/TES"
cd "$BASE_DIR"
ROWS="${BASE_DIR}/02_data/02_processed/callable_bases/rows"
SUMMARY="${BASE_DIR}/02_data/02_processed/callable_summary.tsv"
VALIDATOR="${HF_SCC_ROOT}/TES/03_analysis/03_validate_callable_rows.r"
CONTRACT="${HF_SCC_ROOT}/TES/01_source/tes_callable_contract.tsv"
RUN_CONTRACT="${BASE_DIR}/02_data/02_processed/callable_bases/run_contract.tsv"
R_BIN="${R_BIN:-R}"

HEADER=$("$R_BIN" --vanilla --slave -f "$VALIDATOR" --args \
  "$HF_SCC_ROOT" "$ROWS" "$CONTRACT" "$RUN_CONTRACT")

shopt -s nullglob
rows=( "${ROWS}/"*.tsv )
temporary=$(mktemp "${SUMMARY}.tmp.XXXXXX")
trap 'rm -f "$temporary"' EXIT
{
  printf '%s\n' "$HEADER"
  cat "${rows[@]}" | sort -V
} > "$temporary"
mv "$temporary" "$SUMMARY"
trap - EXIT

echo "[OK] Wrote: $SUMMARY"
