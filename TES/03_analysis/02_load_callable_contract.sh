#!/usr/bin/env bash

load_tes_callable_contract() {
  local root=$1 base=$2 expected declared key expected_hash actual_hash
  expected="$root/TES/01_source/tes_callable_contract.tsv"
  declared="$base/02_data/02_processed/callable_bases/run_contract.tsv"
  [[ -r $declared ]] || { echo "ERROR: missing TES callable run contract" >&2; return 6; }
  cmp -s "$expected" "$declared" || {
    echo "ERROR: TES callable run contract differs from authoritative contract" >&2
    return 6
  }
  contract_value() {
    awk -F '\t' -v key="$1" 'NR==1 {for(i=1;i<=NF;i++) if($i==key) col=i}
      NR==2 && col {print $col}' "$declared"
  }
  ASSAY=$(contract_value assay)
  DEPTH=$(contract_value min_depth)
  MAPQ=$(contract_value mapq)
  BQ=$(contract_value base_quality)
  BED="$root/$(contract_value target_bed)"
  expected_hash=$(contract_value target_bed_sha256)
  [[ $ASSAY == TES && $DEPTH =~ ^[0-9]+$ && $MAPQ =~ ^[0-9]+$ &&
     $BQ =~ ^[0-9]+$ && -r $BED ]] || {
    echo "ERROR: malformed TES callable run contract" >&2
    return 6
  }
  actual_hash=$(sha256sum "$BED" | cut -d' ' -f1)
  [[ $actual_hash == "$expected_hash" ]] || {
    echo "ERROR: TES callable target BED differs from declared contract" >&2
    return 6
  }
}
