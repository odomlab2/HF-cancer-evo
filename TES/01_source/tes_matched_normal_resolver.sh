#!/usr/bin/env bash
tes_validate_matched_normals() {
  local manifest=$1
  awk -F '\t' '
    function fail(message) { print "ERROR: " message > "/dev/stderr"; bad=1 }
    function register_id(id, physical, kind) {
      if (id == "") return
      if (id in id_physical) {
        fail(sprintf("TES normal identifier \"%s\" is assigned more than once (%s); physical controls: \"%s\" and \"%s\"", id, kind, id_physical[id], physical))
      } else id_physical[id]=physical
    }
    NR==1 {
      expected="mouse_group\toperational_normal_id\tphysical_control_id\tsample_role\taliases"
      if ($0 != expected) fail("unexpected TES matched-normal manifest header")
      next
    }
    {
      rows++
      if ($1=="" || $2=="" || $3=="" || $4!="matched_normal")
        fail("invalid TES matched-normal row " NR)
      if ($1 in mouse_physical)
        fail(sprintf("multiple physical controls resolve for mouse_group \"%s\": \"%s\" and \"%s\"", $1, mouse_physical[$1], $3))
      else mouse_physical[$1]=$3
      if ($3 in physical_row)
        fail(sprintf("physical control \"%s\" is represented by canonical rows %s and %s", $3, physical_row[$3], NR))
      else physical_row[$3]=NR
      register_id($2, $3, "operational_normal_id")
      count=split($5, alias, ";")
      for (i=1; i<=count; i++) register_id(alias[i], $3, "alias")
    }
    END {
      if (rows == 0) fail("TES matched-normal manifest has no controls")
      if (bad) exit 2
    }' "$manifest"
}
tes_mouse_group() {
  local sample=${1##*/}
  if [[ "$sample" =~ ^(hairfollicle|skin)([0-9]+) ]]; then
    printf '%s\n' "${BASH_REMATCH[2]}"
  else
    printf "ERROR: invalid TES tumour identifier '%s'\n" "$sample" >&2
    return 2
  fi
}
tes_resolve_matched_normal() {
  local manifest=$1 tumour=$2 mouse
  tes_validate_matched_normals "$manifest" || return
  mouse=$(tes_mouse_group "$tumour") || return
  awk -F '\t' -v mouse="$mouse" -v tumour="$tumour" '
    NR>1 && $1==mouse {normal=$2; found++}
    END {
      if (found != 1) {
        printf "ERROR: TES tumour \"%s\" resolves to %d controls for mouse_group \"%s\"\n",
          tumour, found, mouse > "/dev/stderr"
        exit 3
      }
      print normal
    }' "$manifest"
}
tes_resolve_control_identity() {
  local manifest=$1 identifier=$2
  tes_validate_matched_normals "$manifest" || return
  awk -F '\t' -v id="$identifier" '
    NR>1 {
      hit=($2==id)
      count=split($5, alias, ";")
      for (i=1; i<=count; i++) if (alias[i]==id) hit=1
      if (hit) {printf "%s\t%s\n", $3, $2; found++}
    }
    END {
      if (found != 1) {
        printf "ERROR: TES normal identifier \"%s\" resolves to %d physical controls\n",
          id, found > "/dev/stderr"
        exit 4
      }
    }' "$manifest"
}
tes_validate_tumour_pairings() {
  local manifest=$1 samples=$2 sample
  tes_validate_matched_normals "$manifest" || return
  declare -A seen=()
  while IFS= read -r sample || [[ -n "$sample" ]]; do
    [[ -n "$sample" ]] || continue
    if [[ -n "${seen[$sample]:-}" ]]; then
      printf "ERROR: duplicate TES tumour identifier '%s'\n" "$sample" >&2
      return 5
    fi
    seen[$sample]=1
    if tes_resolve_control_identity "$manifest" "$sample" >/dev/null 2>&1; then
      printf "ERROR: matched-normal identifier '%s' is in the TES tumour universe\n" \
        "$sample" >&2
      return 6
    fi
    tes_resolve_matched_normal "$manifest" "$sample" >/dev/null || return
  done < "$samples"
}
