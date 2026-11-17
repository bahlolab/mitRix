#!/usr/bin/env bash
# aggregate/check.sh — did every sample in the samplesheet finish? (mitrix check)
# Prints: sample  status  detail     status = OK, FAILED, INCOMPLETE or MISSING
# Exit 0 only if every sample is OK.

set -euo pipefail
LOG_PREFIX=check
source "$(dirname "$(readlink -f "$0")")/../lib/common.sh"

SS=""
while (( $# )); do
  case "$1" in
    --samplesheet) SS="${2:?}"; shift 2;;
    --outdir)      OUTDIR="${2:?}"; shift 2;;
    -h|--help)     echo "usage: mitrix check --samplesheet SS --outdir DIR" >&2; exit 0;;
    *)             die "unknown option: $1";;
  esac
done
[[ -s "$SS" ]] || die "--samplesheet not found: $SS"
[[ -d "$OUTDIR" ]] || die "--outdir not found: $OUTDIR"

samplesheet_samples() {  # samplesheet -> sample names
  awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "sample") c = i; next }
              !/^#/ && NF && c { print $c }' "$1"
}

n=0; ok=0
printf 'sample\tstatus\tdetail\n'
while read -r s; do
  n=$((n + 1))
  st="$OUTDIR/$s/$s.status.tsv"
  if [[ ! -s "$st" ]]; then
    printf '%s\tMISSING\t%s\n' "$s" "no $s.status.tsv"; continue
  fi
  failed=$(awk -F'\t' 'NR > 1 && $2 == "FAILED" { printf "%s%s (exit %s)", sep, $1, $3; sep = ", " }' "$st")
  if [[ -n "$failed" ]]; then
    printf '%s\tFAILED\t%s\n' "$s" "$failed"
  elif [[ ! -s "$OUTDIR/$s/$s.tsv" ]]; then
    printf '%s\tINCOMPLETE\t%s\n' "$s" "no $s.tsv (still running or killed?)"
  else
    printf '%s\tOK\t\n' "$s"; ok=$((ok + 1))
  fi
done < <(samplesheet_samples "$SS")

log "$ok/$n samples OK"
(( ok == n ))
