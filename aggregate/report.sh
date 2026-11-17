#!/usr/bin/env bash
# aggregate/report.sh — HTML cohort report from cohort.tsv + metadata (mitrix report).
# Rendering is in aggregate/report.Rmd (R: rmarkdown, knitr, ggplot2; optional glmnet,
# uwot); needs pandoc. Also writes <report>_adjusted.tsv (raw, batch-adjusted, z-scores).

set -euo pipefail
LOG_PREFIX=report
HERE=$(dirname "$(readlink -f "$0")")
source "$HERE/../lib/common.sh"

usage() {
  cat >&2 <<EOF
usage: mitrix report --cohort cohort.tsv --metadata META [-o report.html] [options]

META is a TSV with a 'sample' column plus covariates (e.g. age, sex, batch, disease).

  --subtype-col NAME   metadata column with disease subtypes (default: disease)
  --covariates LIST    covariates to adjust for (default: age,sex,batch; those present)
  --keep-flagged       keep samples whose qc_flag is not PASS (default: exclude them)
  --outdir DIR         mitrix run output directory, for the per-variant heteroplasmy
                       tables (default: the directory of --cohort)
EOF
}

COHORT=""; META=""; OUT=""; RUN_OUTDIR=""; SUBTYPE=disease; COVARS=age,sex,batch; KEEP=FALSE
while (( $# )); do
  case "$1" in
    --cohort)       COHORT="${2:?}"; shift 2;;
    --metadata)     META="${2:?}"; shift 2;;
    -o)             OUT="${2:?}"; shift 2;;
    --subtype-col)  SUBTYPE="${2:?}"; shift 2;;
    --covariates)   COVARS="${2:?}"; shift 2;;
    --keep-flagged) KEEP=TRUE; shift;;
    --outdir)       RUN_OUTDIR="${2:?}"; shift 2;;
    -h|--help)      usage; exit 0;;
    *)              die "unknown option: $1";;
  esac
done
[[ -s "$COHORT" ]] || die "--cohort not found: $COHORT (run mitrix cohort first)"
[[ -s "$META" ]]   || die "--metadata not found: $META"
OUT=$(readlink -f "${OUT:-$(dirname "$COHORT")/report.html}")
RUN_OUTDIR=${RUN_OUTDIR:-$(dirname "$COHORT")}
need_tools Rscript
command -v pandoc >/dev/null 2>&1 || die "pandoc not found in PATH (needed by rmarkdown; e.g. module load pandoc)"
make_tmpdir

Rscript --vanilla -e '
  a <- commandArgs(TRUE)
  rmarkdown::render(a[1], output_file = basename(a[2]), output_dir = dirname(a[2]),
    intermediates_dir = a[3], knit_root_dir = a[3], quiet = TRUE,
    params = list(cohort = a[4], metadata = a[5], subtype_col = a[6], covariates = a[7],
                  keep_flagged = as.logical(a[8]), adjusted_out = a[9], mitrix_version = a[10],
                  outdir = a[11], annotation = a[12]))
' "$HERE/report.Rmd" "$OUT" "$TMP" "$(readlink -f "$COHORT")" "$(readlink -f "$META")" \
  "$SUBTYPE" "$COVARS" "$KEEP" "${OUT%.html}_adjusted.tsv" "$MITRIX_VERSION" \
  "$(readlink -f "$RUN_OUTDIR")" "$MITRIX_HOME/presets/mt_annotation.tsv.gz" \
  || die "rendering failed"
log "wrote $OUT"
