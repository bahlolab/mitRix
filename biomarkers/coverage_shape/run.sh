#!/usr/bin/env bash
# biomarkers/coverage_shape — per-sample model of the chrM coverage shape, to separate the
# technical GC bias from replication intermediates and the common deletion.
# Bins of the corrected (circularised) per-base depth, robust fit (MASS::rlm):
#   log(bin depth / median) ~ GC + GC^2 + e + deletion_span
#   e = L-strand distance from O_L, (pos - O_L) mod 16569 / 16569 (0 just after O_L, 1 just before).
#   Strand-displacement replication: the L strand starts at O_L once the H strand has passed it
#   and runs pos 5.8k -> 16.5k -> 1 -> 5.7k, so the second double-stranded copy appears first just
#   after O_L and last in the minor arc: replication intermediates = coverage falling with e.
# The control region (7S DNA / D-loop excess, d_loop_7s module) is excluded.
set -euo pipefail
MODULE=coverage_shape
LOG_PREFIX=$MODULE
HERE=$(dirname "$(readlink -f "$0")")
source "$HERE/../../lib/common.sh"
source "$MITRIX_HOME/lib/conf.sh"
MODULE_VERSION=$(conf_get "$HERE/module.conf" version)

usage() {
  cat >&2 <<EOF2
$MODULE v$MODULE_VERSION — chrM coverage shape (GC bias, replication gradient, deletion step)

USAGE
  run.sh --sample S --bam F --outdir D [--bin 100] [--exclude 16024-576] [--ol 5761]
         [--deletion-span 8470-13447] [standard options]

METRICS
  gc_slope, gc_curvature     GC bias: log depth per 10% bin GC (linear, quadratic), centred at 44%
  replication_gradient       fitted depth just after O_L / just before (full circle); > 1 = excess of
                             replication intermediates (strand-displacement model)
  deletion_step              depth inside the common-deletion span vs outside, after GC and gradient;
                             < 1 = loss (major-arc deletions)
  shape_r2_gc, shape_r2_full variance of log bin depth explained by GC / by the full model (OLS)
  shape_residual_sd          robust SD (MAD) of residuals — coverage noise
  local_dip, _pos, _z        most negative bin residual (fold), its position and z (local dips / artefacts)
  ND4_ND1_ratio_raw/_gc_adj  ND4/ND1 TaqMan amplicons, raw and divided by the fitted GC effect
  MAJOR_MINOR_ratio_raw/_gc_adj  common-deletion span / ND1-ND2 (3307-5511), raw and GC-adjusted
  shape_bins                 bins used
  rCRS coordinates (1-based); assumes a 16569-bp rCRS mt contig. qc: FAIL if < 100 bins or median
  bin depth < 50; LOW_DEPTH if < 200.
EOF2
}

BIN=100; EXCL=16024-576; OL=5761; DEL=8470-13447
while (( $# )); do
  case "$1" in
    --bin)           BIN="${2:?}"; shift 2;;
    --exclude)       EXCL="${2:?}"; shift 2;;
    --ol)            OL="${2:?}"; shift 2;;
    --deletion-span) DEL="${2:?}"; shift 2;;
    -h|--help) usage; exit 0;;
    *)         std_opt "$@" || die "unknown option: $1"; shift "$STD_SHIFT";;
  esac
done
check_std_args
need_tools Rscript
load_ref_context
make_tmpdir
[[ "${MT_LEN:-16569}" == 16569 ]] || die "mt contig length $MT_LEN; coverage_shape assumes rCRS (16569 bp)"

OUT=$(module_dir "$MODULE"); mkdir -p "$OUT"
CORRECTED="$(core_dir)/${SAMPLE}.chrM.corrected.per-base.bed.gz"
[[ -s "$CORRECTED" ]] || die "corrected per-base missing: $CORRECTED (run core/mt_circular.sh first)"

Rscript --vanilla "$HERE/run.R" "$CORRECTED" "$MITRIX_HOME/refs/chrMC_rCRS.fa" "$TMP/rows.tsv" \
  "$BIN" "$EXCL" "$OL" "$DEL" 16569 || die "coverage-shape fit failed"
{
  echo "$RESULT_HEADER"
  while IFS=$'\t' read -r m v u q; do write_row "$m" "$v" "$u" "$q"; done < "$TMP/rows.tsv"
} > "$TMP/results.tsv"
log "$(awk -F'\t' '$3 ~ /^(gc_slope|replication_gradient|deletion_step|ND4_ND1_ratio_gc_adj)$/ { printf "%s %s; ", $3, $4 }' "$TMP/results.tsv")"
{
  write_param bin "$BIN"; write_param exclude "$EXCL"; write_param O_L "$OL"; write_param deletion_span "$DEL"
  write_param per_base "$(basename "$CORRECTED")"
} > "$TMP/params.tsv"
mv "$TMP/results.tsv" "$OUT/${SAMPLE}.${MODULE}.tsv"
mv "$TMP/params.tsv"  "$OUT/${SAMPLE}.${MODULE}.params.tsv"
