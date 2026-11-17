#!/usr/bin/env bash
# biomarkers/d_loop_7s — 7S DNA / replication proxy: D-loop depth over ND1 depth, and the
# non-coding-region coverage ratios of Gupta et al. (Nature 2023): 7S DNA flap / 7S DNA and
# 7S DNA flap / 7S RNA primer (median depths). All from the corrected chrM per-base
# (core/mt_circular.sh), so reads crossing the 16569→1 seam are counted.

set -euo pipefail
MODULE=d_loop_7s
LOG_PREFIX=$MODULE
HERE=$(dirname "$(readlink -f "$0")")
source "$HERE/../../lib/common.sh"
source "$MITRIX_HOME/lib/conf.sh"
source "$MITRIX_HOME/lib/regions.sh"
source "$MITRIX_HOME/lib/stats.sh"
MODULE_VERSION=$(conf_get "$HERE/module.conf" version)

usage() {
  cat >&2 <<EOF
$MODULE v$MODULE_VERSION — D-loop / ND1 on circularised chrM

USAGE
  run.sh --sample S --bam F --outdir D [standard options]

METRICS (<stat> = mean or median)
  DLOOP_depth_<stat>, ND1_depth_<stat>, DLOOP_ND1_ratio_<stat>
  7S_depth_median       7S DNA, first third (O_H -> 16542; $R7S)
  flap_depth_median     7S DNA flap (before CSBII -> O_H; $RFLAP; excludes the D310 poly-C tract)
  primer_depth_median   7S RNA primer (LSP -> CSBIII; $RPRIMER)
  flap_7S_ratio         flap_depth_median / 7S_depth_median
  flap_primer_ratio     flap_depth_median / primer_depth_median
  (direction as Gupta et al.: flap coverage relative to 7S DNA or RNA primer)
  Regions: Gupta et al. 2023 landmarks (rCRS, 1-based inclusive); always medians, as in
  that paper. Override with --region-7s, --region-flap, --region-primer START-END.
EOF
}

R7S=16542-191; RFLAP=191-298; RPRIMER=346-407
while (( $# )); do
  case "$1" in
    --region-7s)     R7S="${2:?}"; shift 2;;
    --region-flap)   RFLAP="${2:?}"; shift 2;;
    --region-primer) RPRIMER="${2:?}"; shift 2;;
    -h|--help) usage; exit 0;;
    *)         std_opt "$@" || die "unknown option: $1"; shift "$STD_SHIFT";;
  esac
done
check_std_args
need_tools zcat awk sort
load_ref_context
make_tmpdir

OUT=$(module_dir "$MODULE"); mkdir -p "$OUT"
CORRECTED="$(core_dir)/${SAMPLE}.chrM.corrected.per-base.bed.gz"
[[ -s "$CORRECTED" ]] || die "corrected per-base missing: $CORRECTED (run core/mt_circular.sh first)"

resolve_region "$TMP/DLOOP.bed" DLOOP "$MT_CHR"
resolve_region "$TMP/ND1.bed"   ND1   "$MT_CHR"
T=$(region_depth "$STAT" "$TMP/DLOOP.bed" "$CORRECTED")
R=$(region_depth "$STAT" "$TMP/ND1.bed"   "$CORRECTED")
log "D-loop coverage (corrected per-base): $T"
ratio "$T" "$R"
DL_VAL=$RATIO_VAL; DL_QC=$RATIO_QC

# Gupta et al. 2023 non-coding-region medians and ratios
resolve_region "$TMP/7S.bed"     7S     "$MT_CHR" "$R7S"
resolve_region "$TMP/flap.bed"   flap   "$MT_CHR" "$RFLAP"
resolve_region "$TMP/primer.bed" primer "$MT_CHR" "$RPRIMER"
D7S=$(region_depth median "$TMP/7S.bed" "$CORRECTED")
DFL=$(region_depth median "$TMP/flap.bed" "$CORRECTED")
DPR=$(region_depth median "$TMP/primer.bed" "$CORRECTED")
ratio "$DFL" "$D7S"; S7_VAL=$RATIO_VAL; S7_QC=$RATIO_QC
ratio "$DFL" "$DPR"; PR_VAL=$RATIO_VAL; PR_QC=$RATIO_QC
log "medians: 7S $D7S, flap $DFL, primer $DPR; flap/7S $S7_VAL, flap/primer $PR_VAL"

{
  echo "$RESULT_HEADER"
  write_row "DLOOP_depth_${STAT}"     "$T" x "$([[ $T == NA ]] && echo FAIL || echo PASS)"
  write_row "ND1_depth_${STAT}"       "$R" x "$([[ $R == NA ]] && echo FAIL || echo PASS)"
  write_row "DLOOP_ND1_ratio_${STAT}" "$DL_VAL" ratio "$DL_QC"
  write_row 7S_depth_median     "$D7S" x "$([[ $D7S == NA ]] && echo FAIL || echo PASS)"
  write_row flap_depth_median   "$DFL" x "$([[ $DFL == NA ]] && echo FAIL || echo PASS)"
  write_row primer_depth_median "$DPR" x "$([[ $DPR == NA ]] && echo FAIL || echo PASS)"
  write_row flap_7S_ratio       "$S7_VAL" ratio "$S7_QC"
  write_row flap_primer_ratio   "$PR_VAL" ratio "$PR_QC"
} > "$TMP/results.tsv"
{
  write_param stat "$STAT"
  write_param flag "$FLAG"
  write_param region.DLOOP "$(fmt_regions "$TMP/DLOOP.bed")"
  write_param region.ND1   "$(fmt_regions "$TMP/ND1.bed")"
  write_param region.7S     "$(fmt_regions "$TMP/7S.bed")"
  write_param region.flap   "$(fmt_regions "$TMP/flap.bed")"
  write_param region.primer "$(fmt_regions "$TMP/primer.bed")"
} > "$TMP/params.tsv"

mv "$TMP/results.tsv" "$OUT/${SAMPLE}.${MODULE}.tsv"
mv "$TMP/params.tsv"  "$OUT/${SAMPLE}.${MODULE}.params.tsv"
