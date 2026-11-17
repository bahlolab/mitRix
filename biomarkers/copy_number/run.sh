#!/usr/bin/env bash
# biomarkers/copy_number — mtDNA copy number.
# TaqMan mode (as v0.1.0): ND1 depth (chrM per-base) over B2M depth (chr15, mosdepth
# on the B2M region only; uses the index, does not scan chr15). B2M coordinates are
# picked from the detected build. A WGS mode is planned for mitRix Plus.

set -euo pipefail
MODULE=copy_number
LOG_PREFIX=$MODULE
HERE=$(dirname "$(readlink -f "$0")")
source "$HERE/../../lib/common.sh"
source "$MITRIX_HOME/lib/conf.sh"
source "$MITRIX_HOME/lib/regions.sh"
source "$MITRIX_HOME/lib/stats.sh"
MODULE_VERSION=$(conf_get "$HERE/module.conf" version)

usage() {
  cat >&2 <<EOF
$MODULE v$MODULE_VERSION — mtDNA copy number

USAGE
  run.sh --sample S --bam F --outdir D [standard options] [module options]

MODULE OPTIONS
  --mode taqman         ND1/B2M depth ratio (default; the only mode in this version)
  --target START-END    override ND1 on chrM
  --target-bed FILE     override ND1 with a BED
  --reference START-END override B2M on chr15
  --reference-bed FILE  override B2M with a BED (any chromosome)

METRICS (<stat> = mean or median)
  ND1_depth_<stat>, B2M_depth_<stat>, ND1_B2M_ratio_<stat>
  B2M and ratio are SKIPPED when the header has no chr15 and no --reference-bed.
EOF
}

MODE=taqman
T_RANGE=""; R_RANGE=""; T_BED=""; R_BED=""
while (( $# )); do
  case "$1" in
    --mode)          MODE="${2:?}"; shift 2;;
    --target)        T_RANGE="${2:?}"; shift 2;;
    --reference)     R_RANGE="${2:?}"; shift 2;;
    --target-bed)    T_BED="${2:?}"; shift 2;;
    --reference-bed) R_BED="${2:?}"; shift 2;;
    -h|--help)       usage; exit 0;;
    *)               std_opt "$@" || die "unknown option: $1"; shift "$STD_SHIFT";;
  esac
done
[[ "$MODE" == taqman ]] || die "unknown --mode '$MODE' (only taqman for now)"
check_std_args
need_tools mosdepth zcat awk sort
load_ref_context
make_tmpdir

OUT=$(module_dir "$MODULE"); mkdir -p "$OUT"
PERBASE="$(core_dir)/${SAMPLE}.chrM.per-base.bed.gz"

resolve_region "$TMP/ND1.bed" ND1 "$MT_CHR" "$T_RANGE" "$T_BED"
T=$(region_depth "$STAT" "$TMP/ND1.bed" "$PERBASE")

# B2M: user BED, user range on chr15, else preset for the detected build
B2M_QC=PASS
if [[ -n "$R_BED" ]]; then
  resolve_region "$TMP/B2M.bed" B2M "" "" "$R_BED"
elif [[ -z "$CHR15" ]]; then
  log "skipping mtDNA-CN: chr15 not found in BAM header (chrM-only BAM). Supply --reference-bed to use a custom B2M region."
  B2M_QC=SKIPPED
elif [[ -n "$R_RANGE" ]]; then
  resolve_region "$TMP/B2M.bed" B2M "$CHR15" "$R_RANGE"
elif [[ "$BUILD" == unknown ]]; then
  log "could not auto-detect nuclear build from chr15 length=$CHR15_LEN; supply --reference-bed FILE."
  B2M_QC=FAIL
else
  resolve_region "$TMP/B2M.bed" B2M "$CHR15" "" "" "$BUILD"
  log "nuclear build auto-detected: $BUILD (chr15 length=$CHR15_LEN); using B2M as reference"
fi

R=NA
if [[ "$B2M_QC" == PASS ]]; then
  # mosdepth --use-median reports the per-region median in col 4
  MD_OPT=(); [[ "$STAT" == median ]] && MD_OPT=(--use-median)
  mosdepth -n -t "$THREADS" -F "$FLAG" --by "$TMP/B2M.bed" "${MD_OPT[@]}" "${MOSDEPTH_FASTA_OPT[@]}" \
    "$TMP/b2m" "$BAM" 2> >(grep -v "\[W::hts_idx_load3\]" >&2) \
    || die "mosdepth failed for B2M"
  R=$(weighted_mean "$TMP/b2m.regions.bed.gz")
  [[ "$R" == NA ]] && B2M_QC=FAIL
  ratio "$T" "$R"
else
  RATIO_VAL=NA; RATIO_QC=$B2M_QC
fi

{
  echo "$RESULT_HEADER"
  write_row "ND1_depth_${STAT}"     "$T" x "$([[ $T == NA ]] && echo FAIL || echo PASS)"
  write_row "B2M_depth_${STAT}"     "$R" x "$B2M_QC"
  write_row "ND1_B2M_ratio_${STAT}" "$RATIO_VAL" ratio "$RATIO_QC"
} > "$TMP/results.tsv"
{
  write_param mode "$MODE"
  write_param stat "$STAT"
  write_param flag "$FLAG"
  write_param build "$BUILD"
  write_param region.ND1 "$(fmt_regions "$TMP/ND1.bed")"
  if [[ -s "$TMP/B2M.bed" ]]; then write_param region.B2M "$(fmt_regions "$TMP/B2M.bed")"; fi
} > "$TMP/params.tsv"

mv "$TMP/results.tsv" "$OUT/${SAMPLE}.${MODULE}.tsv"
mv "$TMP/params.tsv"  "$OUT/${SAMPLE}.${MODULE}.params.tsv"
