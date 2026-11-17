#!/usr/bin/env bash
# biomarkers/damage — major-arc damage as a coverage ratio (NGS analog of the
# TaqMan Ct ratio). Default ratios ND4/ND1 and ND5/ND1 from chrM per-base coverage.

set -euo pipefail
MODULE=damage
LOG_PREFIX=$MODULE
HERE=$(dirname "$(readlink -f "$0")")
source "$HERE/../../lib/common.sh"
source "$MITRIX_HOME/lib/conf.sh"
source "$MITRIX_HOME/lib/regions.sh"
source "$MITRIX_HOME/lib/stats.sh"
MODULE_VERSION=$(conf_get "$HERE/module.conf" version)

usage() {
  cat >&2 <<EOF
$MODULE v$MODULE_VERSION — major-arc damage ratios

USAGE
  run.sh --sample S --bam F --outdir D [standard options] [module options]

MODULE OPTIONS
  --ratios LIST         comma-separated: ND4_ND1, ND5_ND1, MAJOR_MINOR, custom
                        (default: ND4_ND1,ND5_ND1,MAJOR_MINOR)
                        MAJOR_MINOR = mean depth over the common-deletion span (chrM 8470-13447)
                        / ND1-ND2 on the minor arc (3307-5511): wide-region deletion burden
  --target START-END    override the target region on chrM (single ratio only)
  --reference START-END override the reference region on chrM (single ratio only)
  --target-bed FILE     override the target region with a BED (single ratio only)
  --reference-bed FILE  override the reference region with a BED (single ratio only)

METRICS (<stat> = mean or median)
  <T>_depth_<stat>, <R>_depth_<stat>, <T>_<R>_ratio_<stat>   e.g. ND4_ND1_ratio_mean
  custom: target_depth_<stat>, reference_depth_<stat>, target_reference_ratio_<stat>
EOF
}

RATIOS="ND4_ND1,ND5_ND1,MAJOR_MINOR"
T_RANGE=""; R_RANGE=""; T_BED=""; R_BED=""
while (( $# )); do
  case "$1" in
    --ratios)        RATIOS="${2:?}"; shift 2;;
    --target)        T_RANGE="${2:?}"; shift 2;;
    --reference)     R_RANGE="${2:?}"; shift 2;;
    --target-bed)    T_BED="${2:?}"; shift 2;;
    --reference-bed) R_BED="${2:?}"; shift 2;;
    -h|--help)       usage; exit 0;;
    *)               std_opt "$@" || die "unknown option: $1"; shift "$STD_SHIFT";;
  esac
done
check_std_args
need_tools zcat awk sort
load_ref_context
make_tmpdir

IFS=, read -r -a RATIO_LIST <<< "$RATIOS"
if [[ -n "$T_RANGE$R_RANGE$T_BED$R_BED" ]] && (( ${#RATIO_LIST[@]} > 1 )); then
  die "--target/--reference/--target-bed/--reference-bed need a single ratio (--ratios)"
fi

OUT=$(module_dir "$MODULE"); mkdir -p "$OUT"
PERBASE="$(core_dir)/${SAMPLE}.chrM.per-base.bed.gz"
RESULTS="$TMP/results.tsv"; PARAMS="$TMP/params.tsv"
echo "$RESULT_HEADER" > "$RESULTS"
{ write_param stat "$STAT"; write_param flag "$FLAG"; } > "$PARAMS"

declare -A DONE   # region name -> depth, so a shared reference (ND1) is reported once
depth_row() {  # name bed
  local d
  if [[ -z "${DONE[$1]:-}" ]]; then
    d=$(region_depth "$STAT" "$2" "$PERBASE")
    DONE[$1]=$d
    write_row "${1}_depth_${STAT}" "$d" x "$([[ $d == NA ]] && echo FAIL || echo PASS)" >> "$RESULTS"
    write_param "region.$1" "$(fmt_regions "$2")" >> "$PARAMS"
  fi
}

for r in "${RATIO_LIST[@]}"; do
  case "$r" in
    ND4_ND1|ND5_ND1|MAJOR_MINOR) TN=${r%_*}; RN=${r#*_};;
    custom)
      [[ -n "$T_RANGE$T_BED" && -n "$R_RANGE$R_BED" ]] \
        || die "--ratios custom needs --target/--target-bed and --reference/--reference-bed"
      TN=target; RN=reference;;
    *) die "unknown ratio '$r' (allowed: ND4_ND1, ND5_ND1, MAJOR_MINOR, custom)";;
  esac
  resolve_region "$TMP/$TN.bed" "$TN" "$MT_CHR" "$T_RANGE" "$T_BED"
  resolve_region "$TMP/$RN.bed" "$RN" "$MT_CHR" "$R_RANGE" "$R_BED"
  depth_row "$TN" "$TMP/$TN.bed"
  depth_row "$RN" "$TMP/$RN.bed"
  ratio "${DONE[$TN]}" "${DONE[$RN]}"
  write_row "${TN}_${RN}_ratio_${STAT}" "$RATIO_VAL" ratio "$RATIO_QC" >> "$RESULTS"
done

mv "$RESULTS" "$OUT/${SAMPLE}.${MODULE}.tsv"
mv "$PARAMS"  "$OUT/${SAMPLE}.${MODULE}.params.tsv"
