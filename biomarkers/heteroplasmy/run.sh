#!/usr/bin/env bash
# biomarkers/heteroplasmy — heteroplasmy burden, mutational spectrum and mitochondrial
# local constraint score sum (MSS) from the standard variant table (core/mt_variants.sh).
# Option handling here; the calculation is in run.R.

set -euo pipefail
MODULE=heteroplasmy
LOG_PREFIX=$MODULE
HERE=$(dirname "$(readlink -f "$0")")
source "$HERE/../../lib/common.sh"
source "$MITRIX_HOME/lib/conf.sh"
MODULE_VERSION=$(conf_get "$HERE/module.conf" version)

usage() {
  cat >&2 <<EOF
$MODULE v$MODULE_VERSION — heteroplasmy burden, spectrum and MSS

USAGE
  run.sh --sample S --bam F --outdir D --mt-variants FILE [standard options] [module options]

A heteroplasmy is a variant with HF-MIN <= HF <= HF-MAX (HF = heteroplasmy fraction,
the caller's VAF), supported by at least MIN-ALT-READS reads (HF x depth), not at a
blocklisted position, and whose FILTER has none of the excluded filters.

MODULE OPTIONS
  --hf-min X            lower HF bound (default 0.03)
  --hf-max X            upper HF bound (default 0.95)
  --hf-split X          low/high HF bins: low <= X < high (default 0.2)
  --min-alt-reads N     minimum ALT reads (default 10)
  --exclude-filters L   comma-separated FILTER values to exclude (default: by caller —
                        strand_bias for MitoHPC/Mutect2/mtDNA-Server 2, FAIL for mity)
  --blocklist FILE      positions to exclude: TSV start, end, note (1-based, inclusive)
                        (default presets/heteroplasmy_blocklist.tsv)
  --numt FILE           alleles found in NUMT sequences (default presets/numt_alleles.tsv,
                        from MitoHPC); 'none' to keep them
  --numt-k X            NUMT-copy equivalents (default 1.5): a NUMT allele with
                        HF < X x nuclear depth / chrM depth is excluded as NUMT-like
  --numt-max-hf X       ceiling used when the nuclear depth is unknown (default 0.10)
  Nuclear depth: --nuclear-depth (or samplesheet column nuclear_depth), else the
  copy_number module's B2M depth when it ran, else unknown.

METRICS (counts unless noted)
  heteroplasmy_count          all heteroplasmies
  heteroplasmy_low_count      HF <= hf-split          (likely somatic)
  heteroplasmy_high_count     HF >  hf-split          (often inherited)
  snv_count, indel_count      by variant type
  H_CtoT_count, H_AtoG_count  heavy-strand C>T and A>G (rCRS G>A, T>C)
  L_AtoG_count, L_CtoT_count  light-strand A>G and C>T (rCRS A>G, C>T)
  transversion_count          all 8 transversion types
  age_accumulating_count      H C>T + H A>G + L A>G outside the replication origin
                              (chrM:16172-210) (Gupta et al., Nature 2026)
  age_accumulating_low_count  the same with HF <= hf-split
  MSS                         sum of MLC scores (Lake et al. 2024) (score)
  MSS_low, MSS_high           MSS by HF bin (score)
  MSS_missing_count           heteroplasmies without an MLC score (indels, MNVs)
  nonsynonymous_count         missense, stop_gained, start_lost, stop_lost
  MSS_nonsynonymous           MSS over those (score)
  mt_median_depth             median chrM depth (x)
  hf_detection_limit          max(hf-min, min-alt-reads / mt_median_depth) (fraction)
  numt_excluded_count         NUMT-like calls left out (count)
  numt_hf_threshold           HF below which a NUMT allele is left out (fraction)
  Low-bin metrics get qc_flag low_depth when hf_detection_limit > hf-min.
  All SKIPPED without --mt-variants.

ALSO WRITES
  <sample>.heteroplasmy.variants.tsv   one row per heteroplasmy: HF, depth, ALT reads,
                                       strand, substitution, HF bin, age-accumulating,
                                       gene, region, consequence, MLC score
EOF
}

HF_MIN=0.03; HF_MAX=0.95; HF_SPLIT=0.2; MIN_ALT=10; EXCLUDE=""
BLOCKLIST="$MITRIX_HOME/presets/heteroplasmy_blocklist.tsv"
NUMT="$MITRIX_HOME/presets/numt_alleles.tsv"; NUMT_K=1.5; NUMT_MAX=0.10
ANNOTATION="$MITRIX_HOME/presets/mt_annotation.tsv.gz"
while (( $# )); do
  case "$1" in
    --hf-min)          HF_MIN="${2:?}"; shift 2;;
    --hf-max)          HF_MAX="${2:?}"; shift 2;;
    --hf-split)        HF_SPLIT="${2:?}"; shift 2;;
    --min-alt-reads)   MIN_ALT="${2:?}"; shift 2;;
    --exclude-filters) EXCLUDE="${2:?}"; shift 2;;
    --blocklist)       BLOCKLIST="${2:?}"; shift 2;;
    --numt)            NUMT="${2:?}"; shift 2;;
    --numt-k)          NUMT_K="${2:?}"; shift 2;;
    --numt-max-hf)     NUMT_MAX="${2:?}"; shift 2;;
    -h|--help)         usage; exit 0;;
    *)                 std_opt "$@" || die "unknown option: $1"; shift "$STD_SHIFT";;
  esac
done
check_std_args
awk -v a="$HF_MIN" -v b="$HF_MAX" -v s="$HF_SPLIT" 'BEGIN {
  ok = a ~ /^[0-9.]+$/ && b ~ /^[0-9.]+$/ && s ~ /^[0-9.]+$/ && a+0 >= 0 && b+0 <= 1 && a+0 <= b+0
  exit !ok }' || die "--hf-min/--hf-max/--hf-split must be numbers with 0 <= min <= max <= 1 (got $HF_MIN, $HF_MAX, $HF_SPLIT)"
[[ "$MIN_ALT" =~ ^[0-9]+$ ]] || die "--min-alt-reads must be a non-negative integer: $MIN_ALT"
[[ -s "$BLOCKLIST" ]] || die "blocklist not found: $BLOCKLIST"
[[ -s "$ANNOTATION" ]] || die "annotation not found: $ANNOTATION (tools/make_annotation.sh)"
[[ "$NUMT" == none || -s "$NUMT" ]] || die "NUMT allele list not found: $NUMT (tools/make_numt.sh)"
awk -v k="$NUMT_K" -v m="$NUMT_MAX" 'BEGIN { exit !(k ~ /^[0-9.]+$/ && m ~ /^[0-9.]+$/) }' || die "--numt-k/--numt-max-hf must be numbers"

OUT=$(module_dir "$MODULE"); mkdir -p "$OUT"
METRICS="heteroplasmy_count heteroplasmy_low_count heteroplasmy_high_count snv_count indel_count
  H_CtoT_count H_AtoG_count L_AtoG_count L_CtoT_count transversion_count age_accumulating_count
  age_accumulating_low_count MSS MSS_low MSS_high MSS_missing_count nonsynonymous_count
  MSS_nonsynonymous mt_median_depth hf_detection_limit numt_excluded_count numt_hf_threshold"

if [[ -z "$MT_VARIANTS" ]]; then
  log "no --mt-variants; writing SKIPPED rows"
  {
    echo "$RESULT_HEADER"
    for m in $METRICS; do
      case "$m" in *_count) u=count;; MSS*) u=score;; mt_median_depth) u=x;; hf_detection_limit|numt_hf_threshold) u=fraction;; *) u=count;; esac
      write_row "$m" NA "$u" SKIPPED
    done
  } > "$OUT/${SAMPLE}.${MODULE}.tsv"
  write_param hf_min "$HF_MIN" > "$OUT/${SAMPLE}.${MODULE}.params.tsv"
  exit 0
fi

need_tools Rscript
CORE=$(core_dir)
VARIANTS="$CORE/${SAMPLE}.variants.tsv"
[[ -s "$VARIANTS" ]] || die "variant table missing: $VARIANTS (run core/mt_variants.sh first)"
COVERAGE="$CORE/${SAMPLE}.chrM.per-base.bed.gz"
[[ -s "$COVERAGE" ]] || die "chrM coverage missing: $COVERAGE (run core/mt_coverage.sh first)"
if [[ -z "$EXCLUDE" ]]; then
  EXCLUDE=$(awk -F'\t' '$1 == "exclude_filters" { print $2 }' "$CORE/${SAMPLE}.variants.params.tsv")
  EXCLUDE=${EXCLUDE:-strand_bias}
fi
CALLER=$(awk -F'\t' '$1 == "caller" { print $2 }' "$CORE/${SAMPLE}.variants.params.tsv")

# nuclear depth for the NUMT ceiling: --nuclear-depth, else copy_number's B2M depth
NUC="${NUCLEAR_DEPTH:-}"; NUC_SRC=option
if [[ -z "$NUC" ]]; then
  CN="$(module_dir copy_number)/${SAMPLE}.copy_number.tsv"
  [[ -s "$CN" ]] && NUC=$(awk -F'\t' '$3 ~ /^B2M_depth_(median|mean)$/ && $4 ~ /^[0-9.]+$/ { print $4; exit }' "$CN")
  NUC_SRC=copy_number_B2M
fi
[[ -n "$NUC" ]] || { NUC=NA; NUC_SRC=unknown; }
{
  write_param caller          "$CALLER"
  write_param numt_alleles    "$([[ $NUMT == none ]] && echo none || readlink -f "$NUMT")"
  write_param numt_k          "$NUMT_K"
  write_param numt_max_hf     "$NUMT_MAX"
  write_param nuclear_depth   "$NUC ($NUC_SRC)"
  write_param hf_min          "$HF_MIN"
  write_param hf_max          "$HF_MAX"
  write_param hf_split        "$HF_SPLIT"
  write_param min_alt_reads   "$MIN_ALT"
  write_param exclude_filters "$EXCLUDE"
  write_param blocklist       "$(readlink -f "$BLOCKLIST")"
  write_param annotation      "$(readlink -f "$ANNOTATION")"
} > "$OUT/${SAMPLE}.${MODULE}.params.tsv"

Rscript --vanilla "$HERE/run.R" "$VARIANTS" "$COVERAGE" "$BLOCKLIST" "$ANNOTATION" \
  "$HF_MIN" "$HF_MAX" "$HF_SPLIT" "$MIN_ALT" "$EXCLUDE" "$SAMPLE" "$MODULE_VERSION" \
  "$OUT/${SAMPLE}.${MODULE}.tsv" "$OUT/${SAMPLE}.${MODULE}.variants.tsv" \
  "$NUMT" "$NUMT_K" "$NUMT_MAX" "$NUC" \
  || die "run.R failed"
