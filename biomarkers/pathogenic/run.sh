#!/usr/bin/env bash
# biomarkers/pathogenic — known pathogenic mtDNA variants in a sample: the default list
# (presets/mt_pathogenic.tsv: MITOMAP Cfrm + ClinVar P/LP, tools/make_pathogenic.sh)
# and optional user lists. Option handling here; matching in run.R.

set -euo pipefail
MODULE=pathogenic
LOG_PREFIX=$MODULE
HERE=$(dirname "$(readlink -f "$0")")
source "$HERE/../../lib/common.sh"
source "$MITRIX_HOME/lib/conf.sh"
MODULE_VERSION=$(conf_get "$HERE/module.conf" version)

usage() {
  cat >&2 <<EOF
$MODULE v$MODULE_VERSION — known pathogenic mtDNA variants

USAGE
  run.sh --sample S --bam F --outdir D --mt-variants FILE [standard options] [module options]
  mitrix run ... --opt pathogenic='--list NAME=FILE [--list NAME2=FILE2] ...'

A listed variant is detected when the caller reports it (same pos/REF/ALT after
normalisation against rCRS) with HF >= HF-MIN, at least MIN-ALT-READS ALT reads, and a
FILTER without the excluded filters. Homoplasmic (HF >= 0.95) and heteroplasmic calls
both count.

MODULE OPTIONS
  --pathogenic FILE     default list (default presets/mt_pathogenic.tsv); 'none' to skip
  --list NAME=FILE      extra list, repeatable. TSV with a header and columns pos, ref, alt
                        (rCRS, 1-based, VCF-style), optional gene, disease; '#' lines are
                        comments. NAME: letters, digits, '_' (not 'pathogenic').
  --hf-min X            minimum HF (default 0.03)
  --min-alt-reads N     minimum ALT reads (default 10)
  --exclude-filters L   FILTER values to exclude (default: by caller, as heteroplasmy)

METRICS (per list; LIST = pathogenic or NAME)
  LIST_count            listed variants detected (count)
  LIST_variants         detected variants with HF, e.g. m.3243A>G(0.18);m.11778G>A(1)
                        or 'none' (label)
  LIST_max_hf           highest HF among them, 0 if none (fraction)
  LIST_low_depth_sites  listed positions with depth < MIN-ALT-READS / HF-MIN, where an
                        absent call is not informative (count); LIST_count gets qc_flag
                        low_depth when > 0
  All SKIPPED without --mt-variants.

ALSO WRITES
  <sample>.pathogenic.variants.tsv   one row per detected variant: list, variant, HF,
                                     depth, ALT reads, FILTER, homoplasmic/heteroplasmic,
                                     gene, disease, source
EOF
}

PATHO="$MITRIX_HOME/presets/mt_pathogenic.tsv"; HF_MIN=0.03; MIN_ALT=10; EXCLUDE=""; LISTS=()
while (( $# )); do
  case "$1" in
    --pathogenic)      PATHO="${2:?}"; shift 2;;
    --list)            LISTS+=("${2:?}"); shift 2;;
    --hf-min)          HF_MIN="${2:?}"; shift 2;;
    --min-alt-reads)   MIN_ALT="${2:?}"; shift 2;;
    --exclude-filters) EXCLUDE="${2:?}"; shift 2;;
    -h|--help)         usage; exit 0;;
    *)                 std_opt "$@" || die "unknown option: $1"; shift "$STD_SHIFT";;
  esac
done
check_std_args
awk -v a="$HF_MIN" 'BEGIN { exit !(a ~ /^[0-9.]+$/ && a + 0 > 0 && a + 0 <= 1) }' || die "--hf-min must be in (0, 1]: $HF_MIN"
[[ "$MIN_ALT" =~ ^[0-9]+$ ]] || die "--min-alt-reads must be a non-negative integer: $MIN_ALT"

# list specs: name=file
NAMES=(); FILES=()
if [[ "$PATHO" != none ]]; then
  [[ -s "$PATHO" ]] || die "pathogenic list not found: $PATHO (tools/make_pathogenic.sh)"
  NAMES+=(pathogenic); FILES+=("$PATHO")
fi
for l in "${LISTS[@]}"; do
  [[ "$l" == *=* ]] || die "--list must be NAME=FILE: $l"
  n=${l%%=*}; f=${l#*=}
  [[ "$n" =~ ^[A-Za-z][A-Za-z0-9_]*$ && "$n" != pathogenic ]] || die "--list name must be letters/digits/_ and not 'pathogenic': $n"
  [[ " ${NAMES[*]} " != *" $n "* ]] || die "--list name used twice: $n"
  [[ -s "$f" ]] || die "--list file not found: $f"
  NAMES+=("$n"); FILES+=("$f")
done
(( ${#NAMES[@]} )) || die "no lists: --pathogenic none and no --list"

OUT=$(module_dir "$MODULE"); mkdir -p "$OUT"
if [[ -z "$MT_VARIANTS" ]]; then
  log "no --mt-variants; writing SKIPPED rows"
  {
    echo "$RESULT_HEADER"
    for n in "${NAMES[@]}"; do
      write_row "${n}_count" NA count SKIPPED; write_row "${n}_variants" NA label SKIPPED
      write_row "${n}_max_hf" NA fraction SKIPPED; write_row "${n}_low_depth_sites" NA count SKIPPED
    done
  } > "$OUT/${SAMPLE}.${MODULE}.tsv"
  write_param hf_min "$HF_MIN" > "$OUT/${SAMPLE}.${MODULE}.params.tsv"
  exit 0
fi

need_tools Rscript bcftools
CORE=$(core_dir)
VARIANTS="$CORE/${SAMPLE}.variants.tsv"
[[ -s "$VARIANTS" ]] || die "variant table missing: $VARIANTS (run core/mt_variants.sh first)"
COVERAGE="$CORE/${SAMPLE}.chrM.per-base.bed.gz"
[[ -s "$COVERAGE" ]] || die "chrM coverage missing: $COVERAGE (run core/mt_coverage.sh first)"
if [[ -z "$EXCLUDE" ]]; then
  EXCLUDE=$(awk -F'\t' '$1 == "exclude_filters" { print $2 }' "$CORE/${SAMPLE}.variants.params.tsv")
  EXCLUDE=${EXCLUDE:-strand_bias}
fi
make_tmpdir

# normalise every list against rCRS (left-align, split) as core/mt_variants.sh does for calls
awk '/^>/ { next } { s = s $0 } END { s = substr(s, 1, 16569); print ">chrM"; for (i = 1; i <= length(s); i += 60) print substr(s, i, 60) }' \
  "$MITRIX_HOME/refs/chrMC_rCRS.fa" > "$TMP/rcrs.fa"
printf 'chrM\t16569\t6\t60\t61\n' > "$TMP/rcrs.fa.fai"
SPECS=()
for i in "${!NAMES[@]}"; do
  n=${NAMES[$i]}; f=${FILES[$i]}
  # header line (first non-comment line) must have pos, ref, alt
  awk -F'\t' -v OFS='\t' '
    /^#/ || !NF { next }
    !h { for (j = 1; j <= NF; j++) c[$j] = j; h = 1
         if (!("pos" in c) || !("ref" in c) || !("alt" in c)) { print "missing pos/ref/alt columns" > "/dev/stderr"; exit 3 }
         print "##fileformat=VCFv4.2"; print "##contig=<ID=chrM,length=16569>"
         print "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"; next }
    { g = ("gene" in c) ? $c["gene"] : ""; d = ("disease" in c) ? $c["disease"] : ""
      s = ("source" in c) ? $c["source"] : ""
      ++k
      print "chrM", $c["pos"], k, toupper($c["ref"]), toupper($c["alt"]), ".", ".", "."
      print k "\t" g "\t" d "\t" s > META }' META="$TMP/$n.meta" "$f" > "$TMP/$n.vcf" \
    || die "--list $n: $f must be a TSV with a header containing pos, ref, alt"
  bcftools norm --no-version -f "$TMP/rcrs.fa" -c wx -m -any "$TMP/$n.vcf" 2> "$TMP/$n.norm.log" \
    | bcftools query -f '%ID\t%POS\t%REF\t%ALT\n' > "$TMP/$n.ids" || die "--list $n: bcftools norm failed"
  nm=$(grep -c '^REF_MISMATCH' "$TMP/$n.norm.log" || true)
  (( nm == 0 )) || log "warning: list $n: dropped $nm entr(ies) whose REF does not match rCRS"
  awk -F'\t' -v OFS='\t' 'NR == FNR { m[$1] = $2 "\t" $3 "\t" $4; next }
    FNR == 1 { print "pos", "ref", "alt", "gene", "disease", "source" } { print $2, $3, $4, m[$1] }' \
    "$TMP/$n.meta" "$TMP/$n.ids" > "$TMP/$n.list.tsv"
  log "list $n: $(( $(wc -l < "$TMP/$n.list.tsv") - 1 )) variants ($f)"
  SPECS+=("$n=$TMP/$n.list.tsv")
done

{
  write_param hf_min          "$HF_MIN"
  write_param min_alt_reads   "$MIN_ALT"
  write_param exclude_filters "$EXCLUDE"
  for i in "${!NAMES[@]}"; do write_param "list_${NAMES[$i]}" "$(readlink -f "${FILES[$i]}")"; done
} > "$OUT/${SAMPLE}.${MODULE}.params.tsv"

Rscript --vanilla "$HERE/run.R" "$VARIANTS" "$COVERAGE" "$HF_MIN" "$MIN_ALT" "$EXCLUDE" "$SAMPLE" \
  "$MODULE_VERSION" "$OUT/${SAMPLE}.${MODULE}.tsv" "$OUT/${SAMPLE}.${MODULE}.variants.tsv" "${SPECS[@]}" \
  || die "run.R failed"
