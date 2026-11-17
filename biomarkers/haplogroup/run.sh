#!/usr/bin/env bash
# biomarkers/haplogroup — mtDNA haplogroup from HaploGrep3 (PhyloTree 17 Forensic Update)
# on the clean VCF written by core/mt_variants.sh, plus the macro-haplogroup used for
# grouping in the report.

set -euo pipefail
MODULE=haplogroup
LOG_PREFIX=$MODULE
HERE=$(dirname "$(readlink -f "$0")")
source "$HERE/../../lib/common.sh"
source "$MITRIX_HOME/lib/conf.sh"
MODULE_VERSION=$(conf_get "$HERE/module.conf" version)

usage() {
  cat >&2 <<EOF
$MODULE v$MODULE_VERSION — mtDNA haplogroup (HaploGrep3, tree $HAPLOGREP_TREE)

USAGE
  run.sh --sample S --bam F --outdir D --mt-variants FILE [--het-level X]

  --het-level X   heteroplasmies with HF > X count as alleles of the profile
                  (HaploGrep3 --hetLevel; default 0.9)

Needs HaploGrep3 (tools/install.sh, or MITRIX_HAPLOGREP3=/path/haplogrep3.jar) and
Java >= 11. Input: core/<sample>.mt.vcf.gz.

METRICS (unit label: value is text)
  haplogroup          full haplogroup, e.g. H16a1
  macro_haplogroup    L0-L6 keep the digit (L3); otherwise the leading capitals (H, HV, JT)
  haplogroup_quality  HaploGrep3 quality of the best hit (0-1; score)
  All SKIPPED without --mt-variants. HaploGrep3's report: haplogroup/<sample>.haplogrep3.txt.
EOF
}

HET_LEVEL=0.9
while (( $# )); do
  case "$1" in
    -h|--help)   usage; exit 0;;
    --het-level) HET_LEVEL="${2:?--het-level needs a value}"; shift 2;;
    *)           std_opt "$@" || die "unknown option: $1"; shift "$STD_SHIFT";;
  esac
done
check_std_args

OUT=$(module_dir "$MODULE"); mkdir -p "$OUT"
HG=NA; MACRO=NA; QUAL=NA; QC=PASS
if [[ -z "$MT_VARIANTS" ]]; then
  log "no --mt-variants; writing SKIPPED rows"
  QC=SKIPPED
else
  need_tools java
  [[ -s "$MITRIX_HAPLOGREP3" ]] || die "HaploGrep3 not found: $MITRIX_HAPLOGREP3 (run tools/install.sh)"
  VCF=$(readlink -f "$(core_dir)/${SAMPLE}.mt.vcf.gz")
  [[ -s "$VCF" ]] || die "clean VCF missing: $VCF (run core/mt_variants.sh first)"
  make_tmpdir
  java -Xmx2G -jar "$MITRIX_HAPLOGREP3" classify --tree "$HAPLOGREP_TREE" --hetLevel "$HET_LEVEL" \
      --extend-report --in "$VCF" --out "$TMP/hg.txt" > "$TMP/hg.log" 2>&1 \
    || { cat "$TMP/hg.log" >&2; die "HaploGrep3 failed"; }
  tr -d '"' < "$TMP/hg.txt" > "$OUT/${SAMPLE}.haplogrep3.txt"
  read -r HG QUAL < <(awk -F'\t' '
    NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
    NR == 2 { print $c["Haplogroup"], $c["Quality"] }' "$OUT/${SAMPLE}.haplogrep3.txt")
  if [[ -z "${HG:-}" ]]; then
    log "HaploGrep3 returned no haplogroup"
    HG=NA; QUAL=NA; QC=FAIL
  else
    MACRO=$(awk -v h="$HG" 'BEGIN {
      if (match(h, /^L[0-6]/)) print substr(h, 1, 2)
      else if (match(h, /^[A-Z]+/)) print substr(h, 1, RLENGTH)
      else print "NA" }')
    log "haplogroup $HG (macro $MACRO, quality $QUAL)"
  fi
fi

{
  echo "$RESULT_HEADER"
  write_row haplogroup         "$HG"    label "$QC"
  write_row macro_haplogroup   "$MACRO" label "$QC"
  write_row haplogroup_quality "$QUAL"  score "$QC"
} > "$OUT/${SAMPLE}.${MODULE}.tsv"
{
  write_param jar        "$MITRIX_HAPLOGREP3"
  write_param tree       "$HAPLOGREP_TREE"
  write_param het_level  "$HET_LEVEL"
} > "$OUT/${SAMPLE}.${MODULE}.params.tsv"
