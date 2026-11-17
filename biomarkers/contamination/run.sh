#!/usr/bin/env bash
# biomarkers/contamination — mtDNA contamination from haplocheck on the clean VCF
# written by core/mt_variants.sh. `mitrix cohort` turns contamination_level into a
# sample-level QC flag (contaminated when level > --max-contamination, default 0.02).

set -euo pipefail
MODULE=contamination
LOG_PREFIX=$MODULE
HERE=$(dirname "$(readlink -f "$0")")
source "$HERE/../../lib/common.sh"
source "$MITRIX_HOME/lib/conf.sh"
MODULE_VERSION=$(conf_get "$HERE/module.conf" version)

usage() {
  cat >&2 <<EOF
$MODULE v$MODULE_VERSION — mtDNA contamination (haplocheck)

USAGE
  run.sh --sample S --bam F --outdir D --mt-variants FILE [--min-support N] [--lowconf-depth X]

A haplocheck YES is kept only when its minor haplogroup is supported by at least
MIN-SUPPORT minor heteroplasmies (default 5) and is not rCRS itself (H2a2a1*, the
artefact of near-homoplasmic calls read as heteroplasmic); otherwise the status is NO
and the row is flagged unsupported. Levels below 5% at haplocheck coverage below
LOWCONF-DEPTH (default 2000; Weissensteiner et al. 2021: >600x for 1%, and low-mtCN
samples give unreliable levels) are flagged low_confidence.

Needs haplocheck (tools/install.sh, or MITRIX_HAPLOCHECK=/path/haplocheck.jar) and
Java >= 11. Input: core/<sample>.mt.vcf.gz.

METRICS
  contamination_status   YES, NO or ND (not determined) (label)
  contamination_level    estimated level (fraction); 0 when status is NO
  haplocheck_distance    distance between the major and minor haplogroups (count)
  haplocheck_coverage    sample coverage reported by haplocheck (x)
  haplocheck_status      haplocheck's own status before the support check (label)
  haplocheck_level       haplocheck's own level (fraction)
  minor_haplogroup       haplocheck's minor haplogroup (label)
  minor_heteroplasmies   minor heteroplasmies supporting it (count)
  All SKIPPED without --mt-variants. haplocheck's report: contamination/<sample>.haplocheck.txt.
EOF
}

MIN_SUPPORT=5; LOWCONF_DEPTH=2000
while (( $# )); do
  case "$1" in
    --min-support)   MIN_SUPPORT="${2:?}"; shift 2;;
    --lowconf-depth) LOWCONF_DEPTH="${2:?}"; shift 2;;
    -h|--help) usage; exit 0;;
    *)         std_opt "$@" || die "unknown option: $1"; shift "$STD_SHIFT";;
  esac
done
check_std_args

[[ "$MIN_SUPPORT" =~ ^[0-9]+$ ]] || die "--min-support must be a non-negative integer"
[[ "$LOWCONF_DEPTH" =~ ^[0-9]+$ ]] || die "--lowconf-depth must be a non-negative integer"
OUT=$(module_dir "$MODULE"); mkdir -p "$OUT"
STATUS=NA; LEVEL=NA; DIST=NA; COV=NA; QC=PASS; HSTATUS=NA; HLEVEL=NA; MINOR=NA; NMINOR=NA; LQC=PASS
if [[ -z "$MT_VARIANTS" ]]; then
  log "no --mt-variants; writing SKIPPED rows"
  QC=SKIPPED
else
  need_tools java
  [[ -s "$MITRIX_HAPLOCHECK" ]] || die "haplocheck not found: $MITRIX_HAPLOCHECK (run tools/install.sh)"
  VCF=$(readlink -f "$(core_dir)/${SAMPLE}.mt.vcf.gz")
  [[ -s "$VCF" ]] || die "clean VCF missing: $VCF (run core/mt_variants.sh first)"
  make_tmpdir
  # haplocheck writes <out> plus <out>.html into the output's directory
  ( cd "$TMP" && java -Xmx2G -cp "$MITRIX_HAPLOCHECK" genepi.haplocheck.App --raw --out "$TMP/hc.txt" "$VCF" ) \
      > "$TMP/hc.log" 2>&1 || { cat "$TMP/hc.log" >&2; die "haplocheck failed"; }
  # --raw also writes <out-stem>.raw.txt with the major/minor haplogroups and their support
  RAW="$TMP/hc.raw.txt"; [[ -s "$RAW" ]] || RAW="$TMP/hc.txt"
  tr -d '"' < "$RAW" > "$OUT/${SAMPLE}.haplocheck.txt"
  IFS=$'\t' read -r STATUS LEVEL DIST COV MINOR NMINOR < <(awk -F'\t' -v OFS='\t' '
    NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
    NR == 2 { print $c["Contamination Status"], $c["Contamination Level"], $c["Distance"], $c["Sample Coverage"],
                    ("Minor Haplogroup" in c) ? $c["Minor Haplogroup"] : "NA",
                    ("Minor Heteroplasmies Count" in c) ? $c["Minor Heteroplasmies Count"] : "NA" }' \
    "$OUT/${SAMPLE}.haplocheck.txt")
  HSTATUS=$STATUS; HLEVEL=$LEVEL
  # support check: a YES needs a minor haplogroup other than rCRS with >= MIN_SUPPORT sites
  if [[ "$STATUS" == YES ]]; then
    if [[ "$MINOR" == H2a2a1* || ! "$NMINOR" =~ ^[0-9]+$ || "$NMINOR" -lt "$MIN_SUPPORT" ]]; then
      log "haplocheck YES (level $LEVEL, minor $MINOR, $NMINOR site(s)) not supported: set to NO"
      STATUS=NO; LQC=unsupported
    fi
  fi
  if [[ -z "${STATUS:-}" ]]; then
    log "haplocheck returned no result"
    STATUS=NA; QC=FAIL
  fi
  # haplocheck writes ND (not detected) for the level of uncontaminated samples
  [[ "$STATUS" == NO ]] && LEVEL=0
  for v in LEVEL DIST COV; do
    [[ "${!v}" =~ ^[0-9.eE+-]+$ ]] || printf -v "$v" NA
  done
  if [[ "$STATUS" == YES && "$LEVEL" != NA && "$COV" != NA ]] && \
     awk -v l="$LEVEL" -v c="$COV" -v d="$LOWCONF_DEPTH" 'BEGIN { exit !(l < 0.05 && c < d) }'; then
    LQC=low_confidence
  fi
  [[ "$HLEVEL" =~ ^[0-9.eE+-]+$ ]] || HLEVEL=NA
  log "contamination $STATUS (level $LEVEL, distance $DIST, coverage $COV; haplocheck $HSTATUS, minor $MINOR with $NMINOR site(s))"
fi

{
  echo "$RESULT_HEADER"
  write_row contamination_status "$STATUS" label    "$QC"
  write_row contamination_level  "$LEVEL"  fraction "$([[ $QC == PASS ]] && echo "$LQC" || echo "$QC")"
  write_row haplocheck_distance  "$DIST"   count    "$QC"
  write_row haplocheck_coverage  "$COV"    x        "$QC"
  write_row haplocheck_status    "$HSTATUS" label   "$QC"
  write_row haplocheck_level     "$HLEVEL" fraction "$QC"
  write_row minor_haplogroup     "$MINOR"  label    "$QC"
  write_row minor_heteroplasmies "$NMINOR" count    "$QC"
} > "$OUT/${SAMPLE}.${MODULE}.tsv"
{
  write_param jar  "$MITRIX_HAPLOCHECK"
  write_param min_support   "$MIN_SUPPORT"
  write_param lowconf_depth "$LOWCONF_DEPTH"
} > "$OUT/${SAMPLE}.${MODULE}.params.tsv"
