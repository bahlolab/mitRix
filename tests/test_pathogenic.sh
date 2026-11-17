#!/usr/bin/env bash
# tests/test_pathogenic.sh — unit test of biomarkers/pathogenic on a synthetic variant
# table: default-list matches (homoplasmic/heteroplasmic), HF and ALT-read thresholds,
# excluded filters, a user list with a right-shifted insertion (normalised to match the
# left-aligned call), list names, and low-depth positions. Needs Rscript, bcftools.

set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
TMP=$(mktemp -d -t mitrix_test.XXXXXX); trap 'rm -rf "$TMP"' EXIT
C="$TMP/out/T1/core"; mkdir -p "$C"
touch "$TMP/T1.bam" "$TMP/calls.vcf"

row() { printf 'T1\t%s\t%s\t%s\t%s\t%s\t%s\t%s\ttest\t%s\n' "$@"; }
{
  printf 'sample\tpos\tref\talt\thf\tdepth\talt_reads\tfilter\tcaller\ttype\n'
  row 3243  A G  0.12 2000 240  PASS        SNV    # detected, heteroplasmic
  row 11778 G A  1    2000 2000 PASS        SNV    # detected, homoplasmic
  row 8344  A G  0.02 2000 40   PASS        SNV    # HF < 0.03
  row 3460  G A  0.5  16   8    PASS        SNV    # 8 ALT reads < 10
  row 14484 T C  0.5  2000 1000 strand_bias SNV    # excluded filter
  row 302   A AC 0.81 2000 1620 PASS        INDEL  # user list: 309 C>CC left-aligns to 302 A>AC
  row 263   A G  1    2000 .    PASS        SNV    # no ALT read count: HF x depth
} > "$C/T1.variants.tsv"
printf 'caller\ttest\nexclude_filters\tstrand_bias\n' > "$C/T1.variants.params.tsv"
# depth 2000 everywhere except 1501-1600 at depth 100 (< 10 / 0.03 = 333)
printf 'chrM\t0\t1500\t2000\nchrM\t1500\t1600\t100\nchrM\t1600\t16569\t2000\n' | gzip -c > "$C/T1.chrM.per-base.bed.gz"
printf '# user list\npos\tref\talt\tgene\tdisease\n309\tC\tCC\tMT-CR\tD310 test\n263\tA\tG\tMT-CR\tmarker\n' > "$TMP/user.tsv"

"$HERE/biomarkers/pathogenic/run.sh" --sample T1 --bam "$TMP/T1.bam" --outdir "$TMP/out" \
  --mt-variants "$TMP/calls.vcf" --list PD_test="$TMP/user.tsv" 2> "$TMP/log" || { cat "$TMP/log"; exit 1; }
R="$TMP/out/T1/pathogenic/T1.pathogenic.tsv"; V="$TMP/out/T1/pathogenic/T1.pathogenic.variants.tsv"

fail=0
expect() {  # metric value [qc_flag]
  local got
  got=$(awk -F'\t' -v m="$1" '$3 == m { print $4 "\t" $6 }' "$R")
  if [[ "${got%%$'\t'*}" != "$2" || ( -n "${3:-}" && "${got#*$'\t'}" != "$3" ) ]]; then
    echo "FAIL $1: got '${got//$'\t'/ }', expected '$2${3:+ $3}'"; fail=1
  else echo "ok   $1 = $2${3:+ ($3)}"; fi
}
n_low=$(awk -F'\t' '!/^#/ && $1 != "pos" && $1 > 1500 && $1 <= 1600' "$HERE/presets/mt_pathogenic.tsv" | wc -l)
expect pathogenic_count 2 low_depth
expect pathogenic_variants "m.3243A>G(0.12);m.11778G>A(1)"
expect pathogenic_max_hf 1
expect pathogenic_low_depth_sites "$n_low"
expect PD_test_count 2 PASS
expect PD_test_variants "m.263A>G(1);m.302A>AC(0.81)"
expect PD_test_low_depth_sites 0
p=$(awk -F'\t' '$3 == "m.3243A>G" { print $11 }' "$V"); [[ "$p" == heteroplasmic ]] && echo "ok   3243 heteroplasmic" || { echo "FAIL 3243 plasmy: $p"; fail=1; }
p=$(awk -F'\t' '$3 == "m.11778G>A" { print $11 }' "$V"); [[ "$p" == homoplasmic ]] && echo "ok   11778 homoplasmic" || { echo "FAIL 11778 plasmy: $p"; fail=1; }
d=$(awk -F'\t' '$3 == "m.302A>AC" { print $13 }' "$V"); [[ "$d" == "D310 test" ]] && echo "ok   user-list metadata kept" || { echo "FAIL user disease: $d"; fail=1; }

# bad list name and SKIPPED without calls
"$HERE/biomarkers/pathogenic/run.sh" --sample T1 --bam "$TMP/T1.bam" --outdir "$TMP/out" \
  --mt-variants "$TMP/calls.vcf" --list pathogenic="$TMP/user.tsv" 2> /dev/null \
  && { echo "FAIL list name 'pathogenic' accepted"; fail=1; } || echo "ok   list name 'pathogenic' rejected"
"$HERE/biomarkers/pathogenic/run.sh" --sample T2 --bam "$TMP/T1.bam" --outdir "$TMP/out" 2> /dev/null
s=$(awk -F'\t' 'NR > 1 { print $6 }' "$TMP/out/T2/pathogenic/T2.pathogenic.tsv" | sort -u)
[[ "$s" == SKIPPED ]] && echo "ok   SKIPPED without --mt-variants" || { echo "FAIL skip: $s"; fail=1; }

(( fail )) && { echo FAILED; exit 1; }
echo "all pathogenic tests passed"
