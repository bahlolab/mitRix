#!/usr/bin/env bash
# tests/test_heteroplasmy.sh — unit test of biomarkers/heteroplasmy/run.R on a synthetic
# variant table: strand/substitution mapping for all 12 SNV types, the replication-origin
# boundary (chrM:16172-210), the low/high HF split, --min-alt-reads, blocklist, excluded
# filters, the low_depth flag and the NUMT ceiling (fixed and nuclear-depth scaled). Needs Rscript. Usage: tests/test_heteroplasmy.sh

set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
TMP=$(mktemp -d -t mitrix_test.XXXXXX); trap 'rm -rf "$TMP"' EXIT

row() { printf 'T1\t%s\t%s\t%s\t%s\t%s\t%s\t%s\ttest\t%s\n' "$@"; }
{
  printf 'sample\tpos\tref\talt\thf\tdepth\talt_reads\tfilter\tcaller\ttype\n'
  # the 4 transitions and 8 transversions (rCRS REF>ALT)
  row 1000 C T 0.5 1000 500 PASS SNV     # L C>T
  row 1001 A G 0.5 1000 500 PASS SNV     # L A>G  age
  row 1002 G A 0.5 1000 500 PASS SNV     # H C>T  age
  row 1003 T C 0.5 1000 500 PASS SNV     # H A>G  age
  row 1004 C A 0.5 1000 500 PASS SNV
  row 1005 C G 0.5 1000 500 PASS SNV
  row 1006 A C 0.5 1000 500 PASS SNV
  row 1007 A T 0.5 1000 500 PASS SNV
  row 1008 G T 0.5 1000 500 PASS SNV     # H C>A
  row 1009 G C 0.5 1000 500 PASS SNV     # H C>G
  row 1010 T G 0.5 1000 500 PASS SNV     # H A>C
  row 1011 T A 0.5 1000 500 PASS SNV     # H A>T
  # replication origin 16172-210 (wraps): H A>G inside is not age-accumulating
  row 16171 T C 0.5 1000 500 PASS SNV    # age
  row 16172 T C 0.5 1000 500 PASS SNV
  row 210   T C 0.5 1000 500 PASS SNV
  row 211   T C 0.5 1000 500 PASS SNV    # age
  # HF split at 0.2: low <= 0.2 < high
  row 2000 C T 0.2       1000 200 PASS SNV
  row 2001 C T 0.2000001 1000 200 PASS SNV
  # min ALT reads 10
  row 3000 C T 0.05 199 10 PASS SNV
  row 3001 C T 0.05 180 9  PASS SNV      # removed
  # removed: blocklisted, excluded filter, outside 3-95%
  row 302  A C 0.5  1000 500 PASS SNV
  row 4000 C T 0.5  1000 500 "clustered_events;strand_bias" SNV
  row 4001 C T 0.02 1000 20  PASS SNV
  row 4002 C T 0.96 1000 960 PASS SNV
  # indel; real annotated SNVs for MSS (16360 C>T: 0.13392; 14569 G>A: synonymous, 0)
  row 5000 A AT 0.5 1000 500 PASS INDEL
  row 16360 C T 0.4 1000 400 PASS SNV
  row 3243  A G 0.1 1000 100 PASS SNV   # MT-TL1 tRNA
  row 3308  T C 0.5 1000 500 PASS SNV   # MT-ND1 start codon: start_lost, MLC 0.7
} > "$TMP/variants.tsv"
# coverage: median depth 300 -> detection limit 10/300 = 0.0333 > 0.03 -> low_depth
printf 'chrM\t0\t16569\t300\n' | gzip -c > "$TMP/cov.bed.gz"

Rscript --vanilla "$HERE/biomarkers/heteroplasmy/run.R" "$TMP/variants.tsv" "$TMP/cov.bed.gz" \
  "$HERE/presets/heteroplasmy_blocklist.tsv" "$HERE/presets/mt_annotation.tsv.gz" \
  0.03 0.95 0.2 10 strand_bias T1 test "$TMP/res.tsv" "$TMP/var.tsv" none 2 0.10 NA 2> "$TMP/log"

fail=0
expect() {  # metric value [qc_flag]
  local got
  got=$(awk -F'\t' -v m="$1" '$3 == m { print $4 "\t" $6 }' "$TMP/res.tsv")
  if [[ "${got%%$'\t'*}" != "$2" || ( -n "${3:-}" && "${got#*$'\t'}" != "$3" ) ]]; then
    echo "FAIL $1: expected $2${3:+ ($3)}, got ${got/$'\t'/ (})${got:+)}"; fail=1
  else echo "ok   $1 = $2${3:+ ($3)}"; fi
}
expect_var() {  # pos column value
  local got
  got=$(awk -F'\t' -v p="$1" -v c="$2" 'NR == 1 { for (i = 1; i <= NF; i++) k[$i] = i; next } $2 == p { print $k[c] }' "$TMP/var.tsv")
  if [[ "$got" != "$3" ]]; then echo "FAIL $1 $2: expected $3, got '$got'"; fail=1
  else echo "ok   $1 $2 = $3"; fi
}

# kept: 12 types + 4 origin + 2 split + 3000 + indel + 16360 + 3243 + 3308 = 23
expect heteroplasmy_count 23 PASS
expect heteroplasmy_low_count 3 low_depth          # 2000, 3000, 3243
expect heteroplasmy_high_count 20 PASS
expect snv_count 22
expect indel_count 1
expect L_CtoT_count 5                              # 1000 2000 2001 3000 16360
expect L_AtoG_count 2                              # 1001 3243
expect H_CtoT_count 1                              # 1002
expect H_AtoG_count 6                              # 1003 16171 16172 210 211 3308
expect transversion_count 8
expect age_accumulating_count 7                    # 1001 1002 1003 16171 211 3243 3308
expect age_accumulating_low_count 1 low_depth      # 3243
expect hf_detection_limit 0.03333333333 PASS
expect mt_median_depth 300
# MSS = sum of the variant table's MLC scores (some synthetic SNVs match rCRS by chance);
# MSS_missing_count = rows without one
expect MSS "$(awk -F'\t' 'NR > 1 && $NF != "." { s += $NF } END { printf "%.10g", s }' "$TMP/var.tsv")"
expect MSS_missing_count "$(awk -F'\t' 'NR > 1 && $NF == "." { n++ } END { print n + 0 }' "$TMP/var.tsv")"
expect_var 3243 mlc_score 0.85129
expect_var 3308 mlc_score 0.7
for p in 1000:L:C\>T 1001:L:A\>G 1002:H:C\>T 1003:H:A\>G 1004:L:C\>A 1005:L:C\>G 1006:L:A\>C \
         1007:L:A\>T 1008:H:C\>A 1009:H:C\>G 1010:H:A\>C 1011:H:A\>T; do
  IFS=: read -r pos st sub <<< "$p"
  expect_var "$pos" strand "$st"; expect_var "$pos" substitution "$sub"
done
expect_var 16171 age_accumulating yes
expect_var 16172 age_accumulating no
expect_var 210 age_accumulating no
expect_var 211 age_accumulating yes
expect_var 2000 hf_bin low
expect_var 2001 hf_bin high
expect_var 16360 mlc_score 0.13392
expect_var 3243 gene MT-TL1
expect_var 3243 region tRNA
expect_var 3308 consequence start_lost
expect_var 5000 strand .
for p in 3001 302 4000 4001 4002; do
  if awk -F'\t' -v p="$p" '$2 == p { f = 1 } END { exit !f }' "$TMP/var.tsv"; then
    echo "FAIL $p should be removed"; fail=1
  else echo "ok   $p removed"; fi
done

# NUMT-like calls: toy list with 3000 C>T (HF 0.05) and 16360 C>T (HF 0.4); chrM depth 300
printf '# toy\npos\tref\talt\tnumt\n3000\tC\tT\tX\n16360\tC\tT\tX\n' > "$TMP/numt.tsv"
numt_run() {  # nuclear_depth -> res2.tsv
  Rscript --vanilla "$HERE/biomarkers/heteroplasmy/run.R" "$TMP/variants.tsv" "$TMP/cov.bed.gz" \
    "$HERE/presets/heteroplasmy_blocklist.tsv" "$HERE/presets/mt_annotation.tsv.gz" \
    0.03 0.95 0.2 10 strand_bias T1 test "$TMP/res2.tsv" "$TMP/var2.tsv" "$TMP/numt.tsv" 2 0.10 "$1" 2>> "$TMP/log"
}
expect2() {  # label metric value
  local got; got=$(awk -F'\t' -v m="$2" '$3 == m { print $4 }' "$TMP/res2.tsv")
  [[ "$got" == "$3" ]] && echo "ok   $1 $2 = $3" || { echo "FAIL $1 $2: expected $3, got $got"; fail=1; }
}
numt_run NA     # fixed ceiling 0.10: 3000 (0.05) out, 16360 (0.4) kept
expect2 "NUMT fixed" numt_excluded_count 1; expect2 "NUMT fixed" numt_hf_threshold 0.1; expect2 "NUMT fixed" heteroplasmy_count 22
numt_run 30     # 2 x 30 / 300 = 0.2: same
expect2 "NUMT nuc30" numt_hf_threshold 0.2; expect2 "NUMT nuc30" numt_excluded_count 1
numt_run 100    # 2 x 100 / 300 = 0.667: both out
expect2 "NUMT nuc100" numt_excluded_count 2; expect2 "NUMT nuc100" heteroplasmy_count 21
if (( fail )); then echo "FAILED"; cat "$TMP/log"; exit 1; fi
echo "all heteroplasmy tests passed"
