#!/usr/bin/env bash
# tests/test_carriers.sh — unit test of `mitrix carriers` on a synthetic VCF: het/hom
# dosage, missing genotypes, a multi-allelic record, 'chr' prefix in the list but not in
# the VCF, a right-shifted deletion in the list left-aligned with --reference-fasta,
# per-gene counts, severity sub-columns, two lists, and the cohort columns. Needs bcftools, samtools, Rscript.

set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
TMP=$(mktemp -d -t mitrix_test.XXXXXX); trap 'rm -rf "$TMP"' EXIT
cd "$TMP"
# reference contig '1': ...ACGT then a poly-A run at 21-25 (deleting one A left-aligns to 20)
printf '>1\nACGTACGTACGTACGTACGTAAAAACGTACGTACGTACGTACGT\n' > ref.fa; samtools faidx ref.fa
{
  printf '##fileformat=VCFv4.2\n##contig=<ID=1,length=44>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n'
  printf '1\t5\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t1/1\t0/0\t./.\n'        # gene A: het, hom, ref, missing
  printf '1\t10\t.\tC\tT,G\t.\tPASS\t.\tGT\t0/2\t0/1\t1/2\t0/0\n'     # multi-allelic: list has C>G
  printf '1\t20\t.\tTA\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0|1\t0/0\n'      # deletion (left-aligned)
} | bgzip -c > c.vcf.gz
tabix -p vcf c.vcf.gz
printf '# test list\nchr\tpos\tref\talt\tgene\tlabel\tinheritance\tseverity\nchr1\t5\tA\tG\tGENEA\tp.X1Y\tAD\tpathogenic\nchr1\t10\tC\tG\tGENEB\tp.X2Y\tAR\tsevere\nchr1\t24\tAA\tA\tGENEB\tc.del\tAR\tmild\nchr1\t27\tG\tC\tGENEB\tabsent\tAR\tsevere\n' > l1.tsv
printf 'chr\tpos\tref\talt\tgene\tclassification\tseverity\n1\t10\tC\tT\tGENEC\tpathogenic\tp\n1\t5\tA\tG\tMODG\tdementia risk modifier\te4\n' > l3.tsv
printf 'chr\tpos\tref\talt\tgene\n1\t10\tC\tT\tGENEC\n' > l2.tsv

"$HERE/bin/mitrix" carriers --vcf c.vcf.gz --list T1=l1.tsv --list T2=l2.tsv --list T3=l3.tsv --reference-fasta ref.fa --outdir out 2> log \
  || { cat log; exit 1; }
R=out/carriers/carriers.tsv
fail=0
expect() {  # sample metric value [qc_flag]
  local got
  got=$(awk -F'\t' -v s="$1" -v m="$2" '$1 == s && $3 == m { print $4 "\t" $6 }' "$R")
  if [[ "${got%%$'\t'*}" != "$3" || ( -n "${4:-}" && "${got#*$'\t'}" != "$4" ) ]]; then
    echo "FAIL $1 $2: got '${got//$'\t'/ }', expected '$3${4:+ $4}'"; fail=1
  else echo "ok   $1 $2 = $3${4:+ ($4)}"; fi
}
expect S1 T1_n 2 PASS
expect S1 T1 "GENEA:p.X1Y(het);GENEB:p.X2Y(het)"
expect S2 T1_GENEA 2
expect S2 T1 "GENEA:p.X1Y(hom)"
expect S3 T1_n 2
expect S3 T1_GENEB 2                      # C>G (allele 2 of the multi-allelic) + deletion
expect S3 T1 "GENEB:p.X2Y(het);GENEB:c.del(het)"
expect S3 T1_GENEB_severe 1               # severity sub-columns: GENEB has two classes
expect S3 T1_GENEB_mild 1
n=$(awk -F'\t' '$3 ~ /^T1_GENEA_/' "$R" | wc -l)
(( n == 0 )) && echo "ok   no sub-columns for GENEA (one class)" || { echo "FAIL GENEA sub-columns: $n"; fail=1; }
expect S4 T1_n 0 missing_genotypes
expect S4 T1_missing 1
expect S4 T1 none
expect S2 T2_n 1
expect S3 T2 "GENEC:1:10C>T(het)"
expect S2 T3_n 1                          # modifier (MODG, S2 hom) not counted
expect S2 T3 "GENEC:1:10C>T(het)"
expect S2 T3_MODG 2
expect S1 T3_n 0
expect S1 T3 none
m=$(awk -F'\t' '$1 == "T1" { print $2, $3, $5 }' out/carriers/match.tsv)
[[ "$m" == "4 3 3" ]] && echo "ok   T1 match: 4 listed, 3 found, 3 carriers" || { echo "FAIL T1 match: $m"; fail=1; }

# cohort columns from a minimal run directory
mkdir -p run
printf 'sample\tbam\nS1\tx.bam\nS3\tx.bam\n' > ss.tsv
for s in S1 S3; do
  mkdir -p run/$s
  printf 'step\tstatus\texit_code\tstart\tend\tversion\nref_context\tOK\t0\t.\t.\t1\n' > run/$s/$s.status.tsv
  printf 'sample\tbiomarker\tmetric\tvalue\tunit\tqc_flag\tversion\n' > run/$s/$s.tsv
done
"$HERE/bin/mitrix" cohort --samplesheet ss.tsv --outdir run --carriers "$R" 2> clog || { cat clog; exit 1; }
h=$(head -1 run/cohort.tsv | tr '\t' ' ')
[[ "$h" == *"T1_n T1_GENEA T1_GENEB T1_GENEB_mild T1_GENEB_severe T1"* && "$h" == *"T2_n T2_GENEC T2"* ]] && echo "ok   cohort columns" || { echo "FAIL cohort header: $h"; fail=1; }
v=$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i } $1 == "S3" { print $c["T1_GENEB"] }' run/cohort.tsv)
[[ "$v" == 2 ]] && echo "ok   cohort S3 T1_GENEB = 2" || { echo "FAIL cohort S3 T1_GENEB: $v"; fail=1; }

(( fail )) && { echo FAILED; exit 1; }
echo "all carriers tests passed"
