#!/usr/bin/env bash
# tools/make_numt.sh — build presets/numt_alleles.tsv: mtDNA alleles found in nuclear
# mitochondrial sequences (NUMTs), from MitoHPC's RefSeq/NUMT.vcf.gz (MIT licence,
# https://github.com/dpuiu/MitoHPC; NUMT GenBank accessions in the source column).
# Heteroplasmies at these alleles below the expected NUMT level are excluded by the
# heteroplasmy module. Usage: tools/make_numt.sh [OUT]. Needs bcftools, samtools, curl.
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
OUT="${1:-$HERE/presets/numt_alleles.tsv}"
URL=https://raw.githubusercontent.com/dpuiu/MitoHPC/main/RefSeq/NUMT.vcf.gz
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
curl -sfL "$URL" -o "$TMP/numt.vcf.gz"
awk '/^>/ { next } { s = s $0 } END { s = substr(s, 1, 16569); print ">chrM"; for (i = 1; i <= length(s); i += 60) print substr(s, i, 60) }' \
  "$HERE/refs/chrMC_rCRS.fa" > "$TMP/rcrs.fa"
samtools faidx "$TMP/rcrs.fa"
# chrM (rCRS) records only (the file also has RSRS), left-normalised against rCRS
{ printf '##fileformat=VCFv4.2\n##contig=<ID=chrM,length=16569>\n##INFO=<ID=NUMT,Number=1,Type=String,Description="NUMT">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  bcftools view -H "$TMP/numt.vcf.gz" | awk -F'\t' -v OFS='\t' '$1 == "chrM" { $1 = "chrM"; print }'; } \
  | bcftools norm --no-version -f "$TMP/rcrs.fa" -c x -m -any 2> "$TMP/norm.log" \
  | bcftools query -f '%POS\t%REF\t%ALT\t%INFO/NUMT\n' | sort -k1,1n -u > "$TMP/body.tsv"
{ echo "# mtDNA alleles present in NUMT sequences (mitRix tools/make_numt.sh), from MitoHPC RefSeq/NUMT.vcf.gz"
  echo "# (MIT licence, Copyright (c) MitoHPC authors, https://github.com/dpuiu/MitoHPC). rCRS positions, left-normalised."
  printf 'pos\tref\talt\tnumt\n'; cat "$TMP/body.tsv"; } > "$OUT"
echo "make_numt: $(wc -l < "$TMP/body.tsv") alleles -> $OUT" >&2
