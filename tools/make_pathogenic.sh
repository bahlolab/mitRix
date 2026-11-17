#!/usr/bin/env bash
# tools/make_pathogenic.sh — build presets/mt_pathogenic.tsv: pathogenic mtDNA variants
# from MITOMAP (all "Cfrm" confirmed, incl. Cfrm [VUS*]) and ClinVar (Pathogenic /
# Likely pathogenic, not conflicting, with review criteria: criteria provided, expert
# panel or practice guideline; "no assertion criteria" entries are left out).
#
# Usage: tools/make_pathogenic.sh [OUT]          (default presets/mt_pathogenic.tsv)
#   MITOMAP_CODING=csv MITOMAP_RNA=csv   MITOMAP Foswiki exports (default: the copies in
#                                        tools/pathogenic/; mitomap.org blocks scripted
#                                        downloads, so export them in a browser)
#   CLINVAR=vcf.gz                       ClinVar GRCh38 VCF (default: streamed from NCBI,
#                                        chrMT only)
# Needs bcftools (with remote access for the default ClinVar), samtools, Rscript.
# Both sources are normalised against rCRS (left-aligned) so they match mitRix's
# normalised variant table. MITOMAP: CC BY 3.0; ClinVar: public domain.

set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
OUT="${1:-$HERE/presets/mt_pathogenic.tsv}"
MC="${MITOMAP_CODING:-$HERE/tools/pathogenic/mitomap_coding_control_2026-03-20.csv}"
MR="${MITOMAP_RNA:-$HERE/tools/pathogenic/mitomap_rna_2026-04-16.csv}"
CV="${CLINVAR:-https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz}"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
log() { echo "make_pathogenic: $*" >&2; }

# rCRS as contig "MT" (ClinVar naming)
awk '/^>/ { next } { s = s $0 } END { s = substr(s, 1, 16569); print ">MT"; for (i = 1; i <= length(s); i += 60) print substr(s, i, 60) }' \
  "$HERE/refs/chrMC_rCRS.fa" > "$TMP/rcrs.fa"
samtools faidx "$TMP/rcrs.fa"

# ---- MITOMAP Cfrm
Rscript --vanilla "$HERE/tools/pathogenic/mitomap_cfrm.R" "$MC" "$MR" "$HERE/refs/chrMC_rCRS.fa" "$TMP/mitomap.tsv"

# ---- ClinVar P/LP with review criteria
bcftools view -r MT "$CV" -Oz -o "$TMP/clinvar_mt.vcf.gz" 2> "$TMP/cv.log" || { cat "$TMP/cv.log" >&2; exit 1; }
CV_DATE=$(bcftools view -h "$TMP/clinvar_mt.vcf.gz" | awk -F= '$1 == "##fileDate" { print $2 }')
bcftools query -f '%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/GENEINFO\t%INFO/CLNDN\n' "$TMP/clinvar_mt.vcf.gz" \
  | awk -F'\t' -v OFS='\t' '
      BEGIN { print "pos", "ref", "alt", "gene", "disease", "clinvar_sig", "clinvar_review" }
      $4 ~ /(^|[|_\/])(Pathogenic|Likely_pathogenic)/ && $4 !~ /Conflicting/ &&
      $5 ~ /criteria_provided|reviewed_by_expert_panel|practice_guideline/ && $5 !~ /no_assertion|conflicting/ &&
      $3 != "." {
        g = $6; sub(/:.*/, "", g); if (g == ".") g = ""
        d = $7; gsub(/_/, " ", d); gsub(/\|/, "; ", d)
        print $1, $2, $3, g, d, $4, $5 }' > "$TMP/clinvar.tsv"
log "ClinVar ($CV_DATE): $(( $(wc -l < "$TMP/clinvar.tsv") - 1 )) P/LP chrMT variants with review criteria"

# ---- normalise both against rCRS (left-align), then merge on pos/ref/alt
norm() {  # in.tsv out.tsv
  { printf '##fileformat=VCFv4.2\n##contig=<ID=MT,length=16569>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    awk -F'\t' -v OFS='\t' 'NR > 1 { print "MT", $1, NR - 1, $2, $3, ".", ".", "." }' "$1"; } > "$TMP/n.vcf"
  bcftools norm --no-version -f "$TMP/rcrs.fa" -c x -m -any "$TMP/n.vcf" 2> "$TMP/n.log" \
    | bcftools query -f '%ID\t%POS\t%REF\t%ALT\n' > "$TMP/n.ids"
  grep -q REF_MISMATCH "$TMP/n.log" && log "warning: $(grep -c REF_MISMATCH "$TMP/n.log") REF mismatch(es) dropped from $(basename "$1")"
  awk -F'\t' -v OFS='\t' 'NR == FNR { p[$1] = $2; r[$1] = $3; a[$1] = $4; next }
    FNR == 1 { print; next } ((FNR - 1) in p) { $1 = p[FNR - 1]; $2 = r[FNR - 1]; $3 = a[FNR - 1]; print }' "$TMP/n.ids" "$1" > "$2"
}
norm "$TMP/mitomap.tsv" "$TMP/mitomap.n.tsv"
norm "$TMP/clinvar.tsv" "$TMP/clinvar.n.tsv"

Rscript --vanilla - "$TMP/mitomap.n.tsv" "$TMP/clinvar.n.tsv" "$OUT" "$CV_DATE" <<'EOF'
a <- commandArgs(TRUE)
cc <- c(ref = "character", alt = "character", gene = "character", disease = "character")
m <- read.delim(a[1], colClasses = cc, quote = "", na.strings = character())
c <- read.delim(a[2], colClasses = cc, quote = "", na.strings = character())
k <- function(x) paste(x$pos, x$ref, x$alt)
# one row per variant: MITOMAP duplicates (same variant listed twice) collapsed
m <- do.call(rbind, lapply(split(m, k(m)), function(x) data.frame(pos = x$pos[1], ref = x$ref[1], alt = x$alt[1],
       gene = x$gene[1], disease = paste(unique(x$disease), collapse = " / "), mitomap_status = paste(unique(x$mitomap_status), collapse = "; "))))
c <- c[!duplicated(k(c)), ]
keys <- union(k(m), k(c))
im <- match(keys, k(m)); ic <- match(keys, k(c))
pick <- function(x, y) ifelse(is.na(x) | x == "", y, x)
o <- data.frame(
  pos = as.integer(pick(as.character(m$pos[im]), as.character(c$pos[ic]))),
  ref = pick(m$ref[im], c$ref[ic]), alt = pick(m$alt[im], c$alt[ic]),
  gene = pick(m$gene[im], c$gene[ic]),
  disease = pick(m$disease[im], c$disease[ic]),
  source = ifelse(!is.na(im) & !is.na(ic), "MITOMAP;ClinVar", ifelse(!is.na(im), "MITOMAP", "ClinVar")),
  mitomap_status = ifelse(is.na(im), "", m$mitomap_status[im]),
  clinvar_sig = ifelse(is.na(ic), "", c$clinvar_sig[ic]),
  clinvar_review = ifelse(is.na(ic), "", c$clinvar_review[ic]))
o$variant <- paste0("m.", o$pos, o$ref, ">", o$alt)
o <- o[order(o$pos, o$ref, o$alt), c("pos", "ref", "alt", "variant", "gene", "disease", "source", "mitomap_status", "clinvar_sig", "clinvar_review")]
con <- file(a[3], "w")
writeLines(c("# mitRix pathogenic mtDNA variants (tools/make_pathogenic.sh): MITOMAP Cfrm (CC BY 3.0,",
             "# https://www.mitomap.org, PMID 25489354; exports 2026-03-20 coding/control, 2026-04-16 rRNA/tRNA)",
             paste0("# + ClinVar P/LP with review criteria (", a[4], "). rCRS positions, left-normalised.")), con)
write.table(o, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con)
message(sprintf("make_pathogenic: %d variants (%d both, %d MITOMAP only, %d ClinVar only) -> %s",
  nrow(o), sum(o$source == "MITOMAP;ClinVar"), sum(o$source == "MITOMAP"), sum(o$source == "ClinVar"), a[3]))
EOF
