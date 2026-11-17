#!/usr/bin/env bash
# tools/make_annotation.sh — build presets/mt_annotation.tsv.gz, the rCRS SNV annotation
# used by the heteroplasmy module (same for every variant caller):
#
#   pos  ref  alt  gene  region  consequence  vep_consequence  mlc_score
#
# Sources (downloaded and checksum-verified):
#   - Lake et al. 2024, Nature 635:390, "Quantifying constraint in the human mitochondrial
#     genome", Supplementary Dataset 7: VEP consequence and MLC score for all 49,704
#     possible rCRS SNVs (open access).
#   - leklab/mitochondrial_constraint (MIT): required_files/databases/mitomap_genome_loci.txt
#     (MITOMAP locus boundaries), for gene and region.
#
# gene:        MT-* tRNA, rRNA (MT-RNR1/2) or protein-coding gene; overlapping genes joined
#              with '|' (e.g. MT-ATP8|MT-ATP6); '.' outside genes.
# region:      protein_coding, rRNA, tRNA, control_region (MITOMAP MT-CR, 16024-576) or
#              non_coding (MT-NC1..10).
# consequence: most severe of the VEP terms: stop_gained > start_lost > stop_lost >
#              missense > synonymous (incl. stop_retained) > incomplete_terminal_codon >
#              non_coding (non_coding_transcript_exon) > intergenic.
#
# usage: tools/make_annotation.sh [OUT]   (default: presets/mt_annotation.tsv.gz)

set -euo pipefail
LOG_PREFIX=make_annotation
source "$(dirname "$(readlink -f "$0")")/../lib/common.sh"

OUT=${1:-$MITRIX_HOME/presets/mt_annotation.tsv.gz}
LAKE_URL='https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-08048-x/MediaObjects/41586_2024_8048_MOESM3_ESM.zip'
LAKE_SHA256=c6b338f5fc7226084c03c154d7946a3d3a87a6aa63a54c756ce52948ebe278eb
LOCI_URL=https://raw.githubusercontent.com/leklab/mitochondrial_constraint/main/required_files/databases/mitomap_genome_loci.txt
LOCI_SHA256=72df0617f51d805eb11cc424eff6d1e82d10fb052e1b8f7f0d2a68df9a9f8eec

need_tools curl unzip sha256sum awk gzip
make_tmpdir

fetch() {  # url sha256 out
  log "downloading $1"
  curl -fsSL -A "Mozilla/5.0" -o "$3" "$1" || die "download failed: $1"
  echo "$2  $3" | sha256sum -c --quiet - || die "checksum mismatch: $1"
}
fetch "$LAKE_URL" "$LAKE_SHA256" "$TMP/lake.zip"
fetch "$LOCI_URL" "$LOCI_SHA256" "$TMP/loci.txt"
unzip -p "$TMP/lake.zip" supplementary_dataset_7.tsv > "$TMP/sd7.tsv"

awk -F'\t' -v OFS='\t' '
  # loci: genes and control region; wrap-around loci (end < start) cross position 1
  FNR == NR {
    if (FNR == 1) next
    name = $1; s = $2 + 0; e = $3 + 0
    if (name ~ /^MT-T[A-Z][0-9]?$/)                          type = "tRNA"
    else if (name == "MT-RNR1" || name == "MT-RNR2")          type = "rRNA"
    else if (name ~ /^MT-(ND[1-6]|ND4L|CO[1-3]|ATP[68]|CYB)$/) type = "protein_coding"
    else if (name == "MT-CR")                                 type = "control_region"
    else if (name ~ /^MT-NC[0-9]+$/)                          type = "non_coding"
    else next
    for (p = s; ; p = (p % 16569) + 1) {
      if (type == "control_region" || type == "non_coding") reg0[p] = type
      else { gene[p] = (p in gene) ? gene[p] "|" name : name; reg[p] = type }
      if (p == e) break
    }
    next
  }
  FNR == 1 { print "pos", "ref", "alt", "gene", "region", "consequence", "vep_consequence", "mlc_score"; next }
  {
    p = $1 + 0; v = $4
    if      (v ~ /stop_gained/)                c = "stop_gained"
    else if (v ~ /start_lost/)                 c = "start_lost"
    else if (v ~ /stop_lost/)                  c = "stop_lost"
    else if (v ~ /missense/)                   c = "missense"
    else if (v ~ /synonymous|stop_retained/)   c = "synonymous"
    else if (v ~ /incomplete_terminal_codon/)  c = "incomplete_terminal_codon"
    else if (v ~ /non_coding_transcript_exon/) c = "non_coding"
    else if (v ~ /intergenic/)                 c = "intergenic"
    else { print "unknown consequence: " v > "/dev/stderr"; exit 1 }
    r = (p in reg) ? reg[p] : ((p in reg0) ? reg0[p] : "non_coding")
    print p, $2, $3, ((p in gene) ? gene[p] : "."), r, c, v, $5
  }
' "$TMP/loci.txt" "$TMP/sd7.tsv" > "$TMP/annot.tsv"

n=$(( $(wc -l < "$TMP/annot.tsv") - 1 ))
(( n == 49704 )) || die "expected 49704 SNVs, got $n"
gzip -n -9 -c "$TMP/annot.tsv" > "$OUT"
log "wrote $OUT ($n SNVs)"
