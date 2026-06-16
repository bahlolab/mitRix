#!/usr/bin/env bash
# make_chrMC.sh — one-time setup: create the circularised rCRS reference
#                 used by mitRix --rotate-chrm for accurate D-loop coverage.
#
# Strategy (based on MitoHPC; EXT=500 rather than MitoHPC default 300):
#   chrMC = rCRS (16,569 bp) + first 700 bp of rCRS appended = 17,269 bp
#   The artificial linear break moves from 16569/1 into the coding region,
#   so the D-loop amplicon (orig 16536–16569 + 1–34) becomes the contiguous
#   interval chrMC:16536–16603.
#   EXT=700 (≈ 2 × mean insert) captures seam-crossing read pairs with
#   inserts up to ~850 bp, safely past the 99th percentile for typical WGS.
#
# Run once from the mitRix directory:
#   bash refs/make_chrMC.sh
#
# Requirements: curl, bwa, awk

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OUT="$SCRIPT_DIR/chrMC_rCRS.fa"
EXT=700      # bp to append; rule of thumb: 2 × mean insert size (here ~355 bp)
             # captures inserts up to ~850 bp, safely beyond the 99th percentile
             # for WGS with mean insert ~350 bp (MitoHPC default is 300)
NCBI_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=fasta&retmode=text"

# ----- dependency checks -----
for cmd in curl bwa awk; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "error: $cmd not found in PATH"; exit 2; }
done

# ----- skip if already built -----
if [[ -f "$OUT" && -f "${OUT}.bwt" ]]; then
  echo "Already built: $OUT (with BWA index). Nothing to do."
  exit 0
fi

# ----- download rCRS -----
TMP_FA="$SCRIPT_DIR/.rCRS_tmp_$$.fa"
trap 'rm -f "$TMP_FA"' EXIT

echo "Downloading rCRS (NC_012920.1) from NCBI..."
curl -fsSL "$NCBI_URL" > "$TMP_FA"
[[ -s "$TMP_FA" ]] || { echo "error: download failed or empty"; exit 2; }

# ----- build chrMC.fa -----
echo "Creating circularised reference (appending first ${EXT} bp to rCRS)..."
awk -v ext="$EXT" -v out="$OUT" '
  /^>/ { next }
  { seq = seq $0 }
  END {
    gsub(/[[:space:]]/, "", seq)
    n = length(seq)
    if (n != 16569) {
      printf "warning: expected 16569 bp, got %d bp\n", n > "/dev/stderr"
    }
    full = seq substr(seq, 1, ext)
    print ">chrMC" > out
    for (i = 1; i <= length(full); i += 60)
      print substr(full, i, 60) > out
    printf "rCRS:  %d bp\nchrMC: %d bp (%d bp appended)\n", n, length(full), ext
  }
' "$TMP_FA"

# ----- BWA index -----
echo "Building BWA index for $OUT ..."
bwa index "$OUT" 2>/dev/null

echo ""
echo "Done. Files in $SCRIPT_DIR/:"
ls -lh "$OUT" "${OUT}".{amb,ann,bwt,pac,sa} 2>/dev/null | awk '{print "  "$0}'
echo ""
echo "D-loop amplicon coords in chrMC (used automatically with --ratio DLOOP/ND1):"
echo "  original: chrM:16536-16569 + chrM:1-34  (split across break)"
echo "  chrMC:    chrMC:16536-16603              (contiguous 68 bp)"
