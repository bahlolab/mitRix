#!/usr/bin/env bash
# tools/from_perbase.sh — coverage biomarkers from precomputed mitRix per-base files
# (e.g. an earlier mitRix run on Terra) without the BAM/CRAM.
#
# Input per sample: <S>.chrM.per-base.bed.gz and <S>.chrM.corrected.per-base.bed.gz
# (mosdepth per-base; corrected = circularised chrMC folded back to the mt contig),
# as a directory or a <S>.perbase.tar.gz. Runs the modules that only need these files:
#   damage     ND4/ND1, ND5/ND1 (major-arc deletion burden)   <- chrM.per-base
#   d_loop_7s  D-loop/ND1, 7S, flap/7S, flap/primer           <- chrM.corrected.per-base
#   coverage_shape  GC bias, replication gradient, deletion step  <- chrM.corrected.per-base
# Not available without the BAM: copy_number (needs nuclear B2M depth).
# Assumes an rCRS mt contig (16569 bp); region coordinates are the same in GRCh37/38.
#
# Usage: tools/from_perbase.sh --sample S --perbase DIR|FILE.tar.gz --outdir OUT [--stat mean|median]
#        [--modules damage,d_loop_7s] [--mt-contig chrM] [--build hg38]
set -euo pipefail
MITRIX_HOME="$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)"
SAMPLE=""; PB=""; OUTDIR=""; STAT=mean; MODULES="damage,d_loop_7s"; MTC=chrM; BUILD=hg38
while (( $# )); do
  case "$1" in
    --sample) SAMPLE=$2; shift 2;; --perbase) PB=$2; shift 2;; --outdir) OUTDIR=$2; shift 2;;
    --stat) STAT=$2; shift 2;; --modules) MODULES=$2; shift 2;;
    --mt-contig) MTC=$2; shift 2;; --build) BUILD=$2; shift 2;;
    -h|--help) sed -n '2,16p' "$0"; exit 0;;
    *) echo "unknown option: $1" >&2; exit 1;;
  esac
done
[[ -n "$SAMPLE" && -n "$PB" && -n "$OUTDIR" ]] || { echo "need --sample, --perbase, --outdir" >&2; exit 1; }
CORE="$OUTDIR/$SAMPLE/core"; mkdir -p "$CORE"
if [[ -d "$PB" ]]; then
  cp "$PB"/"$SAMPLE".chrM*.per-base.bed.gz "$CORE"/
else
  tar -xzf "$PB" -C "$CORE" --wildcards --no-anchored "$SAMPLE.chrM.per-base.bed.gz" "$SAMPLE.chrM.corrected.per-base.bed.gz"
fi
for f in chrM.per-base chrM.corrected.per-base; do
  [[ -s "$CORE/$SAMPLE.$f.bed.gz" ]] || { echo "missing $SAMPLE.$f.bed.gz in $PB" >&2; exit 1; }
done
# mt contig length from the per-base file (last interval end)
MTLEN=$(zcat "$CORE/$SAMPLE.chrM.per-base.bed.gz" | awk -v c="$MTC" '$1 == c { e = $3 } END { print e + 0 }')
(( MTLEN > 16000 && MTLEN < 17000 )) || { echo "$SAMPLE: unexpected $MTC length $MTLEN (rCRS expected)" >&2; exit 1; }
printf 'mt_contig\t%s\nmt_length\t%s\nchr15_contig\t\nchr15_length\t\nbuild\t%s\nsource\tfrom_perbase\n' \
  "$MTC" "$MTLEN" "$BUILD" > "$CORE/$SAMPLE.ref_context.tsv"
printf '%s\t%s\n' "$MTC" "$MTLEN" > "$CORE/$SAMPLE.contigs.tsv"
IFS=, read -r -a MODS <<< "$MODULES"
for m in "${MODS[@]}"; do
  case "$m" in damage|d_loop_7s|coverage_shape) ;; *) echo "module $m needs the BAM; only damage, d_loop_7s, coverage_shape" >&2; exit 1;; esac
  # --mt-variants /dev/null satisfies the input check; these modules do not read it
  bash "$MITRIX_HOME/biomarkers/$m/run.sh" --sample "$SAMPLE" --outdir "$OUTDIR" --stat "$STAT" --mt-variants /dev/null
done
