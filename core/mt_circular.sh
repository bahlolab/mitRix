#!/usr/bin/env bash
# core/mt_circular.sh — re-align chrM reads to the circularised reference chrMC
# (rCRS + first 700 bp, 17269 bp) and fold coverage back to the mt contig length.
#
# Reads crossing the linear 16569→1 seam are soft-clipped in a standard BAM, which
# lowers coverage at both ends of the D-loop amplicon. After re-alignment, positions
# 1–700 receive depth from both copies (start and seam) summed.
#
# Writes to <outdir>/<sample>/core/:
#   <sample>.chrMC.per-base.bed.gz            chrMC coverage, one line per position
#   <sample>.chrM.corrected.per-base.bed.gz   folded back to the mt contig length

set -euo pipefail
LOG_PREFIX=mt_circular
source "$(dirname "$(readlink -f "$0")")/../lib/common.sh"

while (( $# )); do
  std_opt "$@" || die "unknown option: $1"
  shift "$STD_SHIFT"
done
check_std_args
need_tools samtools mosdepth bwa zcat awk gzip
load_ref_context
make_tmpdir

CHRMC_REF="$MITRIX_HOME/refs/chrMC_rCRS.fa"
[[ -f "$CHRMC_REF" ]] || die "chrMC reference not found: $CHRMC_REF (run: bash refs/make_chrMC.sh)"
if [[ ! -f "${CHRMC_REF}.bwt" ]]; then
  log "building BWA index for chrMC (one-time)..."
  bwa index "$CHRMC_REF" 2>/dev/null || die "bwa index failed for $CHRMC_REF"
fi

# contig name and length of the circularised reference, from its FASTA
CHRMC_NAME=$(awk '/^>/ { sub(/^>/, ""); print $1; exit }' "$CHRMC_REF")
CHRMC_LEN=$(awk '!/^>/ { n += length($0) } END { print n }' "$CHRMC_REF")
(( CHRMC_LEN > MT_LEN )) \
  || die "$CHRMC_REF ($CHRMC_LEN bp) is not longer than $MT_CHR ($MT_LEN bp); is it built from this mt reference?"
log "chrMC reference: $CHRMC_NAME, $CHRMC_LEN bp ($MT_LEN + $(( CHRMC_LEN - MT_LEN )) bp appended)"

OUT=$(core_dir)
ROTATED_BAM="$TMP/${SAMPLE}.chrMC.bam"
log "re-aligning chrM reads to chrMC (paired-end mode)..."
# Extraction keeps unmapped reads (bit 4 cleared): an unmapped mate placed at its
# chrM partner's position is often a read that crosses the seam. Supplementary
# alignments are excluded on top of FLAG.
EXTRACT_FLAG=$(( (FLAG & ~4) | 2048 ))
MOSDEPTH_FLAG=$(( FLAG | 2048 ))
FQ_R1="$TMP/${SAMPLE}.R1.fq.gz"
FQ_R2="$TMP/${SAMPLE}.R2.fq.gz"
FQ_SE="$TMP/${SAMPLE}.se.fq.gz"
# extract + collate to a temp BAM (process substitution inside a pipeline fails on some HPCs)
COLLATED_BAM="$TMP/${SAMPLE}.collated.bam"
samtools view -b -F "$EXTRACT_FLAG" "${SAMTOOLS_REF_OPT[@]}" "$BAM" "$MT_CHR" \
    2> >(grep -v "^\[W::" >&2) \
  | samtools collate -O -u -@ "$THREADS" - "$TMP/collate" 2>/dev/null \
  > "$COLLATED_BAM" \
  || die "failed to collate chrM reads"
# R1, R2 and singletons (paired read whose mate is not on chrM)
samtools fastq -1 "$FQ_R1" -2 "$FQ_R2" -s "$FQ_SE" -0 /dev/null \
               -@ "$THREADS" "$COLLATED_BAM" 2>/dev/null \
  || die "failed to convert chrM reads to FASTQ"
rm -f "$COLLATED_BAM"
# bwa mem reads input in chunks of 10M bases x threads and estimates the insert
# size per chunk, so results depend on --threads unless the chunk size is fixed.
BWA_CHUNK=100000000
# Paired-end alignment keeps mate information; singletons go single-end.
# Only one SAM header may reach samtools sort.
# An empty gzip is not 0 bytes, so test for content, not size.
has_reads() { [[ -n "$(gzip -cd "$1" 2>/dev/null | head -c1)" ]]; }
{
  if has_reads "$FQ_R1"; then
    bwa mem -t "$THREADS" -K "$BWA_CHUNK" -v 1 -Y "$CHRMC_REF" \
            <(gunzip -c "$FQ_R1") <(gunzip -c "$FQ_R2") 2>/dev/null
    if has_reads "$FQ_SE"; then
      bwa mem -t "$THREADS" -K "$BWA_CHUNK" -v 1 -Y "$CHRMC_REF" \
              <(gunzip -c "$FQ_SE") 2>/dev/null | grep -v '^@'
    fi
  elif has_reads "$FQ_SE"; then
    bwa mem -t "$THREADS" -K "$BWA_CHUNK" -v 1 -Y "$CHRMC_REF" \
            <(gunzip -c "$FQ_SE") 2>/dev/null
  fi
} | samtools sort -T "$TMP/sort" -@ "$THREADS" -o "$ROTATED_BAM" \
  || die "failed to create chrMC BAM"
rm -f "$FQ_R1" "$FQ_R2" "$FQ_SE"
samtools index "$ROTATED_BAM" || die "failed to index chrMC BAM"

# chrMC per-base coverage, expanded from RLE
CHRMC_PREFIX="$OUT/${SAMPLE}.chrMC"
log "computing chrMC per-base coverage → ${CHRMC_PREFIX}.per-base.bed.gz"
mosdepth -t "$THREADS" -F "$MOSDEPTH_FLAG" -c "$CHRMC_NAME" "$CHRMC_PREFIX" "$ROTATED_BAM" 2>/dev/null \
  || die "mosdepth failed for chrMC per-base"
rm -f "${CHRMC_PREFIX}.mosdepth.summary.txt" "${CHRMC_PREFIX}.mosdepth.global.dist.txt"
zcat "${CHRMC_PREFIX}.per-base.bed.gz" \
  | awk '{ for (p=$2; p<$3; p++) print $1"\t"p"\t"(p+1)"\t"$4 }' \
  | gzip -c > "$TMP/chrMC.tmp.gz" \
  || die "RLE expansion of chrMC per-base failed (is /tmp full?)"
mv "$TMP/chrMC.tmp.gz" "${CHRMC_PREFIX}.per-base.bed.gz"

# Fold chrMC back to the mt contig length (modulo-sum)
CORRECTED="$OUT/${SAMPLE}.chrM.corrected.per-base.bed.gz"
log "building corrected chrM per-base ($CHRMC_NAME → $MT_LEN bp)..."
zcat "${CHRMC_PREFIX}.per-base.bed.gz" \
  | awk -v mt_chr="$MT_CHR" -v mt_len="$MT_LEN" '
      {
        pos = $2 + 1
        dep = $4 + 0
        if (pos <= mt_len) { d[pos] += dep }
        else               { d[pos - mt_len] += dep }
      }
      END {
        for (i = 1; i <= mt_len; i++)
          print mt_chr "\t" (i-1) "\t" i "\t" (i in d ? d[i] : 0)
      }' \
  | gzip -c > "$CORRECTED" \
  || die "failed to build corrected chrM per-base"
