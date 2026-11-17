#!/usr/bin/env bash
# core/ref_context.sh — read the BAM/CRAM header once per sample.
# Writes <outdir>/<sample>/core/<sample>.ref_context.tsv (key/value) and
# <sample>.contigs.tsv (name, length).

set -euo pipefail
LOG_PREFIX=ref_context
source "$(dirname "$(readlink -f "$0")")/../lib/common.sh"

while (( $# )); do
  std_opt "$@" || die "unknown option: $1"
  shift "$STD_SHIFT"
done
check_std_args
need_tools samtools awk

OUT=$(core_dir); mkdir -p "$OUT"

case "$BAM" in
  *.cram) [[ -e "${BAM}.crai" || -e "${BAM%.cram}.crai" ]] || log "warning: no .crai index found beside $BAM";;
  *.bam)  [[ -e "${BAM}.bai" || -e "${BAM%.bam}.bai" || -e "${BAM}.csi" ]] \
            || log "warning: no .bai/.csi index found beside $BAM";;
esac

HDR=$(samtools view -H "${SAMTOOLS_REF_OPT[@]}" "$BAM") || die "failed to read header from $BAM"
SQ=$(printf '%s\n' "$HDR" | awk 'BEGIN{OFS="\t"}
  /^@SQ/ {
    sn=""; ln="";
    for (i=2;i<=NF;i++){
      if ($i ~ /^SN:/) sn=substr($i,4);
      else if ($i ~ /^LN:/) ln=substr($i,4);
    }
    if (sn != "") print sn, ln;
  }')
printf '%s\n' "$SQ" > "$OUT/${SAMPLE}.contigs.tsv"

# here-string, not a pipe: awk exits early and would SIGPIPE printf under pipefail
contig_line() { awk -v c="$1" '$1==c{print; exit}' <<< "$SQ"; }

# mtDNA contig: known names, else any contig of rCRS length
MT_CHR=""; MT_LEN=""
for cand in chrM MT chrMT NC_012920.1 NC_012920; do
  line=$(contig_line "$cand")
  if [[ -n "$line" ]]; then MT_CHR=${line%%$'\t'*}; MT_LEN=${line#*$'\t'}; break; fi
done
if [[ -z "$MT_CHR" ]]; then
  line=$(awk '$2==16569{print; exit}' <<< "$SQ")
  [[ -n "$line" ]] && { MT_CHR=${line%%$'\t'*}; MT_LEN=${line#*$'\t'}; }
fi
[[ -n "$MT_CHR" ]] || die "no mtDNA contig (chrM/MT/chrMT/NC_012920.1 or LN=16569) found in BAM header"

if [[ "$MT_LEN" != "16569" ]]; then
  log "warning: mtDNA contig '$MT_CHR' has length $MT_LEN, not 16569 (rCRS). The built-in default coordinates assume rCRS — verify with --target/--reference or --target-bed/--reference-bed."
else
  log "mtDNA contig: $MT_CHR (rCRS, 16569 bp)"
fi

# chr15 (for B2M) and nuclear build from its length
CHR15=""; CHR15_LEN=""
for cand in chr15 15; do
  line=$(contig_line "$cand")
  if [[ -n "$line" ]]; then CHR15=${line%%$'\t'*}; CHR15_LEN=${line#*$'\t'}; break; fi
done
case "$CHR15_LEN" in
  102531392) BUILD=hg19;;
  101991189) BUILD=hg38;;
  *)         BUILD=unknown;;
esac

{
  write_param mt_contig    "$MT_CHR"
  write_param mt_length    "$MT_LEN"
  write_param chr15_contig "$CHR15"
  write_param chr15_length "$CHR15_LEN"
  write_param build        "$BUILD"
} > "$OUT/${SAMPLE}.ref_context.tsv"
