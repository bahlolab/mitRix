#!/usr/bin/env bash
# core/mt_coverage.sh — chrM per-base coverage with mosdepth, expanded from RLE to
# one line per position. Writes <outdir>/<sample>/core/<sample>.chrM.per-base.bed.gz.

set -euo pipefail
LOG_PREFIX=mt_coverage
source "$(dirname "$(readlink -f "$0")")/../lib/common.sh"

while (( $# )); do
  std_opt "$@" || die "unknown option: $1"
  shift "$STD_SHIFT"
done
check_std_args
need_tools mosdepth zcat awk gzip
load_ref_context
make_tmpdir

PREFIX="$(core_dir)/${SAMPLE}.chrM"
log "computing chrM per-base coverage → ${PREFIX}.per-base.bed.gz"
mosdepth -t "$THREADS" -F "$FLAG" -c "$MT_CHR" "${MOSDEPTH_FASTA_OPT[@]}" "$PREFIX" "$BAM" \
  2> >(grep -v "\[W::hts_idx_load3\]" >&2) \
  || die "mosdepth failed for chrM per-base coverage"
[[ -e "${PREFIX}.per-base.bed.gz" ]] || die "chrM per-base output not found: ${PREFIX}.per-base.bed.gz"
rm -f "${PREFIX}.mosdepth.summary.txt" "${PREFIX}.mosdepth.global.dist.txt"

# Expand mosdepth RLE to one line per base position
log "expanding chrM per-base RLE → one line per position"
zcat "${PREFIX}.per-base.bed.gz" \
  | awk '{ for (p=$2; p<$3; p++) print $1"\t"p"\t"(p+1)"\t"$4 }' \
  | gzip -c > "$TMP/perbase.gz" \
  || die "RLE expansion failed (is /tmp full?)"
mv "$TMP/perbase.gz" "${PREFIX}.per-base.bed.gz"
