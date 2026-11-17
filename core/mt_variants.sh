#!/usr/bin/env bash
# core/mt_variants.sh — read mtDNA variant calls from any supported caller, validate and
# normalise them, and write:
#
#   core/<sample>.variants.tsv   standard variant table
#       sample  pos  ref  alt  hf  depth  alt_reads  filter  caller  type
#   core/<sample>.mt.vcf.gz      clean single-sample VCF for HaploGrep3 and haplocheck
#                                (chrM, rCRS; GT 1/1 if HF >= 0.9 else 0/1; FORMAT AF, DP;
#                                records whose FILTER has an excluded filter are left out)
#
# Inputs (--mt-variants; --caller overrides detection):
#   mitohpc / mutect2  VCF; HF = FORMAT/AF (MitoHPC writes homoplasmies as GT 0/1, so GT
#                      is not used); MitoHPC's concatenated VCF is repaired; ALT == REF
#                      rows (reference-allele fraction) dropped
#   mity               VCF; HF = FORMAT/VAF
#   vcf                any other VCF with FORMAT/AF or FORMAT/VAF
#   mtdnaserver2       mtDNA-Server 2 variants.annotated.txt (TSV); HF = VariantLevel,
#                      depth = Coverage
# Multi-sample files: the sample is --variants-sample, else --sample or the BAM/CRAM name
# (a file name such as HG00096.chrM.bam also matches HG00096.chrM).
# hf: heteroplasmy fraction of the ALT allele; type: SNV, MNV or INDEL;
# alt_reads = round(hf * depth). Without --mt-variants, does nothing.

set -euo pipefail
LOG_PREFIX=mt_variants
source "$(dirname "$(readlink -f "$0")")/../lib/common.sh"

while (( $# )); do
  std_opt "$@" || die "unknown option: $1"
  shift "$STD_SHIFT"
done
check_std_args
if [[ -z "$MT_VARIANTS" ]]; then
  log "no --mt-variants; skipping (heteroplasmy, haplogroup and contamination will be SKIPPED)"
  exit 0
fi
need_tools bcftools awk
load_ref_context
make_tmpdir
OUT=$(core_dir); mkdir -p "$OUT"

if gzip -t "$MT_VARIANTS" 2>/dev/null; then CAT=(gzip -cd); else CAT=(cat); fi
"${CAT[@]}" "$MT_VARIANTS" > "$TMP/input"

# ---- sample-name matching
strip_ext() { local n=${1##*/}; n=${n%.bam}; n=${n%.cram}; n=${n%.sam}; printf '%s' "$n"; }
if [[ -n "$VARIANTS_SAMPLE" ]]; then TARGETS=("$VARIANTS_SAMPLE")
else TARGETS=("$SAMPLE" "$(strip_ext "$BAM")"); fi
# exact name or name without .bam/.cram first; then a unique "<target>.<suffix>" name
# (MitoHPC's HG00096.final, mtDNA-Server 2's HG00096.chrM.bam for sample HG00096)
pick_sample() {  # name... -> prints the matching name; 1 if none
  local n t hits
  for t in "${TARGETS[@]}"; do
    for n in "$@"; do
      if [[ "$n" == "$t" || "$(strip_ext "$n")" == "$t" ]]; then printf '%s' "$n"; return 0; fi
    done
  done
  for t in "${TARGETS[@]}"; do
    hits=()
    for n in "$@"; do [[ "$n" == "$t".* ]] && hits+=("$n"); done
    if (( ${#hits[@]} == 1 )); then printf '%s' "${hits[0]}"; return 0; fi
  done
  return 1
}
no_match() {  # names...
  die "sample '${TARGETS[*]}' not found in $MT_VARIANTS (samples: $(printf '%s ' "$@" | cut -c1-200)); set --variants-sample"
}

# ---- detect format and caller
if grep -q '^##fileformat=VCF' "$TMP/input"; then
  FORMAT=vcf
elif head -1 "$TMP/input" | grep -q $'^ID\tFilter\tPos\tRef\tVariant\tVariantLevel'; then
  FORMAT=mtdnaserver2
else
  die "$MT_VARIANTS is neither a VCF nor an mtDNA-Server 2 variants.annotated.txt"
fi

RCRS="$TMP/rcrs.fa"
make_rcrs() {  # contig name
  awk -v name="$1" '
    /^>/ { next }
    { seq = seq $0 }
    END { seq = substr(seq, 1, 16569); print ">" name; for (i = 1; i <= length(seq); i += 60) print substr(seq, i, 60) }
  ' "$MITRIX_HOME/refs/chrMC_rCRS.fa" > "$RCRS"
}

VCF="$TMP/in.vcf"
blocks=1
if [[ "$FORMAT" == vcf ]]; then
  # MitoHPC writes an empty VCF followed by the real one; keep the last header block
  blocks=$(grep -c '^##fileformat' "$TMP/input" || true)
  awk -v n="$blocks" '/^##fileformat/ { k++ } k == n' "$TMP/input" > "$TMP/last.vcf"
  (( blocks == 1 )) || log "warning: $MT_VARIANTS has $blocks header blocks (MitoHPC concatenation); using the last"
  hdr() { grep '^##' "$TMP/last.vcf"; }
  if [[ -z "$CALLER" ]]; then
    if   hdr | grep -q '^##mityCommandline';            then CALLER=mity
    elif hdr | grep -q '^##INFO=<ID=MLC_score,';        then CALLER=mitohpc
    elif hdr | grep -qi '^##source=Mutect2';            then CALLER=mutect2
    else CALLER=vcf; fi
  fi
  [[ "$CALLER" != mtdnaserver2 ]] || die "--caller mtdnaserver2 needs variants.annotated.txt, not a VCF"

  # ---- sample
  read -r -a NAMES <<< "$(grep '^#CHROM' "$TMP/last.vcf" | cut -f10- | tr '\t' ' ')"
  (( ${#NAMES[@]} >= 1 )) || die "VCF has no sample column"
  if ! VCF_SAMPLE=$(pick_sample "${NAMES[@]}"); then
    # MitoHPC: the ##sample line can carry the name when the column does not match
    HDR_SAMPLE=$(hdr | sed -n 's/^##sample=//p' | head -1)
    if (( ${#NAMES[@]} == 1 )) && [[ -n "$HDR_SAMPLE" ]] && pick_sample "$HDR_SAMPLE" >/dev/null; then
      VCF_SAMPLE=${NAMES[0]}
    else
      no_match "${NAMES[@]}"
    fi
  fi
  bcftools view --no-version -s "$VCF_SAMPLE" -o "$VCF" "$TMP/last.vcf" 2> "$TMP/view.log" \
    || die "bcftools view failed: $(head -1 "$TMP/view.log")"

  # ---- one mt contig, rCRS length
  VCF_CHR=$(grep -v '^#' "$VCF" | cut -f1 | sort -u)
  if [[ -z "$VCF_CHR" ]]; then
    VCF_CHR=$(grep '^##contig' "$VCF" | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p' \
      | grep -E '^(chrM|MT|chrMT|rCRS|NC_012920(\.1)?)$' | head -1)
    [[ -n "$VCF_CHR" ]] || VCF_CHR=chrM
  fi
  (( $(printf '%s\n' "$VCF_CHR" | wc -l) == 1 )) || die "VCF has records on several contigs: $(echo $VCF_CHR)"
  case "$VCF_CHR" in chrM|MT|chrMT|rCRS|NC_012920|NC_012920.1) ;; *) die "VCF contig '$VCF_CHR' is not an mtDNA contig";; esac
  VCF_LEN=$(grep '^##contig' "$VCF" | sed -n "s/^##contig=<ID=$VCF_CHR,length=\([0-9]*\).*/\1/p" | head -1)
  if [[ -z "$VCF_LEN" ]]; then
    log "warning: no ##contig length for $VCF_CHR; assuming rCRS"
  elif [[ "$VCF_LEN" != 16569 ]]; then
    die "VCF contig $VCF_CHR has length $VCF_LEN, not 16569 (rCRS)"
  fi
  if   grep -q '^##FORMAT=<ID=AF,'  "$VCF"; then HF_TAG=AF
  elif grep -q '^##FORMAT=<ID=VAF,' "$VCF"; then HF_TAG=VAF
  else die "VCF has no FORMAT/AF or FORMAT/VAF; cannot get heteroplasmy fraction"; fi
  if   grep -q '^##FORMAT=<ID=DP,' "$VCF"; then DP_FMT='[%DP]'
  elif grep -q '^##INFO=<ID=DP,'   "$VCF"; then DP_FMT='%INFO/DP'
  else DP_FMT='.'; log "warning: no DP in the VCF; depth and alt_reads will be '.'"; fi
else
  CALLER=${CALLER:-mtdnaserver2}
  [[ "$CALLER" == mtdnaserver2 ]] || die "a variants.annotated.txt input needs --caller mtdnaserver2"
  read -r -a NAMES <<< "$(tail -n +2 "$TMP/input" | cut -f1 | sort -u | tr '\n' ' ')"
  VCF_SAMPLE=$(pick_sample "${NAMES[@]}") || no_match "${NAMES[@]}"
  VCF_CHR=chrM
  # Convert to a VCF (VariantLevel is the ALT level; with several levels for one ALT,
  # as at 3107N, the first is used)
  {
    printf '##fileformat=VCFv4.2\n##contig=<ID=chrM,length=16569>\n'
    printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    printf '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Variant level">\n'
    printf '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Coverage">\n'
    tail -n +2 "$TMP/input" | awk -F'\t' -v s="$VCF_SAMPLE" '$1 == s { n = split($2, f, ";"); for (i = 1; i <= n; i++) if (f[i] != "PASS") print f[i] }' \
      | sort -u | awk '{ printf "##FILTER=<ID=%s,Description=\"mtDNA-Server 2 filter\">\n", $1 }'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$VCF_SAMPLE"
    tail -n +2 "$TMP/input" | awk -F'\t' -v OFS='\t' -v s="$VCF_SAMPLE" '
      $1 == s { split($6, lv, ","); print "chrM", $3, ".", $4, $5, ".", $2, ".", "GT:AF:DP", (lv[1] >= 0.9 ? "1/1" : "0/1") ":" lv[1] ":" $8 }' \
      | sort -k2,2n
  } > "$VCF"
  HF_TAG=AF; DP_FMT='[%DP]'
fi

# ---- REF vs rCRS, left-align, split multi-allelic records
make_rcrs "$VCF_CHR"
# ALT == REF rows (MitoHPC's reference-allele fraction) are not variants, and bcftools
# norm rejects them as duplicate alleles: drop them first
awk -F'\t' '/^#/ || $4 != $5' "$VCF" > "$TMP/clean.vcf"
n_refalt=$(( $(grep -vc '^#' "$VCF" || true) - $(grep -vc '^#' "$TMP/clean.vcf" || true) ))
NORM="$TMP/norm.vcf"
# records whose REF does not match rCRS are dropped with a warning (MitoHPC writes some
# poly-C indels, e.g. chrM:16180 AC>A, against the sample's own consensus sequence)
bcftools norm --no-version -f "$RCRS" -c wx -m -any -o "$NORM" "$TMP/clean.vcf" 2> "$TMP/norm.log" \
  || die "bcftools norm failed: $(grep -v '^Lines' "$TMP/norm.log" | head -1)"
n_refmis=$(grep -c '^REF_MISMATCH' "$TMP/norm.log" || true)
(( n_refmis == 0 )) || log "warning: dropped $n_refmis record(s) whose REF does not match rCRS: $(awk -F'\t' '$1 == "REF_MISMATCH" { print $3 " REF " $4 " (rCRS " $5 ")" }' "$TMP/norm.log" | head -5 | paste -sd';' -)"
n_norm=$(diff <(grep -v '^#' "$TMP/clean.vcf" | cut -f2,4,5 | sort) <(grep -v '^#' "$NORM" | cut -f2,4,5 | sort) | grep -c '^>' || true)
(( n_norm == 0 )) || log "normalised $n_norm record(s) (left-aligned / split)"

# ---- standard table
bcftools query -f "%POS\t%REF\t%ALT\t%FILTER\t[%$HF_TAG]\t$DP_FMT\n" "$NORM" \
  | awk -F'\t' -v OFS='\t' -v s="$SAMPLE" -v c="$CALLER" -v stats="$TMP/stats" '
      { n_in++ }
      $3 == "*" || $3 == "." { n_star++; next }
      $5 == "." || $5 == "" { n_noaf++; next }
      {
        t = (length($2) == length($3)) ? (length($2) == 1 ? "SNV" : "MNV") : "INDEL"
        dp = ($6 == "" ? "." : $6)
        ar = (dp == "." ? "." : int($5 * dp + 0.5))
        print s, $1, $2, $3, $5 + 0, dp, ar, $4, c, t
        n_out++; if ($4 == "PASS") n_pass++
      }
      END { printf "%d %d %d %d %d\n", n_in, n_star, n_noaf, n_out, n_pass > stats }' \
  > "$TMP/body.tsv"
read -r n_in n_star n_noaf n_out n_pass < "$TMP/stats"
n_in=$(( n_in + n_refalt ))
(( n_noaf == 0 )) || log "warning: dropped $n_noaf record(s) with no HF"
log "$CALLER: $n_out variants ($n_pass PASS) for $VCF_SAMPLE; dropped $n_refalt ALT==REF and $n_star '*' row(s)"
if (( n_out > 0 && n_pass == n_out )); then
  log "warning: every record is PASS — the input looks pre-filtered; heteroplasmy counts will only see PASS calls"
fi
{ printf 'sample\tpos\tref\talt\thf\tdepth\talt_reads\tfilter\tcaller\ttype\n'; cat "$TMP/body.tsv"; } \
  > "$OUT/${SAMPLE}.variants.tsv"

# ---- clean VCF for HaploGrep3 / haplocheck: each caller's artefact filters removed
case "$CALLER" in mity) EXCLUDE=FAIL;; *) EXCLUDE=strand_bias;; esac
{
  printf '##fileformat=VCFv4.2\n##contig=<ID=chrM,length=16569>\n##source=mitRix-%s\n' "$MITRIX_VERSION"
  printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype: 1/1 if AF >= 0.9, else 0/1">\n'
  printf '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Heteroplasmy fraction">\n'
  printf '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$SAMPLE"
  awk -F'\t' -v OFS='\t' -v ex="$EXCLUDE" '
    BEGIN { n = split(ex, e, ","); for (i = 1; i <= n; i++) bad[e[i]] = 1 }
    { k = split($8, f, ";"); for (i = 1; i <= k; i++) if (f[i] in bad) next
      print "chrM", $2, ".", $3, $4, ".", "PASS", ".", "GT:AF:DP", ($5 >= 0.9 ? "1/1" : "0/1") ":" $5 ":" ($6 == "." ? "." : $6) }' \
    "$TMP/body.tsv"
} | bcftools view --no-version -Oz -o "$OUT/${SAMPLE}.mt.vcf.gz"
bcftools index -f -t "$OUT/${SAMPLE}.mt.vcf.gz"
n_clean=$(bcftools view -H "$OUT/${SAMPLE}.mt.vcf.gz" | wc -l)

{
  write_param mt_variants        "$(readlink -f "$MT_VARIANTS")"
  write_param input_format       "$FORMAT"
  write_param caller             "$CALLER"
  write_param variants_sample    "$VCF_SAMPLE"
  write_param input_contig       "$VCF_CHR"
  write_param hf_field           "$HF_TAG"
  write_param header_blocks      "$blocks"
  write_param records_in         "$n_in"
  write_param dropped_ref_eq_alt "$n_refalt"
  write_param dropped_star       "$n_star"
  write_param dropped_no_hf      "$n_noaf"
  write_param normalised         "$n_norm"
  write_param ref_mismatch_dropped "$n_refmis"
  write_param variants           "$n_out"
  write_param variants_pass      "$n_pass"
  write_param exclude_filters    "$EXCLUDE"
  write_param clean_vcf_records  "$n_clean"
} > "$OUT/${SAMPLE}.variants.params.tsv"

# ---- variants-only mode (no BAM): chrM depth from the calls' DP, as a constant
# "per-base" coverage file so heteroplasmy and pathogenic can use it
if [[ -z "$BAM" ]]; then
  dp=$(awk -F'\t' 'NR > 1 && $6 ~ /^[0-9.]+$/ { print $6 }' "$OUT/${SAMPLE}.variants.tsv" | sort -n \
       | awk '{ a[NR] = $1 } END { if (NR) printf "%d", (NR % 2 ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2) }')
  [[ -n "$dp" ]] || { dp=0; log "warning: no DP in the variant calls; depth set to 0 (low-HF metrics flagged low_depth)"; }
  printf 'chrM\t0\t16569\t%s\n' "$dp" | gzip -c > "$OUT/${SAMPLE}.chrM.per-base.bed.gz"
  { write_param depth_source "median DP of variant calls (variants-only mode, no BAM)"
    write_param depth_median "$dp"; } >> "$OUT/${SAMPLE}.variants.params.tsv"
  log "variants-only mode: chrM depth $dp (median DP of the variant calls)"
fi
