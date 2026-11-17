#!/usr/bin/env bash
# aggregate/carriers.sh — carriers of user-listed nuclear variants (mitrix carriers):
# genotypes of each listed variant in a multi-sample nuclear VCF, summarised per gene.

set -euo pipefail
LOG_PREFIX=carriers
HERE=$(dirname "$(readlink -f "$0")")
source "$HERE/../lib/common.sh"

usage() {
  cat >&2 <<EOF
mitrix carriers — carriers of listed nuclear variants (e.g. disease genes)

USAGE
  mitrix carriers --vcf NUCLEAR.vcf.gz --list NAME=FILE [--list NAME2=FILE2] --outdir DIR

OPTIONS
  --vcf FILE           multi-sample nuclear VCF/BCF (bgzipped; indexed is faster).
                       Sample IDs = mitRix sample names. Positions must be on the same
                       build as the lists.
  --list NAME=FILE     variant list, repeatable. TSV with a header; required columns
                       chr, pos, ref, alt, gene; optional label, inheritance (AD, AR,
                       XL, risk, ...), classification, severity, rsid, source. '#'
                       lines are comments. NAME: letters, digits, '_'.
                       Example: presets/examples/pd_nuclear_variants_GRCh38.tsv (PD, ClinVar)
  --reference-fasta FA left-align VCF records and list entries against FA before
                       matching (recommended for indels unless the VCF is already
                       normalised with bcftools norm -f)
  --outdir DIR         writes DIR/carriers/
  --threads N          default 4

Records are split into biallelic records and matched on chr, pos, ref, alt ('chr'
prefix differences are handled). Dosage = copies of the listed ALT allele (0, 1, 2);
missing genotypes are NA, not 0. A listed variant absent from the VCF counts as not
carried, which is right for a joint-called cohort VCF (it lists every variant seen in
any sample) but not for a VCF filtered to some samples or regions; see match.tsv.
Structural variants (e.g. PRKN exon deletions) are not detected.

OUTPUT (DIR/carriers/)
  carriers.tsv   long format (mitrix cohort --carriers), per sample and list:
                   NAME_n        listed variants carried (count); entries whose
                                 classification contains "modifier" (e.g. APOE e4)
                                 count only in their gene columns, not here or in NAME
                   NAME_<GENE>   listed ALT alleles carried in that gene (count; 2 in an
                                 AR gene = biallelic or two variants on one haplotype,
                                 the VCF is unphased)
                   NAME_<GENE>_<SEVERITY>
                                 the same per severity class, for genes whose list
                                 entries have more than one severity value (e.g.
                                 PD_GBA1_severe, PD_GBA1_mild, PD_GBA1_risk)
                   NAME          carried variants, e.g. LRRK2:p.Gly2019Ser(het) or none
                   NAME_missing  listed variants present in the VCF with a missing
                                 genotype for this sample (count)
                 qc_flag: missing_genotypes when NAME_missing > 0
  variants.tsv   one row per carried variant: sample, list, gene, label, chr, pos,
                 ref, alt, dosage, inheritance, classification
  match.tsv      per list: variants listed, found in the VCF, carriers
  logs/
EOF
}

VCF=""; OUTDIR=""; LISTS=(); FASTA=""
while (( $# )); do
  case "$1" in
    --vcf)             VCF="${2:?}"; shift 2;;
    --list)            LISTS+=("${2:?}"); shift 2;;
    --reference-fasta) FASTA="${2:?}"; shift 2;;
    --outdir)          OUTDIR="${2:?}"; shift 2;;
    --threads)         THREADS="${2:?}"; shift 2;;
    -h|--help)         usage; exit 0;;
    *)                 die "unknown option: $1";;
  esac
done
[[ -s "$VCF" ]] || die "--vcf not found: $VCF"
[[ -n "$OUTDIR" ]] || die "--outdir is required"
(( ${#LISTS[@]} )) || die "at least one --list NAME=FILE is required"
[[ -z "$FASTA" || -s "$FASTA" ]] || die "--reference-fasta not found: $FASTA"
THREADS="${THREADS:-4}"
need_tools bcftools Rscript
O="$OUTDIR/carriers"; mkdir -p "$O/logs"
make_tmpdir

# VCF contig naming: 'chr1' or '1'
first=$(bcftools view -h "$VCF" | awk -F'[=,<>]' '/^##contig=<ID=/ && !f++ { print $4 }')
[[ -n "$first" ]] || first=$( (bcftools view -H "$VCF" 2> /dev/null | head -1 | cut -f1) || true )
PREFIX=""; [[ "$first" == chr* ]] && PREFIX=chr
INDEXED=0; [[ -e "$VCF.tbi" || -e "$VCF.csi" ]] && INDEXED=1

NAMES=(); SPECS=()
for l in "${LISTS[@]}"; do
  [[ "$l" == *=* ]] || die "--list must be NAME=FILE: $l"
  n=${l%%=*}; f=${l#*=}
  [[ "$n" =~ ^[A-Za-z][A-Za-z0-9_]*$ ]] || die "--list name must be letters, digits, '_': $n"
  [[ " ${NAMES[*]} " != *" $n "* ]] || die "--list name used twice: $n"
  [[ -s "$f" ]] || die "--list file not found: $f"
  # list -> chr pos ref alt + metadata, contigs renamed to the VCF's style
  awk -F'\t' -v OFS='\t' -v pre="$PREFIX" '
    /^#/ || !NF { next }
    !h { for (j = 1; j <= NF; j++) c[$j] = j; h = 1
         split("chr pos ref alt gene", r, " ")
         for (i = 1; i <= 5; i++) if (!(r[i] in c)) { print "missing column " r[i] > "/dev/stderr"; exit 3 }
         print "id", "chr", "pos", "ref", "alt", "gene", "label", "inheritance", "classification", "severity"; next }
    { ch = $c["chr"]; sub(/^chr/, "", ch); ch = pre ch
      lab = ("label" in c && $c["label"] != "") ? $c["label"] : $c["chr"] ":" $c["pos"] $c["ref"] ">" $c["alt"]
      print ++id, ch, $c["pos"], toupper($c["ref"]), toupper($c["alt"]), $c["gene"], lab,
            ("inheritance" in c) ? $c["inheritance"] : "", ("classification" in c) ? $c["classification"] : "",
            ("severity" in c) ? $c["severity"] : "" }' \
    "$f" > "$TMP/$n.list" || die "--list $n: $f needs a header with chr, pos, ref, alt, gene"
  nl=$(( $(wc -l < "$TMP/$n.list") - 1 )); (( nl > 0 )) || die "--list $n: no variants in $f"
  if [[ -n "$FASTA" ]]; then   # left-align list entries like the VCF records
    { printf '##fileformat=VCFv4.2\n'; awk -F'\t' 'NR > 1 && !s[$2]++ { print "##contig=<ID=" $2 ">" }' "$TMP/$n.list"
      printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
      awk -F'\t' -v OFS='\t' 'NR > 1 { print $2, $3, $1, $4, $5, ".", ".", "." }' "$TMP/$n.list"; } \
      | bcftools norm --no-version -f "$FASTA" -c wx 2> "$O/logs/$n.list_norm.log" \
      | bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\n' > "$TMP/$n.norm" || die "--list $n: bcftools norm failed"
    awk -F'\t' -v OFS='\t' 'NR == FNR { p[$1] = $2 FS $3 FS $4 FS $5; next }
      FNR == 1 { print; next } ($1 in p) { split(p[$1], x, FS); $2 = x[1]; $3 = x[2]; $4 = x[3]; $5 = x[4]; print }' \
      "$TMP/$n.norm" "$TMP/$n.list" > "$TMP/$n.list2" && mv "$TMP/$n.list2" "$TMP/$n.list"
  fi
  awk -F'\t' -v OFS='\t' 'NR > 1 { print $2, $3 }' "$TMP/$n.list" >> "$TMP/regions.tsv"
  log "list $n: $nl variants ($f)"
  NAMES+=("$n"); SPECS+=("$n=$TMP/$n.list")
done
sort -k1,1 -k2,2n -u "$TMP/regions.tsv" > "$TMP/regions.sorted.tsv"

# listed positions only; split multi-allelic records; optionally left-align
if (( INDEXED )); then SEL=(-R "$TMP/regions.sorted.tsv"); else SEL=(-T "$TMP/regions.sorted.tsv"); log "VCF not indexed: streaming the whole file"; fi
NORM=(-m -any); [[ -n "$FASTA" ]] && NORM+=(-f "$FASTA" -c wx)
bcftools view --threads "$THREADS" "${SEL[@]}" "$VCF" -Ou 2> "$O/logs/view.log" \
  | bcftools norm --no-version "${NORM[@]}" -Ou 2> "$O/logs/norm.log" \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > "$TMP/gt.tsv" || die "bcftools failed (see $O/logs/)"
bcftools query -l "$VCF" > "$TMP/samples.txt"
log "$(wc -l < "$TMP/samples.txt") samples, $(wc -l < "$TMP/gt.tsv") VCF records at listed positions"

Rscript --vanilla "$HERE/carriers.R" "$TMP/gt.tsv" "$TMP/samples.txt" "$MITRIX_VERSION" "$O" "${SPECS[@]}" \
  2> "$O/logs/carriers.log" || { cat "$O/logs/carriers.log" >&2; die "carriers.R failed"; }
column -t "$O/match.tsv" | sed 's/^/carriers: /' >&2
log "wrote $O/carriers.tsv"
