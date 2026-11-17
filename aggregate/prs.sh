#!/usr/bin/env bash
# aggregate/prs.sh — polygenic scores from a multi-sample nuclear VCF (mitrix prs).
#
# Needs plink2 and Rscript. Uses the weights, SNP map and ancestry reference shipped in
# presets/prs/ (built by tools/make_prs_weights.sh: PRS-CS-auto on Gupta et al. GWAS).
# Steps: genotype QC -> ancestry (projection onto 1000 Genomes PCs) -> HWE in EUR ->
# match variants to the weights -> plink2 --score -> standardise within EUR.

set -euo pipefail
LOG_PREFIX=prs
HERE=$(dirname "$(readlink -f "$0")")
source "$HERE/../lib/common.sh"
P="${MITRIX_PRS_DIR:-$MITRIX_HOME/presets/prs}"   # override for tests

usage() {
  cat >&2 <<EOF
usage: mitrix prs --vcf FILE --outdir DIR [options]

Scores every sample of a multi-sample nuclear VCF (GRCh38 by default). Sample IDs must
be the mitRix sample names. Writes DIR/prs/: prs.tsv (long format, for mitrix cohort
--prs), match.tsv, qc.tsv, ancestry.tsv, params.tsv, logs/.

  --build 38|37        genome build of the VCF (default 38)
  --scores LIST        default: all in presets/prs/scores.tsv (CN,CN_raw,het,7S,primer)
  --geno X             max variant missingness (default 0.05)
  --mind X             max sample missingness (default 0.05)
  --maf X              min minor-allele frequency (default 0.01)
  --hwe P              min HWE p-value, tested in EUR samples (default 1e-6)
  --min-overlap X      flag scores with fewer weight SNPs found (default 0.8)
  --cohort FILE        cohort.tsv: also writes validation.tsv (each PRS vs its readout)
  --metadata FILE      metadata (sample column), for validation covariates
  --covariates LIST    metadata columns adjusted for in validation.tsv, e.g.
                       age,sex,neutrophil_frac (default: age,sex if present)
Lead-variant scores (PRS_<score>_lead: GWAS betas of independent genome-wide
significant SNPs) are always produced alongside the PRS-CS scores.
  --threads N          default 4
EOF
}

VCF=""; OUTDIR=""; BUILDV=38; SCORES=""; LEAD=1; COVARS=age,sex; GENO=0.05; MIND=0.05; MAF=0.01; HWE=1e-6
MINOV=0.8; COHORT=""; META=""
while (( $# )); do
  case "$1" in
    --vcf)         VCF="${2:?}"; shift 2;;
    --outdir)      OUTDIR="${2:?}"; shift 2;;
    --build)       BUILDV="${2:?}"; shift 2;;
    --scores)      SCORES="${2:?}"; shift 2;;
    --lead)        shift;;   # always on (kept for compatibility)
    --covariates)  COVARS="${2:?}"; shift 2;;
    --geno)        GENO="${2:?}"; shift 2;;
    --mind)        MIND="${2:?}"; shift 2;;
    --maf)         MAF="${2:?}"; shift 2;;
    --hwe)         HWE="${2:?}"; shift 2;;
    --min-overlap) MINOV="${2:?}"; shift 2;;
    --cohort)      COHORT="${2:?}"; shift 2;;
    --metadata)    META="${2:?}"; shift 2;;
    --threads)     THREADS="${2:?}"; shift 2;;
    -h|--help)     usage; exit 0;;
    *)             die "unknown option: $1";;
  esac
done
[[ -s "$VCF" ]] || die "--vcf not found: $VCF"
[[ -n "$OUTDIR" ]] || die "--outdir is required"
[[ "$BUILDV" == 38 || "$BUILDV" == 37 ]] || die "--build must be 38 or 37"
[[ -z "$COHORT" || -s "$COHORT" ]] || die "--cohort not found: $COHORT"
[[ -z "$META" || -s "$META" ]] || die "--metadata not found: $META"
[[ -s "$P/snps.tsv.gz" && -s "$P/scores.tsv" ]] || die "PRS weights missing in $P (tools/make_prs_weights.sh)"
need_tools plink2 Rscript awk
[[ -n "$SCORES" ]] || SCORES=$(awk -F'\t' 'NR > 1 { sub(/^PRS_/, "", $1); print $1 }' "$P/scores.tsv" | paste -sd,)
IFS=, read -r -a SC <<< "$SCORES"
for s in "${SC[@]}"; do [[ -s "$P/PRS_$s.weights.tsv.gz" ]] || die "no weights for score '$s' in $P"; done

O="$OUTDIR/prs"; mkdir -p "$O/logs"
make_tmpdir
G="$TMP/g"
pl() {  # step plink2-args... (log to O/logs/step.log)
  local step=$1; shift
  plink2 --threads "$THREADS" "$@" > "$O/logs/$step.log" 2>&1 || { tail -5 "$O/logs/$step.log" >&2; die "plink2 $step failed"; }
}
count() { wc -l < "$1" | tr -d ' '; }

# ---- 1. import: autosomal biallelic SNVs, IDs chr:pos:ref:alt, call rates
log "importing $VCF (GRCh$BUILDV)"
pl import --vcf "$VCF" --allow-extra-chr --vcf-half-call m --autosome --snps-only just-acgt --max-alleles 2 \
  --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 1 missing --rm-dup exclude-all \
  --output-chr 26 --make-pgen --out "$G.raw"
n_samp_in=$(( $(count "$G.raw.psam") - 1 )); n_var_in=$(grep -vc '^#' "$G.raw.pvar")
pl callrate --pfile "$G.raw" --geno "$GENO" --mind "$MIND" --make-pgen --out "$G.cr"
n_samp=$(( $(count "$G.cr.psam") - 1 )); n_var_cr=$(grep -vc '^#' "$G.cr.pvar")
(( n_samp > 0 )) || die "no samples left after --mind $MIND"
log "$n_samp/$n_samp_in samples, $n_var_cr/$n_var_in variants after call-rate QC"

# target variants keyed by chr:pos -> "ID REF ALT"
awk '!/^#/ { print $1 ":" $2 "\t" $3 "\t" $4 "\t" $5 }' "$G.cr.pvar" > "$TMP/target.tsv"
POSCOL=$([[ "$BUILDV" == 38 ]] && echo 4 || echo 3)

# ---- 2. ancestry: project onto 1000 Genomes PCs
# loadings rows: rsid ref alt alt_freq allele PC1..10 ; map rsid -> chr:pos in this build
zcat "$P/snps.tsv.gz" | awk -v pc="$POSCOL" 'NR > 1 { print $1 "\t" $2 ":" $pc }' > "$TMP/map.pos"
awk -F'\t' -v OFS='\t' -v out="$TMP" '
  FILENAME == ARGV[1] { pos[$1] = $2; next }
  FILENAME == ARGV[2] { id[$1] = $2; ref[$1] = $3; alt[$1] = $4; next }
  FNR == 1 { print "ID", "A1", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10" > out "/proj.score"
             print "#ID", "REF", "ALT", "ALT_FREQS" > out "/proj.afreq"; next }
  !($1 in pos) || !(pos[$1] in id) { next }
  { k = pos[$1]
    if ($2 == ref[k] && $3 == alt[k]) f = $4
    else if ($2 == alt[k] && $3 == ref[k]) f = 1 - $4
    else next
    if (!(k in done)) { print id[k], ref[k], alt[k], f > out "/proj.afreq"; done[k] = 1 }
    printf "%s\t%s", id[k], $5 > out "/proj.score"
    for (i = 6; i <= 15; i++) printf "\t%s", $i > out "/proj.score"
    print "" > out "/proj.score" }' \
  "$TMP/map.pos" "$TMP/target.tsv" <(zcat "$P/ancestry_loadings.tsv.gz")
n_anc_ref=$(( $(zcat "$P/ancestry_loadings.tsv.gz" | awk 'NR > 1 { print $1 }' | sort -u | wc -l) ))
n_anc=$(( $(count "$TMP/proj.afreq") - 1 ))
log "ancestry: $n_anc/$n_anc_ref reference SNPs found"
(( n_anc >= 1000 )) || log "warning: few ancestry SNPs found; ancestry assignment will be unreliable"
pl ancestry --pfile "$G.cr" --read-freq "$TMP/proj.afreq" --score "$TMP/proj.score" 1 2 header-read \
  no-mean-imputation variance-standardize --score-col-nums 3-12 --out "$TMP/proj"
Rscript --vanilla "$HERE/prs_ancestry.R" "$TMP/proj.sscore" "$P/ancestry_reference.tsv.gz" "$O/ancestry.tsv" \
  2> "$O/logs/ancestry_assign.log" || { cat "$O/logs/ancestry_assign.log" >&2; die "ancestry assignment failed"; }
awk -F'\t' 'NR > 1 && $2 == "EUR" { print $1 }' "$O/ancestry.tsv" > "$TMP/eur.ids"
n_eur=$(count "$TMP/eur.ids")
log "$n_eur/$n_samp samples assigned EUR"

# ---- 3. HWE (EUR samples) and MAF
HWE_ARGS=()
if (( n_eur >= 50 )); then
  pl hardy --pfile "$G.cr" --keep "$TMP/eur.ids" --hwe "$HWE" --write-snplist --out "$TMP/hwe"
  HWE_ARGS=(--extract "$TMP/hwe.snplist")
else
  log "warning: $n_eur EUR samples (< 50); HWE filter skipped"
fi
pl final --pfile "$G.cr" "${HWE_ARGS[@]}" --maf "$MAF" --make-pgen --out "$G"
n_var=$(grep -vc '^#' "$G.pvar")
log "$n_var variants after HWE and MAF filters"

# ---- 4. match weights; one plink2 score file with a column per score
awk '!/^#/ { print $1 ":" $2 "\t" $3 "\t" $4 "\t" $5 }' "$G.pvar" > "$TMP/final.tsv"
zcat "$P/snps.tsv.gz" > "$TMP/map.tsv"
wfiles=()
for s in "${SC[@]}"; do zcat "$P/PRS_$s.weights.tsv.gz" > "$TMP/PRS_$s.w"; wfiles+=("$TMP/PRS_$s.w"); done
paste "$TMP/map.tsv" "${wfiles[@]}" > "$TMP/mapw.tsv"
awk -F'\t' -v OFS='\t' -v pc="$POSCOL" -v ns="${#SC[@]}" -v names="${SC[*]}" -v out="$TMP" '
  BEGIN { comp["A"] = "T"; comp["T"] = "A"; comp["C"] = "G"; comp["G"] = "C"; split(names, nm, " ") }
  FILENAME == ARGV[1] { id[$1] = $2; ref[$1] = $3; alt[$1] = $4; next }
  FNR == 1 { printf "ID\tA1" > out "/weights.score"; for (j = 1; j <= ns; j++) printf "\tPRS_%s", nm[j] > out "/weights.score"
             print "" > out "/weights.score"; next }
  {
    for (j = 1; j <= ns; j++) if ($(6 + j) != 0) { nz[j]++ }
    k = $2 ":" $pc
    if (!(k in id)) next
    e = $5; o = $6
    if ((e == ref[k] && o == alt[k]) || (e == alt[k] && o == ref[k])) a = e
    else if ((comp[e] == ref[k] && comp[o] == alt[k]) || (comp[e] == alt[k] && comp[o] == ref[k])) { a = comp[e]; flip++ }
    else { mism++; next }
    printf "%s\t%s", id[k], a > out "/weights.score"
    for (j = 1; j <= ns; j++) { printf "\t%s", $(6 + j) > out "/weights.score"; if ($(6 + j) != 0) found[j]++ }
    print "" > out "/weights.score"
  }
  END {
    print "score\tweight_snps\tfound\tfraction\tstrand_flipped\tallele_mismatch" > out "/match.tsv"
    for (j = 1; j <= ns; j++) printf "PRS_%s\t%d\t%d\t%.4f\t%d\t%d\n", nm[j], nz[j], found[j], (nz[j] ? found[j] / nz[j] : 0), flip, mism > out "/match.tsv"
  }' "$TMP/final.tsv" "$TMP/mapw.tsv"
cp "$TMP/match.tsv" "$O/match.tsv"
column -t "$O/match.tsv" | sed 's/^/prs: /' >&2

# ---- 5. score (sums; missing genotypes mean-imputed from cohort frequencies, passed
# explicitly so plink2 also scores cohorts of fewer than 50 samples)
pl freq --pfile "$G" --freq --out "$TMP/cohort"
(( n_samp >= 50 )) || log "warning: $n_samp samples: missing genotypes imputed from a small cohort's allele frequencies"
pl score --pfile "$G" --read-freq "$TMP/cohort.afreq" --score "$TMP/weights.score" 1 2 header-read cols=+scoresums \
  --score-col-nums 3-$(( 2 + ${#SC[@]} )) --out "$TMP/prs"
LEADF=""
if (( LEAD )); then
  # lead-variant scores: rsid effect_allele weight, matched like the PRS-CS weights
  for s in "${SC[@]}"; do
    [[ -s "$P/PRS_$s.lead.tsv.gz" ]] || { log "PRS_${s}_lead: no lead-variant file"; continue; }
    awk -F'\t' -v OFS='\t' -v pc="$POSCOL" -v s="$s" '
      BEGIN { comp["A"] = "T"; comp["T"] = "A"; comp["C"] = "G"; comp["G"] = "C" }
      FILENAME == ARGV[1] { id[$1] = $2; ref[$1] = $3; alt[$1] = $4; next }
      FILENAME == ARGV[2] { if (FNR > 1) { k[$1] = $2 ":" $pc; e[$1] = $5; o[$1] = $6 }; next }
      FNR == 1 { next }
      ($1 in k) && (k[$1] in id) {
        kk = k[$1]; a = $2; other = (a == e[$1]) ? o[$1] : e[$1]
        if ((a == ref[kk] && other == alt[kk]) || (a == alt[kk] && other == ref[kk])) print id[kk], a, $3
        else if ((comp[a] == ref[kk] && comp[other] == alt[kk]) || (comp[a] == alt[kk] && comp[other] == ref[kk])) print id[kk], comp[a], $3 }' \
      "$TMP/final.tsv" "$TMP/map.tsv" <(zcat "$P/PRS_$s.lead.tsv.gz") > "$TMP/lead_$s.score"
    if [[ -s "$TMP/lead_$s.score" ]]; then
      pl "lead_$s" --pfile "$G" --read-freq "$TMP/cohort.afreq" --score "$TMP/lead_$s.score" 1 2 3 cols=+scoresums --out "$TMP/lead_$s"
      LEADF="$LEADF $s:$TMP/lead_$s.sscore"
    else
      log "PRS_${s}_lead: no lead variants found in the VCF"
    fi
  done
fi

# ---- 6. long table, validation
Rscript --vanilla "$HERE/prs_table.R" "$TMP/prs.sscore" "$O/ancestry.tsv" "$O/match.tsv" "$MINOV" \
  "$MITRIX_VERSION" "$O/prs.tsv" "${COHORT:-}" "${META:-}" "$O/validation.tsv" "${LEADF# }" "$COVARS" \
  2> "$O/logs/table.log" || { cat "$O/logs/table.log" >&2; die "prs table failed"; }

{
  printf 'step\tsamples\tvariants\n'
  printf 'input\t%s\t%s\n' "$n_samp_in" "$n_var_in"
  printf 'call_rate\t%s\t%s\n' "$n_samp" "$n_var_cr"
  printf 'hwe_maf\t%s\t%s\n' "$n_samp" "$n_var"
} > "$O/qc.tsv"
{
  write_param vcf "$(readlink -f "$VCF")"; write_param build "GRCh$BUILDV"
  write_param scores "$SCORES"; write_param validation_covariates "$COVARS"
  write_param geno "$GENO"; write_param mind "$MIND"; write_param maf "$MAF"; write_param hwe "$HWE"
  write_param hwe_samples "$n_eur EUR"; write_param min_overlap "$MINOV"
  write_param ancestry_snps "$n_anc/$n_anc_ref"; write_param weights "$P"
  write_param plink2 "$(plink2 --version | head -1)"; write_param mitrix "$MITRIX_VERSION"
} > "$O/params.tsv"
log "wrote $O/prs.tsv"
