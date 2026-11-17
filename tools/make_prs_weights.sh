#!/usr/bin/env bash
# tools/make_prs_weights.sh — stage A of mitRix PRS (maintainer, run rarely): build the
# PRS weights, SNP map, heritability table and ancestry reference shipped in presets/prs/.
# Users never run this; `mitrix prs` only needs plink2.
#
#   usage: tools/make_prs_weights.sh BUILD_DIR STEP
#
# BUILD_DIR holds the downloads (see BUILD_DIR/download.sh):
#   gwas/<GCST>.h.tsv.gz (+ .meta.yaml)   GWAS Catalog harmonised summary statistics
#   ref/ldblk_1kg_eur/                    PRS-CS 1000G EUR LD panel (HapMap3)
#   src/PRScs, src/ldsc                   PRS-CS and ldsc (git clones)
#   ldsc_ref/LDscore, ldsc_ref/1000G_Phase3_weights_hm3_no_MHC   1000G EUR LD scores, weights
#   kg/all_hg38.{pgen.zst,pvar.zst,psam}  1000 Genomes high-coverage GRCh38 (plink2 resources)
# Python envs: PRSCS_PY (numpy 1.26, scipy 1.11, h5py), LDSC_PY (Python 2.7 for ldsc).
#
# STEP (in order; prscs submits a Slurm array and returns):
#   prep      harmonise summary statistics -> work/<score>.{prscs,ldsc,pos38}.txt
#   prscs     PRS-CS-auto, one Slurm task per score x chromosome -> work/prscs/
#   merge     shared SNP map + weight files -> presets/prs/
#   ldsc      SNP heritability per score -> work/ldsc/
#   ancestry  1000 Genomes PCA reference (loadings, reference PCs) -> presets/prs/
#   lead      lead-variant scores (p < 5e-8, clumped in 1000G EUR) -> presets/prs/
#   table     presets/prs/scores.tsv (sources, N, SNPs, phi, h2, expected accuracy)

set -euo pipefail
LOG_PREFIX=make_prs_weights
source "$(dirname "$(readlink -f "$0")")/../lib/common.sh"

B=$(readlink -f "${1:?usage: make_prs_weights.sh BUILD_DIR STEP}")
STEP=${2:?usage: make_prs_weights.sh BUILD_DIR STEP}
PRSCS_PY=${PRSCS_PY:-$HOME/.conda/envs/prscs/bin/python}
LDSC_PY=${LDSC_PY:-$HOME/.conda/envs/ldsc/bin/python}
OUT="$MITRIX_HOME/presets/prs"
W="$B/work"
LD="$B/ref/ldblk_1kg_eur"
SNPINFO="$LD/snpinfo_1kg_hm3"
mkdir -p "$W" "$OUT"

# score  GWAS Catalog ID  N (sum of sample sizes in the GWAS Catalog metadata)  description
SCORES="CN:GCST90840190:401245:mtDNA copy number, adjusted (Gupta et al. 2026)
CN_raw:GCST90840191:401245:mtDNA copy number, raw (Gupta et al. 2026)
het:GCST90840194:519770:age-related mtSNV accumulation (Gupta et al. 2026)
7S:GCST90268499:130658:coverage discrepancy 7S DNA vs DNA flap (Gupta et al. 2023)
primer:GCST90268500:130658:coverage discrepancy RNA primer vs DNA flap (Gupta et al. 2023)"
# ONLY=CN_raw,... limits prep, prscs and ldsc to some scores (merge, lead, table use all)
each_score() {
  if [[ -n "${ONLY:-}" && "${1:-}" != all ]]; then
    printf '%s\n' "$SCORES" | awk -F: -v only=",$ONLY," 'index(only, "," $1 ",")'
  else printf '%s\n' "$SCORES"; fi
}

case "$STEP" in
prep)
  each_score | while IFS=: read -r s g n d; do
    log "prep $s ($g)"
    "$PRSCS_PY" "$MITRIX_HOME/tools/prs/prep_sumstats.py" "$B/gwas/$g.h.tsv.gz" "$SNPINFO" "$W/$s"
  done
  # PRS-CS restricts to SNPs in a target .bim: use every LD-panel SNP so the weights
  # work for any cohort
  awk 'NR > 1 { print $1 "\t" $2 "\t0\t" $3 "\t" $4 "\t" $5 }' "$SNPINFO" > "$W/hm3.bim"
  ;;
prscs)
  mkdir -p "$W/prscs" "$W/logs"
  each_score | cut -d: -f1,3 > "$W/scores.list"
  cat > "$W/prscs_task.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
i=\$SLURM_ARRAY_TASK_ID
line=\$(sed -n "\$(( i / 22 + 1 ))p" "$W/scores.list"); chr=\$(( i % 22 + 1 ))
s=\${line%%:*}; n=\${line#*:}
"$PRSCS_PY" "$B/src/PRScs/PRScs.py" --ref_dir="$LD" --bim_prefix="$W/hm3" \\
  --sst_file="$W/\$s.prscs.txt" --n_gwas="\$n" --out_dir="$W/prscs/\$s" --chrom="\$chr" --seed=42 \\
  | tee "$W/prscs/\$s.chr\$chr.log"
EOF
  chmod +x "$W/prscs_task.sh"
  n=$(( $(wc -l < "$W/scores.list") * 22 - 1 ))
  sbatch --array=0-"$n" -c 1 --mem 8G -t 8:00:00 -J prscs -o "$W/logs/prscs.%a.log" "$W/prscs_task.sh"
  ;;
merge)
  args=()
  while IFS=: read -r s g n d; do
    args+=("PRS_$s=$W/prscs/${s}_pst_eff_a1_b0.5_phiauto_chr*.txt:$W/$s.pos38.txt")
  done < <(each_score all)
  "$PRSCS_PY" "$MITRIX_HOME/tools/prs/merge_weights.py" "$SNPINFO" "$OUT" "${args[@]}"
  ;;
ldsc)
  mkdir -p "$W/ldsc"
  each_score | while IFS=: read -r s g n d; do
    log "ldsc $s"
    # ldsc .sumstats format (SNP A1 A2 Z N) written directly: the input is already
    # restricted to HapMap3 SNPs with matching alleles (munge_sumstats.py is very slow here)
    awk -v n="$n" 'NR == 1 { print "SNP\tA1\tA2\tZ\tN"; next } { printf "%s\t%s\t%s\t%.6g\t%d\n", $1, $2, $3, $4 / $5, n }' \
      "$W/$s.ldsc.txt" | gzip > "$W/ldsc/$s.sumstats.gz"
    "$LDSC_PY" "$B/src/ldsc/ldsc.py" --h2 "$W/ldsc/$s.sumstats.gz" \
      --ref-ld-chr "$B/ldsc_ref/LDscore/LDscore." --w-ld-chr "$B/ldsc_ref/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC." \
      --out "$W/ldsc/$s.h2" > /dev/null
    grep -E 'Total Observed scale h2' "$W/ldsc/$s.h2.log"
  done
  ;;
ancestry)
  need_tools plink2
  K="$W/kg"; mkdir -p "$K"
  # founders, autosomal biallelic SNVs among the SNP map, MAF >= 5%, outside long-range
  # LD regions (MHC, chr8 inversion), LD-pruned
  awk 'NR > 1 { print $2 }' "$SNPINFO" > "$K/map.ids"   # HapMap3 SNPs (the SNP map is a subset)
  printf '6\t25000000\t35000000\tMHC\n8\t8000000\t12000000\tinv8\n' > "$K/longld.txt"
  plink2 --pfile "$B/kg/all_hg38" vzs --allow-extra-chr --keep-founders --autosome --snps-only just-acgt \
    --max-alleles 2 --rm-dup exclude-all --extract "$K/map.ids" --exclude bed1 "$K/longld.txt" --maf 0.05 \
    --indep-pairwise 1000kb 1 0.1 --out "$K/prune" --threads 8 > /dev/null
  plink2 --pfile "$B/kg/all_hg38" vzs --allow-extra-chr --keep-founders --autosome --snps-only just-acgt --max-alleles 2 --rm-dup exclude-all --extract "$K/prune.prune.in" \
    --freq --pca 10 allele-wts vcols=chrom,ref,alt --out "$K/pca" --threads 8 > /dev/null
  # loadings: one row per variant allele (rsid, REF, ALT, ALT frequency, allele, PC1-10)
  awk 'NR == FNR { if (FNR > 1) fr[$2] = $3 "\t" $4 "\t" $5; next }
       FNR == 1 { printf "rsid\tref\talt\talt_freq\tallele"; for (i = 6; i <= 15; i++) printf "\tPC%d", i - 5; print ""; next }
       { printf "%s\t%s\t%s", $2, fr[$2], $5; for (i = 6; i <= 15; i++) printf "\t%s", $i; print "" }' \
    "$K/pca.afreq" "$K/pca.eigenvec.allele" | gzip -9 > "$OUT/ancestry_loadings.tsv.gz"
  # reference coordinates: the founders projected exactly as `mitrix prs` projects a cohort
  # (plink2 --score variance-standardize), so cohort and reference are on the same scale
  plink2 --pfile "$B/kg/all_hg38" vzs --allow-extra-chr --keep-founders --autosome --snps-only just-acgt --max-alleles 2 --rm-dup exclude-all --extract "$K/prune.prune.in" \
    --read-freq "$K/pca.afreq" --score "$K/pca.eigenvec.allele" 2 5 header-read \
    no-mean-imputation variance-standardize --score-col-nums 6-15 --out "$K/proj" --threads 8 > /dev/null
  awk 'NR == FNR { if (FNR > 1) pop[$1] = $5; next }
       FNR == 1 { c = 0; for (i = 1; i <= NF; i++) if ($i ~ /^PC[0-9]+_AVG$/) col[++c] = i
                  printf "sample\tsuperpop"; for (k = 1; k <= c; k++) printf "\tPC%d", k; print ""; next }
       { printf "%s\t%s", $1, pop[$1]; for (k = 1; k <= c; k++) printf "\t%s", $col[k]; print "" }' \
    "$B/kg/all_hg38.psam" "$K/proj.sscore" | gzip -9 > "$OUT/ancestry_reference.tsv.gz"
  log "ancestry reference: $(zcat "$OUT/ancestry_loadings.tsv.gz" | tail -n +2 | wc -l) SNPs"
  ;;
lead)
  need_tools plink2
  K="$W/kg"; mkdir -p "$K"
  awk 'NR > 1 && $5 == "EUR" && $2 == 0 && $3 == 0 { print $1 }' "$B/kg/all_hg38.psam" > "$K/eur.ids"
  each_score all | while IFS=: read -r s g n d; do
    plink2 --pfile "$B/kg/all_hg38" vzs --allow-extra-chr --keep "$K/eur.ids" --snps-only just-acgt --max-alleles 2 \
      --rm-dup exclude-all --clump "$W/$s.ldsc.txt" --clump-id-field SNP --clump-p-field P \
      --clump-p1 5e-8 --clump-r2 0.1 --clump-kb 1000 --out "$K/$s.lead" --threads 8 > /dev/null 2>&1 || true
    if [[ -s "$K/$s.lead.clumps" ]]; then
      awk 'NR == FNR { if (FNR > 1) lead[$3] = 1; next } FNR == 1 { print "rsid\teffect_allele\tweight"; next }
           ($1 in lead) { print $1 "\t" $2 "\t" $4 }' "$K/$s.lead.clumps" "$W/$s.ldsc.txt" \
        | gzip -9 > "$OUT/PRS_$s.lead.tsv.gz"
    else
      printf 'rsid\teffect_allele\tweight\n' | gzip -9 > "$OUT/PRS_$s.lead.tsv.gz"
    fi
    log "lead $s: $(( $(zcat "$OUT/PRS_$s.lead.tsv.gz" | wc -l) - 1 )) independent genome-wide significant SNPs"
  done
  ;;
table)
  # PRS-CS-auto estimates phi per chromosome: report the median
  phi_median() {
    sed -n 's/.*Estimated global shrinkage parameter: \([0-9.eE+-]*\).*/\1/p' "$W/prscs/$1".chr*.log \
      | sort -g | awk '{ v[NR] = $1 } END { if (!NR) print "NA"; else print (NR % 2 ? v[(NR + 1) / 2] : (v[NR / 2] + v[NR / 2 + 1]) / 2) }'
  }
  # expected accuracy (Daetwyler et al. 2008): R2 = h2 / (1 + Me / (N h2)), Me = 60,000
  # effectively independent segments (European ancestry)
  {
    printf 'score\tgwas_catalog\tdescription\tn_gwas\tn_snps\tn_lead\tphi\th2_snp\th2_se\texpected_r2\n'
    each_score all | while IFS=: read -r s g n d; do
      nsnp=$(zcat "$OUT/PRS_$s.weights.tsv.gz" | awk 'NR > 1 && $1 != 0' | wc -l)
      nlead=$(( $(zcat "$OUT/PRS_$s.lead.tsv.gz" | wc -l) - 1 ))
      read -r h2 se < <(sed -n 's/^Total Observed scale h2: \([-0-9.e]*\) (\([0-9.e]*\)).*/\1 \2/p' "$W/ldsc/$s.h2.log")
      r2=$(awk -v h="$h2" -v n="$n" 'BEGIN { if (h <= 0) print 0; else printf "%.4f", h / (1 + 60000 / (n * h)) }')
      printf 'PRS_%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$s" "$g" "$d" "$n" "$nsnp" "$nlead" \
        "$(phi_median "$s")" "$h2" "$se" "$r2"
    done
  } > "$OUT/scores.tsv"
  column -t -s $'\t' "$OUT/scores.tsv" >&2
  ;;
*)
  die "unknown step: $STEP (prep prscs merge ldsc ancestry lead table)"
  ;;
esac
