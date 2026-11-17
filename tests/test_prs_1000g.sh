#!/usr/bin/env bash
# tests/test_prs_1000g.sh — release test of `mitrix prs` on public 1000 Genomes data
# (Slurm: -c 8 --mem 32G, ~30 min). Needs plink2, Rscript, bcftools, python3.
#
#   KG38=prefix  1000 Genomes high-coverage GRCh38 pfile (pgen + pvar.zst + psam)
#   KG37=prefix  1000 Genomes phase 3 GRCh37 pfile
#   WORK=dir     scratch directory
#
# 1. GRCh38, 3,202 samples: match rate, ancestry vs known super-population, runtime
# 2. GRCh37 (no 'chr' prefix), 2,504 samples: score sums vs the GRCh38 run (r > 0.95)
# 3. robustness (GRCh38): 5% genotypes set missing; 300-sample subset (r > 0.99)
# 4. validation: simulated readouts for EUR samples (PRS_CN-driven r ~ 0.3 and null),
#    with --covariates including a mock cell-composition column

set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
: "${KG38:?set KG38}" "${KG37:?set KG37}" "${WORK:?set WORK}"
mkdir -p "$WORK"; cd "$WORK"
M="$HERE/bin/mitrix"; P="$HERE/presets/prs"
zcat "$P/snps.tsv.gz" | awk 'NR > 1 { print $1 }' > ids
fail=0
ok()  { echo "ok   $*"; }
bad() { echo "FAIL $*"; fail=1; }
corr() {  # file1 file2 metric -> Pearson r of the metric over shared samples
  Rscript --vanilla -e 'a <- commandArgs(TRUE); m <- a[3]
    r <- function(f) { d <- read.delim(f); d <- d[d$metric == m, c("sample", "value")]; d$value <- as.numeric(d$value); d }
    x <- merge(r(a[1]), r(a[2]), by = "sample"); cat(sprintf("%.4f %d\n", cor(x$value.x, x$value.y), nrow(x)))' "$1" "$2" "$3"
}
SCORES=$(awk -F'\t' 'NR > 1 { print $1 }' "$P/scores.tsv")

# ---- export VCFs
[[ -s kg38.vcf.gz ]] || plink2 --pfile "$KG38" vzs --allow-extra-chr --autosome --snps-only just-acgt --max-alleles 2 \
  --rm-dup exclude-all --extract ids --output-chr chr26 --export vcf bgz id-paste=iid --out kg38 --threads 8 > /dev/null
[[ -s kg37.vcf.gz ]] || plink2 --pfile "$KG37" vzs --allow-extra-chr --autosome --snps-only just-acgt --max-alleles 2 \
  --rm-dup exclude-all --extract ids --output-chr 26 --export vcf bgz id-paste=iid --out kg37 --threads 8 > /dev/null

# ---- 1. GRCh38
start=$(date +%s)
"$M" prs --vcf kg38.vcf.gz --outdir t38 --threads 8 2> t38.log
echo "     GRCh38 runtime: $(( $(date +%s) - start )) s"
awk -F'\t' 'NR > 1 { print "     " $1 ": " $3 "/" $2 " (" $4 ")" }' t38/prs/match.tsv
low=$(awk -F'\t' 'NR > 1 && $4 < 0.95' t38/prs/match.tsv | wc -l)
(( low == 0 )) && ok "GRCh38: >= 95% of weight SNPs found for every score" || bad "GRCh38: low overlap"
acc=$(awk -F'\t' -v OFS='\t' 'NR == FNR { if (FNR > 1) sp[$1] = $5; next } FNR > 1 { n++; if ($2 == sp[$1]) k++; if (sp[$1] == "EUR") { e++; if ($2 == "EUR") ek++ } }
  END { printf "%d %d %d %d", k, n, ek, e }' "$KG38.psam" t38/prs/ancestry.tsv)
read -r k n ek e <<< "$acc"
echo "     ancestry: $k/$n correct; EUR $ek/$e"
awk -v k="$k" -v n="$n" -v ek="$ek" -v e="$e" 'BEGIN { exit !(k / n > 0.98 && ek == e) }' \
  && ok "ancestry > 98% correct, all EUR assigned EUR" || bad "ancestry accuracy"
nlead=$(awk -F'\t' '$3 ~ /_lead$/' t38/prs/prs.tsv | cut -f3 | sort -u | wc -l)
(( nlead >= 1 )) && ok "lead-variant scores present ($nlead)" || bad "no lead-variant scores"

# ---- 2. GRCh37 vs GRCh38
"$M" prs --vcf kg37.vcf.gz --build 37 --outdir t37 --threads 8 2> t37.log
for s in $SCORES; do
  read -r r nn <<< "$(corr t38/prs/prs.tsv t37/prs/prs.tsv "${s}_sum")"
  echo "     $s sum GRCh38 vs GRCh37: r = $r (n = $nn)"
  awk -v r="$r" 'BEGIN { exit !(r > 0.95) }' && ok "$s GRCh37 concordant" || bad "$s GRCh37 r = $r"
done

# ---- 3. robustness
[[ -s kg38_miss.vcf.gz ]] || python3 - <<'EOF'
import gzip, random
random.seed(1)
with gzip.open("kg38.vcf.gz", "rt") as f, open("kg38_miss.vcf", "w") as o:
    for line in f:
        if line.startswith("#"):
            o.write(line); continue
        x = line.rstrip("\n").split("\t")
        for i in range(9, len(x)):
            if random.random() < 0.05:
                x[i] = "./."
        o.write("\t".join(x) + "\n")
EOF
[[ -s kg38_miss.vcf.gz ]] || bgzip -f kg38_miss.vcf
"$M" prs --vcf kg38_miss.vcf.gz --outdir tmiss --threads 8 --geno 0.1 --mind 0.1 2> tmiss.log
bcftools query -l kg38.vcf.gz > all.ids; shuf -n 300 --random-source=<(yes) all.ids > sub.ids
bcftools view -S sub.ids -Oz -o kg38_sub.vcf.gz kg38.vcf.gz
"$M" prs --vcf kg38_sub.vcf.gz --outdir tsub --threads 8 2> tsub.log
for t in tmiss tsub; do
  for s in $SCORES; do
    read -r r nn <<< "$(corr t38/prs/prs.tsv $t/prs/prs.tsv "${s}_sum")"
    awk -v r="$r" 'BEGIN { exit !(r > 0.99) }' && ok "$t $s r = $r (n = $nn)" || bad "$t $s r = $r"
  done
done

# ---- 4. validation with simulated readouts (EUR)
Rscript --vanilla - <<'EOF'
set.seed(42)
d <- read.delim("t38/prs/prs.tsv")
z <- d[d$metric == "PRS_CN", c("sample", "value")]; names(z)[2] <- "z"; z$z <- as.numeric(z$z)
anc <- d[d$metric == "ancestry", c("sample", "value")]; names(anc)[2] <- "anc"
x <- merge(z, anc); x <- x[x$anc == "EUR", ]; n <- nrow(x)
age <- round(runif(n, 40, 80)); sex <- sample(c("F", "M"), n, TRUE); neut <- runif(n, 0.4, 0.7)
cohort <- data.frame(sample = x$sample,
  mtDNA_CN = 100 + 10 * (0.3 * x$z + sqrt(1 - 0.09) * rnorm(n)) - 0.2 * age,
  age_accumulating = rpois(n, 2), flap_7S = rnorm(n, 0.9, 0.05), flap_primer = rnorm(n, 1.4, 0.1),
  check.names = FALSE)
write.table(cohort, "sim_cohort.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample = x$sample, age, sex, neutrophil_frac = neut), "sim_meta.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
EOF
"$M" prs --vcf kg38.vcf.gz --outdir tval --threads 8 --cohort sim_cohort.tsv --metadata sim_meta.tsv \
  --covariates age,sex,neutrophil_frac 2> tval.log
column -t tval/prs/validation.tsv | sed 's/^/     /'
r_cn=$(awk -F'\t' '$1 == "PRS_CN" { print $4 }' tval/prs/validation.tsv)
r_7s=$(awk -F'\t' '$1 == "PRS_7S" { print $4 }' tval/prs/validation.tsv)
awk -v r="$r_cn" 'BEGIN { exit !(r > 0.2 && r < 0.4) }' && ok "validation recovers PRS_CN effect (r = $r_cn)" || bad "PRS_CN validation r = $r_cn"
awk -v r="$r_7s" 'BEGIN { exit !(r > -0.12 && r < 0.12) }' && ok "validation null for PRS_7S (r = $r_7s)" || bad "PRS_7S null r = $r_7s"
grep -q 'neutrophil_frac' tval/prs/logs/table.log && ok "validation used neutrophil_frac" || bad "covariate not used"

(( fail )) && { echo FAILED; exit 1; }
echo "all 1000 Genomes PRS tests passed"
