#!/usr/bin/env bash
# tests/test_prs.sh — end-to-end test of `mitrix prs` on synthetic data with known answers:
# 60 samples x 300 SNPs (chr1-2), toy weights for two scores, a toy ancestry reference
# (all EUR). Checks the score sums against values computed here, strand-flipped and
# allele-swapped weights, missing genotypes (mean-imputed), GRCh37 vs GRCh38 positions,
# the overlap fraction, and the long-format table. Needs plink2, Rscript, python3.

set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
TMP=$(mktemp -d -t mitrix_prs_test.XXXXXX); trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/presets"

python3 - "$TMP" <<'EOF'
import gzip, random, sys
t = sys.argv[1]; random.seed(7)
NS, NV = 60, 300
comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
pairs = [("A", "C"), ("A", "G"), ("C", "T"), ("G", "T")]   # no strand-ambiguous pairs
samples = [f"S{i:02d}" for i in range(NS)]
var = []
for j in range(NV):
    chrom = 1 if j < 150 else 2
    p38 = 1_000_000 + j * 5000; p37 = p38 - 12345
    ref, alt = random.choice(pairs)
    f = random.uniform(0.2, 0.8)
    g = [sum(random.random() < f for _ in range(2)) for _ in range(NS)]
    miss = set(random.sample(range(NS), 2)) if j % 10 == 0 else set()   # 2 missing calls per 10th SNP
    var.append(dict(rs=f"rs{100 + j}", chrom=chrom, p37=p37, p38=p38, ref=ref, alt=alt, g=g, miss=miss))

for build, key in (("38", "p38"), ("37", "p37")):
    with open(f"{t}/cohort{build}.vcf", "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for c in (1, 2):
            f.write(f"##contig=<ID=chr{c},length=248956422>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        for v in var:
            gts = ["./." if i in v["miss"] else ["0/0", "0/1", "1/1"][v["g"][i]] for i in range(NS)]
            f.write(f"chr{v['chrom']}\t{v[key]}\t.\t{v['ref']}\t{v['alt']}\t.\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")

# map: first 280 SNPs (20 SNPs absent from the map); effect allele: ALT, REF (swapped) or
# the complement (strand flip) by turn; plus 20 map SNPs missing from the VCF -> overlap 280/300
w = {"A": [], "B": []}
with gzip.open(f"{t}/presets/snps.tsv.gz", "wt") as f:
    f.write("rsid\tchr\tpos_grch37\tpos_grch38\teffect_allele\tother_allele\n")
    for j, v in enumerate(var[:280]):
        mode = j % 3
        e, o = (v["alt"], v["ref"]) if mode == 0 else (v["ref"], v["alt"]) if mode == 1 else (comp[v["alt"]], comp[v["ref"]])
        f.write(f"{v['rs']}\t{v['chrom']}\t{v['p37']}\t{v['p38']}\t{e}\t{o}\n")
        wa, wb = round(random.gauss(0, 1), 3), round(random.gauss(0, 1), 3) if j % 2 else 0.0
        w["A"].append(wa); w["B"].append(wb)
        # expected contribution per sample, in terms of the ALT dosage
        v["eff_alt"] = mode != 1
        v["wA"], v["wB"] = wa, wb
    for k in range(20):
        f.write(f"rs9{k:03d}\t2\t{5_000_000 + k}\t{5_000_000 + k}\tA\tG\n")
        w["A"].append(1.0); w["B"].append(1.0)
for s in ("A", "B"):
    with gzip.open(f"{t}/presets/PRS_{s}.weights.tsv.gz", "wt") as f:
        f.write("weight\n" + "".join(f"{x}\n" for x in w[s]))
with open(f"{t}/presets/scores.tsv", "w") as f:
    f.write("score\tgwas_catalog\n" + "PRS_A\ttoy\nPRS_B\ttoy\n")

# ancestry: loadings on the first 60 SNPs, reference = 25 EUR samples
with gzip.open(f"{t}/presets/ancestry_loadings.tsv.gz", "wt") as f:
    f.write("rsid\tref\talt\talt_freq\tallele" + "".join(f"\tPC{k}" for k in range(1, 11)) + "\n")
    for v in var[:60]:
        for a in (v["ref"], v["alt"]):
            f.write(f"{v['rs']}\t{v['ref']}\t{v['alt']}\t0.5\t{a}" + "".join(f"\t{random.gauss(0, 0.1):.4f}" for _ in range(10)) + "\n")
with gzip.open(f"{t}/presets/ancestry_reference.tsv.gz", "wt") as f:
    f.write("sample\tsuperpop" + "".join(f"\tPC{k}" for k in range(1, 11)) + "\n")
    for i in range(25):
        f.write(f"R{i}\tEUR" + "".join(f"\t{random.gauss(0, 0.01):.4f}" for _ in range(10)) + "\n")

# expected sums: plink2 mean-imputes a missing genotype as 2 x (effect-allele frequency
# among non-missing samples)
with open(f"{t}/expected.tsv", "w") as f:
    f.write("sample\tPRS_A\tPRS_B\n")
    for i, s in enumerate(samples):
        tot = {"A": 0.0, "B": 0.0}
        for v in var[:280]:
            obs = [v["g"][k] for k in range(NS) if k not in v["miss"]]
            alt_d = v["g"][i] if i not in v["miss"] else 2 * sum(obs) / (2 * len(obs))
            d = alt_d if v["eff_alt"] else 2 - alt_d
            tot["A"] += v["wA"] * d; tot["B"] += v["wB"] * d
        f.write(f"{s}\t{tot['A']:.6f}\t{tot['B']:.6f}\n")
EOF

export MITRIX_PRS_DIR="$TMP/presets"
fail=0
for build in 38 37; do
  bgzip -c "$TMP/cohort$build.vcf" > "$TMP/cohort$build.vcf.gz" 2>/dev/null || gzip -c "$TMP/cohort$build.vcf" > "$TMP/cohort$build.vcf.gz"
  "$HERE/bin/mitrix" prs --vcf "$TMP/cohort$build.vcf.gz" --build "$build" --outdir "$TMP/out$build" \
    --maf 0 --hwe 0 --geno 0.1 --threads 1 2> "$TMP/log$build" || { cat "$TMP/log$build"; exit 1; }
  O="$TMP/out$build/prs"
  # score sums vs expected
  bad=$(awk -F'\t' 'NR == FNR { if (FNR > 1) { a[$1] = $2; b[$1] = $3 }; next }
    $3 == "PRS_A_sum" { d = $4 - a[$1]; if (d < 0) d = -d; if (d > 1e-3) print "PRS_A", $1, $4, a[$1] }
    $3 == "PRS_B_sum" { d = $4 - b[$1]; if (d < 0) d = -d; if (d > 1e-3) print "PRS_B", $1, $4, b[$1] }' \
    "$TMP/expected.tsv" "$O/prs.tsv")
  if [[ -n "$bad" ]]; then echo "FAIL GRCh$build scores:"; echo "$bad" | head -5; fail=1
  else echo "ok   GRCh$build: 60 x 2 score sums match"; fi
  # overlap: non-zero weights in the map vs found in the VCF (the 20 extra map SNPs are absent)
  nz=$(zcat "$TMP/presets/PRS_A.weights.tsv.gz" | awk 'NR > 1 && $1 != 0' | wc -l)
  nz_found=$(zcat "$TMP/presets/PRS_A.weights.tsv.gz" | awk 'NR > 1 && NR <= 281 && $1 != 0' | wc -l)
  ov=$(awk -F'\t' '$1 == "PRS_A" { print $2, $3 }' "$O/match.tsv")
  [[ "$ov" == "$nz $nz_found" ]] && echo "ok   GRCh$build: PRS_A overlap $nz_found/$nz" || { echo "FAIL GRCh$build overlap: $ov (expected $nz $nz_found)"; fail=1; }
  fl=$(awk -F'\t' '$1 == "PRS_A" { print $5 }' "$O/match.tsv")
  [[ "$fl" -ge 90 ]] && echo "ok   GRCh$build: $fl strand flips handled" || { echo "FAIL GRCh$build flips: $fl"; fail=1; }
  # 280/300 = 0.933 >= 0.8 -> PASS; all samples EUR
  q=$(awk -F'\t' '$3 == "PRS_A" { print $6 }' "$O/prs.tsv" | sort -u)
  [[ "$q" == PASS ]] && echo "ok   GRCh$build: qc_flag PASS" || { echo "FAIL GRCh$build qc_flag: $q"; fail=1; }
  n=$(awk -F'\t' '$3 == "ancestry" && $4 == "EUR"' "$O/prs.tsv" | wc -l)
  [[ "$n" == 60 ]] && echo "ok   GRCh$build: 60 EUR" || { echo "FAIL GRCh$build ancestry: $n EUR"; fail=1; }
  # z-scores: mean 0 over EUR
  m=$(awk -F'\t' '$3 == "PRS_A" { s += $4; n++ } END { printf "%.6f", s / n }' "$O/prs.tsv")
  awk -v m="$m" 'BEGIN { exit !(m < 1e-4 && m > -1e-4) }' && echo "ok   GRCh$build: z mean 0" || { echo "FAIL z mean $m"; fail=1; }
done
# --min-overlap above 280/300 -> low_overlap flag
"$HERE/bin/mitrix" prs --vcf "$TMP/cohort38.vcf.gz" --outdir "$TMP/outov" --maf 0 --hwe 0 --geno 0.1 \
  --threads 1 --min-overlap 0.95 2> /dev/null
q=$(awk -F'\t' '$3 == "PRS_A" { print $6 }' "$TMP/outov/prs/prs.tsv" | sort -u)
[[ "$q" == low_overlap ]] && echo "ok   low_overlap flag" || { echo "FAIL low_overlap: $q"; fail=1; }

(( fail )) && { echo FAILED; exit 1; }
echo "all prs tests passed"
