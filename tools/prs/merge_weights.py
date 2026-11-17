#!/usr/bin/env python3
"""Merge PRS-CS per-chromosome outputs of several scores into one shared SNP map plus one
weight column file per score (same row order).

usage: merge_weights.py SNPINFO OUTDIR SCORE=PRSCS_GLOB:POS38_FILE [SCORE=...]

SNPINFO   PRS-CS LD panel snpinfo (CHR SNP BP A1 A2 MAF; BP = GRCh37)
PRSCS_GLOB  PRS-CS outputs for the score (CHR SNP BP A1 A2 BETA), e.g. 'dir/CN_pst_*.txt'
POS38_FILE  SNP CHR POS38 from prep_sumstats.py
Writes OUTDIR/snps.tsv.gz (rsid chr pos_grch37 pos_grch38 effect_allele other_allele)
and OUTDIR/<SCORE>.weights.tsv.gz (weight for effect_allele; 0 if the SNP is not in
that score). Map rows: SNPs with a weight in any score, sorted by chr and GRCh37 position.
"""
import glob, gzip, sys

snpinfo, outdir = sys.argv[1:3]
specs = [a.split("=", 1) for a in sys.argv[3:]]

info = {}
with open(snpinfo) as f:
    next(f)
    for line in f:
        chrom, snp, bp, a1, a2, maf = line.split()
        info[snp] = (int(chrom), int(bp), a1, a2)

pos38, weights = {}, {}
for score, spec in specs:
    pattern, pos_file = spec.split(":", 1)
    with open(pos_file) as f:
        next(f)
        for line in f:
            snp, chrom, pos = line.split()
            pos38.setdefault(snp, pos)
    w = {}
    files = sorted(glob.glob(pattern))
    if not files:
        sys.exit(f"{score}: no PRS-CS output matches {pattern}")
    for fn in files:
        with open(fn) as f:
            for line in f:
                chrom, snp, bp, a1, a2, beta = line.split()
                ia1, ia2 = info[snp][2], info[snp][3]
                b = float(beta)
                if (a1, a2) == (ia1, ia2):
                    w[snp] = b
                elif (a1, a2) == (ia2, ia1):
                    w[snp] = -b
                else:
                    sys.exit(f"{score}: {snp} alleles {a1}/{a2} do not match the LD panel {ia1}/{ia2}")
    weights[score] = w
    print(f"{score}: {len(w)} SNPs from {len(files)} files", file=sys.stderr)

snps = sorted({s for w in weights.values() for s in w}, key=lambda s: (info[s][0], info[s][1]))
missing38 = [s for s in snps if s not in pos38]
if missing38:
    sys.exit(f"{len(missing38)} SNPs without a GRCh38 position, e.g. {missing38[:3]}")
with gzip.open(f"{outdir}/snps.tsv.gz", "wt", compresslevel=9) as f:
    f.write("rsid\tchr\tpos_grch37\tpos_grch38\teffect_allele\tother_allele\n")
    for s in snps:
        chrom, bp, a1, a2 = info[s]
        f.write(f"{s}\t{chrom}\t{bp}\t{pos38[s]}\t{a1}\t{a2}\n")
for score, w in weights.items():
    with gzip.open(f"{outdir}/{score}.weights.tsv.gz", "wt", compresslevel=9) as f:
        f.write("weight\n")
        for s in snps:
            f.write(f"{w.get(s, 0.0):.4g}\n")
print(f"map: {len(snps)} SNPs", file=sys.stderr)
