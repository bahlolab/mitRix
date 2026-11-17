#!/usr/bin/env python3
"""Harmonise one GWAS Catalog harmonised summary-statistics file for PRS-CS and LDSC.

usage: prep_sumstats.py GWAS.h.tsv.gz SNPINFO OUT_PREFIX

GWAS: GWAS Catalog harmonised file (GRCh38; columns chromosome, base_pair_location,
      effect_allele, other_allele, beta, standard_error, effect_allele_frequency,
      p_value, rsid).
SNPINFO: PRS-CS LD panel snpinfo (CHR SNP BP A1 A2 MAF; BP in GRCh37), e.g.
      ldblk_1kg_eur/snpinfo_1kg_hm3.

Keeps biallelic SNVs (A/C/G/T) whose rsID is in the LD panel with the same alleles
(either orientation), not strand-ambiguous (A/T, C/G), with effect-allele frequency
0.01-0.99 and finite beta/SE. Writes:
  OUT_PREFIX.prscs.txt   SNP A1 A2 BETA SE        (PRS-CS input; A1 = effect allele)
  OUT_PREFIX.ldsc.txt    SNP A1 A2 BETA SE P      (for ldsc munge_sumstats)
  OUT_PREFIX.pos38.txt   SNP CHR POS38            (GRCh38 positions for the SNP map)
  OUT_PREFIX.prep.log    counts
"""
import gzip, math, sys

gwas_file, snpinfo_file, out = sys.argv[1:4]
AMBIG = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
ACGT = set("ACGT")

ref = {}
with open(snpinfo_file) as f:
    next(f)
    for line in f:
        chrom, snp, bp, a1, a2, maf = line.split()
        ref[snp] = (chrom, a1, a2)

n = {"rows": 0, "no_rsid": 0, "not_snv": 0, "ambiguous": 0, "not_in_ld_panel": 0,
     "allele_mismatch": 0, "chrom_mismatch": 0, "maf": 0, "bad_stats": 0, "duplicate": 0, "kept": 0}
seen = set()
with gzip.open(gwas_file, "rt") as f, open(out + ".prscs.txt", "w") as fp, \
        open(out + ".ldsc.txt", "w") as fl, open(out + ".pos38.txt", "w") as fpos:
    header = f.readline().rstrip("\n").split("\t")
    c = {h: i for i, h in enumerate(header)}
    for need in ("chromosome", "base_pair_location", "effect_allele", "other_allele",
                 "beta", "standard_error", "effect_allele_frequency", "p_value", "rsid"):
        if need not in c:
            sys.exit(f"{gwas_file}: missing column {need}")
    fp.write("SNP\tA1\tA2\tBETA\tSE\n")
    fl.write("SNP\tA1\tA2\tBETA\tSE\tP\n")
    fpos.write("SNP\tCHR\tPOS38\n")
    for line in f:
        n["rows"] += 1
        x = line.rstrip("\n").split("\t")
        rsid = x[c["rsid"]]
        if not rsid.startswith("rs"):
            n["no_rsid"] += 1; continue
        ea, oa = x[c["effect_allele"]].upper(), x[c["other_allele"]].upper()
        if len(ea) != 1 or len(oa) != 1 or ea not in ACGT or oa not in ACGT:
            n["not_snv"] += 1; continue
        if (ea, oa) in AMBIG:
            n["ambiguous"] += 1; continue
        if rsid not in ref:
            n["not_in_ld_panel"] += 1; continue
        rchrom, ra1, ra2 = ref[rsid]
        if {ea, oa} != {ra1, ra2}:
            n["allele_mismatch"] += 1; continue
        chrom = x[c["chromosome"]]
        if chrom != rchrom:
            n["chrom_mismatch"] += 1; continue
        try:
            beta, se = float(x[c["beta"]]), float(x[c["standard_error"]])
            eaf = float(x[c["effect_allele_frequency"]])
            p = float(x[c["p_value"]])
        except ValueError:
            n["bad_stats"] += 1; continue
        if not (math.isfinite(beta) and math.isfinite(se) and se > 0 and 0 <= p <= 1):
            n["bad_stats"] += 1; continue
        if not (0.01 <= eaf <= 0.99):
            n["maf"] += 1; continue
        if rsid in seen:
            n["duplicate"] += 1; continue
        seen.add(rsid)
        fp.write(f"{rsid}\t{ea}\t{oa}\t{beta:.6g}\t{se:.6g}\n")
        fl.write(f"{rsid}\t{ea}\t{oa}\t{beta:.6g}\t{se:.6g}\t{p:.6g}\n")
        fpos.write(f"{rsid}\t{chrom}\t{x[c['base_pair_location']]}\n")
        n["kept"] += 1

with open(out + ".prep.log", "w") as f:
    for k, v in n.items():
        f.write(f"{k}\t{v}\n")
print(" ".join(f"{k}={v}" for k, v in n.items()), file=sys.stderr)
