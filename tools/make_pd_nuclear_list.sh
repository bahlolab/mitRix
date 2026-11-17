#!/usr/bin/env bash
# tools/make_pd_nuclear_list.sh — example nuclear variant list for PD (GRCh38 and GRCh37)
# from ClinVar: presets/examples/pd_nuclear_variants_GRCh{38,37}.tsv
#
# Genes: SNCA, LRRK2, VPS35 (AD); PRKN, PINK1, PARK7 (AR); GBA1 (risk; heterozygous
# Gaucher-causing variants raise PD risk); APOE e4/e2 (rs429358, rs7412; not a PD gene,
# a dementia modifier; severity e4/e2 gives allele-count columns).
# Variants: ClinVar Pathogenic / Likely pathogenic (not conflicting) with review criteria
# (criteria provided, expert panel, practice guideline), plus variants ClinVar labels
# risk_factor with review criteria (GBA1 p.Glu365Lys "E326K", LRRK2 p.Gly2385Arg), plus
# moderate-effect risk variants that ClinVar classifies as conflicting or benign
# (GBA1 p.Thr408Met "T369M", LRRK2 p.Arg1628Pro, p.Asn2081Asp; CURATED_RISK below).
# severity: GBA1 variants graded by Gaucher severity as used for PD risk (Gan-Or et al.
# 2015, Neurology; Parlar et al. 2023, Mov Disord, GBA1-PD browser): severe, mild, risk,
# unclassified (other Gaucher P/LP); other genes: pathogenic or risk. mitrix carriers
# adds a column per severity class for genes with more than one class.
# Structural variants (e.g. PRKN exon deletions) are not in SNV/indel VCFs and are
# not listed. This is an EXAMPLE: curate your own list for analysis.
#
# Usage: tools/make_pd_nuclear_list.sh DIR
#   DIR holds clinvar_GRCh38.vcf.gz, clinvar_GRCh37.vcf.gz and variant_summary.txt.gz
#   from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/ (vcf_GRCh38/, vcf_GRCh37/, tab_delimited/)
# Needs bcftools, Rscript.

set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
D="${1:?usage: tools/make_pd_nuclear_list.sh CLINVAR_DIR}"
OUT="$HERE/presets/examples"; mkdir -p "$OUT"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
GENES="SNCA|LRRK2|VPS35|PRKN|PINK1|PARK7|GBA1|APOE"

for b in 38 37; do
  bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/GENEINFO\t%INFO/RS\n' \
    "$D/clinvar_GRCh$b.vcf.gz" 2> /dev/null \
    | awk -F'\t' -v g="$GENES" '$8 ~ "(^|\\|)(" g "):"' > "$TMP/vcf$b.tsv"
done
zcat "$D/variant_summary.txt.gz" | awk -F'\t' -v OFS='\t' 'NR > 1 && $17 == "GRCh38" { print $31, $3 }' > "$TMP/names.tsv"
DATE=$(bcftools view -h "$D/clinvar_GRCh38.vcf.gz" | awk -F= '$1 == "##fileDate" { print $2 }')

Rscript --vanilla - "$TMP" "$OUT" "$DATE" <<'EOF'
a <- commandArgs(TRUE); tmp <- a[1]; out <- a[2]; date <- a[3]
cols <- c("chr", "pos", "id", "ref", "alt", "clnsig", "review", "geneinfo", "rs")
rd <- function(f) read.delim(f, header = FALSE, col.names = cols, colClasses = "character", quote = "")
inh <- c(SNCA = "AD", LRRK2 = "AD", VPS35 = "AD", PRKN = "AR", PINK1 = "AR", PARK7 = "AR", GBA1 = "risk", APOE = "risk")
nm <- read.delim(file.path(tmp, "names.tsv"), header = FALSE, col.names = c("id", "name"), colClasses = "character", quote = "")
nm <- nm[!duplicated(nm$id), ]
v38 <- rd(file.path(tmp, "vcf38.tsv"))
reviewed <- grepl("criteria_provided|reviewed_by_expert_panel|practice_guideline", v38$review) &
            !grepl("no_assertion|conflicting_classifications", v38$review)
plp <- grepl("^(Pathogenic|Likely_pathogenic|Pathogenic/Likely_pathogenic)(\\||$)", v38$clnsig) & reviewed
risk <- grepl("risk_factor", v38$clnsig) & grepl("criteria_provided|reviewed_by_expert_panel", v38$review) & !grepl("no_assertion", v38$review)
# moderate-effect risk variants not P/LP or risk_factor in ClinVar (variation IDs)
CURATED_RISK <- c("93447",    # GBA1 p.Thr408Met (T369M)
                  "39198",    # LRRK2 p.Arg1628Pro
                  "39221")    # LRRK2 p.Asn2081Asp
cur <- v38$id %in% CURATED_RISK
# APOE e4/e2 (not PD genes; dementia/cognitive-decline modifier)
APOE_IDS <- c(e4 = "17864",   # rs429358 T>C (388T>C): C allele = e4
              e2 = "17848")   # rs7412 C>T (526C>T): T allele = e2
apoe <- v38$id %in% APOE_IDS & v38$alt != "."
keep <- v38[(plp | risk | cur | apoe) & v38$alt != ".", ]
keep$gene <- vapply(strsplit(keep$geneinfo, "|", fixed = TRUE), function(x) {
  g <- sub(":.*", "", x); g <- g[g %in% names(inh)]; if (length(g)) g[1] else NA_character_ }, "")
keep <- keep[!is.na(keep$gene) & (keep$gene != "APOE" | keep$id %in% APOE_IDS), ]
n <- nm$name[match(keep$id, nm$id)]
p <- sub(".*\\((p\\.[^)]*)\\).*", "\\1", n); cdot <- sub("^[^:]*:(c\\.[^ ]*).*", "\\1", n)
keep$label <- ifelse(grepl("\\(p\\.", n), p, ifelse(grepl(":c\\.", n), cdot, paste0("ClinVar:", keep$id)))
ia <- keep$id %in% APOE_IDS; ae <- names(APOE_IDS)[match(keep$id, APOE_IDS)]
keep$label[ia] <- paste0("rs", keep$rs[ia], " (", ae[ia], ")")
keep$classification <- ifelse(ia, "dementia risk modifier (not a PD gene)", ifelse(keep$id %in% CURATED_RISK, "risk_factor (curated)",
                         ifelse(grepl("^(Pathogenic|Likely_pathogenic|Pathogenic/Likely_pathogenic)", keep$clnsig) & !grepl("Conflicting", keep$clnsig),
                              sub("\\|.*", "", keep$clnsig), "risk_factor")))
# GBA1 severity (legacy names in comments; HGVS numbering = legacy + 39)
gba_severe <- c("p.Leu483Pro",  # L444P
                "p.Asp448His",  # D409H
                "p.Val433Leu",  # V394L
                "p.Arg159Trp",  # R120W
                "p.Pro454Arg",  # P415R
                "p.Ala495Pro",  # A456P (RecNciI)
                "p.Leu29fs",    # 84GG (c.84dup)
                "c.115+1G>A")   # IVS2+1G>A
gba_mild <- c("p.Asn409Ser",    # N370S
              "p.Arg535His")    # R496H
keep$severity <- ifelse(ia, ae, ifelse(grepl("^risk_factor", keep$classification), "risk",
                 ifelse(keep$gene != "GBA1", "pathogenic",
                 ifelse(keep$label %in% gba_severe, "severe",
                 ifelse(keep$label %in% gba_mild, "mild", "unclassified")))))
keep$inheritance <- unname(inh[keep$gene])
keep$rsid <- ifelse(keep$rs == ".", "", paste0("rs", keep$rs))
keep$source <- paste0("ClinVar ", date, " VCV", keep$id, " (", gsub("_", " ", keep$review), ")")
keep$source[keep$id %in% CURATED_RISK] <- paste0(keep$source[keep$id %in% CURATED_RISK], "; added as moderate-effect PD risk variant")
hdr <- c("# Example PD nuclear variant list (mitRix tools/make_pd_nuclear_list.sh) from ClinVar " ,
         "# P/LP with review criteria + risk variants; genes SNCA LRRK2 VPS35 (AD), PRKN PINK1 PARK7 (AR), GBA1 (risk; severity graded); APOE e4/e2 (dementia modifier).",
         "# Not included: structural variants (e.g. PRKN exon deletions). Example only: curate your own list.")
write_list <- function(d, build, f, h = hdr) {
  d <- d[order(d$gene, as.integer(d$pos)), ]
  o <- data.frame(chr = paste0("chr", d$chr), pos = d$pos, ref = d$ref, alt = d$alt, gene = d$gene, label = d$label,
                  inheritance = d$inheritance, classification = d$classification, severity = d$severity,
                  rsid = d$rsid, source = d$source)
  con <- file(f, "w"); writeLines(c(h[1], paste0("# ", build, ". ", sub("^# ", "", h[2])), h[3]), con)
  write.table(o, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con)
  message(sprintf("%s: %d variants -> %s", build, nrow(o), f))
}
write_list(keep, "GRCh38", file.path(out, "pd_nuclear_variants_GRCh38.tsv"))
# GRCh37: same ClinVar variation IDs
v37 <- rd(file.path(tmp, "vcf37.tsv"))
k37 <- merge(keep[, c("id", "gene", "label", "classification", "severity", "inheritance", "rsid", "source")],
             v37[, c("id", "chr", "pos", "ref", "alt")], by = "id")
if (nrow(k37) < nrow(keep)) message("GRCh37: ", nrow(keep) - nrow(k37), " variant(s) without GRCh37 position")
write_list(k37, "GRCh37", file.path(out, "pd_nuclear_variants_GRCh37.tsv"))
print(table(keep$gene, keep$severity))
EOF
