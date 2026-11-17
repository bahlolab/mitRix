# Third-party data bundled with mitRix

mitRix code is MIT-licensed (see `LICENSE`). The reference data in `presets/` and
`tools/pathogenic/` were derived from the sources below, which keep their own terms.
Please cite them when you use the corresponding metrics.

| File | Source | Terms |
|---|---|---|
| `presets/mt_pathogenic.tsv`, `tools/pathogenic/mitomap_*.csv` | MITOMAP confirmed pathogenic variants (Lott et al. 2013, PMID 25489354; https://www.mitomap.org) | CC BY 3.0 |
| `presets/mt_pathogenic.tsv` | ClinVar pathogenic / likely pathogenic (NCBI) | public domain |
| `presets/numt_alleles.tsv` | MitoHPC `RefSeq/NUMT.vcf.gz` (https://github.com/dpuiu/MitoHPC) | MIT, Copyright (c) MitoHPC authors |
| `presets/mt_annotation.tsv.gz` | Lake et al. 2024, *Nature* 635:390 (MLC scores, VEP consequences); MITOMAP locus boundaries via leklab/mitochondrial_constraint | open access; MIT |
| `presets/examples/pd_nuclear_variants_*.tsv` | ClinVar (NCBI) | public domain |
| `presets/prs/*` | Weights from GWAS Catalog summary statistics of Gupta et al. 2023 and 2026 (mtDNA copy number, heteroplasmy, 7S DNA traits), built with PRS-CS and the 1000 Genomes EUR LD panel; ancestry reference from 1000 Genomes | per GWAS Catalog / 1000 Genomes terms |
| `refs/chrMC_rCRS.fa` | rCRS, NC_012920.1 (Andrews et al. 1999) | public domain |
