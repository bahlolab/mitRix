# Changelog

## 0.7.1 — 2026-11-17

### New

- **`mitrix cohort --cn-clusters auto|off|K`** (default `auto`): mtDNA copy-number clusters
  across the cohort (possible mixed DNA sources or batches). Gaussian mixture (mclust) on
  log mtDNA-CN, K = 1-4 chosen by ICL with a minimum cluster size; columns `cn_cluster`
  (1 = lowest CN) and `cn_cluster_prob`, soft flag `cn_cluster:uncertain` (prob < 0.9),
  `cohort_cn_clusters.tsv` and `_fit.tsv`. Skipped with a warning without mclust.
  Test: `tests/test_cn_cluster.sh`.

## 0.7.0 — 2026-11-17

### New

- **coverage_shape 0.1.0**: per-sample robust model of the chrM coverage shape (100-bp bins of
  the corrected per-base, control region excluded): `log depth ~ GC + GC^2 + e + deletion span`,
  e = L-strand distance from O_L (strand-displacement replication). Metrics `gc_slope`,
  `replication_gradient` (> 1 = replication-intermediate excess), `deletion_step`,
  GC-adjusted `ND4_ND1_ratio_gc_adj` / `MAJOR_MINOR_ratio_gc_adj`, `local_dip` (+ position, z).
- **damage 0.3.0**: `MAJOR_MINOR` ratio (common-deletion span 8470-13447 / ND1-ND2 3307-5511).
- **tools/from_perbase.sh**: damage, d_loop_7s, coverage_shape from precomputed per-base files
  (no BAM).

### Changed

- Licence: MIT (`LICENSE`); third-party data sources and terms in `NOTICE.md`.
- The v0.1.0 single script (`mitRix`), its WDL, Dockerfile and Docker GitHub Action moved
  to `legacy/v0.1.0/`. Empty `containers/`, `nextflow/`, `wdl/` placeholders removed.
- Removed cohort-specific example lists; PD example list `presets/examples/pd_nuclear_variants_*`
  unchanged apart from source notes.

## 0.6.1 — 2026-11-09

### New

- **Variants-only mode**: `mitrix run --mt-variants FILE` without `--bam` (samplesheet:
  empty or no `bam` column) runs haplogroup, contamination, heteroplasmy and
  pathogenic from the variant calls alone; chrM depth (detection limit, low-depth
  flags) is the median DP of the calls, recorded as `depth_source` in the variant
  params. copy_number, damage and d_loop_7s are SKIPPED (status `SKIPPED`; cohort
  columns NA). HG00096 (MitoHPC, mity, mtDNA-Server 2): all variant-based metrics equal
  the BAM-based run except `mt_median_depth`.
- **contamination 0.3.0**: a haplocheck YES is kept only if its minor haplogroup is
  supported by >= 5 minor heteroplasmies (`--min-support`) and is not rCRS (H2a2a1*);
  otherwise NO, flagged `unsupported`. Levels < 5% at haplocheck coverage < 2000x
  (`--lowconf-depth`) are flagged `low_confidence` (Weissensteiner et al. 2021: >600x
  for 1%; unreliable levels at low mtCN). New metrics: `haplocheck_status`,
  `haplocheck_level`, `minor_haplogroup`, `minor_heteroplasmies`.
- **heteroplasmy 0.3.0, NUMT filter**: a heteroplasmy whose allele is found in NUMT
  sequences (`presets/numt_alleles.tsv`, 379 rCRS alleles from MitoHPC's NUMT.vcf.gz, MIT
  licence; `tools/make_numt.sh`) is left out when HF < k x nuclear depth / chrM depth
  (the HF NUMT reads can give; `--numt-k`, default 1.5). Nuclear depth: `--nuclear-depth`
  or samplesheet column `nuclear_depth`, else the copy_number B2M depth, else a fixed
  ceiling (`--numt-max-hf`, 0.10). New metrics `numt_excluded_count`, `numt_hf_threshold`.
  In low-mtCN samples this removes most NUMT-derived calls, while high-mtCN samples are
  unchanged. Blocklist adds m.567 (567A>AC, 568-573 poly-C, recurrent at ~20% HF).
- **`mitrix cohort --max-heteroplasmy N`** (default 20): `high_heteroplasmy` in
  `qc_flag` for samples with more heteroplasmies (mixtures, NUMT-rich low-mtCN samples).
- Note on MitoHPC input: MitoHPC passes Mutect2's GT 0/1 for every call to haplocheck,
  against haplocheck's input spec (homoplasmies 1/1), so its minor haplogroup collapses
  to rCRS and some contamination is missed. mitRix keeps GT 1/1 at HF >= 0.9 (0.95 and
  0.99 tried: more NUMT-driven calls at low mtCN).

## 0.6.0 — 2026-11-06

### New

- **`pathogenic` module 0.1.0**: known pathogenic mtDNA variants per sample, homoplasmic
  or heteroplasmic (HF >= 0.03, >= 10 ALT reads, caller filters excluded as for
  heteroplasmy). Default list `presets/mt_pathogenic.tsv` (162 variants): MITOMAP
  confirmed (all `Cfrm`, incl. `Cfrm [VUS*]`; exports 2026-03-20/04-16) and ClinVar
  Pathogenic/Likely pathogenic with review criteria (2026-11-01; "no assertion
  criteria" left out), normalised against rCRS; built by `tools/make_pathogenic.sh`
  (MITOMAP allele parsing from the lab's `add_VAR_MitoMAP.R`). User lists with
  `--opt pathogenic='--list NAME=FILE'` (pos, ref, alt; normalised the same way).
  Metrics per list: `_count`, `_variants`, `_max_hf`, `_low_depth_sites`; per-variant
  table `<sample>.pathogenic.variants.tsv`. Cohort columns `pathogenic_n`, `pathogenic`,
  `pathogenic_max_hf` and `NAME_n`, `NAME`, `NAME_max_hf` per user list.
- `tests/test_pathogenic.sh`.
- **`mitrix carriers`**: carriers of user-listed nuclear variants (disease genes) in the
  multi-sample nuclear VCF used by `mitrix prs`: `--list NAME=FILE` (repeatable; chr,
  pos, ref, alt, gene, optional label, inheritance, classification), multi-allelic
  records split, optional left-alignment (`--reference-fasta`), 'chr' naming handled,
  missing genotypes kept as NA. Outputs `DIR/carriers/` (long table, carried variants,
  match summary); `mitrix cohort --carriers` adds `NAME_n`, `NAME_<GENE>`, `NAME`, and
  `NAME_<GENE>_<SEVERITY>` for genes whose list has more than one `severity` class.
  Entries whose `classification` contains "modifier" (APOE e4/e2) count only in their
  gene columns, not in `NAME_n` or `NAME`.
  1000 Genomes (3,202 samples) with the PD example: GBA1 p.Glu365Lys (E326K) in 3.6% of
  EUR, LRRK2 p.Gly2385Arg in 4.3% of EAS. `tests/test_carriers.sh`.
- **Example nuclear PD variant list** for `mitrix carriers`:
  `presets/examples/pd_nuclear_variants_GRCh{38,37}.tsv` (468 variants; built by
  `tools/make_pd_nuclear_list.sh` from ClinVar 2026-11-01): SNCA, LRRK2, VPS35 (AD),
  PRKN, PINK1, PARK7 (AR), GBA1 (risk); P/LP with review criteria plus moderate-effect
  risk variants (GBA1 E326K, T369M; LRRK2 G2385R, R1628P, N2081D). GBA1 `severity` by
  Gaucher severity (Gan-Or 2015; Parlar 2023): severe (L444P, D409H, V394L, R120W,
  84GG, IVS2+1G>A, ...), mild (N370S, R496H), risk, unclassified (330 other Gaucher
  P/LP). Columns chr, pos, ref, alt, gene, label (HGVS protein; N370S = p.Asn409Ser),
  inheritance, classification, severity, rsid, source. APOE rs429358/rs7412 (e4/e2;
  dementia modifier, not a PD gene) give `PD_APOE_e4`/`PD_APOE_e2` allele counts.
  1000 Genomes APOE e4 allele frequency: EUR 0.15, AFR 0.26, EAS 0.085. Common GWAS
  variants are not included.
  Structural variants (PRKN exon deletions) are not included.

## 0.5.1 — 2026-11-06

### Fixes (found on 44 1000 Genomes CRAMs with MitoHPC calls)

- `mt_variants`: records whose REF does not match rCRS are dropped with a warning
  (count in params as `ref_mismatch_dropped`) instead of failing the sample. MitoHPC
  writes some poly-C indels at chrM:16180-16181 (e.g. `AC>A`) against the sample's own
  consensus sequence; 6 of 44 samples had 1-3 such records.
- `prs`: cohorts of fewer than 50 samples failed in plink2 `--score` (no allele
  frequencies to impute missing genotypes); frequencies are now computed from the cohort
  and passed with `--read-freq`. Scores are unchanged for 50 or more samples.
- `cohort`: PRS flags are collapsed into one `qc_flag` entry (`PRS:non_EUR`) instead of
  one per PRS column.

## 0.5.0 — 2026-11-06

### New

- **`mitrix prs`**: polygenic scores `PRS_CN` (blood-adjusted mtCN GWAS, GCST90840190),
  `PRS_CN_raw` (raw mtCN, GCST90840191), `PRS_het`, `PRS_7S`, `PRS_primer` from a
  multi-sample nuclear VCF (GRCh38 or GRCh37) with plink2: genotype QC, ancestry by
  projection onto 1000 Genomes PCs (non-EUR samples flagged), HWE in EUR, variant
  matching with strand flips, z-scores within EUR (raw sums as `PRS_*_sum`),
  lead-variant scores (`PRS_*_lead`, always), validation against the measured readouts
  (`--cohort`, with `--metadata`/`--covariates`). 1000 Genomes
  (3,202 samples, GRCh38): 98.8% of weight SNPs found, 633/633 EUR and 3,186/3,202
  overall assigned to their super-population, about 1 minute on 8 threads.
- **Weights** in `presets/prs/` (31 MB), built once by `tools/make_prs_weights.sh`:
  PRS-CS-auto (1000G EUR LD panel) on Gupta et al. GWAS Catalog summary statistics,
  LDSC SNP-h² and expected R² per score (`scores.tsv`). `PRS_CN_raw` (h² 0.061,
  expected R² ~1.8%, 304 lead SNPs) and `PRS_CN` (h² 0.031, ~0.5%, 139) have clear
  polygenic heritability; `PRS_het`, `PRS_7S` and `PRS_primer` are
  weak (h² < 0.01; 4, 6 and 1 independent lead SNPs).
- **`mitrix cohort --prs`** adds the `PRS_*` and `ancestry` columns.
- **`d_loop_7s` 0.3.0**: `7S_depth_median`, `flap_depth_median`, `primer_depth_median`,
  `flap_7S_ratio` (flap / 7S) and `flap_primer_ratio` (flap / primer), as in Gupta et al.
  2023 (regions: 7S DNA first third 16542–191, flap 191–298, RNA primer 346–407);
  cohort columns `flap_7S` and `flap_primer`. Existing D-loop values are unchanged.
- `tests/test_prs.sh`: exact score sums on synthetic data (GRCh38 and GRCh37, strand
  flips, swapped alleles, missing genotypes, overlap flag); `tests/test_prs_1000g.sh`:
  release test on 1000 Genomes (GRCh38 vs GRCh37 sums r > 0.997, missing genotypes,
  subsets, simulated validation).

## 0.4.0 — 2026-11-05

Any caller's mtDNA variants are accepted; haplogroup and contamination are computed in
mitRix; the heteroplasmy module gains HF bins, a strand-resolved spectrum and functional
counts; new `mitrix check`, `cohort` and `report`. Damage, d_loop_7s and copy_number
values are unchanged from 0.3.0.

### Inputs

- **`--mt-variants FILE`** replaces `--mt-vcf`: a VCF from MitoHPC/Mutect2, mity or any
  caller with `FORMAT/AF` or `VAF`, or mtDNA-Server 2's `variants.annotated.txt`. The
  caller is detected (`--caller` overrides). **`--variants-sample`** replaces
  `--mitohpc-name` and defaults to `--sample`, then the BAM/CRAM name; a unique
  `NAME.<suffix>` also matches. `--mt-vcf` and `--mitohpc-name` still work (deprecated).
- **Removed: `--haplogroup-tab`, `--haplocheck-tab`** (and samplesheet columns). The
  haplogroup and contamination modules now run HaploGrep3 3.3.2 (PhyloTree 17 Forensic
  Update, `phylotree-fu-rcrs@1.3`) and haplocheck 1.3.3 on a clean per-sample VCF
  (`core/<sample>.mt.vcf.gz`). Install both once with `tools/install.sh` (Java ≥ 11).
- Samplesheet columns: `sample bam reference_fasta mt_variants variants_sample caller`.

### Changes that alter values

- **MSS** uses MLC scores from Lake et al. 2024 (Supplementary Dataset 7, rounded to 5
  decimals) for every caller, instead of MitoHPC's `MLC_score` INFO field (an earlier
  release of the scores). Values differ by < 0.007 per variant (HG00096 MSS 0.2226447 →
  0.22264). Consequences are identical.
- **Heteroplasmy** also needs ≥ 10 ALT reads (`--min-alt-reads`). The filters excluded by
  default depend on the caller: `FAIL` for mity, `strand_bias` otherwise.
- **Contamination**: haplocheck now runs on the clean VCF (strand-bias calls excluded),
  so `haplocheck_distance` and `haplocheck_coverage` can differ from MitoHPC's table
  (HG00096: distance 9 → 0; status NO in both). `mitrix cohort` flags `contaminated` only
  when the level is above `--max-contamination` (default 0.02, as in Gupta et al.,
  Nature 2026), not on any YES.
- Core variant table: adds `alt_reads`; `numt` and `mlc_score` are removed (MLC is
  annotated by the heteroplasmy module).

### New

- **`heteroplasmy` 0.2.0** metrics: `heteroplasmy_low_count` / `heteroplasmy_high_count`
  (HF ≤ / > `--hf-split` 0.2), `snv_count`, `indel_count`, `H_CtoT_count`, `H_AtoG_count`,
  `L_AtoG_count`, `L_CtoT_count` (strand: rCRS REF C/A = light, G/T = heavy, written as
  the complement), `transversion_count`, `age_accumulating_count` and
  `age_accumulating_low_count` (H C>T, H A>G, L A>G outside chrM:16172–210; Gupta et al.,
  Nature 2026), `MSS_low`, `MSS_high`, `nonsynonymous_count`, `MSS_nonsynonymous`,
  `mt_median_depth`, `hf_detection_limit` (= max(hf_min, min_alt_reads / median depth);
  low-bin metrics flagged `low_depth` when above hf_min). The per-variant table adds
  `alt_reads strand substitution hf_bin age_accumulating gene region consequence`.
- **`presets/mt_annotation.tsv.gz`**: gene, region, consequence and MLC score for all
  49,704 rCRS SNVs (Lake et al. 2024; MITOMAP loci), built by `tools/make_annotation.sh`.
- **`haplogroup` 0.2.0** adds `haplogroup_quality`; **`contamination` 0.2.0** reports
  `ND` when haplocheck cannot decide.
- **`mitrix check`**, **`mitrix cohort`** (`cohort.tsv`, `cohort_long.tsv`; new columns
  `heteroplasmy_low`, `age_accumulating`) and **`mitrix report`** (HTML with QC, batch
  adjustment, covariate associations, haplogroups, heteroplasmy spectrum and crude
  dN/dS, mito-disease subtype analysis). A `column:low_depth` flag sets only that value
  to `NA` in the report. `<sample>.tsv` lists the shared ND1 depth once.
- `tests/test_heteroplasmy.sh`: strand mapping of all 12 SNV types, origin boundary,
  HF split, ALT-read minimum, filters and blocklist.

HG00096 (MitoHPC, mity and mtDNA-Server 2 calls): 3 heteroplasmies (mtDNA-Server 2: 4,
adding 471 T>G at 3.2%), 2 age-accumulating, MSS 0.22264, H16a1, not contaminated.

## 0.3.0 — 2026-10-23

MitoHPC outputs are now read: heteroplasmy, haplogroup and contamination.
Existing values (damage, d_loop_7s, copy_number) are unchanged from 0.2.0.

### New

- **MitoHPC inputs** (all optional): `--mt-vcf`, `--haplogroup-tab`,
  `--haplocheck-tab`, and `--mitohpc-name` (MitoHPC's sample name; defaults to the
  BAM/CRAM file name without extension, as MitoHPC names samples).
- **Samplesheet**: `mitrix run --samplesheet SS` runs one sample per row (columns
  `sample`, `bam`, optional `reference_fasta`, `mt_vcf`, `mitohpc_name`,
  `haplogroup_tab`, `haplocheck_tab`).
- **`core/mt_variants.sh`**: MitoHPC VCF adapter. Repairs MitoHPC's concatenated
  VCF, validates sample name, contig and REF against rCRS, left-aligns indels, and
  writes the standard variant table (`hf` from `FORMAT/AF`; ALT == REF rows dropped).
- **`heteroplasmy`**: `heteroplasmy_count` (3% ≤ HF ≤ 95%, no `strand_bias`, not
  blocklisted), `MSS` (sum of MLC scores), `MSS_missing_count` (indels, which have no
  MLC score), and a per-variant table with HF. Thresholds, filters and blocklist are
  options (`--opt heteroplasmy='--hf-min 0.01'`).
- **`haplogroup`**: `haplogroup` and `macro_haplogroup` from MitoHPC's table.
- **`contamination`**: haplocheck status, level, distance and coverage from
  MitoHPC's table.
- `mitrix run --opt MODULE='ARGS'` passes options to a module.

Modules whose inputs are not given write `SKIPPED` rows and need none of their tools
(no R without `--mt-vcf`). HG00096: 3 heteroplasmies, MSS 0.2226, H16a1, not
contaminated.

## 0.2.0 — 2026-10-23

mitRix is now modular: `bin/mitrix run` runs core steps and biomarker modules
(`damage`, `d_loop_7s`, `copy_number`) and writes one long-format table per sample.
See `docs/module-contract.md`. `bin/mitrix-legacy` keeps the v0.1.0 command line
and 8-column output.

### Changes that alter values

- **Duplicates are kept by default** (`--flag 2820`, was `1796`). At chrM depth,
  independent fragments often share start and end by chance and are marked
  duplicate: in HG00096, 16.2% of chrM reads vs ~10% of nuclear reads. Excluding
  them lowered chrM depth more than B2M depth, biasing mtDNA-CN down, and more so
  in samples with more mtDNA. Use `--flag 1796` for the old behaviour.
- **Full-precision output.** Values in the long tables are not rounded (`%.10g`);
  v0.1.0 rounded depths and ratios to 2 decimals. `mitrix-legacy` still rounds.

Effect on 1000 Genomes HG00096 (30x WGS, GRCh38), v0.1.0 ratios computed from its
unrounded depths:

| ratio | stat | v0.1.0 | v0.2.0 | change |
|---|---|---|---|---|
| mtDNA-CN (ND1/B2M) | mean | 408.44 | 435.29 | +6.6% |
| ND4/ND1 | mean | 1.0255 | 1.0237 | −0.2% |
| ND5/ND1 | mean | 1.0197 | 1.0160 | −0.4% |
| DLOOP/ND1 | mean | 0.9749 | 0.9944 | +2.0% |
| mtDNA-CN (ND1/B2M) | median | 418.06 | 453.51 | +8.5% |
| ND4/ND1 | median | 1.0232 | 1.0222 | −0.1% |
| ND5/ND1 | median | 1.0140 | 1.0104 | −0.4% |
| DLOOP/ND1 | median | 0.9709 | 0.9943 | +2.4% |

All of the change comes from keeping duplicates: with `--flag 1796` v0.2.0
reproduces v0.1.0 exactly on this sample and on the test slices.

### Fixes (no change on HG00096 at the v0.1.0 flag)

- chrMC realignment keeps unmapped mates (flag 4) placed on chrM, which are often
  seam-crossing reads. HG00096's unmapped chrM reads are all marked duplicate, so
  this shows only with duplicates kept (+1–2 depth at 293 positions on the test slice).
  On the GRCh37 test slice (duplicates re-marked by samtools markdup) it adds +1 depth
  at 9 chrMC positions at `--flag 1796`; ratios unchanged at 2 decimals.
- chrMC realignment no longer depends on `--threads` (`bwa mem -K` fixed).
- DLOOP/ND1 no longer crashes when no chrM read has its mate off chrM (e.g. chrM-only BAMs).
- Fold-back from chrMC uses the detected contig name and length instead of `chrMC`/16569.
- Nothing is written to the working directory; `mitrix run` requires `--outdir`.

### New

- Per-sample `<sample>.tsv`, `<sample>.status.tsv`, `params.tsv` provenance, step logs.

## 0.1.0

First release: single `mitRix` script; ND4/ND1, ND5/ND1, DLOOP/ND1 (chrMC
realignment) and mtDNA-CN (ND1/B2M) from ThermoFisher 4-plex amplicon coordinates.
