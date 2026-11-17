# Module contract

Every core step and biomarker module follows this contract, so that `bin/mitrix`,
Nextflow and WDL can run any of them the same way, and `aggregate/` can merge their
outputs without knowing which module produced them.

Status: draft for ROADMAP p0-4. Module metadata (`module.conf`) is specified in
[module-conf.md](module-conf.md).

## 1. Layout

```
core/<step>.sh                     ref_context, mt_coverage, mt_circular
biomarkers/<name>/module.conf      metadata (see module-conf.md)
biomarkers/<name>/run.sh           entry point (or run.R)
biomarkers/<name>/README.md        method, metrics, options
lib/*.sh                           shared helpers, sourced by modules
```

Outputs:

```
<outdir>/<sample>/<sample>.tsv           all biomarkers for this sample (the main per-sample output)
<outdir>/<sample>/<sample>.status.tsv    one row per core step / module: exit code, times, version
<outdir>/<sample>/logs/<step>.log        stderr of each core step / module
<outdir>/<sample>/core/                  written by core steps, read by modules
<outdir>/<sample>/<module>/              intermediate: written by one module only
<outdir>/cohort.tsv                      from `mitrix cohort` (section 6)
<outdir>/report.html                     from `mitrix report`, when metadata is given
```

A module writes only inside its own directory, never to the working directory.
`bin/mitrix` builds `<sample>.tsv` by concatenating the module tables (section 3).

## 2. Inputs

`bin/mitrix` calls every module with the same flags:

| flag | required | meaning |
|---|---|---|
| `--sample NAME` | yes | sample ID; used in file names and the `sample` column |
| `--bam FILE` | yes | indexed BAM or CRAM |
| `--reference-fasta FA` | CRAM only | reference used to decode the CRAM |
| `--outdir DIR` | yes | top-level output directory; the module reads `<outdir>/<sample>/core/` and writes to `<outdir>/<sample>/<module>/` (creating it) |
| `--threads N` | no (default 1) | |
| `--mt-variants FILE` | no | mtDNA variant calls (samplesheet `mt_variants`): VCF from MitoHPC/Mutect2, mity or another caller with `FORMAT/AF` or `VAF`, or mtDNA-Server 2 `variants.annotated.txt`; read by `core/mt_variants.sh` |
| `--variants-sample NAME` | no (default: `--sample`, then the BAM/CRAM file name without extension) | sample name inside that file; a unique `NAME.<suffix>` also matches (MitoHPC `HG00096.final`, mtDNA-Server 2 `HG00096.chrM.bam`) |
| `--caller NAME` | no (default: detected) | `mitohpc`, `mutect2`, `mity`, `mtdnaserver2` or `vcf` |

`--mt-vcf` and `--mitohpc-name` are deprecated aliases of `--mt-variants` and
`--variants-sample` (0.4.0). `--haplogroup-tab` and `--haplocheck-tab` were removed
in 0.4.0: haplogroup and contamination are computed from the variant calls.

Module-specific options come after these and are documented in the module README
and `--help`. Unknown flags are an error.

### Core outputs a module may read

| file in `<outdir>/<sample>/core/` | written by | content |
|---|---|---|
| `<sample>.ref_context.tsv` | `core/ref_context.sh` | key/value: `mt_contig`, `mt_length`, `chr15_contig` (empty if absent), `chr15_length`, `build` (`hg19`/`hg38`/`unknown`) |
| `<sample>.contigs.tsv` | `core/ref_context.sh` | all `@SQ` contigs: name, length |
| `<sample>.chrM.per-base.bed.gz` | `core/mt_coverage.sh` | one line per position, 0-based BED, column 4 = depth |
| `<sample>.chrM.corrected.per-base.bed.gz` | `core/mt_circular.sh` | same format, after chrMC realignment and fold-back |
| `<sample>.chrMC.per-base.bed.gz` | `core/mt_circular.sh` | chrMC (17269 bp) coverage before fold-back |
| `<sample>.variants.tsv` | `core/mt_variants.sh` | standard variant table: `sample pos ref alt hf depth alt_reads filter caller type`; absent without `--mt-variants` |
| `<sample>.variants.params.tsv` | `core/mt_variants.sh` | adapter provenance: input path, format, caller, sample, HF field, records kept/dropped, `exclude_filters` (caller default) |
| `<sample>.mt.vcf.gz` (+`.tbi`) | `core/mt_variants.sh` | clean single-sample VCF for HaploGrep3 and haplocheck: `chrM`, rCRS, `GT` 1/1 if HF ≥ 0.9 else 0/1, `FORMAT/AF`, `DP`; records with an excluded filter left out |

`core/mt_variants.sh` detects the caller (`##mityCommandline` → mity; `MLC_score`
INFO → MitoHPC; `##source=Mutect2` → mutect2; mtDNA-Server 2's table header), picks
the sample, repairs MitoHPC's concatenated VCF (keeps the last header block), checks
the contig (one mtDNA contig of 16569 bp) and REF against rCRS, left-aligns and splits
multi-allelic records with `bcftools norm`, drops ALT == REF rows, and warns when every
record is `PASS` (pre-filtered input). `hf` is `FORMAT/AF` (MitoHPC, Mutect2, other
VCFs), `FORMAT/VAF` (mity) or `VariantLevel` (mtDNA-Server 2); `alt_reads` =
round(`hf` × `depth`). Filters excluded by default: `FAIL` for mity, `strand_bias`
otherwise.

Core steps run once per sample, before any module. A module lists the core steps it
needs in `module.conf` (`requires_core`).

### Samplesheet

`mitrix run --samplesheet SS --outdir DIR` runs one sample per row. TSV with a header;
required columns `sample`, `bam`; optional `reference_fasta`, `mt_variants`,
`variants_sample`, `caller` (empty or `NA` = not given). Relative paths are
relative to the samplesheet. `mitrix check` / `cohort` (P5) read the same file.

## 3. Outputs

Each module writes two files to `<outdir>/<sample>/<module>/`:

### `<sample>.<module>.tsv` — results

Tab-separated, with a header line and these 7 columns in this order:

| column | type | example | rules |
|---|---|---|---|
| `sample` | text | `HG00096` | as given by `--sample` |
| `biomarker` | text | `damage` | the module name |
| `metric` | text | `ND4_ND1_ratio_mean` | `[A-Za-z0-9_]+`, unique within the file |
| `value` | number, text or `NA` | `1.025496` | full precision (`%.10g`); no rounding. Text only for unit `label` (e.g. haplogroup `H16a1`) |
| `unit` | text | `ratio` | from the unit list below |
| `qc_flag` | text | `PASS` | from the flag list below; several joined with `;` |
| `version` | text | `0.2.0` | `version` from the module's `module.conf` |

- One row per metric. Rows for metrics that could not be computed are kept, with
  `value` = `NA` and a non-`PASS` flag, so every sample has the same rows.
- Metric names put the statistic last: `<region>_depth_mean`, `<target>_<reference>_ratio_median`.

**Units:** `x` (read depth), `ratio`, `copies` (per diploid nuclear genome),
`fraction` (0–1), `count`, `score`, `label` (value is text, e.g. a haplogroup or `YES`/`NO`).

**QC flags:**

| flag | meaning |
|---|---|
| `PASS` | computed, no issues |
| `LOW_DEPTH` | below the module's depth threshold; value still reported |
| `ZERO_DEPTH` | a denominator was 0; `value` = `NA` |
| `SKIPPED` | input not available (e.g. no chr15 contig, no `mt_variants`); `value` = `NA` |
| `FAIL` | computation failed; `value` = `NA` |

Sample-level flags such as contamination (p3-7) are added later by `aggregate/` (p5-3),
not by modules.

### `<sample>.<module>.params.tsv` — provenance

Key/value, tab-separated, no header. Records everything needed to reproduce the
numbers: resolved regions, options, tool versions. For example:

```
stat            mean
flag            2820
region.ND4      chrM:12087-12170
region.ND1      chrM:3485-3553
samtools        1.23.1
mosdepth        0.3.14
```

### `<sample>.tsv` — per-sample result (built by `bin/mitrix`)

After all requested modules for a sample finish, `bin/mitrix` concatenates their
`<sample>.<module>.tsv` files (one header, same 7 columns) into
`<outdir>/<sample>/<sample>.tsv`, and writes `<sample>.status.tsv`:

```
step           status  exit_code  start                 end                   version
ref_context    OK      0          2026-10-23T12:00:01   2026-10-23T12:00:02   0.1.0
damage         OK      0          2026-10-23T12:00:02   2026-10-23T12:00:09   0.1.0
copy_number    FAILED  1          2026-10-23T12:00:02   2026-10-23T12:00:04   0.1.0
```

`<sample>.tsv` is only written when every step is `OK`; a failed module's rows are
never mixed into it.
A metric reported by several modules (ND1 depth: `damage`, `d_loop_7s`,
`copy_number`) appears once, from the first module; the module tables keep their copies.

## 4. Behaviour

- Exit `0` when the results file is written, including when every row is `SKIPPED`.
  Non-zero only when the module could not write a valid results file.
- Log to stderr, one line per message, prefixed `<module>: `. Nothing else on stdout.
- Temporary files go in a directory from `lib/common.sh` (`make_tmpdir`), removed on exit.
- Use the region helpers in `lib/regions.sh` and the statistics in `lib/stats.sh`
  instead of reimplementing them.
- `run.sh --help` prints usage, module-specific options and the metric list.

## 5. v0.1.0 mapping

For p1-12, `bin/mitrix --legacy-tsv` rebuilds the v0.1.0 8-column output
(`sample ratio target reference stat target_cov ref_cov ratio_value`) from the long
tables and `params.tsv`, rounding to 2 decimals as v0.1.0 did. `SKIPPED` rows are
dropped (v0.1.0 omitted them). `baseline_v0.1.0/compare.sh` then runs unchanged.

| v0.1.0 `ratio` | module | metrics (for `--mean`; `_median` for `--median`) |
|---|---|---|
| `ND4/ND1` | `damage` | `ND4_depth_mean`, `ND1_depth_mean`, `ND4_ND1_ratio_mean` |
| `ND5/ND1` | `damage` | `ND5_depth_mean`, `ND1_depth_mean`, `ND5_ND1_ratio_mean` |
| `DLOOP/ND1` | `d_loop_7s` | `DLOOP_depth_mean`, `ND1_depth_mean`, `DLOOP_ND1_ratio_mean` |
| `mtDNA-CN` | `copy_number` | `ND1_depth_mean`, `B2M_depth_mean`, `ND1_B2M_ratio_mean` |

`target` and `reference` come from `region.*` in `params.tsv`. In TaqMan mode,
`copy_number` reports `ND1_B2M_ratio` (unit `ratio`), like v0.1.0. A WGS mode adding
`mtDNA_CN` (unit `copies`) is planned for mitRix Plus.

Example `HG00096.damage.tsv`:

```
sample   biomarker  metric              value     unit   qc_flag  version
HG00096  damage     ND4_depth_mean      13667.33  x      PASS     0.1.0
HG00096  damage     ND1_depth_mean      13327.54  x      PASS     0.1.0
HG00096  damage     ND4_ND1_ratio_mean  1.025496  ratio  PASS     0.1.0
HG00096  damage     ND5_depth_mean      13589.78  x      PASS     0.1.0
HG00096  damage     ND5_ND1_ratio_mean  1.019677  ratio  PASS     0.1.0
```

## 6. Cohort commands

Run after all samples have been processed with `mitrix run`.

### `mitrix check --outdir DIR --samplesheet SS`

Reads every `<sample>.status.tsv` and reports one line per sample in the samplesheet:

```
sample    status      detail
HG00096   OK
HG00097   FAILED      copy_number (exit 1)
HG00099   MISSING     no output directory
```

`INCOMPLETE` means some steps have not finished (e.g. a job still running or killed).
Exits `0` only if every sample is `OK`.

### `mitrix cohort --outdir DIR --samplesheet SS [-o cohort.tsv] [--allow-incomplete] [--max-contamination 0.02] [--cn-clusters auto|off|K]`

Runs `check` first and stops if any sample is not `OK`, unless `--allow-incomplete`
is given (failed samples are then kept as rows of `NA`). Writes:

- `cohort.tsv` — one row per sample, one column per headline readout:

  ```
  sample  mtDNA_CN  deletion  7S     heteroplasmy  heteroplasmy_low  age_accumulating  MSS     haplogroup  macro_haplogroup  contamination  contamination_level  qc_flag
  SIM01   41.47     1.0196    1.032  3             0                 2                 0.223   H1c         H                 NO             0                    PASS
  SIM07   26.55     1.0403    1.005  79            67                48                10.850  V2          V                 YES            0.050                contaminated
  ```

  Each module declares which of its metrics become cohort columns, and under what
  name, with the `cohort` key in `module.conf`; known columns come in the order above.
  A `_mean` metric falls back to `_median` for runs with `--stat median`. `qc_flag` is
  `PASS`, or lists columns whose source row was not `PASS` (e.g. `mtDNA_CN:LOW_DEPTH`),
  `contaminated` when haplocheck's contamination level is above `--max-contamination`
  (default 0.02), and `FAILED`/`MISSING`/
  `INCOMPLETE` for samples kept with `--allow-incomplete` (their values are `NA`).
- `cohort_long.tsv` — every row of every `<sample>.tsv`, for analyses that need more
  than the headline readouts.
- `cn_cluster`, `cn_cluster_prob` (`--cn-clusters`, default `auto`; needs the R package
  mclust): groups of mtDNA copy number in the cohort, e.g. mixed DNA sources (whole blood
  vs cell line / PBMC) or batches. Gaussian mixtures (mclust, univariate equal- and
  unequal-variance models) are fitted to log(`mtDNA_CN`) of samples whose `mtDNA_CN` is not
  flagged. `auto` picks K = 1-4 by ICL, which penalises overlapping components so that one
  skewed group is not split, and accepts K only if every cluster holds at least
  max(10, 5% of samples); a number fixes K. Clusters are numbered by increasing mtDNA-CN
  (1 = lowest) and are not labelled as tissues; with one group every sample is `1`.
  `cn_cluster_prob` is the posterior probability of the assigned cluster; below 0.9
  `qc_flag` gets `cn_cluster:uncertain` (a soft flag: the report only blanks `cn_cluster`).
  `cohort_cn_clusters.tsv` gives per cluster n, median and geometric-mean CN, and the
  chosen K and model; `cohort_cn_clusters_fit.tsv` the ICL and BIC of every K and model.
  Skipped with a warning when mclust is not installed; `--cn-clusters off` disables it.
- With `--carriers DIR/carriers/carriers.tsv` (`mitrix carriers --vcf NUCLEAR.vcf.gz
  --list NAME=FILE`): per user list, `NAME_n` (listed variants carried), `NAME_<GENE>`
  (listed ALT alleles carried per gene) and `NAME` (e.g. `LRRK2:p.Gly2019Ser(het)` or
  `none`); missing genotypes give `NAME:missing_genotypes` in `qc_flag`. Lists are the
  user's own (chr, pos, ref, alt, gene, optional label, inheritance, classification);
  an example for PD from ClinVar is `presets/examples/pd_nuclear_variants_GRCh38.tsv`.
- `pathogenic_n`, `pathogenic`, `pathogenic_max_hf` (`pathogenic` 0.1.0): number of
  known pathogenic mtDNA variants detected (MITOMAP Cfrm + ClinVar P/LP with review
  criteria, `presets/mt_pathogenic.tsv`), the variants with their HF
  (`m.3243A>G(0.18);m.11778G>A(1)` or `none`), and the highest HF. Each user list given
  with `--opt pathogenic='--list NAME=FILE'` adds `NAME_n`, `NAME`, `NAME_max_hf`.
  `pathogenic_n` is flagged `low_depth` when a listed position has depth below
  min-alt-reads / hf-min (an absent call there is not informative).
- `flap_7S` and `flap_primer` (`d_loop_7s` 0.3.0) are median-depth ratios over the
  regions of Gupta et al. 2023: 7S DNA first third 16542–191, 7S DNA flap 191–298,
  7S RNA primer 346–407 (`flap_7S` = flap / 7S, `flap_primer` = flap / primer, the
  direction of their flap-given-7S/primer traits).

### `mitrix prs --vcf NUCLEAR.vcf.gz --outdir DIR [--build 38|37] [--lead] [--cohort cohort.tsv]`

Polygenic scores from one multi-sample nuclear VCF (sample IDs = mitRix sample names).
Needs plink2 and Rscript only; the weights ship in `presets/prs/` (built once by
`tools/make_prs_weights.sh`: PRS-CS-auto with the 1000G EUR LD panel on Gupta et al.
GWAS; LDSC SNP-h² and expected R² in `presets/prs/scores.tsv`).

| score | GWAS (GWAS Catalog) | validation readout |
|---|---|---|
| `PRS_CN` | mtDNA copy number, adjusted (not the "raw" GCST90840191), 2026 (GCST90840190) | `mtDNA_CN` |
| `PRS_CN_raw` | mtDNA copy number, raw, 2026 (GCST90840191) | `mtDNA_CN` |
| `PRS_het` | age-related mtSNV accumulation, 2026 (GCST90840194) | `age_accumulating` |
| `PRS_7S` | coverage discrepancy 7S DNA vs flap, 2023 (GCST90268499) | `flap_7S` |
| `PRS_primer` | coverage discrepancy RNA primer vs flap, 2023 (GCST90268500) | `flap_primer` |

Steps: autosomal biallelic SNVs; variant and sample call rate ≥ 95% (`--geno`, `--mind`);
ancestry by projection onto 1000 Genomes PCs (69k pruned SNPs) and 20-nearest-neighbour
super-population (≥ 80% agreement, else `unassigned`); HWE p ≥ 1e-6 in EUR samples (if
≥ 50); MAF ≥ 0.01; weights matched by chr:pos of the chosen build and alleles (strand
flips handled); `plink2 --score` sums (missing genotypes mean-imputed); z-scores within
EUR samples. `--lead` adds `PRS_<score>_lead` (GWAS betas of independent genome-wide
significant SNPs).

Writes `DIR/prs/`: `prs.tsv` (long format; metrics `PRS_<score>` (z), `PRS_<score>_sum` (raw score sum),
`ancestry`, `ancestry_fraction`, `PC1`–`PC10`; `qc_flag` `non_EUR` and/or `low_overlap`
when fewer than `--min-overlap` (0.8) of a score's SNPs are found), `match.tsv`,
`ancestry.tsv`, `qc.tsv`, `params.tsv`, `logs/`, and with `--cohort` `validation.tsv`
(each PRS vs its readout in EUR samples: Spearman r and incremental R² over PC1–10,
plus age/sex with `--metadata`). `mitrix cohort --prs DIR/prs/prs.tsv` adds the `PRS_*`
and `ancestry` columns to `cohort.tsv`.

### `mitrix report --cohort cohort.tsv --metadata META [-o report.html]`

Only runs when the user supplies metadata: a TSV with a `sample` column plus any
covariates (e.g. age, sex, batch, disease). Options: `--subtype-col` (default
`disease`), `--covariates` (default `age,sex,batch`, those present), `--keep-flagged`,
`--outdir` (run directory with the per-variant heteroplasmy tables; default: the
directory of `--cohort`). A `column:low_depth` flag sets only that column to `NA` for
the sample; any other flag excludes the sample.
Renders `aggregate/report.Rmd` to one self-contained HTML file (R: rmarkdown, knitr,
ggplot2; optional glmnet for the supervised model and uwot for UMAP; needs pandoc),
and writes `<report>_adjusted.tsv` with raw, batch-adjusted (`_adj`) and z-scored
(`_z`) readouts. Sections:

- QC summary: samples in cohort and metadata, flagged samples (excluded unless `--keep-flagged`)
- distribution of each readout, by subtype and by batch; summary statistics
- Spearman correlation between readouts
- batch adjustment (linear model; overall mean kept) and z-scores (p5-2)
- association of each readout with each covariate, adjusted for the others (BH FDR)
- readouts by macro-haplogroup (Kruskal–Wallis)
- heteroplasmy spectrum and function, from the per-variant tables: strand-resolved
  substitution classes by HF bin and subtype, heteroplasmies by region, crude dN/dS by
  HF bin (MT-ATP8/ATP6 overlap excluded) and nonsynonymous heteroplasmies by gene
- **mito-disease subtypes** (ROADMAP p5-6), on batch-adjusted z-scored readouts:
  - *Unsupervised*: PCA (+ loadings), UMAP; k-means and Ward clustering with k =
    number of subtypes, agreement as adjusted Rand index
  - *Supervised*: elastic-net (alpha 0.5) multinomial/binomial model, nested
    cross-validation stratified by subtype (outer ≤5 folds, inner CV for lambda);
    out-of-fold AUC per subtype, accuracy, confusion matrix, coefficients and
    selection frequency across folds. Skipped if any subtype has < 4 samples after QC.

Samples in the cohort table without metadata (and vice versa) are listed in the report.
