# mitRix

**mtDNA biomarkers from whole-genome sequencing.** mitRix turns a WGS BAM/CRAM and, optionally,
mtDNA variant calls into one table of mitochondrial readouts per sample: copy number,
heteroplasmy burden, deletion and 7S DNA ratios, haplogroup, contamination, known pathogenic
variants and coverage-based QC. Cohort commands merge samples, score mtDNA-trait polygenic
scores and flag carriers of listed nuclear variants.

```
mitrix run --bam S1.cram --reference-fasta hg38.fa --mt-variants S1.vcf.gz --sample S1 --outdir results
mitrix cohort --samplesheet samples.tsv --outdir results        # -> results/cohort.tsv
```

## What it measures

| Module | Readouts | Needs |
|---|---|---|
| `copy_number` | mtDNA copy number (ND1/B2M amplicon depth ratio, as the TaqMan assay) | BAM/CRAM |
| `heteroplasmy` | heteroplasmy count (low/high HF), SNV/indel counts, strand-resolved spectrum, MLC score sum (MSS); NUMT-allele filter scaled by nuclear/chrM depth | variant calls |
| `damage` | ND4/ND1, ND5/ND1 (TaqMan amplicons), MAJOR/MINOR arc depth ratio (common-deletion span vs ND1–ND2) | BAM/CRAM |
| `d_loop_7s` | D-loop/ND1 (7S DNA) and 7S flap/7S, flap/primer ratios (Gupta et al. 2023) on circularised chrM | BAM/CRAM |
| `coverage_shape` | per-sample coverage model: GC bias, replication gradient from O_L, common-deletion step; GC-adjusted ND4/ND1 and MAJOR/MINOR | BAM/CRAM |
| `haplogroup` | mtDNA haplogroup and macro-haplogroup (HaploGrep3, PhyloTree 17 FU) | variant calls |
| `contamination` | mtDNA contamination level (haplocheck) with a support check on the minor haplogroup | variant calls |
| `pathogenic` | known pathogenic mtDNA variants (MITOMAP confirmed + ClinVar P/LP) and their HF; your own lists | variant calls |

Cohort commands:

| Command | Output |
|---|---|
| `mitrix check` | did every sample finish? |
| `mitrix cohort` | `cohort.tsv` (one row per sample, headline readouts, `qc_flag`) and `cohort_long.tsv`; flags mixed sample types or batches by clustering mtDNA copy number (`cn_cluster`) |
| `mitrix prs` | polygenic scores for mtDNA copy number, heteroplasmy and 7S DNA traits (Gupta et al. GWAS) from a nuclear multi-sample VCF, with genetic-ancestry PCs (plink2) |
| `mitrix carriers` | carriers of your listed nuclear variants (e.g. PD genes) from a nuclear VCF |

Coverage readouts (`damage`, `d_loop_7s`, `coverage_shape`) are depth-based estimates from
short reads, not deletion calls. In WGS, ND4/ND1 is sensitive to library GC bias; use the
GC-adjusted versions or `gc_slope` as a covariate when comparing groups.

## Input

- **Alignments:** BAM or CRAM (CRAM needs `--reference-fasta`), GRCh38 or GRCh37 (`chrM` or
  `MT`, detected); whole-genome or chrM-only BAMs (copy number needs chr15 for B2M).
- **mtDNA variant calls (optional):** VCF from MitoHPC / Mutect2, mity or another caller with
  FORMAT/AF or VAF, or mtDNA-Server 2 `variants.annotated.txt`. mitRix does not call variants
  itself. Without calls, heteroplasmy, haplogroup, contamination and pathogenic are skipped.
- **Variants-only mode:** `--mt-variants` without `--bam` runs the four variant-based modules.
- **Nuclear VCF (optional):** multi-sample VCF for `mitrix prs` and `mitrix carriers`.

## Install

Requirements: bash, samtools >= 1.10, mosdepth >= 0.3, bwa, bcftools, R >= 4.1 (packages
MASS; mclust for `cn_cluster`; ggplot2, rmarkdown, knitr for the report), Java >= 11 (HaploGrep3, haplocheck);
plink2 for `mitrix prs`.

```bash
git clone https://github.com/bahlolab/mitRix.git
cd mitRix
conda create -n mitrix -c conda-forge -c bioconda samtools mosdepth bwa bcftools plink2 r-base r-mass r-mclust openjdk
conda activate mitrix
tools/install.sh                 # HaploGrep3 + haplocheck jars and the PhyloTree tree (once, needs internet)
export PATH="$PWD/bin:$PATH"
mitrix --version
```

## Usage

One sample:

```bash
mitrix run --bam S1.bam --sample S1 --mt-variants S1.mitohpc.vcf.gz --outdir results
```

Many samples: a tab-separated samplesheet with columns `sample`, `bam` and optionally
`reference_fasta`, `mt_variants`, `variants_sample`, `caller`, `nuclear_depth`
(relative paths are relative to the samplesheet). Each row is independent, so rows can be
split across cluster jobs.

```bash
mitrix run    --samplesheet samples.tsv --outdir results --threads 4
mitrix check  --samplesheet samples.tsv --outdir results
mitrix cohort --samplesheet samples.tsv --outdir results
```

Options and metrics of a module, and module-specific settings:

```bash
mitrix help heteroplasmy
mitrix run ... --biomarkers copy_number,heteroplasmy --opt heteroplasmy='--hf-min 0.01'
```

Polygenic scores and nuclear variant carriers (cohort level):

```bash
mitrix prs      --vcf cohort.nuclear.vcf.gz --outdir results            # results/prs/
mitrix carriers --vcf cohort.nuclear.vcf.gz --list PD=presets/examples/pd_nuclear_variants_GRCh38.tsv --outdir results
mitrix cohort   --samplesheet samples.tsv --outdir results --prs results/prs/prs.tsv --carriers results/carriers/carriers.tsv
```

Coverage modules from precomputed mosdepth per-base files (no BAM):
`tools/from_perbase.sh --sample S1 --perbase DIR --outdir results`.

## Output

```
results/S1/S1.tsv             all readouts of S1 (sample, biomarker, metric, value, unit, qc_flag, version)
results/S1/S1.status.tsv      one row per step (OK / FAILED / SKIPPED)
results/S1/<module>/          per-module results and params (provenance)
results/S1/logs/              one log per step
results/cohort.tsv            one row per sample (headline readouts + qc_flag)
results/cohort_long.tsv       every metric of every sample
```

`qc_flag` marks low-depth readouts, contaminated samples (`--max-contamination`, default 0.02),
samples with unusually many heteroplasmies (`--max-heteroplasmy`, default 20) and failed or
missing samples. Full definitions: [docs/module-contract.md](docs/module-contract.md).

## Tests

```bash
bash tests/test_heteroplasmy.sh
bash tests/test_pathogenic.sh
bash tests/test_carriers.sh
bash tests/test_cn_cluster.sh     # needs R package mclust
bash tests/test_prs.sh           # needs plink2
```

## Documentation

- [docs/module-contract.md](docs/module-contract.md) — inputs, outputs, metrics, cohort commands
- [docs/module-conf.md](docs/module-conf.md) — how modules are declared (adding a module)
- [CHANGELOG.md](CHANGELOG.md) — versions
- [legacy/](legacy/) — the original v0.1.0 single script and its WDL/Docker

## Citation

A manuscript is in preparation. Until then please cite this repository and the data sources
used by the modules you report ([NOTICE.md](NOTICE.md)).

## Licence

MIT ([LICENSE](LICENSE)). Bundled reference data keep their own terms ([NOTICE.md](NOTICE.md)).
