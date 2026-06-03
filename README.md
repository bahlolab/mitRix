# mitRix

A shell tool that estimates mtDNA damage / copy-number proxies as a
**coverage ratio** from a BAM/CRAM. Built to be the NGS analog of a
wet-lab TaqMan Ct ratio so the two read-outs are directly comparable.

```
./mitRix --bam sample.bam --ratio ND4/ND1
# sample  ratio    target            reference       stat  target_cov  ref_cov  ratio_value
# A1      ND4/ND1  chrM:12087-12170  chrM:3441-3504  mean  1284.51     1832.07  0.70
```

## Versions

- **Basic (this directory).** Wet-lab Ct mimic. Defaults to the published
  Thermo Fisher / Picard 4-plex TaqMan amplicon coordinates. Different labs
  can plug in their own targets via `--target` / `--target-bed` (and the
  reference counterparts).
- **Plus (planned).** Dedicated methodology for three damage axes — major-arc
  large deletion burden, 7S DNA / replication activity, and mtDNA copy number.
  Not implemented yet.

## What it does

1. **Header scan** — reads the BAM/CRAM `@SQ` header to auto-detect:
   - mtDNA contig name (`chrM` / `MT` / `chrMT` / `NC_012920.1`)
   - nuclear build (hg19 vs hg38) from `chr15` length, so the correct
     B2M coordinates are used for `mtDNA-CN` without manual input
2. **chrM per-base coverage** — runs `mosdepth -c <mt_contig>` once,
   producing `{sample}.chrM.per-base.bed.gz` (saved to `--outdir`).
   This single run covers all chrM-based ratios (ND4/ND1, ND5/ND1,
   DLOOP/ND1) regardless of whether `--all` or a single `--ratio` is used.
3. **Region extraction** — for each ratio, the target and reference
   intervals are resolved from built-in defaults (or user overrides via
   `--target` / `--reference` / `--target-bed` / `--reference-bed`) and
   extracted from the chrM per-base file.
4. **Coverage statistic** — mean (length-weighted, default) or median
   (`--median`) per-base depth is computed for each side. Read filtering
   is controlled by `--flag` (default 1796: exclude unmapped, non-primary,
   QC-failed, and duplicate reads).
5. **B2M (mtDNA-CN only)** — chr15 coverage for B2M is computed in a
   separate `mosdepth --by` call scoped to the B2M region only (uses the
   BAM index; does not scan all of chr15). Mean uses `-n`, median uses
   `--use-median`.
6. **Output** — one TSV line per ratio written to stdout.

## Install / dependencies

`mitRix` is a bash script. It needs:

- `samtools` (>= 1.10) — reads BAM/CRAM header and extracts reads
- `mosdepth` (>= 0.3) — coverage calculation
- `bwa` — required for `--ratio DLOOP/ND1` only (chrMC re-alignment)
- `awk`, `zcat`, `sort` (standard Unix)

```bash
# Install mosdepth + bwa via conda (recommended)
conda create -n mitrix -c bioconda mosdepth bwa
conda activate mitrix

# Make executable (once only)
chmod +x mitRix

./mitRix --help
# or put on PATH:
mv mitRix ~/.local/bin/
```

**HPC (module system):**
```bash
module load samtools/1.23
module load bwa
module load miniconda3/
conda activate mitrix
```

**One-time setup for `DLOOP/ND1`** — build the chrMC circularised reference:
```bash
bash refs/make_chrMC.sh
```
This downloads rCRS (NC_012920.1), appends 700 bp (≈ 2× mean insert size) to
create `refs/chrMC_rCRS.fa` (17,269 bp), and builds the BWA index. Only needed
once; subsequent runs detect the existing index and skip this step.

For CRAM inputs, pass `--reference-fasta`.

## Ratios and built-in defaults

Default ratio is **`ND4/ND1`**. All default regions follow the published
Thermo Fisher / Picard 4-plex TaqMan amplicons on the rCRS reference
(NC_012920.1, 16,569 bp). Coordinates are 1-based closed.

| `--ratio`             | target (default)                         | reference (default)               | Thermo probe                              |
| --------------------- | ---------------------------------------- | --------------------------------- | ----------------------------------------- |
| `ND4/ND1` *(default)* | ND4, `chrM:12087-12170` (84 bp)          | ND1, `chrM:3441-3504` (64 bp)     | ND4-ABY/QSY (Custom CCU002NR)             |
| `ND5/ND1`             | ND5, `chrM:12874-12969` (96 bp) †        | ND1, `chrM:3441-3504` (64 bp)     | ND5 not in 4-plex; ND1 = APAANDG VIC/MGB |
| `DLOOP/ND1`           | D-loop, `chrM:1-34` + `chrM:16536-16569` | ND1, `chrM:3441-3504` (64 bp)     | DLOOP-JUN/QSY (Custom CCU002NR)           |
| `mtDNA-CN`            | ND1, `chrM:3441-3504` (64 bp)            | B2M exon 2 (chr15, auto-detected) | B2M-FAM/MGB (analog of exon 2 amplicon)   |
| `custom`              | `--target` / `--target-bed`              | `--reference` / `--reference-bed` | bring your own regions                    |

† ND5 is not in the Thermo 4-plex; the default is a reasonable major-arc
window provided so the `ND5/ND1` ratio works out of the box. Override with
`--target` / `--target-bed` if your lab uses a different ND5 amplicon.

B2M defaults (exon 2 of NM_004048, the same exon the B2M-FAM/MGB qPCR probe
sits in):

| Build         | Region (1-based closed)             | Length |
| ------------- | ----------------------------------- | ------ |
| hg38 / GRCh38 | `chr15:44,715,439-44,715,720`       | 282 bp |
| hg19 / GRCh37 | `chr15:45,007,637-45,007,918`       | 282 bp |

## Options

| Option | Default | Description |
| --- | --- | --- |
| `--bam FILE` | — | BAM or CRAM (indexed); **required** |
| `--ratio NAME` | `ND4/ND1` | Ratio to compute (see table above) |
| `--all` | off | Run all four ratios in one pass |
| `--mean` | ✓ | Length-weighted mean per-base coverage |
| `--median` | off | Median per-base coverage |
| `--flag INT` | `1796` | mosdepth `-F`: exclude reads with any of these SAM flag bits set |
| `--sample NAME` | BAM basename | Sample name in output |
| `--outdir DIR` | tmp dir | Save `{sample}.chrM.per-base.bed.gz` here |
| `--reference-fasta FA` | — | Required for CRAM input |
| `--threads N` | `4` | Threads for mosdepth |
| `--no-header` | off | Suppress TSV header (for appending) |
| `--target START-END` | — | Override target interval on chrM (1-based closed) |
| `--reference START-END` | — | Override reference interval on chrM (1-based closed) |
| `--target-bed FILE` | — | Override target as a BED file (any chromosome) |
| `--reference-bed FILE` | — | Override reference as a BED file |
| `--print-bed` | off | Print resolved BED files to stderr |
| `--rotate-perbase` | off | Save `{sample}.chrMC.per-base.bed.gz` and `{sample}.chrM.corrected.per-base.bed.gz` to `--outdir` (DLOOP/ND1 only; useful for QC plots) |
| `--keep-tmp` | off | Do not delete temp dir |

**Common `--flag` values:**

| `--flag` | Excludes |
| --- | --- |
| `1796` | unmapped + not-primary + QC-fail + **duplicate** (default) |
| `2820` | unmapped + not-primary + QC-fail + **supplementary** (keep duplicates) |
| `3844` | unmapped + not-primary + QC-fail + duplicate + supplementary |

## Examples

```bash
# Default: ND4/ND1, mean coverage
./mitRix --bam sample.bam

# Median instead of mean
./mitRix --bam sample.bam --median

# Exclude supplementary reads on top of the default
./mitRix --bam sample.bam --flag 3844

# Keep duplicates (higher coverage for both sides)
./mitRix --bam sample.bam --flag 2820

# Specific ratios
./mitRix --bam sample.bam --ratio ND5/ND1
./mitRix --bam sample.bam --ratio ND5/ND1 --median
./mitRix --bam sample.bam --ratio DLOOP/ND1
./mitRix --bam sample.bam --ratio DLOOP/ND1 --median
./mitRix --bam sample.bam --ratio mtDNA-CN
./mitRix --bam sample.bam --ratio mtDNA-CN --median

# All four ratios in one pass
./mitRix --bam sample.bam --all > results.tsv
./mitRix --bam sample.bam --all --median > results_median.tsv

# CRAM input (reference FASTA required)
./mitRix --bam sample.cram --reference-fasta /refs/hg38.fa --ratio ND4/ND1
./mitRix --bam sample.cram --reference-fasta /refs/hg38.fa --all --median

# Custom target region
./mitRix --bam sample.bam --ratio ND4/ND1 --target 11333-11482 --reference 3441-3504

# Custom B2M amplicon
./mitRix --bam sample.bam --ratio mtDNA-CN --reference-bed my_B2M.bed

# Batch across multiple BAMs (output to v0/, results to TSV)
BAM_DIR=/path/to/bams
first=1
for bam in "$BAM_DIR"/*.bam; do
  sample=$(basename "$bam" .bam | cut -d. -f1)   # e.g. A1.merged.bam → A1
  if (( first )); then
    ./mitRix --bam "$bam" --sample "$sample" --all --outdir v0 > results.tsv
    first=0
  else
    ./mitRix --bam "$bam" --sample "$sample" --all --outdir v0 --no-header >> results.tsv
  fi
done
```

## Output

TSV with 8 columns written to stdout:

```
sample  ratio  target  reference  stat  target_cov  ref_cov  ratio_value
```

- `target` / `reference` encode all intervals, e.g. `chrM:1-34,16536-16569`
  for the D-loop pair.
- `stat` is `mean` or `median` depending on `--mean` / `--median`.
- `target_cov` / `ref_cov` are rounded to 2 decimal places.
- `ratio_value` is `target_cov / ref_cov` (2 decimals; `NA` if either side
  is zero or undefined).

Per-base chrM coverage is saved as `{sample}.chrM.per-base.bed.gz` (and its
`.csi` index) in `--outdir` if given, otherwise the current directory.

## DLOOP/ND1 — chrMC re-alignment

The D-loop amplicon (`chrM:16536-16569` + `chrM:1-34`) spans the linear
reference break. In a standard BAM, reads that physically cross the 16569→1
seam are soft-clipped at the edge, causing artificially low coverage at both
ends of the amplicon and inflating the DLOOP/ND1 ratio by ~1.16×.

`mitRix` fixes this by re-aligning all chrM reads to a circularised reference
(**chrMC**) built by appending the first 700 bp of rCRS to the end of the
genome (16,569 + 700 = 17,269 bp). The 700 bp extension is ≈ 2× the typical
WGS mean insert size, ensuring >99.9% of seam-crossing read pairs are fully
contained within the reference. After re-alignment with BWA in paired-end mode:

- **arm1** (`chrMC:16536-16569`) — single unambiguous copy of the left arm
- **arm2a** (`chrMC:1-34`) — reads BWA placed in the start copy of the right arm
- **arm2b** (`chrMC:16570-16603`) — reads BWA placed in the appended seam copy

D-loop depth = `(arm1_mean + arm2a_mean + arm2b_mean) / 2`

ND1 is extracted from a corrected per-base file that folds the chrMC appended
region back to chrM coordinates (modulo-sum). ND1 is at 3441–3504, well outside
the appended zone (1–700 bp), so it is unaffected.

Use `--rotate-perbase` to save the chrMC and corrected per-base files for QC
inspection, then plot with `Rscript plot_dloop.R <outdir> <sample>`.

## Notes / caveats

- The default `--flag 1796` excludes unmapped, non-primary, QC-failed, and
  duplicate reads — mirroring the deduplication step typical in WGS pipelines.
  Use `--flag 2820` to keep duplicates (closer to raw qPCR signal without
  dedup), or `--flag 3844` to additionally exclude supplementary alignments.
- The D-loop default amplicon wraps the origin of replication (`chrM:16536-16569`
  + `chrM:1-34`). `DLOOP/ND1` automatically re-aligns reads to the chrMC
  circularised reference (see [DLOOP/ND1 — chrMC re-alignment](#dloopnd1--chrmc-re-alignment)
  above) to avoid soft-clipping artefacts at the seam. If you supply a custom
  `--target` that wraps the origin, mitRix will auto-split the interval at the
  seam but will fall back to the original chrM per-base (no chrMC correction)
  for non-default amplicons.
- For a chrM-only BAM (no chr15), `mtDNA-CN` is automatically skipped when
  running `--all`. Supply `--reference-bed` to use a custom nuclear reference.
- For T2T-CHM13 or a non-rCRS mt reference, supply `--target-bed` /
  `--reference-bed` with coordinates appropriate to that assembly.

## Source of default coordinates

- Trumpff C et al. *Stress and circulating cell-free mitochondrial DNA: a
  systematic review of human studies, physiological considerations, and
  technical recommendations.* (PMID: 33839318)
- MitoQuicLy: high-throughput cf-mtDNA quantification (PMC: 9882007)
- Sazonova MA et al. *Simultaneous quantification of mtDNA copy number and
  deletion ratio: a multiplex real-time PCR assay.* Sci Rep 4: 3887 (2014)
- The Picard Lab cell-free mtDNA assay (Columbia University)
