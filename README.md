# mitRix

A small shell tool that estimates mtDNA damage / copy-number proxies as a
**mean-coverage ratio** from a BAM/CRAM. Built to be the NGS analog of a
wet-lab TaqMan Ct ratio so the two read-outs are directly comparable.

```
mitrix --bam sample.bam --ratio ND4/ND1
# sample  ratio    target            reference         target_mean_cov  ref_mean_cov  ratio_value
# sample  ND4/ND1  chrM:12087-12170  chrM:3485-3553    1284.51          1832.07       0.701234
```

## Versions

- **Basic (this directory).** Wet-lab Ct mimic. Defaults to the published
  Thermo Fisher / Picard 4-plex TaqMan amplicon coordinates. Different labs
  can plug in their own targets via `--target` / `--target-bed` (and the
  reference counterparts).
- **Plus (planned).** Same ratios on the same regions, but with methodology
  improvements layered on top — GC-bias correction, NUMT/MAPQ filtering,
  duplicate-aware depth, bootstrap CIs, etc. Not implemented yet.

## What it does

1. Reads the BAM/CRAM `@SQ` header to detect:
   - the mtDNA contig name (`chrM` / `MT` / `chrMT` / `NC_012920.1`)
   - the nuclear build (hg19 vs hg38) from the `chr15` length, so the right
     B2M amplicon coordinates are used for `mtDNA-CN` without intervention
2. Writes per-region BED files for the chosen ratio
3. Runs `mosdepth -n --by <bed>` for both regions
4. Emits length-weighted mean coverage for each side and the ratio, as one
   TSV line on stdout

## Install / dependencies

`mitrix` is a bash script. It needs:

- `samtools` (>= 1.10) — used only to read the BAM/CRAM header
- `mosdepth` (>= 0.3) — does the actual coverage calculation
- `awk`, `zcat`, `sort` (standard Unix)

```
chmod +x mitrix
./mitrix --help
# or put it on PATH:
mv mitrix ~/.local/bin/
```

For CRAM inputs, pass `--reference-fasta`.

## Ratios and built-in defaults

Default ratio is **`ND4/ND1`**. All default regions are the published
Thermo Fisher / Picard 4-plex TaqMan amplicons on the rCRS reference
(NC_012920.1, 16,569 bp). 1-based, closed.

| `--ratio`             | target (default)                              | reference (default)                | Thermo probe                            |
| --------------------- | --------------------------------------------- | ---------------------------------- | --------------------------------------- |
| `ND4/ND1` *(default)* | ND4, `chrM:12087-12170` (84 bp)               | ND1, `chrM:3485-3553` (69 bp)      | ND4-ABY/QSY (Custom CCU002NR)           |
| `ND5/ND1`             | ND5, `chrM:12784-12969` (186 bp) †            | ND1, `chrM:3485-3553` (69 bp)      | ND5 not in 4-plex; ND1 = APAANDG VIC/MGB |
| `DLOOP/ND1`           | D-loop, `chrM:1-34` + `chrM:16536-16569`      | ND1, `chrM:3485-3553` (69 bp)      | DLOOP-JUN/QSY (Custom CCU002NR)         |
| `mtDNA-CN`            | ND1, `chrM:3485-3553` (69 bp)                 | B2M exon 2 (chr15, auto-detected)  | B2M-FAM/MGB (analog of exon 2 amplicon) |
| `custom`              | `--target` / `--target-bed`                   | `--reference` / `--reference-bed`  | bring your own regions                  |

† ND5 is not in the Thermo 4-plex; the default is a reasonable major-arc
window provided so the `ND5/ND1` ratio works out of the box. Override with
`--target` / `--target-bed` if your lab uses a different ND5 amplicon.

B2M defaults (exon 2 of NM_004048, the same exon the B2M-FAM/MGB qPCR probe
sits in):

| Build | Region (1-based closed) | Length |
| --- | --- | --- |
| hg38 / GRCh38 | `chr15:44,715,439-44,715,720` | 282 bp |
| hg19 / GRCh37 | `chr15:45,007,637-45,007,918` | 282 bp |

## Customising for a different probe design

Every default is overridable. Different labs use different amplicons; mitrix
takes a per-side override:

```
# mtDNA range (1-based closed) on the auto-detected mt contig
mitrix --bam s.bam --ratio ND4/ND1 \
       --target 11333-11482 \
       --reference 3441-3504

# arbitrary BED files (any chromosomes; multi-interval allowed)
mitrix --bam s.bam --ratio mtDNA-CN \
       --target-bed my_ND1.bed \
       --reference-bed my_B2M.bed

# fully custom: bring both sides
mitrix --bam s.bam --ratio custom \
       --target-bed my_target.bed \
       --reference-bed my_reference.bed
```

## Examples

```
# Default ND4/ND1 with Thermo amplicons, auto-detect chrM naming
mitrix --bam sample.bam

# mtDNA copy number — auto-detect hg19/hg38 from chr15 length for B2M
mitrix --bam sample.bam --ratio mtDNA-CN

# CRAM input
mitrix --bam sample.cram --reference-fasta /refs/hg38.fa --ratio DLOOP/ND1

# Lab uses a different B2M amplicon
mitrix --bam sample.bam --ratio mtDNA-CN --reference-bed bed/my_B2M.bed

# Batch across BAMs, append to a table
mitrix --bam s1.bam --sample s1 > out.tsv
for b in s*.bam; do
  mitrix --bam "$b" --sample "${b%.bam}" --no-header --ratio ND4/ND1 >> out.tsv
done
```

## Bundled BED files

`bed/` holds ready-to-use BED files (BED is 0-based half-open; the
1-based-closed coordinates are documented above):

```
bed/ucsc/ND1.rCRS.bed       chrM:3485-3553
bed/ucsc/ND4.rCRS.bed       chrM:12087-12170
bed/ucsc/Dloop.rCRS.bed     chrM:1-34 + chrM:16536-16569
bed/ucsc/B2M_hg38.bed       chr15:44,715,439-44,715,720  (exon 2)
bed/ucsc/B2M_hg19.bed       chr15:45,007,637-45,007,918  (exon 2)
```

Ensembl-style copies (`MT` / `15`) sit in `bed/ensembl/`.

## Output

TSV with 7 columns:

```
sample  ratio  target  reference  target_mean_cov  ref_mean_cov  ratio_value
```

- `target` / `reference` encode all intervals used, e.g.
  `chrM:1-34,16536-16569` for the D-loop pair.
- `*_mean_cov` is **length-weighted** across multiple intervals.
- `ratio_value` is `target_mean_cov / ref_mean_cov` (6 decimals; `NA` if
  either side is zero / undefined).

## Notes / caveats

- The Basic version applies **no filtering** — it uses `mosdepth -n --by`
  with default settings to stay as close to the raw qPCR amplification
  signal as possible. NUMTs, PCR duplicates, low-MAPQ reads, and GC bias
  are all left in. The Plus version will address these explicitly.
- The D-loop default amplicon wraps the origin of replication on the
  linear chrM reference (chrM:16536-16569 → chrM:1-34). Pile-up at the
  reference break can be uneven; this is inherent to short-read mt
  mapping. If your aligner used the rotated reference (origin shifted to
  the middle of the linear chrM), provide a `--target-bed` matching that
  layout.
- For T2T-CHM13 or a non-rCRS mt reference, supply `--target-bed` /
  `--reference-bed` with coordinates appropriate to that reference.

## Source of default coordinates

- Trumpff C et al. *Stress and circulating cell-free mitochondrial DNA: a
  systematic review of human studies, physiological considerations, and
  technical recommendations.* (PMID: 33839318)
- MitoQuicLy: high-throughput cf-mtDNA quantification (PMC: 9882007)
- Sazonova MA et al. *Simultaneous quantification of mtDNA copy number and
  deletion ratio: a multiplex real-time PCR assay.* Sci Rep 4: 3887 (2014)
- The Picard Lab cell-free mtDNA assay (Columbia University)
