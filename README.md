# mitRix

Tiny shell tool that estimates mtDNA damage as a **mean-coverage ratio**
between a target region and a reference region in a BAM/CRAM.

```
mitrix --bam sample.bam --ratio ND4/ND1
# sample  ratio    target           reference        target_mean_cov  ref_mean_cov  ratio_value
# sample  ND4/ND1  chrM:11333-11482 chrM:3441-3504   1284.51          1832.07       0.701234
```

## What it does

1. Reads the BAM/CRAM `@SQ` header to detect:
   - the mtDNA contig name (`chrM` / `MT` / `chrMT` / `NC_012920.1`)
   - the nuclear build (hg19 vs hg38) from the `chr15` length, so the right
     B2M coordinates are used for `mtDNA-CN` without you having to think
2. Writes per-region BED files for the chosen ratio
3. Runs `mosdepth -n --by <bed>` for both regions
4. Emits a length-weighted mean coverage for each side and the ratio,
   as one TSV line on stdout

## Install / dependencies

`mitrix` is a bash script. It needs:

- `samtools` (>= 1.10)
- `mosdepth` (>= 0.3)
- `awk`, `zcat`, `sort` (standard Unix)

```
chmod +x mitrix
./mitrix --help
# or put it on PATH:
mv mitrix ~/.local/bin/
```

For CRAM inputs, pass `--reference-fasta`.

## Ratios

Default ratio is **`ND4/ND1`**.

| `--ratio`             | target (default)              | reference (default)                              | meaning                       |
| --------------------- | ----------------------------- | ------------------------------------------------ | ----------------------------- |
| `ND4/ND1` *(default)* | ND4, `chrM:11333-11482`       | ND1, `chrM:3441-3504`                            | major-arc deletion damage     |
| `ND5/ND1`             | ND5, `chrM:12784-12969`       | ND1, `chrM:3441-3504`                            | major-arc deletion damage     |
| `DLOOP/ND1`           | D-loop, `chrM:1-576` + `chrM:16024-16569` | ND1, `chrM:3441-3504`                | control-region coverage       |
| `mtDNA-CN`            | ND1, `chrM:3441-3504`         | B2M, `chr15:44,711,477-44,718,877` (hg38) / `chr15:45,003,675-45,011,087` (hg19) | mtDNA copy number proxy       |
| `custom`              | `--target` / `--target-bed`   | `--reference` / `--reference-bed`                | bring your own regions        |

All defaults are **1-based, closed** intervals on the rCRS mitochondrial
reference (NC_012920.1, 16,569 bp). The tool warns if the mtDNA contig has
a different length.

## Examples

```
# default ND4/ND1, auto-detect chrM
mitrix --bam sample.bam

# mtDNA copy number — auto-detect hg19/hg38 from chr15 length for B2M
mitrix --bam sample.bam --ratio mtDNA-CN

# CRAM input
mitrix --bam sample.cram --reference-fasta /refs/hg38.fa --ratio ND5/ND1

# Custom B2M coords (e.g. T2T-CHM13, or a different transcript)
mitrix --bam sample.bam --ratio mtDNA-CN --reference-bed bed/B2M_hg38.bed

# Pure custom: pick any two mtDNA windows
mitrix --bam sample.bam --ratio custom \
       --target 8470-13447 --reference 3441-3504

# Batch across BAMs, append to a table
mitrix --bam s1.bam --sample s1 > out.tsv
for b in s*.bam; do
  mitrix --bam "$b" --sample "${b%.bam}" --no-header --ratio ND4/ND1 >> out.tsv
done
```

## Bundled BED files

`bed/` holds ready-to-use BED files (1-based-closed coordinates documented in
the filename; BED contents are the standard 0-based half-open):

```
bed/ND1.rCRS.bed         chrM:3441-3504
bed/ND4.rCRS.bed         chrM:11333-11482
bed/ND5.rCRS.bed         chrM:12784-12969
bed/Dloop.rCRS.bed       chrM:1-576 + chrM:16024-16569
bed/B2M_hg38.bed         chr15:44,711,477-44,718,877
bed/B2M_hg19.bed         chr15:45,003,675-45,011,087
```

Each is provided with both `chrM`/`chr15` and `MT`/`15` variants in the
`bed/ucsc/` and `bed/ensembl/` subdirs respectively.

## Output

TSV with 7 columns:

```
sample  ratio  target  reference  target_mean_cov  ref_mean_cov  ratio_value
```

- `target` / `reference` columns encode all intervals used, e.g.
  `chrM:1-576,16024-16569` for the D-loop pair.
- `*_mean_cov` is **length-weighted** across multiple intervals.
- `ratio_value` is `target_mean_cov / ref_mean_cov` (6 decimals, `NA` if
  either side has zero coverage or undefined).

## Notes / caveats

- mosdepth and samtools both apply default filtering. mtDNA reads from short-read
  WGS often have many secondary/supplementary alignments — keep this in mind
  if you want stricter filtering (you can wrap `mitrix` in your own samtools
  pipeline or extend the script with `--min-mapq` / `--exclude-flag` if needed).
- `D-loop / ND1` is sensitive to where the linear chrM reference is broken;
  read pile-up at positions 1 and 16569 can be uneven. The default uses the
  full canonical control region (1-576 + 16024-16569).
- For T2T-CHM13 or a non-rCRS mt reference, supply `--target-bed` /
  `--reference-bed` with coordinates appropriate to that reference.
