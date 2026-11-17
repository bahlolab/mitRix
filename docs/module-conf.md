# module.conf format

Each biomarker module has a `biomarkers/<name>/module.conf`. `bin/mitrix` reads it to
work out run order, check that tools are installed, and list modules in `--help`.
Nextflow and WDL generators (P6) read the same file.

Status: draft for ROADMAP p0-5. The run-time contract is in
[module-contract.md](module-contract.md).

## Syntax

- One `key = value` per line. Spaces around `=` are optional.
- `#` starts a comment, either on its own line or after a value.
- List values are separated by spaces.
- Keys are lower case. Unknown keys are an error, so typos are caught.
- The file is **parsed, never sourced**: it can't run code.

## Keys

| key | required | example | meaning |
|---|---|---|---|
| `name` | yes | `damage` | must match the directory name and the `biomarker` column |
| `version` | yes | `0.1.0` | module version (semver); written to the `version` column |
| `description` | yes | `Major-arc damage ratios (ND4/ND1, ND5/ND1)` | one line, shown in `mitrix --help` |
| `entry` | yes | `run.sh` | executable, relative to the module directory |
| `language` | no | `bash` | `bash` (default) or `R` |
| `requires_core` | no | `ref_context mt_coverage` | core steps from `core/` that must run first |
| `requires_modules` | no | `heteroplasmy` | modules whose outputs this one reads; they run first |
| `requires_inputs` | no | `mt_variants` | inputs needed: `mt_variants` (module writes `SKIPPED` rows without it) or `bam` (module is SKIPPED in variants-only mode) |
| `tools` | no | `samtools>=1.10 mosdepth>=0.3` | external tools, with an optional minimum version |
| `r_packages` | no | `data.table` | R packages (for `language = R`) |
| `options` | no | `--target --reference --stat` | module-specific flags, listed in `--help` |
| `cohort` | no | `deletion:ND4_ND1_ratio_mean` | metrics copied into `cohort.tsv`, as `column:metric` pairs |

## Rules `bin/mitrix` enforces

- `requires_core` names must exist in `core/`.
- `requires_modules` must not form a cycle. Modules run in dependency order, and
  independent modules may run in parallel.
- Before any work starts, every `tools` entry is checked on `PATH` with its version.
  `haplogrep3` and `haplocheck` are Java tools installed by `tools/install.sh`; they are
  checked at `MITRIX_HAPLOGREP3` / `MITRIX_HAPLOCHECK` (default `tools/*/…jar`).
  A missing or too-old tool stops the run with a message naming the module.
- `--biomarkers a,b` also runs anything `a` and `b` require.

## Examples

`biomarkers/damage/module.conf`:

```
name             = damage
version          = 0.1.0
description      = Major-arc damage ratios (ND4/ND1, ND5/ND1, custom)
entry            = run.sh
requires_core    = ref_context mt_coverage
tools            = samtools>=1.10 mosdepth>=0.3
options          = --ratios --target --reference --target-bed --reference-bed --stat
cohort           = deletion:ND4_ND1_ratio_mean
```

`biomarkers/d_loop_7s/module.conf`:

```
name             = d_loop_7s
version          = 0.1.0
description      = 7S DNA / replication proxy (DLOOP/ND1) on circularised chrM
entry            = run.sh
requires_core    = ref_context mt_coverage mt_circular
tools            = samtools>=1.10 mosdepth>=0.3 bwa
options          = --stat
cohort           = 7S:DLOOP_ND1_ratio_mean
```

`biomarkers/haplogroup/module.conf` (P3):

```
name             = haplogroup
version          = 0.1.0
description      = mtDNA haplogroup with HaploGrep3
entry            = run.sh
requires_modules = heteroplasmy   # reads the standard variant table
requires_inputs  = mt_variants
tools            = haplogrep3
cohort           = haplogroup:haplogroup
```

Two modules may not claim the same cohort column name; `bin/mitrix` checks this.
