# legacy/

`v0.1.0/` holds the original single-script mitRix (`mitRix`: coverage ratios ND4/ND1,
ND5/ND1, DLOOP/ND1, mtDNA-CN from one BAM/CRAM), its WDL workflow (`mitRix.wdl`), Docker
image definition (`docker/`) and the GitHub Action that built it (`docker.yml`, disabled
here). They are kept for reproducing results made with v0.1.0; the same ratios are now
produced by the `damage`, `d_loop_7s` and `copy_number` modules of `bin/mitrix`.
A Docker image and WDL workflow for the modular pipeline are planned.
