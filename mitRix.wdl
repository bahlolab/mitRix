version 1.0

struct IndexedCram {
    File cram
    File crai
}

struct IndexedFasta {
    File fasta
    File fai
}

task mitRix {
    input {
        IndexedCram cram
        IndexedFasta reference
        String sample_id

        # Ratio: "ND4/ND1", "ND5/ND1", "DLOOP/ND1", "mtDNA-CN", or "all"
        String ratio = "all"
        # Stat: "mean" or "median"
        String stat  = "mean"
        # mosdepth -F flag (1796 = exclude unmapped/non-primary/QC-fail/dup)
        Int    flag  = 1796
        Int    cpus  = 4

        String memory     = "16G"
        String disk_space = "local-disk ${disk_space_int} HDD"
        Int    preemptible = 2
    }

    Int disk_space_int = ceil(1.2 * size([
        cram.cram,
        cram.crai,
        reference.fasta,
        reference.fai
    ], "G"))

    String ratio_flag = if ratio == "all" then "--all" else "--ratio ~{ratio}"
    String stat_flag  = if stat == "median" then "--median" else "--mean"

    command {
        set -euo pipefail

        source /usr/local/bin/_activate_current_env.sh

        # CRAM index must sit beside the CRAM file
        if [[ ! -f "~{cram.cram}.crai" ]]; then
            ln -sf ~{cram.crai} ~{cram.cram}.crai
        fi

        mkdir -p outdir

        /app/mitRix \
            --bam             ~{cram.cram} \
            --reference-fasta ~{reference.fasta} \
            --sample          ~{sample_id} \
            ~{ratio_flag} \
            ~{stat_flag} \
            --flag            ~{flag} \
            --threads         ~{cpus} \
            --outdir          outdir \
            > ~{sample_id}.mitrix.tsv
    }

    output {
        File tsv              = "~{sample_id}.mitrix.tsv"
        File chrm_perbase     = "outdir/~{sample_id}.chrM.per-base.bed.gz"
        File chrm_perbase_csi = "outdir/~{sample_id}.chrM.per-base.bed.gz.csi"
    }

    runtime {
        docker:      "ghcr.io/bahlolab/mitrix:latest"
        cpu:         cpus
        memory:      memory
        disks:       disk_space
        preemptible: preemptible
    }
}

workflow mitRix_wf {
    call mitRix

    output {
        File tsv              = mitRix.tsv
        File chrm_perbase     = mitRix.chrm_perbase
        File chrm_perbase_csi = mitRix.chrm_perbase_csi
    }
}
