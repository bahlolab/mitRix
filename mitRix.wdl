version 1.0

task mitRix {
    input {
        File   cram
        File   crai
        File   reference_fasta
        File   reference_fasta_fai
        String sample_id

        # Ratio: "ND4/ND1", "ND5/ND1", "DLOOP/ND1", "mtDNA-CN", or "all"
        String ratio = "all"
        # Stat: "mean" or "median"
        String stat  = "mean"
        # mosdepth -F flag
        Int    flag  = 1796
        Int    cpus  = 4

        String memory      = "16G"
        Int    preemptible = 2
    }

    Int disk_space_int = ceil(1.2 * size([cram, crai, reference_fasta, reference_fasta_fai], "G"))

    String ratio_flag = if ratio == "all" then "--all" else "--ratio ~{ratio}"
    String stat_flag  = if stat == "median" then "--median" else "--mean"

    command {
        set -euo pipefail

        source /usr/local/bin/_activate_current_env.sh

        # CRAM index must sit beside the CRAM file
        ln -sf ~{crai} ~{cram}.crai

        mkdir -p outdir

        /app/mitRix \
            --bam             ~{cram} \
            --reference-fasta ~{reference_fasta} \
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
        disks:       "local-disk ~{disk_space_int} HDD"
        preemptible: preemptible
    }
}

workflow mitRix_wf {
    input {
        File   cram
        File   crai
        File   reference_fasta
        File   reference_fasta_fai
        String sample_id
        String ratio = "all"
        String stat  = "mean"
        Int    flag  = 1796
        Int    cpus  = 4
    }

    call mitRix {
        input:
            cram                = cram,
            crai                = crai,
            reference_fasta     = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            sample_id           = sample_id,
            ratio               = ratio,
            stat                = stat,
            flag                = flag,
            cpus                = cpus
    }

    output {
        File tsv              = mitRix.tsv
        File chrm_perbase     = mitRix.chrm_perbase
        File chrm_perbase_csi = mitRix.chrm_perbase_csi
    }
}
