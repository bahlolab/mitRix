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
        Int    preemptible = 0
    }

    Int disk_space_int = ceil(1.2 * size([cram, crai, reference_fasta, reference_fasta_fai], "G"))

    String ratio_flag = if ratio == "all" then "--all" else "--ratio ~{ratio}"
    String stat_flag  = if stat == "median" then "--median" else "--mean"

    command {
        set -eo pipefail

        export MAMBA_SKIP_ACTIVATE=0
        source /usr/local/bin/_activate_current_env.sh

        # CRAM index must sit beside the CRAM file
        if [[ ! -f "~{cram}.crai" ]]; then
            ln -sf ~{crai} ~{cram}.crai
        fi

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
            --rotate-perbase \
            > ~{sample_id}.mitrix.tsv

        tar -czf ~{sample_id}.perbase.tar.gz -C outdir .
    }

    output {
        File tsv       = "~{sample_id}.mitrix.tsv"
        File perbase   = "~{sample_id}.perbase.tar.gz"
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
        File mitRix_result = mitRix.tsv
        File base_cov2 = mitRix.perbase
    }
}
