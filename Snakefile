import os


from utils import ConfigManager

config = ConfigManager(config)


rule all:
    input:
        list(config.get_pipeline_files()),


rule cnasim:
    params:
        b=config.bin_size,
        c=config.coverage,
        n=config.num_cells,
        o=str(config.cnasim_out_dir),
        r=config.read_length,
        cli_args=config.get_cna_cli_args,
    output:
        config.cnasim_reads,
        config.cnasim_profiles,
    conda:
        "envs/cnasim.yaml"
    resources:
        mem="16G",
    threads: config.cnasim_num_threads
    shell:
        "cnasim "
        "-m 1 "
        "-B {params.b} "
        "-C {params.c} "
        "-n {params.n}  "
        "-o {params.o} "
        "-P {threads} "
        "-R {params.r} "
        "-U "
        "{params.cli_args}"


rule snps:
    input:
        config.cnasim_reads,
    output:
        config.snps_file,
    params:
        n=config.num_snps,
        p=config.phase_switch_prob,
    conda:
        "envs/python.yaml"
    resources:
        mem="16G",
    shell:
        "python scripts/simulate_snps.py -i {input} -o {output} -n {params.n} -p {params.p}"


rule simulate_snp_counts:
    input:
        p=config.cnasim_profiles,
        r=config.cnasim_reads,
        s=config.snps_file,
    output:
        config.snp_counts_file,
    params:
        r=config.read_length,
    conda:
        "envs/python.yaml"
    group:
        "cell"
    resources:
        mem="16G",
    shell:
        "python scripts/simulate_snp_counts.py "
        "-c {wildcards.cell} "
        "-p {input.p} "
        "-r {input.r} "
        "-s {input.s} "
        "-o {output} "
        "--read-length {params}"


include: "chisel.smk"
include: "hapclone.smk"
include: "hmmcopy.smk"
include: "signals.smk"
