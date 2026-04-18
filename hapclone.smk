rule format_hapclone_cell:
    input:
        r=config.cnasim_reads,
        s=config.snp_counts_file,
    output:
        config.hapclone_cell_input_file,
    params:
        b=config.bin_size,
        k=config.block_size,
    conda:
        "envs/python.yaml"
    resources:
        mem="8G",
    shell:
        "python scripts/hapclone/format_hapclone_cell.py "
        "-c {wildcards.cell} "
        "-r {input.r} "
        "-s {input.s} "
        "-o {output} "
        "-b {params.b} "
        "-k {params.k} "


rule format_hapclone:
    input:
        expand(
            config.hapclone_cell_input_file, allow_missing=True, cell=config.cell_ids
        ),
    output:
        config.hapclone_input_file,
    params:
        " ".join(config.cell_ids),
    conda:
        "envs/python.yaml"
    resources:
        mem="16G",
    shell:
        "python scripts/hapclone/format_hapclone.py -i {input} -o {output} -c {params}"


rule hapclone_fit:
    input:
        d=config.hapclone_input_file,
    output:
        config.hapclone_fit_file,
    params:
        a=config.get_hapclone_cli_args,
        c=config.max_copy_state,
    # benchmark:
    #     config.get_benchmark_file("", "hapclone_init_fit.txt")
    conda:
        "envs/hapclone.yaml"
    # log:
    #     config.get_log_file(config.init_fit_file),
    resources:
        mem=lambda wildcards, attempt: "{}G".format(32 * attempt),
        runtime=lambda wildcards, attempt: "{}h".format(24 * attempt),
    threads: config.hapclone_num_threads
    shell:
        "hapclone fit "
        "-i {input.d} "
        "-o {output} "
        "-c {params.c} "
        "-t {threads} "
        "{params.a}"


rule hapclone_write_results:
    input:
        d=config.hapclone_input_file,
        f=config.hapclone_fit_file,
    output:
        config.hapclone_results_file,
    conda:
        "envs/hapclone.yaml"
    # log:
    #     config.get_log_file(config.results_file),
    resources:
        mem="16G",
        runtime="8h",
    shell:
        "hapclone write-results "
        "-d {input.d} "
        "-f {input.f} "
        "-o {output}"
