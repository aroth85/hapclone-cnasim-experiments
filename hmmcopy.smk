rule format_hmmcopy:
    conda:
        "envs/python.yaml"
    input:
        config.cnasim_reads,
    output:
        config.hmmcopy_input_template,
    resources:
        mem="32G",
    shell:
        "python scripts/hmmcopy/format_hmmcopy.py -i {input} -o {output} -c {wildcards.cell} -s {wildcards.sim_set}"


rule run_hmmcopy:
    input:
        config.hmmcopy_input_template,
    output:
        config.hmmcopy_output_template,
    params:
        config.max_copy_state,
    conda:
        "envs/hmmcopy.yaml"
    group:
        "cell"
    log:
        config.get_log_file(config.hmmcopy_output_template, config.hmmcopy_log_dir),
    benchmark:
        config.get_benchmark_file(
            config.hmmcopy_output_template, config.hmmcopy_benchmark_dir
        )
    resources:
        mem="16G",
    shell:
        "(python scripts/hmmcopy/run_hmmcopy.py -i {input} -o {output} -m 1 2 3 4 5 6 -s {params}) >{log} 2>&1"


rule write_sample_files:
    input:
        expand(config.hmmcopy_output_template, allow_missing=True, cell=config.cell_ids),
    output:
        p=config.hmmcopy_paths_template,
        m=config.hmmcopy_metrics_template,
        r=config.hmmcopy_reads_template,
        s=config.hmmcopy_segs_template,
    params:
        " ".join(config.cell_ids),
    conda:
        "envs/hmmcopy.yaml"
    log:
        config.get_log_file(config.hmmcopy_metrics_template, config.hmmcopy_log_dir),
    resources:
        mem="32G",
    shell:
        "(python scripts/hmmcopy/write_sample_files.py "
        "-a {wildcards.sim_set} "
        "-c {params} "
        "-i {input} "
        "-m {output.m} "
        "-p {output.p} "
        "-r {output.r} "
        "-s {output.s}) > {log} 2>&1"


rule plot_cnv_profiles:
    input:
        config.hmmcopy_output_template,
    output:
        config.hmmcopy_cnv_profile_template,
    conda:
        "envs/python.yaml"
    group:
        "cell"
    log:
        config.get_log_file(config.hmmcopy_cnv_profile_template, config.hmmcopy_log_dir),
    resources:
        mem="8G",
    shell:
        "(python scripts/hmmcopy/plot_cnv_profile.py -i {input} -o {output}) >{log} 2>&1"
