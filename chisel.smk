# CHISEL
rule format_chisel_baf:
    conda:
        "envs/python.yaml"
    input:
        expand(config.snp_counts_file, allow_missing=True, cell=config.cell_ids),
    output:
        config.chisel_baf_template,
    resources:
        mem="16G",
    shell:
        "python scripts/chisel/format_chisel_baf.py -i {input} -o {output}"


rule format_chisel_rdr:
    conda:
        "envs/python.yaml"
    input:
        config.cnasim_reads,
    resources:
        mem="16G",
    output:
        r=config.chisel_rdr_template,
        t=config.chisel_total_template,
    shell:
        "python scripts/chisel/format_chisel_rdr.py -i {input} -r {output.r} -t {output.t}"


checkpoint split_rdr_by_chunk:
    input:
        config.chisel_rdr_template,
    output:
        directory(config.chisel_working_dir.joinpath("rdr")),
    params:
        config.chisel_chunk_size,
    conda:
        "envs/python.yaml"
    log:
        config.chisel_log_dir.joinpath("split_rdr.log"),
    benchmark:
        config.chisel_benchmark_dir.joinpath("split_rdr.txt")
    resources:
        mem="16G",
    shell:
        "(python scripts/chisel/split_rdr_by_chunk.py "
        "-i {input} "
        "-o {output} "
        "-c {params} ) >{log} 2>&1"


rule split_baf_by_chunk:
    input:
        b=config.chisel_baf_template,
        r=config.chisel_chunk_rdr_template,
    output:
        config.chisel_chunk_baf_template,
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.chisel_chunk_baf_template),
    benchmark:
        config.get_benchmark_file(
            config.chisel_chunk_baf_template, config.chisel_benchmark_dir
        )
    resources:
        mem="16G",
    shell:
        "(python scripts/chisel/split_baf_by_chunk.py "
        "-b {input.b} "
        "-r {input.r} "
        "-o {output} ) >{log} 2>&1"


rule filter_cells:
    input:
        config.chisel_total_template,
    output:
        config.chisel_total_filtered_file,
    params:
        config.chisel_min_cell_reads,
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.chisel_total_filtered_file),
    benchmark:
        config.get_benchmark_file(
            config.chisel_total_filtered_file, config.chisel_benchmark_dir
        )
    resources:
        mem="16G",
    shell:
        "(python scripts/chisel/filter_cells.py "
        "-i {input} "
        "-o {output} "
        "-r {params} ) >{log} 2>&1"


rule run_chisel_combiner:
    input:
        b=config.chisel_chunk_baf_template,
        r=config.chisel_chunk_rdr_template,
        t=config.chisel_total_filtered_file,
    output:
        c=config.chisel_chunk_combo_template,
        s=config.chisel_chunk_swap_template,
    params:
        b=config.block_size,
        n=config.chisel_combo_restarts,
    conda:
        "envs/chisel.yaml"
    log:
        config.get_log_file(config.chisel_chunk_combo_template),
    benchmark:
        config.get_benchmark_file(
            config.chisel_chunk_swap_template, config.chisel_benchmark_dir
        )
    resources:
        mem="16G",
    threads: config.chisel_combo_threads
    shell:
        "(python2.7 scripts/chisel/Combiner.py "
        "-b {input.b} "
        "-l {input.t} "
        "-r {input.r} "
        "-o {output.c} "
        "--allele-swap-file {output.s} "
        "-k {params.b} "
        "-q {params.n} "
        "-j {threads}) >{log} 2>&1"


def get_chunk_combo_files(wildcards):
    checkpoint_output = checkpoints.split_rdr_by_chunk.get(**wildcards).output[0]
    return expand(
        config.chisel_chunk_combo_template,
        sim_set=wildcards.sim_set,
        replicate=wildcards.replicate,
        chunk=glob_wildcards(os.path.join(checkpoint_output, "{chunk}.tsv")).chunk,
    )


rule build_combo_file:
    input:
        get_chunk_combo_files,
    output:
        config.chisel_combo_file,
    log:
        config.get_log_file(config.chisel_combo_file),
    resources:
        mem="16G",
    benchmark:
        config.get_benchmark_file(config.chisel_combo_file, config.chisel_benchmark_dir)
    threads: config.chisel_sort_threads
    shell:
        "(sort -m -k 1,1n -k 2,2n -k 4,4 --parallel {threads} -S 8G -o {output} {input}) >{log} 2>&1"


def get_chunk_swap_files(wildcards):
    checkpoint_output = checkpoints.split_rdr_by_chunk.get(**wildcards).output[0]
    return expand(
        config.chisel_chunk_swap_template,
        sim_set=wildcards.sim_set,
        replicate=wildcards.replicate,
        chunk=glob_wildcards(os.path.join(checkpoint_output, "{chunk}.tsv")).chunk,
    )


rule build_swap_file:
    input:
        get_chunk_swap_files,
    output:
        config.chisel_swap_file,
    log:
        config.get_log_file(config.chisel_swap_file),
    benchmark:
        config.get_benchmark_file(config.chisel_swap_file, config.chisel_benchmark_dir)
    resources:
        mem="16G",
    threads: config.chisel_sort_threads
    shell:
        "(sort -m -k 1,1n -k 2,2n --parallel {threads} -S 8G -o {output} {input}) >{log} 2>&1"


rule run_chisel_call:
    input:
        config.chisel_combo_file,
    output:
        config.chisel_calls_file,
        config.chisel_clones_file,
    params:
        o=directory(config.chisel_results_dir),
        p=config.chisel_max_ploidy,
    conda:
        "envs/chisel.yaml"
    benchmark:
        config.chisel_benchmark_dir.joinpath("chisel_calling.txt")
    log:
        config.chisel_log_dir.joinpath("chisel_calling.log"),
    resources:
        mem="64G",
    threads: config.chisel_call_threads
    shell:
        "(chisel_calling "
        "{input} "
        "-x {params.o} "
        "-P {params.p} "
        "-j {threads} ) >{log} 2>&1"
