rule clustering_benchmark:
    conda:
        "envs/benchmarking.yaml"
    input:
        ha=expand(config.hapclone_results_file, hapclone_run_config=config.hapclone_cli_args, replicate=config.replicate_ids, allow_missing=True),
        t=expand(config.cnasim_tree, replicate=config.replicate_ids, allow_missing=True),
    params:
        n=config.num_replicates,
        hm=expand(config.hmmcopy_reads_template, replicate=config.replicate_ids, allow_missing=True),
        c=expand(config.chisel_clones_file, replicate=config.replicate_ids, allow_missing=True),
        s=expand(config.signals_output_template, replicate=config.replicate_ids, allow_missing=True),
    output:
        rw=config.cluster_within_results,
        ma=config.cluster_max_results,
        mi=config.cluster_min_results,
        rb=config.cluster_between_results
    shell:
        "python scripts/benchmarking/clustering.py "
        "-ha {input.ha} "
        "-c {params.c} "
        "-t {input.t} "
        "-s {params.s} "
        "-hm {params.hm} "
        "-n {params.n} "
        "-rw {output.rw} "
        "-ma {output.ma} "
        "-mi {output.mi} "
        "-rb {output.rb} "

rule copynumber_benchmark:
    conda:
        "envs/benchmarking.yaml"
    input:
        ha=expand(config.hapclone_results_file, hapclone_run_config=config.hapclone_cli_args, replicate=config.replicate_ids, allow_missing=True),
        p=expand(config.cnasim_profiles, replicate=config.replicate_ids, allow_missing=True),
    params:
        n=config.num_replicates,
        hm=expand(config.hmmcopy_reads_template, replicate=config.replicate_ids, allow_missing=True),
        c=expand(config.chisel_calls_file, replicate=config.replicate_ids, allow_missing=True),
        s=expand(config.signals_output_template, replicate=config.replicate_ids, allow_missing=True),
    output:
        po=config.ploidy_results,
        ho=config.hamming_results
    resources:
        mem="8G",
    shell:
        "python scripts/benchmarking/copynumber.py "
        "-ha {input.ha} "
        "-c {params.c} "
        "-p {input.p} "
        "-s {params.s} "
        "-hm {params.hm} "
        "-n {params.n} "
        "-po {output.po} "
        "-ho {output.ho} "



