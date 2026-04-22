rule plot_data:
    conda:
        "envs/benchmarking.yaml"
    input:
        ha=config.hapclone_input_file,
        p=config.cnasim_reads,
        t=config.cnasim_tree
    output:
        config.data_plot
    shell:
       "python scripts/plotting/data_plot.py -ha {input.ha} -p {input.p} -t {input.t} -o {output}"

rule plot_total:
    conda:
        "envs/benchmarking.yaml"
    params:
        c=config.chisel_calls_file,
        s=config.signals_output_template,
        hm=config.hmmcopy_reads_template
    input:
        ha=config.hapclone_default,
        p=config.cnasim_profiles,
        t=config.cnasim_tree
    output:
        config.total_plot
    shell:
        "python scripts/plotting/total_plots.py -ha {input.ha} -p {input.p} -t {input.t} -s {params.s} -hm {params.hm} -c {params.c} -o {output}"

rule plot_baf:
    conda:
        "envs/benchmarking.yaml"
    params:
        c=config.chisel_calls_file,
        s=config.signals_output_template
    input:
        ha=config.hapclone_default,
        p=config.cnasim_profiles,
        t=config.cnasim_tree
    output:
        b=config.baf_plot,
        bm=config.baf_mirror_plot
    shell:
        "python scripts/plotting/baf_plot.py -ha {input.ha} -p {input.p} -t {input.t} -s {params.s} -c {params.c} -b {output.b} -bm {output.bm}"

rule plot_hapclone_adjusted:
    conda:
        "envs/benchmarking.yaml"
    input:
        ha=config.hapclone_default,
        p=config.cnasim_profiles,
        t=config.cnasim_tree
    output:
        b=config.hapclone_baf_adj_plot,
        tp=config.hapclone_total_adj_plot
    shell:
        "python scripts/plotting/hapclone_adj_plot.py -ha {input.ha} -p {input.p} -t {input.t} -b {output.b} -tp {output.tp}"
    
rule plot_phasing:
    conda:
        "envs/benchmarking.yaml"
    input:
        ha=config.hapclone_input_file
    output:
        p=config.phasing_plot
    shell:
        "python scripts/plotting/phasing_plot.py -ha {input.ha} -p {output.p}"

rule plot_hapclone_baf:
    conda:
        "envs/benchmarking.yaml"
    input:
        p=config.cnasim_profiles,
        t=config.cnasim_tree,
        ha=expand(config.hapclone_results_file, hapclone_run_config=config.hapclone_cli_args, allow_missing=True)
    output:
        b=config.hapclone_baf_plot,
        bm=config.hapclone_baf_mirror_plot
    shell:
        "python scripts/plotting/hapclone_baf_plot.py -ha {input.ha} -p {input.p} -t {input.t} -b {output.b} -bm {output.bm}"

rule plot_hapclone_total:
    conda:
        "envs/benchmarking.yaml"
    input:
        p=config.cnasim_profiles,
        t=config.cnasim_tree,
        ha=expand(config.hapclone_results_file, hapclone_run_config=config.hapclone_cli_args, allow_missing=True)
    output:
        o=config.hapclone_total_plot
    shell:
        "python scripts/plotting/hapclone_plots.py -ha {input.ha} -p {input.p} -t {input.t} -o {output.o}"

