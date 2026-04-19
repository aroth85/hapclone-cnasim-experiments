rule format_signals_baf:
    input:
        expand(config.snp_counts_file, allow_missing=True, cell=config.cell_ids),
    output:
        config.signals_baf_template,
    conda:
        "envs/python.yaml"
    resources:
        mem="64G",
    shell:
        "python scripts/signals/format_signals_baf.py -i {input} -o {output}"


rule format_signal_cnv:
    input:
        config.hmmcopy_reads_template,
    output:
        config.signals_cnv_template,
    conda:
        "envs/python.yaml"
    resources:
        mem="16G",
    shell:
        "python scripts/signals/format_signals_cnv.py -i {input} -o {output}"


rule signals:
    input:
        b=config.signals_baf_template,
        c=config.signals_cnv_template,
    output:
        config.signals_output_template,
    conda:
        "envs/signals.yaml"
    benchmark:
        config.get_benchmark_file(config.signals_output_template)
    log:
        config.get_log_file(config.signals_output_template),
    resources:
        mem="64G"
    shell:
        "(Rscript scripts/signals/run.R {input.c} {input.b} {output}) >{log} 2>&1"
