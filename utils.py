import pandas as pd
import pathlib
import os


class ConfigManager(object):
    def __init__(self, config):
        self.config = config

    # Loader methods
    def get_pipeline_files(self):
        for s in self.simulation_set_ids:
            for r in self.replicate_ids:
                for c in self.cell_ids:
                    yield str(self.snp_counts_file).format(cell=c, sim_set=s, replicate=r)
                    #yield str(self.hmmcopy_cnv_profile_template).format(cell=c, sim_set=s, replicate=r)
                    #yield str(self.chisel_calls_file).format(sim_set=s, replicate=r)
                    #yield str(self.chisel_swap_file).format(sim_set=s, replicate=r)
                    yield str(self.hapclone_input_file).format(sim_set=s, replicate=r)
                    #yield str(self.signals_output_template).format(sim_set=s, replicate=r)
                    for h in self.hapclone_cli_args:
                        yield str(self.hapclone_results_file).format(hapclone_run_config=h, sim_set=s, replicate=r)
                    yield str(self.data_plot).format(sim_set=s, replicate=r)
                    yield str(self.total_plot).format(sim_set=s, replicate=r)
                    yield str(self.baf_plot).format(sim_set=s, replicate=r)
                    yield str(self.baf_mirror_plot).format(sim_set=s, replicate=r)
                    yield str(self.hapclone_baf_adj_plot).format(sim_set=s, replicate=r)
                    yield str(self.hapclone_total_adj_plot).format(sim_set=s, replicate=r)
                    yield str(self.phasing_plot).format(sim_set=s, replicate=r)
                    yield str(self.hapclone_baf_plot).format(sim_set=s, replicate=r)
                    yield str(self.hapclone_baf_mirror_plot).format(sim_set=s, replicate=r)
                    yield str(self.hapclone_total_plot).format(sim_set=s, replicate=r)
                    yield str(self.cluster_within_results).format(sim_set=s)
                    yield str(self.cluster_between_results).format(sim_set=s)
                    yield str(self.cluster_max_results).format(sim_set=s)
                    yield str(self.cluster_min_results).format(sim_set=s)
                    yield str(self.ploidy_results).format(sim_set=s)
                    yield str(self.hamming_results).format(sim_set=s)
                    yield str(self.ploidy_plot).format(sim_set=s)
                    yield str(self.hamming_plot).format(sim_set=s)
                    yield str(self.min_plot).format(sim_set=s)
                    yield str(self.max_plot).format(sim_set=s)
                    yield str(self.between_plot).format(sim_set=s)
                    yield str(self.within_plot).format(sim_set=s)

    @property
    def out_dir(self):
        path = pathlib.Path(self.config["outdir"])
        return path.joinpath("sim_{sim_set}", "replicate_{replicate}")

    # Params
    @property
    def bin_size(self):
        return self.config.get("bin_size", int(5e5))

    @property
    def block_size(self):
        return self.config.get("block_size", int(1e4))

    @property
    def cell_ids(self):
        return [f"cell{i + 1}" for i in range(self.num_cells)]

    @property
    def coverage(self):
        return self.config["coverage"]

    @property
    def max_copy_state(self):
        return self.config["max_copystate"]

    @property
    def num_cells(self):
        return self.config["num_cells"]

    @property
    def num_replicates(self):
        return self.config.get("num_replicates", 1)

    @property
    def num_simulations_sets(self):
        return len(self.simulation_set_ids)

    @property
    def num_snps(self):
        return self.config["num_snps"]

    @property
    def phase_switch_prob(self):
        return self.config.get("phase_switch_prob", 0.01)

    @property
    def read_length(self):
        return self.config.get("read_length", 150)

    @property
    def replicate_ids(self):
        return list(range(self.num_replicates))

    @property
    def simulation_set_ids(self):
        return self.cnasim_cli_args.keys()

    # CNASim
    @property
    def cnasim_config(self):
        return self.config.get("cnasim", {})

    @property
    def cnasim_cli_args(self):
        return self.cnasim_config.get("cli_args", {"default": "-M"})

    @property
    def cnasim_num_threads(self):
        return self.cnasim_config.get("num_threads", 1)

    @property
    def cnasim_out_dir(self):
        return self.out_dir.joinpath("cnasim")

    @property
    def cnasim_profiles(self):
        return self.cnasim_out_dir.joinpath("profiles.tsv")

    @property
    def cnasim_reads(self):
        return self.cnasim_out_dir.joinpath("readcounts.tsv")

    @property
    def cnasim_tree(self):
        return self.cnasim_out_dir.joinpath("tree.nwk")

    def get_cna_cli_args(self, wc):
        return self.cnasim_cli_args[wc.sim_set]

    # Phasing files
    @property
    def phasing_dir(self):
        return self.out_dir.joinpath("phasing")

    @property
    def snps_file(self):
        return self.phasing_dir.joinpath("snps.tsv")

    @property
    def snp_counts_file(self):
        return self.phasing_dir.joinpath("snp_counts", "{cell}.tsv.gz")

    # Chisel
    @property
    def chisel_config(self):
        return self.config.get("chisel", {})

    ## Chisel parameters
    @property
    def chisel_call_threads(self):
        return int(self.chisel_config.get("call_threads", 1))

    @property
    def chisel_chunk_size(self):
        return int(self.chisel_config.get("chunk_size", 50))

    @property
    def chisel_combo_restarts(self):
        return int(self.chisel_config.get("combo_restarts", 100))

    @property
    def chisel_combo_threads(self):
        return int(self.chisel_config.get("combo_threads", 1))

    @property
    def chisel_max_ploidy(self):
        return int(self.chisel_config.get("max_ploidy", 4))

    @property
    def chisel_min_cell_reads(self):
        return int(self.chisel_config.get("min_cell_reads", int(1e5)))

    @property
    def chisel_sort_threads(self):
        return int(self.chisel_config.get("sort_threads", 1))

    ## Chisel directories
    @property
    def chisel_dir(self):
        return self.out_dir.joinpath("chisel")

    @property
    def chisel_input_dir(self):
        return self.chisel_dir.joinpath("input")

    @property
    def chisel_results_dir(self):
        return self.chisel_dir.joinpath("results")

    @property
    def chisel_benchmark_dir(self):
        return self.chisel_dir.joinpath("benchmarking")

    @property
    def chisel_log_dir(self):
        return self.chisel_dir.joinpath("log")

    @property
    def chisel_working_dir(self):
        return self.chisel_dir.joinpath("working")

    ## Chisel files
    @property
    def chisel_baf_template(self):
        return self.chisel_input_dir.joinpath("baf.tsv")

    @property
    def chisel_rdr_template(self):
        return self.chisel_input_dir.joinpath("rdr.tsv")

    @property
    def chisel_total_template(self):
        return self.chisel_input_dir.joinpath("total.tsv")

    @property
    def chisel_chunk_baf_template(self):
        return self.chisel_working_dir.joinpath("baf", "{chunk}.tsv")

    @property
    def chisel_chunk_rdr_template(self):
        return self.chisel_working_dir.joinpath("rdr", "{chunk}.tsv")

    @property
    def chisel_chunk_combo_template(self):
        return self.chisel_working_dir.joinpath("combo", "{chunk}.tsv")

    @property
    def chisel_chunk_swap_template(self):
        return self.chisel_working_dir.joinpath("swap", "{chunk}.tsv")

    @property
    def chisel_calls_file(self):
        return self.chisel_results_dir.joinpath("calls", "calls.tsv")

    @property
    def chisel_clones_file(self):
        return self.chisel_results_dir.joinpath("clones", "mapping.tsv")

    @property
    def chisel_combo_file(self):
        return self.chisel_results_dir.joinpath("combo", "combo.tsv")

    @property
    def chisel_swap_file(self):
        return self.chisel_results_dir.joinpath("swap", "swap.tsv")

    @property
    def chisel_total_filtered_file(self):
        return self.chisel_results_dir.joinpath("rdr", "total_filtered.tsv")

    # HMMCopy
    @property
    def hmmcopy_dir(self):
        return self.out_dir.joinpath("hmmcopy")

    @property
    def hmmcopy_input_dir(self):
        return self.hmmcopy_dir.joinpath("input")

    @property
    def hmmcopy_results_dir(self):
        return self.hmmcopy_dir.joinpath("results")

    @property
    def hmmcopy_benchmark_dir(self):
        return self.hmmcopy_dir.joinpath("benchmarking")

    @property
    def hmmcopy_log_dir(self):
        return self.hmmcopy_dir.joinpath("log")

    @property
    def hmmcopy_cnv_profile_template(self):
        return self.hmmcopy_results_dir.joinpath("plots", "{cell}.png")

    @property
    def hmmcopy_input_template(self):
        return self.hmmcopy_input_dir.joinpath("{cell}.tsv.gz")

    @property
    def hmmcopy_metrics_template(self):
        return self.hmmcopy_results_dir.joinpath("metrics.csv.gz")

    @property
    def hmmcopy_output_template(self):
        return self.hmmcopy_results_dir.joinpath("{cell}.h5")

    @property
    def hmmcopy_paths_template(self):
        return self.hmmcopy_results_dir.joinpath("paths.tsv")

    @property
    def hmmcopy_reads_template(self):
        return self.hmmcopy_results_dir.joinpath("reads.csv.gz")

    @property
    def hmmcopy_segs_template(self):
        return self.hmmcopy_results_dir.joinpath("segs.csv.gz")
    
    @property
    def hmmcopy_umap_template(self):
        return self.hmmcopy_results_dir.joinpath("umap.tsv")

    # HapClone
    @property
    def hapclone_config(self):
        return self.config.get("hapclone", {})

    ## Params
    @property
    def hapclone_cli_args(self):
        return self.hapclone_config.get("cli_args", {"default", ""})

    @property
    def hapclone_num_threads(self):
        return self.hapclone_config.get("num_threads", 1)

    ## Directories
    @property
    def hapclone_dir(self):
        return self.out_dir.joinpath("hapclone")

    @property
    def hapclone_input_dir(self):
        return self.hapclone_dir.joinpath("input")

    @property
    def hapclone_results_dir(self):
        return self.hapclone_dir.joinpath("results", "{hapclone_run_config}")

    ## Files
    @property
    def hapclone_cell_input_file(self):
        return self.hapclone_input_dir.joinpath("cells", "{cell}.h5")

    @property
    def hapclone_input_file(self):
        return self.hapclone_input_dir.joinpath("data.h5")

    @property
    def hapclone_fit_file(self):
        return self.hapclone_results_dir.joinpath("fit.h5")

    @property
    def hapclone_results_file(self):
        return self.hapclone_results_dir.joinpath("results.tsv.gz")

    def get_hapclone_cli_args(self, wc):
        return self.hapclone_cli_args[wc.hapclone_run_config]

    # Signals
    ## Directories
    @property
    def signals_dir(self):
        return self.out_dir.joinpath("signals")

    @property
    def signals_input_dir(self):
        return self.signals_dir.joinpath("input")

    @property
    def signals_results_dir(self):
        return self.signals_dir.joinpath("results")

    ## Files
    @property
    def signals_baf_template(self):
        return self.signals_input_dir.joinpath("baf.tsv.gz")

    @property
    def signals_cnv_template(self):
        return self.signals_input_dir.joinpath("hmmcopy.tsv.gz")

    @property
    def signals_output_template(self):
        return self.signals_results_dir.joinpath("signals.tsv.gz")

    @property
    def signals_umap_template(self):
        return self.signals_results_dir.joinpath("umap.tsv")
    
    # Plotting
    @property
    def plot_dir(self):
        path = pathlib.Path(self.config["outdir"])
        return path.joinpath("sim_{sim_set}", "plots")

    @property
    def data_plot(self):
        return self.plot_dir.joinpath("data", "replicate{replicate}.png")

    @property
    def total_plot(self):
        return self.plot_dir.joinpath("total", "replicate{replicate}.png")

    @property
    def baf_plot(self):
        return self.plot_dir.joinpath("baf", "replicate{replicate}.png")
    
    @property
    def baf_mirror_plot(self):
        return self.plot_dir.joinpath("baf_mirror", "replicate{replicate}.png")

    @property
    def hapclone_baf_adj_plot(self):
        return self.plot_dir.joinpath("adjusted", "baf_replicate{replicate}.png")
    
    @property
    def hapclone_total_adj_plot(self):
        return self.plot_dir.joinpath("adjusted", "replicate{replicate}.png")
    
    @property
    def phasing_plot(self):
        return self.plot_dir.joinpath("data", "phasing_{replicate}.png")
    
    @property
    def hapclone_baf_plot(self):
        return self.plot_dir.joinpath("baf", "hapclone", "replicate{replicate}.png")

    @property
    def hapclone_baf_mirror_plot(self):
        return self.plot_dir.joinpath("baf_mirror", "hapclone", "replicate{replicate}.png")

    @property
    def hapclone_total_plot(self):
        return self.plot_dir.joinpath("total", "hapclone", "replicate{replicate}.png")

    @property
    def ploidy_plot(self):
        return self.plot_dir.joinpath("results", "ploidy.png")
    
    @property
    def hamming_plot(self):
        return self.plot_dir.joinpath("results", "hamming.png")

    @property
    def max_plot(self):
        return self.plot_dir.joinpath("results", "max.png")

    @property
    def min_plot(self):
        return self.plot_dir.joinpath("results", "min.png")
    
    @property 
    def between_plot(self):
        return self.plot_dir.joinpath("results", "between.png")
    
    @property
    def within_plot(self):
        return self.plot_dir.joinpath("results", "within.png")
    
    # Benchmarking
    @property
    def benchmark_dir(self):
        path = pathlib.Path(self.config["outdir"])
        return path.joinpath("sim_{sim_set}", "benchmark")

    @property
    def cluster_within_results(self):
        return self.benchmark_dir.joinpath("clustering", "within.tsv")
    
    @property
    def cluster_max_results(self):
        return self.benchmark_dir.joinpath("clustering", "maxs.tsv")
    
    @property
    def cluster_min_results(self):
        return self.benchmark_dir.joinpath("clustering", "mins.tsv")

    @property
    def cluster_between_results(self):
        return self.benchmark_dir.joinpath("clustering", "between.tsv")

    @property
    def ploidy_results(self):
        return self.benchmark_dir.joinpath("copynumber", "ploidy.tsv")

    @property
    def hamming_results(self):
        return self.benchmark_dir.joinpath("copynumber", "hamming.tsv")
    
    @property
    def hapclone_default(self):
        mode = self.config["hapclone_mode"]
        path = self.out_dir.joinpath("hapclone/results/" + mode + "/results.tsv.gz")
        return path

    # Auxiliary files
    def get_log_file(self, template, directory=None):
        file = pathlib.Path(template)
        parent = file.parent
        name = str(file.stem)
        if directory is not None:
            return directory.joinpath(parent.stem, name + ".log")
        else:
            return parent.joinpath("log", name + ".log")

    def get_benchmark_file(self, template, directory=None):
        file = pathlib.Path(template)
        parent = file.parent
        name = str(file.stem)
        if directory is not None:
            return directory.joinpath(parent.stem, name + ".txt")
        return parent.joinpath("benchmark", "benchmark.txt")
