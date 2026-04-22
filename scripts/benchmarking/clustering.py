import numpy as np
from Bio import Phylo
from pathlib import Path
from utils import *


def main(args):

    results_within2 = []
    maxs2 = []
    mins2 = []
    results_between2 = []

    methods, hapclone_files = arrange_hapclone(args.hapclone_files, args.num_replicates)
    methods.extend(["Chisel", "Signals", "HMMCopy"])

    for f in range(args.num_replicates):

        chisel_file = args.chisel_files[f] 
        signals_file = args.signals_files[f]
        hmmcopy_file = args.hmmcopy_files[f]
        hapclone_file = hapclone_files[f]

        tree, leaves = load_cnasim_tree(args.tree_files[f])
        leaves = sorted(leaves, key=lambda x: int(x[4:]))
        distance_matrix = np.zeros((len(leaves), len(leaves)))

        for i, leaf1 in enumerate(leaves):
            for j, leaf2 in enumerate(leaves):
                if i < j:
                    distance = tree.distance(leaf1, leaf2)
                    distance_matrix[i][j] = distance
                    distance_matrix[j][i] = distance

        mins = []
        maxs = []
        results_between = []
        results_within = []

        # HapClone
        for file in hapclone_file:
            hapclone = load_hapclone_results(file)
            hapclone = hapclone.drop_duplicates(subset="cell_id")
            hapclone_results = benchmark_cluster(
                hapclone, distance_matrix, "cluster_id", "cell_id"
            )
            results_within.append(hapclone_results[0])
            maxs.append(hapclone_results[1])
            results_between.append(hapclone_results[2])
            mins.append(hapclone_results[3])

        # Chisel
        if Path(chisel_file).exists():
            chisel = load_chisel_results(chisel_file)
            chisel_results = benchmark_cluster(
                chisel, distance_matrix, "CLUSTER", "#CELL"
            )
            results_within.append(chisel_results[0])
            maxs.append(chisel_results[1])
            results_between.append(chisel_results[2])
            mins.append(chisel_results[3])
        else:
            mins.append([])
            maxs.append([])
            results_between.append([])
            results_within.append([])

        # Signals
        if Path(signals_file).exists():
            signals = load_signals_clones(signals_file)
            signals_results = benchmark_cluster(
                signals, distance_matrix, "clone_id", "cell_id"
            )
            results_within.append(signals_results[0])
            maxs.append(signals_results[1])
            results_between.append(signals_results[2])
            mins.append(signals_results[3])
        else:
            mins.append([])
            maxs.append([])
            results_between.append([])
            results_within.append([])

        # HmmCopy
        if Path(hmmcopy_file).exists():
            hmmcopy = load_hmmcopy_clonse(hmmcopy_file)
            hmmcopy_results = benchmark_cluster(
                hmmcopy, distance_matrix, "clone_id", "cell_id"
            )
            results_within.append(hmmcopy_results[0])
            maxs.append(np.mean(hmmcopy_results[1]))
            results_between.append(hmmcopy_results[2])
            mins.append(np.mean(hmmcopy_results[3]))
        else:
            mins.append([])
            maxs.append([])
            results_between.append([])
            results_within.append([])

        results_within2.append(results_within)
        results_between2.append(results_between)
        mins2.append(mins)
        maxs2.append(maxs)

    save_results(results_within2, methods, args.results_within)
    save_results(maxs2, methods, args.maxs)
    save_results(mins2, methods, args.mins)
    save_results(results_between2, methods, args.results_between)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-ha", "--hapclone-files", nargs="+", required=True)

    parser.add_argument("-c", "--chisel-files", nargs="+", required=True)

    parser.add_argument("-s", "--signals-files", nargs="+", required=True)

    parser.add_argument("-t", "--tree-files", nargs="+", required=True)

    parser.add_argument("-hm", "--hmmcopy-files", nargs="+", required=True)

    parser.add_argument("-rw", "--results-within", required=True)

    parser.add_argument("-ma", "--maxs", required=True)

    parser.add_argument("-mi", "--mins", required=True)

    parser.add_argument("-rb", "--results-between", required=True)

    parser.add_argument("-n", "--num-replicates", required=True, type=int)

    cli_args = parser.parse_args()

    main(cli_args)
