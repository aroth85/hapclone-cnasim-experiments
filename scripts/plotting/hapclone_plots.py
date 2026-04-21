import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import Phylo
from utils import ordered_profiles, plot_hapclone


def main(args):

    print("Plotting")
    hapclone_files = args.hapclone_file
    cnasim_file = args.cnasim_file
    tree_file = args.tree_file

    print("Reading in ground truth data")
    cnasim = pd.read_csv(cnasim_file, sep="\t")
    cnasim[["A", "B"]] = cnasim["CN states"].str.split(r",", expand=True)
    cnasim["total"] = cnasim["A"].astype(int) + cnasim["B"].astype(int)
    cnasim["baf"] = cnasim["A"].astype(int) / cnasim["total"]

    cell = cnasim[cnasim["CELL"] == "cell1"].reset_index(drop=True)
    cell["bin"] = cell.index
    bins = cell.drop_duplicates(subset="chrom")
    cnasim_ticks = bins.bin.values
    tick_labels = [x[3:] for x in bins.chrom.values]

    tree = Phylo.read(tree_file, "newick")
    nodes = tree.get_terminals()
    leaves = []
    for i in nodes:
        leaves.append(i.name)

    profiles = []
    methods = []

    # CNAsim
    cnasim_profile = ordered_profiles("CELL", "total", cnasim, leaves)
    profiles.append(cnasim_profile)
    methods.append("cnasim")

    print("Reading in results")
    for i in range(len(hapclone_files)):
        try:
            hapclone_file = str(hapclone_files[i])
            hapclone = pd.read_csv(hapclone_file, sep="\t", compression="gzip")
            hapclone["total"] = hapclone["cn_A_cell"] + hapclone["cn_B_cell"]

            # plotting info
            hapclone_profile = ordered_profiles("cell_id", "total", hapclone, leaves)
            profiles.append(hapclone_profile)
            name = str(hapclone_files[i]).split("/")[-2]
            methods.append(name)

        except FileNotFoundError:
            hapclone = []
            print("No HapClone results found")

    print("Plotting results")
    fig = plot_hapclone(profiles, methods, cnasim_ticks, tick_labels, "OrRd", 0, 12)
    fig.savefig(args.plot_file)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-ha", "--hapclone-file", nargs="+", required=True)

    parser.add_argument("-p", "--cnasim-file", required=True)

    parser.add_argument("-t", "--tree-file", required=True)

    parser.add_argument("-o", "--plot-file", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
