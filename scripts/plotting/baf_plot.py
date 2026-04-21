import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import Phylo
from utils import ordered_profiles


def main(args):

    print("Plotting")
    hapclone_file = args.hapclone_file
    chisel_file = args.chisel_file
    signals_file = args.signals_file
    cnasim_file = args.cnasim_file
    tree_file = args.tree_file

    print("Reading in ground truth data")
    cnasim = pd.read_csv(cnasim_file, sep="\t")
    cnasim[["A", "B"]] = cnasim["CN states"].str.split(r",", expand=True)
    cnasim["total"] = cnasim["A"].astype(int) + cnasim["B"].astype(int)
    cnasim["baf"] = cnasim["A"].astype(int) / cnasim["total"]

    # Get ticks
    cell = cnasim[cnasim["CELL"] == 0].reset_index(drop=True)
    cell["bin"] = cell.index
    bins = cell.drop_duplicates(subset="chrom")
    cnasim_ticks = bins.bin.values
    tick_labels = [x[3:] for x in bins.chrom.values]

    tree = Phylo.read(tree_file, "newick")
    nodes = tree.get_terminals()
    leaves = []
    for i in nodes:
        leaves.append(i.name)

    ticks = []
    methods = []
    profiles = []
    profiles_mirror = []

    # CNAsim
    cnasim_profile = ordered_profiles("CELL", "baf", cnasim, leaves)
    methods.append("cnasim")
    ticks.append(cnasim_ticks)
    profiles.append(cnasim_profile)
    profiles_mirror.append(cnasim_profile)

    print("Reading in results")
    try:
        hapclone = pd.read_csv(hapclone_file, sep="\t", compression="gzip")
        hapclone["total"] = hapclone["cn_A_cell"] + hapclone["cn_B_cell"]
        hapclone["baf"] = hapclone["cn_A_cell"] / hapclone["total"]

        # ticks
        cell = hapclone[hapclone["cell_id"] == 0].reset_index(drop=True)
        cell["bin"] = cell.index
        bins = cell.drop_duplicates(subset="chrom")
        hapclone_ticks = bins.bin.values

        # plotting info
        hapclone_profile = ordered_profiles("cell_id", "baf", hapclone, leaves)
        profiles.append(hapclone_profile)
        profiles_mirror.append(1 - hapclone_profile)
        methods.append("HapClone")
        ticks.append(hapclone_ticks)

    except FileNotFoundError:
        hapclone = []
        print("No HapClone results found")
    try:
        chisel = pd.read_csv(chisel_file, sep="\t")
        chisel[["A", "B"]] = chisel["HAP_CN"].str.split(r"|", expand=True)
        chisel["total"] = chisel["A"].astype(int) + chisel["B"].astype(int)
        chisel["baf"] = chisel["A"].astype(int) / chisel["total"]
        chisel = chisel.sort_values(["CELL", "#CHR", "START"])

        # Plotting info
        chisel_profile = ordered_profiles("CELL", "baf", chisel, leaves)
        profiles.append(chisel_profile)
        profiles_mirror.append(1 - chisel_profile)
        methods.append("Chisel")
        ticks.append(hapclone_ticks)
    except:  # FileNotFoundError:
        chisel = []
        print("No chisel results found")
    try:
        signals = pd.read_csv(signals_file, sep="\t")

        # Plotting info
        signals_profile = ordered_profiles("cell_id", "BAF", signals, leaves)
        profiles.append(signals_profile)
        profiles_mirror.append(1 - signals_profile)
        methods.append("Signals")
        ticks.append(cnasim_ticks)
    except FileNotFoundError:
        signals = []
        print("No signals results found")

    print("Plotting results")
    # Regular plots
    fig = plot_profiles(profiles, methods, ticks, tick_labels, "GnBu", 0, 1)
    fig.savefig(args.baf_mirror)

    # mirror plots
    fig = plot_profiles(profiles_mirror, methods, ticks, tick_labels, "GnBu", 0, 1)
    fig.savefig(args.baf)

def plot_profiles(profiles, names, ticks, tick_labels, cmap, vmin, vmax):

    fig, ax = plt.subplots(2, 2, figsize=(20, 16))
    order = [(0, 0), (0, 1), (1, 0), (1, 1)]
    for i in range(1, len(profiles)):
        ax_num = order[i]
        im = ax[ax_num].imshow(
            profiles[i], interpolation="none", cmap=cmap, vmin=vmin, vmax=vmax
        )
        ax[ax_num].set_title(names[i])
        ax[ax_num].set_aspect("auto")
        ax[ax_num].set_xlabel("Chromosome")
        ax[ax_num].set_ylabel("Cells")
        ax[ax_num].set_xticks(
            ticks[i], labels=tick_labels, fontsize=6, rotation="vertical"
        )
        cbar = ax[ax_num].figure.colorbar(im)

    # cnasim
    im = ax[0, 0].imshow(profiles[0], interpolation="none", cmap="GnBu", vmin=0, vmax=1)
    ax[0, 0].set_title(names[0])
    ax[0, 0].set_aspect("auto")
    ax[0, 0].set_xlabel("Chromosome")
    ax[0, 0].set_ylabel("Cells")
    ax[0, 0].set_xticks(ticks[0], labels=tick_labels, fontsize=6, rotation="vertical")
    cbar = ax[0, 0].figure.colorbar(im)

    return fig


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-ha", "--hapclone-file", required=True)

    parser.add_argument("-p", "--cnasim-file", required=True)

    parser.add_argument("-t", "--tree-file", required=True)

    parser.add_argument("-s", "--signals-file", required=True)

    parser.add_argument("-c", "--chisel-file", required=True)

    parser.add_argument("-b", "--baf", required=True)

    parser.add_argument("-bm", "--baf-mirror", required=True)

    cli_args = parser.parse_args()

    main(cli_args)