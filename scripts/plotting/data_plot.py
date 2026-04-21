import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import Phylo
import h5py
from utils import ordered_profiles


def main(args):

    cnasim_profile_file = args.cnasim_file
    tree_file = args.tree_file
    hapclone_file = args.hapclone_file

    cnasim = pd.read_csv(cnasim_profile_file, sep="\t")
    cnasim["total"] = cnasim["Acount"].astype(int) + cnasim["Bcount"].astype(int)
    cnasim["baf"] = cnasim["Acount"].astype(int) / cnasim["total"]
    cnasim["CELL"] = [int(x[4:]) - 1 for x in cnasim.CELL.values]

    cell = cnasim[cnasim["CELL"] == 0].reset_index(drop=True)
    cell["bin"] = cell.index
    bins = cell.drop_duplicates(subset="chrom")
    ticks = bins.bin.values
    tick_labels = bins.chrom.values

    with h5py.File(hapclone_file, "r") as fh:
        hapclone_baf = fh["baf"][()]

    a = hapclone_baf[:, :, :, 0].sum(axis=2)
    num_bins = a.shape[1]
    num_cells = a.shape[0]
    hapclone = pd.DataFrame(data=a.flatten(), columns=["acounts"])
    cells = np.repeat(np.arange(num_cells), num_bins)
    hapclone["cells"] = cells
    hapclone["bcounts"] = hapclone_baf[:, :, :, 1].sum(axis=2).flatten()
    hapclone["total"] = hapclone["acounts"] + hapclone["bcounts"]
    hapclone["baf"] = hapclone["acounts"] / hapclone["total"]

    tree = Phylo.read(tree_file, "newick")
    nodes = tree.get_terminals()
    leaves = []
    for i in nodes:
        leaves.append(int(i.name[4:]) - 1)

    # Cnasim reads
    cnasim_profile = ordered_profiles("CELL", "total", cnasim, leaves)
    cnasim_baf = ordered_profiles("CELL", "baf", cnasim, leaves)
    hapclone_profile = ordered_profiles("cells", "total", hapclone, leaves)
    hapclone_baf = ordered_profiles("cells", "baf", hapclone, leaves)
    names = [
        "CNAsim total readcounts",
        "CNAsim baf",
        "Unphased total readcounts",
        "Unphased BAF",
    ]
    profiles = [cnasim_profile, cnasim_baf, hapclone_profile, hapclone_baf]
    cmap = ["OrRd", "GnBu", "OrRd", "GnBu"]

    fig = plot_profiles(profiles, names, ticks, tick_labels, cmap, None, None)
    fig.savefig(args.out_file)


def plot_profiles(profiles, names, ticks, tick_labels, cmap, vmin, vmax):

    fig, ax = plt.subplots(2, 2, figsize=(20, 16))
    order = [(0, 0), (0, 1), (1, 0), (1, 1)]
    for i in range(0, len(profiles)):
        ax_num = order[i]
        im = ax[ax_num].imshow(
            profiles[i], interpolation="none", cmap=cmap[i], vmin=vmin, vmax=vmax
        )
        ax[ax_num].set_title(names[i], fontsize=20)
        ax[ax_num].set_aspect("auto")
        ax[ax_num].set_xlabel("Chromosome")
        ax[ax_num].set_ylabel("Cells")
        ax[ax_num].set_xticks(
            ticks, labels=tick_labels, fontsize=6, rotation="vertical"
        )
        cbar = ax[ax_num].figure.colorbar(im)

    return fig


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-ha", "--hapclone-file", required=True)

    parser.add_argument("-p", "--cnasim-file", required=True)

    parser.add_argument("-t", "--tree-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
