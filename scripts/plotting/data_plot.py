import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import Phylo
import h5py
import utils


def main(args):

    cnasim_profile_file = args.cnasim_file
    tree_file = args.tree_file
    hapclone_file = args.hapclone_file

    _, leaves = load_cnasim_tree(tree_file)

    cnasim = pd.read_csv(cnasim_profile_file, sep="\t")
    cnasim["total"] = cnasim["Acount"].astype(int) + cnasim["Bcount"].astype(int)
    cnasim["baf"] = cnasim["Acount"].astype(int) / cnasim["total"]
    cnasim["CELL"] = [int(x[4:]) - 1 for x in cnasim.CELL.values]

    ticks, tick_labels = get_ticks(cnasim)

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

    fig = plot_profiles(profiles, names, ticks, tick_labels, cmap, None, None, 2)
    fig.savefig(args.out_file)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-ha", "--hapclone-file", required=True)

    parser.add_argument("-p", "--cnasim-file", required=True)

    parser.add_argument("-t", "--tree-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
