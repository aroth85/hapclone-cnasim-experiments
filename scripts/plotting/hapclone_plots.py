import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import Phylo
from utils import *


def main(args):

    # Loading ground truth data
    cnasim = load_cnasim_profile(args.cnasim_file)
    _, leaves = load_cnasim_tree(args.tree_file)
    ticks, tick_labels = get_ticks(cnasim)

    profiles = []
    methods = []

    # CNAsim
    cnasim_profile = ordered_profiles("CELL", "total", cnasim, leaves)
    profiles.append(cnasim_profile)
    methods.append("cnasim")

    hapclone_files = args.hapclone_file
    for i in range(len(hapclone_files)):
        hapclone = load_hapclone_results(str(hapclone_files[i]))
        hapclone_profile = ordered_profiles("cell_id", "total", hapclone, leaves)
        profiles.append(hapclone_profile)
        name = str(hapclone_files[i]).split("/")[-2]
        methods.append(name)

    size = np.sqrt(len(profiles))
    fig = plot_profiles(
        profiles, methods, ticks, tick_labels, ["OrRd"], 0, 12, int(np.ceil(size))
    )
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
