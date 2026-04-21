import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from utils import *


def main(args):

    # Load ground truth data
    cnasim = load_cnasim_profile(args.cnasim_file)
    _, leaves = load_cnasim_tree(args.tree_file)
    ticks, tick_labels = get_ticks(cnasim)

    methods = []
    profs = []
    profs_mirror = []

    # CNAsim
    cnasim_profile = ordered_profiles("CELL", "baf", cnasim, leaves)
    methods.append("cnasim")
    profs.append(cnasim_profile)
    profs_mirror.append(cnasim_profile)

    hapclone_files = args.hapclone_file
    for file in hapclone_files:
        hapclone = load_hapclone_results(str(file), baf=True)
        hapclone_profile = ordered_profiles("cell_id", "baf", hapclone, leaves)
        profs.append(hapclone_profile)
        profs_mirror.append(1 - hapclone_profile)
        name = str(file).split("/")[-2]
        methods.append(name)

    size = np.sqrt(len(profs))
    fig = plot_profiles(
        profs, methods, ticks, tick_labels, ["GnBu"], 0, 1, int(np.ceil(size))
    )
    fig.savefig(args.baf_plot)

    fig = plot_profiles(
        profs_mirror, methods, ticks, tick_labels, ["GnBu"], 0, 1, int(np.ceil(size))
    )
    fig.savefig(args.baf_mirror_plot)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-ha", "--hapclone-file", required=True, nargs="+")

    parser.add_argument("-p", "--cnasim-file", required=True)

    parser.add_argument("-t", "--tree-file", required=True)

    parser.add_argument("-b", "--baf-plot", required=True)

    parser.add_argument("-bm", "--baf-mirror-plot", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
