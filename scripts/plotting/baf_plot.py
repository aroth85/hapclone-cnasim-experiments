import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from utils import *


def main(args):

    cnasim = load_cnasim_profile(args.cnasim_file)
    _, leaves = load_cnasim_tree(args.tree_file)

    # Get ticks
    ticks, tick_labels = get_ticks(cnasim)

    methods = []
    profiles = []
    profiles_mirror = []

    # CNAsim
    cnasim_profile = ordered_profiles("CELL", "baf", cnasim, leaves)
    methods.append("cnasim")
    profiles.append(cnasim_profile)
    profiles_mirror.append(cnasim_profile)

    print("Reading in results")
    # HapClone
    hapclone = load_hapclone_results(args.hapclone_file, baf=True)
    hapclone_profile = ordered_profiles("cell_id", "baf", hapclone, leaves)
    profiles.append(hapclone_profile)
    profiles_mirror.append(1 - hapclone_profile)
    methods.append("HapClone")

    # Chisel
    if Path(args.chisel_file).exists():
        chisel = load_chisel_results(args.chisel_file, baf=True)
        chisel_profile = ordered_profiles("CELL", "baf", chisel,  leaves)
        profiles.append(chisel_profile)
        profiles_mirror.append(1 - chisel_profile)
        methods.append("Chisel")

    # Signals
    if Path(args.signals_file).exists():
        signals = load_signals_results(args.signals_file)
        signals_profile = ordered_profiles("cell_id", "BAF", signals, leaves)
        profiles.append(signals_profile)
        profiles_mirror.append(1 - signals_profile)
        methods.append("Signals")

    print("Plotting results")
    # Regular plots
    fig = plot_profiles(profiles, methods, ticks, tick_labels, ["GnBu"], 0, 1, 2)
    fig.savefig(args.baf_mirror)

    # mirror plots
    fig = plot_profiles(profiles_mirror, methods, ticks, tick_labels, ["GnBu"], 0, 1, 2)
    fig.savefig(args.baf)


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
