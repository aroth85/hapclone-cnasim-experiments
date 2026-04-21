import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from utils import *


def main(args):

    print("Reading in ground truth data")
    cnasim = load_cnasim_profile(args.cnasim_file)
    _, leaves = load_cnasim_tree(args.tree_file)

    ticks, tick_labels = get_ticks(cnasim)

    profiles = []
    methods = []

    print("Reading in results")
    # CNAsim
    cnasim_profile = ordered_profiles("CELL", "total", cnasim, leaves)
    profiles.append(cnasim_profile)
    methods.append("cnasim")

    # HapClone
    hapclone = load_hapclone_results(args.hapclone_file)
    hapclone_profile = ordered_profiles("cell_id", "total", hapclone, leaves)
    profiles.append(hapclone_profile)
    methods.append("HapClone")

    #Chisel
    if Path(args.chisel_file).exists():
        chisel = load_chisel_results(args.chisel_file)
        chisel_profile = ordered_profiles("CELL", "total", chisel, leaves)
        profiles.append(chisel_profile)
        methods.append("Chisel")

    #Signals
    if Path(args.signals_file).exists():
        signals = load_signals_results(args.signals_file)
        signals_profile = ordered_profiles("cell_id", "state", signals, leaves)
        profiles.append(signals_profile)
        methods.append("Signals")

    #HMMCopy
    if Path(args.hmmcopy_file).exists():
        hmmcopy = load_hmmcopy_results(args.hmmcopy_file)
        profiles.append(hmmcopy_profile)
        methods.append("HmmCopy")

    print("Plotting results")
    fig = plot_profiles(profiles, methods, ticks, tick_labels, ["OrRd"], 0, 12, 3)
    fig.savefig(args.out_file)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-ha", "--hapclone-file", required=True)

    parser.add_argument("-p", "--cnasim-file", required=True)

    parser.add_argument("-t", "--tree-file", required=True)

    parser.add_argument("-s", "--signals-file", required=True)

    parser.add_argument("-hm", "--hmmcopy-file", required=True)

    parser.add_argument("-c", "--chisel-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
