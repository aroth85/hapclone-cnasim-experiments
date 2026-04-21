import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from utils import *

def main(args):
 
    hapclone = load_hapclone_results(args.hapclone_file, baf=True, adj=True, baf_adj=True)

    _, leaves = load_cnasim_tree(args.tree_file)
    profile = load_cnasim_profile(args.cnasim_file)

    ticks, tick_labels = get_ticks(profile)

    hapclone_adjusted = ordered_profiles('cell_id', 'cn_cell_adj', hapclone, leaves)
    hapclone_profile = ordered_profiles('cell_id', 'total', hapclone, leaves)
    cnasim_profile = ordered_profiles('CELL', 'total', profile, leaves)
    profiles = [cnasim_profile, hapclone_profile, hapclone_adjusted]
    names = ['CNAsim', 'HapClone - Unadjusted', 'HapClone - Adjusted']

    total = plot_profiles(profiles, names, ticks, tick_labels, ["OrRd"], 0, 12, 2)
    total.savefig(args.total_plot)

    hapclone_adjusted = ordered_profiles('cell_id', 'baf_adjusted', hapclone, leaves)
    hapclone_profile = ordered_profiles('cell_id', 'baf', hapclone, leaves)
    cnasim_profile = ordered_profiles('CELL', 'baf', profile, leaves)
    profiles = [cnasim_profile, hapclone_profile, hapclone_adjusted]

    baf = plot_profiles(profiles, names, ticks, tick_labels, ["GnBu"], 0, 1, 2)
    baf.savefig(args.baf_plot)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-ha", "--hapclone-file", required=True)

    parser.add_argument("-p", "--cnasim-file", required=True)

    parser.add_argument("-t", "--tree-file", required=True)

    parser.add_argument("-b", "--baf-plot", required=True)

    parser.add_argument("-tp", "--total-plot", required=True)

    cli_args = parser.parse_args()

    main(cli_args)