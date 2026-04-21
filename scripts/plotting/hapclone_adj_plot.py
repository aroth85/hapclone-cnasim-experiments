import numpy as np
import pandas as pd
from Bio import Phylo
import matplotlib.pyplot as plt
from utils import ordered_profiles

def main(args):
 
    hapclone = pd.read_csv(args.hapclone_file, sep='\t', compression='gzip')
    tree = Phylo.read(args.tree_file, "newick")

    profile = pd.read_csv(args.cnasim_file, sep='\t')
    profile[['A', 'B']] = profile['CN states'].str.split(r",", expand=True)
    profile['total'] = profile['A'].astype(int) + profile['B'].astype(int)
    profile['baf'] = profile['A'].astype(int)/profile['total']

    nodes = tree.get_terminals()
    leaves = []
    for i in nodes:
        leaves.append(i.name)
    
    cell = profile[profile['CELL'] == 0].reset_index(drop = True)
    cell['bin'] = cell.index
    bins = cell.drop_duplicates(subset='chrom')
    cnasim_ticks = bins.bin.values
    tick_labels = [x[3:] for x in bins.chrom.values]

    cell = hapclone[hapclone['cell_id'] == 0].reset_index(drop = True)
    cell['bin'] = cell.index
    bins = cell.drop_duplicates(subset='chrom')
    hapclone_ticks = bins.bin.values

    hapclone['cn_cell_adj'] = hapclone['cn_A_cell_adj'] + hapclone['cn_B_cell_adj']
    hapclone['cn_total_adj'] = hapclone['cn_A_adj'] + hapclone['cn_B_adj']
    hapclone['total'] = hapclone['cn_A_cell'] + hapclone['cn_B_cell']

    hapclone['baf'] = hapclone['cn_A_cell'] / (hapclone['cn_A_cell'] + hapclone['cn_B_cell'])
    hapclone['baf_adjusted'] = hapclone['cn_A_cell_adj'] / (hapclone['cn_cell_adj'])

    hapclone_adjusted = ordered_profiles('cell_id', 'cn_cell_adj', hapclone, leaves)
    hapclone_profile = ordered_profiles('cell_id', 'total', hapclone, leaves)
    cnasim_profile = ordered_profiles('CELL', 'total', profile, leaves)

    total = plot_profiles(cnasim_profile, hapclone_profile, hapclone_adjusted, cnasim_ticks, hapclone_ticks, tick_labels, 12, 'OrRd')

    total.savefig(args.total_plot)

    hapclone_adjusted = ordered_profiles('cell_id', 'baf_adjusted', hapclone, leaves)
    hapclone_profile = ordered_profiles('cell_id', 'baf', hapclone, leaves)
    cnasim_profile = ordered_profiles('CELL', 'baf', profile, leaves)

    baf = plot_profiles(cnasim_profile, hapclone_profile, hapclone_adjusted, cnasim_ticks, hapclone_ticks, tick_labels, 1, 'GnBu')
    baf.savefig(args.baf_plot)

def plot_profiles(cnasim_profile, hapclone_profile, hapclone_adjusted, cnasim_ticks, hapclone_ticks, tick_labels, vmax, cmap):
        fig, ax = plt.subplots(2, 2, figsize = (16, 12))
        im = ax[0, 0].imshow(cnasim_profile, interpolation= 'none', cmap=cmap, vmin = 0, vmax = vmax)
        ax[0, 0].set_title('CNAsim')
        ax[0, 0].set_aspect('auto')
        ax[0, 0].set_xlabel('Chromosome')
        ax[0, 0].set_ylabel('Cells')
        ax[0, 0].set_xticks(cnasim_ticks, labels = tick_labels, fontsize = 6,  rotation = 'vertical')
        cbar = ax[0, 0].figure.colorbar(im)

        im = ax[0, 1].imshow(hapclone_profile, interpolation= 'none', cmap=cmap, vmin = 0, vmax = vmax)
        ax[0, 1].set_title('HapClone - Unadjusted')
        ax[0, 1].set_aspect('auto')
        ax[0, 1].set_xlabel('Chromosome')
        ax[0, 1].set_ylabel('Cells')
        ax[0, 1].set_xticks(hapclone_ticks, labels = tick_labels, fontsize = 6,  rotation = 'vertical')
        cbar = ax[0, 1].figure.colorbar(im)

        im = ax[1, 0].imshow(hapclone_adjusted, interpolation= 'none', cmap=cmap, vmin = 0, vmax = vmax)
        ax[1, 0].set_title('HapClone - All Adjusted ')
        ax[1, 0].set_aspect('auto')
        ax[1, 0].set_xlabel('Chromosome')
        ax[1, 0].set_ylabel('Cells')
        ax[1, 0].set_xticks(hapclone_ticks, labels = tick_labels, fontsize = 6,  rotation = 'vertical')
        cbar = ax[0, 0].figure.colorbar(im)

        return fig


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