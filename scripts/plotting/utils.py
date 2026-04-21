import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import Phylo


def ordered_profiles(cell_col, total_col, data, order):
    profiles = []
    for j in range(len(order)):
        cell = order[j]
        profiles.append(data.loc[data[cell_col] == cell, total_col].values)
    profiles = np.vstack(profiles)
    return profiles


def load_cnasim_profile(path):
    profile = pd.read_csv(path, sep="\t")
    profile[["A", "B"]] = profile["CN states"].str.split(r",", expand=True)
    profile["total"] = profile["A"].astype(int) + profile["B"].astype(int)
    profile["baf"] = profile["A"].astype(int) / profile["total"]

    return profile


def load_cnasim_tree(path):
    tree = Phylo.read(path, "newick")
    nodes = tree.get_terminals()
    leaves = []
    for i in nodes:
        leaves.append(i.name)

    return tree, leaves


def load_hapclone_results(path, baf=False, adj=False, baf_adj=False):
    hapclone = pd.read_csv(path, sep="\t", compression="gzip")
    hapclone["total"] = hapclone["cn_A_cell"] + hapclone["cn_B_cell"]
    if baf == True:
        hapclone["baf"] = hapclone["cn_A_cell"] / hapclone["total"]
    if adj == True:
        hapclone['cn_cell_adj'] = hapclone['cn_A_cell_adj'] + hapclone['cn_B_cell_adj']
    if baf_adj == True:
        hapclone['baf_adjusted'] = hapclone['cn_A_cell_adj'] / (hapclone['cn_cell_adj'])

    return hapclone


def load_chisel_results(path, baf=False):
    chisel = pd.read_csv(
        path, sep="\t", usecols=["#CHR", "START", "END", "CELL", "HAP_CN"]
    )
    chisel[["A", "B"]] = chisel["HAP_CN"].str.split(r"|", expand=True)
    chisel["total"] = chisel["A"].astype(int) + chisel["B"].astype(int)
    if baf == True:
        chisel["baf"] = chisel["A"].astype(int) / chisel["total"]
    chisel = chisel.sort_values(["CELL", "#CHR", "START"])

    return chisel

def load_hmmcopy_results(path):
    hmmcopy = pd.read_csv(path, sep="\t")

    return hmmcopy

def load_signals_results(path):
    signals = pd.read_csv(
    signals_file, sep="\t", usecols=["cell_id", "chr", "start", "state", "BAF"]
    )
    signals = signals.sort_values(["cell_id", "chr", "start"])

    return signals


def get_ticks(cnasim):
    cell = cnasim[cnasim["CELL"] == "cell1"].reset_index(drop=True)
    cell["bin"] = cell.index
    bins = cell.drop_duplicates(subset="chrom")
    ticks = bins.bin.values
    tick_labels = [x[3:] for x in bins.chrom.values]

    return ticks, tick_labels


def plot_profiles(profiles, names, ticks, tick_labels, cmap, vmin, vmax, size):

    if len(cmap) == 1:
        cmap = cmap * len(profiles)

    fig, ax = plt.subplots(size, size, figsize=(20, 16))

    i = 0
    for j in range(size):
        for k in range(size):
            ax_num = (j, k)
            if i >= len(profiles):
                break
            else:
                im = ax[ax_num].imshow(
                    profiles[i],
                    interpolation="none",
                    cmap=cmap[i],
                    vmin=vmin,
                    vmax=vmax,
                )
                ax[ax_num].set_title(names[i])
                ax[ax_num].set_aspect("auto")
                ax[ax_num].set_xlabel("Chromosome")
                ax[ax_num].set_ylabel("Cells")
                ax[ax_num].set_xticks(
                    ticks, labels=tick_labels, fontsize=6, rotation="vertical"
                )
                cbar = ax[ax_num].figure.colorbar(im)
            i += 1

    return fig
