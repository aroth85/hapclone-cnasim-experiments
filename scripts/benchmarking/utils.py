import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import Phylo
from scipy.spatial import distance


def load_cnasim_tree(path):
    tree = Phylo.read(path, "newick")
    nodes = tree.get_terminals()
    leaves = []
    for i in nodes:
        leaves.append(i.name)

    return tree, leaves


def load_cnasim_profile(path):
    profile = pd.read_csv(path, sep="\t")
    profile[["A", "B"]] = profile["CN states"].str.split(r",", expand=True)
    profile["total"] = profile["A"].astype(int) + profile["B"].astype(int)

    return profile


def load_chisel_results(path):
    chisel = pd.read_csv(
        path, sep="\t", usecols=["#CHR", "START", "END", "CELL", "HAP_CN"]
    )
    chisel[["A", "B"]] = chisel["HAP_CN"].str.split(r"|", expand=True)
    chisel["total"] = chisel["A"].astype(int) + chisel["B"].astype(int)

    return chisel


def load_chisel_clones(path):
    chisel = pd.read_csv(path, sep="\t")
    return chisel


def load_signals_results(path):
    signals = pd.read_csv(
        path, sep="\t", usecols=["cell_id", "chr", "start", "state", "BAF"]
    )
    return signals


def load_signals_clones(path):
    signals = pd.read_csv(path, sep="\t")
    signals = signals[signals["clone_id"] != "0"]
    return signals


def load_hmmcopy_results(path):
    hmmcopy = pd.read_csv(path)
    return hmmcopy


def load_hmmcopy_clones(path):
    hmmcopy = pd.read_csv(hmmcopy_file, sep="\t")
    hmmcopy = hmmcopy[hmmcopy["clone_id"] != "0"]
    return hmmcopy


def load_hapclone_results(path, baf=False, adj=False, baf_adj=False):
    hapclone = pd.read_csv(path, sep="\t", compression="gzip")
    hapclone["total"] = hapclone["cn_A_cell"] + hapclone["cn_B_cell"]
    return hapclone


def arrange_hapclone(hapclone_files, num_replicates):
    cols = int(len(hapclone_files) / num_replicates)
    hapclone_files = sorted(hapclone_files, key=lambda x: int(x.split("/")[-5][-1]))
    hapclone_files = np.array(hapclone_files).reshape(num_replicates, cols)

    methods = []
    for file in hapclone_files[0]:
        name = "HapClone - " + file.split("/")[-2]
        methods.append(name)

    return methods, hapclone_files

def save_results(data, columns, path):
    df = pd.DataFrame(data=data, columns=columns)
    df.to_csv(path, sep="\t")

def benchmark_cluster(data, distance_matrix, cluster_col, cell_col):
    within = [] #Distance between all cells within a clone
    maxs = [] #Maximum distance between cells in a clone
    between = [] #Distances from all cells in one clone to all cells not in the clone 
    mins = [] #Nearest distance between cell in a clone and cell not in a clone 

    clones = data[cluster_col].unique()
    for clone in clones:
        cells = data[data[cluster_col] == clone][cell_col].unique() #Get cells in clone
        not_in_clone = data[data[cluster_col] != clone][cell_col].unique() #Cells not in clone
        
        if len(cells) > 1: #Requires that clone has more than 1 cell in it
            within_clone = []
            for i in range(len(cells)):
                for j in range(len(cells)):
                    if (i < j):
                        within_clone.append(distance_matrix[int(cells[i][4:]) - 1, int(cells[j][4:]) -1])
            maxs.append(np.max(within_clone))
            within.append(np.mean(within_clone)) #mean of distances within a clone
        
        if len(not_in_clone) > 0: #Requires that clone isn't all the cells
            between_clone = []
            for i in range(len(cells)):
                for j in range(len(not_in_clone)):
                        between_clone.append(distance_matrix[int(cells[i][4:]) - 1, int(not_in_clone[j][4:]) - 1])
            between.append(np.mean(between_clone))
            mins.append(np.min(between_clone))
            
    return np.mean(within), np.mean(maxs), np.mean(between), np.mean(mins)

def benchmark_copynumber(data, cnasim, cell_col, total_col, cells):
    hamming = []
    ploidy = []
    for i in range(len(cells)):
        profile_cnasim = cnasim[cnasim['CELL'] == cells[i]].total
        profile = data[data[cell_col] == cells[i]][total_col]
        hamming.append(distance.hamming(profile_cnasim, profile))
        ploidy.append(np.mean(profile_cnasim) - np.mean(profile))
    return np.mean(ploidy), np.mean(hamming)
