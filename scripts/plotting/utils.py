import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def ordered_profiles(cell_col, total_col, data, order):
    profiles = []
    for j in range(len(order)):
        cell = order[j]
        profiles.append(data.loc[data[cell_col] == cell, total_col].values)
    profiles = np.vstack(profiles)
    return profiles


def plot_hapclone(profiles, names, ticks, tick_labels, cmap, vmin, vmax):

    fig, ax = plt.subplots(4, 4, figsize=(24, 20))
    order = [
        (0, 0),
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 0),
        (1, 1),
        (1, 2),
        (1, 3),
        (2, 0),
        (2, 1),
        (2, 2),
        (2, 3),
        (3, 0),
        (3, 1),
        (3, 2),
        (3, 3),
    ]

    for i in range(1, len(profiles)):
        ax_num = order[i]
        im = ax[ax_num].imshow(
            profiles[i], interpolation="none", cmap=cmap, vmin=vmin, vmax=vmax
        )
        ax[ax_num].set_title(names[i])
        ax[ax_num].set_aspect("auto")
        ax[ax_num].set_xlabel("Chromosome")
        ax[ax_num].set_ylabel("Cells")
        ax[ax_num].set_xticks(
            ticks, labels=tick_labels, fontsize=6, rotation="vertical"
        )
        cbar = ax[ax_num].figure.colorbar(im)

    # cnasim
    im = ax[0, 0].imshow(
        profiles[0], interpolation="none", cmap=cmap, vmin=vmin, vmax=vmax
    )
    ax[0, 0].set_title(names[0])
    ax[0, 0].set_aspect("auto")
    ax[0, 0].set_xlabel("Chromosome")
    ax[0, 0].set_ylabel("Cells")
    ax[0, 0].set_xticks(ticks, labels=tick_labels, fontsize=6, rotation="vertical")
    cbar = ax[0, 0].figure.colorbar(im)

    return fig
