import matplotlib.pyplot as pp
import numpy as np
import pandas as pd
import seaborn as sb


def main(args):
    with pd.HDFStore(args.in_file, "r") as store:
        data = store["data"]

        multipliers = store["multipliers"]

        score = []

        for m in multipliers:
            score.append(store["/{}/metrics".format(m)]["scaled_halfiness"])

    idxs = np.argsort(score)

    fig = pp.figure(figsize=(16, 4 * len(multipliers)))

    grid = fig.add_gridspec(len(multipliers), 1, hspace=0.1)

    data["map"] = 1  # ADDED IN

    for i, m in enumerate(idxs):
        with pd.HDFStore(args.in_file, "r") as store:
            bins = store["/{}/bins".format(multipliers[m])]

        bins = pd.merge(bins, data[["chr", "start", "end", "map"]], on=["chr", "start", "end"])

        bins = bins[bins["map"] >= args.min_mappability]

        bins["chr"] = bins["chr"].astype(str).str.replace("chr", "")  # ADDED .astype(str)

        chroms = sort_chroms(bins["chr"].unique())

        chroms_size = bins["chr"].value_counts()

        width_ratios = [chroms_size[x] for x in chroms]

        sub_grid = grid[i].subgridspec(1, len(chroms), width_ratios=width_ratios, wspace=0.05)

        title = "multiplier={0}, scaled_halfiness={1}, ploidy={2}".format(
            multipliers[m], score[m], bins["state"].mean()
        )

        plot_bins(bins, chroms, sub_grid, title=title)

    grid.tight_layout(fig)

    fig.savefig(args.out_file, bbox_inches="tight")


def plot_bins(bins, chroms, grid, y_col="copy", title=None):
    fig = pp.gcf()
    for i, chrom in enumerate(chroms):
        chrom_bins = bins[bins["chr"] == chrom]

        chrom_bins = chrom_bins.sort_values(by=["start"])

        num_bins = chrom_bins.shape[0]

        chrom_bins["idx"] = np.arange(num_bins)

        ax = fig.add_subplot(grid[0, i])

        ax.scatter(np.arange(num_bins), chrom_bins[y_col], s=1)

        sb.despine(ax=ax, offset=10)

        ax.spines["top"].set_visible(False)

        ax.spines["right"].set_visible(False)

        if i != 0:
            ax.spines["left"].set_visible(False)

            ax.set_yticks([])

            ax.set_yticklabels([])

        else:
            ax.tick_params(axis="x", which="major", labelsize=12)

        ax.set_xticks([num_bins / 2])

        ax.set_xticklabels([chrom], fontsize=12)

        ax.set_ylim(0, bins["state"].max() + 2)

        chrom_bins["seg"] = np.cumsum(chrom_bins["state"] != chrom_bins["state"].shift())

        chrom_bins.groupby("seg").apply(
            lambda x: ax.hlines(x["state"].median(), x["idx"].min(), x["idx"].max(), colors="r")
        )

        chrom_bins.groupby("seg").apply(
            lambda x: ax.hlines(x[y_col].median(), x["idx"].min(), x["idx"].max(), colors="k")
        )

    if title is not None:
        ax = fig.add_subplot(grid[:])

        ax.axis("off")

        ax.set_title(title)


def sort_chroms(chroms):
    numeric = []

    string = []

    for c in chroms:
        try:
            numeric.append(int(c))

        except ValueError:
            string.append(c)

    return [str(x) for x in sorted(numeric)] + list(sorted(string))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    parser.add_argument("-m", "--min-mappability", default=0.9, type=float)

    cli_args = parser.parse_args()

    main(cli_args)
