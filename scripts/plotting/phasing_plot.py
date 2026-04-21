import pandas as pd
import numpy as np
import matplotlib.pyplot as pp
import h5py

def main(args):

    with h5py.File(args.hapclone_file, "r") as fh:
        hapclone_baf = fh["baf"][()]

    a = hapclone_baf[1, :, :, 0].sum(axis=1)
    b = hapclone_baf[1, :, :, 1].sum(axis=1)
    d = a / (a + b)
    pp.scatter(np.arange(d.shape[0]), d, s=1)
    pp.ylim(0, 1)
    pp.xlim(0, d.shape[0])

    pp.savefig(args.plot_file)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-ha", "--hapclone-file", required=True)

    parser.add_argument("-p", "--plot-file", required=True)

    cli_args = parser.parse_args()

    main(cli_args)