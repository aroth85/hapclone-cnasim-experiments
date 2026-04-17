import numpy as np
import pandas as pd


def main(args):
    df = pd.read_csv(args.in_file, compression="gzip")

    df = df[["chr", "start", "end", "cell_id", "state", "copy"]]

    df.to_csv(args.out_file, index=False, sep="\t")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
