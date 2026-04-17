import os
import pandas as pd
import pathlib
import warnings

warnings.filterwarnings("ignore")


def main(args):
    reads = pd.read_csv(args.in_file, converters={"chrom": str}, delimiter="\t")

    reads = reads.rename(columns={"CELL": "cell_id", "chrom": "chr"})

    reads = reads[reads["cell_id"] == args.cell_id]

    reads["chr"] = reads["chr"].str.replace("chr", "")

    reads["ideal"] = True

    reads["reads"] = reads["Acount"] + reads["Bcount"]

    reads = reads.drop(columns=["Acount", "Bcount"])

    # Average reads across bin
    mean = reads.reads.mean()

    # Create new column with reads/mean
    reads["cor_gc"] = reads["reads"] / mean

    reads["copy"] = (reads["cor_gc"] - reads["cor_gc"].mean()) + 1

    reads = reads.drop(columns=["cell_id"])

    reads["sample_id"] = str(args.sample_id)

    reads.to_csv(args.out_file, sep="\t", index=False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    parser.add_argument("-c", "--cell-id", required=True)

    parser.add_argument("-s", "--sample-id", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
