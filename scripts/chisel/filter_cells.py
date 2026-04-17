import pandas as pd


def main(args):

    df = pd.read_csv(args.in_file, header=None, sep="\t")

    df.columns = "cell_id", "count"
    df = df[df["count"] >= args.min_reads]

    df.to_csv(args.out_file, header=False, index=False, sep="\t")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    parser.add_argument("-r", "--min-reads", type=int, default=int(1e5))

    cli_args = parser.parse_args()

    main(cli_args)
