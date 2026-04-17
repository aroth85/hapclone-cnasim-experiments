import pandas as pd
import pathlib


def main(args):
    out_dir = pathlib.Path(args.out_dir)

    if not out_dir.exists():
        out_dir.mkdir(parents=True)

    df = pd.read_csv(args.in_file, header=None, sep="\t")

    df.columns = ["chrom", "beg", "end", "cell", "normal_count", "count", "rdr"]

    df_split = []

    i = 0

    j = 0

    for _, x in df.groupby(["chrom", "beg", "end"], as_index=False):
        if i >= args.chunk_size:
            out_file = out_dir.joinpath("{}.tsv".format(j))

            df_split = pd.concat(df_split)
            df_split.to_csv(out_file, header=False, index=False, sep="\t")

            i = 0

            j += 1

            df_split = []

        i += 1

        df_split.append(x)

    if len(df_split) > 0:
        out_file = out_dir.joinpath("{}.tsv".format(j))

        df_split = pd.concat(df_split)

        df_split.to_csv(out_file, header=False, index=False, sep="\t")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-dir", required=True)

    parser.add_argument("-c", "--chunk-size", type=int, default=50)

    cli_args = parser.parse_args()

    main(cli_args)
