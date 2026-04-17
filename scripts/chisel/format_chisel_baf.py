import pandas as pd


def main(args):
    df = []

    for file_name in args.in_files:
        cell_df = pd.read_csv(file_name, sep="\t")

        cell_df["total"] = cell_df["a"] + cell_df["b"]

        cell_df = cell_df[cell_df["total"] == 0]

        df.append(cell_df[["chrom", "pos", "cell_id", "a", "b"]])

    df = pd.concat(df)

    df.to_csv(args.out_file, header=False, index=False, sep="\t")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-files", nargs="+", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
