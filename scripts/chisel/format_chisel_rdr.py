import pandas as pd


def main(args):
    rdr = pd.read_csv(args.in_file, delimiter="\t")

    rdr = rdr.rename(columns={"CELL": "cell"})

    rdr["count"] = rdr["Acount"] + rdr["Bcount"]

    rdr = rdr.drop(columns=["Acount", "Bcount"])

    rdr["cell"] = rdr["cell"].str.replace("cell", "").astype(int)

    rdr["chrom"] = rdr["chrom"].str.replace("chr", "")

    rdr["norm_count"] = 1

    rdr["rdr"] = rdr.groupby("cell")["count"].transform(lambda x: x / x.mean())

    total = rdr.groupby("cell")["count"].sum().reset_index()

    rdr = rdr[["chrom", "start", "end", "cell", "norm_count", "count", "rdr"]]

    rdr = rdr.sort_values(by=["cell", "chrom", "start"])

    rdr.to_csv(args.out_rdr_file, header=False, index=False, sep="\t")

    total.to_csv(args.out_total_file, header=False, index=False, sep="\t")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-r", "--out-rdr-file", required=True)

    parser.add_argument("-t", "--out-total-file", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
