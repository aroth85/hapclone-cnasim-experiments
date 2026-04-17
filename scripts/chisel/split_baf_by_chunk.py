import pandas as pd


def main(args):
    baf = pd.read_csv(args.baf_file, names=["chrom", "coord", "cell", "A", "B"], sep="\t")

    rdr = pd.read_csv(args.rdr_file, names=["chrom", "beg", "end", "cell", "normal_count", "count", "rdr"], sep="\t")
    rdr["beg"] = rdr["beg"].astype(int)
    rdr["end"] = rdr["end"].astype(int)
    # print(rdr.head())
    out_df = []

    for chrom in rdr["chrom"].unique():
        beg_coord = rdr.loc[rdr["chrom"] == chrom, "beg"].min() - args.pad_size
        end_coord = rdr.loc[rdr["chrom"] == chrom, "end"].max() + args.pad_size
        out_df.append(baf[(baf["chrom"] == chrom) & (baf["coord"] >= beg_coord) & (baf["coord"] <= end_coord)])

    out_df = pd.concat(out_df)

    out_df.to_csv(args.out_file, header=False, index=False, sep="\t")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-b", "--baf-file", required=True)

    parser.add_argument("-r", "--rdr-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    parser.add_argument("-p", "--pad-size", type=int, default=1000)

    cli_args = parser.parse_args()

    main(cli_args)
