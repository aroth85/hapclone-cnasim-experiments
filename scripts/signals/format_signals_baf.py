import pandas as pd


def main(args):
    header = True

    mode = "w"

    for file_name in args.in_files:
        cell_df = pd.read_csv(file_name, sep="\t")

        cell_id = cell_df.at[0, "cell_id"]

        cell_df = cell_df.groupby(["chrom", "start", "end", "phase_set"])[["a", "b"]].sum().reset_index()

        cell_df["cell_id"] = cell_id

        cell_df["hap_label"] = range(cell_df.shape[0])

        cell_df["chrom"] = cell_df["chrom"].str.replace("chr", "")

        cell_df = cell_df.rename(columns={"chrom": "chr", "a": "allele0", "b": "allele1"})

        cell_df["totalcounts"] = cell_df["allele0"] + cell_df["allele1"]

        cell_df = cell_df[["chr", "start", "end", "cell_id", "hap_label", "allele1", "allele0", "totalcounts"]]

        cell_df.to_csv(args.out_file, index=False, header=header, mode=mode, sep="\t")

        header = False

        mode = "a"


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-files", nargs="+", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
