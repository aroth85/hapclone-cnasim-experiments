import pandas as pd
import numpy as np
import h5py


def main(args):
    df = load_cell_data(args.cell_id, args.profiles_file, args.reads_file, read_length=args.read_length)

    snps = pd.read_csv(args.snps_file, sep="\t")

    # Compute number of SNPs per bin
    snp_density = (
        snps.groupby(["chrom", "start", "end"])["pos"].count().reset_index().rename(columns={"pos": "num_snps"})
    )

    df = pd.merge(df, snp_density, on=["chrom", "start", "end"])

    df = df.sort_values(by=["chrom_idx", "start", "end"])

    # Allocate SNP reads to bins
    snp_coverage = df["num_snps"] * df["coverage"]

    df["snp_coverage"] = np.random.multinomial(snp_coverage.round().sum(), np.random.dirichlet(snp_coverage + 1e-6))

    # Compute allele coverage per bin
    df["snp_a_coverage"] = np.random.binomial(df["snp_coverage"], df["baf"])

    df["snp_b_coverage"] = df["snp_coverage"] - df["snp_a_coverage"]

    df = pd.merge(df, snps, on=["chrom", "start", "end"], how="left")

    out_df = []

    # Assign per SNP coverage across a bin
    for _, group in df.groupby(["chrom", "start", "end"]):
        group["a_true"] = np.random.multinomial(
            group["snp_a_coverage"].iloc[0], np.random.dirichlet(np.ones(group.shape[0]))
        )

        group["b_true"] = np.random.multinomial(
            group["snp_b_coverage"].iloc[0], np.random.dirichlet(np.ones(group.shape[0]))
        )

        out_df.append(group)

    out_df = pd.concat(out_df)

    # Unphase the SNPs
    out_df["a"] = out_df["phase"] * out_df["a_true"] + (1 - out_df["phase"]) * out_df["b_true"]

    out_df["b"] = out_df["phase"] * out_df["b_true"] + (1 - out_df["phase"]) * out_df["a_true"]

    out_df = out_df[["cell_id", "chrom", "start", "end", "pos", "phase_set", "a", "b"]]

    out_df.to_csv(args.out_file, index=False, sep="\t")


def load_cell_data(cell_id, profiles_file, reads_file, read_length=150):
    df = pd.merge(
        load_profiles(cell_id, profiles_file), load_reads(cell_id, reads_file), on=["cell_id", "chrom", "start", "end"]
    )
    df["chrom_idx"] = df["chrom"].str.replace("chr", "")
    df["chrom_idx"] = df["chrom_idx"].replace({"X": "23", "Y": "24"})
    df["chrom_idx"] = df["chrom_idx"].astype(int)
    bin_len = df.at[0, "end"] - df.at[0, "start"]
    df["coverage"] = (df["reads"] * read_length) / bin_len
    return df


def load_profiles(cell_id, file_name):
    df = pd.read_csv(file_name, sep="\t")
    df = df.rename(columns={"CELL": "cell_id"})
    df = df[df["cell_id"] == cell_id]
    df["a"] = df["CN states"].str.split(",").str.get(0).astype(int)
    df["b"] = df["CN states"].str.split(",").str.get(1).astype(int)
    df["baf"] = df["a"] / (df["a"] + df["b"])
    df["baf"] = df["baf"].fillna(0)
    df = df[["cell_id", "chrom", "start", "end", "baf"]]
    return df


def load_reads(cell_id, file_name):
    df = pd.read_csv(file_name, delimiter="\t")
    df = df.rename(columns={"Acount": "a", "Bcount": "b", "CELL": "cell_id"})
    df = df[df["cell_id"] == cell_id]
    df["reads"] = df["a"] + df["b"]
    df = df[["cell_id", "chrom", "start", "end", "reads"]]
    return df


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--cell-id", required=True)

    parser.add_argument("-p", "--profiles-file", required=True)

    parser.add_argument("-r", "--reads-file", required=True)

    parser.add_argument("-s", "--snps-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    parser.add_argument("--read-length", default=150, type=int)

    cli_args = parser.parse_args()

    main(cli_args)
